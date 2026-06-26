"""Python mirror of src/spf_sampler.hpp, plus a bit-for-bit port of OpenMC's
PCG-RXS-M-XS ``prn`` (src/random_lcg.cpp). Two jobs:

1. Provide an independent oracle for the C++ sampler (tier-1 parity test:
   C++ standalone driver vs this module, identical seed -> identical samples).
2. Let the statistical tests call a sampler directly from Python.

The RNG draw order here MUST match spf_sampler.hpp exactly (selector, shape-u,
phi) or the parity test will (correctly) fail.
"""
from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Exact port of OpenMC's PRN (src/random_lcg.cpp:8-44). 64-bit modular math.
# ---------------------------------------------------------------------------
_MASK64 = (1 << 64) - 1
_PRN_MULT = 6364136223846793005          # prn_mult
_PRN_ADD = 1442695040888963407           # prn_add
_PERM_MULT = 12605985483714917081        # output permutation multiplier
_DEFAULT_SEED = 1                        # master_seed default
_DEFAULT_STRIDE = 152917                 # DEFAULT_STRIDE
STREAM_SOURCE = 1                        # source-sampling stream offset


class Rng:
    """Mutable single-stream RNG state. ``__call__`` returns a double in [0,1),
    matching ``openmc::prn(uint64_t* seed)``."""

    __slots__ = ("state",)

    def __init__(self, seed: int):
        self.state = seed & _MASK64

    def __call__(self) -> float:
        s = (_PRN_MULT * self.state + _PRN_ADD) & _MASK64
        self.state = s
        word = (((s >> ((s >> 59) + 5)) ^ s) * _PERM_MULT) & _MASK64
        result = (word >> 43) ^ word
        return math.ldexp(result, -64)  # result / 2**64  -> [0,1)


def future_seed(n: int, seed: int) -> int:
    """Port of future_seed (F. Brown O(log n) skip-ahead)."""
    g, c = _PRN_MULT, _PRN_ADD
    g_new, c_new = 1, 0
    n &= _MASK64
    while n > 0:
        if n & 1:
            g_new = (g_new * g) & _MASK64
            c_new = (c_new * g + c) & _MASK64
        c = (c * (g + 1)) & _MASK64
        g = (g * g) & _MASK64
        n >>= 1
    return (g_new * seed + c_new) & _MASK64


def init_seed(particle_id: int, offset: int = STREAM_SOURCE,
              master_seed: int = _DEFAULT_SEED, stride: int = _DEFAULT_STRIDE) -> int:
    """Port of init_seed(id, offset): the seed OpenMC uses for source sampling."""
    return future_seed((particle_id * stride) & _MASK64, (master_seed + offset) & _MASK64)


# ---------------------------------------------------------------------------
# Mode weights -- mirror of spf::make_mode_weights.
# ---------------------------------------------------------------------------
PI = math.pi
INT_SIN2 = 8.0 * PI / 3.0   # ∫ sin²θ dΩ
INT_PAR = 2.0 * PI          # ∫ (1/4 + 3/4 cos²θ) dΩ


class ModeWeights:
    __slots__ = ("a", "b", "c", "w_perp", "w_par", "eta", "P_perp")

    def __init__(self, a: float, b: float, c: float):
        if not (a >= 0 and b >= 0 and c >= 0):
            raise ValueError("a,b,c must each be >= 0")
        s = a + b + c
        if not (s > 0):
            raise ValueError("a+b+c must be > 0")
        a, b, c = a / s, b / s, c / s
        self.a, self.b, self.c = a, b, c
        self.w_perp = 0.75 * a
        self.w_par = (2.0 / 3.0) * b + (1.0 / 3.0) * c
        self.eta = a + (2.0 / 3.0) * b + (1.0 / 3.0) * c
        W_perp = self.w_perp * INT_SIN2   # = 2π a
        W_par = self.w_par * INT_PAR      # = 2π w_par
        self.P_perp = W_perp / (W_perp + W_par)
        assert abs(self.P_perp + W_par / (W_perp + W_par) - 1.0) < 1e-12
        assert -1e-15 <= self.P_perp <= 1.0 + 1e-15


def make_mode_weights(a, b, c) -> ModeWeights:
    return ModeWeights(a, b, c)


# ---------------------------------------------------------------------------
# Per-shape inverse-CDF and rejection samplers -- mirror of the header.
# ---------------------------------------------------------------------------
def sample_costheta_perp_invcdf(u: float) -> float:
    # x³ - 3x + (4u-2) = 0 ; in-range trig root (k=1)
    x = 2.0 * math.cos(math.acos(1.0 - 2.0 * u) / 3.0 - 2.0 * PI / 3.0)
    return min(1.0, max(-1.0, x))


def sample_costheta_par_invcdf(u: float) -> float:
    # x³ + x + (2-4u) = 0 ; single real root (Cardano)
    q = 2.0 - 4.0 * u
    disc = q * q / 4.0 + 1.0 / 27.0
    sq = math.sqrt(disc)
    # math.cbrt (py3.11+) -> same libm cbrt as C++ std::cbrt, for bit-for-bit parity.
    x = math.cbrt(-q / 2.0 + sq) + math.cbrt(-q / 2.0 - sq)
    return min(1.0, max(-1.0, x))


def sample_costheta_perp_reject(rng: Rng) -> float:
    while True:
        x = 2.0 * rng() - 1.0
        if rng() <= (1.0 - x * x):
            return x


def sample_costheta_par_reject(rng: Rng) -> float:
    while True:
        x = 2.0 * rng() - 1.0
        if rng() <= (0.25 + 0.75 * x * x):
            return x


def sample_local_direction(m: ModeWeights, rng: Rng, use_rejection: bool = False):
    if rng() < m.P_perp:
        x = (sample_costheta_perp_reject(rng) if use_rejection
             else sample_costheta_perp_invcdf(rng()))
    else:
        x = (sample_costheta_par_reject(rng) if use_rejection
             else sample_costheta_par_invcdf(rng()))
    phi = 2.0 * PI * rng()
    st = math.sqrt(max(0.0, 1.0 - x * x))
    return (st * math.cos(phi), st * math.sin(phi), x)


# ---------------------------------------------------------------------------
# Gram-Schmidt frame + rotation -- mirror of the header.
# ---------------------------------------------------------------------------
def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _normalized(a):
    n = math.sqrt(_dot(a, a))
    if not (n > 0):
        raise ValueError("cannot normalize zero/NaN vector")
    return (a[0] / n, a[1] / n, a[2] / n)


def build_frame(Bhat):
    b = _normalized(Bhat)
    ax, ay, az = abs(b[0]), abs(b[1]), abs(b[2])
    if ax <= ay and ax <= az:
        ref = (1.0, 0.0, 0.0)
    elif ay <= az:
        ref = (0.0, 1.0, 0.0)
    else:
        ref = (0.0, 0.0, 1.0)
    d = _dot(ref, b)
    ex = _normalized((ref[0] - d * b[0], ref[1] - d * b[1], ref[2] - d * b[2]))
    ey = _cross(b, ex)
    return ex, ey, b


def to_global(local, frame):
    ex, ey, ez = frame
    return tuple(local[0] * ex[i] + local[1] * ey[i] + local[2] * ez[i] for i in range(3))


def sample_global_direction(m: ModeWeights, Bhat, rng: Rng, use_rejection: bool = False):
    local = sample_local_direction(m, rng, use_rejection)
    u = to_global(local, build_frame(Bhat))
    assert abs(math.sqrt(_dot(u, u)) - 1.0) < 1e-12
    return u


def pdf_costheta(m: ModeWeights, x: float) -> float:
    """Normalized analytic marginal pdf in x = cosθ (oracle for chi-square)."""
    I = m.w_perp * (1.0 - x * x) + m.w_par * (0.25 + 0.75 * x * x)
    Z = m.w_perp * (4.0 / 3.0) + m.w_par * 1.0
    return I / Z


# ---------------------------------------------------------------------------
# Magnetic-field direction models -- 1:1 mirror of src/spf_field.hpp.
# Parity is checked in tests/test_field_parity.py (C++ vs Python bhat).
# B-hat conventions are documented in spf_field.hpp; the AngledField (alpha,beta)
# convention is pinned to anarrima's angled kernels (see analytic_nwl.angled_*).
# ---------------------------------------------------------------------------
class ConstantField:
    __slots__ = ("b",)

    def __init__(self, b0):
        self.b = _normalized(b0)

    def bhat(self, r):
        return self.b


class ToroidalField:
    __slots__ = ()

    def bhat(self, r):
        x, y, _z = r
        rxy = math.sqrt(x * x + y * y)
        return (-y / rxy, x / rxy, 0.0)  # phi_hat


class AngledField:
    """B-hat = cosβ·phi_hat + sinβ·(cosα·R_hat − sinα·z_hat). Unit; β=0 ⇒ toroidal."""
    __slots__ = ("ca", "sa", "cb", "sb")

    def __init__(self, alpha, beta):
        self.ca, self.sa = math.cos(alpha), math.sin(alpha)
        self.cb, self.sb = math.cos(beta), math.sin(beta)

    def bhat(self, r):
        x, y, _z = r
        rxy = math.sqrt(x * x + y * y)
        cphi, sphi = x / rxy, y / rxy
        return (self.cb * (-sphi) + self.sb * self.ca * cphi,
                self.cb * (cphi) + self.sb * self.ca * sphi,
                -self.sb * self.sa)


def make_field(bmode: str, alpha: float = 0.0, beta: float = 0.0,
               b=(0.0, 0.0, 1.0), path: str = None):
    """Mirror of the C++ plugin's field-selection (bmode parsing)."""
    if bmode == "toroidal":
        return ToroidalField()
    if bmode == "constant":
        return ConstantField(b)
    if bmode == "angled":
        return AngledField(alpha, beta)
    if bmode == "fieldmap":
        from fieldmap import FieldMapField  # lazy: needs numpy
        return FieldMapField(path)
    raise ValueError(f"unknown bmode {bmode!r}")
