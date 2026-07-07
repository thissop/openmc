"""Unit tests for the native SPF addition to ``openmc.TokamakSource``.

Two independent layers:

1. **Sampler-math tests** (always run). A self-contained NumPy mirror of the
   exact math in ``native_spf/source.cpp``
   (``field_direction`` + ``sample_polarized_direction``) reproduces the
   Schwartz P2 birth-direction distribution. These validate the physics
   (isotropic limit, perpendicular vs parallel concentration, (a,b,c)
   normalization, Gram-Schmidt rotation correctness) without building OpenMC.
   The formulas here are line-for-line the ones in source.cpp; only the RNG is
   swapped for NumPy (statistical, not bit-for-bit -- the bit-for-bit port-vs-
   ``spf_sampler.hpp`` parity test is a separate build-time gate, PLAN A §4.1).

2. **Python-API tests** (skipped if ``openmc`` is not importable). Exercise the
   new ``polarization`` / ``field_model`` / ``safety_factor`` / ``field_sign``
   kwargs, the (a,b,c)-vs-spin-fraction setter, and XML round-trip. Once the
   PR + this patch are installed on the cluster these activate automatically.

Tests are written as module-level functions (OpenMC's pytest.ini sets
``python_classes = NoThanks``, i.e. class-based tests are not collected).

Run:  pytest -q native_spf/test_spf_tokamak.py
"""
from __future__ import annotations

import math

import numpy as np
import pytest


# ===========================================================================
# NumPy mirror of native_spf/source.cpp SPF math (physics oracle).
# ===========================================================================

def mode_weights(a, b, c):
    """Mirror of the constructor block in source.cpp: derive (w_perp, w_par,
    eta, p_perp) from (a,b,c). Renormalizes to sum 1. No magic constants:
    coefficients are Schwartz Eq. 2's (3/4), (2/3), (1/3)."""
    if not (a >= 0 and b >= 0 and c >= 0):
        raise ValueError("a,b,c must each be >= 0")
    s = a + b + c
    if not (s > 0):
        raise ValueError("a+b+c must be > 0")
    a, b, c = a / s, b / s, c / s
    w_perp = 0.75 * a
    w_par = (2.0 / 3.0) * b + (1.0 / 3.0) * c
    eta = a + (2.0 / 3.0) * b + (1.0 / 3.0) * c
    p_perp = a / eta
    assert -1e-15 <= p_perp <= 1.0 + 1e-15
    return w_perp, w_par, eta, p_perp


def _costheta_perp(u):
    """sin²θ shape: root in [-1,1] of x³ - 3x + (4u-2) = 0 (source.cpp)."""
    x = 2.0 * np.cos(np.arccos(1.0 - 2.0 * u) / 3.0 - 2.0 * np.pi / 3.0)
    return np.clip(x, -1.0, 1.0)


def _costheta_par(u):
    """1/4+3/4cos²θ shape: single real root of x³ + x + (2-4u) = 0 (Cardano)."""
    q = 2.0 - 4.0 * u
    sq = np.sqrt(q * q / 4.0 + 1.0 / 27.0)
    x = np.cbrt(-q / 2.0 + sq) + np.cbrt(-q / 2.0 - sq)
    return np.clip(x, -1.0, 1.0)


def _build_frame(bhat):
    """Gram-Schmidt with least-aligned reference axis (source.cpp)."""
    b = np.asarray(bhat, dtype=float)
    b = b / np.linalg.norm(b)
    ax, ay, az = abs(b[0]), abs(b[1]), abs(b[2])
    if ax <= ay and ax <= az:
        ref = np.array([1.0, 0.0, 0.0])
    elif ay <= az:
        ref = np.array([0.0, 1.0, 0.0])
    else:
        ref = np.array([0.0, 0.0, 1.0])
    ex = ref - np.dot(ref, b) * b
    ex = ex / np.linalg.norm(ex)
    ey = np.cross(b, ex)
    return ex, ey, b


def sample_polarized_directions(a, b, c, bhat, n, rng):
    """Vectorized mirror of TokamakSource::sample_polarized_direction for a
    single fixed b̂. Returns (u [n,3], x [n]) where x = local cosθ = u·b̂."""
    _, _, _, p_perp = mode_weights(a, b, c)
    sel = rng.random(n)
    u_shape = rng.random(n)
    x = np.where(sel < p_perp, _costheta_perp(u_shape), _costheta_par(u_shape))
    az = 2.0 * np.pi * rng.random(n)
    st = np.sqrt(np.maximum(0.0, 1.0 - x * x))
    local = np.column_stack([st * np.cos(az), st * np.sin(az), x])
    ex, ey, bz = _build_frame(bhat)
    u = (local[:, 0:1] * ex + local[:, 1:2] * ey + local[:, 2:3] * bz)
    return u, x


def field_direction_toroidal(phi):
    """Mirror of field_model==0: b̂ = φ̂ = (-sin φ, cos φ, 0)."""
    return np.array([-math.sin(phi), math.cos(phi), 0.0])


def field_direction_pitched(r, alpha, phi, *, R0, a_minor, kappa, delta, Delta,
                            q, sign):
    """Mirror of field_model==1 (single q value = flat profile) in source.cpp."""
    sphi, cphi = math.sin(phi), math.cos(phi)
    psi = alpha + delta * math.sin(alpha)
    R = R0 + r * math.cos(psi) + Delta * (1.0 - (r * r) / (a_minor * a_minor))
    tR = -r * math.sin(psi) * (1.0 + delta * math.cos(alpha))
    tZ = kappa * r * math.cos(alpha)
    tnorm = math.hypot(tR, tZ)
    if not (tnorm > 0.0):
        return np.array([-sphi, cphi, 0.0])
    pR, pZ = tR / tnorm, tZ / tnorm
    lam = (r / (q * R)) if q > 0.0 else 0.0
    sl = sign * lam
    Bx = -sphi + sl * pR * cphi
    By = cphi + sl * pR * sphi
    Bz = sl * pZ
    return np.array([Bx, By, Bz]) / math.sqrt(Bx * Bx + By * By + Bz * Bz)


# Analytic <cos²θ> of each shape (source of the pass/fail targets):
#   sin²θ shape p(x)∝(1-x²):        <x²> = 1/5
#   1/4+3/4cos²θ shape p=1/4+3/4x²: <x²> = 7/15
COS2_PERP = 1.0 / 5.0
COS2_PAR = 7.0 / 15.0
COS2_ISO = 1.0 / 3.0

N = 400_000


# ===========================================================================
# Layer 1: mode-weight identities.
# ===========================================================================

def test_mode_weights_corner_p_perp():
    # A -> perp shape only; B,C -> par shape only; nonpol -> 1/2.
    assert mode_weights(1, 0, 0)[3] == pytest.approx(1.0)
    assert mode_weights(0, 1, 0)[3] == pytest.approx(0.0)
    assert mode_weights(0, 0, 1)[3] == pytest.approx(0.0)
    assert mode_weights(1, 1, 1)[3] == pytest.approx(0.5)


def test_mode_weights_eta_corners():
    # Total-rate factor eta = a + 2b/3 + c/3.
    assert mode_weights(1, 1, 1)[2] == pytest.approx(2.0 / 3.0)  # nonpol
    assert mode_weights(1, 0, 0)[2] == pytest.approx(1.0)        # A
    assert mode_weights(0, 1, 0)[2] == pytest.approx(2.0 / 3.0)  # B
    assert mode_weights(0, 0, 1)[2] == pytest.approx(1.0 / 3.0)  # C


def test_abc_normalization_scale_invariant():
    # (a,b,c) normalization: only the direction (p_perp) matters, and it is
    # invariant to an overall scale of the input triple.
    base = mode_weights(0.5, 0.3, 0.2)
    scaled = mode_weights(5.0, 3.0, 2.0)
    assert scaled[3] == pytest.approx(base[3])
    assert scaled[0] == pytest.approx(base[0])
    assert scaled[1] == pytest.approx(base[1])


def test_b_and_c_share_shape():
    # B and C differ only in eta (rate), never in angular shape (w_perp == 0 for
    # both) -- the §1.2 correction.
    assert mode_weights(0, 1, 0)[0] == pytest.approx(0.0)  # w_perp
    assert mode_weights(0, 0, 1)[0] == pytest.approx(0.0)  # w_perp


def test_mode_weights_reject_bad_input():
    with pytest.raises(ValueError):
        mode_weights(-0.1, 0.5, 0.6)
    with pytest.raises(ValueError):
        mode_weights(0, 0, 0)


# ===========================================================================
# Layer 1: sampled-direction distribution.
# ===========================================================================

def test_nonpol_is_isotropic():
    # Mode (1/3,1/3,1/3) has a FLAT cosθ pdf -> isotropic. Both the mixture
    # <cos²θ> == 1/3 and a per-bin flatness check.
    rng = np.random.default_rng(12345)
    _, x = sample_polarized_directions(1, 1, 1, [0, 0, 1], N, rng)
    assert np.mean(x) == pytest.approx(0.0, abs=3e-3)
    assert np.mean(x * x) == pytest.approx(COS2_ISO, abs=3e-3)
    counts, _ = np.histogram(x, bins=10, range=(-1, 1))
    expected = N / 10
    # Flat within ~5 sigma of Poisson per bin.
    assert np.max(np.abs(counts - expected)) < 5.0 * math.sqrt(expected)


def test_pure_perp_concentrates_perpendicular():
    # Mode A about a deliberately TILTED b̂: emission concentrates ⊥ b̂, so
    # <(u·b̂)²> -> 1/5 < 1/3, and the density peaks at cosθ = 0.
    rng = np.random.default_rng(7)
    bhat = np.array([1.0, 1.0, 1.0])
    u, _ = sample_polarized_directions(1, 0, 0, bhat, N, rng)
    cos_to_b = u @ (bhat / np.linalg.norm(bhat))
    assert np.mean(cos_to_b ** 2) == pytest.approx(COS2_PERP, abs=3e-3)
    assert np.mean(cos_to_b ** 2) < COS2_ISO
    assert np.allclose(np.linalg.norm(u, axis=1), 1.0, atol=1e-12)
    # Density near the equator (|cos|<0.2) exceeds density near the poles.
    near_eq = np.mean(np.abs(cos_to_b) < 0.2)
    near_pole = np.mean(np.abs(cos_to_b) > 0.8)
    assert near_eq > near_pole


def test_pure_par_concentrates_parallel():
    # Mode B about b̂: concentrates ∥ b̂, <(u·b̂)²> -> 7/15 > 1/3.
    rng = np.random.default_rng(99)
    bhat = np.array([0.3, -0.4, 0.5])
    u, _ = sample_polarized_directions(0, 1, 0, bhat, N, rng)
    cos_to_b = u @ (bhat / np.linalg.norm(bhat))
    assert np.mean(cos_to_b ** 2) == pytest.approx(COS2_PAR, abs=3e-3)
    assert np.mean(cos_to_b ** 2) > COS2_ISO


def test_pure_c_matches_pure_b_shape():
    # Isotropic-vs-C / B-vs-C: pure C must have the SAME angular shape as B
    # (both par-only). <cos²θ> equal within statistics.
    rng = np.random.default_rng(2024)
    _, xb = sample_polarized_directions(0, 1, 0, [0, 0, 1], N, rng)
    _, xc = sample_polarized_directions(0, 0, 1, [0, 0, 1], N, rng)
    assert np.mean(xb ** 2) == pytest.approx(np.mean(xc ** 2), abs=4e-3)


# ===========================================================================
# Layer 1: local -> global rotation correctness.
# ===========================================================================

def test_rotation_preserves_angle_to_bhat():
    # Gram-Schmidt correctness: for an arbitrary tilted b̂, the global u·b̂ must
    # equal the sampled local cosθ to ~machine precision.
    rng = np.random.default_rng(55)
    bhat = np.array([0.2, 0.9, -0.3])
    u, x = sample_polarized_directions(0.5, 0.3, 0.2, bhat, 50_000, rng)
    cos_to_b = u @ (bhat / np.linalg.norm(bhat))
    assert np.max(np.abs(cos_to_b - x)) < 1e-12


def test_isotropic_stays_isotropic_off_axis():
    # Nonpol about a tilted b̂ is isotropic in the GLOBAL frame: each component
    # has zero mean and <u_i²> = 1/3 (no rotation bias).
    rng = np.random.default_rng(321)
    u, _ = sample_polarized_directions(1, 1, 1, [1, 2, 3], N, rng)
    assert np.allclose(np.mean(u, axis=0), 0.0, atol=4e-3)
    assert np.allclose(np.mean(u ** 2, axis=0), 1.0 / 3.0, atol=4e-3)


def test_degenerate_bhat_along_z():
    # b̂ = ±ẑ must not blow up the Gram-Schmidt branch.
    rng = np.random.default_rng(1)
    for bz in (+1.0, -1.0):
        u, x = sample_polarized_directions(1, 0, 0, [0, 0, bz], 20_000, rng)
        # [0,0,bz] is already unit, so u·b̂ == local cosθ == x for either sign.
        cos_to_b = u @ np.array([0.0, 0.0, bz])
        assert np.all(np.isfinite(u))
        assert np.max(np.abs(cos_to_b - x)) < 1e-12


# ===========================================================================
# Layer 1: field-direction models.
# ===========================================================================

def test_toroidal_is_phi_hat():
    for phi in (0.0, 0.7, 2.5, 5.9):
        b = field_direction_toroidal(phi)
        assert b == pytest.approx([-math.sin(phi), math.cos(phi), 0.0], abs=1e-14)
        assert np.linalg.norm(b) == pytest.approx(1.0)


def test_pitched_reduces_to_toroidal_as_q_grows():
    # As q -> ∞ (lambda -> 0) the pitched field converges to φ̂.
    kw = dict(R0=250.0, a_minor=100.0, kappa=1.0, delta=0.0, Delta=0.0, sign=1)
    phi, alpha, r = 1.1, 0.8, 60.0
    b_tor = field_direction_toroidal(phi)
    b_big_q = field_direction_pitched(r, alpha, phi, q=1e6, **kw)
    assert np.allclose(b_big_q, b_tor, atol=1e-4)


def test_pitched_is_unit_and_tilted():
    kw = dict(R0=250.0, a_minor=100.0, kappa=1.8, delta=0.4, Delta=0.0, sign=1)
    b = field_direction_pitched(60.0, 0.8, 1.1, q=2.0, **kw)
    assert np.linalg.norm(b) == pytest.approx(1.0)
    # A finite q introduces a nonzero poloidal (z) component (up-down tilt).
    assert abs(b[2]) > 1e-3


# ===========================================================================
# Layer 2: Python-API tests (skipped if openmc is not importable).
# ===========================================================================

try:
    import openmc
    _HAVE_OPENMC = hasattr(openmc, "TokamakSource")
except Exception:  # pragma: no cover - import env dependent
    openmc = None
    _HAVE_OPENMC = False

requires_openmc = pytest.mark.skipif(
    not _HAVE_OPENMC,
    reason="openmc (with the SPF-patched TokamakSource) not importable here")


def _minimal_tokamak(**extra):
    """A minimal valid TokamakSource for API tests."""
    r = np.linspace(0.0, 1.0, 5)
    return openmc.TokamakSource(
        major_radius=250.0, minor_radius=100.0, elongation=1.0,
        triangularity=0.0, shafranov_shift=0.0, r_over_a=r,
        emission_density=1.0 - r ** 2,
        energy=openmc.stats.Discrete([14.06e6], [1.0]), **extra)


@requires_openmc
def test_polarization_none_backward_compatible():
    src = _minimal_tokamak()
    assert src.polarization is None
    elem = src.to_xml_element()
    # No SPF elements written -> byte-identical XML to the pre-SPF class.
    assert elem.find("polarization") is None
    assert elem.find("field_model") is None


@requires_openmc
def test_abc_setter_normalizes_and_warns():
    src = _minimal_tokamak()
    with pytest.warns(UserWarning):
        src.polarization = (2.0, 0.0, 0.0)
    assert src.polarization == pytest.approx((1.0, 0.0, 0.0))


@requires_openmc
def test_abc_setter_rejects_negative():
    src = _minimal_tokamak()
    with pytest.raises(Exception):
        src.polarization = (-0.1, 0.5, 0.6)


@requires_openmc
def test_spin_fraction_dict_maps_eq1():
    # Fully aligned d+,t+ -> pure A: (a,b,c) = (1,0,0).
    src = _minimal_tokamak(polarization={
        'd_plus': 1.0, 'd_zero': 0.0, 'd_minus': 0.0,
        't_plus': 1.0, 't_minus': 0.0})
    assert src.polarization == pytest.approx((1.0, 0.0, 0.0))


@requires_openmc
def test_spin_fractions_to_abc_helper():
    a, b, c = openmc.source.spin_fractions_to_abc(
        d_plus=0.5, d_zero=0.2, d_minus=0.3, t_plus=0.6, t_minus=0.4)
    assert (a, b, c) == pytest.approx(
        (0.5 * 0.6 + 0.3 * 0.4, 0.2, 0.5 * 0.4 + 0.3 * 0.6))


@requires_openmc
def test_pitched_requires_safety_factor():
    with pytest.raises(ValueError):
        _minimal_tokamak(polarization=(1.0, 0.0, 0.0), field_model='pitched')


@requires_openmc
def test_xml_roundtrip_toroidal():
    src = _minimal_tokamak(polarization=(0.5, 0.3, 0.2))
    elem = src.to_xml_element()
    assert elem.find("polarization") is not None  # element present
    back = openmc.TokamakSource.from_xml_element(elem)
    assert back.polarization == pytest.approx((0.5, 0.3, 0.2))
    assert back.field_model == 'toroidal'


@requires_openmc
def test_xml_roundtrip_pitched():
    rq = np.linspace(0.0, 1.0, 4)
    q = 1.0 + 2.0 * rq ** 2
    src = _minimal_tokamak(polarization=(1.0, 0.0, 0.0), field_model='pitched',
                           safety_factor=(rq, q), field_sign=-1)
    back = openmc.TokamakSource.from_xml_element(src.to_xml_element())
    assert back.field_model == 'pitched'
    assert back.field_sign == -1
    assert np.allclose(back.safety_factor[0], rq)
    assert np.allclose(back.safety_factor[1], q)
