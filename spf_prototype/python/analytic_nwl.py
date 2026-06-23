"""Independent analytic neutron wall loading (NWL) for a filamentary ring source
inside a rectangular-cross-section torus (Schwartz 2025, arXiv:2507.11758, §3).

Two independent engines for the single-ring reduced intensity g:
  * `g_ring_quad`  -- vectorized Gauss-Legendre quadrature of Eq. 3-4 (ours), and
  * `g_ring_anarrima` -- Schwartz's published closed forms (anarrima, the reference).
Plus `g_RiI_closed_eq5` -- our own transcription of the inboard-isotropic closed
form Eq. 5 (a third, independent path for that one case).

Reduced intensity (per unit ring length, dropping E_n/4π) onto a wall target:
    g = ∫_{-φm}^{φm} p·cosν·f(cosθ) / |Δ|² dφ
    |Δ|² = p² + r² + z² − 2 p r cosφ ,  cosθ = r sinφ / |Δ|   (purely toroidal B̂)
geometry-specific cosν and visible half-angle φm:
    inboard  (n̂=+R̂):  cosν=(p cosφ − r)/|Δ|,  φm = arccos(r/p)            [r<p]
    outboard (n̂=−R̂):  cosν=(r − p cosφ)/|Δ|,  φm = arccos(u/p)+arccos(u/r) [r>p]
    floor    (n̂=+ẑ):   cosν= z/|Δ|,            φm = arccos(u/p)+arccos(u/r)
where u is the inboard (center-stack) radius that shadows floor/outboard views.

For a volumetric source Φ(p,z0): Q(target) = ∫∫ Φ(p,z0) g(target;p,z0) dp dz0
(the p of dV=p dp dz0 dφ cancels the 1/p in g/p of Eq. 3 — no extra p weight).
"""
from __future__ import annotations

import numpy as np
from scipy import integrate
from scipy import special

# --------------------------------------------------------------------------
# Angular shape factors f(cosθ)
# --------------------------------------------------------------------------
FACTORS = {
    "iso": lambda ct: np.ones_like(ct),
    "A": lambda ct: 1.0 - ct * ct,          # sin²θ
    "cos2": lambda ct: ct * ct,             # cos²θ
    "BC": lambda ct: 0.25 + 0.75 * ct * ct,  # ¼ + ¾cos²θ
}


def _phi_limit(p, r, wall, u_in):
    """Visible half-angle φm (array-safe). Args clipped to acos domain."""
    if wall == "inboard":
        return np.arccos(np.clip(r / p, -1.0, 1.0))
    # floor / outboard: limited by the center-stack tangent from source and target
    return (np.arccos(np.clip(u_in / p, -1.0, 1.0))
            + np.arccos(np.clip(u_in / r, -1.0, 1.0)))


def _cosnu(p, r, z, cphi, d, wall):
    if wall == "inboard":
        return (p * cphi - r) / d
    if wall == "outboard":
        return (r - p * cphi) / d
    if wall == "floor":
        return z / d
    raise ValueError(wall)


# --------------------------------------------------------------------------
# Engine 1: our vectorized Gauss-Legendre quadrature (p, z may be arrays)
# --------------------------------------------------------------------------
_GL_NODES, _GL_WTS = np.polynomial.legendre.leggauss(256)


def g_ring_quad(p, z, r, wall, factor, u_in):
    """Reduced intensity from ring(s) at (p, z) onto a target of radius r.
    p, z broadcastable arrays; r, u_in scalars. Returns array like p."""
    p = np.asarray(p, dtype=float)
    z = np.asarray(z, dtype=float)
    f = FACTORS[factor]
    phim = _phi_limit(p, r, wall, u_in)                  # shape p
    # nodes: phi = phim * x  (x in [-1,1]); integrate -> phim * sum(w f)
    phi = phim[..., None] * _GL_NODES                    # shape p + (G,)
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    d2 = p[..., None] ** 2 + r ** 2 + z[..., None] ** 2 - 2 * p[..., None] * r * cphi
    d = np.sqrt(d2)
    cosnu = _cosnu(p[..., None], r, z[..., None], cphi, d, wall)
    ct = r * sphi / d
    integ = p[..., None] * cosnu * f(ct) / d2
    return phim * np.sum(integ * _GL_WTS, axis=-1)


def g_ring_quad_scalar(p, z, r, wall, factor, u_in):
    """High-accuracy scalar quadrature (scipy.quad) for the validation table."""
    f = FACTORS[factor]

    def integrand(phi):
        d2 = p * p + r * r + z * z - 2 * p * r * np.cos(phi)
        d = np.sqrt(d2)
        cosnu = _cosnu(p, r, z, np.cos(phi), d, wall)
        ct = r * np.sin(phi) / d
        return p * cosnu * f(ct) / d2

    lim = float(_phi_limit(np.asarray(p), r, wall, u_in))
    val, _ = integrate.quad(integrand, -lim, lim, limit=400)
    return val


# --------------------------------------------------------------------------
# Engine 2: anarrima (Schwartz's published closed forms). Lazy import.
# --------------------------------------------------------------------------
def _anarrima():
    import os
    os.environ.setdefault("JAX_ENABLE_X64", "1")
    from anarrima import reduced_integrals as ri
    return ri


def _an_fn(comp, factor):
    ri = _anarrima()
    return {
        ("H", "iso"): ri.g_HI, ("H", "A"): ri.g_HA, ("H", "cos2"): ri.g_Hc, ("H", "BC"): ri.g_HB,
        ("V", "iso"): ri.g_VI, ("V", "A"): ri.g_VA, ("V", "cos2"): ri.g_Vc, ("V", "BC"): ri.g_VB,
    }[(comp, factor)]


def g_ring_anarrima(p, z, r, wall, factor, u_in):
    comp = "V" if wall == "floor" else "H"
    fn = _an_fn(comp, factor)
    phim = float(_phi_limit(np.asarray(p), r, wall, u_in))
    val = float(fn(p, z, r, phim))
    return -val if wall == "outboard" else val  # outboard normal is -R̂


_AN_VEC = {}


def g_ring_anarrima_vec(p, z, r, wall, factor, u_in):
    """Vectorized anarrima engine (jit+jnp.vectorize) -- engine-compatible with
    g_ring_quad so plasma_patterns can run on the published closed forms."""
    import jax
    import jax.numpy as jnp
    comp = "V" if wall == "floor" else "H"
    key = (comp, factor)
    if key not in _AN_VEC:
        _AN_VEC[key] = jax.jit(jnp.vectorize(_an_fn(comp, factor)))
    p = np.asarray(p, dtype=float)
    z = np.asarray(z, dtype=float)
    phim = _phi_limit(p, r, wall, u_in)
    val = np.asarray(_AN_VEC[key](p, z, r, phim), dtype=float)
    return -val if wall == "outboard" else val


# --------------------------------------------------------------------------
# Engine 3: our transcription of inboard-isotropic closed form (Eq. 5)
# --------------------------------------------------------------------------
def g_RiI_closed_eq5(p, z, u):
    """Eq. 5: inboard isotropic, independent of anarrima. Uses scipy elliptic
    integrals with parameter m (F = ellipkinc, E = ellipeinc)."""
    phi_i = np.arccos(u / p)
    amp = phi_i / 2.0
    m = -4.0 * p * u / ((p - u) ** 2 + z ** 2)
    q = p ** 4 - 2 * p ** 2 * (u ** 2 - z ** 2) + (u ** 2 + z ** 2) ** 2
    s = np.sqrt((p - u) ** 2 + z ** 2)
    F = special.ellipkinc(amp, m)
    E = special.ellipeinc(amp, m)
    term_L = 4 * p * np.sqrt((p ** 2 - u ** 2) * (p ** 2 - u ** 2 + z ** 2)) / q
    term_F = -2 * p / (u * s) * F
    term_E = 2 * p * (p ** 2 - u ** 2 + z ** 2) * s / (q * u) * E
    return term_L + term_F + term_E


# ==========================================================================
# Geometry: exactly matches anarrima's examples/square_torus.py (the script that
# generates Schwartz Fig. 2). Plasma R0=1, a=0.5 (aspect 2.0 -- note the paper
# TEXT says 2.5, but the figure-generating example uses 2.0); walls stand off
# from the plasma by ~0.1: inboard u=0.4, outboard w=1.6, floor/ceiling z=∓0.6.
# ==========================================================================
R0 = 1.0
A_MINOR = 0.5
ASPECT = R0 / A_MINOR        # 2.0 (per the anarrima example; paper text says 2.5)
R_IN = 0.4                   # inboard wall radius u (also the center-stack shadow)
R_OUT = 1.6                  # outboard wall radius w
Z_WALL = 0.6                 # floor/ceiling at z = ∓Z_WALL
P_MIN = R0 - A_MINOR         # 0.5  plasma bounding box (poloidal)
P_MAX = R0 + A_MINOR         # 1.5


def wall_targets(n_per_wall=60):
    """Return list of dicts for points on inboard/outboard/floor/ceiling walls.
    's' is a signed poloidal coordinate (midplane=0 on vertical walls; center=0
    on horizontal walls) used only for plotting."""
    out = []
    zt = np.linspace(-Z_WALL, Z_WALL, n_per_wall)
    rt = np.linspace(R_IN, R_OUT, n_per_wall)
    for z in zt:
        out.append(dict(wall="inboard", R=R_IN, Z=z, s=z))
    for z in zt:
        out.append(dict(wall="outboard", R=R_OUT, Z=z, s=z))
    for r in rt:
        out.append(dict(wall="floor", R=r, Z=-Z_WALL, s=r - R0))
    for r in rt:
        out.append(dict(wall="ceiling", R=r, Z=Z_WALL, s=r - R0))
    return out


def _ring_grid(n=120):
    """Parabolic (1-ρ²) plasma: poloidal-disk grid of rings (p, z0) with weights.
    Grid spans the plasma bounding box [R0-a, R0+a] x [-a, a], masked to the LCFS."""
    ps = np.linspace(P_MIN, P_MAX, n)
    zs = np.linspace(-A_MINOR, A_MINOR, n)
    P, Z = np.meshgrid(ps, zs, indexing="ij")
    rho2 = ((P - R0) ** 2 + Z ** 2) / A_MINOR ** 2
    inside = rho2 <= 1.0
    w = np.where(inside, 1.0 - rho2, 0.0)  # parabolic profile, 0 outside
    return P[inside], Z[inside], w[inside]


def _target_primitives(target, P, Z0, engine=g_ring_quad):
    """g_iso, g_A, g_cos2 summed over rings for one wall target (with relative z)."""
    wall = target["wall"]
    r = target["R"]
    if wall == "ceiling":
        # up-down symmetry: ceiling(Z=+a) == floor(Z=-a) for symmetric plasma
        wall_eff, z_rel = "floor", Z_WALL - Z0   # = (Z0 - (-(Z_WALL)))->mirror
    elif wall == "floor":
        wall_eff, z_rel = "floor", Z0 - (-Z_WALL)
    else:  # inboard / outboard vertical walls
        wall_eff, z_rel = wall, Z0 - target["Z"]
    gi = engine(P, z_rel, r, wall_eff, "iso", R_IN)
    ga = engine(P, z_rel, r, wall_eff, "A", R_IN)
    gc = engine(P, z_rel, r, wall_eff, "cos2", R_IN)
    return gi, ga, gc


# Mode definitions (a,b,c) and the constant-rate factor η = a + 2b/3 + c/3.
MODES = {
    "iso": (1 / 3, 1 / 3, 1 / 3),
    "A": (1.0, 0.0, 0.0),
    "B": (0.0, 1.0, 0.0),
    "C": (0.0, 0.0, 1.0),
    "mixed": (0.5, 0.3, 0.2),
}


def eta(a, b, c):
    return a + 2 * b / 3 + c / 3


def bracket(a, b, c, G_iso, G_A, G_cos2):
    """Physical dσ/dΩ-shape NWL (Eq. 2 bracket, σ0/2π dropped), linear in (a,b,c):
       ¾a·g_A + (⅔b+⅓c)·(¼ g_iso + ¾ g_cos2)."""
    G_BC = 0.25 * G_iso + 0.75 * G_cos2
    return 0.75 * a * G_A + (2 * b / 3 + c / 3) * G_BC


def plasma_patterns(targets, n_grid=120, engine=g_ring_quad):
    """Return dict wall-load arrays {iso,A,cos2} (∫∫(1-ρ²) g dp dz0) over targets,
    plus per-mode bracket and constant-rate-directional patterns."""
    P, Z0, w = _ring_grid(n_grid)
    Gi = np.empty(len(targets)); Ga = np.empty(len(targets)); Gc = np.empty(len(targets))
    for k, t in enumerate(targets):
        gi, ga, gc = _target_primitives(t, P, Z0, engine)
        Gi[k] = np.sum(w * gi); Ga[k] = np.sum(w * ga); Gc[k] = np.sum(w * gc)
    res = {"iso_g": Gi, "A_g": Ga, "cos2_g": Gc}
    for name, (a, b, c) in MODES.items():
        br = bracket(a, b, c, Gi, Ga, Gc)
        res[f"{name}_bracket"] = br
        res[f"{name}_directional"] = br / eta(a, b, c)  # constant total rate
    return res
