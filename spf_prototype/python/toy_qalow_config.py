"""Shared definition of the toy QA-low free-streaming cross-check problem, used by
BOTH the analytic step (anarrima venv) and the OpenMC step (spf venv) so the two
sides solve the IDENTICAL problem. Pure numpy: no anarrima, no OpenMC.

The plasma is the anarrima tier3 toy QA-low spectrum (Nfp=2, eps_eff~0.07 -- the
regime where the analytic series, exact quadrature, and MC should all agree). The
field is the same toy model field anarrima uses. The wall is a square-cross-section
torus (axisymmetric, shared, CSG -> dev-box, no DAGMC); all the non-axisymmetry is
in the SOURCE. We compare DIRECTIONALITY (A/iso, B/iso per wall patch), which
cancels source-position/normalization details and isolates the polarized emission
about the local field -- the physics being validated.

Conventions (match anarrima/VMEC): R(theta,phi)=sum R_mn cos(m th - n Nfp phi),
Z=sum Z_mn sin(...); phi geometric cylindrical; field B0(alpha(phi),beta(phi),phi)
per anarrima Eq.12. cos theta = Bhat . Delta_hat (local field . line of sight).
"""
from __future__ import annotations

import numpy as np

STEM = "toy_qalow"               # /tmp/<STEM>_{analytic,openmc}.npz scratch prefix
TITLE = "Toy QA-low"
NFP = 2
SHAPE = 0.6                       # QA-low non-axisymmetric scale (tier3 qa_spectrum(0.6))
IOTA = 0.42
FIELD_RIPPLE = 0.10
RHO_SURFACES = [0.6, 1.0]
W_RHO = [0.6 ** 2, 1.0 ** 2]      # rho^2 emissivity weight
N_THETA_S = 14                    # source loops per surface (poloidal)

# boundary spectrum rows: (m, n, R_mn, Z_mn). n!=0 rows scaled by SHAPE.
_BASE = [(0, 0, 1.00, 0.00), (1, 0, 0.166, 0.166), (2, 0, 0.012, -0.011)]
_HEL = [(0, 1, 0.030, -0.028), (1, 1, -0.024, 0.021), (2, 1, 0.006, -0.005)]
SPECTRUM = _BASE + [(m, n, R * SHAPE, Z * SHAPE) for (m, n, R, Z) in _HEL]

# polarization modes (a,b,c)
MODES = {"iso": (1 / 3, 1 / 3, 1 / 3), "A": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)}

# square-cross-section torus wall (R0=1 normalized; encloses the toy plasma w/ gap)
R_IN, R_OUT, Z_W = 0.74, 1.30, 0.34


def loop_RZ(theta_s, phi, rho):
    """R(phi), Z(phi) of the loop at poloidal angle theta_s on flux surface rho.
    n!=0 (shaping) and the m>=1 (minor-radius) harmonics scale with rho; the n=0
    major radius R00 (the axis) stays put. phi is an array."""
    R = np.zeros_like(phi); Z = np.zeros_like(phi)
    for m, n, Rmn, Zmn in SPECTRUM:
        s = 1.0 if (m == 0 and n == 0) else rho
        ang = m * theta_s - n * NFP * phi
        R += (Rmn * s) * np.cos(ang)
        Z += (Zmn * s) * np.sin(ang)
    return R, Z


def field_angles(theta_s, phi):
    """Toy field angles alpha(phi), beta(phi) for the loop at theta_s (tier3 model)."""
    a = 0.5 * np.sin(theta_s) + 0.0 * phi
    b = np.arctan(IOTA * 0.6) + 0.0 * phi
    fr = FIELD_RIPPLE
    a1c, a1s = [fr * 0.6, fr * 0.2], [fr * 0.3, -fr * 0.1]
    b1c, b1s = [-fr * 0.5, fr * 0.15], [fr * 0.35, -fr * 0.1]
    for i in range(2):
        k = (i + 1) * NFP
        a = a + a1c[i] * np.cos(k * phi) + a1s[i] * np.sin(k * phi)
        b = b + b1c[i] * np.cos(k * phi) + b1s[i] * np.sin(k * phi)
    return a, b


def field_bhat(theta_s, phi):
    """Unit field B0(alpha,beta,phi) (Cartesian) along the loop (anarrima Eq.12)."""
    a, b = field_angles(theta_s, phi)
    Bx = np.cos(a) * np.sin(b) * np.cos(phi) - np.cos(b) * np.sin(phi)
    By = np.cos(b) * np.cos(phi) + np.cos(a) * np.sin(b) * np.sin(phi)
    Bz = -np.sin(a) * np.sin(b)
    return np.stack([Bx, By, Bz], axis=-1)


def loops():
    """List of (rho, theta_s, weight) source loops."""
    out = []
    for rho, wr in zip(RHO_SURFACES, W_RHO):
        for j in range(N_THETA_S):
            out.append((rho, 2 * np.pi * j / N_THETA_S, wr / N_THETA_S))
    return out


def wall_patches(n_per_wall=16):
    """Wall patch list: dict(R, Z, psi) with psi the INWARD-normal angle
    ((cos psi, sin psi) = inward normal). Square torus: inboard/outboard verticals,
    floor/ceiling horizontals. 's' = signed poloidal coord for plotting."""
    P = []
    zt = np.linspace(-Z_W, Z_W, n_per_wall)
    rt = np.linspace(R_IN, R_OUT, n_per_wall)
    for z in zt:
        P.append(dict(wall="inboard", R=R_IN, Z=z, psi=0.0, s=z))       # inward = +R
    for z in zt:
        P.append(dict(wall="outboard", R=R_OUT, Z=z, psi=np.pi, s=z))   # inward = -R
    for r in rt:
        P.append(dict(wall="floor", R=r, Z=-Z_W, psi=np.pi / 2, s=r))   # inward = +z
    for r in rt:
        P.append(dict(wall="ceiling", R=r, Z=Z_W, psi=-np.pi / 2, s=r))  # inward = -z
    return P


if __name__ == "__main__":  # self-checks
    ph = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    # field is unit
    b = field_bhat(1.0, ph)
    print(f"field |B-hat| max|n-1| = {np.max(np.abs(np.linalg.norm(b, axis=1) - 1)):.2e}")
    # plasma extent vs wall (must be enclosed with a gap)
    Rmin = Rmax = Zmin = Zmax = None
    for rho, th, _ in loops():
        R, Z = loop_RZ(th, ph, rho)
        Rmin = R.min() if Rmin is None else min(Rmin, R.min())
        Rmax = R.max() if Rmax is None else max(Rmax, R.max())
        Zmax = max(abs(Z).max(), Zmax or 0)
    print(f"plasma R in [{Rmin:.3f},{Rmax:.3f}], |Z|<{Zmax:.3f}; "
          f"wall R[{R_IN},{R_OUT}] |Z|<{Z_W} (enclosed: {R_IN < Rmin and Rmax < R_OUT and Zmax < Z_W})")
    # eps_eff (excursion / standoff) sanity ~ small
    eps = []
    for rho, th, _ in loops():
        R, Z = loop_RZ(th, ph, rho)
        p, zc = R.mean(), Z.mean()
        exc = np.max(np.sqrt((R - p) ** 2 + (Z - zc) ** 2))
        r = R_OUT; standoff = np.min(np.sqrt(p ** 2 + r ** 2 + (zc) ** 2 - 2 * p * r * np.cos(ph)))
        eps.append(exc / standoff)
    print(f"eps_eff (outboard-wall) max={max(eps):.3f} (want small, ~0.1)")
