"""Shared problem definition for the GENTLE published-QA free-streaming cross-check,
exposing the same interface as toy_qalow_config but backed by a real QUASR device
(geometry + flux-surface field from quasr_geom, pure numpy). Imported by BOTH the
analytic step (ana_venv) and the OpenMC step (spf_venv) via SPF_CONFIG=quasr_config.

Device 59509: nfp=3 quasi-axisymmetric, aspect ~6.7, the gentlest QA in our QUASR
shortlist (e_geom ~0.29; see data/quasr/eps_eff_results.csv) -- a real published
equilibrium in the analytic series-convergent regime, unlike precise_QA (eps_eff 2.3,
MC-only). The wall is a square-cross-section torus enclosing the boundary with a gap;
all non-axisymmetry is in the SOURCE. We compare DIRECTIONALITY (A/iso, B/iso per wall
patch), which isolates the polarized emission about the local field.
"""
from __future__ import annotations

import numpy as np

import quasr_geom as G

STEM = "quasr59509"
TITLE = "QUASR QA 59509"
DEVICE_ID = 59509
IOTA = 0.1                         # mean_iota from the QUASR catalogue
NFP = 3

RHO_SURFACES = [0.6, 1.0]
W_RHO = [0.6 ** 2, 1.0 ** 2]       # rho^2 emissivity weight
N_THETA_S = 14
GAP_A = 0.40                       # square wall placed GAP_A * a outside the bbox (eps_eff ~0.43)
K_TRIG = 34                        # analytic trig order (>= max |n|*nfp = 30)

MODES = {"iso": (1 / 3, 1 / 3, 1 / 3), "A": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)}

_DEV = G.load_device(DEVICE_ID, IOTA)

# --- square wall from the boundary (rho=1) bounding box + gap ---
_th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
_ph = np.linspace(0, 2 * np.pi, 128, endpoint=False)
_TH, _PH = np.meshgrid(_th, _ph, indexing="ij")
_R, _Z = _DEV.RZ(_TH, _PH, 1.0)
_A = 0.5 * (_R.max() - _R.min())
R_IN = float(_R.min() - GAP_A * _A)
R_OUT = float(_R.max() + GAP_A * _A)
Z_W = float(np.abs(_Z).max() + GAP_A * _A)


def loop_RZ(theta_s, phi, rho):
    """R(phi), Z(phi) of the loop at poloidal angle theta_s on flux surface rho.
    theta_s may be a scalar (one loop) or an array broadcast against phi (vectorized
    presampling, one theta per particle)."""
    phi = np.asarray(phi, dtype=float)
    th = np.broadcast_to(np.asarray(theta_s, dtype=float), phi.shape)
    return _DEV.RZ(th, phi, rho)


def field_bhat(theta_s, phi):
    """Unit flux-surface-tangent field along the loop (boundary field, rho=1),
    Cartesian (...,3) -- identical on both validation sides. theta_s scalar or array."""
    phi = np.asarray(phi, dtype=float)
    th = np.broadcast_to(np.asarray(theta_s, dtype=float), phi.shape)
    return _DEV.bhat(th, phi, 1.0)


def loops():
    out = []
    for rho, wr in zip(RHO_SURFACES, W_RHO):
        for j in range(N_THETA_S):
            out.append((rho, 2 * np.pi * j / N_THETA_S, wr / N_THETA_S))
    return out


def wall_patches(n_per_wall=16):
    P = []
    zt = np.linspace(-Z_W, Z_W, n_per_wall)
    rt = np.linspace(R_IN, R_OUT, n_per_wall)
    for z in zt:
        P.append(dict(wall="inboard", R=R_IN, Z=z, psi=0.0, s=z))
    for z in zt:
        P.append(dict(wall="outboard", R=R_OUT, Z=z, psi=np.pi, s=z))
    for r in rt:
        P.append(dict(wall="floor", R=r, Z=-Z_W, psi=np.pi / 2, s=r))
    for r in rt:
        P.append(dict(wall="ceiling", R=r, Z=Z_W, psi=-np.pi / 2, s=r))
    return P


if __name__ == "__main__":  # self-checks
    eg, ee, a = G.eps_eff(_DEV, gap=GAP_A)
    b = field_bhat(1.0, _ph)
    print(f"device {DEVICE_ID}: nfp {NFP} R0 {_DEV.R0:.3f} a {_A:.3f}")
    print(f"wall R[{R_IN:.3f},{R_OUT:.3f}] |Z|<{Z_W:.3f}")
    print(f"plasma enclosed: {R_IN < _R.min() and _R.max() < R_OUT and np.abs(_Z).max() < Z_W}")
    print(f"e_geom {eg:.3f}, eps_eff(gap={GAP_A}) {ee:.3f}")
    print(f"field |B-hat|-1 max {np.max(np.abs(np.linalg.norm(b, axis=1) - 1)):.1e}")
