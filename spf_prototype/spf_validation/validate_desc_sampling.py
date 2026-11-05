"""End-to-end integration proof: a DESC-built polarized source sampled through the
COMPILED OpenMC engine produces correctly-shaped directions about the DESC field.

Chain exercised: from_desc(eq, polarization) -> XML -> C++ StellaratorSource
constructor -> C++ field_direction (from the DESC-derived lmns/iota) -> C++
sample_polarized_direction. We sample real births with openmc, invert each
Cartesian birth back to flux coordinates (s, theta, zeta) with a vectorized
Newton on the geometry Fourier series, evaluate the local B-hat from the same
field data, and check <cos^2(angle(u, B-hat))> against the analytic mode value
(1/5 perp, 7/15 par, 1/3 unpolarized). If the compiled sampler used the DESC
field correctly, the moment matches.

Run with the descenv312 interpreter (DESC + the built openmc importable).
"""
from __future__ import annotations
import sys
import numpy as np

import openmc
from fluxrep import SpectralRep
from desc.examples import get


def _src_to_wout(src, nfp):
    return dict(xm=np.asarray(src.mode_m, float), xn=np.asarray(src.mode_n, float) * nfp,
                nfp=nfp, ns=len(src.field_rho), s=np.asarray(src.field_rho),
                rmnc=np.asarray(src.rmnc), zmns=np.asarray(src.zmns),
                rmns=None if src.rmns is None else np.asarray(src.rmns),
                zmnc=None if src.zmnc is None else np.asarray(src.zmnc),
                lmns=np.asarray(src.lmns),
                lmnc=None if src.lmnc is None else np.asarray(src.lmnc),
                iotaf=np.asarray(src.iota))


def invert_flux_coords(spec, R_t, Z_t, zeta, iters=60):
    """Vectorized Newton: find (s, theta) with R(s,theta,zeta)=R_t, Z=Z_t (cm)."""
    xm, xn = spec.xm, spec.xn
    # axis (s->0) major radius for the initial guess
    s = np.full_like(R_t, 0.5)
    R0, Z0, _ = spec.eval(np.full_like(R_t, 1e-3), np.zeros_like(R_t), zeta, want_sqrtg=False)
    theta = np.arctan2(Z_t - Z0, R_t - R0)
    for _ in range(iters):
        cos = np.cos(np.multiply.outer(theta, xm) - np.multiply.outer(zeta, xn))
        sin = np.sin(np.multiply.outer(theta, xm) - np.multiply.outer(zeta, xn))
        rc = spec._rmnc(s); zs = spec._zmns(s)
        R = np.einsum("ij,ij->i", rc, cos); Z = np.einsum("ij,ij->i", zs, sin)
        R_th = np.einsum("ij,ij->i", rc, -xm * sin); Z_th = np.einsum("ij,ij->i", zs, xm * cos)
        R_s = np.einsum("ij,ij->i", spec._rmnc(s, 1), cos)
        Z_s = np.einsum("ij,ij->i", spec._zmns(s, 1), sin)
        if spec.asym:
            rs = spec._rmns(s); zc = spec._zmnc(s)
            R += np.einsum("ij,ij->i", rs, sin); Z += np.einsum("ij,ij->i", zc, cos)
            R_th += np.einsum("ij,ij->i", rs, xm * cos); Z_th += np.einsum("ij,ij->i", zc, -xm * sin)
            R_s += np.einsum("ij,ij->i", spec._rmns(s, 1), sin)
            Z_s += np.einsum("ij,ij->i", spec._zmnc(s, 1), cos)
        dR = R - R_t; dZ = Z - Z_t
        det = R_s * Z_th - R_th * Z_s
        det = np.where(np.abs(det) < 1e-12, 1e-12, det)
        ds = -(Z_th * dR - R_th * dZ) / det
        dth = -(-Z_s * dR + R_s * dZ) / det
        s = np.clip(s + ds, 1e-4, 1.0); theta = theta + dth
    resid = np.hypot(dR, dZ)
    return s, np.mod(theta, 2 * np.pi), resid


def run(name, polarization, expected, n_rho=24, n=30000):
    eq = get(name); nfp = eq.NFP
    src = openmc.StellaratorSource.from_desc(
        eq, emission_density=np.ones(n_rho),
        energy=openmc.stats.delta_function(14.07e6),
        n_rho=n_rho, polarization=polarization)
    sphere = openmc.Sphere(r=1.0e4, boundary_type='vacuum')
    model = openmc.Model(
        geometry=openmc.Geometry([openmc.Cell(region=-sphere)]),
        settings=openmc.Settings(particles=100, batches=1,
                                 run_mode='fixed source', source=src))
    sites = model.sample_external_source(n)
    xyz = np.array([s.r for s in sites]); u = np.array([s.u for s in sites])

    spec = SpectralRep(_src_to_wout(src, nfp))
    R_t = np.hypot(xyz[:, 0], xyz[:, 1]); Z_t = xyz[:, 2]
    zeta = np.arctan2(xyz[:, 1], xyz[:, 0]) % (2 * np.pi)
    s, theta, resid = invert_flux_coords(spec, R_t, Z_t, zeta)
    ok = resid < 1e-2                       # keep well-converged inversions (cm)
    bhat = spec.eval_bhat(s[ok], theta[ok], zeta[ok])
    cos2 = np.einsum('ij,ij->i', u[ok], bhat) ** 2
    norm = np.linalg.norm(u, axis=1)
    print(f"DESC {name} pol={polarization}: <cos^2>={cos2.mean():.4f} "
          f"(expect {expected:.4f})  inverted {ok.sum()}/{n}  "
          f"median resid {np.median(resid):.1e} cm  |u|-1 max {np.abs(norm-1).max():.1e}")
    return cos2.mean()


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'precise_QA'
    run(name, (1.0, 0.0, 0.0), 1.0 / 5.0)      # perp
    run(name, (0.0, 1.0, 0.0), 7.0 / 15.0)     # par
    run(name, (1 / 3, 1 / 3, 1 / 3), 1.0 / 3.0)  # unpolarized (isotropic)


if __name__ == "__main__":
    main()
