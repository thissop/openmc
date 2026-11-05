"""Validate a DESC-native field reconstruction against DESC's own compute('B').

Mirrors validate_vmec_bhat.py but for DESC: extract the lambda stream-function
Fourier-Zernike coefficients and the iota profile from a DESC equilibrium, fold
them (Zernike-collapse + Ptolemy) into the same combined VMEC-style form the
geometry uses, reconstruct the field-line direction with the SAME formula as the
compiled/analytic sampler (B^theta = iota - lambda_zeta, B^zeta = 1 + lambda_theta),
and compare its direction to DESC's independently computed B at the same points.

Run with the descenv312 interpreter (DESC + openmc importable). This is the
convention gate: if the angular error is small, the DESC lambda/iota extraction
uses the same theta/zeta convention as the geometry and the sampler is correct.
"""
from __future__ import annotations
import sys
import numpy as np

import openmc
from openmc.source import StellaratorSource as SS
from fluxrep import SpectralRep

from desc.examples import get
from desc.grid import Grid, LinearGrid


def desc_combined(eq, n_rho=16):
    """Return a wout-style dict (xm,xn,rmnc,zmns,rmns,zmnc,lmns,lmnc,iotaf,ns) built
    from a DESC equilibrium on a uniform rho grid, on the UNION of R,Z,L modes."""
    rho = np.linspace(0.0, 1.0, n_rho)
    r_modes, r_lmn, z_modes, z_lmn, nfp = SS._desc_spectral_data(eq)
    l_modes = np.asarray(eq.L_basis.modes, dtype=int)
    l_lmn = np.asarray(eq.L_lmn, dtype=float)

    def collapse(modes, lmn):
        return np.column_stack([x * SS._zernike_radial(rho, l, m)
                                for (l, m, n), x in zip(modes, lmn)])
    radial_r = collapse(r_modes, r_lmn)
    radial_z = collapse(z_modes, z_lmn)
    radial_l = collapse(l_modes, l_lmn)

    surf = []
    for i in range(n_rho):
        tr, tz, tl = {}, {}, {}
        SS._desc_to_combined(r_modes[:, 1], r_modes[:, 2], radial_r[i], tr)
        SS._desc_to_combined(z_modes[:, 1], z_modes[:, 2], radial_z[i], tz)
        SS._desc_to_combined(l_modes[:, 1], l_modes[:, 2], radial_l[i], tl)
        surf.append((tr, tz, tl))
    modes = sorted({k for tr, tz, tl in surf for k in (*tr, *tz, *tl)})
    mode_m = np.array([m for m, _ in modes]); mode_n = np.array([n for _, n in modes])
    nm = len(modes)
    rmnc = np.zeros((n_rho, nm)); rmns = np.zeros((n_rho, nm))
    zmnc = np.zeros((n_rho, nm)); zmns = np.zeros((n_rho, nm))
    lmnc = np.zeros((n_rho, nm)); lmns = np.zeros((n_rho, nm))
    for i, (tr, tz, tl) in enumerate(surf):
        for k, key in enumerate(modes):
            if key in tr: rmnc[i, k], rmns[i, k] = tr[key]
            if key in tz: zmnc[i, k], zmns[i, k] = tz[key]
            if key in tl: lmnc[i, k], lmns[i, k] = tl[key]
    # iota(rho) from DESC
    lg = LinearGrid(rho=rho, M=0, N=0, NFP=nfp)
    di = np.asarray(eq.compute('iota', grid=lg)['iota'])
    rr = np.asarray(lg.nodes)[:, 0]
    iotaf = np.array([di[np.argmin(np.abs(rr - r))] for r in rho])

    sym = not (np.any(rmns) or np.any(zmnc))
    return dict(xm=mode_m.astype(float), xn=(mode_n * nfp).astype(float),
                nfp=nfp, ns=n_rho, s=rho**2,
                rmnc=rmnc * 100.0, zmns=zmns * 100.0,
                rmns=None if sym else rmns * 100.0,
                zmnc=None if sym else zmnc * 100.0,
                lmns=lmns, lmnc=None if not np.any(lmnc) else lmnc,
                iotaf=iotaf), rho, nfp


def validate(name, n_rho, n=4000, seed=3):
    eq = get(name)
    wout, rho, nfp = desc_combined(eq, n_rho=n_rho)
    wout['s'] = rho                          # grid is uniform in rho already
    spec = SpectralRep(wout)
    rng = np.random.default_rng(seed)
    s = rng.uniform(0.15, 0.95, n)
    th = rng.uniform(0, 2*np.pi, n)
    ze = rng.uniform(0, 2*np.pi/nfp, n)
    bhat = spec.eval_bhat(s, th, ze)
    # DESC oracle: B in (R, phi, Z) -> Cartesian via zeta=phi
    B = np.asarray(eq.compute(['B'], grid=Grid(np.column_stack([s, th, ze])))['B'])
    cz, sz = np.cos(ze), np.sin(ze)
    Bxyz = np.column_stack([B[:, 0]*cz - B[:, 1]*sz, B[:, 0]*sz + B[:, 1]*cz, B[:, 2]])
    Bxyz /= np.linalg.norm(Bxyz, axis=1, keepdims=True)
    deg = np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i', bhat, Bxyz)), 0, 1)))
    return dict(example=name, nfp=int(nfp), n_rho=int(n_rho), samples=int(n),
                median_degrees=float(np.median(deg)),
                p95_degrees=float(np.percentile(deg, 95)),
                p99_degrees=float(np.percentile(deg, 99)),
                max_degrees=float(deg.max()))


def main():
    import json
    examples = sys.argv[1:] or ['precise_QA', 'precise_QH']
    out = []
    for name in examples:
        for n_rho in (16, 33, 49):
            r = validate(name, n_rho)
            out.append(r)
            print(f"DESC {name} n_rho={n_rho} (nfp={r['nfp']}): field-dir error vs "
                  f"compute('B') -- median {r['median_degrees']:.3e} p99 "
                  f"{r['p99_degrees']:.3e} max {r['max_degrees']:.3e} deg")
    with open('desc_bhat_validation.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("wrote desc_bhat_validation.json")


if __name__ == "__main__":
    main()
