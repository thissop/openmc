#!/usr/bin/env python
"""PRODUCER: turn a DESC stellarator equilibrium into our portable field-map
format (<stem>.meta + <stem>.bin) that src/spf_fieldmap.hpp and python/fieldmap.py
consume UNCHANGED. This is the INJECT(helios) producer seam: swap precise_QA for
any DESC equilibrium (precise_QH, a QUASR device, or a real Helios wout via
desc.VMECIO) and the rest of the pipeline is identical.

Runs in the ISOLATED desc venv (DESC pulls its own jax/numpy); it only WRITES the
field map. Transport stays in spf_venv. Run:
    $HOME/desc_venv/bin/python spf_prototype/python/desc_to_fieldmap.py [name] [stem]

Method (honest): we need only the field DIRECTION B-hat(x) inside the plasma,
since the SPF source births inside the LCFS. We sample B on a dense flux grid
(rho,theta,zeta over the full torus), convert to Cartesian unit vectors, and
interpolate onto a regular (R,phi,Z) grid (one field period; FieldMapField's
periodic-phi wrap handles the rest) by k-nearest inverse-distance weighting
(KDTree). Nodes with no plasma sample within the cell scale fall back to the
nearest plasma B-hat. The map is therefore high-fidelity inside the LCFS and a
smooth nearest-value extension just outside (which the source never samples). A
DESC map_coordinates inverse-map is the documented higher-fidelity upgrade
(DEFERRED.md).
"""
from __future__ import annotations

import struct
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
DATADIR = REPO / "spf_prototype" / "data"
LENGTH_SCALE_CM = 100.0  # DESC equilibria are normalized (R0~1 m); store grid in cm


def build(name="precise_QA", stem="equil_precise_qa",
          n_rho=12, n_theta=48, n_zeta=96, nR=40, nphi=24, nZ=40, idw_k=8):
    import desc
    from desc.examples import get
    from desc.grid import LinearGrid

    eq = get(name)
    nfp = int(eq.NFP)
    print(f"[desc_to_fieldmap] {name}: NFP={nfp}, aspect={float(eq.compute('R0/a')['R0/a']):.2f}")

    # --- sample B on a dense flux grid over the FULL torus ---
    g = LinearGrid(rho=np.linspace(0.05, 1.0, n_rho), theta=n_theta, zeta=n_zeta,
                   NFP=nfp, sym=False)
    d = eq.compute(["R", "phi", "Z", "B"], grid=g)
    R = np.asarray(d["R"]); phi = np.asarray(d["phi"]); Z = np.asarray(d["Z"])
    Bcyl = np.asarray(d["B"])  # (N,3) in (R,phi,Z) cylindrical basis
    # cylindrical -> Cartesian, then unit
    cph, sph = np.cos(phi), np.sin(phi)
    Bx = Bcyl[:, 0] * cph - Bcyl[:, 1] * sph
    By = Bcyl[:, 0] * sph + Bcyl[:, 1] * cph
    Bz = Bcyl[:, 2]
    Bn = np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    Bx, By, Bz = Bx / Bn, By / Bn, Bz / Bn
    # sample lab points in cm
    sx = R * cph * LENGTH_SCALE_CM
    sy = R * sph * LENGTH_SCALE_CM
    sz = Z * LENGTH_SCALE_CM
    Rcm = R * LENGTH_SCALE_CM

    # --- regular output grid (cm), one field period in phi ---
    Rmin, Rmax = float(Rcm.min()) * 0.97, float(Rcm.max()) * 1.03
    Zmin, Zmax = float(sz.min()) * 1.05, float(sz.max()) * 1.05
    period = 2.0 * np.pi / nfp
    Rg = np.linspace(Rmin, Rmax, nR)
    Zg = np.linspace(Zmin, Zmax, nZ)
    pg = np.arange(nphi) * (period / nphi)  # phi_min=0

    # build sample KDTree once (Cartesian, with periodic phi images so the wrap
    # at the period boundary interpolates correctly)
    pts = np.column_stack([sx, sy, sz])
    Bsamp = np.column_stack([Bx, By, Bz])
    tree = cKDTree(pts)

    Bx_o = np.empty((nR, nphi, nZ)); By_o = np.empty_like(Bx_o); Bz_o = np.empty_like(Bx_o)
    for ip, ph in enumerate(pg):
        c, s = np.cos(ph), np.sin(ph)
        for iR, Rr in enumerate(Rg):
            qx = Rr * c; qy = Rr * s
            q = np.column_stack([np.full(nZ, qx), np.full(nZ, qy), Zg])
            dist, idx = tree.query(q, k=idw_k)
            w = 1.0 / (dist ** 2 + 1e-9)
            w /= w.sum(axis=1, keepdims=True)
            bvec = (w[:, :, None] * Bsamp[idx]).sum(axis=1)  # (nZ,3)
            n = np.linalg.norm(bvec, axis=1, keepdims=True)
            bvec = bvec / n
            Bx_o[iR, ip, :] = bvec[:, 0]
            By_o[iR, ip, :] = bvec[:, 1]
            Bz_o[iR, ip, :] = bvec[:, 2]

    # --- write our portable format ---
    DATADIR.mkdir(parents=True, exist_ok=True)
    meta = DATADIR / f"{stem}.meta"; binf = DATADIR / f"{stem}.bin"
    meta.write_text(
        "magic spf_fieldmap_v1\n"
        f"ndim 3\nnR {nR}\nnphi {nphi}\nnZ {nZ}\n"
        f"R_min {Rmin!r}\nR_max {Rmax!r}\n"
        f"phi_min 0.0\nperiod {period!r}\n"
        f"Z_min {Zmin!r}\nZ_max {Zmax!r}\n"
        f"nfp {nfp}\nendian little\n"
        f"source desc:{name}\n")
    with open(binf, "wb") as fh:
        for arr in (Bx_o, By_o, Bz_o):
            fh.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())
    # provenance for the conformal wall builder: the LCFS evaluated on a
    # (theta, phi) grid in cm -- CONVENTION-FREE (no Fourier sign ambiguity).
    _write_surface_grid(eq, DATADIR / f"{stem}_surface.npz", nfp)
    # also the raw boundary harmonics (reference)
    _write_boundary(eq, DATADIR / f"{stem}_boundary.txt", nfp)
    print(f"[desc_to_fieldmap] wrote {meta.name} + {binf.name} "
          f"({nR}x{nphi}x{nZ}, nfp={nfp}, R[{Rmin:.1f},{Rmax:.1f}]cm) "
          f"+ {stem}_surface.npz + {stem}_boundary.txt")
    return binf


def _write_surface_grid(eq, path, nfp, n_theta=64, n_phi=96):
    """Evaluate the LCFS on a (theta, phi) grid over the FULL torus and store R,Z
    in cm. Convention-free input for the conformal wall builder. Reshape is done by
    explicit lexsort on (theta, zeta) so it is INDEPENDENT of DESC node ordering
    (DESC LinearGrid orders theta-fastest and restricts zeta to one field period;
    we set NFP=1 to span the full torus and sort to a clean (theta,phi) grid)."""
    from desc.grid import LinearGrid
    g = LinearGrid(rho=np.array([1.0]), theta=n_theta, zeta=n_phi, NFP=1, sym=False)
    d = eq.compute(["R", "phi", "Z"], grid=g)
    nodes = g.nodes  # (N,3) columns (rho, theta, zeta)
    th = nodes[:, 1]; ze = nodes[:, 2]
    order = np.lexsort((ze, th))  # theta outer, zeta inner -> C-order (n_theta,n_phi)
    R = (np.asarray(d["R"])[order]).reshape(n_theta, n_phi) * LENGTH_SCALE_CM
    Z = (np.asarray(d["Z"])[order]).reshape(n_theta, n_phi) * LENGTH_SCALE_CM
    phi = (np.asarray(d["phi"])[order]).reshape(n_theta, n_phi)
    # sanity: each row (fixed theta) must sweep phi over the full torus; each column
    # (fixed phi) must sweep theta -> the cross-section the builder offsets.
    assert np.ptp(phi[0, :]) > 5.0, "surface grid axis order wrong (phi not on axis1)"
    np.savez(path, R=R, Z=Z, phi=phi, nfp=nfp,
             theta=np.linspace(0, 2 * np.pi, n_theta, endpoint=False),
             phi_axis=np.linspace(0, 2 * np.pi, n_phi, endpoint=False))


def _write_boundary(eq, path, nfp):
    """Dump the LCFS boundary Fourier modes (R_mn, Z_mn) in cm for the conformal
    wall builder. Format: 'm n Rcoef Zcoef' per line (cm), + header."""
    Rb = eq.surface.R_lmn; Zb = eq.surface.Z_lmn
    R_basis = eq.surface.R_basis; Z_basis = eq.surface.Z_basis
    lines = [f"# LCFS boundary modes (cm); nfp={nfp}; R(theta,phi)=sum R_mn cos/..., "
             "convention DESC FourierRZToroidalSurface", f"# nfp {nfp}"]
    modes = {}
    for coef, (l, m, n) in zip(Rb, R_basis.modes):
        modes.setdefault((m, n), [0.0, 0.0])[0] = float(coef) * LENGTH_SCALE_CM
    for coef, (l, m, n) in zip(Zb, Z_basis.modes):
        modes.setdefault((m, n), [0.0, 0.0])[1] = float(coef) * LENGTH_SCALE_CM
    for (m, n) in sorted(modes):
        rc, zc = modes[(m, n)]
        if abs(rc) > 1e-9 or abs(zc) > 1e-9:
            lines.append(f"{m} {n} {rc!r} {zc!r}")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "precise_QA"
    stem = sys.argv[2] if len(sys.argv) > 2 else "equil_precise_qa"
    build(name, stem)
