#!/usr/bin/env python
"""PRODUCER (DESC; login node / desc venv): build a fixed-boundary DESC equilibrium from
a QUASR device's plasma boundary and write its field map (data/quasr<ID>_equil.meta/.bin)
in the portable spf_fieldmap_v1 format the compiled source, python/fieldmap.py, and
python/equilibrium_field.py all consume UNCHANGED.

WHY (G3): the vacuum COIL field carries each device's real B.n residual on the LCFS, so
field_audit excludes the worst-fit devices -- exactly the high-iota QH configs with
toroidal-sense reversal (e.g. 1190023), where the first-moment C decouples from the
parity-correct nematic order S_phi (docs/THEORY.md, SCOPE_AND_CAVEATS.md). The
EQUILIBRIUM field is a flux function (B.n = 0 on the LCFS by construction), so those
devices pass the audit and can adjudicate C-vs-tensor.

This MIRRORS desc_to_fieldmap.py (the precise_QA producer); the only new part is building
the equilibrium from the QUASR boundary instead of a DESC example. Kept separate so the
working precise_QA producer is untouched.

RUN (DESC pins conflict with OpenMC -> separate venv):
    $HOME/desc_venv/bin/python spf_prototype/python/quasr_equilibrium_field.py <ID>

NOT validatable on the aarch64/osx dev box (no DESC). VERSION-SENSITIVE SPOTS to check on
the DESC machine, flagged inline like the stl_to_h5m seam:
  [C1] QUASR-boundary -> DESC FourierRZToroidalSurface convention (mode signs / NFP);
  [C2] eq.solve() resolution + convergence for aggressive high-iota QH boundaries;
  [C3] B.n on the LCFS of the SOLVED eq must be ~0 -- field_audit re-checks it downstream.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "sweep"))
DATADIR = HERE.parent / "data"
LENGTH_SCALE_CM = 100.0


def _quasr_boundary_modes(ID, ntheta=64, nphi=128, mpol=8, ntor=8, flip_theta=True):
    """Evaluate the QUASR LCFS on a (theta, phi) grid (pure-numpy quasr_geom) and FIT real
    Fourier modes R_mn, Z_mn. Evaluate-then-fit keeps this INDEPENDENT of quasr_geom's
    internal mode convention; the only convention we then own is the fit basis
    cos/sin(m*theta - n*NFP*phi), mapped to DESC in _desc_surface(). [C1]"""
    import quasr_loader as ql
    dev = ql.load_device(ID, coils=False)
    nfp = int(dev.meta["nfp"])
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    # [C1] QUASR's LCFS parameterization winds theta CLOCKWISE (negative poloidal Jacobian -- the m=1
    # R mode comes out negative), so DESC's ensure_positive_jacobian rejects the axisymmetric seed.
    # Evaluate at -theta so the fitted modes wind counter-clockwise (positive Jacobian). This is a pure
    # reparameterization of the SAME physical boundary; gate V4 (rho=1 reconstructs the fitted boundary,
    # in this same convention) verifies it downstream.
    TH_eval = (-TH) % (2.0 * np.pi) if flip_theta else TH
    R, Z = dev.device.RZ(TH_eval, PH, 1.0)                 # LCFS, native units
    Rmodes, Zmodes = {}, {}
    for m in range(mpol + 1):
        for n in range(-ntor, ntor + 1):
            ang = m * TH - n * nfp * PH
            norm = (1.0 if (m == 0 and n == 0) else 2.0) / (ntheta * nphi)
            c, s = np.cos(ang), np.sin(ang)
            Rmodes[(m, n)] = (float((R * c).sum() * norm), float((R * s).sum() * norm))
            Zmodes[(m, n)] = (float((Z * c).sum() * norm), float((Z * s).sum() * norm))
    return nfp, Rmodes, Zmodes, dev.meta


def _desc_surface(nfp, Rmodes, Zmodes):
    """Map fitted (cos, sin) modes to a DESC FourierRZToroidalSurface (stellarator-
    symmetric: R even/cos, Z odd/sin). [C1] Confirm the DESC toroidal-angle sign
    (n*NFP*phi) matches this basis on the DESC machine; a B.n~0 audit downstream (C3)
    catches a mismatch."""
    from desc.geometry import FourierRZToroidalSurface
    R_lmn, R_modes, Z_lmn, Z_modes = [], [], [], []
    for (m, n), (rc, _rs) in Rmodes.items():
        if abs(rc) > 1e-10:
            R_lmn.append(rc); R_modes.append([m, n])          # cos part
    for (m, n), (_zc, zs) in Zmodes.items():
        if abs(zs) > 1e-10:
            Z_lmn.append(zs); Z_modes.append([m, n])          # sin part
    return FourierRZToroidalSurface(
        R_lmn=R_lmn, Z_lmn=Z_lmn,
        modes_R=np.array(R_modes), modes_Z=np.array(Z_modes), NFP=nfp, sym=True)


def build(ID, stem=None, n_rho=12, n_theta=48, n_zeta=144, nR=40, nphi=48, nZ=40, idw_k=8):
    from desc.equilibrium import Equilibrium
    from desc.grid import LinearGrid
    from scipy.spatial import cKDTree

    stem = stem or f"quasr{int(ID)}_equil"
    nfp, Rmodes, Zmodes, meta = _quasr_boundary_modes(ID)
    print(f"[quasr_equil] device {ID}: NFP={nfp}; solving fixed-boundary equilibrium ...")
    surf = _desc_surface(nfp, Rmodes, Zmodes)
    eq = Equilibrium(surface=surf, NFP=nfp, L=6, M=6, N=6, sym=True)    # [C2] resolution
    eq.solve(verbose=1, ftol=1e-2, maxiter=50)                          # [C2] convergence

    # sample b_hat over the FULL torus + IDW onto a regular (R,phi,Z) grid, then write
    # spf_fieldmap_v1 (mirrors desc_to_fieldmap.build; full torus because Cartesian B is
    # NOT field-period periodic -- see that file's docstring).
    g = LinearGrid(rho=np.linspace(0.05, 1.0, n_rho), theta=n_theta, zeta=n_zeta,
                   NFP=1, sym=False)
    d = eq.compute(["R", "phi", "Z", "B"], grid=g, basis="rpz")
    R = np.asarray(d["R"]); phi = np.asarray(d["phi"]); Z = np.asarray(d["Z"])
    Bcyl = np.asarray(d["B"])
    cph, sph = np.cos(phi), np.sin(phi)
    Bx = Bcyl[:, 0] * cph - Bcyl[:, 1] * sph
    By = Bcyl[:, 0] * sph + Bcyl[:, 1] * cph
    Bz = Bcyl[:, 2]
    Bn = np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    Bsamp = np.column_stack([Bx / Bn, By / Bn, Bz / Bn])
    pts = np.column_stack([R * cph * LENGTH_SCALE_CM, R * sph * LENGTH_SCALE_CM,
                           Z * LENGTH_SCALE_CM])
    Rcm = R * LENGTH_SCALE_CM
    Rmin, Rmax = float(Rcm.min()) * 0.97, float(Rcm.max()) * 1.03
    Zmin = float(Z.min()) * LENGTH_SCALE_CM * 1.05
    Zmax = float(Z.max()) * LENGTH_SCALE_CM * 1.05
    period = 2.0 * np.pi
    Rg = np.linspace(Rmin, Rmax, nR); Zg = np.linspace(Zmin, Zmax, nZ)
    pg = np.arange(nphi) * (period / nphi)
    tree = cKDTree(pts)
    Bx_o = np.empty((nR, nphi, nZ)); By_o = np.empty_like(Bx_o); Bz_o = np.empty_like(Bx_o)
    for ip, ph_ in enumerate(pg):
        c, s = np.cos(ph_), np.sin(ph_)
        for iR, Rr in enumerate(Rg):
            q = np.column_stack([np.full(nZ, Rr * c), np.full(nZ, Rr * s), Zg])
            dist, idx = tree.query(q, k=idw_k)
            w = 1.0 / (dist ** 2 + 1e-9); w /= w.sum(axis=1, keepdims=True)
            bvec = (w[:, :, None] * Bsamp[idx]).sum(axis=1)
            nrm = np.linalg.norm(bvec, axis=1, keepdims=True)
            degen = (nrm[:, 0] < 1e-6); bvec[degen] = Bsamp[idx[degen, 0]]
            bvec = bvec / (np.linalg.norm(bvec, axis=1, keepdims=True) + 1e-30)
            Bx_o[iR, ip], By_o[iR, ip], Bz_o[iR, ip] = bvec[:, 0], bvec[:, 1], bvec[:, 2]

    DATADIR.mkdir(parents=True, exist_ok=True)
    meta_p = DATADIR / f"{stem}.meta"; bin_p = DATADIR / f"{stem}.bin"
    meta_p.write_text(
        "magic spf_fieldmap_v1\nndim 3\n"
        f"nR {nR}\nnphi {nphi}\nnZ {nZ}\n"
        f"R_min {Rmin!r}\nR_max {Rmax!r}\n"
        f"phi_min 0.0\nperiod {period!r}\n"
        f"Z_min {Zmin!r}\nZ_max {Zmax!r}\n"
        f"nfp {nfp}\nendian little\nsource quasr_equil:{ID}\n")
    with open(bin_p, "wb") as fh:
        for arr in (Bx_o, By_o, Bz_o):
            fh.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())
    print(f"[quasr_equil] wrote {meta_p.name} + {bin_p.name} ({nR}x{nphi}x{nZ}, nfp={nfp}). "
          "[C3] verify B.n~0 on the LCFS with field_audit before trusting this map.")
    return bin_p


if __name__ == "__main__":
    dev_id = int(sys.argv[1]) if len(sys.argv) > 1 else 1190023
    build(dev_id, sys.argv[2] if len(sys.argv) > 2 else None)
