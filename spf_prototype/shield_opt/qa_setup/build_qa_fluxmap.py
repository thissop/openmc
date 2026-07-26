#!/usr/bin/env python
"""Turn a precise-QA VMEC `wout_*.nc` into our spf_fluxmap_v1 (the format the native
C++ StellaratorSource consumes), matching python/vmec_fluxmap.py's synthesis exactly.

WHY A SEPARATE SCRIPT (vs vmec_fluxmap.py):
  vmec_fluxmap.py couples the VMEC *solve* (vmecpp, from a QUASR boundary) to the
  fluxmap synthesis. For precise-QA we already have the fixed-boundary VMEC INDATA
  (shield_opt/data/raw/equilibria/input.LandremanPaul2021_QA) and only need to run
  VMEC once and read the resulting wout. This script consumes a wout .nc DIRECTLY,
  so it runs anywhere with numpy + a netCDF reader (no vmecpp/DESC needed locally).
  The Fourier synthesis, the sign(sqrt g) gate (V2) and the divergence-theorem volume
  gate (V1b) are byte-for-byte the same physics as vmec_fluxmap.build().

CANONICAL ANCHOR (verified against the repo, 2026-07):
  boundary : Landreman & Paul 2021 precise-QA (PRL 128, 035001), nfp=2, aspect=6,
             iota~0.42. Boundary already vendored at
             shield_opt/data/raw/equilibria/input.LandremanPaul2021_QA
             (VMEC INDATA: NFP=2, MPOL=16, NTOR=12, LASYM=F, unit scale R0~1.01 m,
             a~0.168 m, <B>~1 T).
  coils    : Wechsung 2022 (CoilsForPreciseQS) QA24 == modular coils that approximate
             THIS SAME Landreman-Paul QA boundary (README: "approximate the QA ...
             found in Landreman & Paul 2022"). QA24 = length-bound-24, order-16, 16
             physical coils. So input.LandremanPaul2021_QA and Wechsung_QA24 are a
             MATCHED boundary/coil pair -- no mismatch. Load coils via
             shield_opt/coil_adapter.py load_qa("QA24").

REACTOR SCALE (match the QH device, a=1.7 m):
  Native QA minor radius a=0.1683 m (coil_adapter.NATIVE_A["QA"]); reactor target
  a=1.7 m => scale = 1.7/0.1683 = 10.10 -> R0 ~ 10.2 m, a ~ 1.7 m. coil_adapter.load_qa
  applies this SAME factor to the coils by default, so pass --scale 10.10 here to keep
  the plasma geometry and the coils in one consistent reactor frame. The unit-scale
  map (scale=1.0) is the right input for the analytic anarrima cross-check
  (python/ana_precise_qa.py works in the native frame); the reactor-scaled map is the
  right input for the ParaStell build + OpenMC transport.

  NOTE ON UNITS: VMEC wout is SI (meters, Tesla); this script emits R,Z,phi in the
  native/scaled METERS of the wout times `scale`. The StellaratorSource / ParaStell
  build expect their own length unit (OpenMC = cm) -- apply the m->cm factor in the
  build step, consistently for BOTH the fluxmap and the coils. Only the field
  DIRECTION (b_hat) enters the polarized sampler, so scaling R,Z (and leaving B
  magnitude alone) does not change any sampled direction; scale is pure geometry.

USAGE
  # after VMEC has produced wout_LandremanPaul2021_QA.nc (see QA_STANDUP_PLAN.md):
  python build_qa_fluxmap.py wout_LandremanPaul2021_QA.nc --scale 10.10 \
         --out ../../data/precise_QA_reactor_fluxmap
  # unit-scale (for the analytic anarrima benchmark):
  python build_qa_fluxmap.py wout_LandremanPaul2021_QA.nc --scale 1.0 \
         --out ../../data/precise_QA_vmec_fluxmap

Emits <out>.npz + <out>.meta + <out>.bin (spf_fluxmap_v1).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# reuse the exact, audited writer from the QH pipeline
sys.path.insert(0, str(HERE.parent.parent / "python"))
from quasr_fluxmap import write_fluxmap_bin  # noqa: E402


def _read_wout(path):
    """Read a VMEC wout .nc into a dict of the fields vmec_fluxmap.py uses.
    Tries netCDF4 then scipy.io.netcdf_file (either is fine; wout is SI)."""
    keys = ("nfp", "ns", "rmnc", "zmns", "gmnc", "bsupumnc", "bsupvmnc",
            "xm", "xn", "xm_nyq", "xn_nyq", "volume_p")
    d = {}
    try:
        import netCDF4
        ds = netCDF4.Dataset(path, "r")
        for k in keys:
            if k in ds.variables:
                d[k] = np.asarray(ds.variables[k][:])
        ds.close()
    except Exception:
        from scipy.io import netcdf_file
        ds = netcdf_file(path, "r", mmap=False)
        for k in keys:
            if k in ds.variables:
                d[k] = np.asarray(ds.variables[k][()])
        ds.close()
    # scalars
    d["nfp"] = int(np.ravel(d["nfp"])[0])
    d["ns"] = int(np.ravel(d["ns"])[0])
    d["volume_p"] = float(np.ravel(d["volume_p"])[0])
    # Nyquist fallback (older wout files may omit *_nyq -> gmnc/bsup* share xm/xn)
    d.setdefault("xm_nyq", d["xm"])
    d.setdefault("xn_nyq", d["xn"])
    return d


def build(wout_path, out, scale=1.0, n_rho=12, n_theta=48, n_zeta=144,
          volume_tol=1.5e-2):
    """Synthesize spf_fluxmap_v1 from a VMEC wout. Mirrors vmec_fluxmap.build().

    scale : multiply R,Z (and hence phi-plane geometry) by this factor. Field
            components are NOT rescaled -- only b_hat DIRECTION matters to the
            sampler, and geometry scaling leaves directions invariant. The
            divergence-theorem volume gate is applied on the UNSCALED geometry so
            the tolerance is scale-independent.
    """
    w = _read_wout(wout_path)
    nfp = w["nfp"]
    xm = w["xm"].astype(float)
    xn = w["xn"].astype(float)                      # xn already includes nfp
    xmn = w["xm_nyq"].astype(float)                 # Nyquist set for gmnc/bsup*
    xnn = w["xn_nyq"].astype(float)
    ns = w["ns"]
    s_full = np.linspace(0.0, 1.0, ns)
    s_half = (np.arange(1, ns) - 0.5) / (ns - 1)
    rmnc = w["rmnc"]; zmns = w["zmns"]; gmnc = w["gmnc"]
    busc = w["bsupumnc"]; bvsc = w["bsupvmnc"]

    rho = np.linspace(1.0 / n_rho, 1.0, n_rho)
    s_grid = rho ** 2                               # VMEC s ~ toroidal flux ~ rho^2
    theta = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2 * np.pi, n_zeta, endpoint=False)
    TH, ZE = np.meshgrid(theta, zeta, indexing="ij")
    ang = xm[:, None, None] * TH[None] - xn[:, None, None] * ZE[None]
    cos, sin = np.cos(ang), np.sin(ang)
    angn = xmn[:, None, None] * TH[None] - xnn[:, None, None] * ZE[None]
    cosn = np.cos(angn)

    # VMEC: rmnc/zmns on FULL grid; gmnc/bsup* on HALF grid (skip s=0 column).
    ip_full = lambda C: np.array([np.interp(s_grid, s_full, C[:, k]) for k in range(C.shape[1])])
    ip_half = lambda C: np.array([np.interp(s_grid, s_half, C[1:, k]) for k in range(C.shape[1])])
    rc = ip_full(rmnc); zc = ip_full(zmns); gc = ip_half(gmnc)
    buc = ip_half(busc); bvc = ip_half(bvsc)

    sc = lambda Cs: np.einsum("ks,ktz->stz", Cs, cos)
    ss = lambda Cs: np.einsum("ks,ktz->stz", Cs, sin)
    scn = lambda Cs: np.einsum("ks,ktz->stz", Cs, cosn)
    R = sc(rc); Z = ss(zc); sg = scn(gc)
    Bu = scn(buc); Bv = scn(bvc)
    dRdt = ss(-xm[:, None] * rc); dRdz = ss(xn[:, None] * rc)
    dZdt = sc(xm[:, None] * zc); dZdz = sc(-xn[:, None] * zc)
    BR = Bu * dRdt + Bv * dRdz
    BZ = Bu * dZdt + Bv * dZdz
    Bphi = R * Bv                                   # VMEC zeta = cylindrical phi
    phi = np.broadcast_to(zeta[None, None, :], (n_rho, n_theta, n_zeta)).copy()

    # ---- V2: single-signed sqrt(g) (invariant -> hard assert) ----
    sgn = float(np.sign(np.median(sg)))
    bad = np.sign(sg) != sgn
    assert not np.any(bad) or np.all(np.abs(sg[bad]) < 1e-9 * np.abs(sg).max()), \
        f"V2 FAILED: sqrt(g) changes sign (min {sg.min():.3e}, max {sg.max():.3e})"
    sg = np.abs(sg)

    # ---- V1b: divergence-theorem volume from LCFS vs VMEC volume (unscaled geom) ----
    Rb, Zb, PHb = R[-1], Z[-1], phi[-1]
    Xb = np.stack([Rb * np.cos(PHb), Rb * np.sin(PHb), Zb], axis=-1)
    dth, dze = 2 * np.pi / n_theta, 2 * np.pi / n_zeta
    Xt = np.gradient(Xb, theta, axis=0); Xz = np.gradient(Xb, zeta, axis=1)
    V_div = abs(float((Xb * np.cross(Xt, Xz)).sum() * dth * dze / 3.0))
    V_vmec = w["volume_p"]
    relV = abs(V_div - V_vmec) / abs(V_vmec)
    print(f"[qaflux] V1b: div-theorem vol {V_div:.6g} vs VMEC {V_vmec:.6g}  rel {relV:.2e}",
          flush=True)
    assert relV < volume_tol, f"V1b FAILED rel={relV:.2e} (tol {volume_tol})"

    # ---- geometry scale (pure geometry; b_hat direction invariant) ----
    if scale != 1.0:
        R = R * scale
        Z = Z * scale
        print(f"[qaflux] applied reactor scale x{scale:g}: "
              f"R0~{0.5*(R.max()+R.min()):.3g} m, a~{0.5*(R.max()-R.min()):.3g} m", flush=True)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(str(out) + ".npz", rho=rho, theta=theta, zeta=zeta, nfp=nfp,
             sign_sqrtg=sgn, sqrtg=sg, R=R, Z=Z, phi=phi, BR=BR, Bphi=Bphi, BZ=BZ,
             V_vmec=V_vmec, V_div=V_div, relV=relV, scale=scale,
             source=str(wout_path))
    write_fluxmap_bin(str(out), rho, theta, zeta, sg, R, Z, phi, BR, Bphi, BZ, nfp, sgn)
    print(f"[qaflux] wrote {out.name}.{{npz,meta,bin}} "
          f"({n_rho}x{n_theta}x{n_zeta}, nfp={nfp}, scale={scale:g}) ALL GATES PASSED",
          flush=True)
    return str(out) + ".npz"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wout", help="path to VMEC wout_*.nc (precise-QA)")
    ap.add_argument("--out", default=str(HERE.parent.parent / "data" / "precise_QA_reactor_fluxmap"),
                    help="output stem (.npz/.meta/.bin appended)")
    ap.add_argument("--scale", type=float, default=10.10,
                    help="geometry scale (1.7/0.1683=10.10 => reactor a=1.7 m; "
                         "use 1.0 for the native unit-scale anarrima benchmark)")
    ap.add_argument("--n_rho", type=int, default=12)
    ap.add_argument("--n_theta", type=int, default=48)
    ap.add_argument("--n_zeta", type=int, default=144)
    args = ap.parse_args()
    build(args.wout, args.out, scale=args.scale,
          n_rho=args.n_rho, n_theta=args.n_theta, n_zeta=args.n_zeta)
