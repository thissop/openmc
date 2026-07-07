#!/usr/bin/env python
"""PRODUCER (vmec_venv): VMEC (vmecpp) fixed-boundary equilibrium from a QUASR device boundary ->
export the REAL metric Jacobian sqrt(g), the flux-surface geometry (R,Z,phi) and the unit field b_hat
on ONE (rho,theta,zeta) grid = the spf_fluxmap_v1 the native StellaratorSource consumes. VMEC
converges on QUASR boundaries where DESC's initial guess self-intersects (QUASR is VMEC-native).

Run in the VMEC venv:  $HOME/vmec_venv/bin/python vmec_fluxmap.py <ID>

RIGOR GATES (assert): V1b independent divergence-theorem volume vs VMEC wout volume (<1.5e-2);
V2 sqrt(g) single-signed on rho>0; the reconstructed R,Z reproduce the boundary (implicit in V1b).
b_hat: B = B^u dX/dtheta + B^v dX/dzeta with VMEC zeta = cylindrical phi.
"""
import sys
import glob
import os
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "sweep"))
from quasr_fluxmap import write_fluxmap_bin


def _solve(ID, mpol=6, ntor=6):
    import vmecpp, simsopt
    import quasr_equilibrium_field as qef
    tmpls = sorted(glob.glob(os.path.dirname(simsopt.__file__) + "/**/input.*", recursive=True))
    inp = None
    for c in tmpls:
        try:
            inp = vmecpp.VmecInput.from_file(c); break
        except Exception:
            pass
    assert inp is not None, "no VMEC input template loadable"
    nfp, Rm, Zm, meta = qef._quasr_boundary_modes(ID, flip_theta=True)
    rbc0 = np.asarray(inp.rbc); mt = rbc0.shape[0] - 1; nt = (rbc0.shape[1] - 1) // 2
    rbc = np.zeros_like(rbc0); zbs = np.zeros_like(rbc0)
    mm, nn = min(mpol, mt), min(ntor, nt)
    for (m, n), (rc, rs) in Rm.items():
        if m <= mm and abs(n) <= nn: rbc[m, n + nt] = rc
    for (m, n), (zc, zs) in Zm.items():
        if m <= mm and abs(n) <= nn: zbs[m, n + nt] = zs
    R00 = Rm[(0, 0)][0]; a = abs(Rm.get((1, 0), (0.1, 0))[0]) or 0.1
    rax = np.zeros_like(np.asarray(inp.raxis_c)); rax[0] = R00
    upd = dict(nfp=nfp, lasym=False, rbc=rbc, zbs=zbs, raxis_c=rax,
               zaxis_s=np.zeros_like(np.asarray(inp.zaxis_s)),
               ns_array=np.array([16, 51, 99]), ftol_array=np.array([1e-8, 1e-10, 1e-12]),
               niter_array=np.array([2000, 3000, 6000]), phiedge=float(np.pi * a * a),
               ncurr=1, curtor=0.0, pres_scale=0.0)
    out = vmecpp.run(inp.model_copy(update=upd), verbose=False)
    return out.wout, int(nfp)


def build(ID, n_rho=12, n_theta=48, n_zeta=144, outdir=None):
    w, nfp = _solve(ID)
    xm = np.asarray(w.xm).astype(float); xn = np.asarray(w.xn).astype(float)   # xn already includes nfp
    xmn = np.asarray(getattr(w, "xm_nyq", w.xm)).astype(float)   # Nyquist modes for gmnc/bsup*
    xnn = np.asarray(getattr(w, "xn_nyq", w.xn)).astype(float)
    ns = int(w.ns)
    s_full = np.linspace(0.0, 1.0, ns)
    s_half = (np.arange(1, ns) - 0.5) / (ns - 1)
    rmnc = np.asarray(w.rmnc); zmns = np.asarray(w.zmns); gmnc = np.asarray(w.gmnc)
    busc = np.asarray(w.bsupumnc); bvsc = np.asarray(w.bsupvmnc)

    rho = np.linspace(1.0 / n_rho, 1.0, n_rho)
    s_grid = rho ** 2                                       # VMEC s ~ toroidal flux ~ rho^2
    theta = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2 * np.pi, n_zeta, endpoint=False)
    TH, ZE = np.meshgrid(theta, zeta, indexing="ij")
    ang = xm[:, None, None] * TH[None] - xn[:, None, None] * ZE[None]
    cos, sin = np.cos(ang), np.sin(ang)
    angn = xmn[:, None, None] * TH[None] - xnn[:, None, None] * ZE[None]
    cosn = np.cos(angn)                                    # gmnc/bsup* live on the Nyquist mode set

    ip_full = lambda C: np.array([np.interp(s_grid, s_full, C[k]) for k in range(C.shape[0])])
    ip_half = lambda C: np.array([np.interp(s_grid, s_half, C[k, 1:]) for k in range(C.shape[0])])
    rc = ip_full(rmnc); zc = ip_full(zmns); gc = ip_half(gmnc)
    buc = ip_half(busc); bvc = ip_half(bvsc)

    sc = lambda Cs: np.einsum("ks,ktz->stz", Cs, cos)      # cosine synth on the R,Z mode set
    ss = lambda Cs: np.einsum("ks,ktz->stz", Cs, sin)
    scn = lambda Cs: np.einsum("ks,ktz->stz", Cs, cosn)    # cosine synth on the Nyquist mode set
    R = sc(rc); Z = ss(zc); sg = scn(gc)
    Bu = scn(buc); Bv = scn(bvc)
    dRdt = ss(-xm[:, None] * rc); dRdz = ss(xn[:, None] * rc)
    dZdt = sc(xm[:, None] * zc); dZdz = sc(-xn[:, None] * zc)
    BR = Bu * dRdt + Bv * dRdz
    BZ = Bu * dZdt + Bv * dZdz
    Bphi = R * Bv                                          # VMEC zeta = cylindrical phi
    phi = np.broadcast_to(zeta[None, None, :], (n_rho, n_theta, n_zeta)).copy()

    # ---- V2: single-signed sqrt(g) ----
    sgn = float(np.sign(np.median(sg)))
    bad = np.sign(sg) != sgn
    assert not np.any(bad) or np.all(np.abs(sg[bad]) < 1e-9 * np.abs(sg).max()), \
        f"V2 FAILED: sqrt(g) changes sign (min {sg.min():.3e}, max {sg.max():.3e})"
    sg = np.abs(sg)

    # ---- V1b: independent divergence-theorem volume from the LCFS vs VMEC volume ----
    Rb, Zb, PHb = R[-1], Z[-1], phi[-1]
    Xb = np.stack([Rb * np.cos(PHb), Rb * np.sin(PHb), Zb], axis=-1)
    dth, dze = 2 * np.pi / n_theta, 2 * np.pi / n_zeta
    Xt = np.gradient(Xb, theta, axis=0); Xz = np.gradient(Xb, zeta, axis=1)
    V_div = abs(float((Xb * np.cross(Xt, Xz)).sum() * dth * dze / 3.0))
    V_vmec = float(np.asarray(w.volume_p))
    relV = abs(V_div - V_vmec) / abs(V_vmec)
    print(f"[vmecflux] V1b: div-theorem vol {V_div:.6g} vs VMEC {V_vmec:.6g}  rel {relV:.2e}", flush=True)
    assert relV < 1.5e-2, f"V1b FAILED rel={relV:.2e}"

    out = Path(outdir) if outdir else (HERE.parent / "data")
    out.mkdir(parents=True, exist_ok=True)
    label = f"quasr{int(ID)}_vmec"
    np.savez(out / f"{label}_fluxmap.npz", rho=rho, theta=theta, zeta=zeta, nfp=nfp,
             sign_sqrtg=sgn, sqrtg=sg, R=R, Z=Z, phi=phi, BR=BR, Bphi=Bphi, BZ=BZ,
             V_vmec=V_vmec, V_div=V_div, relV=relV)
    write_fluxmap_bin(out / f"{label}_fluxmap", rho, theta, zeta, sg, R, Z, phi, BR, Bphi, BZ, nfp, sgn)
    print(f"[vmecflux] wrote {label}_fluxmap.{{npz,meta,bin}} ({n_rho}x{n_theta}x{n_zeta}, nfp={nfp}) "
          f"ALL GATES PASSED", flush=True)
    return out / f"{label}_fluxmap.npz"


if __name__ == "__main__":
    build(int(sys.argv[1]) if len(sys.argv) > 1 else 803097)
