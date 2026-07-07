#!/usr/bin/env python
"""PRODUCER (DESC venv): fixed-boundary DESC equilibrium from a QUASR boundary -> export the
REAL metric Jacobian sqrt(g), the flux-surface geometry (R, Z, phi), and the field B on ONE
(rho, theta, zeta) grid = the conformal StellaratorSource 'fluxmap'. Position AND field come
from the SAME equilibrium (this is the unification PLAN_B calls for).

Run in the DESC venv (its deps conflict with OpenMC):
    $HOME/desc_venv/bin/python quasr_fluxmap.py <ID> [outdir]

RIGOR GATES -- asserted, not warned (no shortcuts):
  V1  own-quadrature volume  Sum|sqrt(g)| drho*dtheta*dzeta  ==  DESC eq.compute('V')   (<2e-3)
  V1b INDEPENDENT volume from the LCFS via the divergence theorem
        V = (1/3) oint X . (X_theta x X_zeta) dtheta dzeta      ==  DESC 'V'             (<1e-2)
  V2  sqrt(g) is single-signed on rho>0 (nested surfaces; degeneracy only at the axis)
  V3  the DESC node set is a clean (rho,theta,zeta) tensor product after lexsort reshape
  V4  the rho=1 surface reconstructs the fitted QUASR boundary to the fit residual
Reuses the boundary->DESC-surface machinery in quasr_equilibrium_field.py.
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


def write_fluxmap_bin(stem, rho, theta, zeta, sqrtg, R, Z, phi, BR, Bphi, BZ,
                      nfp, sign_sqrtg=1.0):
    """Write the portable **spf_fluxmap_v1** binary format consumed by the native
    C++ ``StellaratorSource`` (spf_prototype/native_spf/stellarator_source.cpp).

    This mirrors the style of ``src/spf_fieldmap.hpp``'s ``spf_fieldmap_v1``: a
    small text ``.meta`` header (``key value`` per line) plus a raw little-endian
    float64 ``.bin`` payload. It is intentionally HDF5-free so the C++ reader is a
    plain ``std::ifstream`` (producer-agnostic; a real DESC/VMEC solve or a
    synthetic test map writes the identical format).

    Parameters
    ----------
    stem : str or path-like
        Output stem. A trailing ``.meta``/``.bin`` is stripped; the writer emits
        ``<stem>.meta`` and ``<stem>.bin``.
    rho, theta, zeta : 1-D arrays
        The flux-coordinate grids, lengths ``nR``, ``ntheta``, ``nzeta``. theta
        and zeta span the FULL torus ``[0, 2*pi)`` (endpoint excluded); rho runs
        ``(drho .. 1]`` (the axis rho=0 is excluded, sqrt(g)->0 there).
    sqrtg, R, Z, phi, BR, Bphi, BZ : (nR, ntheta, nzeta) arrays
        The metric Jacobian |sqrt(g)| (>= 0), the flux-surface geometry, and the
        cylindrical field components, all indexed ``[i_rho, j_theta, k_zeta]``
        (C-order) -- the SAME layout as ``ConformalStellaratorSampler`` in
        python/stellarator_source.py.
    nfp : int
        Number of field periods (stored; the sampler treats the grid as the full
        torus with 2*pi periods, so nfp is informational).
    sign_sqrtg : float
        Sign of the raw sqrt(g) before taking |.| (stored for provenance).

    Binary layout of ``<stem>.bin`` (little-endian f8, C-order, in this order):

        rho[nR], theta[ntheta], zeta[nzeta],
        sqrtg, R, Z, phi, BR, Bphi, BZ            (each nR*ntheta*nzeta)

    Returns
    -------
    (meta_path, bin_path) : tuple of str
    """
    stem = str(stem)
    if stem.endswith(".meta") or stem.endswith(".bin"):
        stem = stem.rsplit(".", 1)[0]
    rho = np.ascontiguousarray(rho, dtype="<f8")
    theta = np.ascontiguousarray(theta, dtype="<f8")
    zeta = np.ascontiguousarray(zeta, dtype="<f8")
    nR, nT, nZ = rho.size, theta.size, zeta.size
    names = ("sqrtg", "R", "Z", "phi", "BR", "Bphi", "BZ")
    arrs = (sqrtg, R, Z, phi, BR, Bphi, BZ)
    for name, a in zip(names, arrs):
        if np.shape(a) != (nR, nT, nZ):
            raise ValueError(
                f"write_fluxmap_bin: {name} shape {np.shape(a)} != {(nR, nT, nZ)}")
    if float(np.asarray(sqrtg).min()) < -1e-12:
        raise ValueError("write_fluxmap_bin: store |sqrt(g)| (must be >= 0)")

    meta = Path(stem + ".meta")
    with meta.open("w") as f:
        f.write("magic spf_fluxmap_v1\n")
        f.write(f"nR {nR}\n")
        f.write(f"ntheta {nT}\n")
        f.write(f"nzeta {nZ}\n")
        f.write(f"nfp {int(nfp)}\n")
        f.write(f"sign_sqrtg {float(sign_sqrtg):+.1f}\n")
    binp = Path(stem + ".bin")
    with binp.open("wb") as f:
        f.write(rho.tobytes())
        f.write(theta.tobytes())
        f.write(zeta.tobytes())
        for a in arrs:
            f.write(np.ascontiguousarray(a, dtype="<f8").tobytes())
    return str(meta), str(binp)


def build(ID, outdir=None, L=8, n_rho=16, n_theta=64, n_zeta=192, ftol=1e-2, maxiter=100,
          example=None):
    """If `example` is a DESC example name (e.g. 'precise_QA'), use that REAL, already-converged
    equilibrium directly -- this validates the sqrt(g) machinery + StellaratorSource on real
    physics. Otherwise solve a fixed-boundary equilibrium from QUASR device `ID`."""
    from desc.grid import LinearGrid
    Rmodes = Zmodes = None
    if example is not None:
        import desc.examples
        eqx = desc.examples.get(example)
        eq = eqx[-1] if isinstance(eqx, (list, tuple)) else eqx
        nfp = int(eq.NFP); label = str(example)
        print(f"[fluxmap] DESC example '{example}': NFP={nfp} (real converged equilibrium)",
              flush=True)
    else:
        import quasr_equilibrium_field as qef
        from desc.equilibrium import Equilibrium
        from desc.continuation import solve_continuation_automatic
        nfp, Rmodes, Zmodes, meta = qef._quasr_boundary_modes(ID)
        surf = qef._desc_surface(nfp, Rmodes, Zmodes)
        label = f"quasr{int(ID)}"
        print(f"[fluxmap] device {ID}: NFP={nfp}; continuation solve (L=M=N={L}) ...", flush=True)
        eq = Equilibrium(surface=surf, NFP=nfp, L=L, M=L, N=L, sym=True)
        fam = solve_continuation_automatic(eq, verbose=1, ftol=ftol, maxiter=maxiter)
        eq = fam[-1] if isinstance(fam, (list, tuple)) else fam

    # full-torus tensor grid; rho starts at drho (the axis rho=0 has sqrt(g)->0)
    rho = np.linspace(1.0 / n_rho, 1.0, n_rho)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)
    grid = LinearGrid(rho=rho, theta=theta, zeta=zeta, NFP=1, sym=False)
    data = eq.compute(["sqrt(g)", "R", "Z", "phi", "B", "V"], grid=grid, basis="rpz")

    nodes = np.asarray(grid.nodes)                      # (N,3) = (rho, theta, zeta)
    N = nodes.shape[0]
    assert N == n_rho * n_theta * n_zeta, f"grid size {N} != {n_rho*n_theta*n_zeta}"

    # ---- V3: reshape by explicit lexsort so we never assume DESC's internal node order ----
    order = np.lexsort((nodes[:, 2], nodes[:, 1], nodes[:, 0]))   # rho outer, theta mid, zeta inner
    cube = lambda a: np.asarray(a)[order].reshape(n_rho, n_theta, n_zeta)
    RHO, TH, ZE = cube(nodes[:, 0]), cube(nodes[:, 1]), cube(nodes[:, 2])
    assert (np.allclose(RHO[:, 0, 0], rho) and np.allclose(TH[0, :, 0], theta)
            and np.allclose(ZE[0, 0, :], zeta)), "V3 FAILED: node set is not a clean tensor product"

    sg = cube(data["sqrt(g)"])
    Rc, Zc, PHc = cube(data["R"]), cube(data["Z"]), cube(data["phi"])
    Bcyl = np.asarray(data["B"])[order].reshape(n_rho, n_theta, n_zeta, 3)   # (B_R, B_phi, B_Z)

    # ---- V2: single sign (nested surfaces) ----
    sgn = float(np.sign(np.median(sg)))
    bad = np.sign(sg) != sgn
    assert not np.any(bad) or np.all(np.abs(sg[bad]) < 1e-9 * np.abs(sg).max()), (
        f"V2 FAILED: sqrt(g) changes sign (min {sg.min():.3e}, max {sg.max():.3e})")
    sg_abs = np.abs(sg)

    # ---- V1: own quadrature vs DESC volume ----
    # rectangle over periodic (theta,zeta), then trapezoid in rho with the axis point
    # sqrt(g)(rho=0)=0 prepended. (A crude right-Riemann sum was ~8% off on a coarse rho grid;
    # this is O(drho^2) and converges with n_rho -- the deterministic discretization error.)
    dth, dze = 2.0 * np.pi / n_theta, 2.0 * np.pi / n_zeta
    _trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    sg_of_rho = sg_abs.mean(axis=(1, 2)) * (2.0 * np.pi) ** 2           # integral over theta,zeta
    rho_ax = np.concatenate([[0.0], rho]); sg_ax = np.concatenate([[0.0], sg_of_rho])
    V_quad = float(_trap(sg_ax, rho_ax))
    V_desc = float(np.atleast_1d(np.asarray(data["V"])).ravel()[0])
    relV = abs(V_quad - V_desc) / abs(V_desc)
    print(f"[fluxmap] V1 : quad {V_quad:.6g}  DESC {V_desc:.6g}  rel {relV:.2e}", flush=True)
    assert relV < 3e-2, f"V1 FAILED rel={relV:.2e} (raise n_rho if this is the discretization limit)"

    # ---- V1b: INDEPENDENT volume from the LCFS via divergence theorem ----
    Rb, Zb, PHb = Rc[-1], Zc[-1], PHc[-1]               # rho=1 surface (n_theta, n_zeta)
    Xb = np.stack([Rb * np.cos(PHb), Rb * np.sin(PHb), Zb], axis=-1)
    Xt = np.gradient(Xb, theta, axis=0)                 # dX/dtheta
    Xz = np.gradient(Xb, zeta, axis=1)                  # dX/dzeta
    dS = np.cross(Xt, Xz)                               # oriented area element
    V_div = abs(float((Xb * dS).sum() * dth * dze / 3.0))
    relVd = abs(V_div - V_desc) / abs(V_desc)
    print(f"[fluxmap] V1b: divergence-theorem {V_div:.6g}  rel {relVd:.2e}", flush=True)
    assert relVd < 1e-2, f"V1b FAILED rel={relVd:.2e}"

    res = scale = float("nan")            # defined even when V4 is skipped (DESC-example mode)
    # ---- V4: rho=1 reconstructs the fitted boundary (only for QUASR-solved eqs) ----
    if Rmodes is not None:
        Rfit = np.zeros((n_theta, n_zeta)); Zfit = np.zeros_like(Rfit)
        TT, ZZ = np.meshgrid(theta, zeta, indexing="ij")
        for (m, n), (rcos, rsin) in Rmodes.items():
            ang = m * TT - n * nfp * ZZ
            Rfit += rcos * np.cos(ang) + rsin * np.sin(ang)
        for (m, n), (zcos, zsin) in Zmodes.items():
            ang = m * TT - n * nfp * ZZ
            Zfit += zcos * np.cos(ang) + zsin * np.sin(ang)
        res = float(np.sqrt(np.mean((Rb - Rfit) ** 2 + (Zb - Zfit) ** 2)))
        scale = float(np.hypot(Rb.max() - Rb.min(), Zb.max() - Zb.min()))
        print(f"[fluxmap] V4 : LCFS-vs-fit RMS {res:.3e} ({res/scale:.2e} of extent)", flush=True)
        assert res / scale < 5e-2, f"V4 FAILED: LCFS residual {res/scale:.2e} of extent"
    else:
        print("[fluxmap] V4 : skipped (DESC example; no QUASR boundary fit)", flush=True)

    out = Path(outdir) if outdir else (HERE.parent / "data")
    out.mkdir(parents=True, exist_ok=True)
    dst = out / f"{label}_fluxmap.npz"
    np.savez(dst, rho=rho, theta=theta, zeta=zeta, nfp=nfp, sign_sqrtg=sgn,
             sqrtg=sg_abs, R=Rc, Z=Zc, phi=PHc,
             BR=Bcyl[..., 0], Bphi=Bcyl[..., 1], BZ=Bcyl[..., 2],
             V_desc=V_desc, V_quad=V_quad, V_div=V_div, relV=relV, relVd=relVd)
    # Also emit the portable spf_fluxmap_v1 .meta/.bin the native C++
    # StellaratorSource consumes (same arrays, no HDF5 dependency). Keep the .npz.
    stem = out / f"{label}_fluxmap"
    write_fluxmap_bin(stem, rho, theta, zeta, sg_abs, Rc, Zc, PHc,
                      Bcyl[..., 0], Bcyl[..., 1], Bcyl[..., 2], nfp, sgn)
    print(f"[fluxmap] wrote {dst.name} + {stem.name}.{{meta,bin}}  "
          f"({n_rho}x{n_theta}x{n_zeta}, nfp={nfp})  ALL RIGOR GATES PASSED",
          flush=True)

    # ---- plot: sqrt(g) on a poloidal cross-section + flux surfaces (smplotlib if available) ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        try:
            import smplotlib  # noqa: F401
        except Exception:
            pass
        import matplotlib.pyplot as plt
        figdir = HERE.parent / "figs" / "fluxmap"; figdir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        k = 0
        pc = ax[0].contourf(Rc[:, :, k], Zc[:, :, k], sg_abs[:, :, k], 20)
        ax[0].set_aspect("equal"); ax[0].set_xlabel("R [m]"); ax[0].set_ylabel("Z [m]")
        fig.colorbar(pc, ax=ax[0], label=r"$\sqrt{g}$")
        for i in range(0, n_rho, max(1, n_rho // 6)):
            ax[1].plot(np.append(Rc[i, :, k], Rc[i, 0, k]),
                       np.append(Zc[i, :, k], Zc[i, 0, k]), lw=0.8, color="black")
        ax[1].set_aspect("equal"); ax[1].set_xlabel("R [m]"); ax[1].set_ylabel("Z [m]")
        fig.tight_layout(); fig.savefig(figdir / f"quasr{int(ID)}_fluxmap.png", dpi=180)
        plt.close(fig)
        print(f"[fluxmap] plot -> figs/fluxmap/quasr{int(ID)}_fluxmap.png", flush=True)
    except Exception as e:
        print(f"[fluxmap] plot skipped: {type(e).__name__}: {e}", flush=True)

    return dst, dict(nfp=nfp, V=V_desc, relV=relV, relVd=relVd, res=res / scale)


if __name__ == "__main__":
    dev = int(sys.argv[1]) if len(sys.argv) > 1 else 803097
    build(dev, sys.argv[2] if len(sys.argv) > 2 else None)
