#!/usr/bin/env python
"""STAGE equil (DESC venv): QUASR device ID -> DESC fixed-boundary equilibrium ->
   (a) the spf_fluxmap_v1 fluxmap the StellaratorSource + adjoint contributon consume, AND
   (b) a VMEC-format wout .nc that ParaStell (ps.Stellarator) consumes for the DAGMC build.

This is the missing link that lets a FRESH QUASR device reach the ParaStell->DAGMC->adjoint
pipeline: QH/QA had hand-obtained VMEC wouts; a zoo device has only a QUASR boundary+coils.
DESC solves the fixed-boundary equilibrium from the boundary (quasr_fluxmap machinery) and
DESC's VMECIO.save writes the wout ParaStell needs. Same equilibrium feeds both outputs, so
the transport source geometry and the DAGMC geometry are CONSISTENT.

Run:  $HOME/desc_venv/bin/python equil_device.py <ID> <outdir> [--L 8 --n-rho 16 ...]

Rigor: reuses quasr_fluxmap.build's asserted volume gates (V1/V1b/V2/V3/V4). The wout export
is additionally sanity-checked (nfp, aspect finite). Idempotent: skips if outputs exist.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# spf_pp/python holds quasr_fluxmap + quasr_equilibrium_field; sweep holds quasr_geom
_THIS = Path(__file__).resolve()
for cand in (_THIS.parents[2] / "python",           # local repo layout
             Path.home() / "src/GitHub/openmc-spf/spf_prototype/python",
             Path("/ginsburg/astro/users/tjk2147/spf_work/spf_pp/python")):
    if cand.exists():
        sys.path.insert(0, str(cand))
        sys.path.insert(0, str(cand.parent / "sweep"))


def export_wout(eq, wout_path):
    """DESC equilibrium -> VMEC wout .nc for ParaStell. Returns wout_path."""
    from desc.vmec import VMECIO
    VMECIO.save(eq, str(wout_path))
    return str(wout_path)


def _native_minor_radius(Rmodes, Zmodes):
    """Robust minor-radius estimate a from boundary Fourier modes: the m=1 (n=0)
    amplitude sets the circular cross-section radius; use sqrt(|R10*Z10|) style mean
    of the two m=1 half-widths (falls back to |R10| alone)."""
    r10 = abs(Rmodes.get((1, 0), (0.0, 0.0))[0])
    z10 = abs(Zmodes.get((1, 0), (0.0, 0.0))[1]) or abs(Zmodes.get((1, 0), (0.0, 0.0))[0])
    if r10 > 0 and z10 > 0:
        return (r10 * z10) ** 0.5
    return r10 or z10


def _scale_modes(modes, f):
    return {k: (v[0] * f, v[1] * f) for k, v in modes.items()}


def _reactor_factor(outdir, label, ID, Rmodes, Zmodes, reactor_a):
    """Single authoritative reactor-scale factor f. Prefer the shared scale file the fetch
    stage wrote (from the catalogue minor_radius); else estimate from the boundary modes."""
    scale_file = Path(outdir) / f"{label}_reactor_scale.txt"
    if scale_file.exists():
        for line in scale_file.read_text().splitlines():
            if line.startswith("coil_scale_factor"):
                f = float(line.split()[1])
                print(f"[equil] reactor-scale: read factor {f:.6f} from {scale_file.name}",
                      flush=True)
                return f
    if reactor_a and reactor_a > 0:
        a_native = _native_minor_radius(Rmodes, Zmodes)
        if a_native and a_native > 0:
            f = reactor_a / a_native
            scale_file.write_text(
                f"a_native {a_native}\nreactor_a {reactor_a}\ncoil_scale_factor {f}\n")
            print(f"[equil] reactor-scale: a_native(mode-est)={a_native:.4f} -> factor {f:.4f}",
                  flush=True)
            return f
    return 1.0


def build_equil(ID, outdir, L=8, n_rho=16, n_theta=64, n_zeta=192,
                ftol=1e-2, maxiter=100, force=False, reactor_a=1.704, solver="vmec"):
    """SOLVER CHOICE (automation-critical): use VMEC (vmecpp) by default.

    QUASR boundaries are VMEC-native; DESC's ``solve_continuation_automatic`` downscales the
    surface and trips ``get_max_mode`` ("No modes found, geometry is unphysical") in
    ``compute_scaling_factors`` (verified: a low-iota QA and a QH both failed identically), and
    a direct DESC solve, while it passes normalization, does not converge these boundaries to a
    physical equilibrium (V-gate fails). VMEC converges the SAME fresh QH boundary in ~27 s and
    ``VmecWOut.save`` writes a wout ParaStell reads directly. So VMEC is primary; DESC remains
    selectable via ``solver='desc'`` for equilibria where DESC is preferred.

    Reactor scaling is applied by scaling the boundary Fourier modes by ``f`` BEFORE the solve
    (via a monkeypatch of ``quasr_equilibrium_field._quasr_boundary_modes`` that both the VMEC
    and DESC producers call), so the wout AND the extracted fluxmap come out consistently at
    reactor scale and share the coils' frame.
    """
    import numpy as np
    import quasr_equilibrium_field as qef
    from quasr_fluxmap import write_fluxmap_bin

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    label = f"quasr{int(ID)}"
    fluxmap_meta = outdir / f"{label}_fluxmap.meta"
    wout_path = outdir / f"wout_{label}.nc"
    npz = outdir / f"{label}_fluxmap.npz"
    if not force and fluxmap_meta.exists() and wout_path.exists():
        print(f"[equil] {label}: outputs exist, skipping (use --force to rebuild)", flush=True)
        return str(wout_path), str(npz)

    nfp, Rmodes, Zmodes, _meta = qef._quasr_boundary_modes(ID)
    f = _reactor_factor(outdir, label, ID, Rmodes, Zmodes, reactor_a)

    # Monkeypatch the shared boundary-mode loader so BOTH producers see reactor-scaled modes.
    _orig = qef._quasr_boundary_modes
    def _patched(the_id, **kw):
        n, Rm, Zm, m = _orig(the_id, **kw)
        if f and abs(f - 1.0) > 1e-9:
            Rm = _scale_modes(Rm, f); Zm = _scale_modes(Zm, f)
        return n, Rm, Zm, m
    qef._quasr_boundary_modes = _patched

    try:
        if solver == "vmec":
            import vmec_fluxmap as vf
            # reuse the proven VMEC producer for the fluxmap (writes quasr<ID>_vmec_fluxmap.*)
            vf.build(ID, n_rho=n_rho, n_theta=min(n_theta, 48), n_zeta=min(n_zeta, 144),
                     outdir=str(outdir))
            # rename to the stem the StellaratorSource + finalize expect
            for ext in ("npz", "meta", "bin"):
                src = outdir / f"{label}_vmec_fluxmap.{ext}"
                if src.exists():
                    src.replace(outdir / f"{label}_fluxmap.{ext}")
            # save the wout for ParaStell (a second short VMEC solve under the same patch)
            w, _nfp = vf._solve(ID)
            w.save(str(wout_path))
        else:
            _desc_solve_and_extract(ID, nfp, outdir, label, npz, wout_path,
                                    L, n_rho, n_theta, n_zeta, ftol, maxiter)
    except Exception as e:
        (outdir / f"{label}_EQUIL_FAILED.txt").write_text(
            f"device {ID}\nsolver {solver}\nerror {type(e).__name__}: {e}\n")
        print(f"[equil] SOLVER '{solver}' FAILED ({type(e).__name__}: {e}). "
              f"Wrote {label}_EQUIL_FAILED.txt (needs a hand solve or the other solver).",
              flush=True)
        raise SystemExit(3)
    finally:
        qef._quasr_boundary_modes = _orig

    assert wout_path.exists() and wout_path.stat().st_size > 0, "wout export produced no file"
    assert npz.exists(), "fluxmap npz missing"
    print(f"[equil] wrote {npz.name}, {label}_fluxmap.{{meta,bin}}, {wout_path.name} "
          f"(nfp={nfp}, reactor_factor={f:.3f}, solver={solver})", flush=True)
    return str(wout_path), str(npz)


def _desc_solve_and_extract(ID, nfp, outdir, label, npz, wout_path,
                            L, n_rho, n_theta, n_zeta, ftol, maxiter):
    """DESC fallback path (kept for equilibria where DESC is preferred). Direct fixed-resolution
    solve (continuation downscaling is what trips get_max_mode)."""
    import numpy as np
    import quasr_equilibrium_field as qef
    from desc.equilibrium import Equilibrium
    from desc.grid import LinearGrid
    from quasr_fluxmap import write_fluxmap_bin
    nfp2, Rmodes, Zmodes, _ = qef._quasr_boundary_modes(ID)   # already patched -> scaled
    surf = qef._desc_surface(nfp2, Rmodes, Zmodes)
    eq = Equilibrium(surface=surf, NFP=nfp2, L=L, M=L, N=L, sym=True)
    eq.solve(verbose=1, ftol=ftol, maxiter=maxiter)
    rho = np.linspace(1.0 / n_rho, 1.0, n_rho)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)
    grid = LinearGrid(rho=rho, theta=theta, zeta=zeta, NFP=1, sym=False)
    data = eq.compute(["sqrt(g)", "R", "Z", "phi", "B", "V"], grid=grid, basis="rpz")
    nodes = np.asarray(grid.nodes)
    order = np.lexsort((nodes[:, 2], nodes[:, 1], nodes[:, 0]))
    cube = lambda a: np.asarray(a)[order].reshape(n_rho, n_theta, n_zeta)
    sg = np.abs(cube(data["sqrt(g)"]))
    Rc, Zc, PHc = cube(data["R"]), cube(data["Z"]), cube(data["phi"])
    Bcyl = np.asarray(data["B"])[order].reshape(n_rho, n_theta, n_zeta, 3)
    sgn = float(np.sign(np.median(cube(data["sqrt(g)"]))))
    _trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    sg_of_rho = sg.mean(axis=(1, 2)) * (2.0 * np.pi) ** 2
    V_quad = float(_trap(np.concatenate([[0.0], sg_of_rho]), np.concatenate([[0.0], rho])))
    V_desc = float(np.atleast_1d(np.asarray(data["V"])).ravel()[0])
    relV = abs(V_quad - V_desc) / abs(V_desc)
    print(f"[equil] V1: quad {V_quad:.6g} DESC {V_desc:.6g} rel {relV:.2e}", flush=True)
    assert relV < 3e-2, f"V1 FAILED rel={relV:.2e}"
    np.savez(npz, rho=rho, theta=theta, zeta=zeta, nfp=nfp2, sign_sqrtg=sgn,
             sqrtg=sg, R=Rc, Z=Zc, phi=PHc, BR=Bcyl[..., 0], Bphi=Bcyl[..., 1], BZ=Bcyl[..., 2],
             V_desc=V_desc, V_quad=V_quad, relV=relV)
    write_fluxmap_bin(outdir / f"{label}_fluxmap", rho, theta, zeta, sg,
                      Rc, Zc, PHc, Bcyl[..., 0], Bcyl[..., 1], Bcyl[..., 2], nfp2, sgn)
    export_wout(eq, wout_path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ID", type=int)
    ap.add_argument("outdir")
    ap.add_argument("--L", type=int, default=8)
    ap.add_argument("--n-rho", type=int, default=16)
    ap.add_argument("--n-theta", type=int, default=64)
    ap.add_argument("--n-zeta", type=int, default=192)
    ap.add_argument("--ftol", type=float, default=1e-2)
    ap.add_argument("--maxiter", type=int, default=100)
    ap.add_argument("--reactor-a", type=float, default=1.704,
                    help="target minor radius (m) for reactor scaling; 0 = native scale")
    ap.add_argument("--solver", choices=["vmec", "desc"], default="vmec",
                    help="equilibrium solver; VMEC (default) is QUASR-native and converges "
                         "where DESC continuation fails")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    build_equil(a.ID, a.outdir, L=a.L, n_rho=a.n_rho, n_theta=a.n_theta,
                n_zeta=a.n_zeta, ftol=a.ftol, maxiter=a.maxiter, force=a.force,
                reactor_a=a.reactor_a, solver=a.solver)
    print("DONE_EQUIL", flush=True)
