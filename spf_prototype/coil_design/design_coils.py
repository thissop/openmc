#!/usr/bin/env python
"""Stage-2 (fixed-boundary) coil design for a ConStellaration QI boundary.

Given a ConStellaration `boundary.json` (SurfaceRZFourier r_cos/z_sin, nfp,
stellarator-symmetric), design modular filament coils on a winding surface at a
CONTROLLED minimum plasma-coil standoff, minimizing the normal-field error B.n on
the plasma boundary. Emits MAKEGRID `.coils` (cm, x100) + a coil-centroids `.npz`
+ a `_reactor_scale.txt`, matching the format the magnet_suite pipeline consumes
(see shield_opt/magnet_suite/make_coils.py::write_makegrid).

HARD GATE: reports normalized B.n and the achieved min coil-surface distance. The
caller decides pass/fail; a bad coil set is reported, never shipped silently.

Env: the fresh `simsopt` conda env (constellaration bundles simsopt==1.8.1).
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

from simsopt.geo import (SurfaceRZFourier, create_equally_spaced_curves,
                         CurveLength, CurveCurveDistance, CurveSurfaceDistance,
                         LpCurveCurvature, curves_to_vtk)
from simsopt.field import Current, coils_via_symmetries, BiotSavart
from simsopt.objectives import SquaredFlux, QuadraticPenalty
from scipy.optimize import minimize


def load_boundary(path):
    """Boundary -> simsopt SurfaceRZFourier (native scale). Accepts either a
    ConStellaration boundary.json (r_cos/z_sin) or a QUASR simsopt-serial
    (simsopt_objs; LCFS = largest-minor-radius surface, converted to RZFourier)."""
    d = json.load(open(path))
    if "simsopt_objs" in d or "r_cos" not in d:
        import simsopt
        from simsopt.geo import Surface
        acc = []
        def _fl(o):
            [_fl(x) for x in o] if isinstance(o, (list, tuple)) else acc.append(o)
        _fl(simsopt.load(str(path)))
        surfs = [o for o in acc if isinstance(o, Surface)]
        if not surfs:
            raise SystemExit(f"[load] no Surface objects in {path}")
        lcfs = max(surfs, key=lambda s: float(s.minor_radius()))
        rz = lcfs.to_RZFourier()
        return rz, int(rz.nfp)
    rc = np.array(d["r_cos"]); zs = np.array(d["z_sin"])
    nfp = int(d["n_field_periods"])
    stellsym = bool(d.get("is_stellarator_symmetric", True))
    assert stellsym, "this script assumes stellarator symmetry"
    mpol = rc.shape[0] - 1
    ntor = (rc.shape[1] - 1) // 2
    # default quadpoints (one field period grid) -- fine for minor/major radius + scaling
    s = SurfaceRZFourier(nfp=nfp, stellsym=stellsym, mpol=mpol, ntor=ntor)
    # ConStellaration and simsopt share the VMEC angle convention:
    #   R = sum rc[m,n] cos(m*theta - nfp*n*phi),  Z = sum zs[m,n] sin(...)
    for m in range(mpol + 1):
        for n in range(-ntor, ntor + 1):
            s.set_rc(m, n, rc[m, n + ntor])
            s.set_zs(m, n, zs[m, n + ntor])
    return s, nfp


def half_period_surface(s_full, nfp, nphi, ntheta):
    """Rebuild the surface on a half-period grid (stellsym) for the flux objective."""
    qp = np.linspace(0, 0.5 / nfp, nphi, endpoint=False)
    qt = np.linspace(0, 1.0, ntheta, endpoint=False)
    s = SurfaceRZFourier(nfp=nfp, stellsym=True, mpol=s_full.mpol, ntor=s_full.ntor,
                         quadpoints_phi=qp, quadpoints_theta=qt)
    s.x = s_full.x
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("boundary")
    ap.add_argument("outdir")
    ap.add_argument("--label", default=None)
    ap.add_argument("--ncoils", type=int, default=4, help="coils per half field period")
    ap.add_argument("--order", type=int, default=8, help="Fourier order per coil")
    ap.add_argument("--reactor-a", type=float, default=1.704, help="target minor radius (m)")
    ap.add_argument("--standoff", type=float, default=1.30,
                    help="TARGET min plasma-coil distance (m, reactor scale)")
    ap.add_argument("--maxiter", type=int, default=2500)
    ap.add_argument("--nphi", type=int, default=32)
    ap.add_argument("--ntheta", type=int, default=32)
    ap.add_argument("--length-weight", type=float, default=1e-3)
    ap.add_argument("--length-frac", type=float, default=1.15,
                    help="per-coil length ceiling as fraction of 2*pi*R1")
    ap.add_argument("--cc-weight", type=float, default=1.0)
    ap.add_argument("--cc-threshold", type=float, default=0.8,
                    help="coil-coil min distance target, in units of the standoff")
    ap.add_argument("--cs-weight", type=float, default=1.0)
    ap.add_argument("--curv-weight", type=float, default=1e-4)
    ap.add_argument("--curv-threshold", type=float, default=5.0,
                    help="curvature ceiling in units of 1/R1")
    ap.add_argument("--current", type=float, default=1e5)
    ap.add_argument("--init-margin", type=float, default=1.30,
                    help="init winding minor radius = a + margin*standoff (start outside barrier)")
    ap.add_argument("--sym-class", default="QI", help="tag for centroids npz (QI/QA/QH)")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    label = a.label or Path(a.boundary).stem

    # --- load boundary at NATIVE scale (R~1). Optimize here (simsopt weight recipes
    #     are tuned for this scale, and normalized B.n is scale-invariant); apply the
    #     reactor scale factor f ONLY at export + when reporting standoff in meters. ---
    s_full, nfp = load_boundary(a.boundary)
    a_native = float(s_full.minor_radius())
    R_native = float(s_full.major_radius())
    f = a.reactor_a / a_native                          # native -> reactor length factor
    standoff_nat = a.standoff / f                        # target standoff in native units
    print(f"[geom] nfp={nfp} a_native={a_native:.4f} R_native={R_native:.4f} "
          f"aspect={R_native/a_native:.2f} scale_f={f:.4f} "
          f"target_standoff={a.standoff:.3f}m -> {standoff_nat:.4f} native", flush=True)

    s = half_period_surface(s_full, nfp, a.nphi, a.ntheta)

    # --- initialize modular coils on a circular winding surface at the target standoff.
    # R1 (winding minor radius) = plasma minor radius + standoff; CurveSurfaceDistance
    # then holds the true min distance at the target. ---
    # start coils with margin OUTSIDE the target barrier so they relax onto it from
    # outside (initializing exactly at target lets the flux term pull them inside).
    R1 = a_native + a.init_margin * standoff_nat
    base_curves = create_equally_spaced_curves(a.ncoils, nfp, stellsym=True,
                                               R0=R_native, R1=R1, order=a.order)
    base_currents = [Current(a.current) for _ in range(a.ncoils)]
    base_currents[0].fix_all()                 # gauge-fix one current (net poloidal current)
    coils = coils_via_symmetries(base_curves, base_currents, nfp, True)
    bs = BiotSavart(coils)
    bs.set_points(s.gamma().reshape((-1, 3)))

    curves = [c.curve for c in coils]
    # Full-torus surface for the standoff penalty: base coils near the half-period
    # seam are closest to the plasma of the ADJACENT period, which a half-period
    # surface omits -> it under-counts the violation. Enforce distance with the FULL
    # coil set against the full torus so the reported/held min distance is honest.
    s_cs = SurfaceRZFourier(nfp=nfp, stellsym=True, mpol=s_full.mpol, ntor=s_full.ntor,
                            quadpoints_phi=np.linspace(0, 1, 8 * nfp * a.nphi // 4, endpoint=False),
                            quadpoints_theta=np.linspace(0, 1, a.ntheta, endpoint=False))
    s_cs.x = s_full.x
    Jf = SquaredFlux(s, bs, definition="local")          # normalized <(B.n)^2/|B|^2>
    Jls = [CurveLength(c) for c in base_curves]
    Jccdist = CurveCurveDistance(curves, a.cc_threshold * standoff_nat, num_basecurves=a.ncoils)
    Jcsdist = CurveSurfaceDistance(curves, s_cs, standoff_nat)
    Jcs = [LpCurveCurvature(c, 2, a.curv_threshold / R1) for c in base_curves]
    L0 = 2 * np.pi * R1 * a.length_frac        # per-coil length ceiling
    JF = (Jf
          + a.length_weight * sum(QuadraticPenalty(Jl, L0, "max") for Jl in Jls)
          + a.cc_weight * Jccdist
          + a.cs_weight * Jcsdist
          + a.curv_weight * sum(Jcs))

    def fun(x):
        JF.x = x
        return JF.J(), JF.dJ()

    print("[opt] initial: " + diagnostics(Jf, s, bs, base_curves, curves, Jcsdist, Jccdist, f), flush=True)
    res = minimize(fun, JF.x, jac=True, method="L-BFGS-B",
                   options={"maxiter": a.maxiter, "maxcor": 300}, tol=1e-15)
    print(f"[opt] final ({res.nit} its): "
          + diagnostics(Jf, s, bs, base_curves, curves, Jcsdist, Jccdist, f), flush=True)

    # --- honest metrics on the FULL boundary (all field periods), native scale ---
    s_eval = SurfaceRZFourier(nfp=nfp, stellsym=True, mpol=s_full.mpol, ntor=s_full.ntor,
                              quadpoints_phi=np.linspace(0, 1, 4 * a.nphi, endpoint=False),
                              quadpoints_theta=np.linspace(0, 1, a.ntheta, endpoint=False))
    s_eval.x = s_full.x
    bs.set_points(s_eval.gamma().reshape((-1, 3)))
    B = bs.B().reshape(s_eval.gamma().shape)
    n = s_eval.unitnormal()
    Bn = np.sum(B * n, axis=2)
    Bmod = np.linalg.norm(B, axis=2)
    relbn_mean = float(np.mean(np.abs(Bn) / Bmod))       # pointwise |B.n|/|B|, averaged
    relbn_max = float(np.max(np.abs(Bn) / Bmod))
    min_cs_nat = min_coil_surface_dist(curves, s_eval)        # full coil set, native units
    min_cs_m = float(min_cs_nat * f)                     # reactor-scale meters

    metrics = {
        "label": label, "nfp": nfp, "ncoils_per_hp": a.ncoils,
        "ncoils_total": len(coils), "order": a.order, "opt_iters": int(res.nit),
        "a_native": a_native, "R_native": R_native, "aspect": R_native / a_native,
        "scale_f": f, "a_reactor": a.reactor_a, "R_reactor": R_native * f,
        "target_standoff_m": a.standoff, "achieved_min_cs_m": min_cs_m,
        "standoff_ratio": min_cs_m / a.standoff,
        "relBn_mean": relbn_mean, "relBn_max": relbn_max,
        "coil_lengths_m": [float(Jl.J() * f) for Jl in Jls],
        "total_coil_length_m": float(sum(Jl.J() for Jl in Jls) * 2 * nfp * f),
        "min_cc_dist_m": float(Jccdist.shortest_distance() * f),
    }
    (outdir / f"{label}_coil_metrics.json").write_text(json.dumps(metrics, indent=2))
    print("[metrics] " + json.dumps(metrics), flush=True)

    # --- export MAKEGRID coils (cm) + centroids npz (magnet_suite format).
    #     coils are native; export multiplies by f (->reactor m) then x100 (->cm). ---
    write_makegrid_and_centroids(coils, base_curves, nfp, f, a_native,
                                 R_native / a_native, outdir, label, a.reactor_a,
                                 symmetry_class=a.sym_class)
    try:
        curves_to_vtk(curves, str(outdir / f"{label}_coils"))
        s_full.to_vtk(str(outdir / f"{label}_surface"))
    except Exception as e:
        print(f"[vtk] skipped: {e}", flush=True)

    ok_bn = relbn_mean < 0.02
    ok_so = min_cs_m >= 0.95 * a.standoff
    print(f"[GATE] relBn_mean={relbn_mean:.4e} (<0.02:{ok_bn})  "
          f"min_cs={min_cs_m:.3f}m target={a.standoff:.3f}m "
          f"ratio={min_cs_m/a.standoff:.3f} ({'ok' if ok_so else 'SHORT'})",
          flush=True)
    print("DESIGN_DONE", flush=True)


def diagnostics(Jf, s, bs, base_curves, curves, Jcsdist, Jccdist, f=1.0):
    bs.set_points(s.gamma().reshape((-1, 3)))
    B = bs.B().reshape(s.gamma().shape)
    n = s.unitnormal()
    Bn = np.sum(B * n, axis=2); Bmod = np.linalg.norm(B, axis=2)
    rel = np.mean(np.abs(Bn) / Bmod)
    return (f"<|Bn|/|B|>={rel:.3e}  Jf={Jf.J():.3e}  "
            f"min_cs={Jcsdist.shortest_distance()*f:.3f}m  "
            f"min_cc={Jccdist.shortest_distance()*f:.3f}m")


def min_coil_surface_dist(base_curves, surf):
    pts = surf.gamma().reshape((-1, 3))
    dmin = np.inf
    for c in base_curves:
        g = c.gamma()
        d = np.sqrt(((g[:, None, :] - pts[None, :, :]) ** 2).sum(-1)).min()
        dmin = min(dmin, d)
    return dmin


def write_makegrid_and_centroids(coils, base_curves, nfp, f, a_native, aspect,
                                 outdir, label, reactor_a, symmetry_class="QI",
                                 unit_scale=100.0):
    """All symmetry coils -> MAKEGRID (cm) + centroids npz, mirroring make_coils.py.
    Coils are at NATIVE scale; multiply by f (->reactor meters) then unit_scale (->cm)."""
    cm = f * unit_scale                                  # native -> cm at reactor scale
    coils_path = outdir / f"{label}.coils"
    with open(coils_path, "w") as fh:
        fh.write(f"periods {int(nfp)}\n")
        fh.write("begin filament\nmirror NIL\n")
        centroids = []
        currents = []
        for k, c in enumerate(coils):
            P = c.curve.gamma() * cm                    # native -> cm (reactor)
            I = float(c.current.get_value())
            for (x, y, z) in P:
                fh.write(f"{x: .10E} {y: .10E} {z: .10E} {I: .10E}\n")
            x0, y0, z0 = P[0]
            fh.write(f"{x0: .10E} {y0: .10E} {z0: .10E} {0.0: .10E} {k+1} coil{k+1}\n")
            centroids.append(P.mean(axis=0)); currents.append(I)
        fh.write("end\n")
    np.savez(outdir / f"{label}_coil_centroids.npz",
             centroids=np.array(centroids), currents=np.array(currents),
             nfp=nfp, scale_factor=f, symmetry_class=symmetry_class,
             aspect=aspect, minor_radius=a_native)
    (outdir / f"{label}_reactor_scale.txt").write_text(
        f"a_native {a_native}\nreactor_a {reactor_a}\ncoil_scale_factor {f}\n")
    print(f"[export] wrote {coils_path.name} ({len(coils)} coils, cm) + centroids + scale", flush=True)


if __name__ == "__main__":
    main()
