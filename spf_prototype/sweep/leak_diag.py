#!/usr/bin/env python
"""Native (no-transport) DAGMC-leak diagnosis for the conformal-shell sweep.

trimesh is the ground truth. For a device id and a surface-mesh resolution
(ntheta x nph), we rebuild the conformal shell STLs (stellarator_geometry.
build_layers, with the sweep's scale-retry), then reconstruct the EXACT boundary
surfaces the DAGMC .h5m is assembled from (dagmc_writer._read_stl_tris +
_dedup_surface), and trimesh-check every boundary:

  is_watertight, is_winding_consistent, #naked(boundary) edges, #degenerate faces,
  and whether adjacent boundaries INTERSECT (shell crossing / inboard pinch).

This isolates the leak fork:
  - STL/boundary itself non-watertight (naked edges > 0) -> build-resolution problem
  - boundaries watertight but adjacent shells self-intersect -> concave-inboard
    over-offset (intrinsic unless resolution or scale changes the fold)

Private stem/dir so it never touches data/quasr<id>_surface.npz that the running
VM sweep reads. Writes nothing into sweep_out*/.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
import sweep, quasr_loader as ql, stellarator_geometry as sg, dagmc_writer as dw

DATADIR = HERE.parent / "data"
UNIT_CM = 100.0
SCRATCH = HERE / "leak_diag_out"


def stage_surface(dev, stem, nt, nph):
    th = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z = dev.device.RZ(TH, PH, 1.0)
    Rc, Zc = R * UNIT_CM, Z * UNIT_CM
    np.savez(DATADIR / f"{stem}_surface.npz", R=Rc, Z=Zc,
             nfp=np.int64(dev.meta["nfp"]), theta=th, phi_axis=ph,
             phi=np.broadcast_to(ph[None, :], Rc.shape).copy())
    return 0.5 * (R.max() - R.min())


def _boundaries_from_stls(geom_dir):
    """Reproduce dagmc_writer.build_from_stls boundary extraction (dedup path)."""
    man = json.loads((Path(geom_dir) / "manifest.json").read_text())
    layers = [L for L in man["layers"] if L["name"] not in {"sol"}]
    nb = len(layers) + 1
    boundaries = [None] * nb
    names = [None] * nb
    for k, L in enumerate(layers):
        tris = dw._read_stl_tris(Path(geom_dir) / L["stl"])
        m = len(tris) // 2
        outer, inner = tris[:m], tris[m:]
        if k == 0:
            boundaries[0] = dw._orient_outward(dw._dedup_surface(inner), sg)
            names[0] = f"{L['name']}_inner"
        boundaries[k + 1] = dw._orient_outward(dw._dedup_surface(outer), sg)
        names[k + 1] = f"{L['name']}_outer"
    return boundaries, names


def _tm_report(V, T):
    import trimesh
    mesh = trimesh.Trimesh(vertices=np.asarray(V), faces=np.asarray(T), process=False)
    edges = mesh.edges_sorted
    # naked edges = edges used by exactly one face
    from trimesh.grouping import group_rows
    groups = group_rows(edges, require_count=None)
    counts = np.array([len(g) for g in groups])
    naked = int((counts == 1).sum())
    nonman = int((counts > 2).sum())
    # degenerate faces (zero area)
    areas = mesh.area_faces
    degen = int((areas <= 1e-9).sum())
    return dict(
        watertight=bool(mesh.is_watertight),
        winding_consistent=bool(mesh.is_winding_consistent),
        naked_edges=naked, nonmanifold_edges=nonman,
        degenerate_faces=degen, n_faces=int(len(T)), n_verts=int(len(V)),
        euler=int(mesh.euler_number), volume=float(mesh.volume))


def _adjacent_intersections(boundaries, names):
    """Do adjacent boundaries cross? Cheap surrogate: fraction of INNER-boundary
    vertices that lie OUTSIDE the OUTER boundary (a shell that pinches through the
    next one puts inner verts on the wrong side). Uses trimesh contains()."""
    import trimesh
    out = []
    for i in range(len(boundaries) - 1):
        Vi, Ti = boundaries[i]
        Vo, To = boundaries[i + 1]
        outer = trimesh.Trimesh(vertices=np.asarray(Vo), faces=np.asarray(To),
                                process=False)
        if not outer.is_watertight:
            out.append(dict(pair=f"{names[i]}->{names[i+1]}",
                            note="outer not watertight; contains() unreliable"))
            continue
        inside = outer.contains(np.asarray(Vi))
        frac_out = float((~inside).mean())
        out.append(dict(pair=f"{names[i]}->{names[i+1]}",
                        inner_verts_outside_outer=frac_out))
    return out


def diagnose(cid, nt, nph, scale_override=None, check_intersect=True):
    stem = f"quasr{cid}_ld"
    dev = ql.load_device(cid, n_samples=64)
    a_native = stage_surface(dev, stem, nt, nph)
    scale = scale_override or (sweep.TARGET_A_CM / (a_native * UNIT_CM))
    layers = sweep.BLANKETS["baseline"]
    geom_dir = SCRATCH / f"{stem}_{nt}x{nph}_geom"
    geom_dir.mkdir(parents=True, exist_ok=True)

    m_geo, used_scale, folds = None, scale, []
    for _ in range(4):
        m = sg.build_layers(stem, layers=layers, scale=used_scale, outdir=geom_dir)
        bad = [L["name"] for L in m["layers"]
               if not (L["watertight"] and L["simple_cross_section"])]
        if bad:
            folds.append([round(used_scale, 3), bad]); used_scale *= 1.4; continue
        m_geo = m; break
    rec = dict(id=cid, res=[nt, nph], nfp=int(dev.meta["nfp"]),
               a_native_cm=round(a_native * UNIT_CM, 3),
               scale=round(used_scale, 4), scale_retry_folds=folds)
    if m_geo is None:
        rec["verdict"] = "geometry_folded_all_scales"
        return rec

    boundaries, names = _boundaries_from_stls(geom_dir)
    per = []
    any_naked = False
    for nm, (V, T) in zip(names, boundaries):
        r = _tm_report(V, T); r["boundary"] = nm
        any_naked |= (r["naked_edges"] > 0 or not r["watertight"])
        per.append(r)
    rec["boundaries"] = per
    rec["all_boundaries_watertight"] = not any_naked

    if check_intersect:
        rec["adjacent_intersections"] = _adjacent_intersections(boundaries, names)
        worst = max((d.get("inner_verts_outside_outer", 0.0)
                     for d in rec["adjacent_intersections"]), default=0.0)
        rec["worst_shell_crossing_frac"] = round(worst, 5)

    if not rec["all_boundaries_watertight"]:
        rec["verdict"] = "STL_boundary_nonwatertight"
    elif rec.get("worst_shell_crossing_frac", 0.0) > 1e-4:
        rec["verdict"] = "boundaries_wt_but_shells_intersect"
    else:
        rec["verdict"] = "clean"
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", help="comma-separated device ids")
    ap.add_argument("--res", default="64x96,128x192",
                    help="comma list of NTxNPH resolutions to test")
    ap.add_argument("--no-intersect", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ids = [int(x) for x in a.ids.split(",")]
    resolutions = [tuple(int(v) for v in r.split("x")) for r in a.res.split(",")]
    allrecs = []
    for cid in ids:
        for (nt, nph) in resolutions:
            try:
                rec = diagnose(cid, nt, nph, check_intersect=not a.no_intersect)
            except Exception as e:
                import traceback
                rec = dict(id=cid, res=[nt, nph], verdict="EXCEPTION",
                           error=f"{type(e).__name__}: {e}",
                           tb=traceback.format_exc().splitlines()[-3:])
            allrecs.append(rec)
            print(json.dumps({k: rec[k] for k in rec
                              if k not in ("boundaries", "adjacent_intersections")}))
    if a.out:
        Path(a.out).write_text(json.dumps(allrecs, indent=2))
        print(f"[wrote] {a.out}")
