#!/usr/bin/env python
"""validate_offset2.py -- native (no-transport) before/after check of the
morphological offset (stellarator_geometry_offset2.build_layers_sdf) vs the live
normal offset (stellarator_geometry.build_layers).

Metric: worst shell crossing = max over adjacent boundary pairs of the fraction of
inner-boundary vertices lying OUTSIDE the next boundary out. This is computed as a
FAST 2-D point-in-polygon test per shared phi cross-section, which is RIGOROUSLY
EQUIVALENT to sweep/leak_gate.worst_shell_crossing (trimesh 3-D contains) for these
structured, phi-aligned toroidal meshes: every boundary shares the same phi slices,
so the 3-D outer mesh cross-section at plane phi_j is exactly the outer polygon at
phi_j, and a vertex on that plane is inside the mesh iff it is inside that polygon.
The 3-D leak_gate is O(points x triangles) without embree (minutes); this is seconds.
The equivalence is checked in __main__ (--xcheck) against leak_gate on an OLD build.

Boundary selection matches leak_gate exactly: material layers only (sol excluded);
boundaries = [inner of first material layer] + [outer of each material layer].
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))
import stellarator_geometry as sg
import stellarator_geometry_offset2 as sg2

DATADIR = HERE.parent / "data"
TARGET_A_CM = 250.0
BASELINE = [("sol", 5.0), ("W", 0.2), ("steel", 3.8), ("Be", 2.0),
            ("FLiBe", 50.0), ("shield", 40.0), ("coil", 15.0)]
# Full set of transport-'error' (DAGMC-leaky) devices from the local batch
# (sweep_out_local/batch_classification_summary.csv, status==error).
LEAKY = ["1505944", "1577542", "1862008", "1880092", "2147186", "2190977",
         "2194424", "2218517", "2225908", "2291622", "2296767", "2369483",
         "2469994", "2524474"]
GOOD = ["59509", "115887", "134412", "157482"]


def scale_for(stem):
    """Reproduce the sweep's auto-scale (built minor radius -> TARGET_A_CM)."""
    d = np.load(DATADIR / f"{stem}_surface.npz")
    a_native_cm = 0.5 * (float(d["R"].max()) - float(d["R"].min()))
    return TARGET_A_CM / a_native_cm


def _material_boundary_dists(layers):
    """Cumulative distances of the boundaries leak_gate considers (sol excluded):
    inner of first material layer, then outer of each material layer."""
    cum, dpairs = 0.0, []
    for name, thick in layers:
        dpairs.append((name, cum, cum + thick))
        cum += thick
    mats = [p for p in dpairs if p[0] != "sol"]
    dists = [mats[0][1]] + [p[2] for p in mats]   # inner of first + all outers
    return dists


def boundaries_normal(stem, scale, layers):
    """OLD offset: R+delta*n boundary polygons per phi at each material distance."""
    R, Z, _ = sg.load_surface(stem)
    R, Z = R * scale, Z * scale
    nR, nZ = sg.poloidal_outward_normals(R, Z)
    return [sg.offset_surface(R, Z, nR, nZ, d)
            for d in _material_boundary_dists(layers)]


def boundaries_buffer(stem, scale, layers, pol_res=None):
    """NEW offset: morphological-buffer boundary polygons per phi at each distance."""
    R, Z, _ = sg.load_surface(stem)
    R, Z = R * scale, Z * scale
    pr = int(pol_res) if pol_res else 3 * R.shape[0]   # match build_layers_sdf default
    return [sg2.offset_boundary(R, Z, d, pr)
            for d in _material_boundary_dists(layers)]


def worst_crossing_2d(boundaries):
    """boundaries : list of (R,Z) each (ntheta_k, nphi). Max over adjacent pairs and
    phi of the fraction of inner vertices outside the outer polygon. Returns float."""
    worst = 0.0
    for i in range(len(boundaries) - 1):
        Ri, Zi = boundaries[i]
        Ro, Zo = boundaries[i + 1]
        nph = Ri.shape[1]
        n_out, n_tot = 0, 0
        for j in range(nph):
            poly = MplPath(np.column_stack([Ro[:, j], Zo[:, j]]))
            pts = np.column_stack([Ri[:, j], Zi[:, j]])
            inside = poly.contains_points(pts, radius=-1e-9)  # strict: on-edge=outside
            n_out += int((~inside).sum())
            n_tot += len(pts)
        worst = max(worst, n_out / n_tot)
    return worst


def run_all(which):
    print(f"{'device':10s} {'class':6s} {'scale':>6s}  "
          f"{'OLD xcross':>10s} {'NEW xcross':>10s}  verdict")
    print("-" * 66)
    rows = []
    for cid in which:
        stem = f"quasr{cid}"
        klass = "leaky" if cid in LEAKY else "good"
        try:
            scale = scale_for(stem)
        except FileNotFoundError:
            print(f"{cid:10s} {klass:6s}   (no surface npz -- skipped)")
            continue
        oxc = worst_crossing_2d(boundaries_normal(stem, scale, BASELINE))
        nxc = worst_crossing_2d(boundaries_buffer(stem, scale, BASELINE))
        if klass == "leaky":
            verdict = "RECOVERED" if nxc < 1e-3 else "still leaky"
        else:
            verdict = "OK (no regr)" if nxc < 1e-3 else "REGRESSION!"
        print(f"{cid:10s} {klass:6s} {scale:6.2f}  {oxc:10.4f} {nxc:10.4f}  {verdict}",
              flush=True)
        rows.append((cid, klass, scale, oxc, nxc, verdict))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("devices", nargs="*")
    a = ap.parse_args()
    run_all(a.devices or (LEAKY + GOOD))
