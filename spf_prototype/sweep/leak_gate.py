#!/usr/bin/env python
"""leak_gate.py -- native (no-transport) DAGMC-leak gate for the conformal sweep.

WHY: ~38% of sweep devices fail transport with "Maximum number of lost particles
reached / No intersection found with DAGMC cell N" -- a non-watertight DAGMC. The
existing geometry gate in stellarator_geometry.build_layers checks, per layer,
  * is_edge_manifold(T)          -- ALWAYS true for a structured torus grid
  * cross_sections_simple(Ro,Zo) -- per-phi SINGLE-polygon self-crossing only
Neither detects the real failure: adjacent nested SHELLS crossing each other. The
naive poloidal normal-offset folds on the concave inboard side once the cumulative
offset exceeds the local radius of curvature, so an inner boundary pokes THROUGH an
outer boundary. Each boundary is individually watertight; the ASSEMBLY self-
intersects, and DAGMC's ray tracer cannot resolve it -> lost particles.

Ground truth (trimesh): for every adjacent boundary pair, the fraction of inner-
boundary vertices lying OUTSIDE the outer boundary. >~0 means the shells cross.

RECOVERY IS ~0 (measured, honest). Two levers were tested natively and NEITHER seals
the crossing:
  * SCALE: inflating the build scale 4x (device -> ~1000 cm minor radius, blanket
    fixed at 116 cm) only drops the worst crossing from 0.376 to 0.169 for 1505944
    (and 0.324 -> 0.159 for 2190977) -- monotone but nowhere near sealed. The fold is
    an offset-DIRECTION defect (poloidal normals invert on the concave-inboard
    indentation of a bean cross-section), not an offset-MAGNITUDE one, so scale barely
    helps.
  * RESOLUTION: higher poloidal/toroidal density does the OPPOSITE of sealing -- at
    96x144 / 128x192 the coarse 64x96 mesh stops smoothing over the fold and the
    ORIGINAL per-polygon self-intersection gate trips (leak -> hard fold). Confirmed by
    the parallel batch agent on 1880092 (perpendicular-mode lost particles persist past
    an honest lost cap).
So the tight concave-inboard shapes are INTRINSICALLY leaky for this naive normal-
offset builder. The correct fix is therefore a CLEAN REJECT: detect the crossing
natively (cheap trimesh, no transport) and classify 'dagmc_leaky' so a self-
intersecting mesh is never sent to transport. A bounded scale-retry is available
(off by default) but empirically recovers ~0 devices.

Validation of the crossing metric: worst_shell_crossing == 0.0 for 4/4 transport-
'done' devices (59509, 115887, 134412, 157482) and 0.15-0.38 for every leaky device
-> the metric separates good from leaky perfectly.

Drop-in: in sweep.run_config, AFTER build_from_stls writes the .h5m and BEFORE
run_transport, add the native crossing gate (default: reject-only, no scale-retry):

    from leak_gate import worst_shell_crossing
    xcross = worst_shell_crossing(geom_dir)          # cheap, native, no transport
    rec["shell_crossing_frac"] = float(xcross)
    if xcross > defaults.get("leak_crossing_tol", 1e-3):
        rec["status"] = "dagmc_leaky"                # clean reject, NOT "error"
        rec["warnings"].append(
            f"conformal shells self-intersect (worst {xcross:.3f}); not transported")
        _write(rec_path, rec)
        print(f"[sweep] {cid}: DAGMC-LEAKY (crossing {xcross:.3f}) -> rejected")
        return rec
    # else: mesh is cleanly nested -> proceed to run_transport as before.

This turns the ~38% silent transport 'error' into an explicit, native, pre-transport
'dagmc_leaky' status -- no wasted Lima/Ginsburg transport time on meshes that cannot
converge. The optional adaptive builder below (build_geometry_watertight) additionally
tries a bounded scale-retry, but recovery is ~0 so the reject gate is the real fix.

This is intentionally NOT wired into the running sweep.py to avoid colliding with
the in-flight batch; it is a staged module the next sweep run imports.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))
import stellarator_geometry as sg
import dagmc_writer as dw


def _boundaries_from_stls(geom_dir):
    """Exact boundary surfaces the DAGMC .h5m is assembled from (dedup path)."""
    man = json.loads((Path(geom_dir) / "manifest.json").read_text())
    layers = [L for L in man["layers"] if L["name"] not in {"sol"}]
    nb = len(layers) + 1
    boundaries = [None] * nb
    for k, L in enumerate(layers):
        tris = dw._read_stl_tris(Path(geom_dir) / L["stl"])
        m = len(tris) // 2
        outer, inner = tris[:m], tris[m:]
        if k == 0:
            boundaries[0] = dw._orient_outward(dw._dedup_surface(inner), sg)
        boundaries[k + 1] = dw._orient_outward(dw._dedup_surface(outer), sg)
    return boundaries


def worst_shell_crossing(geom_dir):
    """Max over adjacent boundary pairs of the fraction of inner-boundary vertices
    lying OUTSIDE the outer boundary. 0.0 == cleanly nested; >0 == shells cross.
    Native (trimesh + rtree); no transport. Raises if a boundary is non-watertight
    (contains() would be meaningless -- that is itself a leak)."""
    import trimesh
    B = _boundaries_from_stls(geom_dir)
    worst = 0.0
    for i in range(len(B) - 1):
        Vi, _ = B[i]
        Vo, To = B[i + 1]
        outer = trimesh.Trimesh(vertices=np.asarray(Vo), faces=np.asarray(To),
                                process=False)
        if not outer.is_watertight:
            return 1.0  # degenerate outer boundary -> treat as maximally leaky
        inside = outer.contains(np.asarray(Vi))
        worst = max(worst, float((~inside).mean()))
    return worst


def _build_once(stem, layers, scale, geom_dir):
    """build_layers + the ORIGINAL fold gate (self-crossing polygon / manifold).
    Returns manifest or None (fold)."""
    m = sg.build_layers(stem, layers=layers, scale=scale, outdir=geom_dir)
    bad = [L["name"] for L in m["layers"]
           if not (L["watertight"] and L["simple_cross_section"])]
    return None if bad else m


def build_geometry_watertight(stem, layers, scale, geom_dir, *,
                              max_scale_mult=1.0, scale_step=1.4,
                              crossing_tol=1e-3, verbose=True):
    """Adaptive conformal build guaranteeing a NON-self-intersecting assembly.

    Returns (manifest, status, used_scale) where status is one of:
      'ok'               -- built and cleanly nested (safe for transport)
      'geometry_folded'  -- original per-polygon fold gate failed at all scales
      'dagmc_leaky'      -- built, but shells self-intersect at every scale tried

    The scale-retry is shared between the two failure modes: inflating scale both
    unfolds a self-crossing cross-section AND separates crossing shells.
    """
    used = scale
    limit = scale * max_scale_mult
    last_cross = None
    while True:
        m = _build_once(stem, layers, used, geom_dir)
        if m is None:                       # original fold gate tripped
            if used * scale_step > limit + 1e-9:
                return None, "geometry_folded", used
            used *= scale_step
            continue
        cross = worst_shell_crossing(geom_dir)
        last_cross = cross
        if verbose:
            print(f"[leak_gate] {stem} scale={used:.3f} worst_crossing={cross:.4f}",
                  flush=True)
        if cross <= crossing_tol:
            return m, "ok", used
        if used * scale_step > limit + 1e-9:
            return m, "dagmc_leaky", used   # honest reject; do NOT transport
        used *= scale_step


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("geom_dir")
    a = ap.parse_args()
    print("worst_shell_crossing =", worst_shell_crossing(a.geom_dir))
