#!/usr/bin/env python
"""STEP (x86/Ginsburg, conda): convert the conformal per-layer STL shell-solids
(written by stellarator_geometry.py) into a single watertight DAGMC .h5m with
material group tags, ready for OpenMC+DAGMC.

This REQUIRES the compiled DAGMC stack (cad_to_dagmc + moab/dagmc), which is not
available on the aarch64 sandbox -- it is part of the portable package that runs
on Ginsburg via environment.yml. Imports are guarded so running it here gives a
clear message rather than a traceback.

Primary path: cad_to_dagmc (maintained, handles the DAGMC entity-set structure).
The material tag for each layer is its name; OpenMC materials in run_conformal.py
must use the same names. VERIFY the add_stl_file signature against the pinned
cad_to_dagmc version (environment.yml) -- the API has shifted across releases; the
single call site is isolated below for a one-line fix if needed.

Usage (Ginsburg):
    python build_dagmc.py <geom_dir> [out.h5m]
e.g. python build_dagmc.py data/equil_precise_qa_geom stellarator.h5m
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Layers that are real materials (the 'sol' gap is vacuum -> not a tagged volume).
VACUUM_LAYERS = {"sol"}


def build(geom_dir, out_h5m="stellarator.h5m", min_mesh_size=2.0, max_mesh_size=20.0):
    geom_dir = Path(geom_dir)
    manifest = json.loads((geom_dir / "manifest.json").read_text())
    layers = manifest["layers"]
    # hard gate: never feed an invalid (self-intersecting) build to the mesher
    bad = [L["name"] for L in layers
           if not (L["watertight"] and L["simple_cross_section"])]
    if bad:
        raise RuntimeError(
            f"refusing to build DAGMC: layers {bad} are not watertight/simple "
            "(over-thick conformal offset folds inboard). Increase scale or thin "
            "the build (see DEFERRED.md).")

    try:
        from cad_to_dagmc import CadToDagmc
    except ImportError as e:
        raise SystemExit(
            "cad_to_dagmc not available -- this step runs on x86/Ginsburg in the "
            f"conda env (environment.yml). Underlying import error: {e}")

    c = CadToDagmc()
    tagged = []
    for L in layers:
        if L["name"] in VACUUM_LAYERS:
            continue
        stl = geom_dir / L["stl"]
        # >>> single API call site: verify against the pinned cad_to_dagmc version <<<
        c.add_stl_file(str(stl), material_tags=[L["name"]])
        tagged.append(L["name"])
    c.export_dagmc_h5m_file(
        filename=out_h5m,
        min_mesh_size=min_mesh_size, max_mesh_size=max_mesh_size,
    )
    print(f"[build_dagmc] wrote {out_h5m} with material volumes: {tagged}")
    return out_h5m


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    build(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "stellarator.h5m")
