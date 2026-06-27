#!/usr/bin/env python
"""STEP (x86/Ginsburg, conda): convert the conformal per-layer STL shell-solids
(written by stellarator_geometry.py) into a single watertight DAGMC .h5m with
material group tags, ready for OpenMC+DAGMC.

This REQUIRES the compiled DAGMC stack (cad_to_dagmc + moab/dagmc), which is not
available on the aarch64 sandbox -- it is part of the portable package that runs
on Ginsburg via environment.yml. Imports are guarded so running it here gives a
clear message rather than a traceback.

Path: stl_to_h5m (fusion-energy/stl_to_h5m) -- the dedicated STL->DAGMC tool.
NOTE: cad_to_dagmc does NOT read STL (it ingests STEP / CadQuery only), so it is
the wrong tool here; stl_to_h5m is purpose-built and runs make_watertight to merge
coincident surfaces (our nested shells share identical offset surfaces by
construction, so they weld cleanly). The per-layer material tag is the layer name;
run_conformal.py forces each OpenMC Material.name to that same tag (DAGMC matches
volumes to Material.name). VERIFY the stl_to_h5m call signature against the pinned
version (environment.yml) -- isolated at the single call site below.

The SOL/plasma interior is intentionally NOT a tagged volume; it becomes the DAGMC
implicit complement, which run_conformal assigns to vacuum.

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


def build(geom_dir, out_h5m="stellarator.h5m"):
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
        from stl_to_h5m import stl_to_h5m
    except ImportError as e:
        raise SystemExit(
            "stl_to_h5m not available -- this step runs on x86/Ginsburg in the "
            f"conda env (environment.yml). Underlying import error: {e}")

    files_with_tags = [
        {"material_tag": L["name"], "stl_filename": str(geom_dir / L["stl"])}
        for L in layers if L["name"] not in VACUUM_LAYERS
    ]
    tagged = [f["material_tag"] for f in files_with_tags]
    # >>> single API call site: verify against the pinned stl_to_h5m version <<<
    stl_to_h5m(files_with_tags=files_with_tags, h5m_filename=str(out_h5m))
    print(f"[build_dagmc] wrote {out_h5m} with material volumes: {tagged} "
          f"(SOL/plasma interior -> DAGMC implicit complement -> vacuum in run_conformal)")
    return out_h5m


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    build(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "stellarator.h5m")
