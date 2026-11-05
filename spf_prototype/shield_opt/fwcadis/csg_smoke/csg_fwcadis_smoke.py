#!/usr/bin/env python3
"""Minimal CSG smoke test for OpenMC random-ray FW-CADIS.

This isolates the native multigroup/random-ray/weight-window path from DAGMC,
continuous-energy MGXS generation, and the stellarator geometry.  A successful
run must produce a readable, finite, non-empty ``weight_windows.h5`` file.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import openmc
from openmc.examples import random_ray_three_region_cube_with_detectors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def build_model() -> openmc.Model:
    model = random_ray_three_region_cube_with_detectors()

    # Local FW-CADIS should optimize only for the two detector responses.
    background_names = {"Source Tally", "Absorber Tally", "Cavity Tally"}
    detector_tallies = openmc.Tallies([
        tally for tally in model.tallies if tally.name not in background_names
    ])
    assert len(detector_tallies) == 2, [tally.name for tally in model.tallies]
    model.tallies = detector_tallies

    mesh = openmc.RegularMesh()
    mesh.dimension = (7, 7, 7)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (35.0, 35.0, 35.0)

    generator = openmc.WeightWindowGenerator(
        method="fw_cadis",
        targets=model.tallies,
        mesh=mesh,
        max_realizations=model.settings.batches,
    )
    model.settings.weight_window_generators = [generator]
    model.settings.random_ray["volume_estimator"] = "naive"
    return model


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(args.workdir)

    model = build_model()
    preflight = {
        "openmc_version": openmc.__version__,
        "openmc_module": str(Path(openmc.__file__).resolve()),
        "energy_mode": model.settings.energy_mode,
        "run_mode": model.settings.run_mode,
        "random_ray_keys": sorted(model.settings.random_ray),
        "target_tallies": [tally.name for tally in model.tallies],
    }
    print("PREFLIGHT " + json.dumps(preflight, sort_keys=True), flush=True)
    assert model.settings.energy_mode == "multi-group"
    assert model.settings.random_ray

    statepoint = model.run(threads=args.threads, output=True)
    ww_path = Path("weight_windows.h5")
    assert ww_path.is_file() and ww_path.stat().st_size > 0

    windows = openmc.hdf5_to_wws(ww_path)
    assert windows, "weight_windows.h5 contains no weight-window objects"
    lower = np.concatenate([np.asarray(ww.lower_ww_bounds).ravel() for ww in windows])
    upper = np.concatenate([np.asarray(ww.upper_ww_bounds).ravel() for ww in windows])
    assert np.all(np.isfinite(lower)) and np.all(np.isfinite(upper))
    assert np.any(upper > 0.0)
    assert np.all(upper >= lower)

    result = {
        "status": "PASS",
        "statepoint": str(statepoint),
        "weight_windows": str(ww_path.resolve()),
        "weight_windows_bytes": ww_path.stat().st_size,
        "objects": len(windows),
        "positive_upper_bins": int(np.count_nonzero(upper > 0.0)),
        "total_bins": int(upper.size),
    }
    Path("result.json").write_text(json.dumps(result, indent=2) + "\n")
    print("RESULT " + json.dumps(result, sort_keys=True), flush=True)
    print("DONE_CSG_FWCADIS_SMOKE", flush=True)


if __name__ == "__main__":
    main()
