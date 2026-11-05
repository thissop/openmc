#!/usr/bin/env python
"""Benchmark one variance-reduction map on the validated one-period SSR model.

Every invocation uses the same deep-freeze surface source, periodic wedge,
black-inside material overrides, tallies, seed, and split-safety limits.  The
mode label identifies the candidate; it does not silently select an input file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback

import numpy as np
import openmc

from metrics import figures_of_merit, relative_errors


COIL_CELLS = [12, 13, 14, 15, 16]
DEPTH_CELLS = [6, 8, 10, 12, 14]
FIXED_INSIDE_CELLS = [2, 3, 4]
FAST_ENERGY_BOUNDS_EV = [0.1e6, 25.0e6]
MODES = ("analog", "magic", "native", "reciprocity")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--dagmc", type=Path, required=True)
    parser.add_argument("--surface-source", type=Path, required=True)
    parser.add_argument("--weight-windows", type=Path)
    parser.add_argument("--batches", type=int, default=4)
    parser.add_argument("--particles", type=int, default=50_000)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--seed", type=int, default=104_729)
    parser.add_argument("--li6", type=float, default=60.0)
    parser.add_argument("--max-split", type=int, default=4)
    parser.add_argument("--max-history-splits", type=int, default=100)
    parser.add_argument("--max-lower-bound-ratio", type=float, default=5.0)
    # Cell-id overrides: the subdivided-shield geometry shifts coil/depth cell ids.
    # Pass explicitly, or --cellmap-json to read them from cellmap_subdiv.json.
    # Defaults reproduce the original single-shield-cell geometry exactly.
    parser.add_argument("--coil-cells", type=int, nargs="+")
    parser.add_argument("--depth-cells", type=int, nargs="+")
    parser.add_argument("--black-cells", type=int, nargs="+")
    parser.add_argument("--cellmap-json", type=Path)
    args = parser.parse_args()

    if args.cellmap_json is not None:
        cm = json.loads(args.cellmap_json.read_text())
        if args.coil_cells is None:
            args.coil_cells = list(cm["coil_cells"])
        if args.depth_cells is None:
            args.depth_cells = list(cm["depth_cells"])
        if args.black_cells is None:
            args.black_cells = list(cm["black_inside"])
    if args.coil_cells is None:
        args.coil_cells = list(COIL_CELLS)
    if args.depth_cells is None:
        args.depth_cells = list(DEPTH_CELLS)
    if args.black_cells is None:
        args.black_cells = list(FIXED_INSIDE_CELLS)
    if len(args.coil_cells) < 1:
        parser.error(f"need at least one coil cell, got {args.coil_cells}")

    if args.batches < 2:
        parser.error("--batches must be at least 2 for batch statistics")
    if args.particles <= 0 or args.threads <= 0 or args.seed <= 0:
        parser.error("particles, threads, and seed must be positive")
    if (
        args.max_split <= 1
        or args.max_history_splits <= 0
        or args.max_lower_bound_ratio <= 1.0
    ):
        parser.error("split limits must be positive and max-split must exceed 1")
    if args.mode == "analog" and args.weight_windows is not None:
        parser.error("analog mode must not receive --weight-windows")
    if args.mode != "analog" and args.weight_windows is None:
        parser.error(f"{args.mode} mode requires --weight-windows")
    return args


def build_materials(li6: float) -> openmc.Materials:
    tungsten = openmc.Material()
    tungsten.add_element("W", 1.0)
    tungsten.set_density("g/cm3", 19.3)

    rafm = openmc.Material()
    for element, fraction in [
        ("Fe", 0.885), ("Cr", 0.09), ("W", 0.011), ("Mn", 0.006),
        ("V", 0.002), ("C", 0.006),
    ]:
        rafm.add_element(element, fraction, "wo")
    rafm.set_density("g/cm3", 7.9)

    first_wall = openmc.Material.mix_materials(
        [tungsten, rafm], [0.2 / 3.2, 3.0 / 3.2], "vo"
    )
    first_wall.name = "first_wall"

    multiplier = openmc.Material(name="multiplier")
    multiplier.add_element("Be", 1.0)
    multiplier.set_density("g/cm3", 1.85)

    breeder = openmc.Material(name="breeder")
    breeder.add_element(
        "Li", 2.0, enrichment=li6, enrichment_target="Li6", enrichment_type="ao"
    )
    breeder.add_element("Be", 1.0)
    breeder.add_element("F", 4.0)
    breeder.set_density("g/cm3", 1.94)

    back_wall = openmc.Material(name="back_wall")
    for element, fraction in [
        ("Fe", 0.885), ("Cr", 0.09), ("W", 0.011), ("Mn", 0.006),
        ("V", 0.002), ("C", 0.006),
    ]:
        back_wall.add_element(element, fraction, "wo")
    back_wall.set_density("g/cm3", 7.9)

    shield = openmc.Material(name="shield")
    shield.add_element("W", 1.0)
    shield.add_element("C", 1.0)
    shield.set_density("g/cm3", 15.6)

    gap = openmc.Material(name="gap")
    gap.add_nuclide("H1", 1.0)
    gap.set_density("g/cm3", 1e-8)

    vessel = openmc.Material(name="vac_vessel")
    for element, fraction in [
        ("Fe", 0.62), ("Cr", 0.16), ("Ni", 0.11), ("Mo", 0.02),
        ("Mn", 0.015), ("B", 0.03), ("O", 0.03), ("H", 0.005),
    ]:
        vessel.add_element(element, fraction, "wo")
    vessel.set_density("g/cm3", 7.9)

    thermal_shield = openmc.Material(name="thermal_shield")
    thermal_shield.add_element("Fe", 0.89, "wo")
    thermal_shield.add_element("Cr", 0.09, "wo")
    thermal_shield.add_element("Mn", 0.02, "wo")
    thermal_shield.set_density("g/cm3", 7.9)

    magnets = openmc.Material(name="magnets")
    for element, fraction in [
        ("Cu", 0.58), ("Fe", 0.30), ("Cr", 0.04), ("Ni", 0.03),
        ("Y", 0.02), ("Ba", 0.02), ("O", 0.01),
    ]:
        magnets.add_element(element, fraction, "wo")
    magnets.set_density("g/cm3", 8.5)

    vacuum = openmc.Material(name="Vacuum")
    vacuum.add_nuclide("H1", 1.0)
    vacuum.set_density("g/cm3", 1e-10)

    fixed_black = openmc.Material(name="fixed_black")
    fixed_black.add_nuclide("B10", 1.0)
    fixed_black.set_density("g/cm3", 10.0)

    return openmc.Materials([
        first_wall, multiplier, breeder, back_wall, shield, gap, vessel,
        thermal_shield, magnets, vacuum, fixed_black,
    ])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_metadata(path: Path, *, include_sha256: bool) -> dict:
    stat = path.stat()
    result = {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_sha256:
        result["sha256"] = sha256(path)
    return result


def weight_window_metadata(wws, dagmc_bounds: openmc.BoundingBox) -> list[dict]:
    metadata = []
    for ww in wws:
        mesh = ww.mesh
        # RegularMesh (cartesian) maps are validated against the DAGMC bounding box. CylindricalMesh
        # maps (the FW-CADIS / reciprocity-CADIS subdivided-shield maps) are torus-aligned and legitimately
        # do NOT span the cartesian bb -- validate only that they carry active bins.
        is_regular = isinstance(mesh, openmc.RegularMesh)
        if is_regular:
            if not hasattr(mesh, "lower_left") or not hasattr(mesh, "upper_right"):
                raise TypeError("benchmark requires spatially bounded structured weight windows")
            lower_left = np.asarray(mesh.lower_left, dtype=float)
            upper_right = np.asarray(mesh.upper_right, dtype=float)
            if not np.allclose(lower_left, dagmc_bounds.lower_left, rtol=0.0, atol=1e-8):
                raise ValueError("weight-window lower-left bound does not match the DAGMC model")
            if not np.allclose(upper_right, dagmc_bounds.upper_right, rtol=0.0, atol=1e-8):
                raise ValueError("weight-window upper-right bound does not match the DAGMC model")

        lower = np.asarray(ww.lower_ww_bounds)
        upper = np.asarray(ww.upper_ww_bounds)
        active = lower >= 0.0
        if not np.any(active):
            raise ValueError("weight-window object has no active bins")
        if np.any(lower[active] <= 0.0) or np.any(upper[active] <= 0.0):
            raise ValueError("active weight-window bounds must be positive")
        energy_bounds = np.asarray(ww.energy_bounds, dtype=float)
        if energy_bounds[-1] <= FAST_ENERGY_BOUNDS_EV[0]:
            raise ValueError("weight windows do not overlap the coil fast-flux tally")

        metadata.append({
            "particle_type": str(ww.particle_type),
            "mesh_type": type(mesh).__name__,
            "mesh_dimension": [int(value) for value in mesh.dimension],
            "lower_left_cm": (lower_left.tolist() if is_regular else None),
            "upper_right_cm": (upper_right.tolist() if is_regular else None),
            "energy_bounds_ev": energy_bounds.tolist(),
            "active_bins": int(np.count_nonzero(active)),
            "total_bins": int(lower.size),
            "survival_ratio": float(ww.survival_ratio),
            "max_lower_bound_ratio": (
                None if ww.max_lower_bound_ratio is None
                else float(ww.max_lower_bound_ratio)
            ),
            "max_split": int(ww.max_split),
        })
    return metadata


def build_model(args: argparse.Namespace):
    materials = build_materials(args.li6)
    fixed_black = next(mat for mat in materials if mat.name == "fixed_black")

    dagmc = openmc.DAGMCUniverse(args.dagmc, auto_geom_ids=True)
    bounds = dagmc.bounding_box
    extent = float(np.max(np.abs(np.concatenate([bounds.lower_left, bounds.upper_right]))))
    sphere = openmc.Sphere(r=extent * 1.5 + 10.0, boundary_type="vacuum")
    phi_zero = openmc.YPlane(0.0, boundary_type="periodic")
    phi_ninety = openmc.XPlane(0.0, boundary_type="periodic")
    phi_zero.periodic_surface = phi_ninety
    root = openmc.Cell(region=+phi_zero & +phi_ninety & -sphere, fill=dagmc)
    geometry = openmc.Geometry([root])

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.surf_source_read = {"path": args.surface_source}
    settings.batches = args.batches
    settings.particles = args.particles
    settings.seed = args.seed
    settings.photon_transport = False
    settings.max_lost_particles = 500
    settings.rel_max_lost_particles = 2e-3
    settings.output = {"summary": True, "tallies": False}
    settings.verbosity = 6

    ww_metadata = []
    if args.weight_windows is not None:
        wws = openmc.hdf5_to_wws(args.weight_windows)
        if not wws:
            raise ValueError("weight-window file contains no objects")
        for ww in wws:
            ww.max_split = args.max_split
            # Generated bounds can sit many decades below a source-bank
            # particle's starting weight. Apply the same dynamic-rescaling
            # guard to every candidate before allowing secondary cascades.
            ww.max_lower_bound_ratio = args.max_lower_bound_ratio
        ww_metadata = weight_window_metadata(wws, bounds)
        settings.weight_windows = wws
        settings.weight_windows_on = True
        settings.max_history_splits = args.max_history_splits

    cell_filter = openmc.CellFilter(args.coil_cells)
    fast_filter = openmc.EnergyFilter(FAST_ENERGY_BOUNDS_EV)
    coil_flux = openmc.Tally(name="coil_fast_flux")
    coil_flux.filters = [cell_filter, fast_filter]
    coil_flux.scores = ["flux"]
    coil_heating = openmc.Tally(name="coil_heating")
    coil_heating.filters = [cell_filter]
    coil_heating.scores = ["heating-local"]
    depth_flux = openmc.Tally(name="depth_flux")
    depth_flux.filters = [openmc.CellFilter(args.depth_cells)]
    depth_flux.scores = ["flux"]
    tallies = openmc.Tallies([coil_flux, coil_heating, depth_flux])

    model = openmc.Model(
        geometry=geometry, materials=materials, settings=settings, tallies=tallies
    )
    model.init_lib(output=False, threads=args.threads)
    try:
        model.sync_dagmc_universes()
        for cell_id in args.black_cells:
            dagmc.add_material_override(cell_id, fixed_black)
        for cell_id in args.black_cells:
            actual = getattr(dagmc.cells[cell_id].fill, "name", None)
            if actual != "fixed_black":
                raise RuntimeError(
                    f"black-inside override failed for cell {cell_id}: {actual!r}"
                )
    finally:
        model.finalize_lib()

    return model, ww_metadata


def run_openmc(workdir: Path, threads: int) -> tuple[Path, str, float]:
    command = ["openmc", "-s", str(threads)]
    environment = os.environ.copy()
    source_library = environment.get("OPENMC_SOURCE_LIB")
    if source_library:
        environment["LD_PRELOAD"] = source_library
    output_path = workdir / "openmc.out"
    started = time.monotonic()
    lines = []
    with output_path.open("w", encoding="utf-8") as output_file:
        process = subprocess.Popen(
            command,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=environment,
        )
        assert process.stdout is not None
        for line in process.stdout:
            lines.append(line)
            output_file.write(line)
            output_file.flush()
            print(line, end="", flush=True)
        return_code = process.wait()
    wall_seconds = time.monotonic() - started
    output = "".join(lines)
    if return_code != 0:
        raise RuntimeError(f"OpenMC exited with status {return_code}")
    statepoints = sorted(workdir.glob("statepoint.*.h5"), key=lambda path: path.stat().st_mtime_ns)
    if not statepoints:
        raise RuntimeError("OpenMC completed without writing a statepoint")
    return statepoints[-1], output, wall_seconds


def parse_rate(output: str, label: str) -> float | None:
    match = re.search(
        rf"{re.escape(label)}\s*=\s*([0-9.+\-Ee]+)\s+(?:particles|tracks)/second",
        output,
    )
    return float(match.group(1)) if match else None


def tally_payload(tally, simulation_seconds: float) -> dict:
    mean = tally.mean.ravel().astype(float).tolist()
    std_dev = tally.std_dev.ravel().astype(float).tolist()
    errors = relative_errors(mean, std_dev)
    return {
        "mean": mean,
        "std_dev": std_dev,
        "relative_error": [
            error if math.isfinite(error) else None for error in errors
        ],
        "figure_of_merit": figures_of_merit(errors, simulation_seconds),
    }


def execute(args: argparse.Namespace) -> dict:
    args.dagmc = args.dagmc.resolve()
    args.surface_source = args.surface_source.resolve()
    if args.weight_windows is not None:
        args.weight_windows = args.weight_windows.resolve()
    for path in [args.dagmc, args.surface_source, args.weight_windows]:
        if path is not None and not path.is_file():
            raise FileNotFoundError(path)

    args.workdir = args.workdir.resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)
    if (args.workdir / "result.json").exists():
        raise FileExistsError(f"refusing to replace {args.workdir / 'result.json'}")
    os.chdir(args.workdir)

    model, ww_metadata = build_model(args)
    model.export_to_model_xml()
    statepoint, output, wall_seconds = run_openmc(args.workdir, args.threads)

    calculation_rate = parse_rate(output, "Calculation Rate (active)")
    track_rate = parse_rate(output, "Track Rate (active)")
    tracks_per_source = (
        track_rate / calculation_rate
        if track_rate is not None and calculation_rate not in (None, 0.0)
        else (1.0 if args.mode == "analog" else None)
    )

    with openmc.StatePoint(statepoint) as state:
        simulation_seconds = float(state.runtime["simulation"])
        runtime = {key: float(value) for key, value in state.runtime.items()}
        flux = tally_payload(state.get_tally(name="coil_fast_flux"), simulation_seconds)
        heating = tally_payload(state.get_tally(name="coil_heating"), simulation_seconds)
        depth = tally_payload(state.get_tally(name="depth_flux"), simulation_seconds)

    version_environment = os.environ.copy()
    source_library = version_environment.get("OPENMC_SOURCE_LIB")
    if source_library:
        version_environment["LD_PRELOAD"] = source_library
    executable_version = subprocess.run(
        ["openmc", "--version"],
        check=True,
        text=True,
        capture_output=True,
        env=version_environment,
    ).stdout.strip()
    return {
        "status": "PASS",
        "mode": args.mode,
        "openmc_python_version": openmc.__version__,
        "openmc_executable_version": executable_version,
        "execution": {
            "hostname": os.uname().nodename,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        },
        "workdir": str(args.workdir),
        "statepoint": str(statepoint),
        "configuration": {
            "batches": args.batches,
            "particles_per_batch": args.particles,
            "source_histories": args.batches * args.particles,
            "threads": args.threads,
            "seed": args.seed,
            "li6_atom_percent": args.li6,
            "black_inside_cells": args.black_cells,
            "coil_cells": args.coil_cells,
            "depth_cells": args.depth_cells,
            "fast_energy_bounds_ev": FAST_ENERGY_BOUNDS_EV,
            "max_split": None if args.mode == "analog" else args.max_split,
            "max_history_splits": (
                None if args.mode == "analog" else args.max_history_splits
            ),
            "max_lower_bound_ratio": (
                None if args.mode == "analog" else args.max_lower_bound_ratio
            ),
        },
        "inputs": {
            "dagmc": file_metadata(args.dagmc, include_sha256=True),
            "surface_source": file_metadata(args.surface_source, include_sha256=False),
            "weight_windows": (
                None
                if args.weight_windows is None
                else file_metadata(args.weight_windows, include_sha256=True)
            ),
        },
        "weight_window_objects": ww_metadata,
        "runtime": runtime,
        "runner_wall_seconds": wall_seconds,
        "calculation_rate_particles_per_second": calculation_rate,
        "track_rate_per_second": track_rate,
        "estimated_tracks_per_source_history": tracks_per_source,
        "coil_fast_flux": flux,
        "coil_heating_local": heating,
        "depth_flux": depth,
    }


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    result_path = args.workdir.resolve() / "result.json"
    try:
        result = execute(args)
    except Exception as error:
        args.workdir.mkdir(parents=True, exist_ok=True)
        failure = {
            "status": "FAIL",
            "mode": args.mode,
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        }
        if not result_path.exists():
            write_json(result_path, failure)
        raise
    write_json(result_path, result)
    print(f"WROTE {result_path}", flush=True)


if __name__ == "__main__":
    main()
