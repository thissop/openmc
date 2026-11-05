#!/usr/bin/env python3
"""Run native random-ray FW-CADIS on the real one-period QH DAGMC model."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HELPER_DIR = Path(__file__).resolve().parent.parent / "reciprocity_cadis"
sys.path.insert(0, str(HELPER_DIR))
from common import build_materials, stellarator_source_point_cloud  # noqa: E402


def comma_ints(value: str) -> tuple[int, int, int]:
    result = tuple(int(item) for item in value.split(","))
    if len(result) != 3 or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("expected three positive comma-separated integers")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--dagmc", type=Path, required=True)
    parser.add_argument("--fluxmap", type=Path, required=True)
    parser.add_argument("--mgxs", type=Path, required=True)
    parser.add_argument("--rays", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=12)
    parser.add_argument("--inactive", type=int, default=3)
    parser.add_argument("--mgxs-particles", type=int, default=200)
    parser.add_argument("--source-points", type=int, default=2_000)
    parser.add_argument("--mesh", type=comma_ints, default=(10, 10, 10))
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--li6", type=float, default=60.0)
    return parser.parse_args()


def prepare_mgxs_library(openmc, materials, path: Path, particles: int) -> bool:
    """Generate physical-material MGXS and add explicit near-vacuum entries.

    The generic infinite-medium generator is pathological for the two dilute
    hydrogen placeholders: a 1000-m reflective box gives almost every history
    the maximum number of events. Their optical depths in the actual device are
    negligible, so represent them as effectively void two-group absorbers.

    Returns ``True`` when an existing complete library is reused.
    """
    dilute_base_names = {"gap", "Vacuum"}
    runtime_names = {
        material.name: f"{material.name}_{material.id}" for material in materials
    }
    required_names = set(runtime_names.values())
    dilute_names = {runtime_names[name] for name in dilute_base_names}

    def canonicalize_names(library) -> None:
        # The isolated generator appends the material ID to its XSdata names
        # when the same Material objects also belong to the DAGMC model. DAGMC
        # synchronization uses those suffixed names for the runtime macros, so
        # make an older base-name cache conform to that convention.
        for base_name, runtime_name in runtime_names.items():
            xsdata = library.get_by_name(base_name)
            if xsdata is not None and runtime_name not in library.names:
                xsdata.name = runtime_name

    if path.is_file():
        library = openmc.MGXSLibrary.from_hdf5(path)
        canonicalize_names(library)
        missing = required_names.difference(library.names)
        if not missing:
            library.export_to_hdf5(path)
            return True
        if not missing.issubset(dilute_names):
            raise RuntimeError(f"Existing MGXS library is incomplete; missing {sorted(missing)}")
    else:
        physical_materials = openmc.Materials([
            material for material in materials if material.name not in dilute_base_names
        ])
        groups = openmc.mgxs.EnergyGroups("CASMO-2")
        generator_root = openmc.Cell(
            region=-openmc.Sphere(r=1.0, boundary_type="vacuum")
        )
        generator_model = openmc.Model(
            geometry=openmc.Geometry([generator_root]),
            materials=physical_materials,
        )
        original_names = {material.id: material.name for material in materials}
        try:
            generator_model.convert_to_multigroup(
                method="infinite_medium",
                groups=groups,
                nparticles=particles,
                mgxs_path=path,
                source_energy=openmc.stats.Discrete([14.06e6], [1.0]),
            )
        finally:
            # The isolated MGXS workflow may decorate shared Material names
            # with IDs. Restore DAGMC tag names before synchronizing the real
            # model; the library itself intentionally keeps the ID suffixes.
            for material in materials:
                material.name = original_names[material.id]
        library = openmc.MGXSLibrary.from_hdf5(path)
        canonicalize_names(library)

    groups = library.energy_groups
    num_groups = groups.num_groups
    epsilon = np.full(num_groups, 1.0e-12)
    zero_scatter = np.zeros((num_groups, num_groups, 1))
    for name in sorted(required_names.difference(library.names)):
        xsdata = openmc.XSdata(name, groups)
        xsdata.order = 0
        xsdata.set_total(epsilon)
        xsdata.set_absorption(epsilon)
        xsdata.set_scatter_matrix(zero_scatter)
        library.add_xsdata(xsdata)
    library.export_to_hdf5(path)

    missing = required_names.difference(library.names)
    if missing:
        raise RuntimeError(f"Generated MGXS library is incomplete; missing {sorted(missing)}")
    return False


def main() -> None:
    args = parse_args()
    if args.batches <= args.inactive or args.inactive < 0:
        raise ValueError("batches must be greater than inactive")
    if not args.dagmc.is_file():
        raise FileNotFoundError(args.dagmc)

    import openmc

    args.workdir.mkdir(parents=True, exist_ok=True)
    args.mgxs.parent.mkdir(parents=True, exist_ok=True)
    os.chdir(args.workdir)

    materials = build_materials(openmc, args.li6)
    dagmc = openmc.DAGMCUniverse(str(args.dagmc), auto_geom_ids=True)
    dagmc_box = dagmc.bounding_box
    extent = float(np.max(np.abs(np.r_[dagmc_box.lower_left, dagmc_box.upper_right])))
    outer = openmc.Sphere(r=1.5 * extent + 10.0, boundary_type="vacuum")
    root_cell = openmc.Cell(region=-outer, fill=dagmc)
    geometry = openmc.Geometry([root_cell])

    source_points, source_strengths = stellarator_source_point_cloud(
        args.fluxmap, max_points=args.source_points
    )
    source_strengths = source_strengths / np.sum(source_strengths)
    source_angle = openmc.stats.Isotropic()
    source_energy = openmc.stats.Discrete([14.06e6], [1.0])
    sources = [
        openmc.IndependentSource(
            space=openmc.stats.Point(position),
            angle=source_angle,
            energy=source_energy,
            strength=float(strength),
        )
        for position, strength in zip(source_points, source_strengths)
    ]
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.source = sources
    settings.batches = args.batches
    settings.inactive = args.inactive
    settings.particles = args.rays

    fast = openmc.EnergyFilter([0.1e6, 25.0e6])
    coil_cells = openmc.CellFilter([12, 13, 14, 15, 16])
    coil_flux = openmc.Tally(name="coil_fast_flux")
    coil_flux.filters = [coil_cells, fast]
    coil_flux.scores = ["flux"]
    tallies = openmc.Tallies([coil_flux])
    model = openmc.Model(geometry, materials, settings, tallies)

    mgxs_existed = prepare_mgxs_library(
        openmc, materials, args.mgxs, args.mgxs_particles
    )
    print(
        "MGXS " + json.dumps({
            "path": str(args.mgxs),
            "reused": mgxs_existed,
            "particles_per_material_batch": args.mgxs_particles,
            "explicit_dilute_media": ["Vacuum", "gap"],
        }, sort_keys=True),
        flush=True,
    )
    model.convert_to_multigroup(
        method="infinite_medium",
        groups="CASMO-2",
        mgxs_path=args.mgxs,
    )
    model.convert_to_random_ray()

    # convert_to_random_ray supplies a geometric default; the experiment's ray
    # count must be restored after conversion or the CLI argument is ignored.
    model.settings.particles = args.rays
    model.settings.batches = args.batches
    model.settings.inactive = args.inactive

    widths = np.asarray(dagmc_box.upper_right) - np.asarray(dagmc_box.lower_left)
    chord = float(np.linalg.norm(widths))
    model.settings.random_ray["distance_inactive"] = 0.5 * chord
    model.settings.random_ray["distance_active"] = 2.0 * chord
    model.settings.random_ray["source_shape"] = "flat"
    model.settings.random_ray["volume_estimator"] = "naive"

    mesh = openmc.RegularMesh()
    mesh.dimension = args.mesh
    mesh.lower_left = tuple(float(value) for value in dagmc_box.lower_left)
    mesh.upper_right = tuple(float(value) for value in dagmc_box.upper_right)
    model.settings.random_ray["source_region_meshes"] = [
        (mesh, [geometry.root_universe])
    ]

    generator = openmc.WeightWindowGenerator(
        method="fw_cadis",
        targets=model.tallies,
        mesh=mesh,
        max_realizations=args.batches - args.inactive,
    )
    model.settings.weight_window_generators = [generator]

    preflight = {
        "openmc_version": openmc.__version__,
        "openmc_module": str(Path(openmc.__file__).resolve()),
        "energy_mode": model.settings.energy_mode,
        "run_mode": model.settings.run_mode,
        "rays_per_batch": model.settings.particles,
        "batches": model.settings.batches,
        "inactive": model.settings.inactive,
        "mesh": list(args.mesh),
        "target_tallies": [tally.name for tally in model.tallies],
        "source_type": type(model.settings.source[0]).__name__,
        "source_points": int(source_points.shape[0]),
        "random_ray_keys": sorted(model.settings.random_ray),
    }
    print("PREFLIGHT " + json.dumps(preflight, sort_keys=True), flush=True)
    if model.settings.energy_mode != "multi-group" or not generator.targets:
        raise RuntimeError("FW-CADIS preflight did not retain multigroup mode and targets")

    statepoint = model.run(threads=args.threads, output=True)
    ww_path = Path("weight_windows.h5")
    if not ww_path.is_file() or ww_path.stat().st_size == 0:
        raise RuntimeError("OpenMC did not create weight_windows.h5")
    windows = openmc.hdf5_to_wws(ww_path)
    lower = np.concatenate([window.lower_ww_bounds.ravel() for window in windows])
    result = {
        "status": "PASS",
        "statepoint": str(Path(statepoint).resolve()),
        "weight_windows": str(ww_path.resolve()),
        "weight_windows_bytes": ww_path.stat().st_size,
        "mgxs": str(args.mgxs.resolve()),
        "mgxs_reused": mgxs_existed,
        "objects": len(windows),
        "active_bins": int(np.count_nonzero(lower >= 0.0)),
        "total_bins": int(lower.size),
    }
    Path("result.json").write_text(json.dumps(result, indent=2) + "\n")
    print("RESULT " + json.dumps(result, sort_keys=True), flush=True)
    print("DONE_DAGMC_FWCADIS", flush=True)
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
