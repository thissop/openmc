#!/usr/bin/env python3
"""Generate a one-period reciprocal-transport importance map and CADIS windows.

This is the robust Path-B fallback.  An isotropic 14.06-MeV point cloud on the
first five QH coil filaments is transported inward through the continuous-energy
DAGMC model.  By neutron reciprocity, the resulting spatial flux is used as an
approximation to the coil-response adjoint importance.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

from cadis import build_cadis_bounds
from common import build_materials, read_makegrid, sample_coil_points


def comma_ints(value: str) -> tuple[int, int, int]:
    result = tuple(int(item) for item in value.split(","))
    if len(result) != 3 or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("expected three positive comma-separated integers")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--dagmc", type=Path, required=True)
    parser.add_argument("--coils", type=Path, required=True)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--particles", type=int, default=20_000)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--mesh", type=comma_ints, default=(16, 16, 8))
    parser.add_argument("--coil-count", type=int, default=5)
    parser.add_argument("--points-per-coil", type=int, default=64)
    parser.add_argument("--li6", type=float, default=60.0)
    parser.add_argument("--source-energy", type=float, default=14.06e6)
    parser.add_argument("--minimum-energy", type=float, default=0.1e6)
    parser.add_argument("--maximum-energy", type=float, default=25.0e6)
    parser.add_argument("--relative-error", type=float, default=0.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.batches < 2 or args.particles < 1:
        raise ValueError("batches must be at least two and particles must be positive")
    if not args.dagmc.is_file() or not args.coils.is_file():
        raise FileNotFoundError((args.dagmc, args.coils))

    import openmc

    args.workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(args.workdir)

    materials = build_materials(openmc, args.li6)
    dagmc = openmc.DAGMCUniverse(str(args.dagmc), auto_geom_ids=True)
    dagmc_box = dagmc.bounding_box
    extent = float(np.max(np.abs(np.r_[dagmc_box.lower_left, dagmc_box.upper_right])))
    outer = openmc.Sphere(r=1.5 * extent + 10.0, boundary_type="vacuum")
    y_plane = openmc.YPlane(0.0, boundary_type="periodic")
    x_plane = openmc.XPlane(0.0, boundary_type="periodic")
    y_plane.periodic_surface = x_plane
    root_cell = openmc.Cell(region=+y_plane & +x_plane & -outer, fill=dagmc)
    geometry = openmc.Geometry([root_cell])

    coils = read_makegrid(args.coils)[: args.coil_count]
    if len(coils) != args.coil_count:
        raise ValueError(f"Requested {args.coil_count} coils but found {len(coils)}")
    source_points = sample_coil_points(coils, args.points_per_coil)
    source = openmc.IndependentSource(
        space=openmc.stats.PointCloud(source_points),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([args.source_energy], [1.0]),
    )

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.source = [source]
    settings.batches = args.batches
    settings.particles = args.particles
    settings.photon_transport = False
    settings.max_lost_particles = 500
    settings.rel_max_lost_particles = 2.0e-3

    mesh = openmc.RegularMesh()
    mesh.dimension = args.mesh
    mesh.lower_left = tuple(float(value) for value in dagmc_box.lower_left)
    mesh.upper_right = tuple(float(value) for value in dagmc_box.upper_right)
    energy_bounds = np.array([args.minimum_energy, args.maximum_energy])
    importance = openmc.Tally(name="reciprocal_coil_importance")
    importance.filters = [openmc.MeshFilter(mesh), openmc.EnergyFilter(energy_bounds)]
    importance.scores = ["flux"]

    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=openmc.Tallies([importance]),
    )
    preflight = {
        "openmc_version": openmc.__version__,
        "openmc_module": str(Path(openmc.__file__).resolve()),
        "dagmc": str(args.dagmc),
        "coils": str(args.coils),
        "source_points": int(source_points.shape[0]),
        "mesh": list(args.mesh),
        "histories": args.batches * args.particles,
        "energy_bounds_eV": energy_bounds.tolist(),
    }
    print("PREFLIGHT " + json.dumps(preflight, sort_keys=True), flush=True)

    statepoint_path = model.run(threads=args.threads, output=True)
    with openmc.StatePoint(statepoint_path) as statepoint:
        tally = statepoint.get_tally(name=importance.name)
        mean = tally.get_reshaped_data("mean", expand_dims=True)[..., 0, 0]
        std_dev = tally.get_reshaped_data("std_dev", expand_dims=True)[..., 0, 0]

    expected_shape = (*args.mesh, len(energy_bounds) - 1)
    if mean.shape != expected_shape:
        raise RuntimeError(f"Unexpected importance shape {mean.shape}; expected {expected_shape}")
    widths = (np.asarray(mesh.upper_right) - np.asarray(mesh.lower_left)) / np.asarray(args.mesh)
    cell_volume = float(np.prod(widths))
    bounds = build_cadis_bounds(
        mean,
        std_dev,
        cell_volume,
        relative_error_threshold=args.relative_error,
    )

    np.savez_compressed(
        "reciprocal_importance.npz",
        mean=mean,
        std_dev=std_dev,
        relative_error=bounds.relative_error,
        reliable=bounds.reliable,
        lower_ww_bounds=bounds.lower,
        upper_ww_bounds=bounds.upper,
        mesh_dimension=np.asarray(args.mesh),
        mesh_lower_left=np.asarray(mesh.lower_left),
        mesh_upper_right=np.asarray(mesh.upper_right),
        energy_bounds=energy_bounds,
        source_points=source_points,
        cell_volume=cell_volume,
    )

    windows = openmc.WeightWindows(
        mesh=mesh,
        lower_ww_bounds=bounds.lower,
        upper_ww_bounds=bounds.upper,
        energy_bounds=energy_bounds,
        particle_type="neutron",
        survival_ratio=3.0,
        max_lower_bound_ratio=5.0,
        max_split=4,
        weight_cutoff=1.0e-38,
    )
    ww_path = Path("ww_reciprocity_cadis.h5").resolve()
    openmc.WeightWindowsList([windows]).export_to_hdf5(ww_path)
    restored = openmc.hdf5_to_wws(ww_path)
    if len(restored) != 1:
        raise RuntimeError(f"Expected one weight-window object, found {len(restored)}")

    reliable_relerr = bounds.relative_error[bounds.reliable]
    result = {
        "status": "PASS",
        "statepoint": str(Path(statepoint_path).resolve()),
        "importance": str(Path("reciprocal_importance.npz").resolve()),
        "weight_windows": str(ww_path),
        "weight_windows_bytes": ww_path.stat().st_size,
        "reliable_bins": int(np.count_nonzero(bounds.reliable)),
        "positive_sampled_bins": int(np.count_nonzero(mean > 0.0)),
        "total_bins": int(mean.size),
        "reliable_fraction": float(np.mean(bounds.reliable)),
        "median_reliable_relative_error": float(np.median(reliable_relerr)),
        "maximum_reliable_relative_error": float(np.max(reliable_relerr)),
    }
    Path("result.json").write_text(json.dumps(result, indent=2) + "\n")
    print("RESULT " + json.dumps(result, sort_keys=True), flush=True)
    print("DONE_RECIPROCITY_CADIS", flush=True)
    # This source-built Python/OpenMC combination can double-finalize the
    # preloaded library during interpreter teardown after TemporarySession has
    # already finalized it. All outputs are closed at this point; bypassing
    # Python teardown gives Slurm the correct success status.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
