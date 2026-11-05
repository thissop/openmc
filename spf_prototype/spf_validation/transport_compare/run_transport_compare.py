#!/usr/bin/env python3
"""Run one matched-history grid-versus-spectral stellarator transport case.

The geometry, materials, tallies, source emission profile, source energy,
history count, and RNG seed are held fixed.  Only the StellaratorSource
representation changes.  The two samplers consume random numbers differently,
so a shared seed is reproducible but is not a per-particle common-random-number
coupling; uncertainty comparisons therefore treat the runs as independent.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from pathlib import Path

import numpy as np


POLARIZATION = {
    "unpol": None,
    "A": (1.0, 0.0, 0.0),
    "B": (0.0, 1.0, 0.0),
    "C": (0.0, 0.0, 1.0),
}


def parse_cell_spec(spec: str) -> list[int]:
    """Parse comma-separated cell IDs and inclusive ranges (for example 11-26)."""
    cells: list[int] = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            lo_text, hi_text = item.split("-", 1)
            lo, hi = int(lo_text), int(hi_text)
            if hi < lo:
                raise ValueError(f"descending cell range: {item}")
            cells.extend(range(lo, hi + 1))
        else:
            cells.append(int(item))
    if not cells:
        raise ValueError("at least one coil cell is required")
    if len(cells) != len(set(cells)):
        raise ValueError("coil cell list contains duplicates")
    return sorted(cells)


def cells_from_npz(path: Path) -> list[int]:
    """Read coil-cell IDs from the known QH/QA discovery-file schemas."""
    with np.load(path) as data:
        for key in ("magnet_cells", "coil_cells", "cell_ids", "cells"):
            if key in data:
                cells = [int(x) for x in np.asarray(data[key]).ravel()]
                if cells:
                    return sorted(cells)
    raise ValueError(f"no coil-cell key found in {path}")


def _build_material_bank(openmc, li6: float):
    """Return the material names used by the established QA/QH DAGMC builds."""
    w = openmc.Material()
    w.add_element("W", 1.0)
    w.set_density("g/cm3", 19.3)

    rafm = openmc.Material()
    rafm.add_element("Fe", 0.885, "wo")
    rafm.add_element("Cr", 0.09, "wo")
    rafm.add_element("W", 0.011, "wo")
    rafm.add_element("Mn", 0.006, "wo")
    rafm.add_element("V", 0.002, "wo")
    rafm.add_element("C", 0.006, "wo")
    rafm.set_density("g/cm3", 7.9)
    fw = openmc.Material.mix_materials(
        [w, rafm], [0.2 / 3.2, 3.0 / 3.2], "vo")
    fw.name = "first_wall"

    multiplier = openmc.Material(name="multiplier")
    multiplier.add_element("Be", 1.0)
    multiplier.set_density("g/cm3", 1.85)

    breeder = openmc.Material(name="breeder")
    breeder.add_element(
        "Li", 2.0, enrichment=li6, enrichment_target="Li6",
        enrichment_type="ao")
    breeder.add_element("Be", 1.0)
    breeder.add_element("F", 4.0)
    breeder.set_density("g/cm3", 1.94)

    back_wall = openmc.Material(name="back_wall")
    back_wall.add_element("Fe", 0.885, "wo")
    back_wall.add_element("Cr", 0.09, "wo")
    back_wall.add_element("W", 0.011, "wo")
    back_wall.add_element("Mn", 0.006, "wo")
    back_wall.add_element("V", 0.002, "wo")
    back_wall.add_element("C", 0.006, "wo")
    back_wall.set_density("g/cm3", 7.9)

    shield = openmc.Material(name="shield")
    shield.add_element("W", 1.0)
    shield.add_element("C", 1.0)
    shield.set_density("g/cm3", 15.6)

    gap = openmc.Material(name="gap")
    gap.add_nuclide("H1", 1.0)
    gap.set_density("g/cm3", 1e-8)

    vessel = openmc.Material(name="vacuum_vessel")
    vessel.add_element("Fe", 0.62, "wo")
    vessel.add_element("Cr", 0.16, "wo")
    vessel.add_element("Ni", 0.11, "wo")
    vessel.add_element("Mo", 0.02, "wo")
    vessel.add_element("Mn", 0.015, "wo")
    vessel.add_element("B", 0.03, "wo")
    vessel.add_element("O", 0.03, "wo")
    vessel.add_element("H", 0.005, "wo")
    vessel.set_density("g/cm3", 7.9)

    thermal = openmc.Material(name="thermal_shield")
    thermal.add_element("Fe", 0.89, "wo")
    thermal.add_element("Cr", 0.09, "wo")
    thermal.add_element("Mn", 0.02, "wo")
    thermal.set_density("g/cm3", 7.9)

    magnets = openmc.Material(name="magnets")
    magnets.add_element("Cu", 0.58, "wo")
    magnets.add_element("Fe", 0.30, "wo")
    magnets.add_element("Cr", 0.04, "wo")
    magnets.add_element("Ni", 0.03, "wo")
    magnets.add_element("Y", 0.02, "wo")
    magnets.add_element("Ba", 0.02, "wo")
    magnets.add_element("O", 0.01, "wo")
    magnets.set_density("g/cm3", 8.5)

    vacuum = openmc.Material(name="Vacuum")
    vacuum.add_nuclide("H1", 1.0)
    vacuum.set_density("g/cm3", 1e-10)

    return {m.name: m for m in (
        fw, multiplier, breeder, back_wall, shield, gap, vessel, thermal,
        magnets, vacuum)}


def _materials_for_geometry(openmc, dag, li6: float):
    bank = _build_material_bank(openmc, li6)
    tags = list(dag.material_names)
    missing = sorted(set(tags) - set(bank))
    if missing:
        raise ValueError(f"DAGMC contains unsupported material tags: {missing}")
    selected = [bank[tag] for tag in tags]
    return openmc.Materials(selected), bank


def _make_source(openmc, args):
    energy = openmc.stats.Discrete([14.06e6], [1.0])
    polarization = POLARIZATION[args.mode]
    emission_rho = np.linspace(0.0, 1.0, 257)
    emission_density = 1.0 - emission_rho**2

    if args.representation == "grid":
        if args.fluxmap is None:
            raise ValueError("--fluxmap is required for the grid representation")
        return openmc.StellaratorSource(
            fluxmap=str(args.fluxmap),
            polarization=polarization,
            field_model="fluxmap",
            emission=(emission_rho, emission_density),
            energy=energy,
            strength=1.0,
        )

    if args.wout is None:
        raise ValueError("--wout is required for the spectral representation")
    source = openmc.StellaratorSource.from_vmec(
        args.wout,
        emission_density=lambda rho: 1.0 - np.asarray(rho)**2,
        energy=energy,
        polarization=polarization,
        strength=1.0,
    )
    if args.length_scale != 1.0:
        source.rmnc = source.rmnc * args.length_scale
        source.zmns = source.zmns * args.length_scale
        if source.rmns is not None:
            source.rmns = source.rmns * args.length_scale
            source.zmnc = source.zmnc * args.length_scale
    return source


def _build_model(openmc, args, cells):
    openmc.reset_auto_ids()
    dag = openmc.DAGMCUniverse(str(args.dagmc), auto_geom_ids=True)
    materials, bank = _materials_for_geometry(openmc, dag, args.li6)

    box = dag.bounding_box
    lower = np.asarray(box.lower_left)
    upper = np.asarray(box.upper_right)
    extent = float(np.max(np.abs(np.concatenate([lower, upper]))))
    boundary = openmc.Sphere(r=1.5 * extent + 10.0, boundary_type="vacuum")
    geometry = openmc.Geometry([
        openmc.Cell(region=-boundary, fill=dag),
    ])

    fast = openmc.EnergyFilter([0.1e6, 25.0e6])
    coil_filter = openmc.CellFilter(cells)
    coil_flux = openmc.Tally(name="coil_fast_flux")
    coil_flux.filters = [coil_filter, fast]
    coil_flux.scores = ["flux"]
    coil_heating = openmc.Tally(name="coil_heating")
    coil_heating.filters = [coil_filter]
    coil_heating.scores = ["heating-local"]

    wall_flux = openmc.Tally(name="wall_fast_flux")
    wall_flux.filters = [openmc.MaterialFilter([bank["first_wall"]]), fast]
    wall_flux.scores = ["flux"]
    tbr = openmc.Tally(name="tbr")
    tbr.filters = [openmc.MaterialFilter([bank["breeder"]])]
    tbr.scores = ["H3-production"]

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.source = [_make_source(openmc, args)]
    settings.batches = args.batches
    settings.particles = args.particles
    settings.inactive = 0
    settings.seed = args.seed
    settings.photon_transport = False
    settings.max_lost_particles = args.max_lost
    settings.rel_max_lost_particles = args.rel_max_lost
    settings.output = {"summary": True, "tallies": False}

    tallies = openmc.Tallies([coil_flux, coil_heating, wall_flux, tbr])
    return openmc.Model(
        geometry=geometry, materials=materials, settings=settings,
        tallies=tallies), list(dag.material_names)


def _tally_arrays(statepoint, name):
    tally = statepoint.get_tally(name=name)
    return tally.mean.ravel(), tally.std_dev.ravel()


def run(args):
    import openmc

    cells = (cells_from_npz(args.cells_npz) if args.cells_npz is not None
             else parse_cell_spec(args.coil_cells))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    model, material_tags = _build_model(openmc, args, cells)

    manifest = {
        "device": args.device,
        "representation": args.representation,
        "mode": args.mode,
        "dagmc": str(args.dagmc),
        "fluxmap": str(args.fluxmap) if args.fluxmap else None,
        "wout": str(args.wout) if args.wout else None,
        "length_scale": args.length_scale,
        "coil_cells": cells,
        "material_tags": material_tags,
        "batches": args.batches,
        "particles": args.particles,
        "total_histories": args.batches * args.particles,
        "seed": args.seed,
        "li6_percent": args.li6,
        "openmc_version": getattr(openmc, "__version__", "unknown"),
        "uncertainty_note": (
            "shared seed and history count; samplers consume RNG differently, "
            "so comparisons use independent-run standard errors"),
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")

    started = time.time()
    statepoint_path = model.run(
        cwd=args.output_dir, output=True, threads=args.threads)
    elapsed = time.time() - started
    statepoint_copy = args.output_dir / "statepoint.h5"
    if Path(statepoint_path).resolve() != statepoint_copy.resolve():
        shutil.copy2(statepoint_path, statepoint_copy)

    with openmc.StatePoint(statepoint_copy) as statepoint:
        coil_flux, coil_flux_sd = _tally_arrays(statepoint, "coil_fast_flux")
        coil_heat, coil_heat_sd = _tally_arrays(statepoint, "coil_heating")
        wall_flux, wall_flux_sd = _tally_arrays(statepoint, "wall_fast_flux")
        tbr, tbr_sd = _tally_arrays(statepoint, "tbr")

    output = args.output_dir / "result.npz"
    np.savez(
        output,
        device=args.device,
        representation=args.representation,
        mode=args.mode,
        coil_cells=np.asarray(cells),
        coil_fast_flux=coil_flux,
        coil_fast_flux_sd=coil_flux_sd,
        coil_heating=coil_heat,
        coil_heating_sd=coil_heat_sd,
        wall_fast_flux=wall_flux,
        wall_fast_flux_sd=wall_flux_sd,
        tbr=tbr,
        tbr_sd=tbr_sd,
        batches=args.batches,
        particles=args.particles,
        total_histories=args.batches * args.particles,
        seed=args.seed,
        elapsed_seconds=elapsed,
    )
    manifest["elapsed_seconds"] = elapsed
    manifest["result"] = str(output)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", required=True, choices=("qa", "qh"))
    parser.add_argument(
        "--representation", required=True, choices=("grid", "spectral"))
    parser.add_argument("--mode", default="A", choices=tuple(POLARIZATION))
    parser.add_argument("--dagmc", required=True, type=Path)
    cells = parser.add_mutually_exclusive_group(required=True)
    cells.add_argument("--cells-npz", type=Path)
    cells.add_argument("--coil-cells")
    parser.add_argument("--fluxmap", type=Path)
    parser.add_argument("--wout", type=Path)
    parser.add_argument(
        "--length-scale", type=float, default=1.0,
        help="additional scale applied to spectral R/Z after VMEC m-to-cm conversion")
    parser.add_argument("--batches", type=int, default=4)
    parser.add_argument("--particles", type=int, default=25_000)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS", 32)))
    parser.add_argument("--li6", type=float, default=60.0)
    parser.add_argument("--max-lost", type=int, default=100)
    parser.add_argument("--rel-max-lost", type=float, default=1e-4)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.batches < 1 or args.particles < 1 or args.threads < 1:
        parser.error("batches, particles, and threads must be positive")
    if args.length_scale <= 0.0:
        parser.error("length-scale must be positive")
    return args


if __name__ == "__main__":
    run(parse_args())
