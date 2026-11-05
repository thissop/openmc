#!/usr/bin/env python3
"""Compare matched-history grid and spectral transport result files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


SCALAR_METRICS = ("wall_fast_flux", "tbr")
VECTOR_METRICS = ("coil_fast_flux", "coil_heating")


def independent_difference(grid, grid_sd, spectral, spectral_sd):
    """Return spectral-grid difference, independent-run SE, and z score."""
    difference = np.asarray(spectral) - np.asarray(grid)
    standard_error = np.hypot(np.asarray(grid_sd), np.asarray(spectral_sd))
    z_score = np.divide(
        difference, standard_error,
        out=np.full_like(difference, np.nan, dtype=float),
        where=standard_error > 0.0,
    )
    return difference, standard_error, z_score


def _scalar(value):
    array = np.asarray(value).ravel()
    if array.size != 1:
        raise ValueError(f"expected scalar result, got shape {np.asarray(value).shape}")
    return float(array[0])


def _relative(difference, reference):
    return np.divide(
        difference, reference,
        out=np.full_like(np.asarray(difference, dtype=float), np.nan),
        where=np.asarray(reference) != 0.0,
    )


def compare_pair(grid_path: Path, spectral_path: Path):
    with np.load(grid_path) as grid_file, np.load(spectral_path) as spectral_file:
        grid = {key: np.asarray(grid_file[key]) for key in grid_file.files}
        spectral = {key: np.asarray(spectral_file[key]) for key in spectral_file.files}

    identity_keys = ("device", "mode", "seed", "batches", "particles", "coil_cells")
    for key in identity_keys:
        if not np.array_equal(grid[key], spectral[key]):
            raise ValueError(f"grid/spectral mismatch in {key}")

    result = {
        "device": str(grid["device"]),
        "mode": str(grid["mode"]),
        "seed": int(grid["seed"]),
        "histories": int(grid["total_histories"]),
        "grid_result": str(grid_path),
        "spectral_result": str(spectral_path),
        "metrics": {},
    }
    for name in SCALAR_METRICS:
        g, gs = _scalar(grid[name]), _scalar(grid[f"{name}_sd"])
        s, ss = _scalar(spectral[name]), _scalar(spectral[f"{name}_sd"])
        delta, se, z = independent_difference(g, gs, s, ss)
        result["metrics"][name] = {
            "grid": g,
            "grid_sd": gs,
            "spectral": s,
            "spectral_sd": ss,
            "spectral_minus_grid": float(delta),
            "relative_to_grid": float(_relative(delta, g)),
            "difference_se": float(se),
            "z": float(z),
        }

    cells = grid["coil_cells"].astype(int)
    for name in VECTOR_METRICS:
        g, gs = grid[name].astype(float), grid[f"{name}_sd"].astype(float)
        s, ss = spectral[name].astype(float), spectral[f"{name}_sd"].astype(float)
        delta, se, z = independent_difference(g, gs, s, ss)
        rel = _relative(delta, g)
        max_z_index = int(np.nanargmax(np.abs(z))) if np.any(np.isfinite(z)) else 0
        result["metrics"][name] = {
            "cells": cells.tolist(),
            "grid": g.tolist(),
            "grid_sd": gs.tolist(),
            "spectral": s.tolist(),
            "spectral_sd": ss.tolist(),
            "spectral_minus_grid": delta.tolist(),
            "relative_to_grid": rel.tolist(),
            "difference_se": se.tolist(),
            "z": z.tolist(),
            "max_abs_z": float(np.nanmax(np.abs(z))),
            "max_abs_z_cell": int(cells[max_z_index]),
            "grid_peaking": float(np.max(g) / np.mean(g)),
            "spectral_peaking": float(np.max(s) / np.mean(s)),
            "grid_hottest_cell": int(cells[int(np.argmax(g))]),
            "spectral_hottest_cell": int(cells[int(np.argmax(s))]),
        }
    return result


def render_markdown(comparisons):
    lines = [
        "# Grid-versus-spectral transport comparison",
        "",
        "The runs use the same geometry, materials, tallies, emission profile, energy,",
        "history count, and seed. The samplers consume RNG differently, so the quoted",
        "difference z scores conservatively use independent-run standard errors.",
        "",
        "| Device | Mode | Histories | Wall flux rel. diff | Wall z | TBR rel. diff | TBR z | Coil-flux max | Hottest cell (grid/spec) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for comparison in comparisons:
        metrics = comparison["metrics"]
        wall = metrics["wall_fast_flux"]
        tbr = metrics["tbr"]
        coil = metrics["coil_fast_flux"]
        lines.append(
            f"| {comparison['device'].upper()} | {comparison['mode']} | "
            f"{comparison['histories']:,} | {wall['relative_to_grid']:+.2%} | "
            f"{wall['z']:+.2f} | {tbr['relative_to_grid']:+.2%} | {tbr['z']:+.2f} | "
            f"|z|={coil['max_abs_z']:.2f} (cell {coil['max_abs_z_cell']}) | "
            f"{coil['grid_hottest_cell']}/{coil['spectral_hottest_cell']} |")
    lines.extend([
        "",
        "A production claim should require both adequate tally precision and stability",
        "under an independent seed; a smoke result is only a pipeline gate.",
        "",
    ])
    return "\n".join(lines)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pairs", nargs="+", metavar="GRID,SPECTRAL",
        help="comma-separated result.npz pair")
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    comparisons = []
    for pair in args.pairs:
        try:
            grid, spectral = pair.split(",", 1)
        except ValueError as exc:
            raise SystemExit(f"invalid pair {pair!r}; expected GRID,SPECTRAL") from exc
        comparisons.append(compare_pair(Path(grid), Path(spectral)))
    args.json.write_text(json.dumps(comparisons, indent=2) + "\n")
    args.markdown.write_text(render_markdown(comparisons))


if __name__ == "__main__":
    main()
