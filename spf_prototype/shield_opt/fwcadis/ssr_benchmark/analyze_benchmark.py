#!/usr/bin/env python
"""Aggregate SSR benchmark JSON files and compare them with the master run."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from metrics import (
    median_positive_finite,
    normalized_pattern,
    rms_fractional_bias,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="+", help="result.json files or directories")
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--surface-source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def find_results(paths: list[Path]) -> list[Path]:
    found = []
    for path in paths:
        if path.is_dir():
            found.extend(path.rglob("result.json"))
        elif path.name == "result.json":
            found.append(path)
        else:
            raise ValueError(f"not a result.json file or directory: {path}")
    unique = sorted({path.resolve() for path in found})
    if not unique:
        raise ValueError("no result.json files found")
    return unique


def source_bank_metadata(path: Path, chunk_size: int = 1_000_000) -> dict:
    import h5py

    with h5py.File(path, "r") as source_file:
        bank = source_file["source_bank"]
        weight_sum = 0.0
        for start in range(0, len(bank), chunk_size):
            weight_sum += float(np.sum(bank[start : start + chunk_size]["wgt"]))
        return {"sites": len(bank), "weight_sum": weight_sum}


def compatible_bins(
    estimate: np.ndarray,
    estimate_error: np.ndarray,
    reference: np.ndarray,
    reference_error: np.ndarray,
    sigma: float = 2.0,
) -> int:
    uncertainty = np.sqrt(estimate_error**2 + reference_error**2)
    return int(np.count_nonzero(np.abs(estimate - reference) <= sigma * uncertainty))


def finite_or_none(value: float) -> float | None:
    return float(value) if math.isfinite(value) else None


def summarize(result: dict, master: dict, source_factor: float) -> dict:
    raw = np.asarray(result["coil_fast_flux"]["mean"], dtype=float)
    raw_error = np.asarray(result["coil_fast_flux"]["std_dev"], dtype=float)
    reference = np.asarray(master["coil_fast_flux"], dtype=float)
    reference_error = np.asarray(master["coil_fast_flux_err"], dtype=float)
    scaled = raw * source_factor
    scaled_error = raw_error * source_factor
    pattern = normalized_pattern(scaled.tolist())
    reference_pattern = normalized_pattern(reference.tolist())
    pattern_rmse = (
        None
        if any(math.isnan(value) for value in pattern)
        else rms_fractional_bias(pattern, reference_pattern)
    )
    fractional_bias = (scaled - reference) / reference
    relative_error = result["coil_fast_flux"]["relative_error"]
    foms = result["coil_fast_flux"]["figure_of_merit"]
    return {
        "mode": result["mode"],
        "result_path": result.get("workdir"),
        "source_histories": result["configuration"]["source_histories"],
        "simulation_seconds": result["runtime"]["simulation"],
        "runner_wall_seconds": result["runner_wall_seconds"],
        "estimated_tracks_per_source_history": result[
            "estimated_tracks_per_source_history"
        ],
        "positive_coil_bins": int(np.count_nonzero(raw > 0.0)),
        "median_positive_relative_error": finite_or_none(
            median_positive_finite(relative_error)
        ),
        "median_positive_fom": median_positive_finite(foms),
        "source_normalized_coil_flux": scaled.tolist(),
        "source_normalized_coil_std_dev": scaled_error.tolist(),
        "fractional_bias_from_master": fractional_bias.tolist(),
        "rms_fractional_bias_from_master": rms_fractional_bias(
            scaled.tolist(), reference.tolist()
        ),
        "rms_fractional_pattern_error": pattern_rmse,
        "bins_compatible_with_master_at_2sigma": compatible_bins(
            scaled, scaled_error, reference, reference_error
        ),
        "weight_windows_sha256": (
            None
            if result["inputs"]["weight_windows"] is None
            else result["inputs"]["weight_windows"]["sha256"]
        ),
    }


def fmt(value, precision: int = 3) -> str:
    if value is None or not math.isfinite(float(value)):
        return "—"
    return f"{float(value):.{precision}g}"


def markdown(summary: dict) -> str:
    lines = [
        "# SSR variance-reduction benchmark",
        "",
        (
            "Surface-source responses are multiplied by the bank's total starting "
            "weight per master source history. Bias therefore tests the full "
            "deep-freeze replay, including the deliberate black-inside approximation."
        ),
        "",
        f"- Surface-bank sites: {summary['source_bank']['sites']:,}",
        f"- Surface-bank total weight: {summary['source_bank']['weight_sum']:.8g}",
        (
            "- Surface weight per master history: "
            f"{summary['source_weight_per_master_history']:.8g}"
        ),
        "",
        "| mode | coil bins > 0 | median rel. err. | median FOM | transport (s) | tracks/source | RMS absolute bias | RMS pattern error | master-compatible bins |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in summary["cases"]:
        lines.append(
            "| {mode} | {positive}/5 | {rel} | {fom} | {seconds} | {tracks} | {bias} | {pattern} | {compatible}/5 |".format(
                mode=case["mode"],
                positive=case["positive_coil_bins"],
                rel=fmt(case["median_positive_relative_error"]),
                fom=fmt(case["median_positive_fom"]),
                seconds=fmt(case["simulation_seconds"]),
                tracks=fmt(case["estimated_tracks_per_source_history"]),
                bias=fmt(case["rms_fractional_bias_from_master"]),
                pattern=fmt(case["rms_fractional_pattern_error"]),
                compatible=case["bins_compatible_with_master_at_2sigma"],
            )
        )
    lines.extend([
        "",
        "The pilot is a split-safety and wiring check. A production verdict requires all five coil bins to be positive and stable across independent seeds.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    result_paths = find_results(args.results)
    loaded = [json.loads(path.read_text(encoding="utf-8")) for path in result_paths]
    failures = [result for result in loaded if result.get("status") != "PASS"]
    passed = [result for result in loaded if result.get("status") == "PASS"]
    if not passed:
        raise RuntimeError("no successful benchmark results")

    master_npz = np.load(args.master, allow_pickle=False)
    master = {key: master_npz[key].tolist() for key in master_npz.files}
    source_bank = source_bank_metadata(args.surface_source)
    master_histories = int(master["n_source"])
    source_factor = source_bank["weight_sum"] / master_histories
    cases = [summarize(result, master, source_factor) for result in passed]
    mode_order = {mode: index for index, mode in enumerate(("analog", "magic", "native", "reciprocity"))}
    cases.sort(key=lambda case: (mode_order.get(case["mode"], 99), case["result_path"]))

    summary = {
        "master": str(args.master.resolve()),
        "surface_source": str(args.surface_source.resolve()),
        "master_histories": master_histories,
        "source_bank": source_bank,
        "source_weight_per_master_history": source_factor,
        "cases": cases,
        "failures": failures,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "benchmark_summary.json"
    markdown_path = args.output_dir / "BENCHMARK_RESULTS.md"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    markdown_path.write_text(markdown(summary), encoding="utf-8")
    print(markdown(summary))
    print(f"WROTE {json_path}")
    print(f"WROTE {markdown_path}")


if __name__ == "__main__":
    main()

