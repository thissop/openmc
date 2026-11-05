#!/usr/bin/env python3
"""Analyze matched-mass shield-family transport results with uncertainty gates.

The cluster evaluator writes one JSON file per family.  This module keeps the
comparison deliberately small and auditable: it verifies the reported peak,
propagates independent Monte Carlo uncertainties for ratios and differences,
and reports hotspot migration only when the winning coil is separated from the
runner-up by at least the requested number of standard deviations.

No OpenMC or NumPy dependency is required, so the report can be regenerated on
the login node or any stock Python installation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


FAMILY_ORDER = ("famU", "famC", "famR", "famM")
FAMILY_LABELS = {
    "famU": "uniform",
    "famC": "clearance-driven",
    "famR": "response-aware",
    "famM": "linear-minimax stress test",
}


def _as_float_list(value: Any, field: str) -> list[float]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{field} must be a non-empty JSON list")
    result = [float(item) for item in value]
    if not all(math.isfinite(item) for item in result):
        raise ValueError(f"{field} contains a non-finite value")
    return result


def _as_int_list(value: Any, field: str) -> list[int]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{field} must be a non-empty JSON list")
    return [int(item) for item in value]


def load_result(path: str | Path) -> dict[str, Any]:
    """Load and validate one ``forward_eval.py`` result JSON."""

    result_path = Path(path)
    payload = result_path.read_bytes()
    raw = json.loads(payload)

    cells = _as_int_list(raw.get("coil_cells"), "coil_cells")
    labels = _as_int_list(raw.get("coil_labels", cells), "coil_labels")
    flux = _as_float_list(raw.get("coil_fast_flux"), "coil_fast_flux")
    relerr = _as_float_list(
        raw.get("coil_fast_flux_relerr"), "coil_fast_flux_relerr"
    )
    if not (len(cells) == len(labels) == len(flux) == len(relerr)):
        raise ValueError(
            "coil_cells, coil_labels, coil_fast_flux, and relerr lengths differ"
        )
    if len(set(cells)) != len(cells):
        raise ValueError("coil_cells contains duplicates")
    if len(set(labels)) != len(labels):
        raise ValueError("coil_labels contains duplicates")
    if any(item < 0 for item in flux):
        raise ValueError("coil_fast_flux contains a negative value")
    if any(item < 0 for item in relerr):
        raise ValueError("coil_fast_flux_relerr contains a negative value")

    peak_index = max(range(len(flux)), key=flux.__getitem__)
    peak_flux = flux[peak_index]
    peak_cell = cells[peak_index]
    peak_label = labels[peak_index]
    reported_flux = float(raw.get("objective_peak_coil_flux"))
    reported_cell = int(raw.get("peak_coil_cell"))
    if not math.isclose(reported_flux, peak_flux, rel_tol=1e-12, abs_tol=0.0):
        raise ValueError(
            f"reported peak {reported_flux:.8e} != coil-vector max {peak_flux:.8e}"
        )
    if reported_cell != peak_cell:
        raise ValueError(f"reported peak cell {reported_cell} != vector peak {peak_cell}")
    if "peak_coil_label" in raw and int(raw["peak_coil_label"]) != peak_label:
        raise ValueError(
            f"reported peak label {raw['peak_coil_label']} != vector label {peak_label}"
        )

    histories = int(raw.get("histories"))
    if histories <= 0:
        raise ValueError("histories must be positive")

    return {
        "path": str(result_path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "histories": histories,
        "coil_cells": cells,
        "coil_labels": labels,
        "coil_fast_flux": flux,
        "coil_fast_flux_relerr": relerr,
        "peak_index": peak_index,
        "peak_cell": peak_cell,
        "peak_label": peak_label,
        "peak_flux": peak_flux,
        "peak_relerr": relerr[peak_index],
    }


def _peak_margin_z(result: Mapping[str, Any]) -> float:
    """Return separation of the two largest coil estimates in combined sigma."""

    flux = result["coil_fast_flux"]
    relerr = result["coil_fast_flux_relerr"]
    order = sorted(range(len(flux)), key=flux.__getitem__, reverse=True)
    if len(order) < 2:
        return math.inf
    first, second = order[:2]
    denominator = math.hypot(flux[first] * relerr[first], flux[second] * relerr[second])
    return math.inf if denominator == 0 else (flux[first] - flux[second]) / denominator


def analyze_results(
    results: Mapping[str, Mapping[str, Any]],
    *,
    sigma_threshold: float = 2.0,
    peak_relerr_threshold: float = 0.10,
) -> dict[str, Any]:
    """Compare all families with ``famU`` and return a serializable summary."""

    if "famU" not in results:
        raise ValueError("famU is required as the uniform reference")
    if sigma_threshold <= 0:
        raise ValueError("sigma_threshold must be positive")
    if not 0 < peak_relerr_threshold < 1:
        raise ValueError("peak_relerr_threshold must lie between 0 and 1")

    baseline = results["famU"]
    reference_labels = baseline["coil_labels"]
    reference_histories = baseline["histories"]
    base_peak = baseline["peak_flux"]
    base_sigma = base_peak * baseline["peak_relerr"]

    rows: list[dict[str, Any]] = []
    for name in FAMILY_ORDER:
        if name not in results:
            continue
        result = results[name]
        if result["coil_labels"] != reference_labels:
            raise ValueError(f"{name}: coil_labels differ from famU")
        if result["histories"] != reference_histories:
            raise ValueError(f"{name}: histories differ from famU")

        peak = result["peak_flux"]
        peak_sigma = peak * result["peak_relerr"]
        ratio = peak / base_peak
        ratio_sigma = ratio * math.hypot(
            result["peak_relerr"], baseline["peak_relerr"]
        )
        delta_sigma = math.hypot(peak_sigma, base_sigma)
        delta_z = 0.0 if name == "famU" else (peak - base_peak) / delta_sigma
        margin_z = _peak_margin_z(result)

        peak_value_ready = result["peak_relerr"] <= peak_relerr_threshold
        hotspot_ready = margin_z >= sigma_threshold
        if name == "famU":
            verdict = "reference"
        elif abs(delta_z) < sigma_threshold:
            verdict = "statistically unresolved"
        elif peak < base_peak:
            verdict = "improves peak"
        else:
            verdict = "worsens peak"

        rows.append(
            {
                "family": name,
                "label": FAMILY_LABELS.get(name, name),
                "histories": result["histories"],
                "peak_cell": result["peak_label"],
                "tally_peak_cell": result["peak_cell"],
                "peak_flux": peak,
                "peak_relerr": result["peak_relerr"],
                "peak_ratio_vs_uniform": ratio,
                "peak_ratio_sigma": ratio_sigma,
                "peak_change_percent": 100.0 * (ratio - 1.0),
                "peak_change_sigma_percent": 100.0 * ratio_sigma,
                "peak_change_z": delta_z,
                "resolved_coils_relerr_lt_threshold": sum(
                    item < peak_relerr_threshold
                    for item in result["coil_fast_flux_relerr"]
                ),
                "peak_margin_z": margin_z,
                "peak_value_ready": peak_value_ready,
                "hotspot_ready": hotspot_ready,
                "verdict": verdict,
            }
        )

    complete = all(name in results for name in FAMILY_ORDER)
    headline_ready = complete and all(row["peak_value_ready"] for row in rows)
    hotspot_labels_ready = complete and all(row["hotspot_ready"] for row in rows)
    return {
        "reference_family": "famU",
        "sigma_threshold": sigma_threshold,
        "peak_relerr_threshold": peak_relerr_threshold,
        "complete": complete,
        "headline_ready": headline_ready,
        "hotspot_labels_ready": hotspot_labels_ready,
        "families": rows,
    }


def render_markdown(summary: Mapping[str, Any], common_mass_kg: float) -> str:
    """Render the summary as a compact, meeting-ready Markdown report."""

    status = "READY" if summary["headline_ready"] else "NOT READY"
    hotspot_status = "READY" if summary["hotspot_labels_ready"] else "NOT READY"
    lines = [
        "# Equal-mass shield-family transport",
        "",
        f"Common design mass: **{common_mass_kg / 1e6:.3f}×10⁶ kg**. "
        f"Peak-value headline gate: **{status}**. Hotspot-label gate: **{hotspot_status}**.",
        "",
        "| family | peak coil | peak fast flux | rel. err. | Δ vs uniform | z | verdict |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary["families"]:
        lines.append(
            "| {label} | c{peak_cell} | {peak_flux:.4e} | {peak_relerr:.1%} | "
            "{peak_change_percent:+.1f}% ± {peak_change_sigma_percent:.1f}% | "
            "{peak_change_z:+.2f} | {verdict} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## Integrity gates",
            "",
            f"- A peak value is headline-ready only when its relative error is ≤ "
            f"{summary['peak_relerr_threshold']:.0%}.",
            f"- Improvement or degradation requires |z| ≥ {summary['sigma_threshold']:.1f} "
            "after independent-error propagation.",
            f"- A hotspot label is treated as resolved only when the top coil exceeds the "
            f"runner-up by ≥ {summary['sigma_threshold']:.1f} combined σ.",
            "- Coil identity is compared by stable physical labels; each source JSON also "
            "retains the actual DAGMC tally cell IDs when split solids shift geometry IDs.",
            "- The comparison concerns transported peak coil fast flux at matched design mass; "
            "it does not rehabilitate the linear surrogate used to propose the families.",
            "",
        ]
    )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing result_famU/C/R/M.json",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-markdown", type=Path, required=True)
    parser.add_argument("--common-mass-kg", type=float, required=True)
    parser.add_argument("--sigma-threshold", type=float, default=2.0)
    parser.add_argument("--peak-relerr-threshold", type=float, default=0.10)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    results: dict[str, dict[str, Any]] = {}
    for name in FAMILY_ORDER:
        path = args.results_dir / f"result_{name}.json"
        if path.exists():
            results[name] = load_result(path)

    summary = analyze_results(
        results,
        sigma_threshold=args.sigma_threshold,
        peak_relerr_threshold=args.peak_relerr_threshold,
    )
    summary["common_mass_kg"] = args.common_mass_kg
    summary["source_files"] = {
        name: {"path": result["path"], "sha256": result["sha256"]}
        for name, result in results.items()
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    args.output_markdown.write_text(
        render_markdown(summary, args.common_mass_kg), encoding="utf-8"
    )
    print(
        f"families={len(results)}/{len(FAMILY_ORDER)} "
        f"headline_ready={summary['headline_ready']} "
        f"hotspot_labels_ready={summary['hotspot_labels_ready']}"
    )
    print(f"WROTE {args.output_json}")
    print(f"WROTE {args.output_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
