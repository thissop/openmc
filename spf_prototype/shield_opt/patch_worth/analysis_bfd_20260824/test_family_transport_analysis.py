#!/usr/bin/env python3
"""Standard-library tests for family_transport_analysis.py."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from family_transport_analysis import analyze_results, load_result, render_markdown


def _payload(flux, relerr, histories=4_000_000, cells=None, labels=None):
    cells = cells or list(range(19, 19 + len(flux)))
    labels = labels or cells
    peak_index = max(range(len(flux)), key=flux.__getitem__)
    return {
        "histories": histories,
        "coil_cells": cells,
        "coil_labels": labels,
        "coil_fast_flux": flux,
        "coil_fast_flux_relerr": relerr,
        "objective_peak_coil_flux": flux[peak_index],
        "peak_coil_cell": cells[peak_index],
        "peak_coil_label": labels[peak_index],
    }


class FamilyTransportAnalysisTests(unittest.TestCase):
    def _load(self, payload):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "result.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            return load_result(path)

    def test_load_result_rejects_declared_peak_mismatch(self):
        payload = _payload([1.0, 2.0], [0.05, 0.05])
        payload["peak_coil_cell"] = 19
        with self.assertRaisesRegex(ValueError, "reported peak cell"):
            self._load(payload)

    def test_significant_improvement_and_hotspot_migration(self):
        uniform = self._load(_payload([8.0, 10.0, 3.0], [0.02, 0.02, 0.05]))
        improved = self._load(_payload([7.0, 5.0, 3.0], [0.02, 0.02, 0.05]))
        summary = analyze_results({"famU": uniform, "famC": improved})
        row = summary["families"][1]
        self.assertEqual(row["peak_cell"], 19)
        self.assertEqual(row["verdict"], "improves peak")
        self.assertLess(row["peak_change_z"], -2.0)
        self.assertTrue(row["hotspot_ready"])
        self.assertFalse(summary["complete"])

    def test_small_change_is_statistically_unresolved(self):
        uniform = self._load(_payload([5.0, 10.0], [0.05, 0.05]))
        candidate = self._load(_payload([5.0, 9.8], [0.05, 0.05]))
        summary = analyze_results({"famU": uniform, "famR": candidate})
        self.assertEqual(summary["families"][1]["verdict"], "statistically unresolved")

    def test_peak_margin_uses_combined_uncertainty(self):
        uniform = self._load(_payload([9.9, 10.0], [0.10, 0.10]))
        summary = analyze_results({"famU": uniform})
        row = summary["families"][0]
        expected = 0.1 / math.hypot(0.99, 1.0)
        self.assertAlmostEqual(row["peak_margin_z"], expected)
        self.assertFalse(row["hotspot_ready"])

    def test_mismatched_histories_fail_closed(self):
        uniform = self._load(_payload([5.0, 10.0], [0.05, 0.05]))
        candidate = self._load(
            _payload([4.0, 8.0], [0.05, 0.05], histories=2_000_000)
        )
        with self.assertRaisesRegex(ValueError, "histories differ"):
            analyze_results({"famU": uniform, "famM": candidate})

    def test_shifted_geometry_cells_use_stable_coil_labels(self):
        uniform = self._load(_payload([5.0, 10.0], [0.05, 0.05]))
        candidate = self._load(
            _payload(
                [4.0, 8.0],
                [0.05, 0.05],
                cells=[25, 26],
                labels=[19, 20],
            )
        )
        summary = analyze_results({"famU": uniform, "famM": candidate})
        row = summary["families"][1]
        self.assertEqual(row["peak_cell"], 20)
        self.assertEqual(row["tally_peak_cell"], 26)
        self.assertEqual(row["verdict"], "improves peak")

    def test_markdown_states_gate_status(self):
        uniform = self._load(_payload([5.0, 10.0], [0.05, 0.05]))
        summary = analyze_results({"famU": uniform})
        report = render_markdown(summary, 2_163_095.0)
        self.assertIn("Peak-value headline gate: **NOT READY**", report)
        self.assertIn("2.163×10⁶ kg", report)


if __name__ == "__main__":
    unittest.main()
