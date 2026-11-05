#!/usr/bin/env python3
"""Standard-library unit tests for the reciprocity-CADIS helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from cadis import build_cadis_bounds
from common import read_makegrid, sample_coil_points, stellarator_source_point_cloud


class CadisBoundsTests(unittest.TestCase):
    def test_inverts_volume_normalized_importance(self):
        mean = np.array([1.0, 2.0, 4.0])
        result = build_cadis_bounds(mean, np.full(3, 0.1), cell_volume=2.0)
        np.testing.assert_allclose(result.lower, [0.5, 0.25, 0.125])
        np.testing.assert_allclose(result.upper, 5.0 * result.lower)

    def test_ignores_zero_and_unreliable_bins_before_normalization(self):
        mean = np.array([0.0, 1.0e-9, 2.0])
        std_dev = np.array([0.0, 1.0e-8, 0.2])
        result = build_cadis_bounds(mean, std_dev, cell_volume=1.0)
        np.testing.assert_array_equal(result.reliable, [False, False, True])
        np.testing.assert_allclose(result.lower, [-1.0, -1.0, 0.5])
        np.testing.assert_allclose(result.upper, [-5.0, -5.0, 2.5])


class CoilSamplingTests(unittest.TestCase):
    def test_read_and_sample_makegrid(self):
        content = """periods 1\nmirror NIL\n0 0 0 1\n1 0 0 1 end\nend\n"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "coils")
            path.write_text(content)
            coils = read_makegrid(path)
        points = sample_coil_points(coils, points_per_coil=2, scale=1.0)
        np.testing.assert_allclose(points, [[0.25, 0.0, 0.0], [0.75, 0.0, 0.0]])

    def test_fluxmap_cloud_folds_into_one_period(self):
        nr, ntheta, nzeta = 2, 2, 2
        rho = np.array([0.0, 0.5])
        theta = np.array([0.0, np.pi])
        zeta = np.array([0.0, np.pi])
        shape = (nr, ntheta, nzeta)
        sqrtg = np.ones(shape)
        radius = np.full(shape, 10.0)
        vertical = np.zeros(shape)
        phi = np.linspace(0.1, 2.0 * np.pi - 0.1, np.prod(shape)).reshape(shape)
        zeros = np.zeros(shape)
        with tempfile.TemporaryDirectory() as directory:
            stem = Path(directory, "flux")
            stem.with_suffix(".meta").write_text(
                "magic spf_fluxmap_v1\nnR 2\nntheta 2\nnzeta 2\nnfp 4\n"
            )
            np.concatenate([
                rho, theta, zeta, sqrtg.ravel(), radius.ravel(), vertical.ravel(),
                phi.ravel(), zeros.ravel(), zeros.ravel(), zeros.ravel(),
            ]).astype("<f8").tofile(stem.with_suffix(".bin"))
            points, weights = stellarator_source_point_cloud(stem, max_points=4)
        self.assertEqual(points.shape, (4, 3))
        self.assertEqual(weights.shape, (4,))
        angles = np.arctan2(points[:, 1], points[:, 0])
        self.assertTrue(np.all((angles > 0.0) & (angles < np.pi / 2.0)))


if __name__ == "__main__":
    unittest.main()
