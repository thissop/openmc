"""Tests for descriptor_correlation.py -- Chatterjee's xi and the screening matrix.

Key behaviours pinned:
  * xi(y=x) ~ 1                          (Y perfectly a function of X)
  * xi(independent x, y) ~ 0            (over several seeds; xi has small +/- scatter)
  * xi detects a nonmonotonic y=(x-.5)^2 (HIGH) where Pearson ~ 0  (the point of xi)
  * tie-free and tie-corrected forms agree when there are no ties
  * xi handles ties without blowing up
  * correlation_matrix has the right shape/orientation and recovers planted structure
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from descriptor_correlation import (  # noqa: E402
    chatterjee_xi,
    correlation_matrix,
    make_synthetic_dataset,
)


def test_perfect_functional_dependence_is_one():
    x = np.linspace(0, 1, 500)
    # Any strictly monotone Y of X should give xi -> 1.
    assert chatterjee_xi(x, x) == pytest.approx(1.0, abs=0.02)
    assert chatterjee_xi(x, np.exp(3 * x)) == pytest.approx(1.0, abs=0.02)


def test_independent_is_near_zero_over_seeds():
    vals = []
    for s in range(6):
        rng = np.random.default_rng(s)
        x = rng.random(2000)
        y = rng.random(2000)
        vals.append(chatterjee_xi(x, y, seed=s))
    vals = np.array(vals)
    # Independent -> xi centered on 0 with O(1/sqrt(n)) scatter.
    assert np.abs(np.mean(vals)) < 0.05
    assert np.all(np.abs(vals) < 0.12)


def test_detects_nonmonotonic_where_pearson_fails():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1, 4000)
    y = (x - 0.5) ** 2  # symmetric U-shape: Pearson ~ 0, but Y IS a function of X
    pearson = np.corrcoef(x, y)[0, 1]
    xi = chatterjee_xi(x, y, seed=0)
    assert abs(pearson) < 0.1, f"Pearson should be ~0, got {pearson}"
    assert xi > 0.7, f"xi should detect the functional U-shape, got {xi}"


def test_noisy_nonmonotonic_still_beats_pearson():
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1, 3000)
    y = (x - 0.5) ** 2 + rng.normal(0, 0.01, x.size)
    pearson = abs(np.corrcoef(x, y)[0, 1])
    xi = chatterjee_xi(x, y, seed=1)
    assert xi > pearson
    assert xi > 0.5


def test_tie_free_and_corrected_agree_without_ties():
    rng = np.random.default_rng(3)
    x = rng.random(1000)
    y = np.sin(6 * x) + 0.1 * rng.random(1000)  # continuous, no ties
    a = chatterjee_xi(x, y, ties=True, seed=7)
    b = chatterjee_xi(x, y, ties=False, seed=7)
    assert a == pytest.approx(b, abs=1e-9)


def test_ties_do_not_break():
    rng = np.random.default_rng(4)
    # Heavy ties in both X (integer descriptor) and Y (quantized response).
    x = rng.integers(0, 3, 400).astype(float)
    y = rng.integers(0, 4, 400).astype(float)
    xi = chatterjee_xi(x, y, ties=True, seed=0)
    assert np.isfinite(xi)
    assert -0.6 <= xi <= 1.0
    # Functional-through-ties: y = f(x) exactly, even with repeats -> high xi.
    y2 = (x == 2).astype(float)
    assert chatterjee_xi(x, y2, ties=True, seed=0) > 0.9


def test_constant_response_is_nan():
    x = np.arange(50.0)
    y = np.full(50, 3.0)
    assert np.isnan(chatterjee_xi(x, y))


def test_too_few_points_is_nan():
    assert np.isnan(chatterjee_xi([1.0], [2.0]))


def test_mismatched_shapes_raise():
    with pytest.raises(ValueError):
        chatterjee_xi([1, 2, 3], [1, 2])


def test_nan_pairs_dropped():
    x = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    y = np.array([1.0, 4.0, 9.0, np.nan, 25.0])
    # Drops the two bad pairs, leaves a clean monotone triple. n=3 caps the
    # tie-corrected estimator at 0.25 (its max for three points), so just
    # confirm it survived the NaNs, stayed finite, and read as positive/monotone.
    xi = chatterjee_xi(x, y, seed=0)
    assert np.isfinite(xi)
    assert xi == pytest.approx(0.25, abs=1e-9)


def test_reproducible_with_seed_on_ties():
    rng = np.random.default_rng(9)
    x = rng.integers(0, 5, 300).astype(float)
    y = rng.random(300)
    a = chatterjee_xi(x, y, seed=42)
    b = chatterjee_xi(x, y, seed=42)
    assert a == b


def test_correlation_matrix_shape_and_orientation():
    descriptors, responses, _ = make_synthetic_dataset(n=40, seed=1)
    df = correlation_matrix(descriptors, responses, seed=0)
    assert df.values.shape == (len(descriptors), len(responses))
    assert list(df.index) == list(descriptors.keys())
    assert list(df.columns) == list(responses.keys())
    assert np.all(np.isfinite(df.values))


def test_matrix_recovers_planted_structure():
    descriptors, responses, truth = make_synthetic_dataset(n=60, seed=2)
    df = correlation_matrix(descriptors, responses, seed=0)

    def get(row, col):
        # works for pandas DataFrame and the _SimpleMatrix fallback
        return float(df.loc[row, col])

    for resp, spec in truth.items():
        min_strong = min(get(d, resp) for d in spec["strong"])
        max_null = max(get(d, resp) for d in spec["null"])
        assert min_strong > max_null, (
            f"{resp}: strong drivers {spec['strong']} did not beat "
            f"null {spec['null']} (min_strong={min_strong:.3f}, max_null={max_null:.3f})"
        )


def test_iota_drives_heating_more_than_nfp():
    descriptors, responses, _ = make_synthetic_dataset(n=50, seed=5)
    df = correlation_matrix(descriptors, responses, seed=0)
    assert float(df.loc["iota", "peak_coil_heating"]) > float(
        df.loc["nfp", "peak_coil_heating"]
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
