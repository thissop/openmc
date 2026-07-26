"""Tests for the concentration metrics module.

Checks the analytic corner cases (uniform, one-hot) and monotonicity:
concentrating load must raise the Gini coefficient and lower the participation
ratio.
"""

import os
import sys

import numpy as np
import pytest

# Make the parent directory importable so `import concentration` works whether
# pytest is run from shield_opt/ or shield_opt/tests/.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import concentration as C  # noqa: E402


# ---------------------------------------------------------------------------
# Gini corner cases
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [5, 20, 100])
def test_gini_uniform_is_zero(n):
    x = np.ones(n)
    assert C.gini(x) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("n", [2, 5, 20, 100])
def test_gini_one_hot_is_n_minus_1_over_n(n):
    x = np.zeros(n)
    x[0] = 1.0
    # Exact value for a single loaded element is (n-1)/n.
    assert C.gini(x) == pytest.approx((n - 1) / n, abs=1e-12)


def test_gini_between_zero_and_one():
    rng = np.random.default_rng(0)
    for _ in range(20):
        x = rng.random(30)
        g = C.gini(x)
        assert 0.0 <= g < 1.0


def test_gini_all_zero_returns_zero():
    assert C.gini(np.zeros(10)) == 0.0


def test_gini_scale_invariant():
    rng = np.random.default_rng(1)
    x = rng.random(50)
    assert C.gini(x) == pytest.approx(C.gini(1e6 * x), abs=1e-12)


# ---------------------------------------------------------------------------
# Participation ratio corner cases
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [5, 20, 100])
def test_pr_uniform_is_n(n):
    x = np.ones(n)
    pr, pr_over_n = C.participation_ratio(x)
    assert pr == pytest.approx(n, rel=1e-12)
    assert pr_over_n == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("n", [2, 5, 20, 100])
def test_pr_one_hot_is_one(n):
    x = np.zeros(n)
    x[0] = 1.0
    pr, pr_over_n = C.participation_ratio(x)
    assert pr == pytest.approx(1.0, rel=1e-12)
    assert pr_over_n == pytest.approx(1.0 / n, rel=1e-12)


def test_pr_bounds():
    rng = np.random.default_rng(2)
    n = 40
    for _ in range(20):
        x = rng.random(n)
        pr, pr_over_n = C.participation_ratio(x)
        assert 1.0 <= pr <= n + 1e-9
        assert 0.0 < pr_over_n <= 1.0 + 1e-9


# ---------------------------------------------------------------------------
# Monotonicity: concentrating load raises Gini, lowers PR
# ---------------------------------------------------------------------------

def _concentrate(x, alpha):
    """Push load toward the largest element by raising to power alpha >= 1."""
    return x ** alpha


def test_monotonic_gini_up_pr_down():
    rng = np.random.default_rng(3)
    x = rng.random(50) + 0.1  # strictly positive baseline
    prev_g = C.gini(x)
    prev_pr = C.participation_ratio(x)[0]
    for alpha in [1.5, 2.0, 3.0, 5.0]:
        xc = _concentrate(x, alpha)
        g = C.gini(xc)
        pr = C.participation_ratio(xc)[0]
        assert g > prev_g, f"Gini did not increase at alpha={alpha}"
        assert pr < prev_pr, f"PR did not decrease at alpha={alpha}"
        prev_g, prev_pr = g, pr


def test_interpolation_uniform_to_onehot():
    # Moving a uniform load continuously toward one-hot must monotonically
    # raise Gini from 0 toward (n-1)/n.
    n = 20
    uniform = np.ones(n)
    onehot = np.zeros(n)
    onehot[0] = 1.0
    prev = -1.0
    for t in np.linspace(0.0, 1.0, 11):
        x = (1 - t) * (uniform / uniform.sum()) + t * onehot
        g = C.gini(x)
        assert g >= prev - 1e-12
        prev = g
    assert C.gini(onehot) == pytest.approx((n - 1) / n)


# ---------------------------------------------------------------------------
# Lorenz curve
# ---------------------------------------------------------------------------

def test_lorenz_endpoints_and_gini_consistency():
    rng = np.random.default_rng(4)
    x = rng.random(30)
    P, L = C.lorenz_curve(x)
    assert P[0] == 0.0 and L[0] == 0.0
    assert P[-1] == pytest.approx(1.0)
    assert L[-1] == pytest.approx(1.0)
    # Lorenz curve lies on or below the equality line.
    assert np.all(L <= P + 1e-12)
    # Gini = 1 - 2*area under Lorenz curve.
    trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))
    area = trapezoid(L, P)
    assert C.gini(x) == pytest.approx(1.0 - 2.0 * area, abs=1e-12)


# ---------------------------------------------------------------------------
# Peaking
# ---------------------------------------------------------------------------

def test_peaking_uniform_is_one():
    assert C.peaking(np.ones(15)) == pytest.approx(1.0)


def test_peaking_one_hot_is_n():
    n = 8
    x = np.zeros(n)
    x[0] = 1.0
    assert C.peaking(x) == pytest.approx(float(n))


# ---------------------------------------------------------------------------
# top_fraction_for_load
# ---------------------------------------------------------------------------

def test_top_fraction_uniform():
    # For a uniform load, ~load_frac of elements carry load_frac of the load.
    n = 100
    x = np.ones(n)
    assert C.top_fraction_for_load(x, 0.5) == pytest.approx(0.5, abs=0.02)
    assert C.top_fraction_for_load(x, 0.9) == pytest.approx(0.9, abs=0.02)


def test_top_fraction_one_hot():
    n = 50
    x = np.zeros(n)
    x[0] = 1.0
    # A single element carries 100% of the load -> 1/n of elements.
    assert C.top_fraction_for_load(x, 0.5) == pytest.approx(1.0 / n)
    assert C.top_fraction_for_load(x, 0.9) == pytest.approx(1.0 / n)


def test_top_fraction_monotone_in_target():
    rng = np.random.default_rng(5)
    x = rng.random(40)
    f50 = C.top_fraction_for_load(x, 0.5)
    f90 = C.top_fraction_for_load(x, 0.9)
    assert f50 <= f90


def test_top_fraction_invalid():
    with pytest.raises(ValueError):
        C.top_fraction_for_load(np.ones(5), 0.0)
    with pytest.raises(ValueError):
        C.top_fraction_for_load(np.ones(5), 1.5)


# ---------------------------------------------------------------------------
# Weighted variants
# ---------------------------------------------------------------------------

def test_weighted_gini_matches_replication():
    # A weighted Gini with integer weights must equal the unweighted Gini of
    # the replicated array.
    x = np.array([1.0, 2.0, 5.0, 8.0])
    w = np.array([1.0, 2.0, 3.0, 1.0])
    replicated = np.repeat(x, w.astype(int))
    assert C.gini(x, weights=w) == pytest.approx(C.gini(replicated), abs=1e-12)


def test_weighted_uniform_intensity_is_zero_gini():
    # Uniform intensity with non-uniform weights is still perfectly equal.
    x = np.full(10, 3.0)
    w = np.arange(1, 11, dtype=float)
    assert C.gini(x, weights=w) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def test_concentration_report_keys():
    x = np.array([0.1, 0.2, 0.3, 5.0])
    rep = C.concentration_report(x)
    for k in ("n", "n_eff", "max", "mean", "peaking", "gini", "pr",
              "pr_over_n", "top50_frac", "top90_frac", "summary"):
        assert k in rep
    assert rep["n"] == 4
    assert isinstance(rep["summary"], str)


def test_input_validation():
    with pytest.raises(ValueError):
        C.gini(np.array([-1.0, 2.0]))
    with pytest.raises(ValueError):
        C.gini(np.array([]))
    with pytest.raises(ValueError):
        C.gini(np.array([1.0, np.nan]))
