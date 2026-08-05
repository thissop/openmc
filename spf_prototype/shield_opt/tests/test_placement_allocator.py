"""Tests for placement_allocator: multi-coil-aware greedy fixed-envelope shield placement."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import placement_allocator as pa  # noqa: E402
from thickness_field import ThicknessField  # noqa: E402


@pytest.fixture
def tf():
    tor = np.linspace(0.0, 360.0, 24, endpoint=False)
    pol = np.linspace(0.0, 360.0, 18, endpoint=False)
    return ThicknessField(nfp=4, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                          t_breeder0=80.0, t_shield0=20.0, t_breeder_min=15.0)


@pytest.fixture
def priority(tf):
    # a smooth inboard-band priority with distinct values (peaks at theta=pi), + mild toroidal ripple
    TOR, POL = tf.TOR, tf.POL
    P = np.exp(1.5 * (np.cos(POL - np.pi) - 1.0)) * (1.0 + 0.2 * np.cos(4 * TOR))
    return P / P.max()


# ---- combine_priorities -----------------------------------------------------------------------
def test_combine_sum_and_max_and_weighted():
    A = np.array([[1.0, 0.0], [0.0, 0.0]])
    B = np.array([[0.0, 0.0], [0.0, 2.0]])
    s = pa.combine_priorities([A, B], mode="sum")
    assert s.max() == pytest.approx(1.0)                       # normalized
    assert np.argmax(s) == 3                                   # B's cell (value 2) dominates sum
    mx = pa.combine_priorities([A, B], mode="max")
    assert mx[0, 0] > 0 and mx[1, 1] == pytest.approx(1.0)     # per-cell worst-coil envelope
    w = pa.combine_priorities([A, B], mode="weighted", weights=[10.0, 1.0])
    assert np.argmax(w) == 0                                   # A upweighted -> its cell wins


def test_combine_weighted_bad_weights_raises():
    A = np.zeros((2, 2))
    with pytest.raises(ValueError):
        pa.combine_priorities([A, A], mode="weighted", weights=[1.0])


# ---- greedy allocation ------------------------------------------------------------------------
def test_greedy_within_bounds_and_envelope_holds(tf, priority):
    delta = pa.greedy_allocate(priority, tf, budget_frac=0.4)
    assert delta.shape == tf.shape
    assert delta.min() >= 0.0 and delta.max() <= tf.delta_max + 1e-9
    t_s, t_b = tf.trade(delta)                                 # asserts envelope conserved
    assert t_b.min() >= tf.t_b_min - 1e-9                      # breeder floor respected


def test_greedy_budget_is_respected(tf, priority):
    for bf in (0.1, 0.35, 0.75):
        delta = pa.greedy_allocate(priority, tf, budget_frac=bf)
        # uniform cell area -> added volume proxy = mean(delta) = bf * delta_max
        assert tf.added_shield_volume(delta) == pytest.approx(bf * tf.delta_max, rel=1e-6)


def test_greedy_is_monotone_in_priority(tf, priority):
    delta = pa.greedy_allocate(priority, tf, budget_frac=0.5)
    order = np.argsort(priority, axis=None)[::-1]              # high -> low priority
    d_sorted = delta.ravel()[order]
    assert np.all(np.diff(d_sorted) <= 1e-9)                   # non-increasing along priority


def test_greedy_small_budget_concentrates(tf, priority):
    delta = pa.greedy_allocate(priority, tf, budget_frac=0.05)
    assert (delta > 0).sum() < delta.size / 3                  # only the top cells get shield


def test_greedy_full_budget_fills_everywhere(tf, priority):
    delta = pa.greedy_allocate(priority, tf, budget_frac=1.0)
    assert np.allclose(delta, tf.delta_max)


# ---- proportional ------------------------------------------------------------------------------
def test_proportional_monotone_and_peaks_at_delta_max(tf, priority):
    delta = pa.proportional_allocate(priority, tf, gamma=1.0)
    assert delta.max() == pytest.approx(tf.delta_max)         # priority peak -> full trade
    # order-preserving: sorting by priority sorts delta the same way
    o = np.argsort(priority, axis=None)
    assert np.all(np.diff(delta.ravel()[o]) >= -1e-9)


# ---- figure of merit: priority-aware beats uniform at equal budget -----------------------------
def test_greedy_beats_uniform_coverage_at_equal_budget(tf, priority):
    bf = 0.3
    greedy = pa.greedy_allocate(priority, tf, budget_frac=bf)
    uniform = np.full(tf.shape, bf * tf.delta_max)            # same total material, spread evenly
    assert tf.added_shield_volume(greedy) == pytest.approx(tf.added_shield_volume(uniform), rel=1e-6)
    assert pa.coverage(priority, greedy, tf) > pa.coverage(priority, uniform, tf)


# ---- smoothing ---------------------------------------------------------------------------------
def test_smooth_reduces_roughness_and_keeps_bounds(tf, priority):
    delta = pa.greedy_allocate(priority, tf, budget_frac=0.3)
    sm = pa.smooth(delta, tf, passes=2)
    assert sm.min() >= 0.0 and sm.max() <= tf.delta_max + 1e-9
    assert tf.smoothness(sm) < tf.smoothness(delta)           # smoother than the bang-bang field
