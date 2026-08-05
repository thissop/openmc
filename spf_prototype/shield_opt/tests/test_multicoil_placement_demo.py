"""Tests for the end-to-end multi-coil placement demo (combine -> allocate -> radial build)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import attenuation_surrogate as att  # noqa: E402
import multicoil_placement_demo as demo  # noqa: E402
from thickness_field import ThicknessField  # noqa: E402


@pytest.fixture
def tf():
    tor = np.linspace(0, 360, 48, endpoint=False)
    pol = np.linspace(0, 360, 36, endpoint=False)
    t = ThicknessField(nfp=4, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                       t_breeder0=80.0, t_shield0=20.0, t_breeder_min=15.0)
    t.k = att.trade_k()
    return t


def test_fields_shapes_and_normalization(tf):
    P, B, phi = demo.coil_fields(tf, n_coils=5)
    assert P.shape == (5, *tf.shape) and B.shape == (5, *tf.shape) and len(phi) == 5
    assert np.allclose([p.max() for p in P], 1.0)          # each priority normalized


def test_pipeline_emits_valid_radial_build(tf):
    P, B, _ = demo.coil_fields(tf, n_coils=5)
    r = demo.run_pipeline(P, B, tf, budget_frac=0.3, combine_mode="sum",
                          base_dict=demo.base_radial_build(tf))
    assert r["delta"].min() >= 0 and r["delta"].max() <= tf.delta_max + 1e-9
    rb = r["radial_build"]
    t_s, t_b = rb["shield"]["thickness_matrix"], rb["breeder"]["thickness_matrix"]
    assert np.allclose(t_s + t_b, tf.t_s0 + tf.t_b0)       # envelope conserved
    assert t_b.min() >= tf.t_b_min - 1e-9                  # breeder floor


def test_multicoil_beats_single_coil_at_equal_budget(tf):
    P, B, _ = demo.coil_fields(tf, n_coils=5, hotness=[1.35, 1, 1, 1, 1])
    bf = 0.3
    single = demo.run_pipeline(P[[0]], B, tf, budget_frac=bf, combine_mode="sum")
    multi = demo.run_pipeline(P, B, tf, budget_frac=bf, combine_mode="sum")
    # identical material spent...
    assert single["added_volume"] == pytest.approx(multi["added_volume"], rel=1e-6)
    # ...but multi-coil gives a strictly lower peak coil dose
    assert multi["peak_dose"] < single["peak_dose"]


def test_single_coil_placement_migrates_the_peak(tf):
    # coil 0 hottest -> single-coil placement over-protects it and the peak moves elsewhere
    P, B, _ = demo.coil_fields(tf, n_coils=5, hotness=[1.35, 1, 1, 1, 1])
    single = demo.run_pipeline(P[[0]], B, tf, budget_frac=0.3, combine_mode="sum")
    assert int(np.argmax(single["coil_doses"])) != 0      # peak migrated off the protected coil


def test_no_shield_peak_matches_hottest_baseline(tf):
    P, B, _ = demo.coil_fields(tf, n_coils=5, hotness=[1.35, 1, 1, 1, 1])
    peak0 = demo.peak_coil_dose(B, np.zeros(tf.shape), tf.k)
    assert int(np.argmax(demo.coil_doses(B, np.zeros(tf.shape), tf.k))) == 0
    assert peak0 == pytest.approx(float(B[0].sum()))       # exp(0)=1 -> dose = sum of baseline
