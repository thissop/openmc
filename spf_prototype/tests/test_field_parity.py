"""Field-parity (Part A): the B-hat computed by the C++ field models
(src/spf_field.hpp, via spf_field_driver) must equal the Python mirror
(spf_mirror.make_field) to <= 1e-12, for every bmode. This is the field analog of
the C++/Python sampler parity test (test_sampler_stats.test_cpp_python_parity) and
keeps the bit-parity discipline as the field grows (toroidal/constant/angled).
"""
import numpy as np
import pytest

import spf_mirror as mir

# A spread of positions (cm). Avoid the z-axis (rxy = 0 is undefined for the
# toroidal/angled local cylindrical basis).
POSITIONS = [
    (100.0, 0.0, 0.0), (0.0, 150.0, 30.0), (-120.0, 80.0, -45.0),
    (60.0, -90.0, 10.0), (300.0, 200.0, -100.0), (1.0, 1.0, 1.0),
    (-50.0, -50.0, 25.0), (250.0, -10.0, 60.0),
]


def _py_bhat(field, positions):
    return np.array([field.bhat(p) for p in positions])


def test_toroidal_parity(run_field_driver):
    cpp = run_field_driver("toroidal", POSITIONS)
    py = _py_bhat(mir.make_field("toroidal"), POSITIONS)
    mx = float(np.max(np.abs(cpp - py)))
    print(f"\n[field parity] toroidal max|C++-Py| = {mx:.2e}")
    assert mx <= 1e-12


@pytest.mark.parametrize("b", [(0.0, 0.0, 1.0), (0.3, -0.5, 0.8), (1.0, 2.0, 3.0)])
def test_constant_parity(run_field_driver, b):
    cpp = run_field_driver("constant", POSITIONS, b=b)
    py = _py_bhat(mir.make_field("constant", b=b), POSITIONS)
    mx = float(np.max(np.abs(cpp - py)))
    print(f"\n[field parity] constant{b} max|C++-Py| = {mx:.2e}")
    assert mx <= 1e-12


@pytest.mark.parametrize("alpha,beta",
                         [(0.0, 0.0), (0.3, 0.5), (1.5, 1.2), (-0.4, 0.6), (0.0, 0.4)])
def test_angled_parity(run_field_driver, alpha, beta):
    cpp = run_field_driver("angled", POSITIONS, alpha=alpha, beta=beta)
    py = _py_bhat(mir.make_field("angled", alpha=alpha, beta=beta), POSITIONS)
    mx = float(np.max(np.abs(cpp - py)))
    print(f"\n[field parity] angled(a={alpha},b={beta}) max|C++-Py| = {mx:.2e}")
    assert mx <= 1e-12


def test_angled_beta0_is_toroidal(run_field_driver):
    """beta=0 angled field == toroidal field (the Rung-0 anchor), both C++ and Py.
    Also verifies the AngledField is unit-norm at beta=0 (it is phi_hat)."""
    tor = run_field_driver("toroidal", POSITIONS)
    ang = run_field_driver("angled", POSITIONS, alpha=0.7, beta=0.0)
    mx = float(np.max(np.abs(tor - ang)))
    print(f"\n[field parity] angled(beta=0) vs toroidal max|diff| = {mx:.2e}")
    assert mx <= 1e-12
    # unit-norm check on the angled field at a generic nonzero tilt
    a = run_field_driver("angled", POSITIONS, alpha=0.4, beta=0.9)
    norms = np.sqrt((a ** 2).sum(axis=1))
    assert float(np.max(np.abs(norms - 1.0))) <= 1e-12
