"""Field-parity (Part A): the B-hat computed by the C++ field models
(src/spf_field.hpp, via spf_field_driver) must equal the Python mirror
(spf_mirror.make_field) to <= 1e-12, for every bmode. This is the field analog of
the C++/Python sampler parity test (test_sampler_stats.test_cpp_python_parity) and
keeps the bit-parity discipline as the field grows (toroidal/constant/angled).
"""
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import spf_mirror as mir

DATA = Path(__file__).resolve().parents[1] / "data"
PYDIR = Path(__file__).resolve().parents[1] / "python"

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


# ----------------------------------------------------------------------------
# Field-map parity (Part B): C++ FieldMapField (spf_fieldmap.hpp) vs the Python
# mirror (fieldmap.FieldMapField) on the checked-in stand-in maps. Exercises the
# clamped R/Z + periodic-phi trilinear interp; positions sweep phi to hit the wrap.
# ----------------------------------------------------------------------------
def _fieldmap_positions():
    pts = []
    for R in (550.0, 700.0, 850.0, 1050.0):
        for phi in (0.0, 0.3, 1.1, 2.0, 3.0, -0.5, 5.5):  # exercise periodic wrap
            for Z in (-250.0, -50.0, 0.0, 120.0, 280.0):
                pts.append((R * math.cos(phi), R * math.sin(phi), Z))
    return pts


@pytest.mark.parametrize("stem", ["standin_qa", "standin_toroidal"])
def test_fieldmap_parity(run_field_driver, stem):
    from fieldmap import FieldMapField
    path = str(DATA / stem)
    pts = _fieldmap_positions()
    cpp = run_field_driver("fieldmap", pts, path=path)
    fm = FieldMapField(path)
    py = np.array([fm.bhat(p) for p in pts])
    mx = float(np.max(np.abs(cpp - py)))
    print(f"\n[field parity] fieldmap {stem} max|C++-Py| = {mx:.2e}")
    assert mx <= 1e-12
    norms = np.sqrt((cpp ** 2).sum(axis=1))
    assert float(np.max(np.abs(norms - 1.0))) <= 1e-12  # renormalized to unit


def test_fieldmap_generator_reproducible():
    """make_standin_field.py is a pure function of grid indices -> byte-identical
    output on re-run (no RNG/timestamp). Guards the committed .bin."""
    f = DATA / "standin_qa.bin"
    before = f.read_bytes()
    subprocess.run([sys.executable, str(PYDIR / "make_standin_field.py")],
                   check=True, capture_output=True)
    assert f.read_bytes() == before
