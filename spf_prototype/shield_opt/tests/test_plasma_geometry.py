"""Tests for plasma_geometry: geometry-derived R0 / a_minor replacing the voxel heuristics."""
import glob
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import plasma_geometry as pg  # noqa: E402


def _circular_torus(R0=6.2, a=1.7, nth=64, nph=96):
    th = np.linspace(0, 2 * np.pi, nth, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    TH, _PH = np.meshgrid(th, ph, indexing="ij")
    return R0 + a * np.cos(TH), a * np.sin(TH)


def test_circle_recovers_R0_and_a_exactly():
    R, Z = _circular_torus(6.2, 1.7)
    g = pg.axis_minor(R, Z)
    assert g["R0"] == pytest.approx(6.2, abs=1e-9)
    assert g["a_minor"] == pytest.approx(1.7, abs=1e-9)
    assert g["a_max"] == pytest.approx(1.7, abs=1e-9)   # circle: rms == max


def test_a_minor_is_invariant_to_major_radius():
    # moving the torus out in R must not change the minor radius
    _, Z = _circular_torus(6.2, 1.7)
    R_far, _ = _circular_torus(30.0, 1.7)
    assert pg.axis_minor(R_far, Z)["a_minor"] == pytest.approx(1.7, abs=1e-9)


def test_elongation_orders_rms_below_max():
    R, Z0 = _circular_torus(6.2, 1.7)
    Z = 1.8 * Z0                               # kappa = 1.8 vertical elongation
    g = pg.axis_minor(R, Z)
    # a_max reaches the tall extent; a_rms sits strictly between the circular a and kappa*a
    assert g["a_max"] == pytest.approx(1.8 * 1.7, rel=1e-6)
    assert 1.7 < g["a_minor"] < 1.8 * 1.7


def test_scale_converts_metres_to_cm():
    R, Z = _circular_torus(6.2, 1.7)
    g = pg.axis_minor(R, Z, scale=100.0)
    assert g["R0"] == pytest.approx(620.0) and g["a_minor"] == pytest.approx(170.0)


def test_boundary_extracts_lcfs_from_3d_volume():
    # a 3-D (nrho, ntheta, nphi) volume -> boundary is the last rho surface
    R2, Z2 = _circular_torus(6.2, 1.7)
    vol = np.stack([0.5 * R2, R2], axis=0)      # inner (unused) + LCFS
    assert np.allclose(pg._boundary(vol), R2)
    g = pg.axis_minor(vol, np.stack([0.5 * Z2, Z2], axis=0))
    assert g["a_minor"] == pytest.approx(1.7, abs=1e-9)


def test_shield_band_offsets_match_radial_build():
    d_in, d_out = pg.shield_band_offsets()
    # FW 3.2 + mult 2 + breeder 50 + back_wall 4 = 59.2 ; + shield 40 = 99.2
    assert d_in == pytest.approx(59.2) and d_out == pytest.approx(99.2)


@pytest.mark.skipif(not glob.glob(os.path.join(os.path.dirname(__file__), "..", "..",
                    "data", "quasr*_surface.npz")), reason="no local surface file")
def test_real_surface_is_physical():
    f = sorted(glob.glob(os.path.join(os.path.dirname(__file__), "..", "..",
               "data", "quasr*_surface.npz")))[0]
    R0, a, amax = pg.axis_minor_from_surface(f, scale=100.0)
    assert R0 > a > 0 and amax >= a        # major > minor, outer excursion >= rms
    assert 1.5 < R0 / a < 30               # a sane stellarator aspect ratio
