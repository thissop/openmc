"""Unit tests for the G4 near-source directional reducer (pure numpy, no OpenMC).
These verify the reduction MATH so the eta_source observable is trustworthy before the
Ginsburg transport run -- every target is hand-derived."""
import numpy as np
import pytest

import wall_metrics as wm

# nR=4 (r_mid = 15,25,35,45; R0=30 -> inboard = [15,25], outboard = [35,45]),
# nphi=3, nZ=5 (z_mid = -8,-4,0,4,8; midplane band |z|<0.34*8=2.72 -> only z_mid=0).
R_GRID = np.array([10., 20., 30., 40., 50.])
Z_GRID = np.linspace(-10., 10., 6)
R0 = 30.0


def test_known_contrast():
    mean = np.zeros((4, 3, 5)); sd = np.zeros((4, 3, 5))
    mean[0, :, 2] = 3.0; mean[1, :, 2] = 3.0     # inboard midplane: 3 phi * (3+3) = 18
    mean[2, :, 2] = 1.0; mean[3, :, 2] = 1.0     # outboard midplane: 3 phi * (1+1) = 6
    val, s = wm.inboard_outboard_contrast(mean, sd, R_GRID, Z_GRID, R0)
    assert val == pytest.approx((18 - 6) / (18 + 6))         # = 0.5
    assert s == pytest.approx(0.0)


def test_symmetric_is_zero():
    mean = np.ones((4, 3, 5)); sd = np.zeros((4, 3, 5))
    val, _ = wm.inboard_outboard_contrast(mean, sd, R_GRID, Z_GRID, R0)
    assert val == pytest.approx(0.0)                          # 2 inboard r-bins == 2 outboard


def test_only_midplane_band_counts():
    # flux OFF the midplane (z_mid=+-8, indices 0 and 4) must NOT enter the contrast
    mean = np.zeros((4, 3, 5)); sd = np.zeros((4, 3, 5))
    mean[0, :, 0] = 100.0                                     # inboard but off-midplane
    mean[0, :, 2] = 3.0; mean[2, :, 2] = 1.0                  # the only midplane flux
    val, _ = wm.inboard_outboard_contrast(mean, sd, R_GRID, Z_GRID, R0)
    assert val == pytest.approx((9 - 3) / (9 + 3))            # = 0.5, off-midplane ignored


def test_error_propagation():
    mean = np.zeros((4, 3, 5)); sd = np.zeros((4, 3, 5))
    mean[0, :, 2] = 3.0; mean[2, :, 2] = 1.0                  # Iin=9, Iout=3, denom=12
    sd[0, :, 2] = 0.1                                         # Vin = 3*0.01 = 0.03
    _, s = wm.inboard_outboard_contrast(mean, sd, R_GRID, Z_GRID, R0)
    assert s == pytest.approx(np.sqrt((2 * 3 / 12 ** 2) ** 2 * 0.03), rel=1e-9)


def test_empty_returns_nan():
    z = np.zeros((4, 3, 5))
    val, s = wm.inboard_outboard_contrast(z, z, R_GRID, Z_GRID, R0)
    assert np.isnan(val) and np.isnan(s)
