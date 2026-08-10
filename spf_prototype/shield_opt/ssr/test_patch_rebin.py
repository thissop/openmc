"""Tests for the offline re-binner on a synthetic surface source (no cluster / no real h5 needed)."""
import numpy as np
import h5py
import pytest

from patch_rebin import (load_surface_source, source_to_angles, patch_edges,
                         bin_patches, rebin_load)

R0 = 1367.0  # cm, matches the QH build major radius


def _write_source(path, r, E=None, wgt=None):
    n = len(r)
    dt = np.dtype([("r", "<f8", (3,)), ("u", "<f8", (3,)), ("E", "<f8"),
                   ("time", "<f8"), ("wgt", "<f8"), ("delayed_group", "<i4"),
                   ("surf_id", "<i4"), ("particle", "<i4")])
    sb = np.zeros(n, dtype=dt)
    sb["r"] = r
    sb["u"] = np.tile([0, 0, 1.0], (n, 1))
    sb["E"] = 14.06e6 if E is None else E
    sb["wgt"] = 1.0 if wgt is None else wgt
    with h5py.File(path, "w") as f:
        f.create_dataset("source_bank", data=sb)


def _wall_point(phi, theta, R0=R0, a=200.0):
    R = R0 + a * np.cos(theta); Z = a * np.sin(theta)
    return [R * np.cos(phi), R * np.sin(phi), Z]


def test_angles_roundtrip():
    pts = np.array([_wall_point(0.3, 0.5), _wall_point(1.0, -1.2)])
    phi, theta = source_to_angles(pts, R0, nfp=1)
    assert np.allclose(phi, [0.3, 1.0], atol=1e-6)
    assert np.allclose(theta, [0.5, -1.2], atol=1e-6)


def test_fold_and_bin_places_correctly(tmp_path):
    # one particle at (phi=10deg, theta=90deg) -> tor bin 0, pol bin 2 of a 4x4 grid (nfp=4)
    p = np.radians(10.0); t = np.radians(90.0)
    path = tmp_path / "s.h5"; _write_source(str(path), np.array([_wall_point(p, t)]))
    src = load_surface_source(str(path))
    m, _ = rebin_load(src, 4, 4, nfp=4, R0=R0, weight="flux", reduce="sum")
    assert m.sum() == pytest.approx(1.0)
    i, j = np.unravel_index(np.argmax(m), m.shape)
    te, pe = patch_edges(4, 4, nfp=4)
    assert te[i] <= p < te[i + 1] and pe[j] <= t < pe[j + 1]


def test_energy_weight_is_wgt_times_E(tmp_path):
    path = tmp_path / "s.h5"
    _write_source(str(path), np.array([_wall_point(0.1, 0.1)]), E=1.0e6, wgt=2.0)
    src = load_surface_source(str(path))
    m, _ = rebin_load(src, 2, 2, nfp=4, R0=R0, weight="energy", reduce="sum")
    assert m.sum() == pytest.approx(2.0e6)          # wgt*E


def test_max_vs_sum_reduction(tmp_path):
    # two particles in the SAME patch: sum adds, max takes the larger
    pts = np.array([_wall_point(0.1, 0.1), _wall_point(0.12, 0.1)])
    path = tmp_path / "s.h5"; _write_source(str(path), pts, E=np.array([1.0e6, 3.0e6]))
    src = load_surface_source(str(path))
    msum, _ = rebin_load(src, 2, 2, nfp=4, R0=R0, weight="energy", reduce="sum")
    mmax, _ = rebin_load(src, 2, 2, nfp=4, R0=R0, weight="energy", reduce="max")
    assert msum.max() == pytest.approx(4.0e6) and mmax.max() == pytest.approx(3.0e6)


def test_fast_filter(tmp_path):
    pts = np.array([_wall_point(0.1, 0.1), _wall_point(0.1, 0.1)])
    path = tmp_path / "s.h5"; _write_source(str(path), pts, E=np.array([2.0e6, 0.01e6]))  # one fast, one thermal
    src = load_surface_source(str(path))
    m, _ = rebin_load(src, 2, 2, nfp=4, R0=R0, weight="fast", reduce="sum")
    assert m.sum() == pytest.approx(1.0)            # only the fast one counts


def test_area_normalization_runs_if_netcdf(tmp_path):
    wout = "../../talk/geom/wout_qh.nc"
    import os
    if not os.path.exists(wout):
        pytest.skip("wout not present locally")
    from patch_rebin import patch_areas
    te, pe = patch_edges(4, 4, nfp=4)
    A = patch_areas(te, pe, wout)
    assert A.shape == (4, 4) and np.all(A > 0)      # every patch has positive physical area
    assert A.max() / A.min() > 1.05                 # equal-angle patches have UNEQUAL area (the confounder)
