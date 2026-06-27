"""Tier-8 / Part B geometry unit tests: the conformal-wall builder
(stellarator_geometry.py) is pure numpy, so it is tested in the sandbox even
though the DAGMC/transport steps run on x86. Checks the offset/normal math, mesh
watertightness, the self-intersection guard, and outward nesting -- the
properties a valid DAGMC input must have.
"""
import numpy as np
import pytest

import stellarator_geometry as g

STEM = "equil_precise_qa"


def _have_surface():
    from pathlib import Path
    return (Path(g.DATADIR) / f"{STEM}_surface.npz").exists()


pytestmark = pytest.mark.skipif(not _have_surface(),
                                reason="precise_QA surface grid not generated")


def test_surface_grid_axis_order():
    R, Z, nfp = g.load_surface(STEM)
    d = np.load(g.DATADIR / f"{STEM}_surface.npz")
    phi = d["phi"]
    # axis1 = phi (sweeps full torus); axis0 = theta (phi ~const)
    assert np.ptp(phi[0, :]) > 5.0
    assert np.ptp(phi[:, 0]) < 0.5
    assert nfp == 2


def test_normals_outward():
    R, Z, _ = g.load_surface(STEM)
    nR, nZ = g.poloidal_outward_normals(R, Z)
    # outboard (max R) normal points +R; inboard (min R) points -R, at every phi
    for j in range(0, R.shape[1], 8):
        i_out = np.argmax(R[:, j]); i_in = np.argmin(R[:, j])
        assert nR[i_out, j] > 0.5
        assert nR[i_in, j] < -0.5


def test_offset_grows_cross_section():
    R, Z, _ = g.load_surface(STEM)
    nR, nZ = g.poloidal_outward_normals(R, Z)
    Ro, Zo = g.offset_surface(R, Z, nR, nZ, 10.0)
    # bounding box must grow outward at every phi
    for j in range(0, R.shape[1], 8):
        assert Ro[:, j].max() > R[:, j].max()
        assert Ro[:, j].min() < R[:, j].min()


def test_simple_polygon_detector():
    # a square is simple; a bow-tie (figure-8) is not
    sq_r = np.array([0, 1, 1, 0.0]); sq_z = np.array([0, 0, 1, 1.0])
    assert g._poly_is_simple(sq_r, sq_z)
    bt_r = np.array([0, 1, 0, 1.0]); bt_z = np.array([0, 1, 1, 0.0])
    assert not g._poly_is_simple(bt_r, bt_z)


def test_valid_build_watertight_simple_nested(tmp_path):
    """A build with cumulative thickness < minor radius is watertight, simple,
    and nests strictly outward (the DAGMC-ready condition). Writes to an isolated
    tmp dir so it does not clobber the shared data/<stem>_geom manifest."""
    m = g.build_layers(STEM, scale=10.0, outdir=tmp_path)
    vols = []
    for L in m["layers"]:
        assert L["watertight"], f"{L['name']} not watertight"
        assert L["simple_cross_section"], f"{L['name']} self-intersects"
        vols.append(L["enclosed_volume_cm3"])
    assert all(np.diff(vols) > 0), "layers must nest strictly outward"


def test_overthick_build_is_flagged(tmp_path):
    """The self-intersection guard must CATCH an over-thick build (cumulative
    thickness > minor radius folds on the concave inboard side). Isolated outdir."""
    thick = [("blob", 400.0)]  # absurdly thick -> must fold inboard
    m = g.build_layers(STEM, layers=thick, scale=3.0, outdir=tmp_path)
    assert not m["layers"][0]["simple_cross_section"]
