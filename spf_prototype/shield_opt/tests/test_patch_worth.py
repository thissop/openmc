"""Tests for the patch-worth analysis pipeline (the experiment's headline math).

Covers the pure statistics (Chatterjee xi -- which had a transcription bug once -- and the
Spearman wrapper), the patch geometry (folding, patch edges, 2-D binning), and the two field
projections attribution_patches / worth_patches on controlled synthetic adjoint+fluxmap+forward
data, including the shield-band selection that defines W.
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
PW = os.path.abspath(os.path.join(HERE, "..", "patch_worth"))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))   # shield_opt (adjoint_placement, ..)
sys.path.insert(0, PW)
import patch_worth as pw  # noqa: E402


# ---- Chatterjee xi -----------------------------------------------------------------------------
def _ref_xi(x, y):
    """Independent reference (no ties): order y by ascending x, r_i = #{y_j <= y_(i)}."""
    n = len(x)
    yy = np.asarray(y, float)[np.argsort(x, kind="mergesort")]
    r = np.array([int(np.sum(yy <= v)) for v in yy])
    return 1.0 - 3.0 * np.sum(np.abs(np.diff(r))) / (n * n - 1)


def test_chatterjee_matches_reference_on_random():
    rng = np.random.default_rng(0)
    for _ in range(5):
        x = rng.random(150); y = rng.random(150)
        assert pw.chatterjee_xi(x, y) == pytest.approx(_ref_xi(x, y), abs=1e-12)


def test_chatterjee_monotone_gives_textbook_value():
    n = 9
    x = np.arange(n, dtype=float)
    assert pw.chatterjee_xi(x, x) == pytest.approx(1.0 - 3.0 / (n + 1))       # 0.7
    assert pw.chatterjee_xi(x, -x) == pytest.approx(1.0 - 3.0 / (n + 1))      # symmetric


def test_chatterjee_independent_is_near_zero():
    rng = np.random.default_rng(1)
    xi = pw.chatterjee_xi(rng.random(4000), rng.random(4000))
    assert abs(xi) < 0.05


def test_spearman_matches_scipy():
    from scipy.stats import spearmanr
    rng = np.random.default_rng(2)
    x = rng.random(100); y = x + 0.3 * rng.standard_normal(100)
    r, p = pw.spearman(x, y)
    rr, pp = spearmanr(x, y)
    assert r == pytest.approx(rr) and p == pytest.approx(pp)


# ---- patch geometry ----------------------------------------------------------------------------
def test_fold_tor_wraps_into_one_period():
    nfp = 4; period = 2 * np.pi / nfp
    phi = np.array([0.1, period + 0.1, 2 * period + 0.1, -0.1])
    folded = pw._fold_tor(phi, nfp)
    assert np.all((folded >= 0) & (folded < period))
    assert folded[0] == pytest.approx(0.1) and folded[1] == pytest.approx(0.1)
    assert folded[2] == pytest.approx(0.1) and folded[3] == pytest.approx(period - 0.1)


def test_define_patches_edges():
    te, pe = pw.define_patches(4, 4, nfp=4)
    assert len(te) == 5 and len(pe) == 5
    assert te[0] == 0.0 and te[-1] == pytest.approx(np.pi / 2)     # one period, nfp=4
    assert pe[0] == pytest.approx(-np.pi) and pe[-1] == pytest.approx(np.pi)


def test_bin2d_conserves_and_places():
    te, pe = pw.define_patches(4, 4, nfp=4)
    phi = np.array([0.1, 0.1]); theta = np.array([-0.1, -0.1]); w = np.array([2.0, 3.0])
    H = pw._bin2d(phi, theta, w, te, pe)
    assert H.sum() == pytest.approx(5.0)                            # weight conserved
    assert H[0, 1] == pytest.approx(5.0)                           # phi~0 -> tor bin 0; theta~0- -> pol bin 1


# ---- worth_patches: shield-band selection + binning --------------------------------------------
def _mesh():
    return np.array([-200., -200, -120]), np.array([200., 200, 120]), [16, 16, 10]


def _blob_at(phi_deg, theta_deg, rho, R0, ll, ur, dim):
    """Return (ix,iy,iz) of the voxel nearest the (phi,theta,rho) point about axis R0."""
    xc, yc, zc = pw._voxel_centers(ll, ur, dim)
    R = R0 + rho * np.cos(np.radians(theta_deg))
    Z = rho * np.sin(np.radians(theta_deg))
    x, y = R * np.cos(np.radians(phi_deg)), R * np.sin(np.radians(phi_deg))
    return (int(np.argmin(np.abs(xc - x))), int(np.argmin(np.abs(yc - y))),
            int(np.argmin(np.abs(zc - Z))))


def test_worth_selects_in_band_and_bins(tmp_path):
    ll, ur, dim = _mesh(); R0, a, band = 120.0, 40.0, (20.0, 60.0)   # band rho in [60,100]
    ijk = _blob_at(30, 90, 80.0, R0, ll, ur, dim)                    # rho=80 -> inside band
    psi = np.zeros(dim); psi[ijk] = 1.0
    flux = np.zeros(dim); flux[ijk] = 2.0
    adj = tmp_path / "adj.npz"; fwd = tmp_path / "fwd.npz"
    np.savez(adj, importance=psi, lower_left=ll, upper_right=ur,
             dimension=np.array(dim), coil_cells=np.array([1]))
    np.savez(fwd, flux=flux)
    te, pe = pw.define_patches(4, 4, nfp=4)
    W, meta = pw.worth_patches(str(adj), str(fwd), R0, te, pe, a_minor=a, band_cm=band, nfp=4)
    assert W.sum() == pytest.approx(2.0)                             # C = phi*psi = 2 at the blob
    assert np.unravel_index(np.argmax(W), W.shape)[0] == 1           # phi=30deg -> tor bin 1
    assert meta["r_in"] == pytest.approx(60.0) and meta["r_out"] == pytest.approx(100.0)


def test_worth_excludes_out_of_band(tmp_path):
    ll, ur, dim = _mesh(); R0, a, band = 120.0, 40.0, (20.0, 60.0)   # band rho in [60,100]
    ijk = _blob_at(30, 90, 130.0, R0, ll, ur, dim)                   # rho=130 -> OUTSIDE band
    psi = np.zeros(dim); psi[ijk] = 1.0
    flux = np.zeros(dim); flux[ijk] = 2.0
    adj = tmp_path / "adj.npz"; fwd = tmp_path / "fwd.npz"
    np.savez(adj, importance=psi, lower_left=ll, upper_right=ur,
             dimension=np.array(dim), coil_cells=np.array([1]))
    np.savez(fwd, flux=flux)
    te, pe = pw.define_patches(4, 4, nfp=4)
    W, _ = pw.worth_patches(str(adj), str(fwd), R0, te, pe, a_minor=a, band_cm=band, nfp=4)
    assert W.sum() == pytest.approx(0.0)                             # blob beyond the shield band


def test_worth_band_defaults_to_radial_build():
    import plasma_geometry as pg
    # band_cm=None -> the actual build offsets (59.2, 99.2)
    assert pg.shield_band_offsets() == (pytest.approx(59.2), pytest.approx(99.2))


# ---- attribution_patches: localization + conservation (in-memory dicts) ------------------------
def test_attribution_localizes_and_conserves():
    ll, ur, dim = _mesh(); R0 = 120.0
    # single plasma point at (phi=30, theta=60, rho=40) -> S nonzero in one voxel; psi = ones.
    R, Z, phi = 1.4, 0.346, np.radians(30.0)          # metres (scale=100 -> 140 cm, 34.6 cm)
    fluxmap = dict(R=np.array([R]), Z=np.array([Z]), phi=np.array([phi]),
                   sqrtg=np.array([1.0]), rho=np.array([1.0]))
    adjoint = dict(importance=np.ones(dim), lower_left=ll, upper_right=ur,
                   dimension=np.array(dim), coil_cells=np.array([1]))
    te, pe = pw.define_patches(4, 4, nfp=4)
    A, R0out = pw.attribution_patches(adjoint, fluxmap, te, pe, R0=R0,
                                      emissivity="uniform", nfp=4)
    assert A.sum() == pytest.approx(1.0)              # total contributon mass = S*psi at the point
    i, j = np.unravel_index(np.argmax(A), A.shape)
    assert i == 1                                     # phi=30deg -> tor bin 1 (of [0,22.5,45,67.5,90])
    assert R0out == pytest.approx(120.0)              # explicit R0 passed through
