"""Tests for multicoil_importance.combine (map-level coil combination) and patch_worth_analyze
(corner selection + the finite-difference verdict that produces the paper's conclusion)."""
import csv
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
PW = os.path.abspath(os.path.join(HERE, "..", "patch_worth"))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
sys.path.insert(0, PW)
import multicoil_importance as mi  # noqa: E402
import patch_worth_analyze as pa  # noqa: E402


# ---- multicoil_importance.combine --------------------------------------------------------------
def test_parse_coils_ranges_and_lists():
    assert mi._parse_coils("11-14") == [11, 12, 13, 14]
    assert mi._parse_coils("20,25,11") == [11, 20, 25]
    assert mi._parse_coils("11-13,20") == [11, 12, 13, 20]


def test_combine_sum_max_fluxweight():
    cells = [11, 12]
    maps = np.array([np.full((2, 2, 2), 1.0), np.full((2, 2, 2), 3.0)])   # coil 11=1, coil 12=3
    psi_sum, w = mi.combine(cells, maps, "sum")
    assert np.allclose(psi_sum, 4.0) and np.allclose(w, 1.0)
    psi_max, _ = mi.combine(cells, maps, "max")
    assert np.allclose(psi_max, 3.0)                                       # per-voxel worst coil
    psi_fw, wfw = mi.combine(cells, maps, "fluxweight", weights=np.array([2.0, 5.0]))
    assert np.allclose(psi_fw, 2.0 * 1.0 + 5.0 * 3.0)                      # weighted sum = 17
    assert np.allclose(wfw, [2.0, 5.0])


def test_combine_fluxweight_needs_weights():
    with pytest.raises(SystemExit):
        mi.combine([11], np.ones((1, 2, 2, 2)), "fluxweight", weights=None)


# ---- patch_worth_analyze._coil_R ---------------------------------------------------------------
def _coil_npz(path, cell_flux):
    cells = np.arange(11, 31)
    flux = np.array([cell_flux.get(int(c), 1.0) for c in cells], float)
    np.savez(path, coil_cells=cells, coil_fast_flux=flux)


def test_coil_R_reads_target_cell(tmp_path):
    p = tmp_path / "coil_patchbase.npz"
    _coil_npz(p, {20: 3.45e-2, 25: 8.3e-3})
    assert pa._coil_R(str(p), 20) == pytest.approx(3.45e-2)
    assert pa._coil_R(str(p), 25) == pytest.approx(8.3e-3)


# ---- pick_corners: selects the disagreement corners --------------------------------------------
def _write_tier1(path, entries):
    """entries: list of (i,j,A,W)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i_tor", "j_pol", "phi_deg", "theta_deg", "A_attribution", "W_worth"])
        for i, j, A, W in entries:
            w.writerow([i, j, 0.0, 0.0, f"{A:.6e}", f"{W:.6e}"])


def test_pick_corners_includes_disagreement(tmp_path):
    # 4 patches spanning the corners: hiA-hiW, hiA-loW, loA-hiW (the money corner), loA-loW
    tier1 = tmp_path / "t1.csv"
    _write_tier1(tier1, [(0, 0, 100.0, 100.0),   # hiA hiW
                         (0, 1, 100.0, 1.0),     # hiA loW
                         (1, 0, 1.0, 100.0),     # loA hiW  <- attribution says no, worth says yes
                         (1, 1, 1.0, 1.0)])       # loA loW
    out = tmp_path / "patches.txt"
    pa.pick_corners(str(tier1), str(out), n_per_corner=1)
    picked = {tuple(map(int, line.split())) for line in open(out) if line.strip()}
    assert (1, 0) in picked and (0, 1) in picked        # both disagreement corners chosen
    assert (0, 0) in picked and (1, 1) in picked        # and the agreement corners


# ---- finite_diff: the verdict (W_FD vs W vs A) -------------------------------------------------
def test_finite_diff_verdict_attribution_not_actionability(tmp_path):
    # Construct patches where the true finite-difference worth tracks W but NOT A.
    R0 = 1.0e-1
    added = 12 * 20.0                                    # added_cells * delta_cm, same for all
    # W ascending with patch order; A DESCENDING (anti-correlated with W); Rp set so W_FD == W.
    patches = [(2, 1, 1.0, 4.0),                         # (i, j, A, W)
               (2, 3, 2.0, 3.0),
               (1, 0, 3.0, 2.0),
               (1, 2, 4.0, 1.0)]
    tier1 = tmp_path / "t1.csv"
    _write_tier1(tier1, patches)
    base = tmp_path / "coil_patchbase.npz"; _coil_npz(base, {20: R0})
    for i, j, A, W in patches:
        tag = f"t{i}p{j}"
        # W_FD = -(Rp - R0)/added  ==  W   ->   Rp = R0 - W*added
        _coil_npz(tmp_path / f"coil_patch_{tag}.npz", {20: R0 - W * added})
        np.savez(tmp_path / f"delta_qh_patch_{tag}.npz", added_cells=12, delta_cm=20.0)
    res = pa.finite_diff(str(tier1), str(base), str(tmp_path / "coil_patch_*.npz"),
                         cell=20, out_fig=None)
    assert res["rWFD_W"] == pytest.approx(1.0)          # W_FD perfectly tracks adjoint worth W
    assert res["rWFD_A"] == pytest.approx(-1.0)         # ...and is anti-correlated with attribution
    assert res["sign_W"] == pytest.approx(1.0)          # signs all agree (all worths positive)
