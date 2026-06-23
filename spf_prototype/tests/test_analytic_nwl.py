"""Tier-2a gate: the independent analytic NWL must match anarrima (the reference)
and our Eq.5 closed form, satisfy the free identities, and reproduce the Schwartz
Fig.2 directionality (to the anarrima reference, which the paper text rounds)."""
import os
os.environ["JAX_ENABLE_X64"] = "1"

import numpy as np
import pytest

import analytic_nwl as an

P = an.R0
WALLS = [("inboard", an.R_IN, 0.0), ("inboard", an.R_IN, 0.4),
         ("outboard", an.R_OUT, 0.0), ("outboard", an.R_OUT, 0.3),
         ("floor", an.R0, 1.0), ("floor", 0.7, 0.8)]
FACTORS = ["iso", "A", "cos2", "BC"]


@pytest.mark.parametrize("wall,r,z", WALLS)
@pytest.mark.parametrize("fac", FACTORS)
def test_quad_matches_anarrima_single_ring(wall, r, z, fac):
    q = an.g_ring_quad_scalar(P, z, r, wall, fac, an.R_IN)
    a = an.g_ring_anarrima(P, z, r, wall, fac, an.R_IN)
    assert abs(q - a) < 1e-9, f"{wall}/{fac}: quad={q} anarrima={a}"


@pytest.mark.parametrize("z", [0.0, 0.3, 0.7])
def test_eq5_closed_form_inboard_iso(z):
    q = an.g_ring_quad_scalar(P, z, an.R_IN, "inboard", "iso", an.R_IN)
    e = float(an.g_RiI_closed_eq5(P, z, an.R_IN))
    assert abs(q - e) < 1e-9


def test_free_identities():
    gi = an.g_ring_quad_scalar(P, 0.3, an.R_IN, "inboard", "iso", an.R_IN)
    ga = an.g_ring_quad_scalar(P, 0.3, an.R_IN, "inboard", "A", an.R_IN)
    gc = an.g_ring_quad_scalar(P, 0.3, an.R_IN, "inboard", "cos2", an.R_IN)
    gb = an.g_ring_quad_scalar(P, 0.3, an.R_IN, "inboard", "BC", an.R_IN)
    assert abs(gi - (ga + gc)) < 1e-12
    assert abs(gb - (0.25 * gi + 0.75 * gc)) < 1e-12


def test_plasma_quad_vs_anarrima():
    targets = an.wall_targets(n_per_wall=20)
    pq = an.plasma_patterns(targets, n_grid=60, engine=an.g_ring_quad)
    pa = an.plasma_patterns(targets, n_grid=60, engine=an.g_ring_anarrima_vec)
    for k in ("iso_g", "A_g", "cos2_g"):
        rel = np.max(np.abs(pq[k] - pa[k]) / (np.abs(pa[k]) + 1e-300))
        assert rel < 1e-6, f"{k}: max rel diff {rel}"


@pytest.fixture(scope="module")
def oracle():
    ot = [dict(wall="inboard", R=an.R_IN, Z=0.0, s=0.0),
          dict(wall="outboard", R=an.R_OUT, Z=0.0, s=0.0)]
    return an.plasma_patterns(ot, n_grid=200, engine=an.g_ring_quad)


def test_oracle_A_mode_directionality(oracle):
    di, dA = oracle["iso_directional"], oracle["A_directional"]
    inboard = dA[0] / di[0] - 1
    outboard = dA[1] / di[1] - 1
    print(f"\n[oracle] A inboard {100*inboard:+.1f}%  outboard {100*outboard:+.1f}%")
    # anarrima reference ~ +40.6% / -21.2% (paper text rounds to +43%/-22%)
    assert 0.35 < inboard < 0.45
    assert -0.25 < outboard < -0.18


def test_oracle_iso_outboard_hotter(oracle):
    ratio = oracle["iso_g"][1] / oracle["iso_g"][0]
    print(f"\n[oracle] iso outboard/inboard = {ratio:.4f}")
    assert 1.10 < ratio < 1.18  # ~+14.6% (paper text ~12%)


def test_B_equals_C_and_iso_not_C(oracle):
    di, dB, dC = oracle["iso_directional"], oracle["B_directional"], oracle["C_directional"]
    assert np.max(np.abs(dB - dC)) < 1e-9            # B and C share the shape
    assert abs(dC[0] / di[0] - 1) > 0.3              # C is NOT isotropic (spec §4.3 error)


def test_linearity(oracle):
    lin = np.max(np.abs(oracle["mixed_bracket"]
                        - (0.5 * oracle["A_bracket"] + 0.3 * oracle["B_bracket"]
                           + 0.2 * oracle["C_bracket"])))
    assert lin < 1e-9
