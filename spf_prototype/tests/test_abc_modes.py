"""Tier-1 §3.1 -- symbolic + numeric identities for the (a,b,c) mode algebra.

Re-derives the integrated cross section and corner values with sympy (does NOT
trust pasted constants) and checks the numeric weight maps used by the sampler.
"""
import sympy as sp

import abc_modes as ab
from spf_mirror import make_mode_weights


def test_integrated_cross_section_identity():
    """∫ dσ/dΩ dΩ = σ0 (a + 2b/3 + c/3), verified symbolically."""
    total = ab.total_xs(sigma0=1)
    expected = ab.sa + sp.Rational(2, 3) * ab.sb + sp.Rational(1, 3) * ab.sc
    assert sp.simplify(total - expected) == 0


def test_corner_total_cross_sections():
    """Corners: nonpol 2/3, A 1, B 2/3, C 1/3 (units of σ0)."""
    total = ab.total_xs(sigma0=1)
    for name, ((a, b, c), expected) in ab.CORNERS.items():
        val = sp.simplify(total.subs({ab.sa: a, ab.sb: b, ab.sc: c}))
        assert val == expected, f"{name}: got {val}, expected {expected}"


def test_solid_angle_integrals():
    """∫ sin²θ dΩ = 8π/3 and ∫ (1/4+3/4cos²θ) dΩ = 2π."""
    assert sp.simplify(ab.int_sin2_solid_angle() - sp.Rational(8, 3) * sp.pi) == 0
    assert sp.simplify(ab.int_par_solid_angle() - 2 * sp.pi) == 0


def test_B_and_C_share_angular_shape():
    """The b and c terms have identical θ-dependence (1/4+3/4cos²θ); C is NOT
    isotropic. (The correction the spec flagged.)"""
    assert ab.angular_shape_b_minus_c() == 0


def test_weight_maps_match_closed_forms():
    """w_perp=3a/4, w_par=2b/3+c/3, eta=a+2b/3+c/3, P_perp=a/eta."""
    for (a, b, c) in [(1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0.5, 0.3, 0.2)]:
        m = make_mode_weights(a, b, c)
        an, bn, cn = m.a, m.b, m.c  # renormalized
        assert abs(m.w_perp - 0.75 * an) < 1e-15
        assert abs(m.w_par - (2 / 3 * bn + 1 / 3 * cn)) < 1e-15
        assert abs(m.eta - (an + 2 / 3 * bn + 1 / 3 * cn)) < 1e-15
        assert abs(m.P_perp - an / m.eta) < 1e-12


def test_selection_probs_sum_to_one():
    for (a, b, c) in [(1, 1, 1), (0.5, 0.3, 0.2), (0.9, 0.05, 0.05)]:
        m = make_mode_weights(a, b, c)
        P_par = 1.0 - m.P_perp
        assert abs(m.P_perp + P_par - 1.0) < 1e-12
        assert 0.0 <= m.P_perp <= 1.0


def test_collision_mode_fractions_corners():
    """Eq. 1 reproduces the corners and sums to 1 for physical spin fractions."""
    # Fully aligned (d+ = 1, t+ = 1) -> pure A.
    assert ab.collision_mode_fractions(1, 0, 0, 1, 0) == (1, 0, 0)
    # Unpolarized fuel: d+=d0=d-=1/3, t+=t-=1/2 -> (1/3,1/3,1/3).
    a, b, c = ab.collision_mode_fractions(1 / 3, 1 / 3, 1 / 3, 0.5, 0.5)
    assert abs(a - 1 / 3) < 1e-15 and abs(b - 1 / 3) < 1e-15 and abs(c - 1 / 3) < 1e-15
    assert abs((a + b + c) - 1.0) < 1e-15
