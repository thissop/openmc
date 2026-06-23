"""(a,b,c) collision-mode algebra for the spin-polarized DT source.

Numeric maps  (a,b,c) -> (w_perp, w_par, eta, P_perp)  reuse spf_mirror so the
Python/C++ samplers and this module share ONE definition. The symbolic helpers
(sympy) furnish the identities that tests/test_abc_modes.py verifies:

  * ∫ dσ/dΩ dΩ = σ0 (a + 2b/3 + c/3)
  * corner totals: nonpol 2/3, A 1, B 2/3, C 1/3   (in units of σ0)
  * B and C share the same angular shape (1/4 + 3/4 cos²θ)
  * ∫ sin²θ dΩ = 8π/3,  ∫ (1/4+3/4cos²θ) dΩ = 2π
"""
from __future__ import annotations

import sympy as sp

from spf_mirror import make_mode_weights  # single source of truth for numerics

# ---------------------------------------------------------------------------
# Numeric API
# ---------------------------------------------------------------------------


def collision_mode_fractions(d_plus, d0, d_minus, t_plus, t_minus):
    """Schwartz 2025 Eq. 1: deuteron/triton spin fractions -> (a, b, c).

    a = d+ t+ + d- t- ; b = d0 (= d0 t+ + d0 t-) ; c = d+ t- + d- t+.
    """
    a = d_plus * t_plus + d_minus * t_minus
    b = d0 * t_plus + d0 * t_minus  # == d0  when t+ + t- = 1
    c = d_plus * t_minus + d_minus * t_plus
    return a, b, c


def mode_weights(a, b, c) -> dict:
    """(a,b,c) -> dict(a,b,c,w_perp,w_par,eta,P_perp), renormalized to sum 1."""
    m = make_mode_weights(a, b, c)
    return {
        "a": m.a, "b": m.b, "c": m.c,
        "w_perp": m.w_perp, "w_par": m.w_par,
        "eta": m.eta, "P_perp": m.P_perp,
    }


def eta(a, b, c) -> float:
    """Total-rate factor σ_tot/σ0 = a + 2b/3 + c/3 (after renormalization)."""
    return mode_weights(a, b, c)["eta"]


# Corner configurations and their EXPECTED total cross section (units of σ0).
CORNERS = {
    "nonpol": ((sp.Rational(1, 3),) * 3, sp.Rational(2, 3)),
    "A": ((1, 0, 0), sp.Integer(1)),
    "B": ((0, 1, 0), sp.Rational(2, 3)),
    "C": ((0, 0, 1), sp.Rational(1, 3)),
}

# ---------------------------------------------------------------------------
# Symbolic API (sympy)
# ---------------------------------------------------------------------------
sa, sb, sc = sp.symbols("a b c", nonnegative=True)
stheta, sphi = sp.symbols("theta phi", real=True)


def diff_xs(sigma0=1):
    """Schwartz Eq. 2 differential cross section dσ/dΩ (symbolic, in a,b,c,θ)."""
    perp = sp.Rational(3, 4) * sa * sp.sin(stheta) ** 2
    par = (sp.Rational(2, 3) * sb + sp.Rational(1, 3) * sc) * (
        sp.Rational(1, 4) + sp.Rational(3, 4) * sp.cos(stheta) ** 2
    )
    return sigma0 / (2 * sp.pi) * (perp + par)


def total_xs(sigma0=1):
    """∫ dσ/dΩ dΩ over the full sphere (dΩ = sinθ dθ dφ), simplified."""
    integrand = diff_xs(sigma0) * sp.sin(stheta)
    return sp.simplify(sp.integrate(integrand, (stheta, 0, sp.pi), (sphi, 0, 2 * sp.pi)))


def int_sin2_solid_angle():
    """∫ sin²θ dΩ  -> expect 8π/3."""
    return sp.integrate(sp.sin(stheta) ** 2 * sp.sin(stheta),
                        (stheta, 0, sp.pi), (sphi, 0, 2 * sp.pi))


def int_par_solid_angle():
    """∫ (1/4 + 3/4 cos²θ) dΩ  -> expect 2π."""
    expr = (sp.Rational(1, 4) + sp.Rational(3, 4) * sp.cos(stheta) ** 2) * sp.sin(stheta)
    return sp.integrate(expr, (stheta, 0, sp.pi), (sphi, 0, 2 * sp.pi))


def angular_shape_b_minus_c():
    """θ-dependence of dσ/dΩ for pure-B minus (scaled) pure-C; zero ⇒ same shape.

    Pure B (0,1,0): bracket = (2/3)(1/4+3/4cos²θ). Pure C (0,0,1): (1/3)(1/4+3/4cos²θ).
    Their ratio is a constant (2) independent of θ ⇒ identical directionality.
    """
    bracket_B = (sp.Rational(2, 3)) * (sp.Rational(1, 4) + sp.Rational(3, 4) * sp.cos(stheta) ** 2)
    bracket_C = (sp.Rational(1, 3)) * (sp.Rational(1, 4) + sp.Rational(3, 4) * sp.cos(stheta) ** 2)
    return sp.simplify(bracket_B - 2 * bracket_C)  # expect 0
