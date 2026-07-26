"""Tests for the Bayesian shield/standoff co-optimizer.

Covered:
  * objective DECREASES when shield is added over the hot region vs a zero-shield design;
  * the TBR-floor penalty ACTIVATES when the breeder is thinned too far;
  * the optimizer IMPROVES on its initial point on the synthetic problem;
  * bookkeeping helpers (peak dose, magnet lifetime inverse, breeder-sacrificed fraction).
"""
import os
import sys

import numpy as np
import pytest

# make the module importable whether pytest is run from shield_opt/ or tests/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import bayes_optimizer as bo  # noqa: E402


@pytest.fixture(scope="module")
def problem():
    p, _basis, _tf = bo.make_synthetic_problem()
    return p


def _zero_shield_design(p):
    """No standoff, DC bias very negative -> sigmoid ~ 0 -> ~zero shield trade."""
    return p.uniform_field_design(bias=-6.0, standoff=0.0)


def test_objective_decreases_when_shielding_the_hot_region(problem):
    p = problem
    x_zero = _zero_shield_design(p)
    J_zero = bo.objective(x_zero, None, None, None, None, problem=p)

    # A design that raises the uniform shield field -> shield added everywhere, including over the
    # hot inboard band. Modest amount (bias=-1 -> sigmoid~0.27) so the TBR floor is not the story.
    x_shield = p.uniform_field_design(bias=-1.0, standoff=0.0)
    J_shield = bo.objective(x_shield, None, None, None, None, problem=p)

    assert J_shield < J_zero, f"adding shield should lower J: {J_shield} !< {J_zero}"

    # and the actual peak dose must drop too
    assert bo.peak_dose(x_shield, p) < bo.peak_dose(x_zero, p)


def test_tbr_floor_penalty_activates_when_breeder_thinned_too_far(problem):
    p = problem
    # slam the field to saturation everywhere -> delta ~ delta_max -> breeder thinned to the bone.
    x_greedy = p.uniform_field_design(bias=12.0, standoff=0.0)   # sigmoid saturates -> max trade
    tbr = p.tbr_of(x_greedy)
    assert tbr < p.tbr_floor, f"expected TBR below floor, got {tbr}"

    # the penalty term should be a large, positive contribution: compare J with penalty on vs off.
    J_with = bo.objective(x_greedy, None, None, None, None, problem=p)

    p_nopen = bo.OptProblem(p.baseline_coil_dose_field, p.tf, p.basis, p.tbr0,
                            tbr_floor=p.tbr_floor, penalty=0.0)
    J_without = bo.objective(x_greedy, None, None, None, None, problem=p_nopen)

    assert J_with > J_without, "TBR-floor penalty should raise the objective"
    assert (J_with - J_without) > 1.0, "penalty should be a substantial, not token, term"


def test_tbr_floor_respected_for_modest_shield(problem):
    p = problem
    x = p.uniform_field_design(bias=-1.0, standoff=0.0)  # sigmoid~0.27 -> ~18 cm traded
    assert p.tbr_of(x) >= p.tbr_floor  # modest trade keeps TBR healthy


def test_optimizer_improves_on_initial_point(problem):
    p = problem
    x_start = _zero_shield_design(p)
    J_start = bo.objective(x_start, None, None, None, None, problem=p)

    res = bo.optimize(p, n_bootstrap=20, n_verify=6, n_init=6, x0=x_start, seed=3, verbose=False)

    assert res["best_value"] < J_start, (
        f"optimizer failed to beat start: {res['best_value']} !< {J_start}")
    # and it must not have cheated past the TBR floor
    assert res["tbr"] >= p.tbr_floor - 1e-6, f"optimum violates TBR floor: {res['tbr']}"
    # peak dose must be genuinely lower than the starting peak dose
    assert res["peak_dose"] < bo.peak_dose(x_start, p)


def test_magnet_lifetime_is_inverse_of_peak_dose(problem):
    p = problem
    x = _zero_shield_design(p)
    life = bo.magnet_lifetime(x, p)
    assert life == pytest.approx(p.life_ref / bo.peak_dose(x, p), rel=1e-9)


def test_breeder_sacrificed_fraction_monotone(problem):
    p = problem
    x_lo = p.uniform_field_design(bias=-2.0, standoff=0.0)
    x_hi = p.uniform_field_design(bias=2.0, standoff=0.0)
    assert bo.breeder_sacrificed_fraction(x_hi, p) > bo.breeder_sacrificed_fraction(x_lo, p)
    assert 0.0 <= bo.breeder_sacrificed_fraction(x_lo, p) <= 1.0


def test_standoff_reduces_dose(problem):
    p = problem
    x0 = _zero_shield_design(p)
    x1 = x0.copy(); x1[0] = 10.0  # 10 cm standoff
    assert bo.peak_dose(x1, p) < bo.peak_dose(x0, p)


def test_backend_reported(problem):
    p = problem
    res = bo.optimize(p, n_bootstrap=4, n_verify=2, n_init=4, seed=0, verbose=False)
    assert res["backend"] in ("skopt", "botorch", "gp_ei")
