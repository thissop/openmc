"""Phase E scaffolding: the exponential-attenuation surrogate that gives a cheap, noise-free
gradient of coil dose w.r.t. the fixed-envelope shield<->breeder trade.

Why this exists: OpenMC is a stochastic black box you cannot differentiate through, and
finite-difference gradients are wrecked by Monte-Carlo noise (SIMSOPT can't help here). But
shielding is nearly exponential, so the *screening* gradient is analytic:

    coil fast-neutron dose along a sightline  ~  D0 * exp( - sum_i t_i / lambda_i )

Under the fixed-envelope trade (shield += delta, breeder -= delta at the same (theta,phi)):

    D(delta) = D0 * exp( - delta * (1/lambda_shield - 1/lambda_breeder) )
             = D0 * exp( - k * delta ),     k = 1/lambda_shield - 1/lambda_breeder

k > 0 because the WC/steel shield removes fast neutrons in FEWER cm than FLiBe -> trading
breeder for shield over a hot coil always REDUCES its dose (until TBR bites). The gradient
dD/ddelta = -k D is exact for the model. We calibrate lambda_shield, lambda_breeder (and the
per-cell TBR sensitivity) to a handful of real OpenMC runs; the surrogate then optimizes for
free and OpenMC only *verifies/corrects* the candidate.

This module is transport-agnostic: it takes a baseline coil-dose field (from OpenMC, per
(theta,phi) sightline or per coil) and a delta field, and returns dose + gradient + TBR cost.
"""
from __future__ import annotations

import numpy as np

# Effective fast-neutron (E > ~0.1 MeV) removal MFP in cm. Defaults are literature-scale
# starting points; CALIBRATE to OpenMC (fit to dose vs. uniform-thickness scans) before trust.
LAMBDA_DEFAULT = {
    "breeder": 17.0,      # FLiBe: lighter, longer removal MFP
    "shield": 8.0,        # WC/steel: dense, short removal MFP  (< breeder -> trade helps)
    "first_wall": 6.0,
    "back_wall": 6.0,
    "vacuum_vessel": 6.0,
}


def trade_k(lam_shield=LAMBDA_DEFAULT["shield"], lam_breeder=LAMBDA_DEFAULT["breeder"]):
    """The trade attenuation coefficient k = 1/lam_shield - 1/lam_breeder (cm^-1). >0 => helps."""
    return 1.0 / lam_shield - 1.0 / lam_breeder


def coil_dose(baseline_dose, delta, k=None):
    """Surrogate coil dose under the fixed-envelope trade. baseline_dose and delta are
    per-sightline / per-(theta,phi) arrays (same shape). Returns the perturbed dose field."""
    if k is None:
        k = trade_k()
    return np.asarray(baseline_dose) * np.exp(-k * np.asarray(delta))


def dose_gradient(baseline_dose, delta, k=None):
    """Exact model gradient dD/ddelta = -k * D(delta). Noise-free, per cell."""
    if k is None:
        k = trade_k()
    return -k * coil_dose(baseline_dose, delta, k)


def robust_coil_objective(dose_field, rho=8.0, weights=None):
    """Kreisselmeier-Steinhauser 'soft max' of the coil dose -- optimize the hot TAIL, not the
    single noisy hottest cell. rho tunes how peaked. Returns (J, dJ/ddose_field)."""
    d = np.asarray(dose_field, float)
    w = np.ones_like(d) if weights is None else np.asarray(weights, float)
    dref = d.max() + 1e-30
    z = rho * d / dref
    zmax = z.max()
    e = w * np.exp(z - zmax)
    S = e.sum()
    J = dref / rho * (zmax + np.log(S))
    dJ = e / S                       # d J / d(d_i), chain-ruled for the descent
    return float(J), dJ


def tbr_change(delta, breeder_sensitivity, cell_area_frac):
    """Linearized global TBR change from thinning the breeder by `delta` (cm) over each cell.
    breeder_sensitivity = dTBR per cm of local breeder (calibrate to OpenMC; ~1e-3..5e-3 /cm,
    and it SATURATES as breeder thickens, so thinning an already-thick breeder costs little).
    cell_area_frac sums to 1. Returns delta_TBR (negative = breeding lost)."""
    return -float(np.sum(np.asarray(breeder_sensitivity) * np.asarray(delta)
                         * np.asarray(cell_area_frac)))


def calibrate_lambdas(thicknesses, doses, layer):
    """Fit an effective removal MFP for one layer from an OpenMC uniform-thickness scan:
    log(dose) = c - t/lambda  ->  lambda = -1/slope. `thicknesses` (cm), `doses` (coil metric)."""
    t = np.asarray(thicknesses, float)
    y = np.log(np.asarray(doses, float) + 1e-300)
    A = np.vstack([np.ones_like(t), t]).T
    c, slope = np.linalg.lstsq(A, y, rcond=None)[0]
    lam = -1.0 / slope if slope < 0 else np.inf
    return dict(layer=layer, lambda_cm=float(lam), intercept=float(c))


if __name__ == "__main__":  # tiny self-check
    k = trade_k()
    print(f"trade k = {k:.4f} /cm  (>0 means breeder->shield reduces coil dose)")
    base = np.ones((5, 5)) * 100.0
    for dd in (0, 10, 25, 40):
        d = np.full_like(base, float(dd))
        print(f"  delta={dd:2d} cm -> coil dose x{coil_dose(base, d).mean()/100:.3f}, "
              f"grad {dose_gradient(base, d).mean():+.2f}")
