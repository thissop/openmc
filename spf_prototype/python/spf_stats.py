"""Statistical analysis helpers for tier-1 sampler verification.

Shared by tests/test_sampler_stats.py (pass/fail gate) and
python/verify_sampler.py (diagnostics + RESULTS_tier1.md tables) so both use one
definition of every estimator. cosθ bins are equal-width in x = cosθ, i.e.
equal solid angle.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from spf_mirror import ModeWeights


def _antideriv_I(m: ModeWeights, x):
    """Antiderivative of the unnormalized intensity I(x) wrt x.

    I(x) = w_perp (1 - x²) + w_par (1/4 + 3/4 x²)
    ∫I dx = w_perp (x - x³/3) + w_par (x/4 + x³/4)
    """
    return m.w_perp * (x - x ** 3 / 3.0) + m.w_par * (x / 4.0 + x ** 3 / 4.0)


def analytic_bin_probs_costheta(m: ModeWeights, edges: np.ndarray) -> np.ndarray:
    """Exact probability mass of p(cosθ) in each [edges[i], edges[i+1]] bin."""
    Z = m.w_perp * (4.0 / 3.0) + m.w_par * 1.0  # ∫ I dx over [-1,1]
    F = _antideriv_I(m, edges)
    return np.diff(F) / Z


def chi2_costheta(cz: np.ndarray, m: ModeWeights, nbins: int = 40) -> dict:
    """Pearson chi-square of sampled cosθ against the analytic marginal pdf."""
    edges = np.linspace(-1.0, 1.0, nbins + 1)
    observed, _ = np.histogram(cz, bins=edges)
    probs = analytic_bin_probs_costheta(m, edges)
    expected = cz.size * probs
    # All expected counts must be comfortably > 5 for the chi-square to be valid.
    chi2 = float(np.sum((observed - expected) ** 2 / expected))
    dof = nbins - 1
    pval = float(stats.chi2.sf(chi2, dof))
    return {
        "chi2": chi2, "dof": dof, "chi2_per_dof": chi2 / dof, "pval": pval,
        "min_expected": float(expected.min()),
        "observed": observed, "expected": expected, "edges": edges,
    }


def chi2_uniform(values: np.ndarray, lo: float, hi: float, nbins: int = 36) -> dict:
    """Chi-square of values against a uniform distribution on [lo, hi]."""
    edges = np.linspace(lo, hi, nbins + 1)
    observed, _ = np.histogram(values, bins=edges)
    expected = np.full(nbins, values.size / nbins)
    chi2 = float(np.sum((observed - expected) ** 2 / expected))
    dof = nbins - 1
    return {"chi2": chi2, "dof": dof, "chi2_per_dof": chi2 / dof,
            "pval": float(stats.chi2.sf(chi2, dof))}


def sphere_isotropy_chi2(ux, uy, uz, n_cos: int = 12, n_phi: int = 12) -> dict:
    """Chi-square of global directions against isotropy on equal-area sphere cells
    (uniform in cosθ_global and φ_global)."""
    cz = np.asarray(uz)
    phi = np.arctan2(uy, ux)
    cos_edges = np.linspace(-1.0, 1.0, n_cos + 1)
    phi_edges = np.linspace(-np.pi, np.pi, n_phi + 1)
    observed, _, _ = np.histogram2d(cz, phi, bins=[cos_edges, phi_edges])
    n = cz.size
    expected = n / (n_cos * n_phi)
    chi2 = float(np.sum((observed - expected) ** 2 / expected))
    dof = n_cos * n_phi - 1
    return {"chi2": chi2, "dof": dof, "chi2_per_dof": chi2 / dof,
            "pval": float(stats.chi2.sf(chi2, dof))}


def mean_cos2_expected(m: ModeWeights) -> float:
    """Analytic ⟨cos²θ⟩ = ∫ x² p(x) dx.

    ∫ x²(1-x²) dx = 4/15 ; ∫ x²(1/4+3/4 x²) dx = 7/15 (both over [-1,1]).
    Oracle values: A -> 1/5, B/C -> 7/15, nonpol/iso -> 1/3.
    """
    num = m.w_perp * (4.0 / 15.0) + m.w_par * (7.0 / 15.0)
    den = m.w_perp * (4.0 / 3.0) + m.w_par * 1.0
    return num / den


def ks_2samp(a: np.ndarray, b: np.ndarray) -> dict:
    """Two-sample KS statistic + p-value."""
    res = stats.ks_2samp(a, b)
    return {"stat": float(res.statistic), "pval": float(res.pvalue)}
