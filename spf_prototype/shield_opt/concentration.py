"""Concentration metrics for a per-element loading distribution.

This module quantifies how *localized* (peaked on a few elements) versus how
*diffuse* (spread evenly) a non-negative loading array is.  Typical loadings in
this fusion-neutronics project are per-coil fast flux, per-coil nuclear heating,
or per-voxel "culprit" contributions to a hot spot.  The same math applies to
any non-negative vector.

Metrics implemented
--------------------
- Gini coefficient (`gini`):  classic inequality index from economics.
  0 = perfectly uniform loading, -> 1 = all load on a single element.
  Reference: Gini, C. (1912) "Variabilita e mutabilita".  Standard
  mean-absolute-difference form, e.g. Damgaard & Weiner (2000),
  Ecology 81(4):1139-1142, Eq. 1:
      G = ( sum_i sum_j |x_i - x_j| ) / ( 2 n sum_i x_i )
  Weighted generalization (unequal element areas/volumes) follows the
  covariance / weighted-Lorenz form used in survey statistics, e.g.
  Lerman & Yitzhaki (1989), and Handcock & Morris (1999).

- Participation ratio (`participation_ratio`):
      PR = ( sum_i x_i^2 )^2 / sum_i x_i^4
  The "effective number of significantly loaded elements": PR = N for a
  uniform load, PR = 1 for a one-hot (single loaded element) load.  Borrowed
  from the inverse participation ratio (IPR) used to characterize wavefunction
  localization in condensed-matter / Anderson-localization physics, e.g.
  Bell & Dean (1970), Discuss. Faraday Soc. 50:55; Thouless (1974),
  Phys. Rep. 13(3):93.  Here x_i plays the role of an (unnormalized) intensity
  p_i, and PR = 1 / sum(p_i^2) with p_i = x_i^2 / sum(x_j^2).

- Lorenz curve (`lorenz_curve`):  the (cumulative population fraction,
  cumulative load fraction) pairs.  The Gini coefficient is twice the area
  between the Lorenz curve and the 45-degree line of equality.

- Peaking factor (`peaking`):  the conventional max/mean ratio commonly quoted
  in neutronics for first-wall / coil loads.

- Top participation (`top_fraction_for_load`):  the fraction of elements (or of
  weighted area/volume) that carries a target fraction of the total load
  (e.g. the smallest set of coils responsible for 50% or 90% of the flux).
  This is the "participation" number used in the project's culprit maps.

All functions accept an optional `weights` array giving a per-element measure
(area, volume, number of physical sub-elements).  When omitted, every element
counts equally (weight = 1).

The module has no side effects on import.  The `__main__` block loads real
per-coil arrays from this project (if present) and prints a report for each.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "gini",
    "participation_ratio",
    "lorenz_curve",
    "peaking",
    "top_fraction_for_load",
    "concentration_report",
    "plot_lorenz",
]


def _as_loading(x, weights=None):
    """Validate and coerce inputs to (x, weights) float arrays.

    Returns non-negative float64 arrays.  Raises on negatives, NaN/Inf,
    empty input, or shape mismatch.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    if x.size == 0:
        raise ValueError("loading array is empty")
    if not np.all(np.isfinite(x)):
        raise ValueError("loading array contains NaN or Inf")
    if np.any(x < 0):
        raise ValueError("loading array must be non-negative")

    if weights is None:
        w = np.ones_like(x)
    else:
        w = np.asarray(weights, dtype=np.float64).ravel()
        if w.shape != x.shape:
            raise ValueError(
                f"weights shape {w.shape} does not match loading shape {x.shape}"
            )
        if not np.all(np.isfinite(w)):
            raise ValueError("weights contain NaN or Inf")
        if np.any(w < 0):
            raise ValueError("weights must be non-negative")
        if w.sum() <= 0:
            raise ValueError("weights sum to zero")
    return x, w


def gini(x, weights=None):
    """Gini coefficient of a non-negative loading array.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading (e.g. per-coil fast flux).
    weights : array_like, optional
        Per-element measure (area/volume/count).  If given, the Gini is the
        weighted Gini in which each element represents `weights[i]` units of
        "population" and `weights[i] * x[i]` units of load.

    Returns
    -------
    float
        Gini coefficient in [0, 1].  0 = perfectly uniform loading intensity,
        -> 1 = all load concentrated on a single element.  Returns 0.0 for an
        all-zero loading (degenerate but conventionally "perfectly equal").

    Notes
    -----
    Unweighted form (Damgaard & Weiner 2000, Ecology 81:1139, Eq. 1):

        G = sum_i sum_j |x_i - x_j| / (2 n^2 xbar)

    (equivalently divided by 2 n sum(x)).  The weighted form uses the
    weighted-Lorenz / covariance identity: with elements sorted by intensity
    x, cumulative weighted population fraction P and cumulative weighted load
    fraction L,

        G = 1 - sum_k (P_k - P_{k-1}) (L_k + L_{k-1})

    i.e. 1 minus twice the trapezoidal area under the Lorenz curve.  For equal
    weights this reduces exactly to the mean-absolute-difference formula.
    """
    x, w = _as_loading(x, weights)
    total = np.sum(w * x)
    if total <= 0:
        return 0.0

    if weights is None:
        # Direct mean-absolute-difference form (O(n log n) via sorted trick).
        # sum_i sum_j |x_i - x_j| = 2 * sum_k (2k - n - 1) x_(k)  (1-indexed,
        # ascending sort). This is exact and avoids the O(n^2) double loop.
        n = x.size
        xs = np.sort(x)
        idx = np.arange(1, n + 1)
        num = np.sum((2 * idx - n - 1) * xs)
        denom = n * xs.sum()
        return float(num / denom)

    # Weighted Gini via the Lorenz-area identity.
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ws = w[order]
    load = ws * xs
    P = np.concatenate([[0.0], np.cumsum(ws) / ws.sum()])
    L = np.concatenate([[0.0], np.cumsum(load) / load.sum()])
    # Twice the area under the Lorenz curve (trapezoid rule).
    area = np.sum((P[1:] - P[:-1]) * (L[1:] + L[:-1]))
    return float(1.0 - area)


def participation_ratio(x, weights=None):
    """Participation ratio: effective number of significantly loaded elements.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    weights : array_like, optional
        Per-element measure.  Elements are treated as carrying intensity `x`
        with multiplicity `weights`; the PR is computed on the weighted
        moments so that duplicating an element (weight 2) counts as two
        elements.

    Returns
    -------
    (PR, PR_over_N) : tuple of float
        PR : the participation ratio,
            PR = (sum_i x_i^2)^2 / sum_i x_i^4   (unweighted),
            PR = (sum_i w_i x_i^2)^2 / sum_i w_i x_i^4   (weighted).
            PR = 1 for a fully localized (one-hot) load, PR = N (or sum of
            weights) for a uniform load.
        PR_over_N : PR normalized by the effective element count
            (N unweighted, sum(weights) weighted) -> the intuitive fraction of
            elements that are "participating", in (0, 1].

    Notes
    -----
    This is the inverse participation ratio (IPR) construction from
    Anderson-localization physics (Bell & Dean 1970; Thouless 1974).  With
    p_i = x_i^2 / sum_j x_j^2 (a normalized "intensity"), PR = 1 / sum_i p_i^2.
    Returns (1.0, 1.0/N) for an all-zero loading by convention.
    """
    x, w = _as_loading(x, weights)
    n_eff = w.sum()
    s2 = np.sum(w * x**2)
    s4 = np.sum(w * x**4)
    if s4 <= 0:
        # all-zero loading: define PR = 1 (maximally localized, degenerate)
        return 1.0, float(1.0 / n_eff)
    pr = float(s2 * s2 / s4)
    return pr, float(pr / n_eff)


def lorenz_curve(x, weights=None):
    """Lorenz curve of a loading distribution.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    weights : array_like, optional
        Per-element measure (population weight).

    Returns
    -------
    (P, L) : tuple of ndarray
        P : cumulative population fraction, from 0 to 1, length n+1
            (leading 0 included so the curve starts at the origin).
        L : cumulative load fraction, from 0 to 1, length n+1.
        Elements are sorted ascending by intensity `x`, so the curve is convex
        and lies on/below the 45-degree line of equality.  Gini = 1 - 2*area
        under (P, L).
    """
    x, w = _as_loading(x, weights)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ws = w[order]
    load = ws * xs
    total_load = load.sum()
    P = np.concatenate([[0.0], np.cumsum(ws) / ws.sum()])
    if total_load <= 0:
        L = np.linspace(0.0, 0.0, xs.size + 1)
        L[-1] = 0.0
        # degenerate: define load fraction as equal to population (line of eq.)
        L = P.copy()
    else:
        L = np.concatenate([[0.0], np.cumsum(load) / total_load])
    return P, L


def peaking(x, weights=None):
    """Conventional peaking factor = max / (weighted) mean.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    weights : array_like, optional
        Per-element measure; the mean is the weight-weighted mean.

    Returns
    -------
    float
        max(x) / mean(x).  1.0 for a uniform load; large for a peaked load.
        Returns 0.0 for an all-zero loading.
    """
    x, w = _as_loading(x, weights)
    mean = np.sum(w * x) / w.sum()
    if mean <= 0:
        return 0.0
    return float(x.max() / mean)


def top_fraction_for_load(x, load_frac, weights=None):
    """Smallest fraction of elements carrying a target fraction of the load.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    load_frac : float
        Target cumulative load fraction in (0, 1], e.g. 0.5 or 0.9.
    weights : array_like, optional
        Per-element measure.  When given, the returned fraction is a fraction
        of total weighted area/volume rather than of element count.

    Returns
    -------
    float
        The fraction (of elements, or of weighted measure) of the most heavily
        loaded elements whose cumulative load first reaches `load_frac` of the
        total.  This mirrors the "participation" number in the culprit maps:
        e.g. 0.15 means the hottest 15% of the area carries `load_frac` of the
        total load.  Smaller = more concentrated.

    Notes
    -----
    Elements are ranked by intensity `x` descending; we accumulate weighted
    load until it crosses `load_frac * total`, then report the weighted
    population fraction consumed (the crossing element counted in full).
    """
    if not (0.0 < load_frac <= 1.0):
        raise ValueError("load_frac must be in (0, 1]")
    x, w = _as_loading(x, weights)
    total_load = np.sum(w * x)
    if total_load <= 0:
        return 0.0
    order = np.argsort(x, kind="mergesort")[::-1]  # descending by intensity
    ws = w[order]
    load = ws * x[order]
    cum_load = np.cumsum(load) / total_load
    cum_pop = np.cumsum(ws) / ws.sum()
    # first index where cumulative load reaches the target
    k = int(np.searchsorted(cum_load, load_frac * (1.0 - 1e-12)))
    k = min(k, cum_pop.size - 1)
    return float(cum_pop[k])


def concentration_report(x, weights=None):
    """Compute all concentration metrics and a one-line human summary.

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    weights : array_like, optional
        Per-element measure.

    Returns
    -------
    dict
        Keys:
          n            : number of elements
          n_eff        : effective element count (sum of weights)
          max, mean    : loading extremes
          peaking      : max / mean
          gini         : Gini coefficient
          pr           : participation ratio (effective # loaded elements)
          pr_over_n    : pr / n_eff (participating fraction)
          top50_frac   : fraction of elements/area carrying 50% of the load
          top90_frac   : fraction of elements/area carrying 90% of the load
          summary      : one-line human-readable string
    """
    x, w = _as_loading(x, weights)
    n = int(x.size)
    n_eff = float(w.sum())
    mean = float(np.sum(w * x) / w.sum())
    g = gini(x, weights)
    pr, pr_over_n = participation_ratio(x, weights)
    pk = peaking(x, weights)
    top50 = top_fraction_for_load(x, 0.5, weights)
    top90 = top_fraction_for_load(x, 0.9, weights)

    summary = (
        f"n={n}: Gini={g:.3f}, PR={pr:.1f}/{n} ({pr_over_n*100:.0f}% of elements "
        f"effectively loaded), peaking={pk:.2f}x; 50% of load on hottest "
        f"{top50*100:.0f}% of elements, 90% on hottest {top90*100:.0f}%."
    )

    return {
        "n": n,
        "n_eff": n_eff,
        "max": float(x.max()),
        "mean": mean,
        "peaking": pk,
        "gini": g,
        "pr": pr,
        "pr_over_n": pr_over_n,
        "top50_frac": top50,
        "top90_frac": top90,
        "summary": summary,
    }


def plot_lorenz(x, weights=None, ax=None, label=None, title=None,
                savepath=None):
    """Plot a Lorenz curve with the Gini area shaded.

    Uses the matplotlib Agg backend (no display required).

    Parameters
    ----------
    x : array_like
        Non-negative per-element loading.
    weights : array_like, optional
        Per-element measure.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure/axis is created if omitted.
    label : str, optional
        Legend label for the Lorenz curve.
    title : str, optional
        Axis title.
    savepath : str, optional
        If given, the figure is saved to this path (PNG) at 150 dpi.

    Returns
    -------
    matplotlib.axes.Axes
        The axis containing the plot.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    P, L = lorenz_curve(x, weights)
    g = gini(x, weights)

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    # Line of equality
    ax.plot([0, 1], [0, 1], "--", color="0.5", lw=1, label="Line of Equality")
    # Lorenz curve
    ax.plot(P, L, "-", lw=2, color="C0", label=label or "Lorenz Curve")
    # Shade the Gini area (between equality line and Lorenz curve)
    ax.fill_between(P, L, P, color="C0", alpha=0.2)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("Cumulative Fraction of Elements")
    ax.set_ylabel("Cumulative Fraction of Load")
    ax.set_title(title or f"Lorenz Curve (Gini = {g:.3f})")
    ax.legend(loc="upper left", fontsize=8)

    if savepath is not None:
        ax.figure.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax


def _load_array(path, key):
    """Load a single named array from an .npz, or return None if unavailable."""
    import os
    if not os.path.exists(path):
        return None
    try:
        with np.load(path, allow_pickle=True) as d:
            if key not in d:
                return None
            return np.asarray(d[key], dtype=np.float64).ravel()
    except Exception as exc:  # pragma: no cover - diagnostic path only
        print(f"  [warn] could not load {path}: {exc}")
        return None


if __name__ == "__main__":
    sources = [
        (
            "coil_base_unpol.npz [coil_fast_flux]",
            "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/shield_opt/"
            "data/processed/coil_base_unpol.npz",
            "coil_fast_flux",
        ),
        (
            "percoil_unpol_corr_breed.npz [flux]",
            "/Users/tkiker/.claude/jobs/876e30be/tmp/"
            "percoil_unpol_corr_breed.npz",
            "flux",
        ),
    ]

    for name, path, key in sources:
        print("=" * 72)
        print(name)
        print("-" * 72)
        arr = _load_array(path, key)
        if arr is None:
            print(f"  (not found: {path} key={key})")
            continue
        rep = concentration_report(arr)
        for k in ("n", "max", "mean", "peaking", "gini", "pr", "pr_over_n",
                  "top50_frac", "top90_frac"):
            v = rep[k]
            if isinstance(v, float):
                print(f"  {k:12s}= {v:.4g}")
            else:
                print(f"  {k:12s}= {v}")
        print(f"  summary: {rep['summary']}")
    print("=" * 72)
