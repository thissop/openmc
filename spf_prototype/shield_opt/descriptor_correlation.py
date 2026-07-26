"""Chatterjee's xi rank correlation + a descriptor -> response sensitivity framework
for stellarator shield/magnet neutronics.

Why this exists
---------------
Miralles-Dolz et al. (IEEE Trans. Plasma Sci. 2026, "Stellarator Design Exploration
Using Symbolic-Regression Neutronics Surrogates") screen how *global* neutronics
responses (TBR, integral heating, dpa) depend on stellarator *descriptors* using
Chatterjee's xi correlation coefficient (Chatterjee 2021, JASA 116:2009). xi is a
rank statistic that is 0 under independence, -> 1 when Y is a measurable function of
X, and -- unlike Pearson/Spearman -- picks up NONLINEAR, non-monotonic dependence.
That last property is the whole point: a U-shaped heating response vs. rotational
transform reads as ~0 Pearson correlation but a large xi.

We reproduce that method and EXTEND it to the MAGNET/SPATIAL response quantities the
paper explicitly deferred (they screened only global/integral quantities). The new
responses -- coil peaking, coil-current concentration (Gini / participation ratio),
peak coil fast flux/heating, min magnet lifetime -- are exactly the fusion-magnet
figures of merit our shield-optimization pipeline cares about.

Chatterjee's xi (Chatterjee 2021)
---------------------------------
Given paired (X_i, Y_i), i = 1..n:
  1. Sort the pairs by X (break X-ties at random).
  2. Let r_i = #{ j : Y_j <= Y_(i) } be the "rank" of the Y that sits at sorted
     position i (rank = number of Y's less-than-or-equal, i.e. the max/right rank).

  Tie-free form:
      xi_n = 1 - 3 * sum_{i=1}^{n-1} |r_{i+1} - r_i| / (n^2 - 1)

  General (tie-corrected) form, with l_i = #{ j : Y_j >= Y_i }:
      xi_n = 1 - n * sum_{i=1}^{n-1} |r_{i+1} - r_i| / ( 2 * sum_{i=1}^{n} l_i (n - l_i) )

With no ties the two coincide: sum l_i(n-l_i) = n(n^2-1)/6, so the denominator
2*sum = n(n^2-1)/3 and the n on top collapses the general form to the tie-free one.
We default to the tie-corrected form (it is the correct estimator when Y has ties,
e.g. a QA/QH indicator response or quantized lifetimes).

Note on asymmetry: xi is NOT symmetric -- xi(X,Y) measures "is Y a function of X",
which is the direction we want for descriptor -> response screening.

Graceful degradation: pandas/matplotlib are optional. Without pandas,
`correlation_matrix` returns a small dict-backed result that still supports
`.values` / row/col access and prints as a table. Without matplotlib,
`plot_correlation_heatmap` is a no-op that warns.
"""
from __future__ import annotations

import warnings
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

# ---- optional deps -------------------------------------------------------------
try:
    import pandas as pd  # type: ignore

    _HAVE_PANDAS = True
except Exception:  # pragma: no cover - environment dependent
    pd = None  # type: ignore
    _HAVE_PANDAS = False

try:
    import matplotlib  # type: ignore

    matplotlib.use("Agg")  # headless / non-interactive: required on run nodes
    import matplotlib.pyplot as plt  # type: ignore

    _HAVE_MPL = True
except Exception:  # pragma: no cover - environment dependent
    plt = None  # type: ignore
    _HAVE_MPL = False


# ================================================================================
# Core statistic
# ================================================================================
def _max_rank(y: np.ndarray) -> np.ndarray:
    """r_i = #{ j : y_j <= y_i }  (the 'max'/right rank; ties share the top rank)."""
    order = np.argsort(y, kind="mergesort")
    y_sorted = y[order]
    # position of the last element <= y_i, +1  -> count of elements <= y_i
    counts_sorted = np.searchsorted(y_sorted, y_sorted, side="right")
    r = np.empty_like(counts_sorted)
    r[order] = counts_sorted
    return r


def _ge_rank(y: np.ndarray) -> np.ndarray:
    """l_i = #{ j : y_j >= y_i }  (max rank of -y)."""
    return _max_rank(-y)


def chatterjee_xi(
    x: Sequence[float],
    y: Sequence[float],
    ties: bool = True,
    seed: Optional[int] = None,
) -> float:
    """Chatterjee's xi correlation, measuring how strongly Y is a function of X.

    Parameters
    ----------
    x, y : array-like, same length n >= 2
        Paired samples. NaNs in either coordinate drop the pair.
    ties : bool, default True
        If True use the tie-corrected general estimator (correct when Y has ties).
        If False use the classic tie-free 1 - 3*sum|dr|/(n^2-1) form.
    seed : int, optional
        Seed for the random tie-break of the X-ordering (only matters when X has
        repeated values). Fixing it makes the statistic reproducible.

    Returns
    -------
    float
        xi in roughly [-1/2, 1]; ~0 = independent, ->1 = Y a measurable fn of X.
        Returns nan if n < 2 or Y is constant (undefined).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same shape, got {x.shape} vs {y.shape}")

    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = x.size
    if n < 2:
        return float("nan")

    rng = np.random.default_rng(seed)
    # Order by X; break X-ties at random so the estimator is well defined (Chatterjee 2021).
    jitter = rng.random(n)
    order = np.lexsort((jitter, x))  # primary key x, secondary key random
    y_ord = y[order]

    r = _max_rank(y_ord).astype(float)          # r_i along the X-ordering
    abs_dr = np.abs(np.diff(r)).sum()

    if ties:
        l = _ge_rank(y_ord).astype(float)
        denom = 2.0 * np.sum(l * (n - l))
        if denom == 0.0:
            return float("nan")  # Y constant -> dependence undefined
        return float(1.0 - n * abs_dr / denom)

    if n == 1:
        return float("nan")
    return float(1.0 - 3.0 * abs_dr / (n * n - 1.0))


# ================================================================================
# Descriptor / response definitions (the screening design)
# ================================================================================
# Stellarator descriptors. Reproduces the Miralles-Dolz global-design descriptor
# family and keeps our own naming so real QUASR/precise-QA-QH data drops straight in.
DESCRIPTORS: Dict[str, str] = {
    "nfp":    "Number of field periods (integer symmetry count)",
    "iota":   "Rotational transform iota-bar (edge or volume-averaged)",
    "A":      "Aspect ratio R0 / a",
    "kappa":  "Maximum plasma cross-section elongation",
    "qs_type": "Quasi-symmetry type indicator (QA=0, QH=1)",
    "Pf":     "Fusion power [MW]",
    "S_A":    "Plasma surface area [m^2]",
}

# Responses. These are the MAGNET/SPATIAL quantities beyond Miralles-Dolz's global
# TBR/heating/dpa set (they deferred magnets + spatial tallies). 'coil_gini' and
# 'coil_participation' come from concentration.py once it exists (soft import below);
# TBR is kept as the one shared global anchor for cross-comparison to the paper.
RESPONSES: Dict[str, str] = {
    "coil_peaking_maxmean": "Coil load peaking factor: max/mean over coil surface",
    "coil_gini":            "Coil-current / coil-load concentration (Gini coefficient)",
    "coil_participation":   "Coil load participation ratio (inverse Simpson concentration)",
    "peak_coil_fast_flux":  "Peak fast (E>0.1 MeV) neutron flux on any coil [n/cm^2/s]",
    "peak_coil_heating":    "Peak nuclear heating density in the coil pack [W/cm^3]",
    "min_magnet_lifetime":  "Minimum magnet lifetime over all coils [full-power-years]",
    "TBR":                  "Tritium breeding ratio (global anchor vs. Miralles-Dolz)",
}


# Soft import of concentration.py (may not exist yet). We only need it to *label*
# which responses are concentration-derived; the numeric matrix works regardless.
def _load_concentration():
    try:
        import concentration  # type: ignore  # noqa: F401

        return concentration
    except Exception:
        return None


_CONCENTRATION = _load_concentration()
if _CONCENTRATION is None:
    warnings.warn(
        "concentration.py not importable yet; 'coil_gini'/'coil_participation' "
        "responses fall back to caller-supplied arrays. This is expected until "
        "concentration.py lands.",
        stacklevel=2,
    )


# ================================================================================
# Correlation matrix
# ================================================================================
class _SimpleMatrix:
    """Minimal DataFrame stand-in used when pandas is unavailable.

    Supports `.values`, `.index`, `.columns`, `.loc[row, col]`, and pretty-printing,
    which is all the rest of this module and the tests rely on.
    """

    def __init__(self, values: np.ndarray, index: Sequence[str], columns: Sequence[str]):
        self.values = np.asarray(values, dtype=float)
        self.index = list(index)
        self.columns = list(columns)

    def _at(self, row: str, col: str) -> float:
        return float(self.values[self.index.index(row), self.columns.index(col)])

    class _Loc:
        def __init__(self, parent):
            self._p = parent

        def __getitem__(self, key):
            row, col = key
            return self._p._at(row, col)

    @property
    def loc(self):
        return _SimpleMatrix._Loc(self)

    def __repr__(self) -> str:
        colw = max([len(c) for c in self.columns] + [8])
        rw = max(len(r) for r in self.index)
        head = " " * rw + "  " + "  ".join(c.rjust(colw) for c in self.columns)
        lines = [head]
        for i, r in enumerate(self.index):
            row = "  ".join(f"{v:{colw}.3f}" for v in self.values[i])
            lines.append(f"{r.ljust(rw)}  {row}")
        return "\n".join(lines)


def correlation_matrix(
    descriptors: Mapping[str, Sequence[float]],
    responses: Mapping[str, Sequence[float]],
    ties: bool = True,
    seed: Optional[int] = 0,
):
    """xi correlation matrix, rows = descriptors, cols = responses.

    Each cell is xi(descriptor, response) = "how strongly is this response a
    (possibly nonlinear) function of this descriptor". Reproduces the
    Miralles-Dolz Fig. 7/8-style descriptor-vs-response heatmap.

    Returns a pandas.DataFrame if pandas is available, else a `_SimpleMatrix`
    with the same `.values` / `.loc[row, col]` access.
    """
    dnames = list(descriptors.keys())
    rnames = list(responses.keys())
    n = len(next(iter(descriptors.values())))

    mat = np.full((len(dnames), len(rnames)), np.nan)
    for i, dn in enumerate(dnames):
        xd = np.asarray(descriptors[dn], dtype=float).ravel()
        if xd.size != n:
            raise ValueError(f"descriptor '{dn}' length {xd.size} != {n}")
        for j, rn in enumerate(rnames):
            yr = np.asarray(responses[rn], dtype=float).ravel()
            if yr.size != n:
                raise ValueError(f"response '{rn}' length {yr.size} != {n}")
            mat[i, j] = chatterjee_xi(xd, yr, ties=ties, seed=seed)

    if _HAVE_PANDAS:
        return pd.DataFrame(mat, index=dnames, columns=rnames)
    return _SimpleMatrix(mat, dnames, rnames)


# ================================================================================
# Plot
# ================================================================================
def plot_correlation_heatmap(
    df,
    ax=None,
    cmap: str = "Reds",
    vmin: float = 0.0,
    vmax: float = 1.0,
    annotate: bool = True,
    title: str = "Chatterjee $\\xi$: descriptor $\\rightarrow$ response",
    savepath: Optional[str] = None,
):
    """Red heatmap of the xi matrix, styled after Miralles-Dolz Fig. 7/8.

    Returns the matplotlib Axes, or None if matplotlib is unavailable.
    """
    if not _HAVE_MPL:
        warnings.warn("matplotlib unavailable; plot_correlation_heatmap is a no-op.", stacklevel=2)
        return None

    values = np.asarray(getattr(df, "values"))
    index = list(getattr(df, "index"))
    columns = list(getattr(df, "columns"))

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(1.1 * len(columns) + 2.5, 0.7 * len(index) + 1.5)
        )
    else:
        fig = ax.figure

    im = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(columns)))
    ax.set_yticks(range(len(index)))
    ax.set_xticklabels(columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(index, fontsize=9)
    ax.set_xlabel("Response")
    ax.set_ylabel("Descriptor")
    ax.set_title(title, fontsize=11)

    if annotate:
        for i in range(len(index)):
            for j in range(len(columns)):
                v = values[i, j]
                if not np.isfinite(v):
                    continue
                # white text on dark (high-xi) cells for contrast
                color = "white" if v > 0.55 * vmax + 0.45 * vmin else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color=color, fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("$\\xi$ (0 = independent, 1 = functional)")
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax


# ================================================================================
# Synthetic demonstration / self-validation
# ================================================================================
def make_synthetic_dataset(
    n: int = 30, seed: int = 12345
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, list]]:
    """Build a fake 30-config stellarator dataset with KNOWN dependence structure.

    Mirrors the qualitative Miralles-Dolz finding that deep-component responses
    (heating in shielded regions) are driven by rotational transform and elongation.
    Returns (descriptors, responses, truth) where `truth` lists, per response, the
    descriptors that genuinely drive it ('strong') vs. those that are null.
    """
    rng = np.random.default_rng(seed)

    # ---- descriptors (independent draws over plausible ranges) ----
    nfp = rng.integers(2, 6, size=n).astype(float)          # 2..5 field periods
    iota = rng.uniform(0.2, 1.1, size=n)                    # rotational transform
    A = rng.uniform(4.0, 9.0, size=n)                       # aspect ratio
    kappa = rng.uniform(1.0, 2.2, size=n)                   # max elongation
    qs_type = rng.integers(0, 2, size=n).astype(float)      # QA=0 / QH=1
    Pf = rng.uniform(500.0, 2500.0, size=n)                 # fusion power [MW]
    S_A = rng.uniform(200.0, 900.0, size=n)                 # surface area [m^2]

    descriptors = {
        "nfp": nfp, "iota": iota, "A": A, "kappa": kappa,
        "qs_type": qs_type, "Pf": Pf, "S_A": S_A,
    }

    def z(v):  # standardize for building clean response relationships
        return (v - v.mean()) / v.std()

    zi, zk, zA, zP = z(iota), z(kappa), z(A), z(Pf)
    noise = lambda s: rng.normal(0.0, s, size=n)

    # ---- responses with planted, mostly-nonlinear structure ----
    # Peak coil heating: driven by iota (nonmonotonic, U-shaped) and elongation.
    peak_coil_heating = 3.0 * (zi ** 2) + 1.6 * zk + noise(0.25)
    # Peak coil fast flux: tracks heating physics -> iota + elongation.
    peak_coil_fast_flux = 2.2 * np.abs(zi) + 1.3 * zk + noise(0.30)
    # Coil peaking factor: driven by aspect ratio (tighter -> sharper peaking).
    coil_peaking_maxmean = 1.5 + 0.6 * (zA ** 2) + noise(0.20)
    # Concentration metrics: driven by nfp (more periods -> more even sharing).
    coil_gini = 0.4 - 0.05 * nfp + noise(0.03)
    coil_participation = 0.3 + 0.08 * nfp + noise(0.04)
    # Min magnet lifetime: inverse of heating (more heating -> shorter life).
    min_magnet_lifetime = 20.0 - 2.5 * (zi ** 2) - 1.0 * zk + noise(0.5)
    # TBR: driven by iota (as the paper found for the global anchor), mild.
    TBR = 1.05 + 0.10 * np.sin(2.5 * zi) + noise(0.02)

    responses = {
        "coil_peaking_maxmean": coil_peaking_maxmean,
        "coil_gini": coil_gini,
        "coil_participation": coil_participation,
        "peak_coil_fast_flux": peak_coil_fast_flux,
        "peak_coil_heating": peak_coil_heating,
        "min_magnet_lifetime": min_magnet_lifetime,
        "TBR": TBR,
    }

    # Ground-truth drivers we planted (used for the self-assert).
    truth = {
        "coil_peaking_maxmean": {"strong": ["A"], "null": ["qs_type", "S_A"]},
        "coil_gini": {"strong": ["nfp"], "null": ["iota", "Pf"]},
        "coil_participation": {"strong": ["nfp"], "null": ["A", "S_A"]},
        "peak_coil_fast_flux": {"strong": ["iota", "kappa"], "null": ["nfp", "S_A"]},
        "peak_coil_heating": {"strong": ["iota", "kappa"], "null": ["nfp", "Pf", "S_A"]},
        "min_magnet_lifetime": {"strong": ["iota"], "null": ["nfp", "S_A"]},
        "TBR": {"strong": ["iota"], "null": ["nfp", "A", "S_A"]},
    }
    return descriptors, responses, truth


def _demo() -> int:
    descriptors, responses, truth = make_synthetic_dataset()
    df = correlation_matrix(descriptors, responses, ties=True, seed=0)

    print("Chatterjee xi matrix  (rows = descriptors, cols = responses)")
    print("=" * 70)
    if _HAVE_PANDAS:
        with pd.option_context("display.float_format", lambda v: f"{v:6.3f}",
                               "display.width", 200):
            print(df)
    else:
        print(df)
    print()

    # Self-validation: every planted 'strong' driver must out-score every 'null'
    # driver for the same response. This is the check that the estimator recovers
    # the known dependence structure (and that we didn't transpose the matrix).
    def get(row, col):
        return float(df.loc[row, col]) if _HAVE_PANDAS else df.loc[row, col]

    failures = []
    print("Self-validation (strong drivers must beat null drivers):")
    for resp, spec in truth.items():
        strong_scores = {d: get(d, resp) for d in spec["strong"]}
        null_scores = {d: get(d, resp) for d in spec["null"]}
        min_strong = min(strong_scores.values())
        max_null = max(null_scores.values())
        ok = min_strong > max_null
        flag = "OK " if ok else "FAIL"
        print(f"  [{flag}] {resp:22s} min(strong)={min_strong:5.3f} "
              f"> max(null)={max_null:5.3f}   "
              f"strong={ {k: round(v,3) for k,v in strong_scores.items()} } "
              f"null={ {k: round(v,3) for k,v in null_scores.items()} }")
        if not ok:
            failures.append(resp)

    # Hard-assert the headline dependences the way Miralles-Dolz report them.
    assert get("iota", "peak_coil_heating") > get("nfp", "peak_coil_heating"), \
        "iota should drive coil heating more than nfp (null)"
    assert get("kappa", "peak_coil_heating") > get("Pf", "peak_coil_heating"), \
        "elongation should drive coil heating more than fusion power (null)"
    assert not failures, f"planted dependences not recovered for: {failures}"

    print("\nAll planted dependences recovered; xi framework validated.")

    if _HAVE_MPL:
        import os
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, "descriptor_xi_heatmap_synthetic.png")
        plot_correlation_heatmap(df, savepath=out,
                                 title="Chatterjee $\\xi$ (synthetic 30-config demo)")
        print(f"Heatmap written to {out}")
    else:
        print("matplotlib unavailable; skipped heatmap.")
    return 0


if __name__ == "__main__":
    raise SystemExit(_demo())
