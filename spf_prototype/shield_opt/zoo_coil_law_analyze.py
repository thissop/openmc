#!/usr/bin/env python
"""Analyze zoo_coil_law.csv: engineering cuts, Chatterjee-xi descriptor->response
table with bootstrap CIs, and the standoff-confound decomposition (the whole point:
separate geometry->concentration from standoff->everything).

Outputs: prints tables; writes figs to spf_prototype/figs/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from descriptor_correlation import chatterjee_xi  # noqa: E402

FIGS = ROOT / "figs"
FIGS.mkdir(exist_ok=True)

# Reactor-relevance cuts (SCIENCE_STRATEGY sec 4a/5; coil_standoff.py convention).
REACTOR_R0_CM = 1000.0        # scale QUASR (R0~1) to a 10 m reactor
BLANKET_CM = 101.0            # sol5+W0.2+steel3.8+Be2+FLiBe50+shield40 stack to coil
ASPECT_LO, ASPECT_HI = 4.0, 12.0   # reactor-plausible stellarator aspect band

DESCRIPTORS = ["d_min_over_a", "nonplanarity", "nonplanarity_max", "curvature",
               "elongation", "iota", "nfp", "qs_type"]
RESPONSES = ["point_gini", "coil_peak", "coil_gini"]


def load():
    df = pd.read_csv(HERE / "data" / "zoo_coil_law.csv")
    df["gap_cm"] = df["d_min_over_R0"] * REACTOR_R0_CM
    return df


def apply_cuts(df):
    m_gap = df["gap_cm"] >= BLANKET_CM
    m_asp = (df["aspect"] >= ASPECT_LO) & (df["aspect"] <= ASPECT_HI)
    eng = df[m_gap & m_asp].copy()
    return eng, dict(n_total=len(df), n_gap=int(m_gap.sum()),
                     n_aspect=int(m_asp.sum()), n_eng=len(eng))


def xi_ci(x, y, nperm=2000, seed=0):
    """xi point estimate + jackknife SE + permutation p-value + analytic null SE.

    NOTE: we deliberately do NOT bootstrap with replacement. Resampling rows creates
    duplicate X-values; Chatterjee's xi with random X-tie-breaking is spuriously
    inflated by ties, so a with-replacement bootstrap over-reports xi (the corrupted
    'CI' fails to even contain the point estimate). Jackknife (leave-one-out) induces
    no ties; permutation gives an exact null.

    Returns (xi, se_jack, p_perm, n). Under independence xi ~ N(0, 2/(5n))
    (Chatterjee 2021), so the null 95% one-sided threshold is 1.645*sqrt(2/(5n))."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = x.size
    if n < 5:
        return np.nan, np.nan, np.nan, n
    pt = chatterjee_xi(x, y, seed=seed)
    # jackknife SE
    jk = np.array([chatterjee_xi(np.delete(x, i), np.delete(y, i), seed=seed)
                   for i in range(n)])
    se = float(np.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2)))
    # permutation p-value (shuffle y)
    rng = np.random.default_rng(seed)
    ge = 0
    for b in range(nperm):
        if chatterjee_xi(x, rng.permutation(y), seed=b) >= pt:
            ge += 1
    p = (ge + 1) / (nperm + 1)
    return pt, se, p, n


def null_thresh(n):
    return 1.645 * np.sqrt(2.0 / (5.0 * n))


def xi_table(df, descriptors=DESCRIPTORS, responses=RESPONSES, label=""):
    thr = null_thresh(len(df))
    print(f"\n=== Chatterjee xi  descriptor -> response   {label} (n={len(df)}) ===")
    print(f"(xi +/- jackknife SE, *=perm p<0.05;  null 95% threshold xi~{thr:.2f})")
    rows = []
    for d in descriptors:
        cells = []
        for r in responses:
            sub = df
            if r == "coil_gini":       # degenerate for nc_per_hp=1
                sub = df[df["nc_per_hp"] >= 2]
            pt, se, p, n = xi_ci(sub[d], sub[r])
            cells.append((pt, se, p, n))
            rows.append(dict(descriptor=d, response=r, xi=pt, se=se, p=p, n=n))
        line = f"  {d:16s} " + "  ".join(
            f"{r}={c[0]:+.2f}+/-{c[1]:.2f}{'*' if c[2] < 0.05 else ' '}"
            for r, c in zip(responses, cells))
        print(line)
    return pd.DataFrame(rows)


def standoff_confound(df):
    """Is concentration just a standoff proxy? Three tests."""
    print("\n=== STANDOFF-CONFOUND DECOMPOSITION ===")
    resp = "point_gini"
    s = df["d_min_over_a"].values
    y = df[resp].values

    # (1) baseline: how strongly does standoff alone drive concentration?
    pt, se, p, n = xi_ci(s, y)
    print(f"(1) xi(standoff -> {resp}) = {pt:+.2f} +/-{se:.2f} (perm p={p:.3f}, n={n})")
    print(f"    -> standoff is {'NOT ' if p>=0.05 else ''}a significant driver of concentration")

    # (2) standoff-stratified xi for each descriptor: does the descriptor still drive
    #     concentration WITHIN narrow standoff bands? (breaks collinearity)
    print(f"(2) standoff-stratified xi(descriptor -> {resp})  (tertiles of standoff):")
    q = np.nanquantile(s, [1/3, 2/3])
    strat = np.digitize(s, q)
    for d in DESCRIPTORS:
        if d == "d_min_over_a":
            continue
        vals = []
        for g in range(3):
            mask = strat == g
            if mask.sum() >= 8:
                vals.append(chatterjee_xi(df[d].values[mask], y[mask], seed=0))
        if vals:
            full = chatterjee_xi(df[d].values, y, seed=0)
            print(f"    {d:16s} full={full:+.2f}   within-standoff-tertile "
                  f"mean={np.mean(vals):+.2f}  [{', '.join(f'{v:+.2f}' for v in vals)}]")

    # (3) residualize concentration on standoff (rank-isotonic-ish via binned mean),
    #     then xi(descriptor -> residual). Surviving xi = genuine, non-standoff signal.
    order = np.argsort(s)
    # smooth standoff->y with a moving average over rank (window ~ n/6)
    k = max(5, len(s) // 8)
    ys = y[order].astype(float)
    kern = np.ones(k) / k
    trend = np.convolve(np.pad(ys, (k//2, k-1-k//2), mode="edge"), kern, mode="valid")
    resid = np.empty_like(y)
    resid[order] = ys - trend
    print(f"(3) xi(descriptor -> standoff-residualized {resp})  (surviving xi = non-standoff signal):")
    for d in DESCRIPTORS:
        if d == "d_min_over_a":
            continue
        pt, se, p, n = xi_ci(df[d].values, resid)
        print(f"    {d:16s} xi={pt:+.2f} +/-{se:.2f}  (perm p={p:.3f})")

    # (4) QS-type: binary descriptor -> use Mann-Whitney U (proper for 2 groups),
    #     within full and standoff-tertiles (xi is unreliable for binary X).
    from scipy.stats import mannwhitneyu
    print("(4) QS-type effect on concentration (Mann-Whitney QA vs QH):")
    qa = y[df["qs_type"].values == 0]; qh = y[df["qs_type"].values == 1]
    if len(qa) >= 5 and len(qh) >= 5:
        U, pmw = mannwhitneyu(qa, qh, alternative="two-sided")
        print(f"    full: QA med={np.median(qa):.3f} (n={len(qa)}) vs "
              f"QH med={np.median(qh):.3f} (n={len(qh)})  p={pmw:.3f}")
    return resid


def figures(full, eng):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Fig 1: xi heatmap (engineering subset)
    from descriptor_correlation import correlation_matrix, plot_correlation_heatmap
    desc = {d: eng[d].values for d in DESCRIPTORS}
    sub2 = eng[eng["nc_per_hp"] >= 2]
    resp = {r: (eng[r].values if r != "coil_gini" else eng[r].values) for r in RESPONSES}
    df_xi = correlation_matrix(desc, resp, seed=0)
    plot_correlation_heatmap(df_xi, savepath=str(FIGS / "zoo_coil_xi_heatmap.png"),
                             title="Chatterjee xi: geometry -> coil-load concentration (free-streaming)")
    plt.close("all")

    # Fig 2: concentration vs standoff, colored by QS type
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for label, sub, c in [("full zoo", full, "0.7"), ("engineering", eng, "C3")]:
        pass
    qa = full[full["qs_type"] == 0]; qh = full[full["qs_type"] == 1]
    ax[0].scatter(qa["d_min_over_a"], qa["point_gini"], s=14, alpha=0.6, label="QA", color="C0")
    ax[0].scatter(qh["d_min_over_a"], qh["point_gini"], s=14, alpha=0.6, label="QH", color="C1")
    ax[0].axvspan(0, 0, alpha=0)  # noop
    ax[0].set_xlabel("Coil-Plasma Standoff  d_min / a")
    ax[0].set_ylabel("Coil-Load Concentration (Point Gini)")
    ax[0].set_title("Concentration vs Standoff")
    ax[0].legend()
    # Fig 2b: concentration vs elongation
    ax[1].scatter(qa["elongation"], qa["point_gini"], s=14, alpha=0.6, label="QA", color="C0")
    ax[1].scatter(qh["elongation"], qh["point_gini"], s=14, alpha=0.6, label="QH", color="C1")
    ax[1].set_xlabel("Max Plasma Elongation")
    ax[1].set_ylabel("Coil-Load Concentration (Point Gini)")
    ax[1].set_title("Concentration vs Elongation")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(FIGS / "zoo_coil_concentration_scatter.png", dpi=150)
    plt.close("all")
    print(f"\nwrote figures to {FIGS}/zoo_coil_xi_heatmap.png, zoo_coil_concentration_scatter.png")


if __name__ == "__main__":
    df = load()
    eng, cuts = apply_cuts(df)
    print("=== ENGINEERING-RELEVANT SUBSET ===")
    print(f"  total cached devices        : {cuts['n_total']}")
    print(f"  pass standoff gap>={BLANKET_CM:.0f}cm@R0=10m : {cuts['n_gap']}")
    print(f"  pass aspect in [{ASPECT_LO},{ASPECT_HI}]     : {cuts['n_aspect']}")
    print(f"  pass BOTH (engineering set) : {cuts['n_eng']}")
    print(f"  QA/QH in eng set            : {int((eng['qs_type']==0).sum())}/"
          f"{int((eng['qs_type']==1).sum())}")

    xi_table(df, label="FULL ZOO")
    xi_table(eng, label="ENGINEERING SUBSET")
    standoff_confound(eng)
    figures(df, eng)
