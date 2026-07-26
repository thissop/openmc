#!/usr/bin/env python
"""Regenerate data-driven talk figures for spf_petmtg.tex, house style (smplotlib).

RULE (per T.K.): every figure imports smplotlib BEFORE any plotting so the whole
deck shares the serif Super-Mongo house style. Titles/labels in Title Case.

Figures written to talk/figs/:
  zoo_descriptor_xi.png   Chatterjee-xi heatmap, geometry descriptor -> wall-load
                          response, over the 83-device free-streaming scale-up set.

Run:  python talk/make_talk_figs.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import smplotlib  # noqa: F401  house serif style — MUST precede pyplot use
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PROTO = HERE.parent
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)
sys.path.insert(0, str(PROTO / "shield_opt"))

from descriptor_correlation import correlation_matrix, plot_correlation_heatmap  # noqa: E402


def _load_zoo():
    """The 83-device free-streaming geometry/peaking scale-up set."""
    p = PROTO / "sweep" / "sweep_out" / "geometry_peaking_scaleup.json"
    rows = json.load(open(p))
    return rows


def fig_descriptor_xi():
    rows = _load_zoo()
    # Geometry descriptors (independent inputs of the configuration) ...
    descriptors = {
        "Aspect Ratio":     [r["aspect_boundary"] for r in rows],
        "Field Periods":    [r["nfp"] for r in rows],
        "Elongation":       [r["circular_std"] for r in rows],   # poloidal shape spread
        "Shape Anisotropy": [r["anisotropy"] for r in rows],
        "Field Reversal":   [r["reversal_frac"] for r in rows],
    }
    # ... vs free-streaming wall-load responses (what we care about).
    responses = {
        "Peak Load":      [r["PF_unpol"] for r in rows],
        "Concentration":  [r["C"] for r in rows],
        "Perp. Lever":    [r["PF_perp_over_unpol"] for r in rows],
        "Par. Lever":     [r["PF_par_over_unpol"] for r in rows],
        "In/Out Asym.":   [r["io_asym_unpol"] for r in rows],
    }
    n = len(rows)
    df = correlation_matrix(descriptors, responses)
    ax = plot_correlation_heatmap(df, title="")   # title lives in the LaTeX frame
    ax.set_title("")
    ax.set_xlabel("Response", fontsize=18)
    ax.set_ylabel("Descriptor", fontsize=18)
    ax.set_xticklabels(list(responses), rotation=45, ha="right", fontsize=14)
    ax.set_yticklabels(list(descriptors), fontsize=14)
    ax.figure.savefig(FIGS / "zoo_descriptor_xi.png", dpi=200, bbox_inches="tight")
    plt.close(ax.figure)
    # Honest console dump of the matrix (no cherry-picking).
    vals = df.values
    print(f"[zoo_descriptor_xi] n={n} devices")
    print("rows(descriptors):", list(descriptors))
    print("cols(responses):  ", list(responses))
    for i, dn in enumerate(descriptors):
        print(f"  {dn:16s}", "  ".join(f"{vals[i, j]:.2f}" for j in range(vals.shape[1])))
    return FIGS / "zoo_descriptor_xi.png"


def fig_magnet_matrix():
    """Magnet-side SPF matrix, smplotlib restyle of the overnight result.

    Source: docs/OVERNIGHT_2026-07-22.md 'COMPLETE MAGNET-SIDE SPF MATRIX'
    (both devices, identical radial build; QH is the golden flux reference).
    """
    modes = ["Unpol", "A\n(Perp)", "B/C\n(Par)"]
    peaking = {"QH": [2.906, 2.878, 2.808], "QA": [2.185, 2.410, 1.805]}
    ratio = {d: [p / v[0] for p in v] for d, v in peaking.items()}
    # Black-and-white encoding, consistent across both panels:
    # QH = black outline / white fill; QA = black outline / hatched.
    style = {
        "QH": dict(facecolor="white", edgecolor="black", linewidth=1.3),
        "QA": dict(facecolor="white", edgecolor="black", linewidth=1.3, hatch="////"),
    }
    import numpy as np
    x = np.arange(len(modes)); w = 0.38

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    for i, d in enumerate(["QH", "QA"]):
        axL.bar(x + (i - 0.5) * w, peaking[d], w, label=d, **style[d])
        axR.bar(x + (i - 0.5) * w, ratio[d], w, label=d, **style[d])
    axL.set_ylabel("Coil Peaking (Max/Mean)")
    axL.set_title("Magnet Peaking by Polarization Mode")
    axR.axhline(1.0, ls="--", color="k", lw=1)
    axR.set_ylabel("Peaking / Unpolarized")
    axR.set_title("Polarization Effect (Within Device)")
    axR.set_ylim(0, 1.35)
    axL.set_xticks(x); axL.set_xticklabels(modes)          # no legend on left (no room)
    axR.set_xticks(x); axR.set_xticklabels(modes)
    axR.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "magnet_polarization_matrix.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return FIGS / "magnet_polarization_matrix.png"


def fig_standoff():
    """Coil-plasma standoff, smplotlib restyle.

    Source: docs/OVERNIGHT_2026-07-22.md follow-up (LCFS-to-coil-filament).
    QH (Wiedman): min 1.63 m, mean-nearest 1.85 m.
    QA (Wechsung): min 3.11 m, mean-nearest 3.47 m.
    Bar = global-min standoff; whisker up to the mean-nearest standoff.
    """
    dev = ["QH\n(Wiedman)", "QA\n(Wechsung)"]
    mn = [1.63, 3.11]; mean_near = [1.85, 3.47]
    err = [[0, 0], [mean_near[0] - mn[0], mean_near[1] - mn[1]]]
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    bars = ax.bar(dev, mn, 0.6, yerr=err, facecolor="white", edgecolor="black",
                  linewidth=1.3, capsize=6, error_kw=dict(lw=1.4))
    bars[1].set_hatch("////")                              # QA hatched (match magnet chart)
    ax.set_ylabel("Min Plasma-Coil Standoff (m)")
    ax.set_title("QA Coils Sit $\\sim$1.9$\\times$ Further From Plasma")
    fig.tight_layout()
    fig.savefig(FIGS / "standoff_qa_qh.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return FIGS / "standoff_qa_qh.png"


def fig_decomp_iota():
    """Single-panel: rotational transform drives the first-wall peak.

    Reconstructed with the committed free-streaming machinery
    (sweep/geometry_peaking.analyze_device, PF = max/area-weighted-mean) over the
    QUASR devices that have BOTH a committed vmec fluxmap and a tabulated iota
    (data/quasr_metadata.json). xi = Chatterjee(iota, peaking).
    """
    import numpy as np
    from descriptor_correlation import chatterjee_xi
    d = json.load(open(FIGS / "decomp_iota_peak.json"))
    iota = np.array([r["iota"] for r in d])
    pf = np.array([r["PF_unpol"] for r in d])
    xi = chatterjee_xi(iota, pf)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(iota, pf, s=48, facecolor="white", edgecolor="black",
               linewidth=1.2, zorder=3)
    if pf.max() / pf.min() > 20:        # peak spans ~1..50 across the zoo
        ax.set_yscale("log")
    ax.set_xlabel("Rotational Transform $\\iota$", fontsize=14)
    ax.set_ylabel("First-Wall Peaking (Max/Mean)", fontsize=14)
    ax.set_title(f"Rotational Transform Drives The Peak "
                 f"($\\xi = {xi:.2f}$, $n = {len(d)}$)", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGS / "decomp_iota_peak.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[decomp_iota] n={len(d)} xi(iota,peaking)={xi:.3f}")
    return FIGS / "decomp_iota_peak.png"


if __name__ == "__main__":
    print("wrote", fig_descriptor_xi())
    print("wrote", fig_magnet_matrix())
    print("wrote", fig_standoff())
    print("wrote", fig_decomp_iota())

