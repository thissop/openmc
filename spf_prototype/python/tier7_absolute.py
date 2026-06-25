#!/usr/bin/env python
"""Tier 7: the ABSOLUTE (un-normalized) rate effect. Tiers 2b/3/6 compared modes
at *constant fusion rate* to isolate directionality. Here we restore the rate:
at fixed fuel density a fuel mix produces neutrons in proportion to the total-rate
factor eta = a + 2b/3 + c/3, so pure-A fuel makes +50% more neutrons and pure-C
makes -50%. We scale the (existing) per-mode wall load + TBR by eta -- pure
postprocessing of the Tier 5/6 statepoints, no new transport runs.

Run:  $HOME/spf_venv/bin/python spf_prototype/python/tier7_absolute.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
FIGDIR = REPO / "spf_prototype" / "figs"

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import smplotlib  # noqa: E402,F401  (scientific style; sets rcParams on import)
import matplotlib.pyplot as plt  # noqa: E402

import reactor_model as rm  # noqa: E402

MODES = {"nonpol": (1/3, 1/3, 1/3), "A": (1, 0, 0), "B": (0, 1, 0), "C": (0, 0, 1)}
SP = {"nonpol": "iso", "A": "A", "B": "B", "C": "C"}  # statepoint dir tags from Tier 5/6
FLIBE_ID = 5


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    eta = {m: rm.an_eta(*abc) if hasattr(rm, "an_eta") else (abc[0] + 2*abc[1]/3 + abc[2]/3)
           for m, abc in MODES.items()}
    eta_ref = eta["nonpol"]

    # per-neutron quantities from the existing Tier 5/6 runs (scattering on)
    per = {}
    for m, tag in SP.items():
        sp = f"/tmp/spf_t56_{tag}/statepoint.10.h5"
        tbr_n, _ = rm.tbr(sp, FLIBE_ID)
        inb, outb = rm.midplane_intensities(sp, 40, 40)
        per[m] = dict(tbr=tbr_n, inb=inb, outb=outb)

    # absolute = (relative neutron yield = eta/eta_ref) x per-neutron, normalized to nonpol=1
    rel = {m: eta[m] / eta_ref for m in MODES}
    inb0 = per["nonpol"]["inb"]; outb0 = per["nonpol"]["outb"]; tbr0 = per["nonpol"]["tbr"]
    table = {}
    for m in MODES:
        table[m] = dict(
            rate=rel[m],                                   # fusion power / total wall load
            tritium=rel[m] * per[m]["tbr"] / tbr0,         # absolute tritium production
            inboard=rel[m] * per[m]["inb"] / inb0,         # absolute center-stack load
            outboard=rel[m] * per[m]["outb"] / outb0,      # absolute outboard load
            tbr_per_n=per[m]["tbr"],                       # the self-sufficiency ratio (unchanged)
        )

    # ---- figure: grouped bars vs unpolarized = 1 ----
    ms = ["nonpol", "A", "B", "C"]
    mlabel = {"nonpol": "Unpolarized", "A": "A", "B": "B", "C": "C"}
    metrics = [("rate", "Fusion Rate / Total Wall Load"), ("tritium", "Absolute Tritium Production"),
               ("inboard", "Center-Stack (Inboard) Load")]
    x = np.arange(len(ms)); wbar = 0.25
    hatches = ["///", "...", "xxx"]   # distinguish series in black & white
    fig, ax = plt.subplots(figsize=(9.5, 5))
    for i, (k, lab) in enumerate(metrics):
        ax.bar(x + (i - 1) * wbar, [table[m][k] for m in ms], wbar, label=lab,
               facecolor="white", edgecolor="black", linewidth=1.1, hatch=hatches[i])
    ax.axhline(1.0, color="black", ls="dashed", lw=1)
    ax.set_xticks(x); ax.set_xticklabels([mlabel[m] for m in ms], fontsize=13)
    ax.set_ylabel("Relative to Unpolarized Fuel (= 1)", fontsize=14)
    ax.set_title("Absolute Effect at Fixed Fuel Density (Rate Restored)", fontsize=16)
    ax.legend(fontsize=12)
    for i, m in enumerate(ms):  # annotate rate
        ax.text(x[i] - wbar, table[m]["rate"] + 0.02, f"{table[m]['rate']:.2f}×",
                ha="center", fontsize=10)
    fig.tight_layout(); fig.savefig(FIGDIR / "tier7_absolute.png", dpi=130); plt.close(fig)

    # ---- RESULTS_tier7.md ----
    L = ["# RESULTS - Tier 7 (absolute rate effect: rate restored)\n",
         "Tiers 2b/3/6 held the fusion rate constant to isolate *directionality*. Here we "
         "restore it: at fixed fuel density the neutron yield scales with the total-rate "
         "factor eta = a + 2b/3 + c/3. Pure postprocessing of the Tier 5/6 statepoints "
         "(scattering on); all values relative to unpolarized fuel = 1.\n",
         "## The rate factor eta (the +50% / -50% effect)\n",
         "| mode | eta | fusion rate vs unpolarized |", "|---|---|---|"]
    for m in ms:
        L.append(f"| {m} | {eta[m]:.3f} | {table[m]['rate']:.2f}x ({100*(table[m]['rate']-1):+.0f}%) |")
    L += ["\n**A mode = +50% fusion power; C mode = -50%** -- the cross-section enhancement "
          "(Kulsrud 1982), which our eta carried all along but which the constant-rate "
          "comparisons normalized away.\n",
          "## Absolute quantities at fixed fuel density (vs unpolarized = 1)\n",
          "| mode | fusion rate | abs. tritium prod. | center-stack load | outboard load | TBR per neutron |",
          "|---|---|---|---|---|---|"]
    for m in ms:
        t = table[m]
        L.append(f"| {m} | {t['rate']:.2f}x | {t['tritium']:.2f}x | {t['inboard']:.2f}x | "
                 f"{t['outboard']:.2f}x | {t['tbr_per_n']:.3f} |")
    L += ["\n![absolute](figs/tier7_absolute.png)\n",
          "## What this shows (the rate vs steering trade-off)\n",
          f"- **A (rate play):** +50% power and +{100*(table['A']['tritium']-1):.0f}% absolute "
          f"tritium -- but it steers neutrons *inboard*, so the center-stack load compounds to "
          f"**{table['A']['inboard']:.2f}x** (rate x inboard steering). Great for power, hard on the center stack.",
          f"- **B (free steering):** same power, same breeding, center-stack load **{table['B']['inboard']:.2f}x** "
          f"-- the sweet spot: directional benefit at no rate penalty.",
          f"- **C (steering at a cost):** halves power and breeding ({table['C']['rate']:.2f}x), giving the "
          f"lowest absolute center-stack load ({table['C']['inboard']:.2f}x) -- only worth it if center-stack "
          f"load is the binding constraint.",
          f"- **Self-sufficiency is preserved:** the TBR *per neutron* is ~mode-independent "
          f"({min(table[m]['tbr_per_n'] for m in ms):.3f}-{max(table[m]['tbr_per_n'] for m in ms):.3f}), so the "
          f"rate boost breeds proportionally more tritium -- the rate gain does not cost breeding adequacy.\n",
          "## Scope: where the *other* headline SPF benefits live (NOT modeled here)\n",
          "- The cross-section **+50%** is in our model (eta, above).",
          "- **Net-electricity / lower-density gains** (disproportionate Q, recirculating power) "
          "are a plasma + **systems** power-balance result -- upstream of neutron transport.",
          "- **~10x tritium burn efficiency / inventory** is a plasma **fuel-cycle** (particle-balance) "
          "effect -- distinct from our blanket TBR (breeding supply vs burn demand). Both feed tritium "
          "self-sufficiency; we model only the blanket-breeding half.\n",
          "**Tier-7 takeaway: the rate boost is real (A +50%) but lives in a different spin state than "
          "the center-stack steering (B/C); B uniquely buys steering at no rate cost.**"]
    (REPO / "spf_prototype" / "RESULTS_tier7.md").write_text("\n".join(L) + "\n")

    print("wrote RESULTS_tier7.md, figs/tier7_absolute.png")
    for m in ms:
        t = table[m]
        print(f"  {m:7s} rate={t['rate']:.2f}x tritium={t['tritium']:.2f}x "
              f"inboard={t['inboard']:.2f}x TBR/n={t['tbr_per_n']:.3f}")


if __name__ == "__main__":
    main()
