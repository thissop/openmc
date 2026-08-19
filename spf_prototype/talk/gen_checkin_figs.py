"""Figures for the 2026-08-19 Peterson check-in deck (smplotlib style, matches prior decks)."""
import numpy as np
import smplotlib  # noqa: F401
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = "/Users/tkiker/.claude/jobs/876e30be/tmp"
FIG = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/figs"
KEY = "#1a599f"; WARN = "#b2182b"; GREY = "#7a7a7a"; GOOD = "#2a7f4f"

M   = np.load(f"{T}/coil_real_master.npz")
P   = np.load(f"{T}/coil_real_pristine.npz")
VK  = np.load(f"{T}/coil_real_gate2_vk.npz")
cells = M["coil_cells"]; full = M["coil_fast_flux"]
r_pri = full / P["coil_fast_flux"]; r_vk = full / VK["coil_fast_flux"]
hot = np.array([c in range(16, 22) for c in cells])

# ---- Figure 1: per-coil ratio, void-kill bowl vs pristine (hot band collapses to ~1) ----
fig, ax = plt.subplots(figsize=(8.2, 4.6))
ax.axhline(1.0, color=GREY, ls="--", lw=1.2, zorder=1)
ax.text(30.3, 1.0, "ideal\n(k=1)", va="center", fontsize=8, color=GREY)
ax.plot(cells, r_vk, "-o", color=GREY, ms=5, lw=1.4, label="void-kill B-10 (CV 47%)", zorder=2)
ax.plot(cells, r_pri, "-o", color=KEY, ms=6, lw=1.8, label="pristine + vacuum void (CV 47%)", zorder=3)
ax.axvspan(15.5, 21.5, color=GOOD, alpha=0.10, zorder=0)
ax.text(18.5, 5.6, "inboard hot band\n(peak-load coils)\npristine CV = 4.1%, k = 1.26",
        ha="center", va="top", fontsize=8.5, color=GOOD)
ax.annotate("far-side residual\n(void-transit coils)\nstill 3.5-5x low", xy=(29, 4.5), xytext=(24.5, 5.9),
            fontsize=8.5, color=WARN, ha="center",
            arrowprops=dict(arrowstyle="->", color=WARN, lw=1.2))
ax.set_xlabel("coil cell ID"); ax.set_ylabel(r"ratio  full / SSR  (1 = faithful)")
ax.set_title("Pristine fix flattens the peak-load coils toward 1; far side still low")
ax.set_xticks(range(11, 31, 2)); ax.set_ylim(0.7, 6.3)
ax.legend(loc="upper center", fontsize=8.5, frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/ci_ratio.pdf"); plt.close(fig)

# ---- Figure 2: the frozen-source energy breakdown (58% backscatter) ----
labels = ["~14.06 MeV\n(pristine\nplasma current)", "1-14 MeV", "0.1-1 MeV", "1-100 keV", "<1 keV"]
frac = [42.0, 13.2, 22.0, 20.3, 2.6]
cols = [GOOD, WARN, WARN, WARN, WARN]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
bars = ax.bar(range(len(frac)), frac, color=cols, edgecolor="k", linewidth=1.0, width=0.7)
for b, v in zip(bars, frac):
    ax.text(b.get_x()+b.get_width()/2, v+0.6, f"{v:.0f}%", ha="center", fontsize=9)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("fraction of 56.7 M banked crossings")
ax.set_title("Naive frozen source is 58% collided backscatter (the double-count)")
ax.text(2.5, 36, "58% collided re-exits\n(dropped in the pristine fix)", ha="center",
        color=WARN, fontsize=9)
ax.text(0.0, 46, "42% = true\n1 per birth", ha="center", color=GOOD, fontsize=8.5)
ax.set_ylim(0, 52)
fig.tight_layout(); fig.savefig(f"{FIG}/ci_backscatter.pdf"); plt.close(fig)

# ---- Figure 3: full vs pristine-SSR per-coil flux (track on hot band, diverge far) ----
fig, ax = plt.subplots(figsize=(8.2, 4.4))
x = np.arange(len(cells)); w = 0.4
ax.bar(x-w/2, full*1e3, w, color=KEY, edgecolor="k", lw=0.6, label="full run (truth)")
ax.bar(x+w/2, P["coil_fast_flux"]*1e3, w, color="#8fb2d4", edgecolor="k", lw=0.6, label="pristine SSR")
ax.set_xticks(x); ax.set_xticklabels([int(c) for c in cells], fontsize=7)
ax.set_xlabel("coil cell ID"); ax.set_ylabel(r"coil fast flux ($\times10^{-3}$ per source)")
ax.set_title("Full vs pristine SSR: match at the peak, SSR low on the far side")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/ci_fluxbars.pdf"); plt.close(fig)
print("WROTE ci_ratio.pdf, ci_backscatter.pdf, ci_fluxbars.pdf")
