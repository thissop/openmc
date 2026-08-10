import os, sys, glob
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "..")
import smplotlib  # noqa
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import patch_worth_analyze as pwa
from scipy.stats import spearmanr

P = "cluster_pull"
FIG = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/figs"
KEY = "#1a599f"; WARN = "#b2182b"; GREY = "#7a7a7a"

def coilmax(f):
    d = np.load(os.path.join(P, f)); return float(np.max(d["coil_fast_flux"]))

# ---- Figure 1: peak coil fast flux, three placements ----
labels = ["Uniform", "Single-Coil\nPlacement", "Distributed\nPlacement"]
vals = [coilmax("coil_step1b_uniform.npz"), coilmax("coil_step1b_placed.npz"),
        coilmax("coil_step1b_multicoil.npz")]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
bars = ax.bar(labels, vals, color=[GREY, KEY, KEY], width=0.6, edgecolor="k", linewidth=1.1)
u = vals[0]
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v*1.02, f"{v/u:.2f}x", ha="center", va="bottom")
ax.set_ylabel("Peak Coil Fast Flux (n/cm$^2$ per source)")
ax.set_title("Peak Coil Flux at Equal Shield Material")
fig.tight_layout(); fig.savefig(f"{FIG}/peak_flux.pdf"); plt.close(fig)

# ---- worth data: W_FD, W scalar, A attribution per patch ----
key = pwa.load_tier1(os.path.join(P, "patch_worth_tier1.csv"))
R0 = pwa._coil_R(os.path.join(P, "coil_step1b_uniform.npz"), 20)
rows = []
for p in sorted(glob.glob(os.path.join(P, "coil_patch_*.npz"))):
    tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
    i = int(tag.split("t")[1].split("p")[0]); j = int(tag.split("p")[1])
    Rp = pwa._coil_R(p, 20)
    dd = np.load(os.path.join(P, f"delta_qh_patch_{tag}.npz"))
    added = float(dd["added_cells"]) * float(dd["delta_cm"])
    A, W = key[(i, j)]
    rows.append(dict(W_FD=-(Rp-R0)/added, A=A, W=W))
wfd = np.array([r["W_FD"] for r in rows])
Wsc = np.array([r["W"] for r in rows]); A = np.array([r["A"] for r in rows])

# ---- Figure 2: proxy vs finite-difference truth (two panels) ----
fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.2))
def panel(a, x, y, xl, rho):
    back = y < 0
    a.scatter(x[~back], y[~back]*1e5, s=55, color=KEY, edgecolor="k", linewidth=0.5, zorder=3)
    a.scatter(x[back], y[back]*1e5, s=90, color=WARN, edgecolor="k", linewidth=0.6, marker="v",
              zorder=4, label="Shielding backfires")
    a.axhline(0, color=GREY, lw=1.0, ls="--")
    a.set_xlabel(xl); a.set_ylabel(r"Measured Worth $-\Delta$dose/material ($\times10^{-5}$)")
    a.set_title(f"Spearman $\\rho={rho:+.2f}$")
    a.legend(loc="lower left", fontsize=8, frameon=False)
panel(ax[0], Wsc/1e6, wfd, r"Predicted Scalar Worth $\phi\psi^\dagger$ ($\times10^6$)",
      spearmanr(wfd, Wsc)[0])
panel(ax[1], A/1e12, wfd, r"Attribution $S\psi^\dagger$ ($\times10^{12}$)", spearmanr(wfd, A)[0])
fig.suptitle("Cheap Adjoint Proxy Versus Monte Carlo Truth (Twelve Patches)")
fig.tight_layout(); fig.savefig(f"{FIG}/worth_vs_truth.pdf"); plt.close(fig)

# ---- Figure 3: finite-difference worth per patch, sorted ----
order = np.argsort(wfd)[::-1]
w = wfd[order]*1e5
cols = [WARN if v < 0 else KEY for v in w]
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.bar(range(len(w)), w, color=cols, edgecolor="k", linewidth=0.9)
ax.axhline(0, color="k", lw=1.0)
ax.set_xlabel("Shield Patch (Sorted)")
ax.set_ylabel(r"Measured Worth ($\times10^{-5}$)")
ax.set_title("Most Patches Help About Equally; One Backfires")
ax.set_xticks([])
fig.tight_layout(); fig.savefig(f"{FIG}/fd_worth_bars.pdf"); plt.close(fig)

print("wrote 3 figures to", FIG)
print(f"peak flux ratios: uniform=1.00 single={vals[1]/u:.2f} multi={vals[2]/u:.2f}")
print(f"Spearman W_FD vs W={spearmanr(wfd,Wsc)[0]:+.2f}  vs A={spearmanr(wfd,A)[0]:+.2f}")
