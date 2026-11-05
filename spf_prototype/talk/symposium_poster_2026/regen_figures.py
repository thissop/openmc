"""Regenerate the two non-conforming HBFS figures clean (smplotlib serif,
legend-only, no title, no in-plot annotation) for the poster. Real data only."""
import os, glob, re
import numpy as np
import smplotlib  # noqa: F401  (SM serif house style)
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype"
CP = f"{ROOT}/shield_opt/patch_worth/cluster_pull"
OUT = "/private/tmp/claude-501/-Users-tkiker-Documents-GitHub-openmc/638db22c-edea-421d-88fe-e7f5d5ecfd10/scratchpad/poster_figs"
os.makedirs(OUT, exist_ok=True)

# ---------- HARDER: per-coil fast-neutron flux (uniform shield), 13.2x swing ----------
base = np.load(f"{CP}/coil_step1b_uniform.npz")
flux = np.ravel(base["coil_fast_flux"]).astype(float)
ferr = np.ravel(base["coil_fast_flux_err"]).astype(float)
idx = np.arange(11, 11 + len(flux))           # coil indices 11..30
hot = int(np.argmax(flux))
peaking = flux.max() / flux.min()

fig, ax = plt.subplots(figsize=(6.6, 4.6))
colors = ["0.72"] * len(flux)
colors[hot] = "black"
ax.bar(idx, flux, yerr=ferr, color=colors, edgecolor="black", linewidth=0.7,
       error_kw=dict(ecolor="black", elinewidth=0.8, capsize=2))
ax.set_xlabel("Coil Index")
ax.set_ylabel("Fast-Neutron Flux per Coil")
ax.set_xticks(idx[::2])
ax.margins(x=0.01)
ax.legend(handles=[Patch(facecolor="black", edgecolor="black", label="Hottest Coil"),
                   Patch(facecolor="0.72", edgecolor="black", label="Other Coils")],
          loc="upper right", frameon=False)
fig.tight_layout()
fig.savefig(f"{OUT}/harder_coil_peaking.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"HARDER: {len(flux)} coils, peaking {peaking:.1f}x, hottest coil {idx[hot]}")

# ---------- BETTER: peak coil load vs uniform baseline (sign is the robust finding) ----------
patches = ["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
peak_base = flux.max()
ratios = []
for p in patches:
    R = np.ravel(np.load(f"{CP}/coil_patch_{p}.npz")["coil_fast_flux"]).astype(float)
    ratios.append(R.max() / peak_base)
ratios = np.sort(np.array(ratios))

fig, ax = plt.subplots(figsize=(6.6, 4.6))
x = np.arange(len(ratios))
for xi, r in zip(x, ratios):
    if r < 1.0:
        ax.bar(xi, r, color="0.72", edgecolor="black", linewidth=0.7)
    else:
        ax.bar(xi, r, color="0.30", edgecolor="black", linewidth=0.7, hatch="//")
ax.axhline(1.0, color="black", linestyle="--", linewidth=1.2)
ax.set_ylabel("Peak Coil Load  (× Uniform Baseline)")
ax.set_xlabel("Shield Placement Candidate")
ax.set_xticks([])
ax.set_ylim(0, 1.65)
ax.legend(handles=[Patch(facecolor="0.72", edgecolor="black", label="Reduces Peak"),
                   Patch(facecolor="0.30", edgecolor="black", hatch="//", label="Backfires"),
                   Line2D([0],[0], color="black", linestyle="--", label="Uniform Baseline")],
          loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(f"{OUT}/better_placed_vs_uniform.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"BETTER: {len(ratios)} candidates, {int((ratios<1).sum())} reduce "
      f"({ratios.min():.2f}-{ratios[ratios<1].max():.2f}x), "
      f"{int((ratios>=1).sum())} backfire ({ratios[ratios>=1].min():.2f}-{ratios.max():.2f}x)")
print("wrote", OUT)
