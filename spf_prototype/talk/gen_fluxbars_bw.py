"""Black-and-white Full-vs-SSR coil-load figure (the only figure in the deck).
Full Run = white bars, black outline. SSR = diagonal hatch. Each series scaled to its own peak coil
(shape comparison; the benign ~2.4x normalization divides out). No shading, no in-plot annotations."""
import numpy as np
import smplotlib  # noqa: F401
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = "/Users/tkiker/.claude/jobs/876e30be/tmp"
FIG = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/figs"
M = np.load(f"{T}/coil_real_master.npz"); P = np.load(f"{T}/coil_real_pristine.npz")
cells = M["coil_cells"]
fn = M["coil_fast_flux"] / M["coil_fast_flux"].max()
sn = P["coil_fast_flux"] / P["coil_fast_flux"].max()

fig, ax = plt.subplots(figsize=(8.6, 4.5))
x = np.arange(len(cells)); w = 0.42
ax.bar(x - w/2, fn, w, facecolor="white", edgecolor="black", linewidth=1.2, label="Full Run")
ax.bar(x + w/2, sn, w, facecolor="white", edgecolor="black", linewidth=1.2, hatch="////", label="SSR")
ax.set_xticks(x); ax.set_xticklabels([int(c) for c in cells], fontsize=8)
ax.set_xlabel("Coil Cell ID"); ax.set_ylabel("Coil Load (Relative To Peak Coil)")
ax.set_title("Coil Load Shape: Full Run vs SSR")
ax.legend(fontsize=14, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(f"{FIG}/ci_fluxbars.pdf"); plt.close(fig)
print("WROTE ci_fluxbars.pdf")
