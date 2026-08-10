"""Color the real QH patch grid by the MEASURED finite-difference worth from our runs.
value[pp,pt] = W_FD[i_tor=pt%4, j_pol=(pp+2)%4]; missing patches (undosed) stay grey.
Diverging scale centered at 0: red = shielding backfires (dose up), blue = dose down."""
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa
sys.path.insert(0, os.path.join("..", "shield_opt", "patch_worth"))
import patch_worth_analyze as pwa
from gen_patch_fig import draw

P = "../shield_opt/patch_worth/cluster_pull"
R0 = pwa._coil_R(os.path.join(P, "coil_step1b_uniform.npz"), 20)
wfd = np.full((4, 4), np.nan)            # [i_tor, j_pol]
for p in sorted(glob.glob(os.path.join(P, "coil_patch_*.npz"))):
    tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
    i = int(tag.split("t")[1].split("p")[0]); j = int(tag.split("p")[1])
    dd = np.load(os.path.join(P, f"delta_qh_patch_{tag}.npz"))
    added = float(dd["added_cells"]) * float(dd["delta_cm"])
    wfd[i, j] = -(pwa._coil_R(p, 20) - R0) / added

# map experiment (i_tor, j_pol) -> viz grid (pp in 0..3 poloidal, pt in 0..15 toroidal)
NPOL, NTOR = 4, 16
val = np.array([[wfd[pt % 4, (pp + 2) % 4] for pt in range(NTOR)] for pp in range(NPOL)])
print("W_FD grid [i_tor rows, j_pol cols] (x1e-5):")
print(np.array2string(wfd*1e5, precision=1))

vmax = np.nanmax(np.abs(wfd))
norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

fig = plt.figure(figsize=(9.2, 9))
ax = fig.add_axes([0, 0.06, 1, 0.94], projection="3d")
draw(ax, patch_values=val, cmap="RdBu", norm=norm)   # RdBu: negative=red (backfire), positive=blue
ax.view_init(elev=34, azim=-58)
# slim colorbar
sm = plt.cm.ScalarMappable(norm=norm, cmap="RdBu"); sm.set_array([])
cax = fig.add_axes([0.30, 0.05, 0.40, 0.022])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("measured worth  (red = shielding backfires)", fontsize=10)
cb.ax.tick_params(labelsize=8)
fig.savefig("figs/patch_geometry_worth.png", dpi=150, bbox_inches="tight", pad_inches=0.05)
fig.savefig("figs/patch_geometry_worth.pdf", bbox_inches="tight", pad_inches=0.05)
print("saved figs/patch_geometry_worth.png/.pdf")
