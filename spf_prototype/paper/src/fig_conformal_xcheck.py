"""Fig: free-streaming cross-check ON the conformal wall (device 59509). Analytic point-patch
NWL directionality vs the ray-traced OpenMC StellaratorSource, per wall bin, on the real
conformal geometry -- they agree to 0.5-2.3% by region (per-bin corr ~0.91)."""
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

f = np.load(DATA / "conformal_maps" / "quasr59509_conformal_xcheck.npz")
ok = f["ok"]
A_an, A_mc = f["A_iso_an"][ok], f["A_iso_mc"][ok]
B_an, B_mc = f["B_iso_an"][ok], f["B_iso_mc"][ok]

fig, ax = plt.subplots(figsize=(4.5, 4.3))
lim = [0.6, 1.5]
ax.plot(lim, lim, "-", color="0.7", lw=0.9, zorder=0)
ax.scatter(A_an, A_mc, s=10, facecolor="none", edgecolor="0.15", lw=0.5, label="A (perp)")
ax.scatter(B_an, B_mc, s=10, marker="s", facecolor="none", edgecolor="0.55", lw=0.5, label="B/C (par)")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("analytic  load / isotropic")
ax.set_ylabel("OpenMC  load / isotropic")
ax.set_title("Conformal-wall free-streaming check (59509)", fontsize=9.5)
ax.legend(frameon=False, loc="upper left", fontsize=8)
ax.text(0.97, 0.06, "per-bin r = 0.91\nregions agree 0.5-2.3%", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=7.5, color="0.35")
finish(fig, "fig_conformal_xcheck")
