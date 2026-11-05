"""(phi, theta) grid figures for the 'meshing capacity' question.
LEFT/MID: the REAL 4x4 attribution A(phi,theta) and worth W(phi,theta) we currently work with
          (from patch_worth_tier1.csv) -- capped near 4x4 by the coarse 30x30x20 transport mesh
          (120 cm voxels vs a 40 cm shield band).
RIGHT:    the resolution ladder the FROZEN SURFACE SOURCE unlocks -- exact neutron positions bin to
          any (phi,theta) grid, statistics-limited (56.7 M crossings), not mesh-limited.
"""
import csv, numpy as np
import smplotlib  # noqa: F401
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/shield_opt/patch_worth/cluster_pull/patch_worth_tier1.csv"
FIG = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/figs"

A = np.zeros((4, 4)); W = np.zeros((4, 4))
with open(CSV) as f:
    for r in csv.DictReader(f):
        i = int(r["i_tor"]); j = int(r["j_pol"])
        A[i, j] = float(r["A_attribution"]); W[i, j] = float(r["W_worth"])

phi_edges = np.linspace(0, 90, 5)          # one field period, nfp=4
th_edges = np.linspace(-180, 180, 5)

fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
for ax, M, ttl, cmap in [(axes[0], A.T, "Attribution  $A=S\\,\\psi^\\dagger$  (where coil-loading neutrons are born)", "viridis"),
                         (axes[1], W.T, "Worth  $W=\\sum_{\\rm shield}\\phi\\,\\psi^\\dagger$  (where shielding blocks flux)", "magma")]:
    pc = ax.pcolormesh(phi_edges, th_edges, M, cmap=cmap, shading="flat", edgecolors="k", linewidth=0.5)
    ax.set_xlabel(r"toroidal $\phi$ (deg, one period)"); ax.set_ylabel(r"poloidal $\theta$ (deg)")
    ax.set_title(ttl, fontsize=8.5)
    ax.set_xticks([0, 45, 90]); ax.set_yticks([-180, -90, 0, 90, 180])
    fig.colorbar(pc, ax=ax, fraction=0.046)

# right: resolution ladder (empty grids at increasing density over the same domain)
ax = axes[2]
for (nt, npol, c, lw, lab) in [(4, 4, "#b2182b", 2.2, "4x4  (current, mesh-limited)"),
                               (13, 19, "#1a599f", 1.1, "13x19"),
                               (60, 90, "#2a7f4f", 0.35, "60x90  (frozen source)")]:
    for x in np.linspace(0, 90, nt+1):
        ax.plot([x, x], [-180, 180], color=c, lw=lw, alpha=0.9)
    for y in np.linspace(-180, 180, npol+1):
        ax.plot([0, 90], [y, y], color=c, lw=lw, alpha=0.9)
ax.set_xlim(0, 90); ax.set_ylim(-180, 180)
ax.set_xlabel(r"toroidal $\phi$ (deg)"); ax.set_ylabel(r"poloidal $\theta$ (deg)")
ax.set_title("Resolution unlocked by the frozen surface source", fontsize=8.5)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],color="#b2182b",lw=2.2,label="4x4  current (mesh-limited)"),
                   Line2D([0],[0],color="#1a599f",lw=1.1,label="13x19"),
                   Line2D([0],[0],color="#2a7f4f",lw=0.6,label="60x90  frozen source")],
          fontsize=7.5, frameon=False, loc="upper right")
fig.suptitle("(phi, theta) shielding grids: what we work with now, and what surface-freezing unlocks", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/grids_phitheta.pdf"); plt.close(fig)
print("WROTE grids_phitheta.pdf  |  A range", A.min(), A.max(), " W range", W.min(), W.max())
