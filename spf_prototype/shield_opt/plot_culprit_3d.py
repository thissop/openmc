"""3D version of the culprit map: the plasma LCFS colored by each region's CONTRIBUTION to
the peak first-wall NWL, for unpolarized vs SPF A-mode (device 59509). Complements the 2D
(theta,zeta) culprit maps -- same data (culprit_59509.npz), rendered on the boundary surface."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

f = np.load("data/processed/culprit_59509.npz", allow_pickle=True)
nrho, nt, nz = int(f["nrho"]), int(f["nt"]), int(f["nz"])
Sv = f["Sv"].reshape(nrho, nt, nz, 3)          # source voxels (x,y,z)
lcfs = Sv[-1]                                    # outermost rho = plasma boundary (nt,nz,3)
X, Y, Z = lcfs[..., 0], lcfs[..., 1], lcfs[..., 2]

fig = plt.figure(figsize=(12, 5.4))
for k, (key, label) in enumerate([("Cmap_unpol", "Unpolarized"),
                                  ("Cmap_A", "SPF A-Mode (Perp)")]):
    C = np.asarray(f[key])                       # (nt, nz) contribution to peak NWL
    # normalize both panels to a common scale for honest comparison
    if k == 0:
        vmax = max(np.asarray(f["Cmap_unpol"]).max(), np.asarray(f["Cmap_A"]).max())
    norm = colors.Normalize(0, vmax)
    fc = cm.inferno(norm(C))
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    ax.plot_surface(X, Y, Z, facecolors=fc, rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False)
    ax.set_title(f"{label}: Plasma Contribution To Peak NWL", fontsize=10)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.view_init(elev=26, azim=-60)
    try: ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
    except Exception: pass

m = cm.ScalarMappable(norm=norm, cmap="inferno")
m.set_array([])
cb = fig.colorbar(m, ax=fig.axes, fraction=0.025, pad=0.04)
cb.set_label("Contribution To Peak NWL  (bright = culprit region)")
out = "figs/culprit_3d_59509"
fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{out}.pdf", bbox_inches="tight")
print(f"wrote {out}.png  | LCFS {X.shape}, vmax {vmax:.2e}")
