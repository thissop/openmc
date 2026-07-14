"""3D overview of the QH baseline device: the 20 real Wiedman coils around the
LandremanPaul-QH plasma (reactor scale). Diagnostic, not a paper figure."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import coil_adapter as ca

WOUT = HERE / "data" / "raw" / "QH_unzipped" / "zenodo" / "configurations" / "LandremanPaulQH_coils" / "wout_vmec.nc"

# --- coils (reactor scale, already ~14 m) ---
cs = ca.load_qh()
print(f"{cs.n_coils} coils, R0~{cs.extent()['R0']:.1f} m")

# --- plasma LCFS from the VMEC wout ---
from simsopt.geo import SurfaceRZFourier
surf = SurfaceRZFourier.from_wout(str(WOUT), s=1.0,
                                  quadpoints_phi=np.linspace(0, 1, 160),
                                  quadpoints_theta=np.linspace(0, 1, 80))
g = surf.gamma()          # (nphi, ntheta, 3), meters
X, Y, Z = g[..., 0], g[..., 1], g[..., 2]

fig = plt.figure(figsize=(11, 5.2))

# left: full device
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_surface(X, Y, Z, color="#6a9fd8", alpha=0.35, linewidth=0,
                rstride=2, cstride=2, shade=True)
for c, cur in zip(cs.centerlines, cs.currents):
    P = np.vstack([c, c[:1]])   # close the loop
    ax.plot(P[:, 0], P[:, 1], P[:, 2], color="#b5622c", lw=1.6)
ax.set_title("QH Device: 20 Coils + Plasma", fontsize=11)
ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]"); ax.set_zlabel("Z [m]")
ax.view_init(elev=32, azim=-52)
try: ax.set_box_aspect((1, 1, 0.55))
except Exception: pass

# right: one field period, top-down, showing coil-plasma standoff
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
nph = X.shape[0]; q = nph // 4          # one of nfp=4 periods
ax2.plot_surface(X[:q+1], Y[:q+1], Z[:q+1], color="#6a9fd8", alpha=0.45,
                 linewidth=0, rstride=2, cstride=2, shade=True)
phi_c = np.array([np.arctan2(c[:, 1], c[:, 0]).mean() for c in cs.centerlines])
for c in [cs.centerlines[i] for i in range(cs.n_coils) if 0 <= phi_c[i] <= np.pi/2 + 0.3]:
    P = np.vstack([c, c[:1]])
    ax2.plot(P[:, 0], P[:, 1], P[:, 2], color="#b5622c", lw=1.8)
ax2.set_title("One Field Period (nfp=4)", fontsize=11)
ax2.set_xlabel("X [m]"); ax2.set_ylabel("Y [m]"); ax2.set_zlabel("Z [m]")
ax2.view_init(elev=68, azim=-60)
try: ax2.set_box_aspect((1, 1, 0.5))
except Exception: pass

fig.tight_layout()
out = HERE / "figs" / "qh_device_3d"
fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{out}.pdf", bbox_inches="tight")
print(f"wrote {out}.png")
