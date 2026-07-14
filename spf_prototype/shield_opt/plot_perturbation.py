"""Intuitive visualization of the fixed-envelope shield<->breeder trade, on the real QH
plasma. Left: a poloidal cross-section showing the shield growing INWARD into the FLiBe over
a targeted region while the OUTER envelope (and coils) stay put. Right: the (theta,phi)
control field delta -- where and how much we thicken shield / thin breeder."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import thickness_field as tf
import attenuation_surrogate as asur

WOUT = HERE / "data" / "raw" / "QH_unzipped" / "zenodo" / "configurations" / "LandremanPaulQH_coils" / "wout_vmec.nc"

# baseline radial build (cm), plasma-facing -> out
LAYERS = [("first_wall", 4.0), ("breeder", 50.0), ("back_wall", 3.0),
          ("shield", 40.0), ("vacuum_vessel", 10.0)]
T_BREEDER0, T_SHIELD0 = 50.0, 40.0

# --- LCFS cross-section at one toroidal angle from the VMEC wout ---
from simsopt.geo import SurfaceRZFourier
nph, nth = 64, 160
surf = SurfaceRZFourier.from_wout(str(WOUT), s=1.0,
                                  quadpoints_phi=np.linspace(0, 1, nph),
                                  quadpoints_theta=np.linspace(0, 1, nth))
g = surf.gamma()                       # (nph, nth, 3) meters
jphi = 0                               # toroidal cut
R = np.hypot(g[jphi, :, 0], g[jphi, :, 1]) * 100.0   # -> cm
Z = g[jphi, :, 2] * 100.0
theta = np.linspace(0, 2 * np.pi, nth, endpoint=False)

# outward poloidal normals for this single loop
def outward_normals(R, Z):
    tR = np.gradient(R); tZ = np.gradient(Z)
    nR, nZ = tZ.copy(), -tR.copy()                    # rotate tangent -90
    L = np.hypot(nR, nZ) + 1e-30; nR /= L; nZ /= L
    cx, cy = R.mean(), Z.mean()
    flip = ((R - cx) * nR + (Z - cy) * nZ) < 0        # orient away from centroid
    nR[flip] *= -1; nZ[flip] *= -1
    return nR, nZ
nR, nZ = outward_normals(R, Z)

# --- the control field: a smooth bump on the inboard region (theta~180 deg) ---
field = tf.ThicknessField(nfp=4, toroidal_angles_deg=np.linspace(0, 90, 9),
                          poloidal_angles_deg=np.degrees(theta),
                          t_breeder0=T_BREEDER0, t_shield0=T_SHIELD0, t_breeder_min=15.0)
delta_grid = field.gaussian_bump(theta0_deg=180.0, phi0_deg=0.0, amp=30.0, width_deg=30.0)
delta_theta = delta_grid[0]            # the phi=0 poloidal slice, delta(theta)

def boundaries(delta_th):
    """Return list of (R,Z) closed loops for each layer interface, plasma-out."""
    tb = T_BREEDER0 - delta_th
    ts = T_SHIELD0 + delta_th
    thick = {"first_wall": np.full(nth, 4.0), "breeder": tb,
             "back_wall": np.full(nth, 3.0), "shield": ts,
             "vacuum_vessel": np.full(nth, 10.0)}
    cum = np.zeros(nth); loops = [(R.copy(), Z.copy())]
    for name, _ in LAYERS:
        cum = cum + thick[name]
        loops.append((R + nR * cum, Z + nZ * cum))
    return loops

base = boundaries(np.zeros(nth))
pert = boundaries(delta_theta)

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(11.4, 5.0))

# LEFT: stacked radial build vs poloidal angle (no offset-folding). Depth from the
# plasma-facing surface outward; each band a layer. FLiBe pinches, shield bulges, at
# fixed total depth (the outer envelope is a flat line).
COL = {"first_wall": "#888888", "breeder": "#38b8bf", "back_wall": "#bfbfbf",
       "shield": "#e8973a", "vacuum_vessel": "#7f9fc4"}
thd = np.degrees(theta)
order = np.argsort(thd); thd_s = thd[order]
tb = (T_BREEDER0 - delta_theta)[order]
ts = (T_SHIELD0 + delta_theta)[order]
bands = [("first_wall", np.full(nth, 4.0)[order]), ("breeder", tb),
         ("back_wall", np.full(nth, 3.0)[order]), ("shield", ts),
         ("vacuum_vessel", np.full(nth, 10.0)[order])]
low = np.zeros(nth)
for name, th in bands:
    up = low + th
    ax[0].fill_between(thd_s, low, up, color=COL[name], alpha=0.9, label=name)
    low = up
# baseline breeder/shield interface (fixed layers) as a dashed guide
base_interface = np.full(nth, 4.0 + T_BREEDER0)   # fw + baseline breeder
ax[0].plot(thd_s, base_interface, "k--", lw=1.1, label="baseline FLiBe/shield interface")
ax[0].plot(thd_s, low, color="0.2", lw=1.3, label="fixed outer envelope")
ax[0].set_xlim(0, 360); ax[0].set_ylim(0, 112)
ax[0].set_xlabel("Poloidal Angle $\\theta$ [deg]")
ax[0].set_ylabel("Radial Depth From First Wall [cm]")
ax[0].set_title("Radial Build At $\\phi=0$: FLiBe$\\to$Shield Trade, Envelope Fixed", fontsize=10)
ax[0].legend(loc="lower center", fontsize=6.5, ncol=3, frameon=False)

# RIGHT: (theta,phi) control field delta
im = ax[1].imshow(delta_grid.T, origin="lower", aspect="auto", cmap="RdBu_r",
                  vmin=-delta_grid.max(), vmax=delta_grid.max(),
                  extent=[0, 90, 0, 360])
ax[1].set_xlabel("Toroidal Angle $\\phi$ [deg]")
ax[1].set_ylabel("Poloidal Angle $\\theta$ [deg]")
ax[1].set_title("Control Field: +Shield / -FLiBe [cm]", fontsize=10)
cb = fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.02)
cb.set_label("$\\delta t$ [cm]  (red = +shield, -FLiBe)")

# surrogate estimate of the coil-dose benefit in the bump + TBR note
k = asur.trade_k()
dose_ratio = np.exp(-k * delta_theta.max())
fig.text(0.5, 0.005,
         f"Surrogate at peak $\\delta t$={delta_theta.max():.0f} cm: coil fast-dose "
         f"$\\times${dose_ratio:.2f}  (k={k:.3f}/cm).  Envelope conserved; "
         f"FLiBe locally {T_BREEDER0:.0f}$\\to${(T_BREEDER0-delta_theta.max()):.0f} cm "
         f"(TBR cost -- see note).", ha="center", fontsize=8)

fig.tight_layout(rect=(0, 0.03, 1, 1))
out = HERE / "figs" / "perturbation_demo"
fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{out}.pdf", bbox_inches="tight")
print(f"wrote {out}.png  | peak delta {delta_theta.max():.1f} cm, coil dose x{dose_ratio:.2f}")
