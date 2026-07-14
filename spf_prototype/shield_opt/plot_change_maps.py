"""Two (theta,phi) change maps for the fixed-envelope trade, with the requested convention:
white = no change, BLUE = material ADDED, RED = material REMOVED. Left = neutron shield,
right = FLiBe breeder (mirror image, since the trade is volume-neutral). Uses a signed
control field so both add and remove show -- the optimizer can thicken shield over a hot
region AND recover breeder elsewhere."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nfp = 4
tor = np.linspace(0, 90, 91)          # one field period, deg
pol = np.linspace(0, 360, 181)        # deg
TOR, POL = np.meshgrid(np.radians(tor), np.radians(pol), indexing="ij")

def bump(theta0, phi0, amp, w_deg=28.0):
    dth = (POL - np.radians(theta0) + np.pi) % (2*np.pi) - np.pi
    dph = ((nfp*TOR - nfp*np.radians(phi0) + np.pi) % (2*np.pi) - np.pi) / nfp
    return amp * np.exp(-(dth**2 + dph**2) / (2*np.radians(w_deg)**2))

# signed control field: +shield at inboard hotspot (theta~180), recover a bit at outboard
delta_shield = bump(180, 0, +25.0) - bump(0, 45, 12.0, w_deg=32)   # cm added to shield
delta_flibe = -delta_shield                                         # removed from FLiBe

vmax = np.abs(delta_shield).max()
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6), sharey=True)
ext = [0, 90, 0, 360]
for a, field, title in [(ax[0], delta_shield, "Neutron Shield"),
                        (ax[1], delta_flibe, "FLiBe Breeder")]:
    im = a.imshow(field.T, origin="lower", aspect="auto", cmap="RdBu",
                  vmin=-vmax, vmax=vmax, extent=ext)
    a.set_xlabel("Toroidal Angle $\\phi$ [deg]")
    a.set_title(f"{title}: Added (Blue) / Removed (Red)", fontsize=11)
ax[0].set_ylabel("Poloidal Angle $\\theta$ [deg]")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cb.set_label("$\\Delta$ Thickness [cm]   (blue = added, red = removed)")
fig.suptitle("Fixed-Envelope Trade: Where Shield Is Added, FLiBe Is Removed (And Vice Versa)",
             fontsize=11)
out = "figs/shield_flibe_change_maps"
fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{out}.pdf", bbox_inches="tight")
print(f"wrote {out}.png  | shield add max {delta_shield.max():.0f} cm, remove {delta_shield.min():.0f} cm")
