"""Fig: volumetric flux tally on the conformal blanket (QUASR 59509, scattering on).
Left: poloidal fast-flux map (phi-averaged), grayscale. Right: heating per blanket layer.
Native StellaratorSource through the multi-layer DAGMC blanket."""
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

f = np.load(DATA / "conformal_maps" / "quasr59509_conformal_flux.npz", allow_pickle=True)
r, z = np.asarray(f["r_grid"]), np.asarray(f["z_grid"])
layers = [str(x) for x in f["layers"]]
fmap = f["flux_unpolarized"].mean(axis=1)          # (nr, nz), phi-averaged fast flux
heat = f["heat_unpolarized"]

fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.5), gridspec_kw={"width_ratios": [1.15, 1]})

im = ax[0].imshow(fmap.T, origin="lower", aspect="auto", cmap="gray_r",
                  extent=[r[0], r[-1], z[0], z[-1]])
ax[0].set_xlabel("Major Radius R [cm]")
ax[0].set_ylabel("Height Z [cm]")
cb = fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.02)
cb.set_label("Fast Flux [/cm$^2$ per Source]")

xb = np.arange(len(layers))
ax[1].bar(xb, heat, color="none", edgecolor="k", width=0.7)   # black outline, no fill
ax[1].set_yscale("log")
_LBL = {"W": "W", "steel": "Steel", "Be": "Be", "FLiBe": "FLiBe", "shield": "Shield", "coil": "Coil"}
ax[1].set_xticks(xb)
ax[1].set_xticklabels([_LBL.get(l, l) for l in layers], rotation=30, ha="right")
ax[1].set_ylabel("Heating [eV per Source]")

finish(fig, "fig_conformal_flux")
