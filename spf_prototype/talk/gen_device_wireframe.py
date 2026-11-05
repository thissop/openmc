"""Black-and-white wireframe of the QH (Landreman-Paul) stellarator plasma boundary (LCFS),
from data/wout_qh.nc. Full torus, nfp=4. For the title-area device slide."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
WOUT = os.path.join(HERE, "..", "data", "wout_qh.nc")

import netCDF4
d = netCDF4.Dataset(WOUT)
rmnc = d.variables["rmnc"][-1, :]; zmns = d.variables["zmns"][-1, :]
xm = d.variables["xm"][:]; xn = d.variables["xn"][:]
nfp = int(d.variables["nfp"][:])

theta = np.linspace(0, 2*np.pi, 61)          # poloidal
phi = np.linspace(0, 2*np.pi, 361)           # toroidal, full torus
TH, PH = np.meshgrid(theta, phi, indexing="ij")
ang = xm[:, None, None]*TH[None] - xn[:, None, None]*PH[None]
R = np.tensordot(rmnc, np.cos(ang), axes=(0, 0))
Z = np.tensordot(zmns, np.sin(ang), axes=(0, 0))
X = R*np.cos(PH); Y = R*np.sin(PH)

fig = plt.figure(figsize=(9.5, 5.2))
ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(X, Y, Z, rstride=4, cstride=9, color="black", linewidth=0.35)
# equal aspect
mx = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()/2
cx, cy, cz = X.mean(), Y.mean(), Z.mean()
ax.set_xlim(cx-mx, cx+mx); ax.set_ylim(cy-mx, cy+mx); ax.set_zlim(cz-mx, cz+mx)
ax.set_box_aspect((1, 1, 1))
ax.view_init(elev=32, azim=-52)
ax.set_axis_off()
for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
    pane.pane.set_visible(False)
fig.tight_layout(pad=0)
fig.savefig(os.path.join(FIG, "device_wireframe.pdf"), bbox_inches="tight", pad_inches=0)
plt.close(fig)
print(f"WROTE device_wireframe.pdf  nfp={nfp}  R0~{rmnc[np.argmax((xm==0)&(xn==0))]:.2f} m")
