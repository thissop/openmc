"""Real QH stellarator: plasma boundary tiled into the 16x4 patch grid (inboard red, hottest
patch yellow) with the real coil filaments. draw(ax) builds it; run this file to save PNG/PDF,
or run view_patch_geometry.py for a draggable interactive window."""
import sys, os
import numpy as np
import netCDF4
import smplotlib  # noqa
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "coil_design"))
from read_makegrid import read_makegrid

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = "#1a599f"; RED = "#c0392b"; YEL = "#f4d03f"; GREY = "#bdbdbd"; EDGE = "#555555"

def draw(ax, patch_values=None, cmap="inferno", norm=None):
    """patch_values: optional (NPOL, NTOR) array of real per-patch data to color by (NaN = uncolored).
    None -> no coloring (just the grey wall, patch grid, and coils). norm: optional matplotlib norm
    (e.g. TwoSlopeNorm centered at 0 for a signed worth); default is linear min..max."""
    d = netCDF4.Dataset(os.path.join(HERE, "geom/wout_qh.nc"))
    xm = d.variables["xm"][:]; xn = d.variables["xn"][:]
    rmnc = d.variables["rmnc"][-1, :]; zmns = d.variables["zmns"][-1, :]
    nfp = int(d.variables["nfp"][:]); R0 = float(rmnc[np.argmax((xm == 0) & (xn == 0))])

    def boundary(theta, phi):
        ang = np.multiply.outer(xm, theta)[..., None] - np.multiply.outer(xn, phi)[:, None, :]
        return (np.tensordot(rmnc, np.cos(ang), axes=(0, 0)),
                np.tensordot(zmns, np.sin(ang), axes=(0, 0)))

    NTOR, NPOL, FT, FP = 4*nfp, 4, 4, 8
    Nphi = NTOR*FT + 1; Ntheta = NPOL*FP + 1
    phi = np.linspace(0, 2*np.pi, Nphi); theta = np.linspace(0, 2*np.pi, Ntheta)
    R, Z = boundary(theta, phi)
    X = R*np.cos(phi)[None, :]; Y = R*np.sin(phi)[None, :]

    ax.plot_surface(X, Y, Z, color=GREY, alpha=0.12, linewidth=0, antialiased=True, shade=False)

    def fine_quads(pt, pp):
        q = []
        for i in range(pp*FP, pp*FP+FP):
            for j in range(pt*FT, pt*FT+FT):
                q.append([(X[i, j], Y[i, j], Z[i, j]), (X[i, j+1], Y[i, j+1], Z[i, j+1]),
                          (X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]), (X[i+1, j], Y[i+1, j], Z[i+1, j])])
        return q

    if patch_values is not None:
        import matplotlib.cm as cm, matplotlib.colors as mcolors
        vals = np.asarray(patch_values, float)
        finite = vals[np.isfinite(vals)]
        if norm is None:
            norm = mcolors.Normalize(finite.min(), finite.max())
        m = cm.get_cmap(cmap)
        polys, cols = [], []
        for pp in range(NPOL):
            for pt in range(NTOR):
                if np.isfinite(vals[pp, pt]):
                    qs = fine_quads(pt, pp); polys += qs; cols += [m(norm(vals[pp, pt]))]*len(qs)
        ax.add_collection3d(Poly3DCollection(polys, facecolor=cols, alpha=0.75, edgecolor="none"))

    segs = []
    for jt in range(0, Nphi, FT):
        segs += [[(X[i, jt], Y[i, jt], Z[i, jt]), (X[i+1, jt], Y[i+1, jt], Z[i+1, jt])] for i in range(Ntheta-1)]
    for it in range(0, Ntheta, FP):
        segs += [[(X[it, j], Y[it, j], Z[it, j]), (X[it, j+1], Y[it, j+1], Z[it, j+1])] for j in range(Nphi-1)]
    ax.add_collection3d(Line3DCollection(segs, colors=EDGE, linewidths=0.7))

    _, coils = read_makegrid(os.path.join(HERE, "geom/coils_qh"))
    for arr, _ in coils:
        ax.plot(arr[:, 0], arr[:, 1], arr[:, 2], color=KEY, lw=1.0, alpha=0.8)

    xr = (X.min(), X.max()); yr = (Y.min(), Y.max()); zr = (Z.min(), Z.max())
    ax.set_box_aspect((xr[1]-xr[0], yr[1]-yr[0], (zr[1]-zr[0])*1.7))
    ax.set_xlim(xr); ax.set_ylim(yr); ax.set_zlim(zr[0]*1.7, zr[1]*1.7)
    ax.set_axis_off()

if __name__ == "__main__":
    plt.switch_backend("Agg")
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_axes([0, 0, 1, 1], projection="3d")   # fill the whole figure
    draw(ax); ax.view_init(elev=34, azim=-58)
    fig.savefig("figs/patch_geometry.png", dpi=150, bbox_inches="tight", pad_inches=0)
    fig.savefig("figs/patch_geometry.pdf", bbox_inches="tight", pad_inches=0)
    print("saved figs/patch_geometry.png/.pdf")
