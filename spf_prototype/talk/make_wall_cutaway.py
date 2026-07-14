#!/usr/bin/env python
"""Presentation-quality two-panel figure of the multi-layer conformal
stellarator blanket wall (for a beamer talk).

LEFT  : flat poloidal (R,Z) cross-section at a single toroidal angle phi --
        nested closed polygons, one filled band per material layer, legend.
RIGHT : 3D 3/4-torus wedge cutaway exposing the same colored radial build.

Reuses the pure-numpy geometry API in stellarator_geometry.py (load_surface,
poloidal_outward_normals, offset_surface, DEFAULT_LAYERS). Does NOT modify any
existing file. Regenerable: `python make_wall_cutaway.py`.

Outputs: talk/figs/wall_cutaway.pdf and .png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# House style (SM / AASTeX look). smplotlib registers its rcParams on import.
try:
    import smplotlib  # noqa: F401,E402
    _SMPLOTLIB = True
except Exception as _e:  # pragma: no cover - reported at runtime
    _SMPLOTLIB = False
    print(f"WARNING: smplotlib import failed ({_e!r}); "
          "falling back to plain matplotlib.")

from matplotlib.patches import Polygon as MplPolygon  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

CM_PER_M = 100.0  # divide cm coordinates by this to get meters

# import sibling geometry module (run-dir-independent)
HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent / "python"
sys.path.insert(0, str(PYDIR))
import stellarator_geometry as sg  # noqa: E402

OUTDIR = HERE / "figs"
OUTDIR.mkdir(parents=True, exist_ok=True)

STEM = "equil_precise_qa"
SCALE = 10.0  # reactor-scale (matches plot_stellarator_3d.py)

# Distinct, professional categorical colors (talk style, NOT the grayscale paper
# style). Each material is its own category.
LAYER_COLORS = {
    "plasma": "#f2f2f5",  # near-white interior
    "sol":    "#d9dce1",  # light gray vacuum gap
    "W":      "#3a4a63",  # tungsten: dark slate blue-gray
    "steel":  "#7f9fc4",  # RAFM steel: muted blue
    "Be":     "#6bbf73",  # beryllium: green
    "FLiBe":  "#38b8bf",  # breeder: teal/cyan
    "shield": "#e8973a",  # shield: orange
    "coil":   "#b5622c",  # coil: copper/brown
}
LAYER_LABEL = {
    "plasma": "Plasma / SOL interior",
    "sol":    "SOL (vacuum gap)",
    "W":      "W first-wall armor",
    "steel":  "RAFM steel first wall",
    "Be":     "Be multiplier",
    "FLiBe":  "FLiBe breeder",
    "shield": "WC/steel shield",
    "coil":   "Coil block",
}

LABELSIZE = 14


def build_boundaries():
    """Return (phi_axis, boundaries) where boundaries is a list of
    (name, R(theta), Z(theta), thickness) from plasma outward. R,Z are per-node
    arrays of shape (ntheta, nphi)."""
    R, Z, nfp = sg.load_surface(STEM)
    R = R * SCALE
    Z = Z * SCALE
    nR, nZ = sg.poloidal_outward_normals(R, Z)
    nt, nph = R.shape
    phi_axis = np.linspace(0.0, 2 * np.pi, nph, endpoint=False)

    boundaries = [("plasma", R, Z, 0.0)]  # innermost boundary = plasma edge
    cum = 0.0
    for name, thick in sg.DEFAULT_LAYERS:
        cum += thick
        Ro, Zo = sg.offset_surface(R, Z, nR, nZ, cum)
        boundaries.append((name, Ro, Zo, thick))
    return phi_axis, boundaries, nfp


def close(a):
    return np.append(a, a[:1])


# --------------------------------------------------------------------------- #
# LEFT panel: flat poloidal cross-section (onion paint, outermost first)
# --------------------------------------------------------------------------- #
def draw_flat(ax, phi_axis, boundaries, jcol):
    # Paint outermost boundary first, then each inner boundary on top, so the
    # visible remainder of each fill is exactly that layer's annular band.
    # Coordinates converted cm -> m.
    for k in range(len(boundaries) - 1, -1, -1):
        name, R, Z, _ = boundaries[k]
        rr = close(R[:, jcol]) / CM_PER_M
        zz = close(Z[:, jcol]) / CM_PER_M
        # fill color: boundary[k] is the OUTER edge of layer `name` (for k>=1);
        # boundary[0] is the plasma interior.
        color = LAYER_COLORS[name]
        # innermost boundary must sit ON TOP so each inner fill leaves only its
        # own annular band visible -> zorder increases inward.
        poly = MplPolygon(np.column_stack([rr, zz]), closed=True,
                          facecolor=color, edgecolor="0.25", linewidth=0.7,
                          zorder=(len(boundaries) - k))
        ax.add_patch(poly)

    ax.set_aspect("equal")
    ax.autoscale_view()
    allR = np.concatenate([b[1][:, jcol] for b in boundaries]) / CM_PER_M
    allZ = np.concatenate([b[2][:, jcol] for b in boundaries]) / CM_PER_M
    mR = 0.05 * (allR.max() - allR.min())
    mZ = 0.05 * (allZ.max() - allZ.min())
    ax.set_xlim(allR.min() - mR, allR.max() + mR)
    ax.set_ylim(allZ.min() - mZ, allZ.max() + mZ)
    ax.set_xlabel("R [m]", fontsize=LABELSIZE)
    ax.set_ylabel("Z [m]", fontsize=LABELSIZE)
    ax.tick_params(labelsize=LABELSIZE - 2)

    phi_deg = np.degrees(phi_axis[jcol])
    ax.set_title(rf"Poloidal Section ($\phi = {phi_deg:.0f}^\circ$)",
                 fontsize=LABELSIZE)


# --------------------------------------------------------------------------- #
# RIGHT panel: 3D 3/4-torus wedge cutaway (adapted from fig_conformal_build_3d)
# --------------------------------------------------------------------------- #
def torus_xyz(R, Z, phi_axis, wrap_pol=True):
    nt, nph = R.shape
    Rg, Zg = R.copy(), Z.copy()
    Pg = np.broadcast_to(phi_axis, (nt, nph)).copy()
    if wrap_pol:
        Rg = np.concatenate([Rg, Rg[:1, :]], axis=0)
        Zg = np.concatenate([Zg, Zg[:1, :]], axis=0)
        Pg = np.concatenate([Pg, Pg[:1, :]], axis=0)
    # convert cm -> m
    return (Rg * np.cos(Pg) / CM_PER_M,
            Rg * np.sin(Pg) / CM_PER_M,
            Zg / CM_PER_M)


def draw_wedge(ax, phi_axis, boundaries):
    nt, nph = boundaries[0][1].shape
    jcut = int(round(0.75 * nph))  # wedge phi in [0, 1.5 pi]
    sl = slice(0, jcut + 1)

    # faint outer (coil) surface over the wedge for 3D context
    _, Rco, Zco, _ = boundaries[-1]
    Xo, Yo, Zo3 = torus_xyz(Rco[:, sl], Zco[:, sl], phi_axis[sl])
    ax.plot_surface(Xo, Yo, Zo3, color=LAYER_COLORS["coil"], alpha=0.12,
                    rstride=2, cstride=2, linewidth=0, shade=False)

    def ring_poly(Rin, Zin, Rout, Zout, ph):
        # coordinates converted cm -> m
        outer = np.column_stack([Rout * np.cos(ph), Rout * np.sin(ph),
                                 Zout]) / CM_PER_M
        inner = np.column_stack([Rin * np.cos(ph), Rin * np.sin(ph),
                                 Zin]) / CM_PER_M
        return np.vstack([outer, inner[::-1]])

    # nested filled annuli at the two cut faces
    for jend in (0, jcut):
        ph = phi_axis[jend]
        # innermost filled plasma bean (cm -> m)
        _, Rp, Zp, _ = boundaries[0]
        cap = np.column_stack([Rp[:, jend] * np.cos(ph),
                               Rp[:, jend] * np.sin(ph),
                               Zp[:, jend]]) / CM_PER_M
        ax.add_collection3d(Poly3DCollection(
            [cap], facecolor=LAYER_COLORS["plasma"], edgecolor="0.3",
            alpha=0.95, linewidths=0.4))
        for k in range(1, len(boundaries)):
            nm = boundaries[k][0]
            Ri, Zi = boundaries[k - 1][1][:, jend], boundaries[k - 1][2][:, jend]
            Ro_, Zo_ = boundaries[k][1][:, jend], boundaries[k][2][:, jend]
            poly = ring_poly(Ri, Zi, Ro_, Zo_, ph)
            ax.add_collection3d(Poly3DCollection(
                [poly], facecolor=LAYER_COLORS[nm], edgecolor="0.3",
                alpha=0.95, linewidths=0.3))

    ax.set_xlabel("X [m]", fontsize=LABELSIZE, labelpad=12)
    ax.set_ylabel("Y [m]", fontsize=LABELSIZE, labelpad=12)
    ax.set_zlabel("Z [m]", fontsize=LABELSIZE, labelpad=12)
    ax.tick_params(labelsize=LABELSIZE - 4)
    ax.view_init(elev=28, azim=-55)
    dx, dy, dz = np.ptp(Xo), np.ptp(Yo), np.ptp(Zo3)
    ax.set_box_aspect((dx, dy, max(dz, 0.6 * max(dx, dy))))


# --------------------------------------------------------------------------- #
def main():
    phi_axis, boundaries, nfp = build_boundaries()
    print(f"loaded {STEM} (NFP={nfp}), SCALE={SCALE}, "
          f"{len(boundaries)-1} layers, grid {boundaries[0][1].shape}")

    fig = plt.figure(figsize=(16, 7.5))
    axL = fig.add_subplot(1, 2, 1)
    axR = fig.add_subplot(1, 2, 2, projection="3d")

    jcol = 0  # toroidal column for the flat section (phi = 0)
    draw_flat(axL, phi_axis, boundaries, jcol)
    draw_wedge(axR, phi_axis, boundaries)

    # No shared legend (materials listed in the slide text). Leave generous
    # margins so the 3D panel's Z-axis label is not clipped on save (a known
    # mpl 3D issue: bbox_inches='tight' can crop the outermost axis label).
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.08, top=0.92,
                        wspace=0.15)
    pdf = OUTDIR / "wall_cutaway.pdf"
    png = OUTDIR / "wall_cutaway.png"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.4)
    fig.savefig(png, dpi=200, bbox_inches="tight", pad_inches=0.4)
    plt.close(fig)
    for p in (pdf, png):
        print(f"  wrote {p}  ({p.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
