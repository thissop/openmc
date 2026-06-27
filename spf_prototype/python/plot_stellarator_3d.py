#!/usr/bin/env python
"""Publication-quality 3D + supporting visualizations of the SPF stellarator work.

PURE matplotlib visualization from the committed DESC / field-map / conformal-build
data under spf_prototype/data/. No OpenMC, no DAGMC, no transport. Only the
pure-numpy helpers (fieldmap.py, stellarator_geometry.py) are imported.

Run from the repo root:
    /home/tkiker.guest/spf_venv/bin/python spf_prototype/python/plot_stellarator_3d.py

Figures written to spf_prototype/figs/:
  1. stellarator_lcfs_3d.png            full-torus LCFS + 3/4 cutaway, colored by B_z
  2. stellarator_cross_sections.png     poloidal cross-sections vs toroidal angle
  3. stellarator_conformal_build_3d.png nested conformal radial build (cutaway)
  4. stellarator_field_3d.png           3D quiver of B-hat on the LCFS (iota)
  5. stellarator_field_pitch.png        field pitch vs poloidal angle (3D variation)
  6. stellarator_spf_lobes_3d.png       SPF emission lobes about the local B-hat
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import smplotlib  # noqa: E402,F401  house serif style
import matplotlib.cm as cm  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "spf_prototype" / "data"
FIGS = REPO / "spf_prototype" / "figs"
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))

from fieldmap import FieldMapField  # noqa: E402
import stellarator_geometry as sg  # noqa: E402

STEM = "equil_precise_qa"
DPI = 130

# smplotlib renders fonts small; force readable sizes everywhere.
plt.rcParams.update({
    "axes.titlesize": 18,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 20,
})

TS, LS = 17, 14  # explicit title / label fontsizes for 3D axes


# --------------------------------------------------------------------------- #
# Data loading (cached at module import via load_all())
# --------------------------------------------------------------------------- #
def load_all():
    """Load the LCFS surface grid and sample B-hat on every surface node."""
    d = np.load(DATA / f"{STEM}_surface.npz")
    R = np.asarray(d["R"]); Z = np.asarray(d["Z"])
    phi = np.asarray(d["phi"]); nfp = int(d["nfp"])
    theta = np.asarray(d["theta"])
    nt, nph = R.shape
    phi_axis = phi[0, :].copy()  # phi is constant along theta

    fm = FieldMapField(str(DATA / STEM))
    X = R * np.cos(phi); Y = R * np.sin(phi)
    B = np.empty((nt, nph, 3))
    for i in range(nt):
        for j in range(nph):
            B[i, j] = fm.bhat((X[i, j], Y[i, j], Z[i, j]))
    return dict(R=R, Z=Z, phi=phi, phi_axis=phi_axis, theta=theta,
                nfp=nfp, nt=nt, nph=nph, X=X, Y=Y, B=B, fm=fm)


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def torus_xyz(R, Z, phi, wrap_tor=True, wrap_pol=True):
    """Build closed Cartesian grids X,Y,Z for plot_surface. Optionally wrap the
    seam in the toroidal (phi) and/or poloidal (theta) directions so the surface
    has no visible gap."""
    Rg, Zg, Pg = R.copy(), Z.copy(), phi.copy()
    if wrap_tor:
        Rg = np.concatenate([Rg, Rg[:, :1]], axis=1)
        Zg = np.concatenate([Zg, Zg[:, :1]], axis=1)
        Pg = np.concatenate([Pg, Pg[:, :1] + 2 * np.pi], axis=1)
    if wrap_pol:
        Rg = np.concatenate([Rg, Rg[:1, :]], axis=0)
        Zg = np.concatenate([Zg, Zg[:1, :]], axis=0)
        Pg = np.concatenate([Pg, Pg[:1, :]], axis=0)
    Xc = Rg * np.cos(Pg); Yc = Rg * np.sin(Pg)
    return Xc, Yc, Zg


def quad_centers(A):
    """Vertex field -> quad-centered field (for plot_surface facecolors)."""
    return 0.25 * (A[:-1, :-1] + A[1:, :-1] + A[:-1, 1:] + A[1:, 1:])


def equalize_3d(ax, X, Y, Z, zfrac=0.0):
    """Box aspect from the data extents. This is a high-aspect (~6) torus, so a
    faithful aspect renders as a pancake; `zfrac` sets a floor on the z box height
    (as a fraction of the xy footprint) to gently exaggerate the tube for
    readability without changing the X/Y proportions."""
    dx = np.ptp(X); dy = np.ptp(Y); dz = np.ptp(Z)
    dz_eff = max(dz, zfrac * max(dx, dy))
    ax.set_box_aspect((dx, dy, dz_eff))


def local_frame(bhat):
    """Gram-Schmidt local frame with z_local = bhat (least-aligned ref axis)."""
    bhat = np.asarray(bhat, float)
    bhat /= np.linalg.norm(bhat)
    e = np.eye(3)[np.argmin(np.abs(bhat))]
    x_l = e - np.dot(e, bhat) * bhat
    x_l /= np.linalg.norm(x_l)
    y_l = np.cross(bhat, x_l)
    return x_l, y_l, bhat


# --------------------------------------------------------------------------- #
# Figure 1 -- LCFS full torus + 3/4 cutaway, colored by B_z
# --------------------------------------------------------------------------- #
def fig_lcfs_3d(D):
    R, Z, phi, B = D["R"], D["Z"], D["phi"], D["B"]
    cmap = cm.coolwarm
    norm = Normalize(vmin=B[:, :, 2].min(), vmax=B[:, :, 2].max())

    fig = plt.figure(figsize=(15, 7))

    # --- left: full torus ---
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    Xf, Yf, Zf = torus_xyz(R, Z, phi, wrap_tor=True, wrap_pol=True)
    Bzf = np.concatenate([B[:, :, 2], B[:, :1, 2]], axis=1)
    Bzf = np.concatenate([Bzf, Bzf[:1, :]], axis=0)
    fc = cmap(norm(quad_centers(Bzf)))
    ax.plot_surface(Xf, Yf, Zf, facecolors=fc, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False)
    ax.set_title("precise_QA LCFS (NFP=2 Quasi-Axisymmetric)", fontsize=TS, pad=12)
    ax.set_xlabel("X [cm]", fontsize=LS); ax.set_ylabel("Y [cm]", fontsize=LS)
    ax.set_zlabel("Z [cm]", fontsize=LS)
    ax.view_init(elev=30, azim=-60)
    equalize_3d(ax, Xf, Yf, Zf, zfrac=0.6)

    # --- right: 3/4 cutaway (phi in [0, 1.5 pi]) with filled bean caps ---
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    jcut = int(round(0.75 * D["nph"]))  # column nearest 1.5*pi
    sl = slice(0, jcut + 1)
    Xc, Yc, Zc = torus_xyz(R[:, sl], Z[:, sl], phi[:, sl],
                           wrap_tor=False, wrap_pol=True)
    Bzc = np.concatenate([B[:, sl, 2], B[:1, sl, 2]], axis=0)
    fc2 = cmap(norm(quad_centers(Bzc)))
    ax2.plot_surface(Xc, Yc, Zc, facecolors=fc2, rstride=1, cstride=1,
                     linewidth=0, antialiased=True, shade=False)
    # filled cross-section caps at the two cut faces
    for jend, ec in ((0, "0.25"), (jcut, "0.25")):
        ph = phi[0, jend]
        cap = np.column_stack([R[:, jend] * np.cos(ph),
                               R[:, jend] * np.sin(ph), Z[:, jend]])
        pc = Poly3DCollection([cap], facecolor="#9ec6e0", edgecolor=ec,
                              alpha=0.85, linewidths=0.6)
        ax2.add_collection3d(pc)
    ax2.set_title("LCFS 3/4 Cutaway (Rotating Bean Cross-Section)",
                  fontsize=TS, pad=12)
    ax2.set_xlabel("X [cm]", fontsize=LS); ax2.set_ylabel("Y [cm]", fontsize=LS)
    ax2.set_zlabel("Z [cm]", fontsize=LS)
    ax2.view_init(elev=30, azim=-60)
    equalize_3d(ax2, Xf, Yf, Zf, zfrac=0.6)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=[ax, ax2], shrink=0.6, pad=0.02)
    cb.set_label(r"$\hat{B}_z$ on LCFS (poloidal tilt)", fontsize=LS)
    out = FIGS / "stellarator_lcfs_3d.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Figure 2 -- poloidal cross-sections vs toroidal angle
# --------------------------------------------------------------------------- #
def fig_cross_sections(D):
    R, Z, phi_axis = D["R"], D["Z"], D["phi_axis"]
    targets = [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]
    labels = [r"$\phi = 0$", r"$\phi = \pi/4$", r"$\phi = \pi/2$",
              r"$\phi = 3\pi/4$"]
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

    fig, ax = plt.subplots(figsize=(8, 8))
    for tgt, lab, col in zip(targets, labels, colors):
        j = int(np.argmin(np.abs(phi_axis - tgt)))
        rr = np.append(R[:, j], R[0, j]); zz = np.append(Z[:, j], Z[0, j])
        ax.plot(rr, zz, "-", color=col, lw=2.0,
                label=f"{lab}  (col {j}, {np.degrees(phi_axis[j]):.0f}$^\\circ$)")
    ax.set_aspect("equal")
    ax.set_xlabel("R [cm]", fontsize=15); ax.set_ylabel("Z [cm]", fontsize=15)
    ax.set_title("LCFS Cross-Sections vs Toroidal Angle", fontsize=18, pad=10)
    ax.legend(fontsize=12, loc="upper right")
    ax.grid(True, alpha=0.25)
    out = FIGS / "stellarator_cross_sections.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Figure 3 -- nested conformal radial build (cutaway)
# --------------------------------------------------------------------------- #
LAYER_COLORS = {
    "plasma": "#e9a3c9",
    "sol": "#cccccc",
    "W": "#4d4d4d",
    "steel": "#5a8fbf",
    "Be": "#7bc47f",
    "FLiBe": "#f0c040",
    "shield": "#b07aa1",
    "coil": "#c0612a",
}


def fig_conformal_build_3d(D):
    SCALE = 10.0
    Rs, Zs, nfp = sg.load_surface(STEM)
    Rs = Rs * SCALE; Zs = Zs * SCALE
    nRn, nZn = sg.poloidal_outward_normals(Rs, Zs)
    nt, nph = Rs.shape
    phi_axis = np.linspace(0.0, 2 * np.pi, nph, endpoint=False)
    phi2d = np.broadcast_to(phi_axis, (nt, nph))

    cum = 0.0
    boundaries = [("plasma", Rs, Zs)]  # innermost = plasma boundary
    for name, t in sg.DEFAULT_LAYERS:
        cum += t
        Ro, Zo = sg.offset_surface(Rs, Zs, nRn, nZn, cum)
        boundaries.append((name, Ro, Zo))

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    jcut = int(round(0.75 * nph))  # wedge phi in [0, 1.5 pi]
    sl = slice(0, jcut + 1)

    # faint outer (coil) surface over the wedge for 3D context
    _, Rco, Zco = boundaries[-1]
    Xo, Yo, Zo3 = torus_xyz(Rco[:, sl], Zco[:, sl], phi2d[:, sl],
                            wrap_tor=False, wrap_pol=True)
    ax.plot_surface(Xo, Yo, Zo3, color=LAYER_COLORS["coil"], alpha=0.12,
                    rstride=2, cstride=2, linewidth=0, shade=False)

    # nested filled annuli at the two cut faces -> shows the radial build
    def ring_poly(Rin, Zin, Rout, Zout, ph):
        outer = np.column_stack([Rout * np.cos(ph), Rout * np.sin(ph), Zout])
        inner = np.column_stack([Rin * np.cos(ph), Rin * np.sin(ph), Zin])
        return np.vstack([outer, inner[::-1]])

    layer_names = ["plasma"] + [n for n, _ in sg.DEFAULT_LAYERS]
    for jend in (0, jcut):
        ph = phi_axis[jend]
        # innermost: filled plasma bean
        _, Rp, Zp = boundaries[0]
        cap = np.column_stack([Rp[:, jend] * np.cos(ph),
                               Rp[:, jend] * np.sin(ph), Zp[:, jend]])
        ax.add_collection3d(Poly3DCollection(
            [cap], facecolor=LAYER_COLORS["plasma"], edgecolor="0.3",
            alpha=0.9, linewidths=0.4))
        # each layer annulus
        for k in range(1, len(boundaries)):
            nm = boundaries[k][0]
            Ri, Zi = boundaries[k - 1][1][:, jend], boundaries[k - 1][2][:, jend]
            Ro_, Zo_ = boundaries[k][1][:, jend], boundaries[k][2][:, jend]
            poly = ring_poly(Ri, Zi, Ro_, Zo_, ph)
            ax.add_collection3d(Poly3DCollection(
                [poly], facecolor=LAYER_COLORS[nm], edgecolor="0.3",
                alpha=0.92, linewidths=0.3))

    handles = [Patch(facecolor=LAYER_COLORS[n], edgecolor="0.3",
                     label=n) for n in layer_names]
    ax.legend(handles=handles, fontsize=11, loc="upper left",
              title="Radial build", title_fontsize=12)
    ax.set_title("Conformal Radial Build on a QA Stellarator (Machinery Test)",
                 fontsize=TS, pad=12)
    ax.set_xlabel("X [cm]", fontsize=LS); ax.set_ylabel("Y [cm]", fontsize=LS)
    ax.set_zlabel("Z [cm]", fontsize=LS)
    ax.view_init(elev=28, azim=-55)
    equalize_3d(ax, Xo, Yo, Zo3, zfrac=0.6)
    fig.text(0.5, 0.03,
             "Layer thicknesses are reactor-scaled placeholders "
             "(SCALE=10, plasma->W->steel->Be->FLiBe->shield->coil); "
             "geometry-machinery test, not an engineered build.",
             ha="center", fontsize=10, style="italic", color="0.3")
    out = FIGS / "stellarator_conformal_build_3d.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Figure 4 -- 3D quiver of B-hat on the LCFS
# --------------------------------------------------------------------------- #
def fig_field_3d(D):
    R, Z, phi, B = D["R"], D["Z"], D["phi"], D["B"]
    X, Y = D["X"], D["Y"]
    cmap = cm.coolwarm
    norm = Normalize(vmin=B[:, :, 2].min(), vmax=B[:, :, 2].max())

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")

    # faint LCFS underneath
    Xf, Yf, Zf = torus_xyz(R, Z, phi, wrap_tor=True, wrap_pol=True)
    ax.plot_surface(Xf, Yf, Zf, color="0.7", alpha=0.10, rstride=2, cstride=2,
                    linewidth=0, shade=False)

    ti = slice(None, None, 6)   # every 6th theta
    pj = slice(None, None, 4)   # every 4th phi
    xs = X[ti, pj].ravel(); ys = Y[ti, pj].ravel(); zs = Z[ti, pj].ravel()
    us = B[ti, pj, 0].ravel(); vs = B[ti, pj, 1].ravel(); ws = B[ti, pj, 2].ravel()
    bz = B[ti, pj, 2].ravel()

    c = cmap(norm(bz))
    colors = np.concatenate((c, np.repeat(c, 2, axis=0)))  # shafts, then heads
    ax.quiver(xs, ys, zs, us, vs, ws, colors=colors, length=18.0,
              normalize=True, linewidth=1.6, arrow_length_ratio=0.4)

    ax.set_title(r"Magnetic Field Direction $\hat{B}(x)$ on the LCFS "
                 r"(Rotational Transform)", fontsize=TS, pad=12)
    ax.set_xlabel("X [cm]", fontsize=LS); ax.set_ylabel("Y [cm]", fontsize=LS)
    ax.set_zlabel("Z [cm]", fontsize=LS)
    ax.view_init(elev=32, azim=-60)
    equalize_3d(ax, Xf, Yf, Zf, zfrac=0.6)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label(r"$\hat{B}_z$ (poloidal pitch component)", fontsize=LS)
    out = FIGS / "stellarator_field_3d.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Figure 5 -- field pitch vs poloidal angle
# --------------------------------------------------------------------------- #
def fig_field_pitch(D):
    B, phi, theta, phi_axis = D["B"], D["phi"], D["theta"], D["phi_axis"]
    targets = [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]
    labels = [r"$\phi = 0$", r"$\phi = \pi/4$", r"$\phi = \pi/2$",
              r"$\phi = 3\pi/4$"]
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    theta_deg = np.degrees(theta)

    fig, ax = plt.subplots(figsize=(10, 6.5))
    for tgt, lab, col in zip(targets, labels, colors):
        j = int(np.argmin(np.abs(phi_axis - tgt)))
        ph = phi[0, j]
        phihat = np.array([-np.sin(ph), np.cos(ph), 0.0])
        dotp = np.abs(B[:, j, :] @ phihat)
        pitch = np.degrees(np.arccos(np.clip(dotp, 0.0, 1.0)))
        ax.plot(theta_deg, pitch, "-", color=col, lw=2.0, label=lab)
    ax.set_xlabel(r"Poloidal Angle $\theta$ [deg]", fontsize=15)
    ax.set_ylabel(r"Field Pitch $\arccos|\hat{B}\cdot\hat{\phi}|$ [deg]",
                  fontsize=15)
    ax.set_title("Field Pitch vs Poloidal Angle (3D Field-Period Variation)",
                 fontsize=18, pad=10)
    ax.set_xlim(0, 360); ax.set_xticks(np.arange(0, 361, 60))
    ax.legend(fontsize=12, loc="best")
    ax.grid(True, alpha=0.25)
    out = FIGS / "stellarator_field_pitch.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Figure 6 -- SPF emission lobes about the local B-hat (optional)
# --------------------------------------------------------------------------- #
def _lobe_surface(P, bhat, radial_fn, scale, n_th=60, n_az=60):
    """Surface of revolution r(theta_l) about bhat, anchored at point P (cm)."""
    x_l, y_l, z_l = local_frame(bhat)
    th = np.linspace(0, np.pi, n_th)
    az = np.linspace(0, 2 * np.pi, n_az)
    TH, AZ = np.meshgrid(th, az, indexing="ij")
    r = radial_fn(TH)
    dir_ = (np.sin(TH) * np.cos(AZ))[..., None] * x_l \
        + (np.sin(TH) * np.sin(AZ))[..., None] * y_l \
        + (np.cos(TH))[..., None] * z_l
    pts = P + scale * r[..., None] * dir_
    return pts[..., 0], pts[..., 1], pts[..., 2], r


def fig_spf_lobes_3d(D):
    R, Z, phi, B = D["R"], D["Z"], D["phi"], D["B"]
    X, Y = D["X"], D["Y"]
    # outboard-midplane-ish anchor: theta=0 row, phi=pi/4 column
    j = int(np.argmin(np.abs(D["phi_axis"] - np.pi / 4)))
    i = 0
    P = np.array([X[i, j], Y[i, j], Z[i, j]])
    bhat = B[i, j].copy()

    fig = plt.figure(figsize=(13, 7))

    a_fn = lambda t: np.sin(t) ** 2                       # A mode  (perp lobe)
    bc_fn = lambda t: 0.25 + 0.75 * np.cos(t) ** 2        # B/C mode (parallel)
    scale = 45.0

    for k, (fn, name, cmap_, ttl) in enumerate([
        (a_fn, "A", cm.Reds, r"A mode  $\propto \sin^2\theta$  ($\perp\,\hat{B}$)"),
        (bc_fn, "BC", cm.Blues,
         r"B/C mode  $\propto \frac{1}{4}+\frac{3}{4}\cos^2\theta$  "
         r"($\parallel\,\hat{B}$)")]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        # faint local LCFS patch for context
        ii = slice(max(0, i - 8), i + 9) if i - 8 >= 0 else slice(0, 17)
        jj = slice(max(0, j - 8), j + 9)
        ax.plot_surface(X[:, jj], Y[:, jj], Z[:, jj], color="0.7", alpha=0.10,
                        rstride=2, cstride=2, linewidth=0, shade=False)
        Lx, Ly, Lz, rr = _lobe_surface(P, bhat, fn, scale)
        norm = Normalize(rr.min(), rr.max())
        fc = cmap_(0.35 + 0.6 * norm(quad_centers(rr)))
        ax.plot_surface(Lx, Ly, Lz, facecolors=fc, rstride=1, cstride=1,
                        linewidth=0, antialiased=True, shade=False, alpha=0.95)
        # local B-hat arrow
        ax.quiver(*P, *bhat, length=70.0, color="k", linewidth=2.2,
                  arrow_length_ratio=0.25)
        ax.text(*(P + 75 * bhat), r"$\hat{B}$", fontsize=14)
        ax.set_title(ttl, fontsize=15, pad=10)
        ax.set_xlabel("X [cm]", fontsize=12); ax.set_ylabel("Y [cm]", fontsize=12)
        ax.set_zlabel("Z [cm]", fontsize=12)
        ax.view_init(elev=22, azim=-70)
        allx = np.concatenate([Lx.ravel(), [P[0]]])
        ally = np.concatenate([Ly.ravel(), [P[1]]])
        allz = np.concatenate([Lz.ravel(), [P[2]]])
        ax.set_box_aspect((np.ptp(allx), np.ptp(ally), np.ptp(allz)))

    fig.suptitle("SPF Neutron Emission Lobes about the Local Field "
                 r"$\hat{B}$ on the LCFS", fontsize=18, y=0.98)
    out = FIGS / "stellarator_spf_lobes_3d.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight"); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    print("Loading data + sampling B-hat on LCFS ...")
    D = load_all()
    print(f"  LCFS grid {D['nt']}x{D['nph']} (theta x phi), NFP={D['nfp']}, "
          f"B_z in [{D['B'][:,:,2].min():.3f}, {D['B'][:,:,2].max():.3f}]")

    makers = [
        ("Figure 1 LCFS 3D", fig_lcfs_3d),
        ("Figure 2 cross-sections", fig_cross_sections),
        ("Figure 3 conformal build", fig_conformal_build_3d),
        ("Figure 4 field quiver", fig_field_3d),
        ("Figure 5 field pitch", fig_field_pitch),
        ("Figure 6 SPF lobes", fig_spf_lobes_3d),
    ]
    for name, fn in makers:
        out = fn(D)
        kb = out.stat().st_size / 1024
        flag = "OK " if kb > 20 else "!! "
        print(f"  {flag}{name:28s} -> {out.name}  ({kb:.0f} KB)")


if __name__ == "__main__":
    main()
