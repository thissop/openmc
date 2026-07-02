#!/usr/bin/env python3
"""
figures.py -- regenerate every figure for the SPF conformal-stellarator deck.

Reads ONLY ./spf_export (plain data: npz / csv / json -- no OpenMC, no DAGMC,
no cluster) and writes 200-dpi PNGs into ./figs, plus ./figs/derived_numbers.json
(the exact headline numbers, so the slides / talk / Q&A quote computed values).

Run:  python3 figures.py

Data model notes (established by inspecting the bundle):
  * results/tallies_long.csv is machine-written and NOT RFC-clean: several
    fields contain unescaped commas -- the mesh-tally filter tuples
    ('1','1','1','x-min out'), the reaction name '(n,Xt)', and the coil energy
    filter '6;(100000.0, 20000000.0)'. A naive csv reader mis-columns those
    rows. We parse tally-by-tally against the known layout instead. We only
    need the well-behaved-once-realigned cell_response (108 rows) and coil_fast
    (6 rows); the 1.9M noisy wall_current_phi mesh rows are not used here.
  * summary.json in this bundle is empty ({}), so all headline metrics are
    recomputed from tallies_long.csv and cross-checked against
    results/RESULTS_tier8_conformal.md (they agree to the printed precision).
  * Geometry is exported at the BASE equilibrium scale (LCFS R ~ 63-133 cm),
    not multiplied by the run's device scale of 15. We plot it as-is, in cm.
"""

from __future__ import annotations
import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------- #
# Paths and global style
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent
SRC = ROOT / "spf_export"
GEO = SRC / "geometry"
RES = SRC / "results"
OUT = ROOT / "figs"
OUT.mkdir(exist_ok=True)

DPI = 200

# Mode colors (validated categorical palette; blue<->orange is the strongest
# CVD-safe pair, gray is a designated neutral baseline -- every bar is also
# direct-labeled so identity is never carried by color alone).
MODE_COLOR = {
    "unpolarized": "#8a8884",   # neutral gray  -- reference / baseline
    "perpendicular": "#2a78d6",  # blue (cool)   -- A mode, sin^2 th, reduces coil flux
    "parallel": "#eb6834",       # orange (warm) -- B/C mode, 1+3cos^2 th, raises coil flux
}
MODE_ORDER = ["unpolarized", "perpendicular", "parallel"]
MODE_LABEL = {
    "unpolarized": "unpolarized",
    "perpendicular": "perpendicular (A)",
    "parallel": "parallel (B/C)",
}
STREAM_LABEL = {"free": "free-streaming", "scatter": "with scattering"}

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.edgecolor": INK2,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "axes.grid": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
})


def _save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {p.relative_to(ROOT)}")


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_layers():
    with open(GEO / "layers.json") as f:
        return json.load(f)


def load_boundaries():
    return np.load(GEO / "boundaries_RZ.npz")


def load_shells():
    return np.load(GEO / "shells_xyz.npz")


def parse_results():
    """Robustly parse the well-formed physics rows out of tallies_long.csv.

    Returns two lists of dicts:
      cell -- keys: stream, mode, score, mat, mean, std, relerr
              score in {heating, damage-energy, (n,Xt)};  mat in
              {W, steel, Be, FLiBe, shield, coil}
      coil -- keys: stream, mode, mean, std, relerr   (coil fast flux >0.1 MeV)
    We dispatch on the tally name and reassemble fields the comma-split broke.
    """
    cell, coil = [], []
    with open(RES / "tallies_long.csv") as f:
        header = f.readline()  # discard
        for line in f:
            f_ = line.rstrip("\n").split(",")
            tally = f_[2]
            stream, mode = f_[0], f_[1]
            if tally == "cell_response":
                if f_[4] == "(n":
                    # score '(n,Xt)' was split into f_[4]='(n', f_[5]='Xt)'
                    score = "(n,Xt)"
                    mat = f_[8]
                    mean, std, relerr = f_[9], f_[10], f_[11]
                else:
                    # 'heating' or 'damage-energy' -- row is aligned
                    score = f_[4]
                    mat = f_[7]
                    mean, std, relerr = f_[8], f_[9], f_[10]
                cell.append(dict(
                    stream=stream, mode=mode, score=score, mat=mat,
                    mean=float(mean), std=float(std),
                    relerr=float("nan") if relerr == "nan" else float(relerr),
                ))
            elif tally == "coil_fast":
                # filter '6;(100000.0, 20000000.0)' split f_[6]/f_[7];
                # material col f_[8] empty; mean/std/relerr at 9/10/11.
                coil.append(dict(
                    stream=stream, mode=mode,
                    mean=float(f_[9]), std=float(f_[10]), relerr=float(f_[11]),
                ))
            # wall_current_phi: skip (noisy mesh, not used for these figures)
    return cell, coil


def get_cell(cell, stream, score, mat):
    """Return (mean, std) for one cell_response entry, keyed by mode."""
    out = {}
    for r in cell:
        if r["stream"] == stream and r["score"] == score and r["mat"] == mat:
            out[r["mode"]] = (r["mean"], r["std"])
    return out


def get_coil(coil, stream):
    return {r["mode"]: (r["mean"], r["std"]) for r in coil if r["stream"] == stream}


# --------------------------------------------------------------------------- #
# Figure 1 -- 3D LCFS wireframe
# --------------------------------------------------------------------------- #
def fig1_wireframe():
    import pandas as pd
    df = pd.read_csv(GEO / "lcfs_wireframe.csv")
    nt = df.theta_idx.max() + 1
    nph = df.phi_idx.max() + 1
    X = np.full((nt, nph), np.nan)
    Y = np.full((nt, nph), np.nan)
    Z = np.full((nt, nph), np.nan)
    for _, r in df.iterrows():
        X[int(r.theta_idx), int(r.phi_idx)] = r.x_cm
        Y[int(r.theta_idx), int(r.phi_idx)] = r.y_cm
        Z[int(r.theta_idx), int(r.phi_idx)] = r.z_cm
    # close the torus in phi so the wireframe wraps
    X = np.concatenate([X, X[:, :1]], axis=1)
    Y = np.concatenate([Y, Y[:, :1]], axis=1)
    Z = np.concatenate([Z, Z[:, :1]], axis=1)
    nfp = int(load_boundaries()["nfp"])

    fig = plt.figure(figsize=(8, 6.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(X, Y, Z, rstride=3, cstride=2,
                      color="#2a78d6", linewidth=0.5, alpha=0.85)
    ax.set_xlabel("x [cm]", labelpad=10)
    ax.set_ylabel("y [cm]", labelpad=10)
    ax.set_zlabel("z [cm]", labelpad=6)
    # equal aspect from data extents
    xr = np.nanmax(X) - np.nanmin(X)
    yr = np.nanmax(Y) - np.nanmin(Y)
    zr = np.nanmax(Z) - np.nanmin(Z)
    ax.set_box_aspect((xr, yr, zr))
    ax.view_init(elev=32, azim=-58)
    ax.set_title(f"Plasma boundary (LCFS), precise_QA, {nfp} field periods",
                 pad=6, fontsize=13)
    fig.tight_layout()
    _save(fig, "fig1_lcfs_wireframe.png")


# --------------------------------------------------------------------------- #
# Figure 2 -- poloidal cross sections at two toroidal angles (non-axisymmetry)
# --------------------------------------------------------------------------- #
# gap (inner boundary, outer boundary) -> material filling that shell
GAPS = [
    ("LCFS", "plasma_edge", "sol"),
    ("plasma_edge", "W_outer", "W"),
    ("W_outer", "steel_outer", "steel"),
    ("steel_outer", "Be_outer", "Be"),
    ("Be_outer", "FLiBe_outer", "FLiBe"),
    ("FLiBe_outer", "shield_outer", "shield"),
    ("shield_outer", "coil_outer", "coil"),
]


def _ring(df_phi, boundary):
    s = df_phi[df_phi.layer == boundary].sort_values("theta_idx")
    R = s.R_cm.to_numpy()
    Z = s.Z_cm.to_numpy()
    return np.column_stack([R, Z])


def _annulus_patch(outer, inner, color, alpha=1.0, lw=0.0):
    """Compound-path annulus (outer ring minus inner ring), even-odd fill.
    Robust for non-convex / axis-crossing offset surfaces."""
    outer = np.vstack([outer, outer[:1]])
    inner = np.vstack([inner, inner[:1]])
    verts = np.vstack([outer, inner[::-1]])
    codes = ([MplPath.MOVETO] + [MplPath.LINETO] * (len(outer) - 2) + [MplPath.CLOSEPOLY]
             + [MplPath.MOVETO] + [MplPath.LINETO] * (len(inner) - 2) + [MplPath.CLOSEPOLY])
    return PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="white",
                     lw=lw, alpha=alpha, joinstyle="round")


def _draw_cross_section(ax, df_phi, layers, colors):
    # plasma interior (inside LCFS): faint fill
    lcfs = _ring(df_phi, "LCFS")
    ax.fill(lcfs[:, 0], lcfs[:, 1], color="#f3ecff", zorder=0)
    # shells outward
    for inner_b, outer_b, mat in GAPS:
        inner = _ring(df_phi, inner_b)
        outer = _ring(df_phi, outer_b)
        col = colors[mat]
        ax.add_patch(_annulus_patch(outer, inner, col, alpha=1.0))
    # boundary outlines (thin) for legibility
    for b in ["LCFS", "FLiBe_outer", "coil_outer"]:
        r = _ring(df_phi, b)
        r = np.vstack([r, r[:1]])
        ax.plot(r[:, 0], r[:, 1], color="white", lw=0.8, zorder=5)
    ax.set_aspect("equal")
    ax.set_xlabel("R [cm]")


def fig2_cross_section():
    import pandas as pd
    df = pd.read_csv(GEO / "cross_sections.csv")
    layers = {l["name"]: l for l in load_layers()["layers"]}
    colors = {name: layers[name]["color"] for name in layers}
    colors["sol"] = layers["sol"]["color"]

    phis = [0.0, 90.0]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
    # shared limits from the outermost boundary across both phis
    allR, allZ = [], []
    for ph in phis:
        d = df[df.phi_deg == ph]
        co = _ring(d, "coil_outer")
        allR.append(co[:, 0]); allZ.append(co[:, 1])
    allR = np.concatenate(allR); allZ = np.concatenate(allZ)
    pad = 12
    xlim = (allR.min() - pad, allR.max() + pad)
    ylim = (allZ.min() - pad, allZ.max() + pad)

    for ax, ph in zip(axes, phis):
        d = df[df.phi_deg == ph]
        _draw_cross_section(ax, d, layers, colors)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_title(rf"$\phi = {int(ph)}^\circ$", fontsize=12)
    axes[0].set_ylabel("Z [cm]")

    # legend from GAPS materials + plasma
    handles = [Patch(facecolor="#f3ecff", edgecolor=MUTED, label="plasma (void)")]
    seen = set()
    for _, _, mat in GAPS:
        if mat in seen:
            continue
        seen.add(mat)
        lab = f"{layers[mat]['material']}" + (" (vacuum)" if layers[mat].get("is_vacuum") else "")
        handles.append(Patch(facecolor=colors[mat], edgecolor="white", label=lab))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.995, 0.5),
               frameon=False, fontsize=10, title="layer")
    fig.suptitle("Conformal blanket cross-section, precise_QA (nfp = 2)",
                 fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0, 0.9, 0.97))
    _save(fig, "fig2_cross_section.png")


# --------------------------------------------------------------------------- #
# Figure 3 -- 3D nested translucent shells (conformal wrap)
# --------------------------------------------------------------------------- #
def fig3_nested_shells():
    d = load_shells()
    layers = {l["name"]: l for l in load_layers()["layers"]}
    # boundary -> (color, alpha, material label)
    shows = [
        ("coil_outer", layers["coil"]["color"], 0.16, "coil case (outer)"),
        ("FLiBe_outer", layers["FLiBe"]["color"], 0.28, "FLiBe breeder (outer)"),
        ("LCFS", "#7a4fd0", 0.9, "plasma (LCFS)"),
    ]
    # toroidal cutaway: keep a 300-deg wedge so the nesting is visible at the cut
    nph = d["LCFS_x"].shape[1]
    keep = slice(0, int(nph * 0.80))

    fig = plt.figure(figsize=(8.2, 6.6))
    ax = fig.add_subplot(111, projection="3d")
    for b, col, al, _lab in shows:
        X = d[f"{b}_x"][:, keep]
        Y = d[f"{b}_y"][:, keep]
        Z = d[f"{b}_z"][:, keep]
        ax.plot_surface(X, Y, Z, color=col, alpha=al, linewidth=0,
                        antialiased=True, shade=True, rstride=1, cstride=1)
    ax.set_xlabel("x [cm]", labelpad=10)
    ax.set_ylabel("y [cm]", labelpad=10)
    ax.set_zlabel("z [cm]", labelpad=6)
    Xall = d["coil_outer_x"]; Yall = d["coil_outer_y"]; Zall = d["coil_outer_z"]
    ax.set_box_aspect((np.ptp(Xall), np.ptp(Yall), np.ptp(Zall)))
    ax.view_init(elev=30, azim=-62)
    ax.set_title("Conformal blanket wraps the plasma (300° cutaway)",
                 pad=6, fontsize=13)
    handles = [Patch(facecolor=c, alpha=min(a + 0.25, 1.0), label=l)
               for _, c, a, l in shows]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, "fig3_nested_shells.png")


# --------------------------------------------------------------------------- #
# Figure 4 -- radial layer schematic (to scale)
# --------------------------------------------------------------------------- #
def _draw_stack(ax, layers, y0, h):
    """Draw each layer as a to-scale rectangle at its true radial offset."""
    for L in layers:
        ax.add_patch(plt.Rectangle((L["inner_offset_cm"], y0), L["thickness_cm"], h,
                                   facecolor=L["color"], edgecolor="white",
                                   lw=1.0, zorder=2))


def fig4_radial_schematic():
    layers = load_layers()["layers"]
    total = sum(L["thickness_cm"] for L in layers)
    fig, ax = plt.subplots(figsize=(11, 3.7))
    y0, h = 0.0, 1.0
    _draw_stack(ax, layers, y0, h)

    # inline labels for the thick outer layers (legible at full scale)
    for L in layers:
        if L["thickness_cm"] >= 8:
            cx = L["inner_offset_cm"] + L["thickness_cm"] / 2
            ax.text(cx, y0 + h / 2, f"{L['material']}\n{L['thickness_cm']:g} cm",
                    ha="center", va="center", fontsize=11, color="white",
                    fontweight="bold", zorder=3)

    ax.set_xlim(-3, total + 3)
    ax.set_ylim(-0.15, 2.35)
    ax.set_yticks([])
    ax.set_xlabel("radial offset from plasma surface [cm]")
    ax.spines[["left", "top", "right"]].set_visible(False)

    # inset zoom for the thin plasma-side layers (sol / W / steel / Be, 0-11 cm)
    inner = [L for L in layers if L["inner_offset_cm"] < 11.0]
    axz = ax.inset_axes([0.04, 1.02, 0.44, 0.82], transform=ax.transAxes)
    _draw_stack(axz, inner, 0.0, 1.0)
    heights = {"sol": 1.30, "W": 2.00, "steel": 1.30, "Be": 2.00}
    for L in inner:
        cx = L["inner_offset_cm"] + L["thickness_cm"] / 2
        name = L["material"] + (" (vac)" if L.get("is_vacuum") else "")
        axz.annotate(f"{name}  {L['thickness_cm']:g} cm", xy=(cx, 1.0),
                     xytext=(cx, heights.get(L["material"], 1.6)),
                     ha="center", va="bottom", fontsize=8.5, color=INK,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    axz.set_xlim(-0.3, 11.0)
    axz.set_ylim(0, 2.9)
    axz.set_yticks([])
    axz.set_xticks([0, 5, 11])
    axz.tick_params(labelsize=8)
    axz.set_title("inner layers, zoomed (0 to 11 cm)", fontsize=9, color=INK2)
    for s in ["left", "top", "right"]:
        axz.spines[s].set_visible(False)
    ax.indicate_inset_zoom(axz, edgecolor=MUTED, alpha=0.6)

    ax.set_title(f"Radial blanket stack (plasma surface to coil), total {total:g} cm",
                 fontsize=13, pad=8, loc="right")
    fig.subplots_adjust(left=0.04, right=0.98, bottom=0.16, top=0.62)
    _save(fig, "fig4_radial_schematic.png")


# --------------------------------------------------------------------------- #
# Helpers for result bar figures
# --------------------------------------------------------------------------- #
def _bar_labels(ax, bars, vals, fmt, dy=0.0, fontsize=8.5, rotation=0):
    for b, v in zip(bars, vals):
        ax.annotate(fmt(v), (b.get_x() + b.get_width() / 2, b.get_height() + dy),
                    ha="center", va="bottom", fontsize=fontsize, color=INK,
                    rotation=rotation)


def _clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------- #
# Figure 5a -- TBR (FLiBe H3 production) vs mode, scattering
# --------------------------------------------------------------------------- #
def fig5a_tbr(cell):
    sc = get_cell(cell, "scatter", "(n,Xt)", "FLiBe")
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    xs = np.arange(3)
    means = [sc[m][0] for m in MODE_ORDER]
    stds = [sc[m][1] for m in MODE_ORDER]
    cols = [MODE_COLOR[m] for m in MODE_ORDER]
    bars = ax.bar(xs, means, yerr=stds, capsize=5, color=cols,
                  edgecolor="white", width=0.62, zorder=3,
                  error_kw=dict(ecolor=INK2, lw=1.1))
    for x, m, s in zip(xs, means, stds):
        ax.annotate(f"{m:.3f}\n±{s:.3f}", (x, m + s), ha="center", va="bottom",
                    fontsize=9.5, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([MODE_LABEL[m] for m in MODE_ORDER])
    ax.set_ylabel("TBR  (tritium atoms / source neutron)")
    ax.set_ylim(0, max(means) * 1.28)
    _clean(ax)
    ax.set_title("Tritium breeding ratio in FLiBe, with scattering", fontsize=12.5)
    ax.text(0.5, -0.22, "Free-streaming TBR ≈ 0 by construction "
            "(no scattering, so no moderation and no capture). 10k-history smoke run.",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color=INK2)
    fig.tight_layout()
    _save(fig, "fig5a_tbr.png")


# --------------------------------------------------------------------------- #
# Figure 5b -- coil fast flux vs mode, free vs scatter
# --------------------------------------------------------------------------- #
def fig5b_coil_flux(coil):
    free = get_coil(coil, "free")
    sca = get_coil(coil, "scatter")
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    groups = ["free", "scatter"]
    data = {"free": free, "scatter": sca}
    gx = np.arange(2) * 1.6
    w = 0.42
    offs = np.array([-w, 0, w])
    for j, m in enumerate(MODE_ORDER):
        means = [data[g][m][0] for g in groups]
        stds = [data[g][m][1] for g in groups]
        xs = gx + offs[j]
        bars = ax.bar(xs, means, yerr=stds, width=w, color=MODE_COLOR[m],
                      edgecolor="white", capsize=4, zorder=3,
                      label=MODE_LABEL[m], error_kw=dict(ecolor=INK2, lw=1.0))
        for x, v in zip(xs, means):
            ax.annotate(f"{v:.1f}", (x, v), ha="center", va="bottom",
                        fontsize=8.5, color=INK, xytext=(0, 2),
                        textcoords="offset points")
    ax.set_xticks(gx)
    ax.set_xticklabels([STREAM_LABEL[g] for g in groups])
    ax.set_ylabel("coil fast flux  >0.1 MeV  [/cm$^2$/src]")
    ax.set_ylim(0, 42)
    _clean(ax)
    ax.legend(frameon=False, fontsize=9.5, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.14))
    ax.set_title("Coil fast flux -- magnet-lifetime proxy", fontsize=12.5, pad=26)
    ax.text(0.5, -0.20, "Parallel raises, perpendicular lowers the coil load; "
            "the ordering survives scattering. 10k-history smoke run "
            "(± shown; large MC error at depth).",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color=INK2)
    fig.tight_layout()
    _save(fig, "fig5b_coil_flux.png")


# --------------------------------------------------------------------------- #
# Figure 5c -- per-material heating (scattering) by layer and mode
# --------------------------------------------------------------------------- #
def fig5c_heating(cell):
    mats = ["W", "steel", "Be", "FLiBe", "shield", "coil"]
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    gx = np.arange(len(mats)) * 1.5
    w = 0.40
    offs = np.array([-w, 0, w])
    for j, m in enumerate(MODE_ORDER):
        means, stds = [], []
        for mat in mats:
            d = get_cell(cell, "scatter", "heating", mat)
            means.append(d[m][0]); stds.append(d[m][1])
        xs = gx + offs[j]
        ax.bar(xs, means, yerr=stds, width=w, color=MODE_COLOR[m],
               edgecolor="white", capsize=3, zorder=3, label=MODE_LABEL[m],
               error_kw=dict(ecolor=INK2, lw=0.9))
    ax.set_yscale("log")
    ax.set_xticks(gx)
    ax.set_xticklabels(mats)
    ax.set_ylabel("heating  [eV / source neutron]")
    ax.set_ylim(3e3, 2e7)
    ax.grid(axis="y", which="both", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.09))
    fig.suptitle("Per-layer nuclear heating, with scattering", fontsize=12.5, y=1.01)
    ax.text(0.5, -0.16, "FLiBe (50 cm breeder) absorbs most of the energy. At the "
            "coil, parallel > unpolarized > perpendicular -- same directional "
            "ordering as the coil flux. 10k-history smoke run.",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color=INK2)
    fig.tight_layout()
    _save(fig, "fig5c_heating.png")


# --------------------------------------------------------------------------- #
# Figure 5d -- directional efficiency eta
# --------------------------------------------------------------------------- #
def _frac_change(p, u):
    """(p-u)/u with 1-sigma error, p=(mean,std), u=(mean,std)."""
    pm, ps = p
    um, us = u
    r = pm / um
    sig = r * math.sqrt((ps / pm) ** 2 + (us / um) ** 2)
    return (r - 1.0), sig  # fractional change, its sigma


def fig5d_eta(coil, numbers):
    free = get_coil(coil, "free")
    sca = get_coil(coil, "scatter")
    # parallel vs unpolarized fractional change, free and scatter
    df_free, sf_free = _frac_change(free["parallel"], free["unpolarized"])
    df_sca, sf_sca = _frac_change(sca["parallel"], sca["unpolarized"])
    eta = df_sca / df_free
    eta_err = abs(eta) * math.sqrt((sf_free / df_free) ** 2 + (sf_sca / df_sca) ** 2)
    # perpendicular (context)
    pf_free, _ = _frac_change(free["perpendicular"], free["unpolarized"])
    pf_sca, _ = _frac_change(sca["perpendicular"], sca["unpolarized"])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.2, 4.7),
                                   gridspec_kw=dict(width_ratios=[2.1, 1]))
    # Left: fractional change of the two polarized modes vs unpolarized
    xs = np.array([0, 1])
    w = 0.36
    par = [df_free * 100, df_sca * 100]
    per = [pf_free * 100, pf_sca * 100]
    par_err = [sf_free * 100, sf_sca * 100]
    b1 = axL.bar(xs - w / 2, par, w, yerr=par_err, capsize=5,
                 color=MODE_COLOR["parallel"], edgecolor="white",
                 label="parallel (B/C)", zorder=3, error_kw=dict(ecolor=INK2, lw=1.1))
    b2 = axL.bar(xs + w / 2, per, w, color=MODE_COLOR["perpendicular"],
                 edgecolor="white", label="perpendicular (A)", zorder=3)
    axL.axhline(0, color=INK2, lw=1.0)
    for x, v in zip(xs - w / 2, par):
        axL.annotate(f"{v:+.1f}%", (x, v), ha="center",
                     va="bottom" if v >= 0 else "top",
                     xytext=(0, 4 if v >= 0 else -4), textcoords="offset points",
                     fontsize=9.5, color=INK)
    for x, v in zip(xs + w / 2, per):
        axL.annotate(f"{v:+.1f}%", (x, v), ha="center",
                     va="top" if v < 0 else "bottom",
                     xytext=(0, -4 if v < 0 else 4), textcoords="offset points",
                     fontsize=9.5, color=INK)
    axL.set_xticks(xs)
    axL.set_xticklabels(["free-streaming", "with scattering"])
    axL.set_ylabel("coil fast-flux change vs unpolarized [%]")
    axL.set_ylim(-24, 22)
    axL.spines[["top", "right"]].set_visible(False)
    axL.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    axL.set_axisbelow(True)
    axL.legend(frameon=False, fontsize=9.5, loc="lower left")
    axL.set_title("Directional signature vs unpolarized", fontsize=12)

    # Right: hero eta
    axR.axis("off")
    axR.text(0.5, 0.74, "$\\eta$", ha="center", fontsize=30, color=INK)
    axR.text(0.5, 0.50, f"{eta:.2f}", ha="center", fontsize=46, color=MODE_COLOR["parallel"],
             fontweight="bold")
    axR.text(0.5, 0.34, f"± {eta_err:.2f}", ha="center", fontsize=14, color=INK2)
    axR.text(0.5, 0.17, "fraction of the free-streaming\nsteering signal that survives\n"
             "conformal walls + scattering", ha="center", va="top",
             fontsize=9.5, color=INK2)
    axR.text(0.5, 0.98, f"$\\eta = \\dfrac{{{df_sca*100:+.1f}\\%}}{{{df_free*100:+.1f}\\%}}$",
             ha="center", va="top", fontsize=13, color=INK)

    fig.suptitle("Does SPF steering survive a real stellarator wall?",
                 fontsize=13.5, y=1.02)
    fig.text(0.5, -0.02, "Large Monte-Carlo error on a 10k-history smoke run "
             "(η = 0.63 ± 0.22); indicative, needs variance reduction "
             "(weight windows / FW-CADIS) for a precise value.",
             ha="center", fontsize=8.5, color=INK2)
    fig.tight_layout()
    _save(fig, "fig5d_eta.png")

    numbers["eta"] = dict(value=eta, error=eta_err,
                          frac_change_free=df_free, frac_change_free_err=sf_free,
                          frac_change_scatter=df_sca, frac_change_scatter_err=sf_sca,
                          perp_frac_change_free=pf_free, perp_frac_change_scatter=pf_sca)


# --------------------------------------------------------------------------- #
# Derived numbers dump (grounds the slides / talk / Q&A)
# --------------------------------------------------------------------------- #
def dump_numbers(cell, coil, numbers):
    numbers["tbr_scatter"] = {m: get_cell(cell, "scatter", "(n,Xt)", "FLiBe")[m]
                              for m in MODE_ORDER}
    numbers["coil_fast_free"] = {m: get_coil(coil, "free")[m] for m in MODE_ORDER}
    numbers["coil_fast_scatter"] = {m: get_coil(coil, "scatter")[m] for m in MODE_ORDER}
    numbers["heating_scatter"] = {
        mat: {m: get_cell(cell, "scatter", "heating", mat)[m] for m in MODE_ORDER}
        for mat in ["W", "steel", "Be", "FLiBe", "shield", "coil"]
    }
    # parallel-vs-perpendicular contrast (the "directional signature" magnitude)
    cf = numbers["coil_fast_free"]; cs = numbers["coil_fast_scatter"]
    numbers["par_vs_perp_free_pct"] = (cf["parallel"][0] / cf["perpendicular"][0] - 1) * 100
    numbers["par_vs_perp_scatter_pct"] = (cs["parallel"][0] / cs["perpendicular"][0] - 1) * 100
    with open(OUT / "derived_numbers.json", "w") as f:
        json.dump(numbers, f, indent=2)
    print(f"  wrote {(OUT / 'derived_numbers.json').relative_to(ROOT)}")


# --------------------------------------------------------------------------- #
def main():
    print("Reading spf_export, writing figs/ ...")
    cell, coil = parse_results()
    numbers = {}
    fig1_wireframe()
    fig2_cross_section()
    fig3_nested_shells()
    fig4_radial_schematic()
    fig5a_tbr(cell)
    fig5b_coil_flux(coil)
    fig5c_heating(cell)
    fig5d_eta(coil, numbers)
    dump_numbers(cell, coil, numbers)
    print("Done.")


if __name__ == "__main__":
    main()
