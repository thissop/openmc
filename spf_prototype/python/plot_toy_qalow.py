#!/usr/bin/env python
"""Paper figures for the toy QA-low free-streaming OpenMC<->analytic cross-check.

Reads the two products of the validation (run those first):
  /tmp/toy_qalow_analytic.npz  (anarrima exact quadrature, per wall patch; fine grid)
  /tmp/toy_qalow_openmc.npz    (OpenMC free-streaming, per mesh bin)
and writes, to spf_prototype/figs/:
  toy_qalow_directionality.pdf/.png  A/iso, B/iso poloidal profiles, MC vs analytic
  toy_qalow_agreement.pdf/.png       MC vs analytic scatter (y=x), all walls/modes
  toy_qalow_geometry.pdf/.png        rotating plasma cross-sections in the square wall

House style: smplotlib serif, black-and-white (line style, not colour, separates the
modes), open white-faced markers for the Monte-Carlo points. Fixed seed upstream.

Run: $HOME/spf_venv/bin/python spf_prototype/python/plot_toy_qalow.py
"""
import importlib
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import smplotlib  # noqa: E402,F401  house serif style

REPO = Path(__file__).resolve().parents[2]
FIGS = REPO / "spf_prototype" / "figs"
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
# config-driven: SPF_CONFIG=quasr_config plots the QUASR cross-check instead
C = importlib.import_module(os.environ.get("SPF_CONFIG", "toy_qalow_config"))

plt.rcParams.update({
    "axes.titlesize": 16, "axes.labelsize": 14,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11, "figure.titlesize": 18,
})

WALLS = ("inboard", "outboard", "floor", "ceiling")
WALL_TITLE = {"inboard": "Inboard wall", "outboard": "Outboard wall",
              "floor": "Floor", "ceiling": "Ceiling"}
XLAB = {"inboard": "Height $z$ (norm.)", "outboard": "Height $z$ (norm.)",
        "floor": "Major radius $R$ (norm.)", "ceiling": "Major radius $R$ (norm.)"}
TITLE = getattr(C, "TITLE", C.STEM)
DPI = 150


def load():
    an = np.load(f"/tmp/{C.STEM}_analytic.npz")
    mc = np.load(f"/tmp/{C.STEM}_openmc.npz")
    return an, mc


def mc_ratio(mc, wall, mode, floor_frac=0.10):
    """Per-bin directionality mode/iso for one wall; mask low-statistics corner bins."""
    s = mc[f"s_{wall}"]; iso = mc[f"iso_{wall}"]; num = mc[f"{mode}_{wall}"]
    good = iso > floor_frac * iso.max()
    r = np.where(good, num / np.where(good, iso, 1.0), np.nan)
    return s, r, iso


# --------------------------------------------------------------------------- #
def fig_directionality(an, mc):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, wall in zip(axes.ravel(), WALLS):
        m = an["wall"] == wall
        sa = an["s"][m]; order = np.argsort(sa)
        sa = sa[order]
        Aa = an["A_iso"][m][order]; Ba = an["B_iso"][m][order]
        # analytic (smooth lines)
        ax.plot(sa, Aa, "k-", lw=2.0, label="A  analytic")
        ax.plot(sa, Ba, "k--", lw=2.0, label="B/C  analytic")
        # MC per-bin points + wall mean (the validated table number)
        for mode, mk in (("A", "^"), ("B", "s")):
            s, r, iso = mc_ratio(mc, wall, mode)
            ax.plot(s, r, mk, ms=6, mfc="white", mec="0.35", mew=1.1, ls="none",
                    label=f"{('A' if mode=='A' else 'B/C')}  OpenMC")
            wm = np.nansum(mc[f"{mode}_{wall}"]) / np.nansum(iso)
            ax.axhline(wm, color="0.55", lw=0.8, ls=":")
        ax.axhline(1.0, color="0.7", lw=1.0, ls=(0, (1, 2)))
        ax.set_title(WALL_TITLE[wall])
        ax.set_xlabel(XLAB[wall]); ax.set_ylabel("Directionality (NWL / isotropic)")
        ax.set_ylim(0.6, 1.45)
    axes[0, 0].legend(loc="upper center", ncol=2, frameon=True, handlelength=2.4)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{C.STEM}_directionality.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {C.STEM}_directionality.pdf/.png")


def fig_agreement(an, mc):
    """MC vs analytic (interpolated to MC bin centres), every wall and both modes."""
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    pts = {"A": ([], []), "B": ([], [])}
    for wall in WALLS:
        m = an["wall"] == wall
        sa = an["s"][m]; o = np.argsort(sa); sa = sa[o]
        for mode in ("A", "B"):
            aa = an[f"{mode}_iso"][m][o]
            s, r, _ = mc_ratio(mc, wall, mode)
            an_at = np.interp(s, sa, aa)
            ok = np.isfinite(r)
            pts[mode][0].extend(an_at[ok]); pts[mode][1].extend(r[ok])
    lo, hi = 0.7, 1.35
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=1.2, ls="--", label="y = x")
    ax.plot(pts["A"][0], pts["A"][1], "^", ms=6, mfc="white", mec="k", mew=1.0,
            ls="none", label="A mode")
    ax.plot(pts["B"][0], pts["B"][1], "s", ms=6, mfc="0.7", mec="k", mew=1.0,
            ls="none", label="B/C mode")
    allx = np.array(pts["A"][0] + pts["B"][0]); ally = np.array(pts["A"][1] + pts["B"][1])
    rel = np.abs(ally - allx) / np.abs(allx)
    ax.text(0.04, 0.93, f"per-bin |rel. err|: mean {rel.mean():.1%}, "
            f"median {np.median(rel):.1%}", transform=ax.transAxes, fontsize=11)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("Analytic directionality (anarrima)")
    ax.set_ylabel("OpenMC directionality (free-streaming)")
    ax.set_title(f"Per-bin agreement, {TITLE}")
    ax.legend(loc="lower right", frameon=True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{C.STEM}_agreement.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {C.STEM}_agreement.pdf/.png")


def fig_geometry():
    """Rotating plasma cross-sections (several toroidal cuts) inside the square wall."""
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    th = np.linspace(0, 2 * np.pi, 200)
    cuts = [(0.0, "k-", r"$\phi=0$"),
            (np.pi / (2 * C.NFP), "k--", r"$\phi=\pi/2N_{\rm fp}$"),
            (np.pi / C.NFP, "k:", r"$\phi=\pi/N_{\rm fp}$")]
    for phi, sty, lab in cuts:
        for rho, lw in ((1.0, 1.8), (0.6, 1.0)):
            R = np.array([C.loop_RZ(t, np.array([phi]), rho)[0][0] for t in th])
            Z = np.array([C.loop_RZ(t, np.array([phi]), rho)[1][0] for t in th])
            ax.plot(R, Z, sty, lw=lw, label=(lab if rho == 1.0 else None))
    # square wall
    ax.plot([C.R_IN, C.R_OUT, C.R_OUT, C.R_IN, C.R_IN],
            [-C.Z_W, -C.Z_W, C.Z_W, C.Z_W, -C.Z_W], color="0.45", lw=2.5)
    ax.text(C.R_OUT - 0.02, C.Z_W - 0.03, "wall", ha="right", va="top", color="0.45")
    ax.set_xlabel("Major radius $R$ (norm.)"); ax.set_ylabel("Height $z$ (norm.)")
    ax.set_title(f"{TITLE} source loops (LCFS + $\\rho=0.6$) in the square wall")
    ax.set_aspect("equal"); ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{C.STEM}_geometry.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {C.STEM}_geometry.pdf/.png")


if __name__ == "__main__":
    an, mc = load()
    fig_directionality(an, mc)
    fig_agreement(an, mc)
    fig_geometry()
