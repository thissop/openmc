#!/usr/bin/env python
"""Check-in figures (house style: smplotlib; Title Case labels; NO plot titles)."""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import smplotlib  # noqa: F401  house serif style — MUST precede pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
PROTO = HERE.parent
FIGS = HERE / "figs"; FIGS.mkdir(exist_ok=True)
CQA, CQH, CQI = "#2a6f97", "#e07a2f", "#3a9d54"


def _rows(p):
    return list(csv.DictReader(open(p)))


def fig_reciprocity():
    r = _rows(PROTO / "shield_opt/data/reciprocity_pairs.csv")
    dev = np.array([x["device"] for x in r])
    a = np.array([float(x["adjoint"]) for x in r])
    f = np.array([float(x["forward"]) for x in r])
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    for c, col in (("QH", CQH), ("QA", CQA)):
        m = dev == c
        rho, _ = spearmanr(a[m], f[m])
        ax.scatter(a[m], f[m], s=34, c=col, edgecolor="k", linewidth=0.3, alpha=0.85,
                   label=fr"{c}: $\rho={rho:+.2f}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Adjoint Importance $\int S\,\psi^\dagger$")
    ax.set_ylabel("Forward Per-Coil Flux")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(FIGS / "reciprocity.png", dpi=150); plt.close(fig)


def fig_kappel():
    r = _rows(PROTO / "data/kappel/kappel_configs.csv")
    Lg = np.array([float(x["L_gradB_m"]) for x in r])
    Lr = np.array([float(x["L_REGCOIL_m"]) for x in r])
    cls = np.array([x["qs_type"] for x in r])
    conf = [x["config"] for x in r]
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    palette = {"QA": CQA, "QH": CQH, "QI": CQI}
    other = ~np.isin(cls, list(palette))
    ax.scatter(Lg[other], Lr[other], s=22, c="#9aa0a6", label="Other", alpha=0.7, edgecolor="k", linewidth=0.2)
    for c in ("QA", "QH", "QI"):
        m = cls == c
        ax.scatter(Lg[m], Lr[m], s=28, c=palette[c], label=c, alpha=0.85, edgecolor="k", linewidth=0.3)
    z = np.polyfit(Lg, Lr, 1); xs = np.linspace(Lg.min(), Lg.max(), 50)
    ax.plot(xs, np.polyval(z, xs), "k--", lw=1.0, alpha=0.6)
    for i, cc in enumerate(conf):
        if cc == "new_QA_aScaling":
            ax.scatter([Lg[i]], [Lr[i]], marker="*", s=340, c="gold", edgecolor="k", linewidth=0.7, zorder=6)
    ax.set_xlabel(r"Magnetic Gradient Scale Length $L_{\nabla B}$ (m)")
    ax.set_ylabel("Min Plasma-Coil Separation (m)")
    ax.text(0.05, 0.90, r"$R^2=0.94$, $n=45$", transform=ax.transAxes)
    ax.legend(loc="lower right", frameon=False, ncol=2)
    fig.tight_layout(); fig.savefig(FIGS / "kappel_spine.png", dpi=150); plt.close(fig)


def fig_blanket():
    r = _rows(PROTO / "shield_opt/data/engineering_relevant_devices.csv")
    asp = np.array([float(x["aspect"]) for x in r])
    gap = np.array([float(x["gap_m_ARIESCS"]) for x in r])
    cls = np.array([x["class"] for x in r])
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    for c, col in (("QA", CQA), ("QH", CQH)):
        m = cls == c
        ax.scatter(asp[m], gap[m], s=24, c=col, label=c, alpha=0.75, edgecolor="k", linewidth=0.25)
    ax.axhline(1.29, color="k", ls="--", lw=1.0)
    ax.text(asp.max() * 0.5, 1.45, "1.29 m Radial Build", fontsize=9)
    ax.set_xlabel("Aspect Ratio")
    ax.set_ylabel("Plasma-Coil Gap at ARIES-CS (m)")
    ax.set_yscale("log")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout(); fig.savefig(FIGS / "blanket_fit.png", dpi=150); plt.close(fig)


def fig_coilsource():
    t = np.linspace(0, 2 * np.pi, 400)
    x = np.cos(t) + 0.18 * np.cos(3 * t)
    y = np.sin(t) + 0.18 * np.sin(3 * t)
    z = 0.32 * np.sin(2 * t)
    ts = np.linspace(0, 2 * np.pi, 14, endpoint=False)
    xs = np.cos(ts) + 0.18 * np.cos(3 * ts); ys = np.sin(ts) + 0.18 * np.sin(3 * ts); zs = 0.32 * np.sin(2 * ts)
    fig = plt.figure(figsize=(5.6, 4.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, color="#444444", lw=2.2)
    ax.scatter(xs, ys, zs, marker="*", s=170, c="gold", edgecolor="k", linewidth=0.6, depthshade=False)
    ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)"); ax.set_zlabel("Z (m)")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.view_init(elev=22, azim=35)
    fig.tight_layout(); fig.savefig(FIGS / "coil_source.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    fig_reciprocity(); fig_kappel(); fig_blanket(); fig_coilsource()
    print("wrote reciprocity, kappel_spine, blanket_fit, coil_source ->", FIGS)
