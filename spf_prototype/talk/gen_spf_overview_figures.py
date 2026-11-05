"""Generate the added SPF-overview and StellaratorSource-validation figures.

The figures are sized for full-width placement in the 16:9 Beamer deck and use
the repository plotting protocol's monochrome technical-illustration language.
"""
import smplotlib  # noqa: F401 -- plotting protocol: first import statement

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIGS = HERE / "figs" / "summer_recap_2026"
FIGS.mkdir(parents=True, exist_ok=True)


plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "axes.grid": False,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.9,
    "axes.labelsize": 12,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "lines.linewidth": 1.35,
    "mathtext.fontset": "stix",
    "text.color": "black",
})


def finish(fig, stem, pad=0.02):
    fig.savefig(FIGS / f"{stem}.pdf", bbox_inches="tight", pad_inches=pad)
    fig.savefig(FIGS / f"{stem}.png", bbox_inches="tight", pad_inches=pad)
    plt.close(fig)


def _spin_arrow(ax, x, projection, label):
    """Draw a compact glyph for a nuclear-spin projection on the field axis."""
    if projection == 0:
        ax.plot([x - 0.10, x + 0.10], [0, 0], color="black", linewidth=1.35)
        ax.scatter([x], [0], s=25, facecolor="white", edgecolor="black",
                   linewidth=1.0, zorder=3)
    else:
        y0 = -0.08 * projection
        y1 = 0.78 * projection
        ax.annotate("", xy=(x, y1), xytext=(x, y0),
                    arrowprops=dict(arrowstyle="-|>", color="black",
                                    linewidth=1.35, mutation_scale=11))
    ax.text(x, -1.04, label, ha="center", va="center", fontsize=11)


def spin_substates_bw():
    """Examples of the A, B, and C population groups about local B-hat."""
    cases = [
        ("A-State Example", r"$d_{+}t_{+}$", +1, +1,
         r"$m_D=+1,\ m_T=+\frac{1}{2}$"),
        ("B-State Example", r"$d_{0}t_{+}$", 0, +1,
         r"$m_D=0,\ m_T=+\frac{1}{2}$"),
        ("C-State Example", r"$d_{-}t_{+}$", -1, +1,
         r"$m_D=-1,\ m_T=+\frac{1}{2}$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.35))
    for ax, (title, state, d_proj, t_proj, quantum_label) in zip(axes, cases, strict=True):
        ax.annotate("", xy=(0, 1.13), xytext=(0, -0.88),
                    arrowprops=dict(arrowstyle="-|>", color="0.55",
                                    linewidth=1.15, mutation_scale=10))
        ax.text(0.04, 1.11, r"$\widehat{\mathbf{B}}$", ha="left", va="center",
                fontsize=12)
        _spin_arrow(ax, -0.34, d_proj, r"D")
        _spin_arrow(ax, +0.34, t_proj, r"T")
        ax.text(0, 0.98, state, ha="center", va="top", fontsize=13)
        ax.text(0, -1.27, quantum_label, ha="center", va="center", fontsize=10)
        ax.set_title(title, fontsize=14, pad=16)
        ax.set_xlim(-0.82, 0.82)
        ax.set_ylim(-1.44, 1.22)
        ax.set_aspect("equal")
        ax.set_axis_off()
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.08, top=0.82, wspace=0.06)
    finish(fig, "spin_substates_bw", pad=0.0)


def quasr59509_cross_sections_bw():
    """Overlay the QUASR 59509 LCFS at three toroidal phases."""
    fluxmap = dict(np.load(DATA / "quasr59509_vmec_fluxmap.npz"))
    radius = np.asarray(fluxmap["R"])[-1]
    height = np.asarray(fluxmap["Z"])[-1]
    zeta = np.asarray(fluxmap["zeta"])
    nfp = int(fluxmap["nfp"])

    phases = np.array([0.0, np.pi / (2 * nfp), np.pi / nfp])
    indices = [int(np.argmin(np.abs(zeta - phase))) for phase in phases]
    labels = [
        r"$\zeta=0$",
        r"$\zeta=\frac{\pi}{2N}$",
        r"$\zeta=\frac{\pi}{N}$",
    ]
    styles = [
        dict(color="black", linestyle="-", linewidth=1.55),
        dict(color="0.30", linestyle="--", linewidth=1.45),
        dict(color="0.55", linestyle=":", linewidth=1.65),
    ]

    fig, ax = plt.subplots(figsize=(7.6, 4.25))
    plotted_r = []
    plotted_z = []
    for index, label, style in zip(indices, labels, styles, strict=True):
        r_curve = np.append(radius[:, index], radius[0, index])
        z_curve = np.append(height[:, index], height[0, index])
        ax.plot(r_curve, z_curve, label=label, **style)
        plotted_r.append(r_curve)
        plotted_z.append(z_curve)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Major Radius $R$ (m)", fontsize=15)
    ax.set_ylabel(r"Vertical Position $Z$ (m)", fontsize=15)
    ax.tick_params(axis="both", which="major", labelsize=13.5)

    # Reserve explicit white space at right for the legend so no legend sample
    # or label crosses the rotating LCFS curves.
    r_all = np.concatenate(plotted_r)
    z_all = np.concatenate(plotted_z)
    r_span = np.ptp(r_all)
    z_span = np.ptp(z_all)
    ax.set_xlim(np.min(r_all) - 0.10 * r_span, np.max(r_all) + 0.68 * r_span)
    ax.set_ylim(np.min(z_all) - 0.13 * z_span, np.max(z_all) + 0.13 * z_span)
    ax.legend(loc="center right", bbox_to_anchor=(0.985, 0.52),
              fontsize=14, handlelength=2.2, labelspacing=0.62)
    fig.tight_layout(pad=0.55)
    finish(fig, "quasr59509_cross_sections_bw")


if __name__ == "__main__":
    spin_substates_bw()
    quasr59509_cross_sections_bw()
    print(f"Wrote full-width figures to {FIGS}")
