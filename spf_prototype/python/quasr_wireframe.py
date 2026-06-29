#!/usr/bin/env python
"""Bare 3-D wireframe renders of QUASR plasma boundaries (no axes, labels, or grid)
for figures/. Geometry from quasr_geom (VMEC namelist RBC/ZBS, pure numpy).

Run: $HOME/spf_venv/bin/python spf_prototype/python/quasr_wireframe.py [ID ...]
Default: the gentle shortlist devices used in the validation.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402

REPO = Path(__file__).resolve().parents[2]
FIGS = REPO / "spf_prototype" / "figs"
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
import quasr_geom as G  # noqa: E402

DEFAULT_IDS = [59509, 802624, 932746, 399015]
NPHI, NTHETA = 160, 56


def wireframe(ID, iota=0.1):
    dev = G.load_device(ID, iota)
    ph = np.linspace(0, 2 * np.pi, NPHI)
    th = np.linspace(0, 2 * np.pi, NTHETA)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z = dev.RZ(TH, PH, 1.0)
    X = R * np.cos(PH); Y = R * np.sin(PH)

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=3, color="k", linewidth=0.4)
    # bare: no panes, axes, grid, ticks, labels
    ax.set_axis_off()
    ax.grid(False)
    ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))   # true proportions
    ax.view_init(elev=52, azim=35)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"quasr{ID}_wireframe.{ext}", dpi=200,
                    bbox_inches="tight", pad_inches=0, transparent=True)
    plt.close(fig)
    print(f"wrote quasr{ID}_wireframe.png/.pdf  (nfp {dev.nfp}, R0 {dev.R0:.3f})")


if __name__ == "__main__":
    ids = [int(a) for a in sys.argv[1:]] or DEFAULT_IDS
    for ID in ids:
        wireframe(ID)
