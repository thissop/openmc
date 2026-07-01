#!/usr/bin/env python
"""Figure for the bare-cylinder demonstration (run bare_cylinder_demo.py first).

Left: in-band directionality A/iso from the analytic point-patch flux integral vs the
exact ray-trace, for the inner (convex-to-source) cylinder, against shaping amplitude
s -- coincident at s=0, diverging as s grows. Right: the analytic<->ray-trace
discrepancy vs s for the inner and outer cylinder, with an s^2 guide; zero (to MC
statistics) without shaping, growing with shaping and largest on the convex inner wall.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import smplotlib  # noqa: E402,F401

REPO = Path(__file__).resolve().parents[2]
FIGS = REPO / "spf_prototype" / "figs"
plt.rcParams.update({"axes.titlesize": 15, "axes.labelsize": 13,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11})

d = np.load("/tmp/bare_cylinder_demo.npz")["rows"]
s = d[:, 0]; ai, ri, di = d[:, 1], d[:, 2], d[:, 3]; ao, ro, do = d[:, 4], d[:, 5], d[:, 6]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.6))

axL.plot(s, ai, "k-o", mfc="white", mec="k", lw=1.8, label="Analytic (point-patch flux)")
axL.plot(s, ri, "k--s", mfc="0.6", mec="k", lw=1.8, label="Ray-trace (exact transport)")
axL.set_xlabel("Shaping amplitude $s$")
axL.set_ylabel("Inner-wall directionality $A/\\mathrm{iso}$")
axL.set_title("Inner cylinder (convex toward source)")
axL.legend(loc="upper left", frameon=True)

axR.plot(s, di, "k-o", mfc="white", mec="k", lw=1.8, label="Inner (convex, `inboard`)")
axR.plot(s, do, "k--^", mfc="0.6", mec="k", lw=1.8, label="Outer (concave, `outboard`)")
sg = np.linspace(0, 1, 50)
axR.plot(sg, di[-1] * sg ** 2, color="0.6", lw=1.0, ls=":", label="$\\propto s^2$ guide")
axR.axhline(0.003, color="0.75", lw=0.9, ls=(0, (1, 2)))
axR.text(0.03, 0.006, "MC statistical floor", fontsize=9, color="0.4")
axR.set_xlabel("Shaping amplitude $s$")
axR.set_ylabel("$|A/\\mathrm{iso}_{\\rm analytic} - A/\\mathrm{iso}_{\\rm raytrace}|$")
axR.set_title("Point-patch error vs shaping")
axR.legend(loc="upper left", frameon=True)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIGS / f"bare_cylinder_demo.{ext}", dpi=150, bbox_inches="tight")
print("wrote bare_cylinder_demo.pdf/.png")
