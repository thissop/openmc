"""Fig: the filamentary point-patch analytic error grows ~s^2 with shaping and is
~3x worse on the convex inboard wall -- where Monte Carlo is required."""
import csv
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

rows = list(csv.DictReader(r for r in open(DATA / "bare_cylinder.csv") if not r.startswith("#")))
s = np.array([float(r["s"]) for r in rows])
inner = np.array([float(r["inner_absdiff"]) for r in rows])
outer = np.array([float(r["outer_absdiff"]) for r in rows])

fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.plot(s, inner, "o-", color="k", label="Inboard (Convex)")
ax.plot(s, outer, "s--", color="k", label="Outboard (Concave)")
ss = np.linspace(0, 1, 50)
ax.plot(ss, inner[-1] * ss**2, ":", color="k", lw=1.0, label="~ s^2")
ax.set_xlabel("Plasma Shaping Amplitude $s$")
ax.set_ylabel("Analytic Point-Patch Error |A/iso - ray|")
ax.legend(frameon=False, loc="upper left", fontsize=8)
ax.annotate("3x Worse\nInboard", (1.0, inner[-1]), xytext=(0.62, 0.075),
            fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
finish(fig, "fig_bare_cylinder")
