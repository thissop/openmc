"""Fig: the rate-vs-steering trade-off. A buys +50% power but compounds the center-stack
load to 1.70x; B is the sweet spot (steering at no rate cost); C halves everything."""
import csv
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

rows = {r["mode"]: r for r in csv.DictReader(
    x for x in open(DATA / "rate_steering.csv") if not x.startswith("#"))}
modes = ["nonpol", "A", "B", "C"]
metrics = [("fusion_rate", "fusion rate"), ("center_stack_load", "center-stack load"),
           ("outboard_load", "outboard load")]
x = np.arange(len(modes)); w = 0.26
shades = ["0.20", "0.50", "0.75"]

fig, ax = plt.subplots(figsize=(5.6, 3.5))
for i, (key, lab) in enumerate(metrics):
    vals = [float(rows[m][key]) for m in modes]
    ax.bar(x + (i - 1) * w, vals, w, label=lab, color=shades[i])
ax.axhline(1.0, color="k", lw=0.8, ls=":")
ax.set_xticks(x); ax.set_xticklabels(modes)
ax.set_ylabel("relative to unpolarized = 1")
ax.legend(frameon=False, loc="upper left", fontsize=8, ncol=1)
ax.annotate("1.70x", (1 + 0*w, 1.70), ha="center", va="bottom", fontsize=8)
ax.annotate("0.85x", (2 + 0*w, 0.85), ha="center", va="bottom", fontsize=8)
finish(fig, "fig_rate_steering")
