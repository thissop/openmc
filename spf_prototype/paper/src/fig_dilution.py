"""Fig: scattering dilutes the free-streaming steering ~3x (the headline honesty result)."""
import csv
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

rows = list(csv.DictReader(r for r in open(DATA / "steering_dilution.csv") if not r.startswith("#")))
labels = [f"{r['mode']}\n{r['wall'].title()}" for r in rows]
fs = [float(r["free_streaming_pct"]) for r in rows]
sc = [float(r["with_scattering_pct"]) for r in rows]
x = np.arange(len(rows)); w = 0.38

fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.bar(x - w/2, fs, w, label="Free-Streaming", color="0.30")
ax.bar(x + w/2, sc, w, label="Full Transport", color="0.72")
ax.axhline(0, color="k", lw=0.8)
for xi, a, b in zip(x, fs, sc):
    ax.annotate(f"{a:+.0f}", (xi - w/2, a), ha="center",
                va="bottom" if a >= 0 else "top", fontsize=8)
    ax.annotate(f"{b:+.0f}", (xi + w/2, b), ha="center",
                va="bottom" if b >= 0 else "top", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Steering, Load vs Isotropic [%]")
ax.legend(frameon=False, loc="lower right", fontsize=8)
finish(fig, "fig_dilution")
