"""Fig: free-streaming analytic vs OpenMC on real QUASR QA 59509, per wall.
Outboard/floor/ceiling agree <=1.8%; the strongly-curved inboard shows the bounded 4-6% residual."""
import csv
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

rows = list(csv.DictReader(r for r in open(DATA / "quasr59509_validation.csv") if not r.startswith("#")))
walls = [r["wall"] for r in rows]
A_mc = np.array([float(r["A_iso_MC"]) for r in rows])
A_an = np.array([float(r["A_iso_ana"]) for r in rows])
B_mc = np.array([float(r["B_iso_MC"]) for r in rows])
B_an = np.array([float(r["B_iso_ana"]) for r in rows])

fig, ax = plt.subplots(figsize=(4.4, 4.2))
lim = [0.70, 1.30]
ax.plot(lim, lim, "-", color="0.7", lw=0.9, zorder=0)
for mc, an, mk, lab in [(A_mc, A_an, "o", "A (perp)"), (B_mc, B_an, "s", "B/C (par)")]:
    ax.scatter(an, mc, marker=mk, s=45, facecolor="none", edgecolor="0.15", label=lab)
for w, a, m in zip(walls, A_an, A_mc):
    ax.annotate(w, (a, m), (4, -2), textcoords="offset points", fontsize=7, color="0.4")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("analytic  load / isotropic")
ax.set_ylabel("OpenMC  load / isotropic")
ax.set_title("Analytic vs Monte Carlo: QUASR QA 59509")
ax.legend(frameon=False, loc="upper left", fontsize=8)
finish(fig, "fig_quasr_validation")
