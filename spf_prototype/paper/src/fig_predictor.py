"""Fig: across 15 conformal-wall devices, SPF global benefit Y is weak everywhere (0-16%,
median 4%) and not robustly predictable from geometry: the transport-free headroom proxy
(PF_geo) is the best predictor (rho=+0.53, borderline at n=15) while pure-plasma shape
(nematic Q) does not predict it (rho=+0.40, not significant). Data: predictor_n15.csv."""
import csv
import numpy as np
import matplotlib.pyplot as plt
from _style import DATA, finish

rows = list(csv.DictReader(r for r in open(DATA / "predictor_n15.csv") if not r.startswith("#")))
Y = np.array([float(r["Y"]) * 100 for r in rows])
PFg = np.array([float(r["PF_geo"]) for r in rows])
Q = np.array([float(r["nematic_Q"]) for r in rows])

fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4))
ax[0].scatter(PFg, Y, s=42, facecolor="none", edgecolor="0.15")
ax[0].set_xlabel("headroom proxy PF$_\\mathrm{geo}$ (geometry only)")
ax[0].set_ylabel("SPF benefit Y  [% peaking reduction]")
ax[1].scatter(Q, Y, s=42, marker="s", facecolor="none", edgecolor="0.4")
ax[1].set_xlabel("plasma shape (nematic Q)")
for a in ax:
    a.axhline(0, color="0.8", lw=0.8)
    a.set_ylim(-1, 18)
finish(fig, "fig_predictor")
