#!/usr/bin/env python
"""Zoo-wide polarization law: does the cheap field-coherence C predict the transport
polarization steering effect eta (free-streaming coil-flux change under polarized emission)?
Reads shield_opt/data/eta_C_zoo.csv (64 QA/QH devices from the Ginsburg conformal sweep).
"""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "data/eta_C_zoo.csv")))
C = np.array([float(r["C"]) for r in rows])
perp = np.array([float(r["delta_free_perp"]) for r in rows])
par = np.array([float(r["delta_free_par"]) for r in rows])
cls = np.array([r["symmetry_class"] for r in rows])
nfp = np.array([int(r["nfp"]) for r in rows])
eff = np.abs(perp)

print(f"n={len(rows)}  QA={np.sum(cls=='QA')} QH={np.sum(cls=='QH')}  C[{C.min():.3f},{C.max():.3f}]")
rho, p = spearmanr(C, eff)
rp, pp = pearsonr(C, eff)
print(f"|eta| vs C: Spearman rho={rho:+.3f} (p={p:.2e})  Pearson r={rp:+.3f} (p={pp:.2e})")
lo, hi = C < np.median(C), C >= np.median(C)
print(f"low-C |eta| {eff[lo].mean():.3f}+/-{eff[lo].std():.3f}  high-C {eff[hi].mean():.3f}+/-{eff[hi].std():.3f}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for c, col in (("QA", "#2a6f97"), ("QH", "#e07a2f")):
        m = cls == c
        ax.scatter(C[m], eff[m], s=34, c=col, label=f"{c} (n={m.sum()})", alpha=0.85, edgecolor="k", linewidth=0.3)
    # trend
    z = np.polyfit(C, eff, 1)
    xs = np.linspace(C.min(), C.max(), 50)
    ax.plot(xs, np.polyval(z, xs), "k--", lw=1, alpha=0.6)
    ax.set_xlabel("Field Coherence C (Cheap, Coils-Only Predictor)")
    ax.set_ylabel("Polarization Steering Effect |Δφ| (Free-Streaming Coil Flux)")
    ax.set_title(f"Zoo-Wide Polarization Law (n={len(rows)}): C Predicts η\nSpearman ρ={rho:+.2f}, p={p:.1e}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fp = HERE.parent / "figs/eta_C_zoo.png"
    fp.parent.mkdir(exist_ok=True)
    fig.savefig(fp, dpi=120)
    print("wrote", fp)
except Exception as e:
    print("fig skipped:", e)
