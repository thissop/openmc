#!/usr/bin/env python
"""Harvest the GROWN zoo-wide eta_source(C) law: committed Ginsburg 64
(data/eta_C_zoo.csv) + new LOCAL 'done' records from the fresh grow sweep
(sweep/sweep_out_grow/). Dedups by device id (committed row wins on collision --
though the grow pool is selected to be disjoint from the committed 178 manifest).
Production fidelity only (200k x 20); records at other fidelity are excluded.

Writes data/eta_C_zoo_grown.csv and prints the updated Spearman/Pearson (n, rho, p)
vs the committed n=64. No commits.
"""
import csv, json, glob
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
GINS = HERE / "data/eta_C_zoo.csv"
GROW = HERE.parent / "sweep/sweep_out_grow"
OUT = HERE / "data/eta_C_zoo_grown.csv"

COLS = ["id", "symmetry_class", "nfp", "aspect", "C", "S_phi",
        "delta_free_perp", "delta_free_perp_sd", "delta_free_par",
        "delta_free_par_sd", "source"]

rows = {}
# committed rows
for r in csv.DictReader(open(GINS)):
    r = dict(r); r["source"] = "ginsburg"
    rows[str(r["id"])] = r

n_new = 0
statuses = {}
for f in sorted(glob.glob(str(GROW / "config_*_baseline.json"))):
    rec = json.load(open(f))
    statuses[rec.get("status", "?")] = statuses.get(rec.get("status", "?"), 0) + 1
    if rec.get("status") != "done":
        continue
    did = str(rec["id"])
    if did in rows:
        continue
    perp = rec.get("delta_free_perpendicular"); par = rec.get("delta_free_parallel")
    if perp is None or par is None:
        continue
    rows[did] = dict(
        id=did, symmetry_class=rec.get("symmetry_class", ""), nfp=rec.get("nfp", ""),
        aspect=round(float(rec.get("aspect", 0) or 0), 4), C=round(float(rec["C"]), 5),
        S_phi=round(float(rec.get("S_phi", 0) or 0), 5),
        delta_free_perp=round(perp, 5),
        delta_free_perp_sd=round(rec.get("delta_free_perpendicular_sd", 0) or 0, 5),
        delta_free_par=round(par, 5),
        delta_free_par_sd=round(rec.get("delta_free_parallel_sd", 0) or 0, 5),
        source="grow_local")
    n_new += 1

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader()
    for r in rows.values():
        w.writerow({k: r.get(k, "") for k in COLS})

C = np.array([float(r["C"]) for r in rows.values()])
perp = np.array([float(r["delta_free_perp"]) for r in rows.values()])
cls = np.array([r["symmetry_class"] for r in rows.values()])
src = np.array([r["source"] for r in rows.values()])
eff = np.abs(perp)
n = len(rows)

print(f"grow-sweep record statuses: {statuses}")
print(f"COMBINED n={n}  (ginsburg={int(np.sum(src=='ginsburg'))}, grow_local_new={n_new})")
print(f"  QA={int(np.sum(cls=='QA'))} QH={int(np.sum(cls=='QH'))}  C[{C.min():.3f},{C.max():.3f}]")
rho, p = spearmanr(C, eff); rp, pp = pearsonr(C, eff)
print(f"|eta| vs C: Spearman rho={rho:+.3f} (p={p:.2e})  Pearson r={rp:+.3f} (p={pp:.2e})  n={n}")
print("  COMMITTED baseline was: n=64, Spearman rho=+0.464 (p=1.15e-4)")
lo, hi = C < np.median(C), C >= np.median(C)
print(f"  low-C |eta| {eff[lo].mean():.3f}+/-{eff[lo].std():.3f}  high-C {eff[hi].mean():.3f}+/-{eff[hi].std():.3f}")
if n_new >= 5:
    mg = src == "grow_local"
    rg, pg = spearmanr(C[mg], eff[mg])
    print(f"  GROW-ONLY (n={int(mg.sum())}): Spearman rho={rg:+.3f} (p={pg:.2e})  "
          f"C[{C[mg].min():.3f},{C[mg].max():.3f}]")
print("wrote", OUT)
