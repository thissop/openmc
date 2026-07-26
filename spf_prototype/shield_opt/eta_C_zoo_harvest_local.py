#!/usr/bin/env python
"""Re-harvest the zoo-wide eta_source(C) law on the COMBINED set:
Ginsburg 64 (data/eta_C_zoo.csv) + new LOCAL done records (sweep_out_local/).
Dedups by device id (Ginsburg row wins if a device appears in both). Writes the
combined CSV and prints the updated Spearman/Pearson n, rho, p. No commits.
"""
import csv, json, glob
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
GINS = HERE / "data/eta_C_zoo.csv"
LOCAL = HERE.parent / "sweep/sweep_out_local"
OUT = HERE / "data/eta_C_zoo_combined_local.csv"

COLS = ["id","symmetry_class","nfp","aspect","C","S_phi",
        "delta_free_perp","delta_free_perp_sd","delta_free_par","delta_free_par_sd","source"]

rows = {}
# 1) Ginsburg committed rows
for r in csv.DictReader(open(GINS)):
    r = dict(r); r["source"] = "ginsburg"
    rows[str(r["id"])] = r

# devices that were run at SMOKE fidelity (not the 200kx20 production law fidelity) and
# must be excluded from the law (their deltas are far noisier).
SMOKE_IDS = {"59509", "932746"}

# 2) new local done records
n_local_new = 0
for f in sorted(glob.glob(str(LOCAL / "config_*_baseline.json"))):
    rec = json.load(open(f))
    if rec.get("status") != "done":
        continue
    did = str(rec["id"])
    if did in SMOKE_IDS:
        continue  # low-fidelity smoke run, not law-grade
    if did in rows:
        continue  # Ginsburg wins on dedup
    perp = rec.get("delta_free_perpendicular"); par = rec.get("delta_free_parallel")
    if perp is None or par is None:
        continue
    rows[did] = dict(
        id=did, symmetry_class=rec.get("symmetry_class",""), nfp=rec.get("nfp",""),
        aspect=round(float(rec.get("aspect",0)),4), C=round(float(rec["C"]),5),
        S_phi=round(float(rec.get("S_phi",0)),5),
        delta_free_perp=round(perp,5), delta_free_perp_sd=round(rec.get("delta_free_perpendicular_sd",0),5),
        delta_free_par=round(par,5), delta_free_par_sd=round(rec.get("delta_free_parallel_sd",0),5),
        source="local")
    n_local_new += 1

# write combined
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
print(f"COMBINED n={n}  (ginsburg={np.sum(src=='ginsburg')}, local_new={n_local_new})")
print(f"  QA={np.sum(cls=='QA')} QH={np.sum(cls=='QH')}  C[{C.min():.3f},{C.max():.3f}]")
rho, p = spearmanr(C, eff); rp, pp = pearsonr(C, eff)
print(f"|eta| vs C: Spearman rho={rho:+.3f} (p={p:.2e})  Pearson r={rp:+.3f} (p={pp:.2e})  n={n}")
lo, hi = C < np.median(C), C >= np.median(C)
print(f"  low-C |eta| {eff[lo].mean():.3f}+/-{eff[lo].std():.3f}  high-C {eff[hi].mean():.3f}+/-{eff[hi].std():.3f}")
# local-only law for a clean cross-check
if n_local_new >= 5:
    ml = src == "local"
    rl, pl = spearmanr(C[ml], eff[ml])
    print(f"  LOCAL-ONLY (n={ml.sum()}): Spearman rho={rl:+.3f} (p={pl:.2e})")
print("wrote", OUT)
