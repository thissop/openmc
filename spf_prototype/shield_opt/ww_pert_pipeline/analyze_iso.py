#!/usr/bin/env python
"""Compare baseline vs perturbed per-coil results, write perturbation_result.npz + figures."""
import os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WD = os.environ.get("WORKDIR", "/burg-archive/home/tjk2147/pstl_test/kiso")
os.makedirs(f"{WD}/figs", exist_ok=True)

def load(tag):
    d = np.load(f"{WD}/coil_{tag}.npz", allow_pickle=True)
    return d

cen = np.load(f"{WD}/coil_centroids.npz")
cell_ids = cen["cell_ids"]; phi = cen["phi"]; Rc = cen["R"]; zc = cen["cz"]

def relerr(d, key, ekey):
    return d[ekey] / np.maximum(d[key], 1e-30)

pairs = []  # (label, base_tag, pert_tag)
# figure which tags exist
import os
def have(tag): return os.path.exists(f"{WD}/coil_{tag}.npz")

configs = []
for mode in ["unpol", "A"]:
    bt, pt = f"base_{mode}", f"pert_{mode}"
    if have(bt) and have(pt):
        configs.append((mode, bt, pt))

print("Configs:", [c[0] for c in configs])
results = {}
for mode, bt, pt in configs:
    b = load(bt); p = load(pt)
    bf = b["coil_fast_flux"]; be = b["coil_fast_flux_err"]
    pf = p["coil_fast_flux"]; pe = p["coil_fast_flux_err"]
    bpk = bf.max()/bf.mean(); ppk = pf.max()/pf.mean()
    ihot_b = int(np.argmax(bf))
    # track same hottest cells in perturbed
    frac = pf/bf; frac_err = frac*np.sqrt((pe/bf)**2 + (be*pf/bf**2)**2)
    print(f"\n===== MODE {mode} =====")
    print(f"  peaking baseline={bpk:.3f}  perturbed={ppk:.3f}")
    print(f"  hottest baseline cell={cell_ids[ihot_b]} flux={bf[ihot_b]:.4e}+/-{be[ihot_b]:.1e} "
          f"-> perturbed {pf[ihot_b]:.4e}+/-{pe[ihot_b]:.1e}  ratio={frac[ihot_b]:.3f}+/-{frac_err[ihot_b]:.3f}")
    print(f"  TBR baseline={float(b['tbr']):.4f}+/-{float(b['tbr_err']):.4f}  "
          f"perturbed={float(p['tbr']):.4f}+/-{float(p['tbr_err']):.4f}  "
          f"dTBR={float(p['tbr'])-float(b['tbr']):+.4f}")
    # per-coil table
    print("  cell phi_deg  base_flux  pert_flux  ratio")
    order = np.argsort(-bf)
    for i in order:
        print(f"   {cell_ids[i]:3d} {phi[i]:6.1f}  {bf[i]:.3e} {pf[i]:.3e}  {frac[i]:.3f}+/-{frac_err[i]:.3f}")
    results[mode] = dict(
        cell_ids=cell_ids, phi=phi, Rc=Rc, zc=zc,
        base_flux=bf, base_flux_err=be, pert_flux=pf, pert_flux_err=pe,
        base_heat=b["coil_heating"], base_heat_err=b["coil_heating_err"],
        pert_heat=p["coil_heating"], pert_heat_err=p["coil_heating_err"],
        base_peaking=bpk, pert_peaking=ppk,
        base_tbr=float(b["tbr"]), base_tbr_err=float(b["tbr_err"]),
        pert_tbr=float(p["tbr"]), pert_tbr_err=float(p["tbr_err"]),
        base_hist=int(b["batches"])*int(b["particles"]),
        pert_hist=int(p["batches"])*int(p["particles"]))

# ---- figure: per-coil flux baseline vs perturbed (unpol) ----
if "unpol" in results:
    r = results["unpol"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    x = np.arange(len(cell_ids))
    ax.bar(x-0.2, r["base_flux"], 0.4, yerr=r["base_flux_err"], label="Baseline", color="#4C78A8")
    ax.bar(x+0.2, r["pert_flux"], 0.4, yerr=r["pert_flux_err"], label="Perturbed", color="#F58518")
    ax.set_xticks(x); ax.set_xticklabels(cell_ids, fontsize=7)
    ax.set_xlabel("Coil Cell ID"); ax.set_ylabel("Fast Flux (>0.1 MeV) per Source Neutron")
    ax.set_title("Per-Coil Fast Flux: Baseline vs Perturbed (Unpolarized)")
    ax.legend()
    fig.tight_layout(); fig.savefig(f"{WD}/figs/coil_flux_compare.png", dpi=130)
    print("\nSAVED figs/coil_flux_compare.png")

np.savez(f"{WD}/perturbation_result.npz",
         **{f"{m}_{k}": v for m, d in results.items() for k, v in d.items()},
         configs=np.array([c[0] for c in configs]))
print("SAVED perturbation_result.npz")
