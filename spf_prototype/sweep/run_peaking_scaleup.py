#!/usr/bin/env python
"""Run the analytic peaking-factor study on the flat-S_phi SELECTED device set (~120), in parallel,
and recompute the S_phi -> peaking-modulation correlations at scale. Refreshes the T3 result with
~10x the devices spanning the coherence range.

In : sweep/sweep_out/selected_devices.csv
Out: sweep/sweep_out/geometry_peaking_scaleup.json + figs/geometry_peaking/scaleup_*.png + a printed
     correlation report.
"""
import csv
import sys
import json
from pathlib import Path
from multiprocessing import Pool
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
SEL = HERE / "sweep_out" / "selected_devices.csv"
OUT = HERE / "sweep_out" / "geometry_peaking_scaleup.json"
FIGD = HERE.parent / "figs" / "geometry_peaking"


def one(ID):
    try:
        import geometry_peaking as gp
        rec, _ = gp.analyze_device(ID, keep_maps=False)
        return rec
    except Exception as e:  # noqa
        return {"ID": ID, "error": f"{type(e).__name__}: {e}"}


def report(recs):
    from scipy.stats import spearmanr, mannwhitneyu
    S = np.array([r["S_phi"] for r in recs])
    pperp = np.array([r["PF_perp_over_unpol"] for r in recs])
    ppar = np.array([r["PF_par_over_unpol"] for r in recs])
    nfp = np.array([r["nfp"] for r in recs])
    cls = np.array([r["qs_class"] for r in recs])
    print(f"\n=== SCALE-UP correlations (n={len(recs)}) ===")
    for nm, y in [("PF_perp/unpol", pperp), ("PF_par/unpol", ppar)]:
        rS, pS = spearmanr(S, y); rN, pN = spearmanr(nfp, y)
        print(f"{nm}: vs S_phi rho={rS:+.3f} (p={pS:.1e}) | vs nfp rho={rN:+.3f} (p={pN:.1e})")
    qa, qh = cls == "QA", cls == "QH"
    for nm, y in [("PF_perp/unpol", pperp), ("PF_par/unpol", ppar)]:
        u, p = mannwhitneyu(y[qa], y[qh]) if qa.sum() and qh.sum() else (np.nan, np.nan)
        print(f"{nm}: QA {y[qa].mean():.3f}+/-{y[qa].std():.3f}  QH {y[qh].mean():.3f}+/-{y[qh].std():.3f}"
              f"  (Mann-Whitney p={p:.1e})")
    return S, pperp, ppar, cls


def plot(S, pperp, ppar, cls):
    try:
        import matplotlib; matplotlib.use("Agg")
        try:
            import smplotlib  # noqa
        except Exception:
            pass
        import matplotlib.pyplot as plt
        FIGD.mkdir(parents=True, exist_ok=True)
        col = np.where(cls == "QA", "#2a6fdb", "#e2711d")
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for a, y, t in [(ax[0], pperp, "perpendicular / A"), (ax[1], ppar, "parallel / B,C")]:
            a.scatter(S, y, s=22, c=col, edgecolor="k", linewidths=0.3, alpha=0.85)
            a.axhline(1.0, ls=":", color="grey", lw=0.8)
            a.set_xlabel(r"nematic order $S_\phi$"); a.set_ylabel(f"PF({t.split()[0]}) / PF(unpol)")
            a.set_title(t)
        fig.suptitle(f"SPF peaking modulation vs field coherence (n={S.size}; blue QA, orange QH)")
        fig.tight_layout(); fig.savefig(FIGD / "scaleup_pf_vs_Sphi.png", dpi=180); plt.close(fig)
        print(f"plot -> figs/geometry_peaking/scaleup_pf_vs_Sphi.png")
    except Exception as e:
        print(f"plot skipped: {type(e).__name__}: {e}")


def main(nproc=6):
    ids = [int(r["ID"]) for r in csv.DictReader(open(SEL))]
    print(f"peaking scale-up: {len(ids)} devices on {nproc} procs ...", flush=True)
    with Pool(nproc) as p:
        out = p.map(one, ids)
    recs = [r for r in out if "error" not in r]
    print(f"ok {len(recs)}/{len(out)} (failed {len(out)-len(recs)})")
    OUT.write_text(json.dumps(recs, indent=2))
    S, pperp, ppar, cls = report(recs)
    plot(S, pperp, ppar, cls)
    print("wrote", OUT)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
