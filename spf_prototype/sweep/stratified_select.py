"""Down-select ~TARGET devices from the S_phi pool so the S_phi histogram is roughly FLAT across the
range within each class -- the anti-bias step. QA clusters near S_phi~1, so we cap the dense bins and
keep all of the rare low-S_phi devices; the result spans the coherence range instead of piling up at
high S_phi.

In : sweep/sweep_out/sphi_pool.csv
Out: sweep/sweep_out/selected_devices.csv  (ID, class, nfp, S_phi, C, lambda_phi, reversal_frac)
"""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
POOL = HERE / "sweep_out" / "sphi_pool.csv"
OUT = HERE / "sweep_out" / "selected_devices.csv"
TARGET = 120
NBINS = 8
RNG = np.random.default_rng(7)


def main():
    recs = [r for r in csv.DictReader(open(POOL)) if r["ok"] == "1" and r["S_phi"]]
    for r in recs:
        r["S_phi"] = float(r["S_phi"])
    s = np.array([r["S_phi"] for r in recs])
    lo, hi = s.min(), min(1.0, s.max())
    edges = np.linspace(lo, hi + 1e-9, NBINS + 1)
    # per (class, S_phi-bin) cap so the joint histogram is ~flat; ~TARGET/(2*NBINS) each
    cap = max(1, round(TARGET / (2 * NBINS)))
    sel = []
    for cls in ("QA", "QH"):
        for b in range(NBINS):
            grp = [r for r in recs if r["class"] == cls
                   and edges[b] <= r["S_phi"] < edges[b + 1]]
            if not grp:
                continue
            idx = RNG.choice(len(grp), size=min(cap, len(grp)), replace=False)
            sel.extend(grp[i] for i in idx)

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "S_phi", "C", "lambda_phi", "reversal_frac"])
        for r in sorted(sel, key=lambda r: r["S_phi"]):
            w.writerow([r["ID"], r["class"], r["nfp"], f"{r['S_phi']:.4f}",
                        r["C"], r["lambda_phi"], r["reversal_frac"]])

    ss = np.array([r["S_phi"] for r in sel])
    cls = [r["class"] for r in sel]
    print(f"selected {len(sel)} devices (QA {cls.count('QA')}, QH {cls.count('QH')}) from {len(recs)} pool")
    print("S_phi histogram of the SELECTED set (should be ~flat):")
    hcounts, _ = np.histogram(ss, bins=edges)
    for b in range(NBINS):
        print(f"  [{edges[b]:.2f},{edges[b+1]:.2f}): {hcounts[b]}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
