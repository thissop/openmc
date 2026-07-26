"""Robust cached-only S_phi computation. Processes the candidate pool in batches with
a fresh Pool per batch and INCREMENTAL writes, so a worker segfault (which deadlocks a
single big Pool.map) can at worst lose one batch, which is then retried serially.
No network: only devices whose serial is already cached are processed."""
import csv, sys, os
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, "."); sys.path.insert(0, "../python")
SER = Path("../data/quasr/serials")
OUT = Path("sweep_out/sphi_pool.csv")
HDR = ["ID", "class", "nfp", "S_phi", "C", "lambda_phi", "reversal_frac", "ok"]
BATCH = 64


def _row(o):
    return [o["ID"], o["cls"], o["nfp"],
            f"{o['S_phi']:.4f}" if o["ok"] else "",
            f"{o['C']:.4f}" if o["ok"] else "",
            f"{o['lambda_phi']:.4f}" if o["ok"] else "",
            f"{o['reversal_frac']:.4f}" if o["ok"] else "", o["ok"]]


def main():
    from compute_sphi_pool import one
    pool = [r for r in csv.DictReader(open("sweep_out/candidate_pool.csv"))]
    recs = [(int(r["ID"]), r["class"], int(r["nfp"])) for r in pool
            if (SER / f"serial{int(r['ID']):07d}.json").exists()]
    print(f"computing S_phi for {len(recs)} CACHED pool devices in batches of {BATCH}...",
          flush=True)
    with open(OUT, "w", newline="") as f:
        csv.writer(f).writerow(HDR)
    done = 0
    for i in range(0, len(recs), BATCH):
        batch = recs[i:i + BATCH]
        try:
            with Pool(8) as p:
                out = p.map(one, batch)
        except Exception as e:  # pool broke (worker crash) -> serial fallback
            print(f"  batch {i//BATCH}: pool failed ({type(e).__name__}); serial fallback",
                  flush=True)
            out = [one(r) for r in batch]
        with open(OUT, "a", newline="") as f:
            w = csv.writer(f)
            for o in out:
                w.writerow(_row(o))
        done += len(out)
        print(f"  {done}/{len(recs)} done", flush=True)
    # summary
    rows = [r for r in csv.DictReader(open(OUT)) if r["ok"] == "1" and r["S_phi"]]
    import numpy as np
    a = np.array(sorted(float(r["S_phi"]) for r in rows))
    print(f"ok {len(rows)}/{len(recs)}  S_phi span [{a.min():.3f},{a.max():.3f}] "
          f"med {np.median(a):.3f}", flush=True)
    for lo, hi in [(-1, 0), (0, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
        print(f"  S_phi in [{lo},{hi}): {int(((a >= lo) & (a < hi)).sum())}")


if __name__ == "__main__":
    main()
