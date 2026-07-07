#!/usr/bin/env python
"""Compute the cheap field-coherence stats (S_phi, C, lambda_phi, reversal_frac) for every device in
the candidate pool -- b_hat from the REAL coils (Biot-Savart), auto-downloading serials on demand.
Parallel. This is the CHEAP predictor pass (no wall load); the flat-S_phi down-select and the
expensive peaking run follow.

In : sweep/sweep_out/candidate_pool.csv
Out: sweep/sweep_out/sphi_pool.csv  (ID, class, nfp, S_phi, C, lambda_phi, reversal_frac, ok)
"""
import csv
import sys
from pathlib import Path
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))
POOL = HERE / "sweep_out" / "candidate_pool.csv"
OUT = HERE / "sweep_out" / "sphi_pool.csv"


def one(rec):
    ID, cls, nfp = rec
    try:
        import quasr_loader as ql
        import geometry_peaking as gp
        dev = ql.load_device(ID)
        xyz, w, rho, bhat = gp.sample_source(dev)
        s = gp.cheap_stats(xyz, w, bhat, dev.meta, dev)
        return dict(ID=ID, cls=cls, nfp=nfp, S_phi=s["S_phi"], C=s["C"],
                    lambda_phi=s["lambda_phi"], reversal_frac=s.get("reversal_frac", 0.0), ok=1)
    except Exception as e:  # noqa
        return dict(ID=ID, cls=cls, nfp=nfp, S_phi="", C="", lambda_phi="",
                    reversal_frac="", ok=0, err=f"{type(e).__name__}:{e}")


def main(nproc=8):
    recs = [(int(r["ID"]), r["class"], int(r["nfp"]))
            for r in csv.DictReader(open(POOL))]
    print(f"computing S_phi for {len(recs)} candidates on {nproc} procs ...", flush=True)
    with Pool(nproc) as p:
        out = p.map(one, recs)
    ok = [o for o in out if o["ok"]]
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "S_phi", "C", "lambda_phi", "reversal_frac", "ok"])
        for o in out:
            w.writerow([o["ID"], o["cls"], o["nfp"],
                        f"{o['S_phi']:.4f}" if o["ok"] else "", f"{o['C']:.4f}" if o["ok"] else "",
                        f"{o['lambda_phi']:.4f}" if o["ok"] else "",
                        f"{o['reversal_frac']:.4f}" if o["ok"] else "", o["ok"]])
    sp = sorted(o["S_phi"] for o in ok)
    print(f"ok {len(ok)}/{len(out)}  (failed {len(out)-len(ok)})")
    if ok:
        import numpy as np
        a = np.array(sp)
        print(f"S_phi span [{a.min():.3f}, {a.max():.3f}]  "
              f"quartiles {np.percentile(a,25):.3f}/{np.percentile(a,50):.3f}/{np.percentile(a,75):.3f}")
        for lo, hi in [(0, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
            print(f"  S_phi in [{lo:.2f},{hi:.2f}): {int(((a>=lo)&(a<hi)).sum())}")
    print("wrote", OUT)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
