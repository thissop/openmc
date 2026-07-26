#!/usr/bin/env python
"""Stratified candidate pool from the full QUASR catalogue (371k devices), chosen to SPAN the
field-coherence range rather than cluster. We can't compute S_phi for all 371k, so we stratify on
the catalogue columns that drive S_phi -- class (helicity), nfp, quasisymmetry error, aspect ratio --
and sample evenly across strata. S_phi is then computed for the pool (compute_sphi_pool.py) and the
final ~120 are down-selected to a FLAT S_phi histogram (stratified_select.py).

Out: sweep/sweep_out/candidate_pool.csv  (ID, class, nfp, qs_error, aspect_ratio, mean_iota)
"""
import gzip
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CAT = HERE.parent / "data" / "quasr" / "catalogue.csv.gz"
OUT = HERE / "sweep_out" / "candidate_pool.csv"
N_PER_STRATUM = 25          # sampled per (class x nfp x qs_error-quintile) stratum
RNG = np.random.default_rng(20260707)


def _f(r, k):
    try:
        return float(r[k])
    except Exception:
        return None


def main():
    rows = [r for r in csv.DictReader(gzip.open(CAT, "rt"))]
    for r in rows:
        r["_cls"] = "QA" if int(float(r["helicity"])) == 0 else "QH"
        r["_nfp"] = int(float(r["nfp"]))
        r["_qs"] = _f(r, "qs_error")
        r["_asp"] = _f(r, "aspect_ratio")
    rows = [r for r in rows if r["_qs"] is not None and r["_asp"] is not None]
    # reactor-relevance: cap aspect ratio (ARIES-CS ~4.5, HELIAS ~10); drop compact
    # research configs with extreme aspect that are not engineering-relevant.
    rows = [r for r in rows if 2.5 <= r["_asp"] <= 12.0]

    # qs_error quintile edges (global) -> a coherence proxy axis
    qs = np.array([r["_qs"] for r in rows])
    edges = np.quantile(qs, np.linspace(0, 1, 6))
    edges[-1] += 1e-9

    def qbin(q):
        return int(np.clip(np.searchsorted(edges, q, side="right") - 1, 0, 4))

    strata = {}
    for r in rows:
        key = (r["_cls"], r["_nfp"], qbin(r["_qs"]))
        strata.setdefault(key, []).append(r)

    picked = []
    for key, group in sorted(strata.items()):
        idx = RNG.choice(len(group), size=min(N_PER_STRATUM, len(group)), replace=False)
        picked.extend(group[i] for i in idx)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "qs_error", "aspect_ratio", "mean_iota"])
        for r in picked:
            w.writerow([int(float(r["ID"])), r["_cls"], r["_nfp"],
                        f"{r['_qs']:.4f}", f"{r['_asp']:.3f}", r.get("mean_iota", "")])

    cls = [r["_cls"] for r in picked]
    print(f"candidate pool: {len(picked)} devices "
          f"(QA {cls.count('QA')}, QH {cls.count('QH')}) across {len(strata)} strata")
    import collections
    print("by nfp:", dict(sorted(collections.Counter(r["_nfp"] for r in picked).items())))
    print(f"qs_error span: [{min(r['_qs'] for r in picked):.2f}, {max(r['_qs'] for r in picked):.2f}]  "
          f"aspect span: [{min(r['_asp'] for r in picked):.1f}, {max(r['_asp'] for r in picked):.1f}]")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
