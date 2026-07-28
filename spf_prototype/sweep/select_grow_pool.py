#!/usr/bin/env python
"""Fresh candidate pool to GROW the eta(C) zoo law past n=64.

Selects engineering-relevant QUASR devices from the full catalogue (371k) that are
NOT already in the exhausted 178-device production_scaleup.json manifest (= 64 done +
114 known-intrinsic failures) and NOT in the committed eta_C_zoo.csv (64). Same
reactor-relevance cuts as select_candidate_pool.py (aspect in [2.5,12]) plus a light
coil-complexity relevance signal, stratified across QS-class x nfp x coherence
(qs_error quintile, the cheap coherence proxy) so the fresh sample SPANS C.

Out: sweep/sweep_out_grow/grow_candidate_pool.csv
"""
import gzip
import csv
import json
from pathlib import Path
import collections
import numpy as np

HERE = Path(__file__).resolve().parent
CAT = HERE.parent / "data" / "quasr" / "catalogue.csv.gz"
PROD = HERE / "configs" / "production_scaleup.json"
ZOO = HERE.parent / "shield_opt" / "data" / "eta_C_zoo.csv"
OUTDIR = HERE / "sweep_out_grow"
OUT = OUTDIR / "grow_candidate_pool.csv"

N_PER_STRATUM = 14          # per (class x nfp x qs_error-quintile) stratum
TARGET = 220
RNG = np.random.default_rng(20260726)


def _f(r, k):
    try:
        return float(r[k])
    except Exception:
        return None


def excluded_ids():
    prod = json.loads(PROD.read_text())
    ids = {int(c["id"]) for c in prod["configs"]}
    for r in csv.DictReader(open(ZOO)):
        ids.add(int(r["id"]))
    return ids


def main():
    excl = excluded_ids()
    rows = [r for r in csv.DictReader(gzip.open(CAT, "rt"))]
    for r in rows:
        try:
            r["_id"] = int(float(r["ID"]))
            r["_cls"] = "QA" if int(float(r["helicity"])) == 0 else "QH"
            r["_nfp"] = int(float(r["nfp"]))
        except Exception:
            r["_id"] = None
        r["_qs"] = _f(r, "qs_error")
        r["_asp"] = _f(r, "aspect_ratio")
        r["_nchp"] = _f(r, "nc_per_hp")
        r["_tcl"] = _f(r, "total_coil_length")
        r["_kappa"] = _f(r, "max_kappa")
        r["_c2s"] = _f(r, "min_coil2surface_dist")
    rows = [r for r in rows if r["_id"] is not None and r["_qs"] is not None
            and r["_asp"] is not None and r["_id"] not in excl]
    # reactor-relevance: aspect band (ARIES-CS ~4.5, HELIAS ~10) as in select_candidate_pool
    rows = [r for r in rows if 2.5 <= r["_asp"] <= 12.0]
    # engineering-relevant nfp band (production coverage spanned nfp 1-6)
    rows = [r for r in rows if 1 <= r["_nfp"] <= 6]
    # light coil-complexity relevance: buildable coils (>=1 coil/hp; finite standoff).
    # This keeps the "does a blanket fit / are the coils real" flavour without needing
    # the full Kappel standoff computed per device (RESULTS_engineering_relevance.md).
    rows = [r for r in rows if (r["_nchp"] or 0) >= 1
            and (r["_c2s"] is None or r["_c2s"] > 0.05)]
    print(f"catalogue after cuts + exclusion: {len(rows)} devices")

    qs = np.array([r["_qs"] for r in rows])
    edges = np.quantile(qs, np.linspace(0, 1, 6))
    edges[-1] += 1e-9

    def qbin(q):
        return int(np.clip(np.searchsorted(edges, q, side="right") - 1, 0, 4))

    strata = {}
    for r in rows:
        strata.setdefault((r["_cls"], r["_nfp"], qbin(r["_qs"])), []).append(r)

    picked = []
    for key, group in sorted(strata.items()):
        idx = RNG.choice(len(group), size=min(N_PER_STRATUM, len(group)), replace=False)
        picked.extend(group[i] for i in idx)
    # if we overshoot TARGET, thin uniformly at random preserving stratification
    if len(picked) > TARGET:
        keep = RNG.choice(len(picked), size=TARGET, replace=False)
        picked = [picked[i] for i in sorted(keep)]

    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "qs_error", "aspect_ratio",
                    "nc_per_hp", "total_coil_length", "min_coil2surface_dist"])
        for r in picked:
            w.writerow([r["_id"], r["_cls"], r["_nfp"], f"{r['_qs']:.4f}",
                        f"{r['_asp']:.3f}", int(r["_nchp"] or 0),
                        f"{r['_tcl']:.3f}" if r['_tcl'] else "",
                        f"{r['_c2s']:.4f}" if r['_c2s'] else ""])

    cls = [r["_cls"] for r in picked]
    print(f"grow pool: {len(picked)} devices (QA {cls.count('QA')}, QH {cls.count('QH')}) "
          f"across {len(strata)} strata")
    print("by nfp:", dict(sorted(collections.Counter(r["_nfp"] for r in picked).items())))
    print(f"qs_error span [{min(r['_qs'] for r in picked):.2f},{max(r['_qs'] for r in picked):.2f}]  "
          f"aspect span [{min(r['_asp'] for r in picked):.1f},{max(r['_asp'] for r in picked):.1f}]")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
