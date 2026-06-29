#!/usr/bin/env python
"""Build a local CSV catalogue of the QUASR stellarator database and shortlist
gentle (low-eps_eff) quasi-axisymmetric (QA) candidates for the free-streaming
validation -- WITHOUT the 13 GB bulk download.

The QUASR navigator (quasr.flatironinstitute.org) is a SPA backed by a single
gzipped table `database.json.gz` (~37 MB, 371k devices) and per-device simsopt
serial files at  simsopt_serials/<ID[:4]>/serial<ID>.json  (~0.2 MB each).  We
pull the table once -> our own CSV; the metadata has no boundary harmonics, so it
cannot give the geometric eps_eff directly, but it is plenty to RANK and SHORTLIST.
The shortlisted boundaries are downloaded and eps_eff computed exactly elsewhere
(quasr_eps_eff.py).

Columns of interest: helicity (0=QA, 1=QH), nfp, aspect_ratio, mean/max_elongation,
qs_error (log10 of the root QS error; lower is better), mean_iota, minor_radius.

Run: $HOME/desc_venv/bin/python spf_prototype/python/quasr_catalogue.py
"""
import gzip
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BASE = "https://quasr.flatironinstitute.org/"
DATA = Path(__file__).resolve().parents[1] / "data" / "quasr"
GZ = DATA / "database.json.gz"
CSV = DATA / "catalogue.csv.gz"          # our full catalogue (compressed)
SHORT = DATA / "shortlist.csv"           # candidates to download + eps_eff


def fetch_catalogue():
    DATA.mkdir(parents=True, exist_ok=True)
    if not GZ.exists():
        print(f"downloading {BASE}database.json.gz ...")
        urllib.request.urlretrieve(BASE + "database.json.gz", GZ)
    d = json.load(gzip.open(GZ))
    df = pd.DataFrame(d["data"], columns=d["columns"])
    print(f"catalogue: {len(df)} devices, {len(df.columns)} columns")
    return df


def shortlist(df, n_top=40, n_rand=20, seed=0):
    """QA, nfp in {2,3}, reactor-ish aspect, good QS, multi-surface. Rank by
    max_elongation (roundest cross-sections == gentlest shaping proxy); also keep a
    random sample of the filtered pool to MEASURE how well the proxy predicts the
    geometric eps_eff once we download boundaries."""
    qa = df[(df.helicity == 0)
            & df.nfp.isin([2, 3])
            & df.aspect_ratio.between(4.0, 12.0)
            & (df.qs_error < -2.0)            # root QS error < 1e-2
            & (df.Nsurfaces >= 2)].copy()
    print(f"QA filtered pool: {len(qa)} devices "
          f"(nfp2={int((qa.nfp==2).sum())}, nfp3={int((qa.nfp==3).sum())})")
    qa = qa.sort_values("max_elongation")
    top = qa.head(n_top)
    rng = np.random.default_rng(seed)
    rest = qa.iloc[n_top:]
    rand = rest.iloc[rng.choice(len(rest), min(n_rand, len(rest)), replace=False)] \
        if len(rest) else rest
    sel = pd.concat([top, rand]).drop_duplicates("ID")
    sel = sel[["ID", "nfp", "aspect_ratio", "mean_elongation", "max_elongation",
               "qs_error", "mean_iota", "minor_radius", "volume", "Nsurfaces"]]
    sel["pick_reason"] = ["low_elong"] * len(top) + ["random"] * (len(sel) - len(top))
    return sel


if __name__ == "__main__":
    df = fetch_catalogue()
    df.to_csv(CSV, index=False)
    print(f"wrote {CSV}  ({CSV.stat().st_size/1e6:.1f} MB)")
    sel = shortlist(df)
    sel.to_csv(SHORT, index=False)
    print(f"\nshortlist: {len(sel)} devices -> {SHORT}")
    print(sel.head(12).to_string(index=False))
