#!/usr/bin/env python
"""Download the shortlisted QUASR boundaries (VMEC namelists) and compute the
geometric eps_eff EXACTLY, to pick a genuinely gentle (low-eps_eff) published QA
config for the free-streaming validation -- without the bulk database and without
simsopt's compiled core (which SIGILLs on this CPU).

Per device we fetch only the small VMEC input namelist
  nml/<zpad7(ID)[:4]>/input.<zpad7(ID)>   (~18 kB)
parse its RBC/ZBS boundary spectrum, and measure the helical excursion of each
constant-theta filament (see quasr_geom.eps_eff):

  e_geom        = max_theta (max_phi |loop-<loop>|) / a   intrinsic, aspect-invariant
  eps_eff(gap)  = max_theta excursion / standoff to a square wall gap*a outside bbox

The metadata proxy (max_elongation) is reported alongside to show how well it
predicted gentleness.

Run: $HOME/desc_venv/bin/python spf_prototype/python/quasr_eps_eff.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import quasr_geom as G  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "data" / "quasr"


def main():
    sl = pd.read_csv(DATA / "shortlist.csv")
    rows = []
    for _, r in sl.iterrows():
        ID = int(r.ID)
        try:
            dev = G.load_device(ID, r.mean_iota)
            eg, ee, a = G.eps_eff(dev)
            rows.append(dict(ID=ID, nfp=int(r.nfp), aspect=r.aspect_ratio,
                             max_elong=r.max_elongation, qs_log10=r.qs_error,
                             iota=r.mean_iota, e_geom=eg, eps_eff_gap25=ee,
                             R0=dev.R0, a=a, reason=r.pick_reason))
            print(f"{ID:7d} nfp{int(r.nfp)} asp{r.aspect_ratio:5.1f} "
                  f"elong{r.max_elongation:5.2f} -> e_geom {eg:.3f}  eps_eff {ee:.3f}",
                  flush=True)
        except Exception as e:
            print(f"{ID:7d} SKIP: {str(e)[:70]}", flush=True)
    out = pd.DataFrame(rows).sort_values("e_geom").reset_index(drop=True)
    out.to_csv(DATA / "eps_eff_results.csv", index=False)
    print(f"\nwrote {DATA/'eps_eff_results.csv'}  ({len(out)} devices)")
    c = np.corrcoef(out.max_elong, out.e_geom)[0, 1]
    print(f"proxy corr(max_elongation, e_geom) = {c:.2f}")
    print("\ngentlest 8 by e_geom:")
    print(out.head(8)[["ID", "nfp", "aspect", "max_elong", "iota", "e_geom",
                       "eps_eff_gap25", "reason"]].to_string(index=False))


if __name__ == "__main__":
    main()
