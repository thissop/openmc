#!/usr/bin/env python
"""Gate 2: does the SSR run (frozen first-wall source through the full build) reproduce the full-run
coil dose? Compares per-coil fast flux from the master (full run) vs the SSR run.

Verdict logic: if the per-coil RATIO full/ssr is CONSTANT across coils, SSR reproduces the spatial
pattern and the decoupling is valid -- the constant is just the source normalization k (crossings per
birth). A non-constant ratio (high coefficient of variation) means SSR distorts the pattern and the
decoupling is not yet trustworthy (likely the interior-surface backscatter double-count; refine by
killing plasma-void re-entry).

Usage: gate2_analyze.py coil_real_master.npz coil_real_gate2.npz
"""
import sys
import numpy as np
from scipy.stats import spearmanr

full = np.load(sys.argv[1]); ssr = np.load(sys.argv[2])
cf = full["coil_fast_flux"]; cs = ssr["coil_fast_flux"]
cells = full["coil_cells"]
assert np.array_equal(cells, ssr["coil_cells"]), "coil cell mismatch"

ratio = np.where(cs > 0, cf / cs, np.nan)
rmean = np.nanmean(ratio); rcv = np.nanstd(ratio) / rmean
rho, _ = spearmanr(cf, cs)
pk_full = int(np.argmax(cf)); pk_ssr = int(np.argmax(cs))

print(f"{'coil':>6} {'full':>12} {'ssr':>12} {'ratio':>8}")
for i, c in enumerate(cells):
    print(f"{int(c):6d} {cf[i]:12.4e} {cs[i]:12.4e} {ratio[i]:8.3f}")
print(f"\nper-coil ratio full/ssr: mean k={rmean:.3f}  CV={rcv*100:.1f}%")
print(f"Spearman(full, ssr) per-coil = {rho:+.3f}")
print(f"hottest coil  full=cell {int(cells[pk_full])}  ssr=cell {int(cells[pk_ssr])}  "
      f"{'MATCH' if pk_full == pk_ssr else 'MISMATCH'}")
print("\nVERDICT:")
if rcv < 0.05 and pk_full == pk_ssr:
    print("  PASS -- SSR reproduces the coil-dose pattern (CV<5%, same peak). Decoupling valid;")
    print(f"          proceed to high-resolution breeder optimization. Normalization k={rmean:.3f}.")
elif rcv < 0.12:
    print("  MARGINAL -- pattern mostly reproduced but CV is elevated; inspect the interior-surface")
    print("          backscatter double-count before trusting fine rankings.")
else:
    print("  FAIL -- SSR distorts the pattern (CV high). Refine the capture (kill plasma-void re-entry)")
    print("          before any optimization; the frozen source is not yet a faithful stand-in.")
