#!/usr/bin/env python
"""Regenerate the Step-1 placed/uniform shield-thickening delta field NATIVELY on the
CORRECTED 8-layer build grid, from the adjoint coil-20 (QH adjoint-peak) contributon.

Why regenerate (not reuse ~/qh_step1_deltas.npz):
  - The old (9,9) field was designed for the SIMPLIFIED 5-layer build (breeder0=60, so peaks
    to 45 cm). The literature-anchored corrected build has breeder0=50, delta_max=35,
    breeder_min=15 -> a 45 cm swing drives breeder to 5 cm and TRIPS build_corrected's
    assertion. So the old field is INFEASIBLE on the validated radial build.
  - placement_priority() regenerates the field on ANY (theta,phi) grid straight from the
    adjoint contributon, so we build it native to the corrected (13,19) grid: no (9,9)->(13,19)
    interpolation artifact.

Consistency choice: emissivity="uniform" for the plasma source S(r), because coil_run_v2's
forward StellaratorSource samples births ~|sqrt(g)| (uniform emissivity) -- the same source the
reciprocity study matched. The placed field is thus optimized for the SAME source the dose run
transports. (Bosch-Hale-weighted source + matching placement is a documented follow-up.)

Design: delta_placed = DELTA_MAX * P(theta,phi), P in [0,1] from the adjoint priority, so the
single hottest cell uses the full breeder->shield swing (breeder 50->15 exactly at the TBR floor).
delta_uniform = mean(delta_placed) constant -> EQUAL material budget (area-integral of delta).
The experiment: is concentrating that budget at the adjoint hotspot (placed) better than
spreading it evenly (uniform), in FULL TRANSPORT coil dose?
"""
import os, sys
import numpy as np

HOME = os.path.expanduser("~")
sys.path.insert(0, os.path.join(HOME, "src/GitHub/openmc-spf/spf_prototype/shield_opt"))
from adjoint_placement import contributon, placement_priority

MAP = os.path.join(HOME, "adjoint_test/recip_w35f_coil20_flat/adjoint_importance_flat_P0.npz")
FLM = os.path.join(HOME, "qh_freestream_fluxmap.npz")
OUT = os.path.join(HOME, "pstl_test/corrected/qh_step1_deltas_corr.npz")

# grid EXACTLY as build_corrected.py
TOR = list(np.linspace(0.0, 90.0, 13))    # 7.5 deg, one nfp=4 period
POL = list(np.linspace(0.0, 360.0, 19))   # 20 deg
# corrected-build constraints
T_BR0, T_SH0, DELTA_MAX, T_BR_MIN = 50.0, 40.0, 35.0, 15.0

print("=== Step-1 delta regeneration (corrected 8-layer grid) ===", flush=True)
print(f"MAP  = {MAP}", flush=True)
print(f"FLM  = {FLM}", flush=True)

field = contributon(MAP, FLM, scale=100.0, emissivity="uniform")  # match forward StellaratorSource
res = placement_priority(field, TOR, POL)                          # width_deg=15 default splat
P = res["priority"]                                                # (13,19) in [0,1]
assert P.shape == (13, 19), P.shape
assert abs(P.max() - 1.0) < 1e-9 and P.min() >= 0.0

delta_placed = DELTA_MAX * P
delta_uniform = np.full_like(delta_placed, float(delta_placed.mean()))

# feasibility on the corrected build
br_placed = T_BR0 - delta_placed
assert np.all(br_placed >= T_BR_MIN - 1e-9), f"breeder below floor: min={br_placed.min():.2f}"
assert np.all(delta_placed <= DELTA_MAX + 1e-9)
assert np.allclose(delta_placed.sum(), delta_uniform.sum())   # equal budget

# concentration diagnostics
flat = np.sort(delta_placed.ravel())[::-1]
n = flat.size
top25 = flat[:max(1, n // 4)].sum() / flat.sum()
print(f"placement target: phi0={res['phi0_deg']:.1f} deg  theta0={res['theta0_deg']:.1f} deg  "
      f"R0={res['R0']:.1f} cm", flush=True)
print(f"delta_placed : min {delta_placed.min():.3f}  max {delta_placed.max():.3f}  "
      f"mean {delta_placed.mean():.4f}  sum {delta_placed.sum():.3f} cm", flush=True)
print(f"delta_uniform: const {delta_uniform.flat[0]:.4f} cm (equal budget)", flush=True)
print(f"breeder placed: min {br_placed.min():.2f} cm (floor {T_BR_MIN})  "
      f"shield placed max {(T_SH0+delta_placed).max():.2f} cm", flush=True)
print(f"CONCENTRATION: top-25% cells hold {100*top25:.1f}% of placed shield material", flush=True)

np.savez(OUT, delta_placed=delta_placed, delta_uniform=delta_uniform,
         priority=P, phi0_deg=res["phi0_deg"], theta0_deg=res["theta0_deg"], R0=res["R0"],
         toroidal_angles=np.array(TOR), poloidal_angles=np.array(POL),
         t_breeder0=T_BR0, t_shield0=T_SH0, delta_max=DELTA_MAX, t_breeder_min=T_BR_MIN,
         emissivity="uniform", coil="cell20_adjoint_peak", build="corrected_w35f")
print(f"WROTE {OUT}", flush=True)
