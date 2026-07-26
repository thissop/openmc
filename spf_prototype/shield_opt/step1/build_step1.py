#!/usr/bin/env python
"""Step-1 kill-shot builds on the CORRECTED literature-anchored 8-layer radial build.

Identical to build_corrected.py EXCEPT the breeder<->shield delta comes from the regenerated,
feasibility-checked Step-1 field (qh_step1_deltas_corr.npz) instead of a gaussian bump:
  VARIANT=placed  -> delta = delta_placed  (adjoint-concentrated at coil-20 hotspot)
  VARIANT=uniform -> delta = delta_uniform (equal material budget, spread evenly)
Everything else (8 layers, breeder0=50, shield0=40, material tags, grid, magnets) is byte-for-byte
the corrected build so coil_run_v2.py's material set matches exactly (this is the fix for the
job-9191387 5-layer/9-material mismatch).

Usage: python build_step1.py {placed|uniform} [BE_CM]   ->  dagmc_qh_step1_<VARIANT>.h5m
"""
import sys
import numpy as np
import parastell.parastell as ps

VARIANT = sys.argv[1]
assert VARIANT in ("placed", "uniform"), VARIANT
BE_CM = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
OUTNAME = f"qh_step1_{VARIANT}"

export_dir = "/burg-archive/home/tjk2147/pstl_test/corrected"

# ---- Fine angular grid (MUST match build_corrected.py AND the regenerated delta grid) ----
toroidal_angles = list(np.linspace(0.0, 90.0, 13))    # one period, 7.5 deg
poloidal_angles = list(np.linspace(0.0, 360.0, 19))   # 20 deg
repeat = 3
nfp = 4
wall_s = 1.08
ones = np.ones((len(toroidal_angles), len(poloidal_angles)))

# ---- corrected baseline thicknesses (cm), same as build_corrected.py ----
T_FW, T_MULT, T_BR0, T_BW, T_SH0, T_GAP, T_VV, T_TS = \
    3.2, BE_CM, 50.0, 4.0, 40.0, 2.0, 25.0, 3.0

# ---- load the regenerated, feasibility-checked Step-1 delta ----
dd = np.load(f"{export_dir}/qh_step1_deltas_corr.npz")
delta = dd["delta_placed"] if VARIANT == "placed" else dd["delta_uniform"]
assert delta.shape == ones.shape, (delta.shape, ones.shape)
# re-assert the corrected-build invariants HERE (hard gate; the gen script also checked)
t_shield = T_SH0 + delta
t_breeder = T_BR0 - delta
assert np.all(t_breeder >= 15.0 - 1e-9), f"breeder below floor: {t_breeder.min():.3f}"
assert np.allclose(t_shield + t_breeder, T_SH0 + T_BR0), "envelope not conserved"

total = T_FW + T_MULT + T_BR0 + T_BW + T_SH0 + T_GAP + T_VV + T_TS
print(f"VARIANT={VARIANT} OUTNAME={OUTNAME} BE_CM={BE_CM}", flush=True)
print(f"total plasma->coil standoff (uniform envelope) = {total:.1f} cm", flush=True)
print(f"delta: max={delta.max():.2f} mean={delta.mean():.3f} cm  budget(sum)={delta.sum():.3f}  "
      f"shield[{t_shield.min():.1f},{t_shield.max():.1f}] breeder[{t_breeder.min():.1f},{t_breeder.max():.1f}]",
      flush=True)

radial_build_dict = {
    "first_wall":    {"thickness_matrix": ones * T_FW},
    "multiplier":    {"thickness_matrix": ones * T_MULT},
    "breeder":       {"thickness_matrix": t_breeder},
    "back_wall":     {"thickness_matrix": ones * T_BW},
    "shield":        {"thickness_matrix": t_shield},
    "gap":           {"thickness_matrix": ones * T_GAP},
    "vacuum_vessel": {"thickness_matrix": ones * T_VV, "mat_tag": "vac_vessel"},
    "thermal_shield":{"thickness_matrix": ones * T_TS},
}

stellarator = ps.Stellarator("/burg-archive/home/tjk2147/pstl_test/wout_qh.nc")
print("=== Constructing in-vessel build...", flush=True)
stellarator.construct_invessel_build(
    toroidal_angles, poloidal_angles, wall_s, radial_build_dict, repeat=repeat)
print("=== Exporting in-vessel STEP...", flush=True)
stellarator.export_invessel_build_step(export_dir=export_dir)

print("=== Constructing magnets from filaments...", flush=True)
stellarator.construct_magnets_from_filaments(
    "/burg-archive/home/tjk2147/pstl_test/coils_qh", 40.0, 50.0, 360.0, sample_mod=1)
print("=== Exporting magnets STEP...", flush=True)
stellarator.export_magnets_step(export_dir=export_dir)
print("=== n coil solids:", len(stellarator.magnet_set.all_coil_solids), flush=True)

print("=== Building cad_to_dagmc model...", flush=True)
stellarator.build_cad_to_dagmc_model()
print("=== Material tags:", stellarator._material_tags, flush=True)
print(f"=== Exporting DAGMC H5M (dagmc_{OUTNAME})...", flush=True)
stellarator.export_cad_to_dagmc(filename=f"dagmc_{OUTNAME}", export_dir=export_dir)
np.savez(f"{export_dir}/delta_{OUTNAME}.npz",
         delta=delta, t_shield=t_shield, t_breeder=t_breeder, variant=VARIANT, be_cm=BE_CM,
         toroidal_angles=np.array(toroidal_angles), poloidal_angles=np.array(poloidal_angles))
print(f"DONE build_{OUTNAME}", flush=True)
