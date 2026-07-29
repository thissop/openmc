#!/usr/bin/env python
"""Tier-2 finite-difference-worth build: the corrected 8-layer QH build with extra
shield DELTA_CM concentrated in a SINGLE (theta,phi) patch (breeder traded locally,
envelope held), everything else baseline. Byte-for-byte identical to step1/build_step1.py
except the delta field is a single-patch indicator instead of the placed/uniform field,
so coil_run's material set matches exactly.

The delta is applied on the per-period (13 toroidal x 19 poloidal) grid, so it repeats in
all nfp field periods (the shield design is stellarator-symmetric). The peak-coil (cell 20)
response is dominated by the image in coil-20's own period; the finite-difference worth
W_FD = -(R_patch - R_baseline)/Delta_mass is the worth of that period-patch -- exactly the
per-period design control variable. Choose patches from patch_worth_tier1.csv (high/low A,
high/low W corners).

Usage: python build_patch.py <I_TOR> <J_POL> [N_TOR N_POL DELTA_CM]
       -> dagmc_qh_patch_t<I>p<J>.h5m   (baseline: build_step1.py with a zero delta)
"""
import sys

import numpy as np
import parastell.parastell as ps

I_TOR = int(sys.argv[1])
J_POL = int(sys.argv[2])
N_TOR = int(sys.argv[3]) if len(sys.argv) > 3 else 4
N_POL = int(sys.argv[4]) if len(sys.argv) > 4 else 4
DELTA_CM = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0   # extra shield in the patch
BE_CM = 2.0
OUTNAME = f"qh_patch_t{I_TOR}p{J_POL}"

# export_dir: 6th arg (per-task RUNDIR) so PARALLEL builds don't clobber each other's
# intermediate STEP/gmsh files in a shared dir. Default = the shared corrected dir (serial).
export_dir = sys.argv[6] if len(sys.argv) > 6 else "/burg-archive/home/tjk2147/pstl_test/corrected"

# ---- grid (MUST match build_step1.py) ----
toroidal_angles = list(np.linspace(0.0, 90.0, 13))    # one period, nfp=4
poloidal_angles = list(np.linspace(0.0, 360.0, 19))
repeat = 3
nfp = 4
wall_s = 1.08
ones = np.ones((len(toroidal_angles), len(poloidal_angles)))

T_FW, T_MULT, T_BR0, T_BW, T_SH0, T_GAP, T_VV, T_TS = \
    3.2, BE_CM, 50.0, 4.0, 40.0, 2.0, 25.0, 3.0

# ---- single-patch delta: DELTA_CM inside patch (I_TOR,J_POL) of the N_TOR x N_POL
#      coarse grid over one period (toroidal 0..90) x poloidal (0..360), else 0. ----
tor_c = np.asarray(toroidal_angles)                    # grid-cell centers (deg)
pol_c = np.asarray(poloidal_angles)
tor_edges = np.linspace(0.0, 90.0, N_TOR + 1)
pol_edges = np.linspace(0.0, 360.0, N_POL + 1)
i_of = np.clip(np.digitize(tor_c, tor_edges) - 1, 0, N_TOR - 1)   # which patch each grid col
j_of = np.clip(np.digitize(pol_c, pol_edges) - 1, 0, N_POL - 1)   # which patch each grid row
delta = np.zeros_like(ones)
sel_i = np.where(i_of == I_TOR)[0]
sel_j = np.where(j_of == J_POL)[0]
delta[np.ix_(sel_i, sel_j)] = DELTA_CM
assert sel_i.size and sel_j.size, f"empty patch ({I_TOR},{J_POL}) for {N_TOR}x{N_POL}"

t_shield = T_SH0 + delta
t_breeder = T_BR0 - delta
assert np.all(t_breeder >= 15.0 - 1e-9), f"breeder below floor: {t_breeder.min():.3f}"
assert np.allclose(t_shield + t_breeder, T_SH0 + T_BR0), "envelope not conserved"

# added-mass proxy for W_FD normalization: delta-volume-weighted (grid cells * DELTA_CM).
# Exact mass needs cell areas; the per-patch added shield VOLUME is proportional to
# (n grid cells in patch) * DELTA_CM at fixed poloidal/toroidal cell size -> report it.
added_cells = int(sel_i.size * sel_j.size)
print(f"VARIANT=patch({I_TOR},{J_POL}) of {N_TOR}x{N_POL}  OUTNAME={OUTNAME}", flush=True)
print(f"delta=DELTA_CM={DELTA_CM} on {added_cells} grid cells "
      f"(tor cols {sel_i.tolist()}, pol rows {sel_j.tolist()})", flush=True)
print(f"shield[{t_shield.min():.1f},{t_shield.max():.1f}] "
      f"breeder[{t_breeder.min():.1f},{t_breeder.max():.1f}]", flush=True)

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
         delta=delta, t_shield=t_shield, t_breeder=t_breeder,
         i_tor=I_TOR, j_pol=J_POL, n_tor=N_TOR, n_pol=N_POL, delta_cm=DELTA_CM,
         added_cells=added_cells,
         toroidal_angles=tor_c, poloidal_angles=pol_c)
print(f"DONE build_{OUTNAME}", flush=True)
