#!/usr/bin/env python
"""PHASE 2: build a fixed-envelope shield-thickening / breeder-thinning perturbation
targeting the hot lower coils, OR a matched-grid uniform baseline (AMP=0).

Usage: python build_perturbed.py AMP OUTNAME
  AMP=0    -> uniform baseline on the fine grid (dagmc_<OUTNAME>.h5m)
  AMP>0    -> gaussian shield bump (peak=AMP cm) at (theta0,phi0); breeder thinned by same.

Reuses build_baseline.py structure: same wout, repeat=3, magnets from coils_qh.
Fixed envelope: shield=40+delta, breeder=50-delta, breeder_min=15, first/back_wall/vv fixed.
Bump centered on the hot lower-coil angular cluster: poloidal theta0=300, toroidal phi0=30
(relative within the 90-deg field period). nfp=4.
"""
import sys
import numpy as np
import parastell.parastell as ps
from thickness_field import ThicknessField

AMP     = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else ("perturbed" if AMP > 0 else "baseline_fine")

export_dir = "/burg-archive/home/tjk2147/pstl_test/ww_pert"

# ---- Fine angular grid (resolves the ~35-deg bump; uniform regions unchanged) ----
toroidal_angles = list(np.linspace(0.0, 90.0, 13))    # one period, 7.5 deg
poloidal_angles = list(np.linspace(0.0, 360.0, 19))   # 20 deg
repeat = 3
nfp = 4
wall_s = 1.08
ones = np.ones((len(toroidal_angles), len(poloidal_angles)))

# ---- baseline thicknesses (cm) ----
T_FW, T_BR0, T_BW, T_SH0, T_VV = 4.0, 50.0, 3.0, 40.0, 10.0

# ---- fixed-envelope breeder<->shield trade via smooth gaussian bump ----
tf = ThicknessField(nfp, toroidal_angles, poloidal_angles,
                    t_breeder0=T_BR0, t_shield0=T_SH0,
                    delta_max=35.0, t_breeder_min=15.0)
if AMP > 0:
    delta = tf.gaussian_bump(theta0_deg=300.0, phi0_deg=30.0, amp=AMP, width_deg=35.0)
else:
    delta = np.zeros_like(ones)
t_shield = T_SH0 + delta
t_breeder = T_BR0 - delta
assert np.all(t_breeder >= 15.0 - 1e-9)
assert np.allclose(t_shield + t_breeder, T_SH0 + T_BR0)
print(f"AMP={AMP} OUTNAME={OUTNAME}", flush=True)
print(f"delta: max={delta.max():.2f} mean={delta.mean():.3f} cm  "
      f"shield[{t_shield.min():.1f},{t_shield.max():.1f}] "
      f"breeder[{t_breeder.min():.1f},{t_breeder.max():.1f}]", flush=True)
# where is the peak?
it, ip = np.unravel_index(np.argmax(delta), delta.shape)
print(f"peak at toroidal={toroidal_angles[it]:.1f}deg poloidal={poloidal_angles[ip]:.1f}deg", flush=True)

radial_build_dict = {
    "first_wall":    {"thickness_matrix": ones * T_FW},
    "breeder":       {"thickness_matrix": t_breeder},
    "back_wall":     {"thickness_matrix": ones * T_BW},
    "shield":        {"thickness_matrix": t_shield},
    "vacuum_vessel": {"thickness_matrix": ones * T_VV, "mat_tag": "vac_vessel"},
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
         toroidal_angles=np.array(toroidal_angles), poloidal_angles=np.array(poloidal_angles),
         amp=AMP, theta0=300.0, phi0=30.0, width=35.0)
print(f"DONE build_{OUTNAME}", flush=True)
