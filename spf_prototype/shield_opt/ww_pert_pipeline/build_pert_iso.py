#!/usr/bin/env python
"""PHASE 2 (isolated): build the delta=25 perturbed geometry ALONE in $WORKDIR.
No concurrent build shares this export_dir, so a wedge here => geometry pathology (not a race).
Usage: python build_pert_iso.py AMP OUTNAME
"""
import os, sys
import numpy as np
import parastell.parastell as ps
from thickness_field import ThicknessField

AMP     = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "perturbed"
WORKDIR = os.environ.get("WORKDIR", "/burg-archive/home/tjk2147/pstl_test/kiso")
export_dir = WORKDIR
BASE = "/burg-archive/home/tjk2147/pstl_test"

toroidal_angles = list(np.linspace(0.0, 90.0, 13))
poloidal_angles = list(np.linspace(0.0, 360.0, 19))
repeat = 3; nfp = 4; wall_s = 1.08
ones = np.ones((len(toroidal_angles), len(poloidal_angles)))
T_FW, T_BR0, T_BW, T_SH0, T_VV = 4.0, 50.0, 3.0, 40.0, 10.0

tf = ThicknessField(nfp, toroidal_angles, poloidal_angles,
                    t_breeder0=T_BR0, t_shield0=T_SH0, delta_max=35.0, t_breeder_min=15.0)
delta = tf.gaussian_bump(theta0_deg=300.0, phi0_deg=30.0, amp=AMP, width_deg=35.0)
t_shield = T_SH0 + delta; t_breeder = T_BR0 - delta
assert np.all(t_breeder >= 15.0 - 1e-9)
assert np.allclose(t_shield + t_breeder, T_SH0 + T_BR0)
print(f"AMP={AMP} OUTNAME={OUTNAME} WORKDIR={WORKDIR}", flush=True)
print(f"delta max={delta.max():.2f} mean={delta.mean():.3f}  shield[{t_shield.min():.1f},{t_shield.max():.1f}] "
      f"breeder[{t_breeder.min():.1f},{t_breeder.max():.1f}]", flush=True)
it, ip = np.unravel_index(np.argmax(delta), delta.shape)
print(f"peak toroidal={toroidal_angles[it]:.1f} poloidal={poloidal_angles[ip]:.1f}", flush=True)

radial_build_dict = {
    "first_wall":    {"thickness_matrix": ones * T_FW},
    "breeder":       {"thickness_matrix": t_breeder},
    "back_wall":     {"thickness_matrix": ones * T_BW},
    "shield":        {"thickness_matrix": t_shield},
    "vacuum_vessel": {"thickness_matrix": ones * T_VV, "mat_tag": "vac_vessel"},
}
stellarator = ps.Stellarator(f"{BASE}/wout_qh.nc")
print("=== construct_invessel_build...", flush=True)
stellarator.construct_invessel_build(toroidal_angles, poloidal_angles, wall_s, radial_build_dict, repeat=repeat)
print("=== export invessel STEP...", flush=True)
stellarator.export_invessel_build_step(export_dir=export_dir)
print("=== magnets...", flush=True)
stellarator.construct_magnets_from_filaments(f"{BASE}/coils_qh", 40.0, 50.0, 360.0, sample_mod=1)
stellarator.export_magnets_step(export_dir=export_dir)
print("=== n coils:", len(stellarator.magnet_set.all_coil_solids), flush=True)
print("=== build_cad_to_dagmc_model...", flush=True)
stellarator.build_cad_to_dagmc_model()
stellarator.export_cad_to_dagmc(filename=f"dagmc_{OUTNAME}", export_dir=export_dir)
np.savez(f"{export_dir}/delta_{OUTNAME}.npz", delta=delta, t_shield=t_shield, t_breeder=t_breeder,
         toroidal_angles=np.array(toroidal_angles), poloidal_angles=np.array(poloidal_angles),
         amp=AMP, theta0=300.0, phi0=30.0, width=35.0)
print(f"DONE build_{OUTNAME}", flush=True)
