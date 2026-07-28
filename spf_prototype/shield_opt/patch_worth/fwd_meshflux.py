#!/usr/bin/env python
"""Tier-0 prerequisite: forward flux phi(x) on the SAME mesh as the adjoint psi_dagger,
so patch_worth.py can form the contributon C = phi*psi_dagger in the shield band.

REWRITTEN to match the REAL cluster recipe (coil_run_v3.py), NOT the CompiledSource
assumption of the first draft:
  * source   = openmc.StellaratorSource(fluxmap=..., polarization=unpol, field_model="fluxmap")
  * materials = built IN-CODE (build_materials, copied verbatim from coil_run_v3.py) so the
               forward flux is physically identical to the per-coil dose runs
  * geometry = DAGMC baseline (uniform) wrapped in a vacuum sphere (same as coil_run_v3)
  * mesh     = RegularMesh read straight from the adjoint npz (lower_left/upper_right/dimension)
               -> forward and adjoint fields are co-located voxel-for-voxel
  * WW       = optional FW-CADIS weight windows (fwcadis/weight_windows.h5)

Usage (positional, coil_run_v3 style):
  python fwd_meshflux.py BATCHES PARTICLES DAGMC ADJOINT_NPZ OUT_NPZ LI6 [WW_FILE]

Writes OUT_NPZ with key 'flux' (shape = adjoint 'dimension', C-order to match psi.reshape).
Env: conda activate spf-stellarator; PYTHONPATH=$SRC; LD_PRELOAD=$SRC/build/lib/libopenmc.so.
"""
import os
import sys
import time

import numpy as np
import openmc

BASE = "/burg-archive/home/tjk2147/pstl_test"
os.chdir(os.path.join(BASE, "corrected"))

BATCHES   = int(sys.argv[1])
PARTICLES = int(sys.argv[2])
DAGMC     = sys.argv[3]
ADJOINT   = sys.argv[4]
OUT_NPZ   = sys.argv[5]
LI6       = float(sys.argv[6]) if len(sys.argv) > 6 else 60.0
WW_FILE   = sys.argv[7] if len(sys.argv) > 7 else None

# clear stale XML so OpenMC loads THIS geometry, not a leftover model.xml
for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml"):
    try: os.remove(f)
    except FileNotFoundError: pass


def build_materials(li6):
    """VERBATIM copy of coil_run_v3.build_materials -> forward flux matches the dose runs."""
    w = openmc.Material(); w.add_element("W", 1.0); w.set_density("g/cm3", 19.3)
    rafm = openmc.Material()
    rafm.add_element("Fe", 0.885, "wo"); rafm.add_element("Cr", 0.09, "wo")
    rafm.add_element("W", 0.011, "wo"); rafm.add_element("Mn", 0.006, "wo")
    rafm.add_element("V", 0.002, "wo"); rafm.add_element("C", 0.006, "wo")
    rafm.set_density("g/cm3", 7.9)
    fW, fR = 0.2 / 3.2, 3.0 / 3.2
    fw = openmc.Material.mix_materials([w, rafm], [fW, fR], "vo"); fw.name = "first_wall"

    mult = openmc.Material(name="multiplier")
    mult.add_element("Be", 1.0); mult.set_density("g/cm3", 1.85)

    br = openmc.Material(name="breeder")
    br.add_element("Li", 2.0, enrichment=li6, enrichment_target="Li6", enrichment_type="ao")
    br.add_element("Be", 1.0); br.add_element("F", 4.0); br.set_density("g/cm3", 1.94)

    bw = openmc.Material(name="back_wall")
    bw.add_element("Fe", 0.885, "wo"); bw.add_element("Cr", 0.09, "wo")
    bw.add_element("W", 0.011, "wo"); bw.add_element("Mn", 0.006, "wo")
    bw.add_element("V", 0.002, "wo"); bw.add_element("C", 0.006, "wo")
    bw.set_density("g/cm3", 7.9)

    sh = openmc.Material(name="shield")
    sh.add_element("W", 1.0); sh.add_element("C", 1.0); sh.set_density("g/cm3", 15.6)

    gap = openmc.Material(name="gap")
    gap.add_nuclide("H1", 1.0); gap.set_density("g/cm3", 1e-8)

    vv = openmc.Material(name="vac_vessel")
    vv.add_element("Fe", 0.62, "wo"); vv.add_element("Cr", 0.16, "wo")
    vv.add_element("Ni", 0.11, "wo"); vv.add_element("Mo", 0.02, "wo")
    vv.add_element("Mn", 0.015, "wo"); vv.add_element("B", 0.03, "wo")
    vv.add_element("O", 0.03, "wo"); vv.add_element("H", 0.005, "wo")
    vv.set_density("g/cm3", 7.9)

    ts = openmc.Material(name="thermal_shield")
    ts.add_element("Fe", 0.89, "wo"); ts.add_element("Cr", 0.09, "wo")
    ts.add_element("Mn", 0.02, "wo"); ts.set_density("g/cm3", 7.9)

    mg = openmc.Material(name="magnets")
    mg.add_element("Cu", 0.58, "wo"); mg.add_element("Fe", 0.30, "wo")
    mg.add_element("Cr", 0.04, "wo"); mg.add_element("Ni", 0.03, "wo")
    mg.add_element("Y", 0.02, "wo"); mg.add_element("Ba", 0.02, "wo")
    mg.add_element("O", 0.01, "wo"); mg.set_density("g/cm3", 8.5)

    vac = openmc.Material(name="Vacuum")
    vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-10)

    return openmc.Materials([fw, mult, br, bw, sh, gap, vv, ts, mg, vac])


# --- mesh from adjoint npz (co-located with psi_dagger) ---
da = np.load(ADJOINT)
ll = np.asarray(da["lower_left"], float)
ur = np.asarray(da["upper_right"], float)
dim = [int(v) for v in np.asarray(da["dimension"])]

mats = build_materials(LI6)

dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
lo, hi = np.array(bb.lower_left), np.array(bb.upper_right)
ext = float(np.max(np.abs(np.concatenate([lo, hi]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
root = openmc.Cell(region=-sph, fill=dag)
geom = openmc.Geometry([root])

# --- source: SAME as coil_run_v3 (unpol, fluxmap field model, 14.06 MeV) ---
src = openmc.StellaratorSource(
    fluxmap=f"{BASE}/qh_fluxmap", polarization=(1/3, 1/3, 1/3), field_model="fluxmap",
    energy=openmc.stats.Discrete([14.06e6], [1.0]), strength=1.0)

s = openmc.Settings()
s.run_mode = "fixed source"; s.source = [src]
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False
s.max_lost_particles = 100
s.rel_max_lost_particles = 1e-4
s.output = {"summary": True, "tallies": False}
if WW_FILE and os.path.exists(WW_FILE):
    s.weight_windows = openmc.hdf5_to_wws(WW_FILE)
    s.weight_windows_on = True
    print("Using weight windows from", WW_FILE, flush=True)

mesh = openmc.RegularMesh()
mesh.lower_left = ll.tolist()
mesh.upper_right = ur.tolist()
mesh.dimension = dim
t = openmc.Tally(name="fwd_meshflux")
t.filters = [openmc.MeshFilter(mesh)]
t.scores = ["flux"]

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=openmc.Tallies([t]))
print(f"=== fwd_meshflux DAGMC={os.path.basename(DAGMC)} mesh={dim} "
      f"batches={BATCHES} particles={PARTICLES} LI6={LI6} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)

with openmc.StatePoint(sp_path) as sp:
    tf = sp.get_tally(name="fwd_meshflux")
    flux = tf.mean.ravel()          # C-order, matches psi.reshape(dim) in patch_worth
    err = tf.std_dev.ravel()
np.savez(OUT_NPZ, flux=flux.reshape(dim), err=err.reshape(dim), dimension=np.array(dim))
rel = err[flux > 0] / flux[flux > 0]
print(f"WROTE {OUT_NPZ}  mesh={dim}  nonzero voxels={int((flux>0).sum())}  "
      f"median relerr={np.median(rel):.3f}  p90 relerr={np.percentile(rel,90):.3f}", flush=True)
print("DONE_fwd_meshflux", flush=True)
