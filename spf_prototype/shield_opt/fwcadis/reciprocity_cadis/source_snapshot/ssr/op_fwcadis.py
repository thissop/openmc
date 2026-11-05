#!/usr/bin/env python
"""FW-CADIS weight-window generation for the one-period magnet-shielding problem, using OpenMC's
native random-ray multigroup solver (no external Denovo/ADVANTG). Adjoint-sourced from the coil
tally so the windows drive neutrons THROUGH the shield TO the coils -- the deep-penetration gold
standard, and far more efficient than MAGIC (which over-splits).

Workflow (OpenMC automates forward+adjoint random-ray, then writes weight_windows.h5):
  convert_to_multigroup -> convert_to_random_ray -> WeightWindowGenerator(method='fw_cadis',
  targets=[coil_tally]) -> run.

Usage: op_fwcadis.py PARTICLES [LI6]     (PARTICLES = random-ray rays/batch)
Env: RUNDIR. Writes ww_op_fwcadis.h5.
NOTE: uses a VACUUM-sphere boundary (not periodic) for the random-ray solve -- the importance map
through the shield is a local/radial quantity, insensitive to the toroidal BC, so the windows are
valid to apply to the periodic-BC continuous-energy SSR run.
"""
import os, sys, shutil
import numpy as np
import openmc

BASE = "/burg-archive/home/tjk2147/pstl_test"
WORK = f"{BASE}/corrected"
RUNDIR = os.environ.get("RUNDIR", WORK); os.makedirs(RUNDIR, exist_ok=True); os.chdir(RUNDIR)
PARTICLES = int(sys.argv[1]); LI6 = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
DAGMC = f"{WORK}/dagmc_qh_op_uniform.h5m"
MAGCELLS = [12, 13, 14, 15, 16]

for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml"):
    try: os.remove(f)
    except FileNotFoundError: pass


def build_materials(li6):
    w = openmc.Material(); w.add_element("W", 1.0); w.set_density("g/cm3", 19.3)
    rafm = openmc.Material()
    for el, wo in [("Fe", .885), ("Cr", .09), ("W", .011), ("Mn", .006), ("V", .002), ("C", .006)]:
        rafm.add_element(el, wo, "wo")
    rafm.set_density("g/cm3", 7.9)
    fw = openmc.Material.mix_materials([w, rafm], [0.2/3.2, 3.0/3.2], "vo"); fw.name = "first_wall"
    mult = openmc.Material(name="multiplier"); mult.add_element("Be", 1.0); mult.set_density("g/cm3", 1.85)
    br = openmc.Material(name="breeder")
    br.add_element("Li", 2.0, enrichment=li6, enrichment_target="Li6", enrichment_type="ao")
    br.add_element("Be", 1.0); br.add_element("F", 4.0); br.set_density("g/cm3", 1.94)
    bw = openmc.Material(name="back_wall")
    for el, wo in [("Fe", .885), ("Cr", .09), ("W", .011), ("Mn", .006), ("V", .002), ("C", .006)]:
        bw.add_element(el, wo, "wo")
    bw.set_density("g/cm3", 7.9)
    sh = openmc.Material(name="shield"); sh.add_element("W", 1.0); sh.add_element("C", 1.0); sh.set_density("g/cm3", 15.6)
    gap = openmc.Material(name="gap"); gap.add_nuclide("H1", 1.0); gap.set_density("g/cm3", 1e-8)
    vv = openmc.Material(name="vac_vessel")
    for el, wo in [("Fe", .62), ("Cr", .16), ("Ni", .11), ("Mo", .02), ("Mn", .015), ("B", .03), ("O", .03), ("H", .005)]:
        vv.add_element(el, wo, "wo")
    vv.set_density("g/cm3", 7.9)
    ts = openmc.Material(name="thermal_shield")
    ts.add_element("Fe", .89, "wo"); ts.add_element("Cr", .09, "wo"); ts.add_element("Mn", .02, "wo"); ts.set_density("g/cm3", 7.9)
    mg = openmc.Material(name="magnets")
    for el, wo in [("Cu", .58), ("Fe", .30), ("Cr", .04), ("Ni", .03), ("Y", .02), ("Ba", .02), ("O", .01)]:
        mg.add_element(el, wo, "wo")
    mg.set_density("g/cm3", 8.5)
    vac = openmc.Material(name="Vacuum"); vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-10)
    return openmc.Materials([fw, mult, br, bw, sh, gap, vv, ts, mg, vac])


mats = build_materials(LI6)
dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
ext = float(np.max(np.abs(np.concatenate([np.array(bb.lower_left), np.array(bb.upper_right)]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
root = openmc.Cell(region=-sph, fill=dag)
geom = openmc.Geometry([root])

# forward source for the random-ray solve: the plasma (14.06 MeV). The FW-CADIS adjoint source is the
# coil tally, so the windows are built to reach the coils regardless of the CE run's frozen source.
src = openmc.StellaratorSource(fluxmap=f"{BASE}/qh_fluxmap", polarization=(1/3, 1/3, 1/3),
                               field_model="fluxmap", energy=openmc.stats.Discrete([14.06e6], [1.0]),
                               strength=1.0)

s = openmc.Settings()
s.run_mode = "fixed source"; s.source = [src]
s.particles = PARTICLES
s.batches = 20; s.inactive = 5

# target tally = coil fast flux (drives the FW-CADIS adjoint source at the coils)
fast = openmc.EnergyFilter([0.1e6, 25.0e6]); cf = openmc.CellFilter(MAGCELLS)
t_coil = openmc.Tally(name="coil_fast_flux"); t_coil.filters = [cf, fast]; t_coil.scores = ["flux"]
tals = openmc.Tallies([t_coil])

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tals)

# --- convert to multigroup + random ray ---
print("=== convert_to_multigroup ===", flush=True)
model.convert_to_multigroup(method="infinite_medium")
print("=== convert_to_random_ray ===", flush=True)
model.convert_to_random_ray()
ml = 2.0*float(np.max(np.array(bb.upper_right)-np.array(bb.lower_left)))
model.settings.random_ray['distance_inactive'] = 0.2*ml
model.settings.random_ray['distance_active'] = ml

# source-region mesh (also the WW mesh) over the wedge bounding box
n = 10
mesh = openmc.RegularMesh(); mesh.dimension = (n, n, n)
mesh.lower_left = [float(v) for v in bb.lower_left]
mesh.upper_right = [float(v) for v in bb.upper_right]
model.settings.random_ray['source_shape'] = 'flat'

# FW-CADIS weight window generator targeting the coil tally
wwg = openmc.WeightWindowGenerator(method='fw_cadis', mesh=mesh)
try:
    wwg.targets = [t_coil]
except Exception as e:
    print("targets set failed:", e, flush=True)
model.settings.weight_window_generators = wwg

print(f"=== FW-CADIS random-ray solve rays={PARTICLES} mesh={n}^3 cwd={RUNDIR} ===", flush=True)
model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))

for cand in ("weight_windows.h5", os.path.join(RUNDIR, "weight_windows.h5")):
    if os.path.exists(cand):
        shutil.move(cand, f"{WORK}/ww_op_fwcadis.h5"); print(f"WROTE {WORK}/ww_op_fwcadis.h5", flush=True); break
else:
    print("WARNING: weight_windows.h5 not found", flush=True)
print("DONE_op_fwcadis", flush=True)
