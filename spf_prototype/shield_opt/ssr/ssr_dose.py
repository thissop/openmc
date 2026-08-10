#!/usr/bin/env python
"""SSR shielding run: read a frozen first-wall surface source and transport it through a given radial
build to the magnets, tallying the coil dose. Skips the expensive plasma+source sampling.

For Gate 2, run this on the SAME (uniform) build the master run used and compare per-coil dose:
if the per-coil RATIO full/ssr is constant across coils, SSR reproduces the spatial pattern and the
decoupling is valid (the constant is just the source-normalization factor k = crossings/birth).

Usage: ssr_dose.py BATCHES PARTICLES DAGMC OUT_TAG LI6 CELLS_NPZ SURF_SRC
Env: RUNDIR, OMP_NUM_THREADS.
"""
import os, sys, shutil, time
import numpy as np
import openmc

BASE = "/burg-archive/home/tjk2147/pstl_test"
WORK = f"{BASE}/corrected"
RUNDIR = os.environ.get("RUNDIR", WORK)
os.makedirs(RUNDIR, exist_ok=True); os.chdir(RUNDIR)

BATCHES   = int(sys.argv[1]); PARTICLES = int(sys.argv[2])
DAGMC     = sys.argv[3];      OUT_TAG   = sys.argv[4]
LI6       = float(sys.argv[5]); CELLS_NPZ = sys.argv[6]
SURF_SRC  = sys.argv[7]

for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml"):
    try: os.remove(f)
    except FileNotFoundError: pass
assert os.path.exists(SURF_SRC), f"surface source not found: {SURF_SRC}"

# build_materials identical to ssw_master (kept in sync deliberately -- same physics)
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
    return openmc.Materials([fw, mult, br, bw, sh, gap, vv, ts, mg, vac]), br


mats, breeder_mat = build_materials(LI6)
dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
ext = float(np.max(np.abs(np.concatenate([np.array(bb.lower_left), np.array(bb.upper_right)]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
geom = openmc.Geometry([openmc.Cell(region=-sph, fill=dag)])

COIL_CELLS = sorted(int(x) for x in np.load(CELLS_NPZ)["magnet_cells"])
fast = openmc.EnergyFilter([0.1e6, 25.0e6]); cf = openmc.CellFilter(COIL_CELLS)
t_flux = openmc.Tally(name="coil_fast_flux"); t_flux.filters = [cf, fast]; t_flux.scores = ["flux"]
t_heat = openmc.Tally(name="coil_heating"); t_heat.filters = [cf]; t_heat.scores = ["heating-local"]
tallies = openmc.Tallies([t_flux, t_heat])

s = openmc.Settings()
s.run_mode = "fixed source"
s.surf_source_read = {"path": SURF_SRC}          # <-- the frozen first-wall source becomes the source
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False; s.max_lost_particles = 100; s.rel_max_lost_particles = 1e-4
s.output = {"summary": True, "tallies": False}

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tallies)
print(f"=== SSR TAG={OUT_TAG} DAGMC={os.path.basename(DAGMC)} src={os.path.basename(SURF_SRC)} "
      f"batches={BATCHES} particles={PARTICLES} cwd={RUNDIR} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)

with openmc.StatePoint(sp_path) as sp:
    flux = sp.get_tally(name="coil_fast_flux").mean.ravel(); ferr = sp.get_tally(name="coil_fast_flux").std_dev.ravel()
    heat = sp.get_tally(name="coil_heating").mean.ravel()

np.savez(f"{WORK}/coil_{OUT_TAG}.npz", tag=OUT_TAG, kind="ssr", coil_cells=np.array(COIL_CELLS),
         coil_fast_flux=flux, coil_fast_flux_err=ferr, coil_heating=heat,
         batches=BATCHES, particles=PARTICLES, n_source=BATCHES*PARTICLES, surf_src=SURF_SRC)
print(f"COIL peak/mean={flux.max()/flux.mean():.3f} hottest_idx={int(np.argmax(flux))}", flush=True)
print(f"WROTE {WORK}/coil_{OUT_TAG}.npz", flush=True)
print("DONE_" + OUT_TAG, flush=True)
