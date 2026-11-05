#!/usr/bin/env python
"""One-field-period master: periodic-BC wedge + REAL symmetric source restricted to the wedge, with a
DEEP freeze at the fixed-breeder boundary (cellfrom=5 breeder_fixed -> cellto=6 breeder_var, so ONLY
the outward current at that surface is banked -- no inner backscatter, no non-convex double-count).
Also tallies the reference coil dose (cells 12-16) for the Gate-2 comparison.

Usage: op_master.py BATCHES PARTICLES [LI6]
"""
import os, sys, time, shutil
import numpy as np
import openmc

BASE = "/burg-archive/home/tjk2147/pstl_test"
WORK = f"{BASE}/corrected"
RUNDIR = os.environ.get("RUNDIR", WORK); os.makedirs(RUNDIR, exist_ok=True); os.chdir(RUNDIR)
BATCHES = int(sys.argv[1]); PARTICLES = int(sys.argv[2])
LI6 = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
DAGMC = f"{WORK}/dagmc_qh_op_uniform.h5m"
MAGCELLS = [12, 13, 14, 15, 16]
CF_FROM, CF_TO = 5, 6          # freeze at breeder_fixed -> breeder_var boundary (outward only)

for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml", "surface_source.h5"):
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

# rotational-periodic wedge: phi=0 face is the y=0 plane, phi=90 face is the x=0 plane.
py = openmc.YPlane(0.0, boundary_type="periodic")     # wedge on +y side
px = openmc.XPlane(0.0, boundary_type="periodic")     # wedge on +x side
py.periodic_surface = px
root = openmc.Cell(region=+py & +px & -sph, fill=dag)
geom = openmc.Geometry([root])

# REAL symmetric source, restricted to the wedge by domain rejection (only births inside the DAGMC
# wedge are kept; the other 3/4 of the full-torus fluxmap are resampled away).
src = openmc.StellaratorSource(
    fluxmap=f"{BASE}/qh_fluxmap", polarization=(1/3, 1/3, 1/3), field_model="fluxmap",
    energy=openmc.stats.Discrete([14.06e6], [1.0]), strength=1.0,
    constraints={"domains": [root], "rejection_strategy": "resample"})

s = openmc.Settings()
s.run_mode = "fixed source"; s.source = [src]
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False
s.max_lost_particles = 500; s.rel_max_lost_particles = 2e-3
s.output = {"summary": True, "tallies": False}
s.surf_source_write = {"cellfrom": CF_FROM, "max_particles": 80_000_000}  # bank all crossings leaving fixed zone; inner ones die on black-inside replay

fast = openmc.EnergyFilter([0.1e6, 25.0e6]); cf = openmc.CellFilter(MAGCELLS)
t_flux = openmc.Tally(name="coil_fast_flux"); t_flux.filters = [cf, fast]; t_flux.scores = ["flux"]
t_heat = openmc.Tally(name="coil_heating"); t_heat.filters = [cf]; t_heat.scores = ["heating-local"]
tals = openmc.Tallies([t_flux, t_heat])

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tals)
print(f"=== OP MASTER (1 period, periodic) freeze {CF_FROM}->{CF_TO} batches={BATCHES} "
      f"particles={PARTICLES} coils={MAGCELLS} cwd={RUNDIR} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)

with openmc.StatePoint(sp_path) as sp:
    fl = sp.get_tally(name="coil_fast_flux"); flux = fl.mean.ravel(); ferr = fl.std_dev.ravel()
    heat = sp.get_tally(name="coil_heating").mean.ravel()

src_out = f"{WORK}/surf_op_master.h5"
for cand in ("surface_source.h5", os.path.join(RUNDIR, "surface_source.h5")):
    if os.path.exists(cand):
        shutil.move(cand, src_out); break
else:
    print("WARNING surface_source.h5 not found", flush=True); src_out = ""

np.savez(f"{WORK}/coil_op_master.npz", tag="op_master", coil_cells=np.array(MAGCELLS),
         coil_fast_flux=flux, coil_fast_flux_err=ferr, coil_heating=heat,
         batches=BATCHES, particles=PARTICLES, n_source=BATCHES*PARTICLES,
         surf_file=src_out, cellfrom=CF_FROM, cellto=CF_TO)
print("COIL flux:", np.array2string(flux, precision=4), flush=True)
print(f"WROTE {WORK}/coil_op_master.npz and {src_out}", flush=True)
print("DONE_op_master", flush=True)
