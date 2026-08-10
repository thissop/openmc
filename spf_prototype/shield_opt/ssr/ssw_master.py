#!/usr/bin/env python
"""Master surface-source-write run: full transport from the plasma source through the radial build,
freezing the INWARD first-wall neutron source (cellfrom = plasma void, cell 2) so the shield can be
re-optimized cheaply later without re-running the plasma+geometry sampling.

This is a full run (identical physics to coil_run_v3), so it ALSO produces the coil dose -- which is
the FULL-RUN REFERENCE for Gate 2 (compare against the SSR run reading this frozen source).

Usage: ssw_master.py BATCHES PARTICLES DAGMC OUT_TAG LI6 CELLS_NPZ MAXPART CELLFROM
Writes: coil_<TAG>.npz (full-run coil dose) and surf_<TAG>.h5 (the frozen surface source).
Env: RUNDIR (isolated cwd, avoids statepoint/model.xml races); OMP_NUM_THREADS.
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
MAXPART   = int(sys.argv[7]); CELLFROM  = int(sys.argv[8])

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
    return openmc.Materials([fw, mult, br, bw, sh, gap, vv, ts, mg, vac]), br


mats, breeder_mat = build_materials(LI6)
dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
ext = float(np.max(np.abs(np.concatenate([np.array(bb.lower_left), np.array(bb.upper_right)]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
geom = openmc.Geometry([openmc.Cell(region=-sph, fill=dag)])

COIL_CELLS = sorted(int(x) for x in np.load(CELLS_NPZ)["magnet_cells"])
assert COIL_CELLS, "no magnet cells"

fast = openmc.EnergyFilter([0.1e6, 25.0e6])
cf = openmc.CellFilter(COIL_CELLS)
t_flux = openmc.Tally(name="coil_fast_flux"); t_flux.filters = [cf, fast]; t_flux.scores = ["flux"]
t_heat = openmc.Tally(name="coil_heating"); t_heat.filters = [cf]; t_heat.scores = ["heating-local"]
t_fw = openmc.Tally(name="fw_fast_flux"); t_fw.filters = [openmc.MaterialFilter([mats[0]]), fast]; t_fw.scores = ["flux"]
t_tbr = openmc.Tally(name="tbr"); t_tbr.filters = [openmc.MaterialFilter([breeder_mat])]; t_tbr.scores = ["H3-production"]
tallies = openmc.Tallies([t_flux, t_heat, t_fw, t_tbr])

src = openmc.StellaratorSource(fluxmap=f"{BASE}/qh_fluxmap", polarization=(1/3, 1/3, 1/3),
                               field_model="fluxmap", energy=openmc.stats.Discrete([14.06e6], [1.0]), strength=1.0)
s = openmc.Settings()
s.run_mode = "fixed source"; s.source = [src]
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False; s.max_lost_particles = 100; s.rel_max_lost_particles = 1e-4
s.output = {"summary": True, "tallies": False}
s.surf_source_write = {"cellfrom": CELLFROM, "max_particles": MAXPART}   # FREEZE the inward FW source

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tallies)
print(f"=== SSW MASTER TAG={OUT_TAG} DAGMC={os.path.basename(DAGMC)} batches={BATCHES} "
      f"particles={PARTICLES} cellfrom={CELLFROM} maxpart={MAXPART} cwd={RUNDIR} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)

with openmc.StatePoint(sp_path) as sp:
    flux = sp.get_tally(name="coil_fast_flux").mean.ravel(); ferr = sp.get_tally(name="coil_fast_flux").std_dev.ravel()
    heat = sp.get_tally(name="coil_heating").mean.ravel()
    fwf = sp.get_tally(name="fw_fast_flux").mean.ravel()[0]
    tbr = sp.get_tally(name="tbr").mean.ravel()[0]

# move the frozen surface source next to WORK under a tagged name
src_out = f"{WORK}/surf_{OUT_TAG}.h5"
for cand in ("surface_source.h5", os.path.join(RUNDIR, "surface_source.h5")):
    if os.path.exists(cand):
        shutil.move(cand, src_out); break
else:
    print("WARNING: surface_source.h5 not found", flush=True); src_out = ""

np.savez(f"{WORK}/coil_{OUT_TAG}.npz", tag=OUT_TAG, kind="ssw_master", coil_cells=np.array(COIL_CELLS),
         coil_fast_flux=flux, coil_fast_flux_err=ferr, coil_heating=heat, fw_fast_flux=fwf, tbr=tbr,
         batches=BATCHES, particles=PARTICLES, n_source=BATCHES*PARTICLES, surf_file=src_out, cellfrom=CELLFROM)
print(f"COIL peak/mean={flux.max()/flux.mean():.3f} hottest_idx={int(np.argmax(flux))} "
      f"FW_fast={fwf:.4e} TBR={tbr:.4f}", flush=True)
print(f"WROTE {WORK}/coil_{OUT_TAG}.npz  and  {src_out}", flush=True)
print("DONE_" + OUT_TAG, flush=True)
