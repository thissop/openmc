#!/usr/bin/env python
"""One-field-period SSR replay: read the deep-freeze surface source (banked at the fixed-breeder
boundary in op_master), transport it through the variable zone to the coils on the SAME periodic
wedge, and tally the coil dose for Gate 2.

BLACK_INSIDE=1 (default) overrides breeder_fixed (cell 5) with a strong B-10 absorber so any replayed
neutron that backscatters inward across the freeze surface DIES instead of streaming back across the
plasma void -- the deep-freeze analogue of the pristine fix. BLACK_INSIDE=0 leaves the fixed zone
as-is (void vacuum) for comparison.

Usage: op_ssr.py BATCHES PARTICLES SURF_SRC [LI6]
"""
import os, sys, time
import numpy as np
import openmc

BASE = "/burg-archive/home/tjk2147/pstl_test"
WORK = f"{BASE}/corrected"
RUNDIR = os.environ.get("RUNDIR", WORK); os.makedirs(RUNDIR, exist_ok=True); os.chdir(RUNDIR)
BATCHES = int(sys.argv[1]); PARTICLES = int(sys.argv[2]); SURF = sys.argv[3]
LI6 = float(sys.argv[4]) if len(sys.argv) > 4 else 60.0
DAGMC = f"{WORK}/dagmc_qh_op_uniform.h5m"
MAGCELLS = [12, 13, 14, 15, 16]
BLACK = os.environ.get("BLACK_INSIDE", "1") == "1"
FIXED_CELLS = [2, 3, 4]   # void, first wall, multiplier -> black (INSIDE the freeze; kills inward neutrons before the void, leaves cell 5/6 normal so outward source transports out)

for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml"):
    try: os.remove(f)
    except FileNotFoundError: pass
assert os.path.exists(SURF), f"surface source not found: {SURF}"


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
void_abs = None
if BLACK:
    void_abs = openmc.Material(name="fixed_black")
    void_abs.add_nuclide("B10", 1.0); void_abs.set_density("g/cm3", 10.0)   # ~black over the 15 cm fixed zone
    mats.append(void_abs)

dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
ext = float(np.max(np.abs(np.concatenate([np.array(bb.lower_left), np.array(bb.upper_right)]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
py = openmc.YPlane(0.0, boundary_type="periodic")
px = openmc.XPlane(0.0, boundary_type="periodic")
py.periodic_surface = px
root = openmc.Cell(region=+py & +px & -sph, fill=dag)
geom = openmc.Geometry([root])

s = openmc.Settings()
s.run_mode = "fixed source"
s.surf_source_read = {"path": SURF}
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False
s.max_lost_particles = 500; s.rel_max_lost_particles = 2e-3
s.output = {"summary": True, "tallies": False}

fast = openmc.EnergyFilter([0.1e6, 25.0e6]); cf = openmc.CellFilter(MAGCELLS)
t_flux = openmc.Tally(name="coil_fast_flux"); t_flux.filters = [cf, fast]; t_flux.scores = ["flux"]
t_heat = openmc.Tally(name="coil_heating"); t_heat.filters = [cf]; t_heat.scores = ["heating-local"]
# DIAGNOSTIC: TOTAL flux (no energy cut) at increasing radial depth -> where does the flux die?
t_diag = openmc.Tally(name="depth_flux")
t_diag.filters = [openmc.CellFilter([6, 8, 10, 12, 14])]; t_diag.scores = ["flux"]   # breeder_var, shield, vv, coils
tals = openmc.Tallies([t_flux, t_heat, t_diag])

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tals)

if BLACK:
    # override breeder_fixed (cell 5) -> black, via the synced-DAGMC path (as in ssr_dose void-kill)
    model.init_lib(output=False)
    model.sync_dagmc_universes()
    for fc in FIXED_CELLS:
        dag.add_material_override(fc, void_abs)
    model.finalize_lib()
    for fc in FIXED_CELLS:
        got = getattr(dag.cells[fc].fill, "name", None)
        assert got == "fixed_black", f"BLACK override failed: cell {fc} fill is {got!r}"
    print(f"BLACK_INSIDE: cells {FIXED_CELLS} -> B-10 black", flush=True)
else:
    print("BLACK_INSIDE=0: fixed zone left as-is", flush=True)

print(f"=== OP SSR replay src={os.path.basename(SURF)} black={BLACK} batches={BATCHES} "
      f"particles={PARTICLES} cwd={RUNDIR} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)

with openmc.StatePoint(sp_path) as sp:
    fl = sp.get_tally(name="coil_fast_flux"); flux = fl.mean.ravel(); ferr = fl.std_dev.ravel()
    heat = sp.get_tally(name="coil_heating").mean.ravel()

tag = "op_ssr_black" if BLACK else "op_ssr_vac"
np.savez(f"{WORK}/coil_{tag}.npz", tag=tag, coil_cells=np.array(MAGCELLS),
         coil_fast_flux=flux, coil_fast_flux_err=ferr, coil_heating=heat,
         batches=BATCHES, particles=PARTICLES, surf_src=SURF, black=BLACK)
with openmc.StatePoint(sp_path) as sp:
    dg = sp.get_tally(name="depth_flux").mean.ravel()
print("DEPTH total flux [breeder_var6, shield8, vv10, coil12, coil14]:", np.array2string(dg, precision=4), flush=True)
print("COIL flux:", np.array2string(flux, precision=4), flush=True)
print(f"WROTE {WORK}/coil_{tag}.npz", flush=True)
print("DONE_op_ssr", flush=True)
