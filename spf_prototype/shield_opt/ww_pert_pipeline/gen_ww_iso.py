#!/usr/bin/env python
"""PHASE 1 (isolated): right-sized MAGIC weight-window generation.
cwd = $WORKDIR. Uses the CLEAN original 9x9 baseline geometry. Writes weight_windows.h5.
"""
import os, sys, time
import numpy as np
import openmc

WORKDIR = os.environ.get("WORKDIR", os.path.expanduser("~/pstl_test/kiso"))
os.chdir(WORKDIR)
BASE = "/burg-archive/home/tjk2147/pstl_test"

BATCHES   = int(sys.argv[1]) if len(sys.argv) > 1 else 8
PARTICLES = int(sys.argv[2]) if len(sys.argv) > 2 else 25000
MAXREAL   = int(sys.argv[3]) if len(sys.argv) > 3 else 8

def build_materials():
    fw = openmc.Material(name="first_wall"); fw.add_element("W", 1.0); fw.set_density("g/cm3", 19.3)
    br = openmc.Material(name="breeder")
    br.add_element("Li", 2.0, enrichment=90.0, enrichment_target="Li6", enrichment_type="ao")
    br.add_element("Be", 1.0); br.add_element("F", 4.0); br.set_density("g/cm3", 1.94)
    bw = openmc.Material(name="back_wall")
    bw.add_element("Fe", 0.89, "wo"); bw.add_element("Cr", 0.09, "wo"); bw.add_element("Mn", 0.02, "wo"); bw.set_density("g/cm3", 7.9)
    sh = openmc.Material(name="shield"); sh.add_element("W", 1.0); sh.add_element("C", 1.0); sh.set_density("g/cm3", 15.6)
    vv = openmc.Material(name="vac_vessel")
    vv.add_element("Fe", 0.89, "wo"); vv.add_element("Cr", 0.09, "wo"); vv.add_element("Mn", 0.02, "wo"); vv.set_density("g/cm3", 7.9)
    mg = openmc.Material(name="magnets")
    mg.add_element("Cu", 0.58, "wo"); mg.add_element("Fe", 0.30, "wo"); mg.add_element("Cr", 0.04, "wo")
    mg.add_element("Ni", 0.03, "wo"); mg.add_element("Y", 0.02, "wo"); mg.add_element("Ba", 0.02, "wo"); mg.add_element("O", 0.01, "wo")
    mg.set_density("g/cm3", 8.5)
    vac = openmc.Material(name="Vacuum"); vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-10)
    return openmc.Materials([fw, br, bw, sh, vv, mg, vac])

mats = build_materials()
dag = openmc.DAGMCUniverse(f"{BASE}/dagmc_baseline.h5m", auto_geom_ids=True)
bb = dag.bounding_box
ll, ur = np.array(bb.lower_left), np.array(bb.upper_right)
ext = float(np.max(np.abs(np.concatenate([ll, ur]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
root = openmc.Cell(region=-sph, fill=dag)
geom = openmc.Geometry([root])

mesh = openmc.CylindricalMesh(
    r_grid=np.linspace(0.0, ext, 26),
    phi_grid=np.linspace(0.0, 2*np.pi, 17),
    z_grid=np.linspace(ll[2]-1, ur[2]+1, 26),
    origin=(0.0, 0.0, 0.0))
wwg = openmc.WeightWindowGenerator(
    mesh=mesh, energy_bounds=[1.0, 1.0e5, 2.5e7], particle_type="neutron",
    method="magic", max_realizations=MAXREAL, update_interval=1, on_the_fly=True)

src = openmc.StellaratorSource(
    fluxmap=f"{BASE}/qh_fluxmap", polarization=(1/3, 1/3, 1/3),
    field_model="fluxmap", energy=openmc.stats.Discrete([14.06e6], [1.0]), strength=1.0)

s = openmc.Settings()
s.run_mode = "fixed source"; s.source = [src]
s.batches = BATCHES; s.particles = PARTICLES
s.photon_transport = False
s.max_lost_particles = 100000; s.rel_max_lost_particles = 0.99
s.weight_window_generators = [wwg]
s.output = {"summary": False, "tallies": False}
model = openmc.Model(geometry=geom, materials=mats, settings=s)
print(f"=== WWGEN iso batches={BATCHES} particles={PARTICLES} maxreal={MAXREAL} ext={ext:.0f} ===", flush=True)
t0 = time.time()
model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"WW_GEN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)
print("weight_windows.h5 exists:", os.path.exists("weight_windows.h5"), flush=True)
print("DONE_WWGEN_ISO", flush=True)
