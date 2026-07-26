#!/usr/bin/env python
"""CORRECTED-build transport run (v3): identical physics to coil_run_v2, but the coil-cell list
is LOADED from a precomputed pymoab npz instead of the openmc.lib c.fill discovery -- because
this checkout's openmc.lib.Cell.fill raises NotImplementedError on DAGMC (segfaults at teardown).
Also removes stale XML in cwd so the run loads the DAGMC we pass, not a leftover model.xml.

Usage: python coil_run_v3.py MODE BATCHES PARTICLES DAGMC_FILE OUT_TAG LI6 CELLS_NPZ [WW_FILE]
  MODE in {unpol, A, B, C};  LI6 = Li-6 enrichment percent;  CELLS_NPZ has magnet_cells (from
  step1_cells.py). Writes coil_<OUT_TAG>.npz + prints per-coil table + TBR.
"""
import os, sys, shutil, time
import numpy as np
import openmc

os.chdir(os.path.expanduser("~/pstl_test/corrected"))
BASE = "/burg-archive/home/tjk2147/pstl_test"

MODE      = sys.argv[1]
BATCHES   = int(sys.argv[2])
PARTICLES = int(sys.argv[3])
DAGMC     = sys.argv[4]
OUT_TAG   = sys.argv[5]
LI6       = float(sys.argv[6]) if len(sys.argv) > 6 else 40.0
CELLS_NPZ = sys.argv[7]
WW_FILE   = sys.argv[8] if len(sys.argv) > 8 else None

# --- clear stale XML so OpenMC does not prefer a leftover model.xml over this geometry ---
for f in ("model.xml", "geometry.xml", "materials.xml", "settings.xml", "tallies.xml"):
    try: os.remove(f)
    except FileNotFoundError: pass

POL = {"unpol": (1/3, 1/3, 1/3), "A": (1.0, 0.0, 0.0),
       "B": (0.0, 1.0, 0.0), "C": (0.0, 0.0, 1.0)}[MODE]


def build_materials(li6):
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

    return openmc.Materials([fw, mult, br, bw, sh, gap, vv, ts, mg, vac]), br


mats, breeder_mat = build_materials(LI6)

dag = openmc.DAGMCUniverse(DAGMC, auto_geom_ids=True)
bb = dag.bounding_box
ll, ur = np.array(bb.lower_left), np.array(bb.upper_right)
ext = float(np.max(np.abs(np.concatenate([ll, ur]))))
sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
root = openmc.Cell(region=-sph, fill=dag)
geom = openmc.Geometry([root])

# ---------------- coil cells: LOADED (pymoab precompute), not c.fill ----------------
cd = np.load(CELLS_NPZ)
COIL_CELLS = sorted(int(x) for x in cd["magnet_cells"])
assert COIL_CELLS, "no magnet cells in " + CELLS_NPZ
print(f"LOADED {len(COIL_CELLS)} coil cells from {CELLS_NPZ}: {COIL_CELLS}", flush=True)

# ---------------- tallies ----------------
fast = openmc.EnergyFilter([0.1e6, 25.0e6])
coil_filter = openmc.CellFilter(COIL_CELLS)
t_flux = openmc.Tally(name="coil_fast_flux"); t_flux.filters = [coil_filter, fast]; t_flux.scores = ["flux"]
t_heat = openmc.Tally(name="coil_heating"); t_heat.filters = [coil_filter]; t_heat.scores = ["heating-local"]
t_fw = openmc.Tally(name="fw_fast_flux")
t_fw.filters = [openmc.MaterialFilter([mats[0]]), fast]; t_fw.scores = ["flux"]
t_tbr = openmc.Tally(name="tbr")
t_tbr.filters = [openmc.MaterialFilter([breeder_mat])]; t_tbr.scores = ["H3-production"]
t_tbr_nuc = openmc.Tally(name="tbr_nuc")
t_tbr_nuc.filters = [openmc.MaterialFilter([breeder_mat])]
t_tbr_nuc.nuclides = ["Li6", "Li7"]; t_tbr_nuc.scores = ["H3-production"]
tallies = openmc.Tallies([t_flux, t_heat, t_fw, t_tbr, t_tbr_nuc])

src = openmc.StellaratorSource(
    fluxmap=f"{BASE}/qh_fluxmap", polarization=POL, field_model="fluxmap",
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

model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tallies)
print(f"=== MODE={MODE} TAG={OUT_TAG} DAGMC={os.path.basename(DAGMC)} LI6={LI6} "
      f"batches={BATCHES} particles={PARTICLES} ===", flush=True)
t0 = time.time()
sp_path = model.run(output=True, threads=int(os.environ.get("OMP_NUM_THREADS", 32)))
print(f"RUN_WALL_SECONDS {time.time()-t0:.1f}", flush=True)
dst = f"sp_{OUT_TAG}.h5"; shutil.copy(sp_path, dst)

with openmc.StatePoint(dst) as sp:
    tf = sp.get_tally(name="coil_fast_flux"); flux = tf.mean.ravel(); ferr = tf.std_dev.ravel()
    th = sp.get_tally(name="coil_heating"); heat = th.mean.ravel(); herr = th.std_dev.ravel()
    fw = sp.get_tally(name="fw_fast_flux"); fwf = fw.mean.ravel()[0]; fwe = fw.std_dev.ravel()[0]
    tbr = sp.get_tally(name="tbr").mean.ravel()[0]; tbre = sp.get_tally(name="tbr").std_dev.ravel()[0]
    tn = sp.get_tally(name="tbr_nuc"); tbr_nuc = tn.mean.ravel(); tbr_nuc_e = tn.std_dev.ravel()

rel = np.where(flux > 0, ferr / np.maximum(flux, 1e-30), np.nan)
print("\nCOIL_CELL  FAST_FLUX      RELERR   HEATING        RELERR")
for i, cid in enumerate(COIL_CELLS):
    hr = herr[i] / heat[i] if heat[i] > 0 else np.nan
    print(f"{cid:8d}  {flux[i]:.4e}   {rel[i]:6.4f}  {heat[i]:.4e}   {hr:6.4f}")
peak = flux.max() / flux.mean() if flux.mean() > 0 else np.nan
ihot = int(np.argmax(flux))
print(f"\nCOIL_PEAKING(max/mean) = {peak:.4f}")
print(f"HOTTEST_COIL_CELL = {COIL_CELLS[ihot]}  flux={flux.max():.4e}  relerr={rel[ihot]:.4f}")
print(f"FW_FAST_FLUX = {fwf:.4e} +/- {fwe:.4e}")
print(f"TBR = {tbr:.4f} +/- {tbre:.4f}   (rel {tbre/tbr*100:.2f}%)")
print(f"  Li6 H3-prod = {tbr_nuc[0]:.4f} +/- {tbr_nuc_e[0]:.4f}")
print(f"  Li7 H3-prod = {tbr_nuc[1]:.4f} +/- {tbr_nuc_e[1]:.4f}")
np.savez(f"coil_{OUT_TAG}.npz",
    mode=MODE, tag=OUT_TAG, pol=np.array(POL), li6=LI6, coil_cells=np.array(COIL_CELLS),
    coil_fast_flux=flux, coil_fast_flux_err=ferr, coil_heating=heat, coil_heating_err=herr,
    fw_fast_flux=fwf, fw_fast_flux_err=fwe, tbr=tbr, tbr_err=tbre,
    tbr_li6=tbr_nuc[0], tbr_li7=tbr_nuc[1],
    peaking=peak, hottest_cell=COIL_CELLS[ihot], batches=BATCHES, particles=PARTICLES)
print("NPZ coil_" + OUT_TAG + ".npz", flush=True)
print("DONE_" + OUT_TAG, flush=True)
