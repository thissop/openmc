#!/usr/bin/env python
"""FW-CADIS (random ray) weight-window generation for the coil->plasma RECIPROCITY run.

Replaces the failed MAGIC attempt (particle-splitting explosion, ~1h/batch).

Physics: emit neutrons FROM a coil (DAGMC cell, isotropic 14.06 MeV, box source
constrained to that cell), and tally flux on a RegularMesh over the plasma
region. That mesh flux is the "plasma importance map" -- which plasma regions can
reach the coil through ~1 m of shield+blanket. Because the flux drops ~6 decades,
we build FW-CADIS weight windows with the random-ray solver (forward + adjoint,
run automatically), targeting the plasma mesh tally, then apply those weight
windows in a continuous-energy fixed-source run.

WORKFLOW (see FWCADIS_PLAN.md for the verified API and failure points):
  1. Build the CE model (materials + DAGMC geometry + coil source + plasma tally).
  2. deepcopy -> convert_to_multigroup() -> convert_to_random_ray().
  3. Overlay the WW mesh as the random-ray source-region mesh (same mesh!).
  4. WeightWindowGenerator(method='fw_cadis', targets=[plasma tally]); model.run()
     -> forward+adjoint random-ray solve -> weight_windows.h5.
  5. Apply weight_windows.h5 to the CE model; run the reciprocity flux tally.
  6. Structure gate: assert the importance map is non-empty, localized, structured.

This is a runnable TEMPLATE. The geometry/material/source setup is copied from the
working MAGIC deck (../ww_pert_pipeline/coil_run.py); adjust paths, coil cell IDs,
plasma-mesh extent, and group structure to your device. Places that genuinely
depend on cluster data or a runtime API check are flagged with # >>> CONFIRM.

Anchored to OpenMC 0.15.4.dev (tag 0.15.3+182); Ginsburg fork is 0.15.1-dev.
"""
import argparse
import copy
import os
import sys
import time

import numpy as np
import openmc


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dagmc", required=True, help="Path to DAGMC .h5m file")
    p.add_argument("--workdir", default=os.path.expanduser("~/pstl_test/fwcadis"))
    # Coil cells: the reciprocity SOURCE is emitted from these DAGMC cells.
    # coil_run.py uses COIL_CELLS = range(8, 28); pick the coil(s) of interest.
    p.add_argument("--coil-cells", type=int, nargs="+", default=[8],
                   help="DAGMC cell ID(s) to emit the reciprocity source from")
    # MGXS generation
    p.add_argument("--mgxs-method", default="stochastic_slab",
                   choices=["material_wise", "stochastic_slab", "infinite_medium"],
                   help="stochastic_slab is safest for deep-penetration (no zero-XS "
                        "for far materials); material_wise is higher fidelity.")
    p.add_argument("--mgxs-particles", type=int, default=20000)
    p.add_argument("--mgxs-path", default="mgxs.h5")
    # WW / random-ray mesh resolution (same mesh used for both, per doc warning)
    p.add_argument("--nr", type=int, default=30)
    p.add_argument("--nphi", type=int, default=16)
    p.add_argument("--nz", type=int, default=30)
    # Random-ray WW-gen run
    p.add_argument("--rr-batches", type=int, default=25)
    p.add_argument("--rr-inactive", type=int, default=5)
    p.add_argument("--rr-particles", type=int, default=2000, help="rays per batch")
    # CE reciprocity run
    p.add_argument("--ce-batches", type=int, default=20)
    p.add_argument("--ce-particles", type=int, default=100000)
    p.add_argument("--threads", type=int,
                   default=int(os.environ.get("OMP_NUM_THREADS", 32)))
    p.add_argument("--skip-wwgen", action="store_true",
                   help="Reuse an existing weight_windows.h5 and only run step 5-6.")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# Materials  (copied verbatim from ../ww_pert_pipeline/coil_run.py)
# --------------------------------------------------------------------------- #
def build_materials(vacuum=True):
    # EXACT tag set of dagmc_corr_breed (matches coil_percoil5.py): 9 in-vessel + Vacuum.
    w = openmc.Material(); w.add_element("W", 1.0); w.set_density("g/cm3", 19.3)
    rafm = openmc.Material(); rafm.add_element("Fe", 1.0); rafm.set_density("g/cm3", 7.9)
    fw = openmc.Material.mix_materials([w, rafm], [0.2/3.2, 3.0/3.2], "vo"); fw.name = "first_wall"
    mult = openmc.Material(name="multiplier"); mult.add_element("Be", 1.0); mult.set_density("g/cm3", 1.85)
    br = openmc.Material(name="breeder")
    br.add_element("Li", 2.0, enrichment=60.0, enrichment_target="Li6", enrichment_type="ao")
    br.add_element("Be", 1.0); br.add_element("F", 4.0); br.set_density("g/cm3", 1.94)
    bw = openmc.Material(name="back_wall"); bw.add_element("Fe", 1.0); bw.set_density("g/cm3", 7.9)
    sh = openmc.Material(name="shield"); sh.add_element("W", 1.0); sh.add_element("C", 1.0); sh.set_density("g/cm3", 15.6)
    # Near-void densities floored at 1e-6 g/cm3 (MFP ~ 4e5 cm >> device ~4e3 cm, so still
    # physically void) so their MC-tallied random-ray total XS is DETERMINISTICALLY positive.
    # At 1e-8/1e-10 the total sits at the tally-noise floor and can come out <=0, which random
    # ray rejects ("No zero or negative total macroscopic cross sections") -- geometry/seed
    # dependent (w20 happened to pass, w35f/thicker builds fail). VERIFIED fix 2026-07-25.
    gap = openmc.Material(name="gap"); gap.add_element("H", 1.0); gap.set_density("g/cm3", 1e-6)
    vv = openmc.Material(name="vac_vessel")
    vv.add_element("Fe", 0.7, "wo"); vv.add_element("Cr", 0.2, "wo"); vv.add_element("Ni", 0.1, "wo"); vv.set_density("g/cm3", 7.9)
    ts = openmc.Material(name="thermal_shield"); ts.add_element("Fe", 0.9, "wo"); ts.add_element("Cr", 0.1, "wo"); ts.set_density("g/cm3", 7.9)
    mg = openmc.Material(name="magnets"); mg.add_element("Cu", 1.0); mg.set_density("g/cm3", 8.5)
    vac = openmc.Material(name="Vacuum"); vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-6)
    base = [fw, mult, br, bw, sh, gap, vv, ts, mg]
    # QA DAGMC has a 9-tag set (no Vacuum); QH has 10. Drop Vacuum for QA.
    return openmc.Materials(base + ([vac] if vacuum else []))


# --------------------------------------------------------------------------- #
# Geometry + source + tallies (the CE reciprocity model)
# --------------------------------------------------------------------------- #
def build_ce_model(args):
    mats = build_materials()

    dag = openmc.DAGMCUniverse(args.dagmc, auto_geom_ids=True)
    bb = dag.bounding_box
    ll, ur = np.array(bb.lower_left), np.array(bb.upper_right)
    ext = float(np.max(np.abs(np.concatenate([ll, ur]))))
    sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
    root = openmc.Cell(region=-sph, fill=dag)
    geom = openmc.Geometry([root])

    # ---- Reciprocity SOURCE: isotropic 14.06 MeV emitted FROM the coil cell(s).
    # Random-ray fixed source = IndependentSource with a domain constraint + MG
    # energy (verified: examples.py L1271-1276, random_ray_fixed_source_domain).
    # A Box spatial over the coil bbox is combined with the cell-domain constraint
    # so rejection sampling only keeps sites inside the coil cell.
    # tighten the source box to the coil vicinity using discovered centroids (efficient rejection)
    _cc = np.load("/burg-archive/home/tjk2147/pstl_test/corrected/coil_centroids_corr.npz")
    _cids = [int(x) for x in _cc["cell_ids"]]; _cents = _cc["centroids"]
    _sel = np.array([_cents[_cids.index(c)] for c in args.coil_cells])
    _cll = _sel.min(0) - 60.0; _cur = _sel.max(0) + 60.0   # ~coil cross-section
    _pt = _sel.mean(0)   # coil centroid (representative coil location)
    print("POINT source at coil centroid", args.coil_cells, "->", _pt, flush=True)
    # Random-ray mode requires a POINT source (or domain-constrained); DAGMC cells aren't
    # enumerable pre-run, so use a point at the coil centroid (coil is tiny vs coil->plasma dist).
    src = openmc.IndependentSource(
        space=openmc.stats.Point(tuple(float(v) for v in _pt)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([14.06e6], [1.0]),
        strength=1.0,
    )

    # ---- Plasma importance map: RegularMesh flux tally = FW-CADIS target + deliverable.
    # >>> CONFIRM the plasma envelope; here we use a Cartesian box around the
    # magnetic axis. For a device with major radius R0 and plasma minor radius a,
    # set lower_left/upper_right to bound the plasma region only (NOT the whole
    # device -- the target is the plasma).
    plasma_mesh = openmc.RegularMesh()
    plasma_mesh.dimension = (30, 30, 20)                 # QH plasma bore
    plasma_mesh.lower_left = (-1800.0, -1800.0, -550.0)  # QH plasma envelope (cm)
    plasma_mesh.upper_right = (1800.0, 1800.0, 550.0)

    plasma_filter = openmc.MeshFilter(plasma_mesh)
    plasma_tally = openmc.Tally(name="plasma_importance")
    plasma_tally.filters = [plasma_filter]
    plasma_tally.scores = ["flux"]
    tallies = openmc.Tallies([plasma_tally])

    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.source = [src]
    s.batches = args.ce_batches
    s.particles = args.ce_particles
    s.photon_transport = False
    s.max_lost_particles = 100000
    s.rel_max_lost_particles = 0.99
    s.output = {"summary": True, "tallies": False}

    model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tallies)
    return model, plasma_mesh, ext, (ll, ur)


# --------------------------------------------------------------------------- #
# Step 2-4: FW-CADIS weight-window generation via random ray
# --------------------------------------------------------------------------- #
def generate_weight_windows(ce_model, plasma_mesh, ext, bounds, args):
    ll, ur = bounds

    # Deep copy so the CE model stays pristine for the final reciprocity run.
    mg_model = copy.deepcopy(ce_model)

    # --- WW / source-region mesh. A CylindricalMesh spanning the whole device
    # (plasma out past the coils), mirroring the MAGIC gen_ww.py. This SAME mesh
    # is used for source_region_meshes AND the WW generator (doc warning:
    # variance_reduction.rst L150-155 -- WW mesh must not subdivide source regions).
    ww_mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0.0, ext, args.nr + 1),
        phi_grid=np.linspace(0.0, 2 * np.pi, args.nphi + 1),
        z_grid=np.linspace(ll[2] - 1, ur[2] + 1, args.nz + 1),
        origin=(0.0, 0.0, 0.0),
    )

    # --- (a) CE -> multigroup. DAGMC-aware (model.py L2758-2773). Writes mgxs.h5.
    print("=== [1] convert_to_multigroup (%s) ===" % args.mgxs_method, flush=True)
    t0 = time.time()
    mg_model.convert_to_multigroup(
        method=args.mgxs_method,
        groups=openmc.mgxs.EnergyGroups([1e-5,0.025,1.0,100.0,1e3,1e4,1e5,5e5,1e6,2e6,5e6,1e7,1.4e7,2e7]),  # fusion, fast detail
        # fusion-appropriate multi-group structure with FAST detail (>0.1 MeV).
        # e.g. groups=openmc.mgxs.EnergyGroups([...eV edges...]) or a named struct
        # from openmc.mgxs.GROUP_STRUCTURES. CASMO-2 is a coarse placeholder.
        nparticles=args.mgxs_particles,
        mgxs_path=args.mgxs_path,
        source_energy=openmc.stats.delta_function(14.06e6),  # DT collapse spectrum
        overwrite_mgxs_library=False,
    )
    print("MGXS_WALL_SECONDS %.1f" % (time.time() - t0), flush=True)

    # Guard: material_wise can leave zero total XS for materials the CE gen source
    # never reached (deep-penetration hazard). Flag it loudly.
    _check_mgxs_nonzero(args.mgxs_path)

    # --- (b) multigroup -> random ray. Auto-fills ray_source, distances, particles.
    print("=== [2] convert_to_random_ray ===", flush=True)
    mg_model.convert_to_random_ray()

    # Overlay the WW mesh as the source-region decomposition mesh on the root universe.
    root_universe = mg_model.geometry.root_universe
    mg_model.settings.random_ray["source_region_meshes"] = [(ww_mesh, [root_universe])]
    mg_model.settings.random_ray["volume_estimator"] = "naive"   # per fw_cadis tests
    # Override the ray count / batches (the auto guess is just a starting point).
    mg_model.settings.particles = args.rr_particles
    mg_model.settings.batches = args.rr_batches
    mg_model.settings.inactive = args.rr_inactive

    # --- (c) FW-CADIS generator targeting the plasma mesh tally (local VR).
    # The target tally object must be the one already inside mg_model.tallies
    # (deepcopy preserved it) so the export-time registration check passes.
    plasma_tally = None
    for t in mg_model.tallies:
        if t.name == "plasma_importance":
            plasma_tally = t
            break
    assert plasma_tally is not None, "plasma_importance tally missing in MG model"

    wwg = openmc.WeightWindowGenerator(
        method="fw_cadis",
        mesh=ww_mesh,
        targets=[plasma_tally.id],   # this fork: tally IDs (ints), not objects
        max_realizations=mg_model.settings.batches,
    )
    mg_model.settings.weight_window_generators = [wwg]

    print("=== [3] random-ray FW-CADIS solve (forward + adjoint) ===", flush=True)
    t0 = time.time()
    # Model.run() (not Settings.export_to_xml) so the FW-CADIS target registration
    # check runs (variance_reduction.rst L191-201).
    mg_model.run(output=True, threads=args.threads)
    print("WWGEN_WALL_SECONDS %.1f" % (time.time() - t0), flush=True)

    assert os.path.exists("weight_windows.h5"), \
        "FW-CADIS run did not produce weight_windows.h5"
    print("weight_windows.h5 written", flush=True)


def _check_mgxs_nonzero(mgxs_path):
    """Warn loudly if any material has all-zero total XS (deep-penetration hazard)."""
    try:
        import h5py
    except ImportError:
        return
    zero_mats = []
    with h5py.File(mgxs_path, "r") as f:
        for name in f.keys():
            grp = f[name]
            # MGXSLibrary layout: /<material>/<temperature>/total etc.
            found_total = False
            for key in grp:
                sub = grp[key]
                if hasattr(sub, "keys") and "total" in sub:
                    tot = np.array(sub["total"])
                    found_total = True
                    if not np.any(tot > 0):
                        zero_mats.append(name)
            if not found_total:
                # layout differs across versions; skip silently
                pass
    if zero_mats:
        print("WARNING: zero total MGXS for materials %s -- CE gen source likely "
              "never reached them. Switch --mgxs-method stochastic_slab." % zero_mats,
              file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# Step 5: CE reciprocity run with weight windows applied
# --------------------------------------------------------------------------- #
def run_reciprocity(ce_model, args):
    print("=== [4] CE reciprocity run with weight windows ===", flush=True)
    assert os.path.exists("weight_windows.h5"), "weight_windows.h5 not found"
    ce_model.settings.weight_windows = openmc.hdf5_to_wws("weight_windows.h5")
    ce_model.settings.weight_windows_on = True
    ce_model.settings.weight_window_checkpoints = {"collision": True, "surface": True}
    ce_model.settings.survival_biasing = False

    t0 = time.time()
    sp_path = ce_model.run(output=True, threads=args.threads)
    print("RECIP_WALL_SECONDS %.1f" % (time.time() - t0), flush=True)
    return sp_path


# --------------------------------------------------------------------------- #
# Step 6: structure gate on the plasma importance map
# --------------------------------------------------------------------------- #
def structure_gate(sp_path, plasma_mesh):
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="plasma_importance")
        flux = t.mean.ravel()
        err = t.std_dev.ravel()

    imp = flux.reshape(plasma_mesh.dimension)
    total = float(flux.sum())
    nnz = int(np.count_nonzero(flux))
    frac_nnz = nnz / flux.size
    peak = float(flux.max())
    mean_nz = float(flux[flux > 0].mean()) if nnz else 0.0
    peak_to_mean = peak / mean_nz if mean_nz > 0 else np.nan
    rel = np.where(flux > 0, err / np.maximum(flux, 1e-300), np.nan)
    med_rel = float(np.nanmedian(rel))

    print("\n--- STRUCTURE GATE (plasma importance map) ---")
    print(f"  mesh dims        = {plasma_mesh.dimension}")
    print(f"  total flux       = {total:.4e}")
    print(f"  nonzero voxels   = {nnz}/{flux.size}  ({frac_nnz:.1%})")
    print(f"  peak/mean(nz)    = {peak_to_mean:.2f}")
    print(f"  median rel.err   = {med_rel:.3f}")

    # (1) non-empty
    assert total > 0.0, "GATE FAIL: importance map is entirely zero (WW/transport failed)"
    assert nnz > 0, "GATE FAIL: no nonzero voxels"
    # (2) localized: a reciprocity importance map should NOT light up every voxel
    #     uniformly, nor be a single spike. Expect a structured, partial map.
    assert frac_nnz < 0.98, "GATE FAIL: map fills ~all voxels -- not localized"
    assert peak_to_mean > 1.5, "GATE FAIL: map is flat -- no importance structure"
    # (3) structured / converged enough to be meaningful
    assert med_rel < 0.5, ("GATE FAIL: median rel. error %.2f too high -- WW not "
                           "effective or too few histories" % med_rel)
    print("  GATE PASS: importance map is non-empty, localized, and structured.\n")

    np.savez("plasma_importance.npz",
             importance=imp, flux=flux, err=err,
             dimension=np.array(plasma_mesh.dimension),
             lower_left=np.array(plasma_mesh.lower_left),
             upper_right=np.array(plasma_mesh.upper_right))
    print("saved plasma_importance.npz", flush=True)


# --------------------------------------------------------------------------- #
def main():
    args = parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    os.chdir(args.workdir)
    print("workdir:", args.workdir, flush=True)

    ce_model, plasma_mesh, ext, bounds = build_ce_model(args)

    if not args.skip_wwgen:
        generate_weight_windows(ce_model, plasma_mesh, ext, bounds, args)
    else:
        print("skip-wwgen: reusing existing weight_windows.h5", flush=True)

    sp_path = run_reciprocity(ce_model, args)
    structure_gate(sp_path, plasma_mesh)
    print("DONE_FWCADIS_RECIPROCITY", flush=True)


if __name__ == "__main__":
    main()
