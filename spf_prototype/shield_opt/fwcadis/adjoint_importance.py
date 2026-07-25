#!/usr/bin/env python
"""Adjoint magnet-importance map via OpenMC random-ray (Peterson recipe, 2026-07-23).

This is the PRINCIPLED version of the coil->plasma attribution. Instead of the
reciprocity shortcut (emit isotropic 14 MeV FROM the coil, tally plasma flux --
see reciprocity_fwcadis.py), we localize an ADJOINT SOURCE = the magnet response
function inside the coil cell and run the random-ray ADJOINT solver. The adjoint
flux everywhere IS the importance map: psi_dagger(r) = how much a neutron born at
plasma point r contributes to the chosen magnet response (heating / fast flux),
THROUGH the ~1 m of shield+blanket. Energy-resolved and exact for the response,
with no isotropic-monoenergetic emission assumption.

Adjoint source energy dependence (Peterson gave two options):
  --response flat   : uniform over [E_lo, E_hi]  (simpler starting point)
  --response kerma  : weighted by the Cu (magnet) KERMA(E) -- nuclear-heating
                      response function; the physically correct magnet-heating
                      adjoint source.

Why this runs unchanged on the existing Ginsburg build: the localized adjoint
source is already in mainline (`Local adjoint source for Random Ray`, PR #3717,
+ stabilization PRs #3476/#3449/#3325). FlatSourceDomain::set_local_adjoint_sources()
consumes settings.random_ray['adjoint_source']. No rebuild against a feature
branch is needed. (Jack Fletcher's local_adjoint branch adds only convenience UI:
CADIS-on-WWG linking + worked examples.)

WORKFLOW:
  1. Build the CE model: materials + DAGMC geometry + plasma-mesh flux tally.
  2. convert_to_multigroup() -> convert_to_random_ray()  (same as reciprocity).
  3. settings.random_ray['adjoint'] = True
     settings.random_ray['adjoint_source'] = [ localized magnet-response source ].
  4. model.run() -> adjoint random-ray solve -> plasma-mesh tally = importance map.
  5. Structure gate + save adjoint_importance.npz  (drop-in for the shield optimizer,
     replacing the free-streaming culprit map with a through-shield importance map).

RUN ENVIRONMENT: Ginsburg (no runnable OpenMC locally; transport is x86/cluster).
Places that depend on cluster data / a device envelope are flagged  # >>> CONFIRM.

Anchored to OpenMC 0.15.x (localized-adjoint API). Sibling of reciprocity_fwcadis.py.
"""
import argparse
import copy
import os
import sys
import time

import numpy as np
import openmc

# Reuse the exact material stack from the reciprocity deck (magnets = Cu).
from reciprocity_fwcadis import build_materials, _check_mgxs_nonzero


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dagmc", required=True, help="Path to DAGMC .h5m file")
    p.add_argument("--workdir", default=os.path.expanduser("~/pstl_test/adjoint"))
    p.add_argument("--coil-cells", type=int, nargs="+", default=[8],
                   help="DAGMC cell ID(s) whose coil guide curve hosts the adjoint source")
    p.add_argument("--coils-file",
                   default="/burg-archive/home/tjk2147/pstl_test/coils_qh",
                   help="MAKEGRID coils file with the coil guide curves (filaments)")
    p.add_argument("--coils-scale", type=float, default=100.0,
                   help="scale to bring the coils file into the DAGMC frame (m->cm = 100)")
    p.add_argument("--coil-centroids",
                   default="/burg-archive/home/tjk2147/pstl_test/corrected/coil_centroids_corr.npz",
                   help="npz (cell_ids, centroids) mapping DAGMC coil cells -> filaments")
    p.add_argument("--filament-points", type=int, default=40,
                   help="number of point sources sampled along each coil guide curve")
    p.add_argument("--response", choices=["flat", "kerma"], default="flat",
                   help="Adjoint-source energy dependence: flat over a band, or "
                        "Cu-KERMA-weighted (magnet nuclear-heating response).")
    p.add_argument("--e-lo", type=float, default=1.0e5,
                   help="Low edge of the flat adjoint-source band (eV). Default 0.1 MeV "
                        "(fast neutrons dominate Cu heating & damage).")
    p.add_argument("--e-hi", type=float, default=1.5e7,
                   help="High edge of the flat adjoint-source band (eV). Default 15 MeV.")
    # MGXS generation (identical knobs to reciprocity)
    p.add_argument("--mgxs-method", default="stochastic_slab",
                   choices=["material_wise", "stochastic_slab", "infinite_medium"])
    p.add_argument("--mgxs-particles", type=int, default=20000)
    p.add_argument("--mgxs-path", default="mgxs.h5")
    # Plasma importance mesh (the deliverable) + random-ray source-region mesh
    p.add_argument("--nr", type=int, default=30)
    p.add_argument("--nphi", type=int, default=16)
    p.add_argument("--nz", type=int, default=30)
    # Random-ray adjoint solve
    p.add_argument("--rr-batches", type=int, default=40)
    p.add_argument("--rr-inactive", type=int, default=10)
    p.add_argument("--rr-particles", type=int, default=3000, help="rays per batch")
    p.add_argument("--threads", type=int,
                   default=int(os.environ.get("OMP_NUM_THREADS", 32)))
    return p.parse_args()


# Fusion-appropriate group structure with FAST detail (>0.1 MeV), shared with the
# reciprocity deck so the two importance maps are directly comparable.
GROUP_EDGES = [1e-5, 0.025, 1.0, 100.0, 1e3, 1e4, 1e5, 5e5, 1e6, 2e6, 5e6, 1e7, 1.4e7, 2e7]


# --------------------------------------------------------------------------- #
# The localized magnet-response adjoint source
# --------------------------------------------------------------------------- #
def _parse_coils(path, scale):
    """MAKEGRID coils file -> list of (Npts,3) filament polylines, scaled (e.g. m->cm).

    Each coil is a run of 'x y z current' rows terminated by a row with current == 0.
    """
    coils, cur = [], []
    for line in open(path):
        p = line.split()
        if len(p) >= 4:
            try:
                x, y, z, c = map(float, p[:4])
            except ValueError:
                continue                                   # header lines (periods/begin/mirror)
            cur.append((x, y, z))
            if c == 0.0:
                coils.append(np.array(cur) * scale)
                cur = []
    return coils


def _response_energy(args):
    """DISCRETE (multigroup) adjoint-source energy = the response weighting at group midpoints.

    Random-ray adjoint sources REQUIRE a discrete distribution (VERIFIED 2026-07-24: a
    Uniform/Tabular energy aborts with "Only discrete (multigroup) energy distributions
    are allowed for adjoint sources in random ray mode"). flat = equal weight to groups
    in [e_lo, e_hi]; kerma = Cu KERMA(E) per group (magnet nuclear-heating response).
    """
    edges = np.asarray(GROUP_EDGES, dtype=float)
    mids = np.sqrt(edges[:-1] * edges[1:])
    if args.response == "flat":
        w = ((mids >= args.e_lo) & (mids <= args.e_hi)).astype(float)
        if w.sum() == 0:
            w[:] = 1.0
    else:
        w = _cu_kerma_weights(mids)
    p = w / w.sum()
    nz = p > 0
    return openmc.stats.Discrete(mids[nz], p[nz])


def build_adjoint_source(args):
    """Adjoint source = magnet response as a FILAMENTARY source ON the coil guide curve.

    A QH coil is a non-planar curved tube; its heating response lives ALONG the winding.
    A point/box at the coil-cell CENTROID is unreliable (the centroid of a curved coil's
    volume can sit off the conductor -> VERIFIED 2026-07-24 the cell-15 centroid map
    pointed AWAY from the coil), and DAGMC cells are not enumerable pre-run so a per-cell
    domain constraint is impossible (dag.cells is empty). Instead we localize the response
    ON the coil's guide curve:
      1. parse the coil filament (MAKEGRID coils file) and align it to the DAGMC frame
         (scale m->cm; VERIFIED 2026-07-24 the scaled filaments match the DAGMC coil
         centroids BIJECTIVELY to ~35-50 cm, i.e. one filament per coil cell);
      2. pick each target DAGMC coil cell's filament by nearest centroid;
      3. subsample the guide curve to N points and emit a LIST of Point adjoint sources.
    Random ray SUMS multiple localized adjoint sources (set_local_adjoint_sources over
    model::adjoint_sources), so the list == the filamentary coil. Each point carries the
    response spectrum. This puts the source exactly on the conductor path everywhere it
    runs -- no off-conductor centroid, no axis-aligned box that could catch a neighbor.
    """
    filaments = _parse_coils(args.coils_file, args.coils_scale)
    fcen = np.array([f.mean(0) for f in filaments])
    cc = np.load(args.coil_centroids)
    cids = [int(x) for x in cc["cell_ids"]]
    dcen = cc["centroids"]
    energy = _response_energy(args)
    print(f"adjoint response = {args.response}", flush=True)

    srcs = []
    for cell in args.coil_cells:
        dci = dcen[cids.index(cell)]
        j = int(np.linalg.norm(fcen - dci, axis=1).argmin())
        fil = filaments[j]
        idx = np.unique(np.linspace(0, len(fil) - 1,
                                    min(args.filament_points, len(fil))).round().astype(int))
        pts = fil[idx]
        print(f"  coil cell {cell} -> filament {j}: {len(pts)} pts on guide curve, "
              f"centroid offset {np.linalg.norm(fcen[j] - dci):.0f} cm", flush=True)
        for pt in pts:
            srcs.append(openmc.IndependentSource(
                space=openmc.stats.Point(tuple(float(v) for v in pt)),
                angle=openmc.stats.Isotropic(), energy=energy, strength=1.0))
    print(f"adjoint source = FILAMENT, {len(srcs)} point sources total", flush=True)
    return srcs


def _cu_kerma_weights(mids):
    """Per-group Cu KERMA (heating) weights evaluated at group midpoints.

    Uses the CE library reachable via OPENMC_CROSS_SECTIONS on the cluster. The
    magnet-heating adjoint response is proportional to the material KERMA, so we
    weight each group's adjoint-source probability by the mid-group Cu kerma
    coefficient. Falls back to a flat weighting (with a loud warning) if kerma data
    is absent, so a run never silently misrepresents the response.
    """
    try:
        import openmc.data
        lib = openmc.data.DataLibrary.from_xml()
        cu = openmc.data.IncidentNeutron.from_hdf5(lib.get_by_material("Cu63")["path"])
        kerma = cu[301].xs["294K"](mids)                 # MT=301 = total heating (KERMA)
        w = np.clip(kerma, 0.0, None)
        if not np.any(w > 0):
            raise ValueError("Cu kerma evaluated to all-zero")
        return w
    except Exception as e:  # pragma: no cover -- data-dependent, cluster only
        print(f"WARNING: Cu kerma unavailable ({e!r}); falling back to FLAT "
              f"response weighting. Magnet-heating weighting is NOT applied.",
              file=sys.stderr, flush=True)
        return np.ones_like(mids)


# --------------------------------------------------------------------------- #
# CE base model: geometry + plasma-mesh flux tally (the importance deliverable)
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

    # Plasma importance mesh = where we read the adjoint flux (the attribution map).
    # >>> CONFIRM the plasma envelope per device (default: QH plasma bore, cm).
    plasma_mesh = openmc.RegularMesh()
    plasma_mesh.dimension = (30, 30, 20)
    plasma_mesh.lower_left = (-1800.0, -1800.0, -550.0)
    plasma_mesh.upper_right = (1800.0, 1800.0, 550.0)

    plasma_tally = openmc.Tally(name="magnet_importance")
    plasma_tally.filters = [openmc.MeshFilter(plasma_mesh)]
    plasma_tally.scores = ["flux"]        # in adjoint mode this scores the ADJOINT flux
    tallies = openmc.Tallies([plasma_tally])

    # A forward placeholder source is required for model construction but is
    # unused once adjoint mode + a localized adjoint_source are set.
    placeholder = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([14.06e6], [1.0]), strength=1.0)

    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.source = [placeholder]
    s.batches = args.rr_batches
    s.particles = args.rr_particles
    s.photon_transport = False
    s.output = {"summary": True, "tallies": False}

    model = openmc.Model(geometry=geom, materials=mats, settings=s, tallies=tallies)
    return model, plasma_mesh, ext, (ll, ur)


# --------------------------------------------------------------------------- #
# MG -> random ray -> ADJOINT solve
# --------------------------------------------------------------------------- #
def run_adjoint(ce_model, plasma_mesh, ext, bounds, args):
    ll, ur = bounds
    mg_model = copy.deepcopy(ce_model)

    ww_mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0.0, ext, args.nr + 1),
        phi_grid=np.linspace(0.0, 2 * np.pi, args.nphi + 1),
        z_grid=np.linspace(ll[2] - 1, ur[2] + 1, args.nz + 1),
        origin=(0.0, 0.0, 0.0),
    )

    print("=== [1] convert_to_multigroup (%s) ===" % args.mgxs_method, flush=True)
    t0 = time.time()
    mg_model.convert_to_multigroup(
        method=args.mgxs_method,
        groups=openmc.mgxs.EnergyGroups(GROUP_EDGES),
        nparticles=args.mgxs_particles,
        mgxs_path=args.mgxs_path,
        source_energy=openmc.stats.delta_function(14.06e6),
        overwrite_mgxs_library=False,
    )
    print("MGXS_WALL_SECONDS %.1f" % (time.time() - t0), flush=True)
    _check_mgxs_nonzero(args.mgxs_path)

    print("=== [2] convert_to_random_ray ===", flush=True)
    mg_model.convert_to_random_ray()
    rr = mg_model.settings.random_ray
    rr["source_region_meshes"] = [(ww_mesh, [mg_model.geometry.root_universe])]
    rr["volume_estimator"] = "naive"

    # --- THE ADJOINT RECIPE ---------------------------------------------------
    rr["adjoint"] = True
    rr["adjoint_source"] = build_adjoint_source(args)
    # -------------------------------------------------------------------------
    mg_model.settings.particles = args.rr_particles
    mg_model.settings.batches = args.rr_batches
    mg_model.settings.inactive = args.rr_inactive

    print("=== [3] random-ray ADJOINT solve ===", flush=True)
    t0 = time.time()
    sp_path = mg_model.run(output=True, threads=args.threads)
    print("ADJOINT_WALL_SECONDS %.1f" % (time.time() - t0), flush=True)
    return sp_path


# --------------------------------------------------------------------------- #
# Structure gate + save (drop-in for the shield optimizer)
# --------------------------------------------------------------------------- #
def structure_gate(sp_path, plasma_mesh, args):
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="magnet_importance")
        flux = t.mean.ravel()
        err = t.std_dev.ravel()

    imp = flux.reshape(plasma_mesh.dimension)
    total = float(flux.sum())
    nnz = int(np.count_nonzero(flux))
    frac_nnz = nnz / flux.size
    mean_nz = float(flux[flux > 0].mean()) if nnz else 0.0
    peak_to_mean = float(flux.max()) / mean_nz if mean_nz > 0 else np.nan
    rel = np.where(flux > 0, err / np.maximum(flux, 1e-300), np.nan)
    med_rel = float(np.nanmedian(rel))

    print("\n--- STRUCTURE GATE (adjoint magnet-importance map) ---")
    print(f"  response         = {args.response}")
    print(f"  mesh dims        = {plasma_mesh.dimension}")
    print(f"  total adj. flux  = {total:.4e}")
    print(f"  nonzero voxels   = {nnz}/{flux.size}  ({frac_nnz:.1%})")
    print(f"  peak/mean(nz)    = {peak_to_mean:.2f}")
    print(f"  median rel.err   = {med_rel:.3f}")

    assert total > 0.0, "GATE FAIL: adjoint importance map is entirely zero"
    assert frac_nnz < 0.98, "GATE FAIL: map fills ~all voxels -- not localized"
    assert peak_to_mean > 1.5, "GATE FAIL: map is flat -- no importance structure"
    assert med_rel < 0.5, ("GATE FAIL: median rel. error %.2f too high -- more "
                           "rays/batches needed" % med_rel)
    print("  GATE PASS: adjoint importance map is non-empty, localized, structured.\n")

    out = f"adjoint_importance_{args.response}.npz"
    np.savez(out, importance=imp, flux=flux, err=err,
             dimension=np.array(plasma_mesh.dimension),
             lower_left=np.array(plasma_mesh.lower_left),
             upper_right=np.array(plasma_mesh.upper_right),
             response=args.response, coil_cells=np.array(args.coil_cells))
    print("saved", out, flush=True)


def main():
    args = parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    os.chdir(args.workdir)
    print("workdir:", args.workdir, flush=True)

    ce_model, plasma_mesh, ext, bounds = build_ce_model(args)
    sp_path = run_adjoint(ce_model, plasma_mesh, ext, bounds, args)
    structure_gate(sp_path, plasma_mesh, args)
    print("DONE_ADJOINT_IMPORTANCE", flush=True)


if __name__ == "__main__":
    main()
