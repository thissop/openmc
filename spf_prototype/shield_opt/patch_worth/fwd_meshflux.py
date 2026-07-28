#!/usr/bin/env python
"""Tier-0 prerequisite: forward flux phi(x) on the SAME mesh as the adjoint psi_dagger,
so patch_worth.py can form the contributon C = phi*psi_dagger in the shield.

Runs the QH StellaratorSource (uniform emissivity, to match the forward per-coil runs and
the plasma source S used for attribution) on the BASELINE 8-layer DAGMC, tallying neutron
flux on a RegularMesh whose extent/dimension are read straight from the adjoint npz -> the
forward and adjoint fields are co-located voxel-for-voxel. Uses the FW-CADIS weight windows
so the deep-shield flux is resolved cheaply. Writes qh_fwd_meshflux.npz with key 'flux'
(shape = adjoint 'dimension', flattened C-order to match psi.reshape in patch_worth.py).

NOTE: the SOURCE block below must match the existing forward per-coil runs
(percoil_unpol_corr_w35f) -- same StellaratorSource params, same energy, uniform emissivity.
Edit SOURCE_LIB/params to point at spf_work/openmc_src's StellaratorSource if the API differs.
"""
import argparse
import os

import numpy as np
import openmc


def build(args):
    da = np.load(args.adjoint)
    ll = np.asarray(da["lower_left"], float)
    ur = np.asarray(da["upper_right"], float)
    dim = [int(v) for v in np.asarray(da["dimension"])]

    model = openmc.Model()
    # --- geometry: baseline DAGMC (same materials as step1 build) ---
    dag = openmc.DAGMCUniverse(args.dagmc).bounded_universe()
    model.geometry = openmc.Geometry(dag)
    # materials.xml is reused from the baseline run dir (must match the DAGMC tags)
    model.materials = openmc.Materials.from_xml(args.materials)

    # --- source: MATCH the forward per-coil runs (uniform emissivity DT source) ---
    src = openmc.CompiledSource(library=args.source_lib, parameters=args.source_params)
    settings = openmc.Settings()
    settings.source = src
    settings.run_mode = "fixed source"
    settings.particles = args.particles
    settings.batches = args.batches
    settings.seed = args.seed
    if args.weight_windows and os.path.exists(args.weight_windows):
        ww = openmc.hdf5_to_wws(args.weight_windows)
        settings.weight_windows = ww
        settings.weight_windows_on = True
    model.settings = settings

    # --- mesh flux tally on the ADJOINT mesh (co-located with psi_dagger) ---
    mesh = openmc.RegularMesh()
    mesh.lower_left = ll.tolist()
    mesh.upper_right = ur.tolist()
    mesh.dimension = dim
    t = openmc.Tally(name="fwd_meshflux")
    t.filters = [openmc.MeshFilter(mesh)]
    t.scores = ["flux"]
    model.tallies = openmc.Tallies([t])
    return model, dim


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adjoint", required=True, help="adjoint npz (for the mesh extent)")
    ap.add_argument("--dagmc", required=True, help="baseline dagmc_qh_step1_*.h5m")
    ap.add_argument("--materials", required=True, help="baseline materials.xml")
    ap.add_argument("--source-lib", required=True, help="libstellarator_source.so path")
    ap.add_argument("--source-params", default="emissivity=uniform",
                    help="StellaratorSource params (MATCH the forward per-coil runs)")
    ap.add_argument("--weight-windows", default=None, help="FW-CADIS weight_windows.h5")
    ap.add_argument("--particles", type=int, default=4_000_000)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="qh_fwd_meshflux.npz")
    args = ap.parse_args()

    model, dim = build(args)
    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="fwd_meshflux")
        flux = t.mean.ravel()            # C-order, matches psi.reshape(dim) in patch_worth
        err = t.std_dev.ravel()
    np.savez(args.out, flux=flux.reshape(dim), err=err.reshape(dim),
             dimension=np.array(dim))
    rel = err[flux > 0] / flux[flux > 0]
    print(f"wrote {args.out}  mesh={dim}  nonzero voxels={int((flux>0).sum())}  "
          f"median relerr={np.median(rel):.3f}  p90 relerr={np.percentile(rel,90):.3f}")


if __name__ == "__main__":
    main()
