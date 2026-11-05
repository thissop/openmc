#!/usr/bin/env python
"""Phase-1 head-to-head: grid (Thaddaeus) vs spectral (Ethan) flux-surface representations,
on a real VMEC equilibrium, with NO compiled OpenMC required.

For a given wout, we build the exact SpectralRep, then build GridRep at a sweep of grid
resolutions and report -- per resolution -- the accuracy (grid-vs-spectral RMS/max relative
deviation in R, Z, sqrt(g)), the evaluation throughput of each representation, and the stored
memory footprint. This quantifies the tradeoff the two sources make: the grid is O(1) per draw
but pays interpolation error + memory; the spectral form is exact + compact but costs a Fourier
sum per draw.

Usage:
    python compare_representations.py <wout.nc> [--out results.json] [--nquery 20000]

Phase 2 (compare_samplers.py, needs both builds) draws actual births from each compiled source
and KS-tests the (R, Z, phi) distributions across the QUASR QA/QH set.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np

from fluxrep import read_wout, SpectralRep, GridRep
import metrics


# grid resolutions to sweep (nrho, ntheta, nzeta); the last matches Thaddaeus' production fluxmap
RESOLUTIONS = [(8, 32, 96), (12, 48, 144), (16, 64, 192), (24, 96, 288)]


def run(wout_path, nquery=20000, resolutions=RESOLUTIONS):
    wout = read_wout(wout_path)
    nfp = wout["nfp"]
    spec = SpectralRep(wout)
    spec_speed = metrics.throughput(spec, n=20000, nfp=nfp)
    spec_bhat_speed = (metrics.bhat_throughput(spec, n=20000, nfp=nfp)
                       if spec.has_field else None)
    spec_mem = metrics.footprint(spec)

    rows = []
    for (nr, nt, nz) in resolutions:
        grid = GridRep(spec, nrho=nr, ntheta=nt, nzeta=nz, nfp=nfp)
        acc = metrics.accuracy(spec, grid, n=nquery, nfp=nfp)
        spd = metrics.throughput(grid, n=100000, nfp=nfp)
        bhat_spd = (metrics.bhat_throughput(grid, n=100000, nfp=nfp)
                    if grid.has_field else None)
        mem = metrics.footprint(grid)
        rows.append({
            "resolution": [nr, nt, nz],
            "accuracy": acc,
            "grid_points_per_sec": spd["points_per_sec"],
            "grid_bhat_points_per_sec": (
                bhat_spd["points_per_sec"] if bhat_spd else None),
            "grid_MiB": mem["MiB"],
        })
    return {
        "wout": str(wout_path),
        "nfp": nfp,
        "ns_surfaces": wout["ns"],
        "n_modes": int(len(wout["xm"])),
        "spectral": {
            "points_per_sec": spec_speed["points_per_sec"],
            "bhat_points_per_sec": (
                spec_bhat_speed["points_per_sec"] if spec_bhat_speed else None),
            "MiB": spec_mem["MiB"],
            "note": "exact to VMEC truncation; coefficient storage only",
        },
        "grid_sweep": rows,
    }


def print_table(res):
    s = res["spectral"]
    print(f"\nDevice: {res['wout']}  (nfp={res['nfp']}, ns={res['ns_surfaces']}, "
          f"{res['n_modes']} Fourier modes)")
    print(f"SPECTRAL (Ethan):  {s['points_per_sec']:,.0f} pts/s   {s['MiB']:.3f} MiB   exact")
    print(f"{'grid (nrho,nt,nz)':>20} {'R rms':>9} {'sqrtg rms':>10} "
          f"{'source TV':>9} {'B_w p99':>9} {'pts/s':>12} {'MiB':>8}")
    for r in res["grid_sweep"]:
        a = r["accuracy"]
        print(f"{str(tuple(r['resolution'])):>20} {a['R']['rms']:>9.2e} "
              f"{a['sqrtg']['rms']:>10.2e} "
              f"{a['source_pdf']['total_variation']:>9.3g} "
              f"{a.get('bhat_degrees', {}).get('source_weighted_p99', float('nan')):>9.3g} "
              f"{r['grid_points_per_sec']:>12,.0f} {r['grid_MiB']:>8.3f}")
    print("\nRead: grid buys throughput and pays it back in interpolation error + memory as "
          "resolution rises; spectral is exact at fixed tiny memory.\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("wout")
    ap.add_argument("--out", default=None)
    ap.add_argument("--nquery", type=int, default=20000)
    args = ap.parse_args(argv)
    res = run(args.wout, nquery=args.nquery)
    print_table(res)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=2)
        print(f"wrote {args.out}")
    return res


if __name__ == "__main__":
    main()
