#!/usr/bin/env python
"""Multi-device aggregate of the grid-vs-spectral representation comparison across a set of
QUASR QA/QH VMEC wout files. Reports MEDIAN and IQR (25th-75th percentile) per grid resolution
for accuracy (grid-vs-spectral RMS deviation in R, Z, sqrt(g)), grid throughput, and grid memory
-- plus the spectral baseline. Median/IQR chosen over mean/std because the sqrt(g) interpolation
error is axis-dominated and skewed across devices.

Usage:
    python compare_multidevice.py wout_A.nc wout_B.nc ...  [--out agg.json] [--nquery 20000]
    python compare_multidevice.py --manifest wouts.txt      # one wout path per line
"""
from __future__ import annotations
import argparse, glob, json
import numpy as np

from compare_representations import run, RESOLUTIONS


def _iqr(x):
    x = np.asarray(x, float)
    q25, q50, q75 = np.percentile(x, [25, 50, 75])
    return {"median": float(q50), "iqr": [float(q25), float(q75)], "n": int(x.size)}


def aggregate(wout_paths, nquery=20000, resolutions=RESOLUTIONS):
    per_device = []
    for p in wout_paths:
        try:
            per_device.append(run(p, nquery=nquery, resolutions=resolutions))
        except Exception as e:                       # keep going; report failures honestly
            per_device.append({"wout": str(p), "error": repr(e)})
    ok = [d for d in per_device if "error" not in d]
    fail = [d for d in per_device if "error" in d]

    # spectral baseline aggregated across devices
    spec = {
        "points_per_sec": _iqr([d["spectral"]["points_per_sec"] for d in ok]),
        "MiB": _iqr([d["spectral"]["MiB"] for d in ok]),
    }
    # per-resolution aggregation. Use whichever resolution list the successful
    # devices were actually run on, so a custom --resolutions stays consistent.
    resolutions = [tuple(r["resolution"]) for r in ok[0]["grid_sweep"]] if ok \
        else RESOLUTIONS
    res_rows = []
    for i, resn in enumerate(resolutions):
        R_rms = [d["grid_sweep"][i]["accuracy"]["R"]["rms"] for d in ok]
        Z_rms = [d["grid_sweep"][i]["accuracy"]["Z"]["rms"] for d in ok]
        G_rms = [d["grid_sweep"][i]["accuracy"]["sqrtg"]["rms"] for d in ok]
        G_max = [d["grid_sweep"][i]["accuracy"]["sqrtg"]["max"] for d in ok]
        source_tv = [d["grid_sweep"][i]["accuracy"]["source_pdf"]["total_variation"]
                     for d in ok]
        B_p99 = [d["grid_sweep"][i]["accuracy"]["bhat_degrees"]["source_weighted_p99"]
                 for d in ok]
        pps = [d["grid_sweep"][i]["grid_points_per_sec"] for d in ok]
        mib = [d["grid_sweep"][i]["grid_MiB"] for d in ok]
        res_rows.append({
            "resolution": list(resn),
            "R_rms": _iqr(R_rms), "Z_rms": _iqr(Z_rms),
            "sqrtg_rms": _iqr(G_rms), "sqrtg_max": _iqr(G_max),
            "source_total_variation": _iqr(source_tv),
            "source_weighted_bhat_p99_degrees": _iqr(B_p99),
            "grid_points_per_sec": _iqr(pps), "grid_MiB": _iqr(mib),
        })
    return {"n_devices": len(ok), "n_failed": len(fail),
            "failed": [d["wout"] for d in fail],
            "spectral": spec, "grid_sweep": res_rows, "per_device": per_device}


def print_table(agg):
    print(f"\n=== grid vs spectral across {agg['n_devices']} QUASR devices "
          f"({agg['n_failed']} failed) ===")
    s = agg["spectral"]
    print(f"SPECTRAL: {s['points_per_sec']['median']:,.0f} pts/s "
          f"[{s['points_per_sec']['iqr'][0]:,.0f}-{s['points_per_sec']['iqr'][1]:,.0f}]   "
          f"{s['MiB']['median']:.3f} MiB   exact\n")
    hdr = (f"{'grid':>16} {'R rms (med)':>14} {'sqrtg rms (med)':>18} "
           f"{'source TV':>10} {'B_w p99':>10} {'MiB (med)':>10}")
    print(hdr)
    for r in agg["grid_sweep"]:
        print(f"{str(tuple(r['resolution'])):>16} "
              f"{r['R_rms']['median']:>14.2e} {r['sqrtg_rms']['median']:>18.2e} "
              f"{r['source_total_variation']['median']:>10.3g} "
              f"{r['source_weighted_bhat_p99_degrees']['median']:>10.3g} "
              f"{r['grid_MiB']['median']:>10.3f}")
    print("\n(medians shown; full IQR + per-device in the JSON)\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("wouts", nargs="*")
    ap.add_argument("--manifest", default=None, help="file with one wout path per line")
    ap.add_argument("--out", default="agg_representations.json")
    ap.add_argument("--nquery", type=int, default=20000)
    ap.add_argument("--resolutions", default=None,
                    help="semicolon-separated nr,nt,nz triples, "
                         "e.g. '24,96,288;32,128,384;48,192,576'")
    args = ap.parse_args(argv)
    if args.resolutions:
        resolutions = [tuple(int(v) for v in tri.split(","))
                       for tri in args.resolutions.split(";")]
    else:
        resolutions = RESOLUTIONS
    paths = list(args.wouts)
    if args.manifest:
        paths += [ln.strip() for ln in open(args.manifest) if ln.strip() and not ln.startswith("#")]
    paths = [p for pat in paths for p in (glob.glob(pat) or [pat])]
    if not paths:
        ap.error("no wout paths given")
    agg = aggregate(paths, nquery=args.nquery, resolutions=resolutions)
    print_table(agg)
    with open(args.out, "w") as f:
        json.dump(agg, f, indent=2)
    print(f"wrote {args.out}")
    return agg


if __name__ == "__main__":
    main()
