#!/usr/bin/env python
"""recover_transport.py -- ISOLATED transport validation of the morphological offset.

Runs the sweep's own transport pipeline on a previously-'error' (DAGMC-leaky) device
but builds the conformal shells with stellarator_geometry_offset2.build_layers_sdf
(the non-inverting Minkowski-buffer offset) instead of the live normal offset. Writes
records to a PRIVATE out dir (default sweep_out_recover/) and NEVER touches the live
sweep dirs. Reuses sweep.run_transport / dagmc_writer / deltas_from_vals verbatim so
the numbers are directly comparable to a normal sweep 'done' record.

Success criterion: a device that errored under the normal offset now runs to a clean
'done' with lost_particles == 0 and finite delta_free_* (steering) values.

Run ON the Lima 'spf' worker (needs openmc+DAGMC + compiled source .so):
  limactl shell spf
  source spf_prototype/local_worker/activate_worker.sh
  cd spf_prototype/sweep
  OMP_NUM_THREADS=2 python -u recover_transport.py 1505944 2190977 2296767 \
      --out sweep_out_recover --particles 20000 --batches 5
"""
from __future__ import annotations
import argparse
import json
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))
import sweep
import quasr_loader as ql
import dagmc_writer
import stellarator_geometry_offset2 as sg2

DATADIR = HERE.parent / "data"


def recover_one(cid, out_dir, particles, batches, pol_res=None, force=False,
                field_res=None):
    blanket = "baseline"
    stem = f"quasr{cid}"
    tag = f"{cid}_{blanket}"
    rec_path = Path(out_dir) / f"config_{tag}.json"
    if rec_path.exists() and not force:
        print(f"[recover] {tag}: record exists, skipping (use --force)", flush=True)
        return json.loads(rec_path.read_text())

    rec = dict(id=cid, stem=stem, blanket=blanket, offset="morphological_buffer",
               status="started", warnings=[])
    try:
        dev = ql.load_device(cid, n_samples=256)
        rec.update(nfp=dev.meta["nfp"])
        field_native = sweep._make_field(dev, {}, {}, stem)
        layers = sweep.BLANKETS[blanket]

        # Field grid resolution: lost-particle counting is PURE GEOMETRY (a leaky
        # DAGMC loses particles regardless of the source direction), so this leak
        # test uses a coarse field map by default to make Biot-Savart staging fast.
        # delta_free is then INDICATIVE (coarse field); the leak verdict is exact.
        fr = field_res or (8, 12, 8)
        a_native = sweep.stage_surface_and_field(dev, field_native, stem,
                                                 nR=fr[0], nphi=fr[1], nZ=fr[2])
        rec["field_res"] = list(fr)
        scale = sweep.TARGET_A_CM / (a_native * sweep.UNIT_CM)
        rec["scale"] = float(scale)
        rec["built_minor_cm"] = float(a_native * sweep.UNIT_CM * scale)

        # ---- NON-INVERTING morphological offset (the fix under test) ----------
        geom_dir = DATADIR / f"{stem}_{blanket}_recover_geom"
        m_geo = sg2.build_layers_sdf(stem, layers=layers, scale=scale,
                                     outdir=str(geom_dir), pol_res=pol_res)
        bad = [L["name"] for L in m_geo["layers"]
               if not (L["watertight"] and L["simple_cross_section"])]
        if bad:
            rec["status"] = "geometry_folded"
            rec["warnings"].append(f"sdf offset non-simple layers: {bad}")
            sweep._write(rec_path, rec)
            return rec
        rec["pol_res"] = m_geo["pol_res"]

        # ---- DAGMC (pymoab-free) ---------------------------------------------
        h5m = str(geom_dir / f"{stem}.h5m")
        info = dagmc_writer.build_from_stls(str(geom_dir), h5m)
        errs = dagmc_writer.verify_structure(h5m)
        if errs:
            rec["status"] = "dagmc_invalid"
            rec["warnings"] += errs
            sweep._write(rec_path, rec)
            return rec
        rec["dagmc"] = dict(volumes=info["volumes"], tris=info["tris"])

        if not sweep.have_openmc():
            rec["status"] = "staged_no_transport"
            sweep._write(rec_path, rec)
            print(f"[recover] {cid}: staged (no openmc here). Geometry watertight, "
                  "shells non-crossing; run this on the Lima worker to transport.",
                  flush=True)
            return rec

        # ---- transport (free-stream, near-void) reusing the sweep verbatim ----
        results_dir = Path(out_dir) / f"transport_{tag}"
        vals, lost = sweep.run_transport(
            stem, scale, h5m, layers, a_native, results_dir,
            particles, batches, "free")
        rec["lost_particles"] = int(lost)
        rec.update(sweep.deltas_from_vals(vals))
        rec["status"] = "done"
        sweep._write(rec_path, rec)
        print(f"[recover] {cid}: DONE  lost={lost}  "
              f"delta_free_perp={rec.get('delta_free_perpendicular')}  "
              f"delta_free_par={rec.get('delta_free_parallel')}", flush=True)
        return rec
    except Exception as e:  # noqa
        rec["status"] = "error"
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc()
        sweep._write(rec_path, rec)
        print(f"[recover] {cid}: ERROR {e}", flush=True)
        return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("devices", nargs="+")
    ap.add_argument("--out", default="sweep_out_recover")
    ap.add_argument("--particles", type=int, default=20000)
    ap.add_argument("--batches", type=int, default=5)
    ap.add_argument("--pol-res", type=int, default=None)
    ap.add_argument("--field-res", type=int, nargs=3, default=None,
                    metavar=("NR", "NPHI", "NZ"),
                    help="Biot-Savart field grid (default 8 12 8; coarse=fast, "
                         "leak verdict is field-independent)")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    summary = []
    for cid in a.devices:
        r = recover_one(cid, a.out, a.particles, a.batches, a.pol_res, a.force,
                        tuple(a.field_res) if a.field_res else None)
        summary.append((cid, r.get("status"), r.get("lost_particles"),
                        r.get("delta_free_perpendicular"),
                        r.get("delta_free_parallel")))
    print("\n== recover_transport summary ==")
    print(f"{'device':10s} {'status':10s} {'lost':>5s} {'d_perp':>10s} {'d_par':>10s}")
    for cid, st, lost, dp, dq in summary:
        dp = f"{dp:.4f}" if isinstance(dp, float) else "-"
        dq = f"{dq:.4f}" if isinstance(dq, float) else "-"
        print(f"{cid:10s} {str(st):10s} {str(lost):>5s} {dp:>10s} {dq:>10s}")


if __name__ == "__main__":
    main()
