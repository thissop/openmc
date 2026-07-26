#!/usr/bin/env python
"""Leak-recovery probe for ERROR (DAGMC-leaky) devices.

For a given device id, rebuild the conformal shell at a chosen surface-mesh
resolution (ntheta x nph -> STL facet density) and run the CHEAP watertightness
probe (magnet_suite/dagmc_watertight_check.py: 1 batch, isotropic 14 MeV point
source, low-density filler) to measure the lost-particle fraction. Compare the
baseline resolution (64x96, what the sweep used) against a higher one.

Honest verdict: leak FIXABLE if lost-fraction drops to ~0 at higher res;
INTRINSIC if it persists (concave-inboard offset self-near-intersection).

Usage: python recover_leaky.py <id> --nt 128 --nph 192 [--probe-particles 20000]
Writes data/<stem>_baseline_geom/<stem>.h5m and a small JSON verdict to stdout.
Does NOT touch sweep records. Reuses the exact pipeline modules.
"""
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
import sweep, quasr_loader as ql, stellarator_geometry as sg, dagmc_writer
import biotsavart_field as bsf

DATADIR = HERE.parent / "data"
UNIT_CM = 100.0
PROBE = HERE.parent / "shield_opt/magnet_suite/dagmc_watertight_check.py"


def stage_surface_hires(dev, stem, ntheta, nph):
    """Re-write data/<stem>_surface.npz at higher poloidal/toroidal resolution."""
    th = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2*np.pi, nph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z = dev.device.RZ(TH, PH, 1.0)
    Rc, Zc = R*UNIT_CM, Z*UNIT_CM
    np.savez(DATADIR / f"{stem}_surface.npz", R=Rc, Z=Zc, nfp=np.int64(dev.meta["nfp"]),
             theta=th, phi_axis=ph, phi=np.broadcast_to(ph[None,:], Rc.shape).copy())
    return 0.5*(R.max()-R.min())


def build_and_probe(cid, nt, nph, probe_particles, max_lost, rel_max_lost, scale_override=None):
    stem = f"quasr{cid}"
    dev = ql.load_device(cid, n_samples=256)
    a_native = stage_surface_hires(dev, stem, nt, nph)
    scale = scale_override or (sweep.TARGET_A_CM / (a_native*UNIT_CM))
    layers = sweep.BLANKETS["baseline"]
    geom_dir = DATADIR / f"{stem}_baseline_geom"
    # scale-retry like the sweep
    m_geo, used_scale, folds = None, scale, []
    for _ in range(4):
        try:
            m_geo = sg.build_layers(stem, layers=layers, scale=used_scale, outdir=geom_dir)
            bad = [L["name"] for L in m_geo["layers"]
                   if not (L["watertight"] and L["simple_cross_section"])]
            if bad:
                folds.append((used_scale, bad)); m_geo=None; used_scale*=1.4; continue
            break
        except Exception as e:
            folds.append((used_scale, str(e))); m_geo=None; used_scale*=1.4
    if m_geo is None:
        return dict(id=cid, res=[nt,nph], status="geometry_folded", folds=folds)
    h5m = str(geom_dir / f"{stem}.h5m")
    info = dagmc_writer.build_from_stls(str(geom_dir), h5m)
    tris = info["tris"]
    # cheap watertight probe as subprocess (contained abort)
    out_json = geom_dir / f"{stem}_wt_{nt}x{nph}.json"
    cmd = [sys.executable, str(PROBE), h5m, "--particles", str(probe_particles),
           "--max-lost", str(max_lost), "--rel-max-lost", str(rel_max_lost),
           "--out", str(out_json)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    verdict = {}
    if out_json.exists():
        verdict = json.load(open(out_json))
    return dict(id=cid, res=[nt,nph], scale=round(used_scale,3), tris=tris,
                probe_rc=p.returncode, probe=verdict,
                stderr_tail=p.stderr.strip().splitlines()[-3:] if p.stderr else [])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("id", type=int)
    ap.add_argument("--nt", type=int, default=128)
    ap.add_argument("--nph", type=int, default=192)
    ap.add_argument("--probe-particles", type=int, default=20000)
    ap.add_argument("--max-lost", type=int, default=200)
    ap.add_argument("--rel-max-lost", type=float, default=1e-2)
    ap.add_argument("--baseline-too", action="store_true",
                    help="also probe at the sweep's 64x96 resolution for A/B comparison")
    a = ap.parse_args()
    res = []
    if a.baseline_too:
        res.append(build_and_probe(a.id, 64, 96, a.probe_particles, a.max_lost, a.rel_max_lost))
    res.append(build_and_probe(a.id, a.nt, a.nph, a.probe_particles, a.max_lost, a.rel_max_lost))
    print(json.dumps(res, indent=2))
