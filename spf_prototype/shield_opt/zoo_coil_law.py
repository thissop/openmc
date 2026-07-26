#!/usr/bin/env python
"""Step 2 (SCIENCE_STRATEGY.md sec 5): zoo-wide FREE-STREAMING coil-load
concentration -> geometry law.

For every cached QUASR device we compute a geometric, transport-free neutron load
on the real coil filaments from the (1-rho^2)-weighted plasma fusion source, then
concentration metrics (Gini / participation ratio, concentration.py) and correlate
them against coil + plasma geometry descriptors with Chatterjee xi
(descriptor_correlation.py). Cheap, no DAGMC, no OpenMC transport.

Coil-side free-streaming load (convex free-streaming limit, no occlusion / no
scattering -- matching the first-wall culprit-map caveats in geometry_peaking.py):
    phi(p) = sum_x  s(x) / |p - x|^2          (isotropic point emission, 4pi dropped)
for each coil-filament point p, source point x with (1-rho^2) emissivity weight s(x).
This is exactly the solid-angle / line-of-sight neutron load on the conductor.

Two concentration views (both reported; they answer different questions):
  * POINT-pooled (primary): Gini/PR over phi at ALL filament points across all coils.
    Non-degenerate; measures the spatial hot-spot structure of the coil-side load
    (the reactor-relevant "are the leak paths few and targetable" question).
  * PER-COIL (strategy-literal): Gini/PR over the per-coil MEAN load. NOTE: coils in a
    stellarator-symmetry orbit see identical load, so this collapses to the nc_per_hp
    distinct base coils and is exactly 0 for nc_per_hp=1 devices. Reported with caveat.

No transport. Labelled a geometric proxy throughout.
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "sweep"))
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(HERE))

import quasr_loader as ql          # noqa: E402
import concentration as conc       # noqa: E402

DATA = ROOT / "data" / "quasr"
SERIAL_DIR = DATA / "serials"

# Free-streaming resolution. Source: (1-rho^2)-weighted plasma volume; coil: filament
# samples per base curve. Modest but convergence-checked on a handful of devices.
SRC_RES = dict(n_rho=6, n_theta=16, n_phi=48, rho_max=0.85)
COIL_NSAMP = 160


def cached_ids():
    return sorted(int(p.stem.replace("serial", "")) for p in SERIAL_DIR.glob("serial*.json"))


def load_catalogue_index():
    idx = DATA / "_catalogue_index.json"
    if idx.exists():
        return {int(k): v for k, v in json.loads(idx.read_text()).items()}
    ids = set(cached_ids())
    cat = {}
    import csv
    with gzip.open(DATA / "catalogue.csv.gz", "rt") as f:
        for row in csv.DictReader(f):
            try:
                I = int(float(row["ID"]))
            except (TypeError, ValueError):
                continue
            if I in ids:
                cat[I] = row
    idx.write_text(json.dumps({str(i): cat[i] for i in cat}))
    return cat


def load_standoff():
    import csv
    out = {}
    with open(DATA / "coil_standoff.csv") as f:
        for r in csv.DictReader(f):
            out[int(r["ID"])] = r
    return out


# --------------------------------------------------------------------------- #
# Geometry descriptors from filaments
# --------------------------------------------------------------------------- #
def coil_nonplanarity(poly):
    """Out-of-plane fraction of a closed coil polyline via SVD of centered points.
    s0>=s1>=s2 singular values; np = s2 / sqrt(s0^2 + s1^2)  (0 = planar coil)."""
    P = poly - poly.mean(axis=0, keepdims=True)
    s = np.linalg.svd(P, compute_uv=False)
    return float(s[2] / np.sqrt(s[0] ** 2 + s[1] ** 2 + 1e-30))


# --------------------------------------------------------------------------- #
# Free-streaming per-coil / per-point flux
# --------------------------------------------------------------------------- #
def coil_free_stream(dev, src_res=None, coil_nsamp=COIL_NSAMP, block=4096, eps=1e-30):
    """Return per-point flux, per-point coil-id, per-coil mean/peak, non-planarity.

    phi(p) = sum_x w_x / |p-x|^2 over the (1-rho^2) plasma source. No occlusion."""
    src_res = src_res or SRC_RES
    xyz, w, _ = dev.source_sample(**src_res)                      # (Nsrc,3),(Nsrc,)
    polys = [poly for poly, _ in dev.coils]
    pts = np.vstack(polys)                                        # (Npts,3)
    coil_id = np.concatenate([np.full(p.shape[0], k) for k, p in enumerate(polys)])
    ncoil = len(polys)

    phi = np.empty(pts.shape[0])
    for i0 in range(0, pts.shape[0], block):
        i1 = min(i0 + block, pts.shape[0])
        d = pts[i0:i1, None, :] - xyz[None, :, :]                # (m,Nsrc,3)
        r2 = np.sum(d * d, axis=2)                               # (m,Nsrc)
        phi[i0:i1] = np.sum(w[None, :] / (r2 + eps), axis=1)

    per_coil_mean = np.array([phi[coil_id == k].mean() for k in range(ncoil)])
    per_coil_peak = np.array([phi[coil_id == k].max() for k in range(ncoil)])
    nonplanar = np.array([coil_nonplanarity(p) for p in polys])
    return dict(phi=phi, coil_id=coil_id, ncoil=ncoil,
                per_coil_mean=per_coil_mean, per_coil_peak=per_coil_peak,
                nonplanarity=float(nonplanar.mean()),
                nonplanarity_max=float(nonplanar.max()))


# --------------------------------------------------------------------------- #
# Per-device record
# --------------------------------------------------------------------------- #
def analyze_device(ID, cat, standoff, src_res=None):
    dev = ql.load_device(ID, n_samples=COIL_NSAMP)
    fs = coil_free_stream(dev, src_res=src_res)

    phi = fs["phi"]
    # POINT-pooled concentration (primary, non-degenerate).
    pt = conc.concentration_report(phi)
    # PER-COIL concentration (strategy-literal; degenerate for nc_per_hp=1).
    pc_mean = fs["per_coil_mean"]
    coil_gini = conc.gini(pc_mean)
    coil_pr, coil_pr_over_n = conc.participation_ratio(pc_mean)

    row = cat.get(ID, {})
    so = standoff.get(ID, {})

    def cf(d, k, default=np.nan):
        v = d.get(k)
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    hel = cf(row, "helicity")
    qs_type = 0.0 if (np.isfinite(hel) and round(hel) == 0) else 1.0  # QA=0, QH=1

    return dict(
        ID=int(ID),
        nfp=cf(row, "nfp"),
        iota=cf(row, "mean_iota"),
        aspect=cf(row, "aspect_ratio"),
        elongation=cf(row, "max_elongation"),
        mean_elongation=cf(row, "mean_elongation"),
        curvature=cf(row, "max_kappa"),
        msc=cf(row, "max_msc"),
        qs_type=qs_type,
        nc_per_hp=cf(row, "nc_per_hp"),
        n_coils=fs["ncoil"],
        # standoff (scale-free + reactor-cm proxy)
        d_min_over_a=cf(so, "d_min_over_a"),
        d_min_over_R0=cf(so, "d_min_over_R0"),
        # coil geometry
        nonplanarity=fs["nonplanarity"],
        nonplanarity_max=fs["nonplanarity_max"],
        # ---- responses ----
        point_gini=pt["gini"],
        point_pr_over_n=pt["pr_over_n"],
        coil_peak=pt["peaking"],           # max/mean over all filament points
        coil_gini=coil_gini,               # across per-coil means (degenerate if nc_per_hp=1)
        coil_pr_over_n=coil_pr_over_n,
    )


def run(ids=None, src_res=None, verbose=True):
    cat = load_catalogue_index()
    standoff = load_standoff()
    ids = ids or cached_ids()
    recs, fails = [], []
    t0 = time.time()
    for k, ID in enumerate(ids):
        try:
            recs.append(analyze_device(ID, cat, standoff, src_res=src_res))
        except Exception as e:  # noqa: BLE001
            fails.append((ID, f"{type(e).__name__}: {e}"))
            if verbose:
                print(f"[skip] {ID}: {type(e).__name__}: {e}")
            continue
        if verbose and (k % 25 == 0 or k == len(ids) - 1):
            print(f"[{k+1}/{len(ids)}] {ID} pt_gini={recs[-1]['point_gini']:.3f} "
                  f"coil_peak={recs[-1]['coil_peak']:.2f} "
                  f"({time.time()-t0:.0f}s)")
    return recs, fails


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None, help="limit to first n devices (debug)")
    ap.add_argument("--out", default=str(HERE / "data" / "zoo_coil_law.csv"))
    args = ap.parse_args()
    ids = cached_ids()
    if args.n:
        ids = ids[:args.n]
    recs, fails = run(ids)
    import csv
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(recs[0].keys()))
        wr.writeheader()
        wr.writerows(recs)
    print(f"\nwrote {outp}  ({len(recs)} devices, {len(fails)} failed)")
    if fails:
        print("failures (first 8):", fails[:8])
