#!/usr/bin/env python
"""Reactor-relevance filter, step 1: min coil-to-plasma-surface standoff per QUASR device.

Ethan (2026-07-13): scrape QUASR for *reactor-relevant* devices BEFORE building any
magnet geometry -- filter by the minimum coil-to-plasma-surface distance. A device whose
modular coils hug the plasma cannot fit a ~1 m blanket+shield stack when scaled to reactor
size, so it is not reactor-relevant no matter how good its quasisymmetry.

For each device we load the real coil filaments (sweep/quasr_loader.load_coils) and the
LCFS boundary grid (quasr_geom via load_device), then compute the minimum Euclidean
distance from any coil filament point to any LCFS point. We report it three ways:
  d_min        raw, in QUASR native length units (R0 ~ O(1))
  d_min / a    normalized by the plasma minor radius  <- the scale-free reactor-room metric
  d_min / R0   normalized by the major radius
The blanket+shield stack to the coil is ~1.01 m (sol5+W0.2+steel3.8+Be2+FLiBe50+shield40),
so "reactor-relevant" means d_min/a large enough that, scaled to a reactor minor radius,
the physical gap clears that stack. The threshold itself is a tunable (see --min_gap_over_a).

Offline: uses only cached serials in data/quasr/serials/ and boundaries in data/quasr/nml/.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                       # sweep/ (quasr_loader)
sys.path.insert(0, str(HERE.parent / "python"))     # python/ (quasr_geom)
import quasr_loader as ql                 # noqa: E402

DATADIR = HERE.parent / "data" / "quasr"
SERIAL_DIR = DATADIR / "serials"


def _min_point_to_set(P, Q, chunk=2048):
    """min over p in P of min over q in Q of |p-q|, and the index in P achieving it.
    Chunked so we never materialize the full |P|x|Q| matrix. Uses cKDTree if available."""
    try:
        from scipy.spatial import cKDTree
        d, _ = cKDTree(Q).query(P, k=1)
        return float(d.min()), int(d.argmin()), d
    except Exception:
        dmin_per_p = np.empty(P.shape[0])
        for i in range(0, P.shape[0], chunk):
            blk = P[i:i + chunk]
            d2 = ((blk[:, None, :] - Q[None, :, :]) ** 2).sum(-1)
            dmin_per_p[i:i + chunk] = np.sqrt(d2.min(1))
        return float(dmin_per_p.min()), int(dmin_per_p.argmin()), dmin_per_p


def coil_standoff(ID, n_samples=256, ntheta=64, nphi=192):
    """Return a dict of standoff metrics for one device (loads coils + LCFS)."""
    d = ql.load_device(ID, n_samples=n_samples, coils=True)
    coil_pts = np.vstack([poly for poly, _ in d.coils])          # (Ncoilpts, 3)
    bg = d.boundary_grid(ntheta=ntheta, nphi=nphi, rho=1.0)
    lcfs_pts = bg["xyz"].reshape(-1, 3)                           # (Nlcfs, 3)

    d_min, imin, dper = _min_point_to_set(coil_pts, lcfs_pts)
    # per-coil closest approach (useful: which coil is tightest)
    per_coil = []
    off = 0
    for poly, cur in d.coils:
        n = poly.shape[0]
        per_coil.append(float(dper[off:off + n].min()))
        off += n

    a = d.meta.get("minor_radius") or np.nan
    R0 = d.meta.get("R0") or np.nan
    return dict(
        ID=int(ID), nfp=d.meta.get("nfp"), symmetry_class=d.meta.get("symmetry_class"),
        aspect=d.meta.get("aspect"), a=a, R0=R0,
        n_coils=len(d.coils), nc_per_hp=d.meta.get("nc_per_hp"),
        qs_error=d.meta.get("qs_error"),
        d_min=d_min,
        d_min_over_a=(d_min / a) if a and a == a else np.nan,
        d_min_over_R0=(d_min / R0) if R0 and R0 == R0 else np.nan,
        d_min_per_coil_min=min(per_coil), d_min_per_coil_max=max(per_coil),
    )


def _cached_ids():
    ids = []
    for fn in sorted(SERIAL_DIR.glob("serial*.json")):
        try:
            ids.append(int(fn.stem.replace("serial", "")))
        except ValueError:
            pass
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*", type=int,
                    help="device IDs; default = all cached serials")
    ap.add_argument("--out", default=str(DATADIR / "coil_standoff.csv"))
    ap.add_argument("--from_csv", default=None,
                    help="skip recompute; re-summarize/filter an existing coil_standoff.csv")
    ap.add_argument("--reactor_R0_cm", type=float, default=1000.0,
                    help="target reactor major radius (cm) for the absolute-gap filter")
    ap.add_argument("--blanket_cm", type=float, default=101.0,
                    help="required physical coil-plasma gap (cm): blanket+shield stack to the coil")
    args = ap.parse_args()

    if args.from_csv:
        rows = list(csv.DictReader(open(args.from_csv)))
        for r in rows:  # cast the numeric columns we use
            for k in ("d_min", "d_min_over_a", "d_min_over_R0", "aspect"):
                r[k] = float(r[k]) if r.get(k) not in (None, "", "nan", "None") else np.nan
        fails = []
    else:
        ids = args.ids or _cached_ids()
        rows, fails = [], []
        for k, ID in enumerate(ids):
            try:
                r = coil_standoff(ID)
                rows.append(r)
                print(f"[{k+1}/{len(ids)}] {ID}: d_min={r['d_min']:.4f}  "
                      f"d_min/R0={r['d_min_over_R0']:.3f}  n_coils={r['n_coils']}  "
                      f"class={r['symmetry_class']} nfp={r['nfp']} aspect={r['aspect']}")
            except Exception as e:
                fails.append((ID, repr(e)))
                print(f"[{k+1}/{len(ids)}] {ID}: FAIL {e!r}")
        if rows:
            with open(args.out, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print(f"\nwrote {args.out}  ({len(rows)} devices, {len(fails)} failed)")

    # ---- reactor-relevance by ABSOLUTE physical gap at a target reactor scale ----
    # gap_cm = (d_min/R0) * reactor_R0_cm. This is aspect-independent, unlike d_min/a.
    for r in rows:
        r["gap_cm"] = (float(r["d_min_over_R0"]) * args.reactor_R0_cm
                       if np.isfinite(r["d_min_over_R0"]) else np.nan)
    dovr = np.array([r["d_min_over_R0"] for r in rows], float); dovr = dovr[np.isfinite(dovr)]
    gap = np.array([r["gap_cm"] for r in rows], float); gap = gap[np.isfinite(gap)]
    print(f"\nd_min/R0 (aspect-independent): min {dovr.min():.3f}  median {np.median(dovr):.3f}  "
          f"max {dovr.max():.3f}")
    print(f"physical gap at R0={args.reactor_R0_cm/100:.0f} m (cm): "
          f"min {gap.min():.0f}  median {np.median(gap):.0f}  max {gap.max():.0f}")
    keep = [r for r in rows if np.isfinite(r["gap_cm"]) and r["gap_cm"] >= args.blanket_cm]
    keep.sort(key=lambda r: -r["gap_cm"])
    print(f"\n{len(keep)}/{len(rows)} devices clear a {args.blanket_cm:.0f} cm blanket+shield "
          f"stack at R0={args.reactor_R0_cm/100:.0f} m (reactor-relevant). Top 15 by gap:")
    for r in keep[:15]:
        print(f"  {r['ID']}  gap={r['gap_cm']:.0f} cm  d_min/R0={float(r['d_min_over_R0']):.3f}  "
              f"class={r['symmetry_class']} nfp={r['nfp']} aspect={float(r['aspect']):.1f}")
    if fails:
        print(f"\n{len(fails)} failures (first 5): {fails[:5]}")


if __name__ == "__main__":
    main()
