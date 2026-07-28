#!/usr/bin/env python
"""Combine the 20 per-coil adjoint importance maps (recip_w35f_coilNN_flat, coils 11..30) into a
single multi-coil importance psi_dagger_multi -- WITHOUT any new transport. All 20 maps already
exist and are co-located on the SAME mesh ([30,30,20], identical lower_left/upper_right), so the
combination is pure numpy post-processing.

MOTIVATION (kill-shot 9193680): attribution-guided placement using the single coil-20 adjoint
crushed coil 20 (4.7x) but the peak MIGRATED to coil 25 -- single-coil importance over-protects
its own target and ignores the coil that then sets the peak. A multi-coil importance protects the
whole coil set, so the placement fed by it (via the EXISTING build_step1 'placed' allocator and
patch_worth) addresses the migration directly. This is needed regardless of how A-vs-W resolves.

Combine modes (design objective -> weighting of the per-coil response R_k):
  sum         psi_multi = sum_k psi_k            total coil dose; CADIS response = sum over coils
  fluxweight  psi_multi = sum_k (R_k^p) psi_k    R_k = uniform coil-k fast flux; hot coils dominate
                                                 -> surrogate for reducing the PEAK (default, p=1;
                                                 raise p or iterate to concentrate on the hottest)
  max         psi_multi = max_k psi_k            voxel-wise worst coil; protect the envelope

Output: combined npz in the SAME format as a single-coil adjoint npz (importance/dimension/
lower_left/upper_right/coil_cells + carried 'flux'/'err'/'response'/'scatter_order' template),
so it drops straight into  patch_worth.py --adjoint  and  adjoint_placement.contributon().
Adds combine_mode / component_cells / weights for provenance.

Usage:
  python multicoil_importance.py --adjoint-dir DIR --out OUT.npz
        [--mode sum|fluxweight|max] [--flux-weights coil_run.npz] [--power P] [--coils 11-30]
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np

DEF_DIR = "/burg-archive/home/tjk2147/adjoint_test"


def _coil_num(path):
    m = re.search(r"coil(\d+)_flat", os.path.basename(path))
    return int(m.group(1)) if m else -1


def _parse_coils(spec):
    """'11-30' or '20,25,28' -> sorted list of ints."""
    out = set()
    for tok in str(spec).split(","):
        tok = tok.strip()
        if "-" in tok:
            a, b = tok.split("-"); out.update(range(int(a), int(b) + 1))
        elif tok:
            out.add(int(tok))
    return sorted(out)


def load_maps(adjoint_dir, coils):
    """Load per-coil importance maps; assert identical mesh; return (cells, maps, template npz)."""
    dirs = sorted(glob.glob(os.path.join(adjoint_dir, "recip_w35f_coil*_flat")), key=_coil_num)
    cells, maps, template = [], [], None
    ref = None
    for d in dirs:
        k = _coil_num(d)
        if k not in coils:
            continue
        f = os.path.join(d, "adjoint_importance_flat_P0.npz")
        if not os.path.exists(f):
            raise SystemExit(f"missing importance map: {f}")
        z = np.load(f)
        dim = tuple(int(x) for x in z["dimension"])
        ll = tuple(round(float(x), 6) for x in z["lower_left"])
        ur = tuple(round(float(x), 6) for x in z["upper_right"])
        key = (dim, ll, ur)
        if ref is None:
            ref, template = key, dict(z)  # first map is the format template
        elif key != ref:
            raise SystemExit(f"MESH MISMATCH coil {k}: {key} != {ref} -- maps are not co-located")
        cc = [int(x) for x in z["coil_cells"]] if "coil_cells" in z else [k]
        assert cc == [k], f"coil dir {d} has coil_cells={cc}, expected [{k}]"
        cells.append(k)
        maps.append(z["importance"].astype(float))
    if not cells:
        raise SystemExit(f"no coil maps found in {adjoint_dir} for coils {coils}")
    return cells, np.stack(maps, axis=0), template


def flux_weights(cells, coil_run_npz, power=1.0):
    """R_k^power for each coil cell from a coil_run_v3 output npz (uniform baseline recommended)."""
    z = np.load(coil_run_npz)
    cc = [int(x) for x in z["coil_cells"]]
    flux = np.asarray(z["coil_fast_flux"], float)
    lut = {c: flux[i] for i, c in enumerate(cc)}
    w = np.array([lut[k] for k in cells], float)
    return w ** power, w


def combine(cells, maps, mode, weights=None):
    if mode == "sum":
        w = np.ones(len(cells))
        psi = (w[:, None, None, None] * maps).sum(axis=0)
    elif mode == "fluxweight":
        if weights is None:
            raise SystemExit("mode=fluxweight needs --flux-weights coil_run.npz")
        w = weights
        psi = (w[:, None, None, None] * maps).sum(axis=0)
    elif mode == "max":
        w = np.ones(len(cells))
        psi = maps.max(axis=0)
    else:
        raise SystemExit(f"unknown mode {mode}")
    return psi, w


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adjoint-dir", default=DEF_DIR)
    ap.add_argument("--coils", default="11-30", help="coil cells to combine, e.g. 11-30 or 20,25")
    ap.add_argument("--mode", default="fluxweight", choices=["sum", "fluxweight", "max"])
    ap.add_argument("--flux-weights", default=None, help="coil_run npz for fluxweight R_k")
    ap.add_argument("--power", type=float, default=1.0, help="R_k exponent (fluxweight); >1 = peakier")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    coils = _parse_coils(args.coils)
    cells, maps, template = load_maps(args.adjoint_dir, coils)
    print(f"loaded {len(cells)} coil maps: cells={cells}  mesh={maps.shape[1:]}")

    weights = None; Rk = None
    if args.mode == "fluxweight":
        weights, Rk = flux_weights(cells, args.flux_weights, power=args.power)

    psi, w = combine(cells, maps, args.mode, weights)

    out = dict(template)                         # inherit format (dimension/ll/ur/flux/err/...)
    out["importance"] = psi
    out["coil_cells"] = np.array(cells)          # all combined cells
    out["combine_mode"] = args.mode
    out["component_cells"] = np.array(cells)
    out["weights"] = np.asarray(w, float)
    if Rk is not None:
        out["R_k_flux"] = Rk
        out["power"] = args.power
    np.savez(args.out, **out)

    # provenance: which coils dominate the combined importance mass
    mass = np.array([float((wk * maps[i]).sum()) for i, wk in enumerate(w)])
    order = np.argsort(mass)[::-1]
    tot = mass.sum()
    print(f"mode={args.mode}  nonzero voxels={int((psi>0).sum())}  max={psi.max():.3e}")
    print("top coils by importance mass:")
    for idx in order[:6]:
        print(f"  coil {cells[idx]:2d}: weight={w[idx]:.3e}  mass_frac={mass[idx]/tot:.3f}")
    print(f"WROTE {args.out}")


if __name__ == "__main__":
    main()
