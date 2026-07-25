#!/usr/bin/env python
"""M7 payoff: does direction-resolved (P1) adjoint change the REAL coil importance map?

Compares the scalar (P0) and angular (P1) adjoint magnet-importance maps for the same
coil on the real DAGMC geometry (openmc-m7 build). M6 predicted, on a controlled slab,
that P1/P0 importance grows with penetration depth into scattering shield (up to ~1.7x
at 5 mfp) but is ~flat in clean gaps. Here we test it on the coil:

  * P1/P0 importance ratio overall, and vs distance from the coil (does it grow deep?)
  * contributon (S*psi_dagger) centroid direction for both -- does P1 re-point it?
  * placement target (theta0, phi0) for both -- does angular MOVE where shield goes?

Run:  python compare_p0_p1.py --p0 <P0.npz> --p1 <P1.npz> --fluxmap <qh_fluxmap.npz> \
          --coil-centroid X Y Z
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # shield_opt/
import adjoint_placement as ap  # noqa: E402


def _centers(d):
    ll, ur, dim = d["lower_left"], d["upper_right"], d["dimension"]
    return ap._voxel_centers(ll, ur, dim)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--p0", required=True)
    p.add_argument("--p1", required=True)
    p.add_argument("--fluxmap", required=True)
    p.add_argument("--coil-centroid", type=float, nargs=3, required=True)
    p.add_argument("--fluxmap-scale", type=float, default=100.0)
    args = p.parse_args()

    d0 = np.load(args.p0); d1 = np.load(args.p1)
    i0 = d0["importance"].astype(float); i1 = d1["importance"].astype(float)
    xc, yc, zc = _centers(d0)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    coil = np.array(args.coil_centroid)
    dist = np.sqrt((X - coil[0])**2 + (Y - coil[1])**2 + (Z - coil[2])**2)

    # --- (1) where BOTH are nonzero: the P1/P0 importance ratio ---
    both = (i0 > 0) & (i1 > 0)
    ratio = np.where(both, i1 / np.maximum(i0, 1e-300), np.nan)
    print("=== P1 vs P0 adjoint importance (real DAGMC coil) ===")
    print(f"  voxels compared: {int(both.sum())}")
    print(f"  P1/P0 ratio: median {np.nanmedian(ratio):.3f}  "
          f"mean {np.nanmean(ratio):.3f}  "
          f"5-95%% [{np.nanpercentile(ratio,5):.2f}, {np.nanpercentile(ratio,95):.2f}]")
    frac = np.nanmean(np.abs(i1[both] - i0[both]) / i0[both])
    print(f"  mean |P1-P0|/P0 = {frac:.1%}")

    # --- (2) ratio vs distance from the coil (M6: grows with depth) ---
    print("\n  P1/P0 ratio vs distance from coil (M6 predicts growth with depth):")
    dd = dist[both]; rr = ratio[both]
    edges = np.percentile(dd, [0, 20, 40, 60, 80, 100])
    for a, b in zip(edges[:-1], edges[1:]):
        m = (dd >= a) & (dd < b)
        if m.sum():
            print(f"    d in [{a/100:5.1f}, {b/100:5.1f}] m : median P1/P0 = "
                  f"{np.nanmedian(rr[m]):.3f}  (n={int(m.sum())})")

    # --- (3) contributon direction + placement target, P0 vs P1 ---
    print("\n=== contributon direction + shield-placement target ===")
    tor = np.linspace(0, 360, 48, endpoint=False)
    pol = np.linspace(0, 360, 48, endpoint=False)
    u = coil / np.linalg.norm(coil)
    for name, npz in [("P0 (scalar)", args.p0), ("P1 (angular)", args.p1)]:
        C = ap.contributon(npz, args.fluxmap, scale=args.fluxmap_scale)
        Cimp = C["importance"]
        cen = np.array([(X*Cimp).sum(), (Y*Cimp).sum(), (Z*Cimp).sum()]) / Cimp.sum()
        cos = float(cen @ coil / (np.linalg.norm(cen) * np.linalg.norm(coil)))
        res = ap.placement_priority(C, tor, pol)
        cphi, _ = ap.coil_angles(args.coil_centroid, res["R0"])
        dphi = ((res["phi0_deg"] - cphi + 180) % 360) - 180
        print(f"  {name:13s} centroid cos(coil)={cos:+.2f}  "
              f"placement phi0={res['phi0_deg']:.0f} (|dphi_coil|={abs(dphi):.0f})  "
              f"theta0={res['theta0_deg']:.0f}")

    print("\nInterpretation: if P1/P0 grows with depth and the placement target shifts, "
          "angular treatment materially changes coil attribution (M6 regime confirmed on "
          "the real coil). If ratios ~1 and the target is unchanged, scalar suffices here.")


if __name__ == "__main__":
    main()
