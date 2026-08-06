#!/usr/bin/env python
"""Test the SIGNED trade worth against the Tier-2 finite differences.

Tier-2 showed the scalar worth W = phi*psi_dagger (shield band only) is ANTI-correlated with the
true finite-difference worth W_FD (Spearman -1.0 on 4 corners): it ranks highest the patch (2,3)
where adding shield actually RAISED coil-20 dose, because the fixed-envelope trade also REMOVES
breeder and the scalar worth ignores that. The signed worth debits the breeder band:

    W_signed = sigma_sh * sum_[shield band] C  -  sigma_br * sum_[breeder band] C ,   C = phi*psi_dagger

This script computes A, W (scalar), and W_signed on the (n_tor x n_pol) grid, reads W_FD for the
built+dosed patches, and reports Spearman(W_FD, .) for each -- the test of whether the signed
worth is the quantity that actually predicts actionability.
"""
import argparse
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))     # shield_opt
sys.path.insert(0, HERE)                       # patch_worth
import patch_worth as pw            # noqa: E402
import patch_worth_analyze as pwa   # noqa: E402
import plasma_geometry as pg        # noqa: E402
from scipy.stats import spearmanr   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adjoint", required=True)
    ap.add_argument("--fluxmap", required=True)
    ap.add_argument("--fwd", required=True)
    ap.add_argument("--baseline", required=True, help="uniform coil dose npz (delta=0)")
    ap.add_argument("--patch-glob", required=True, help="coil_patch_t*.npz")
    ap.add_argument("--delta-dir", default=None, help="dir with delta_qh_patch_*.npz (default: patch dir)")
    ap.add_argument("--n-tor", type=int, default=4)
    ap.add_argument("--n-pol", type=int, default=4)
    ap.add_argument("--nfp", type=int, default=4)
    ap.add_argument("--cell", type=int, default=20)
    ap.add_argument("--sigma-shield", type=float, default=1.0 / 8.0)
    ap.add_argument("--sigma-breeder", type=float, default=1.0 / 17.0)
    args = ap.parse_args()

    R0, a_minor, _ = pg.axis_minor_from_fluxmap(args.fluxmap)
    print(f"geometry: R0={R0:.1f} cm  a_minor={a_minor:.1f} cm")
    te, pe = pw.define_patches(args.n_tor, args.n_pol, nfp=args.nfp)

    A, _ = pw.attribution_patches(args.adjoint, args.fluxmap, te, pe, R0=R0, nfp=args.nfp)
    W, _ = pw.worth_patches(args.adjoint, args.fwd, R0, te, pe, a_minor=a_minor, nfp=args.nfp)
    Wsig, meta = pw.signed_worth_patches(args.adjoint, args.fwd, R0, te, pe, a_minor=a_minor,
                                         sigma_shield=args.sigma_shield,
                                         sigma_breeder=args.sigma_breeder, nfp=args.nfp)
    print(f"bands: shield {pg.shield_band_offsets()} ({meta['n_sh']} vox)  "
          f"breeder {pg.breeder_band_offsets()} ({meta['n_br']} vox)")

    # finite-difference worth per built+dosed patch
    R0dose = pwa._coil_R(args.baseline, args.cell)
    ddir = args.delta_dir or os.path.dirname(os.path.abspath(args.patch_glob))
    rows = []
    for p in sorted(glob.glob(args.patch_glob)):
        tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
        i = int(tag.split("t")[1].split("p")[0]); j = int(tag.split("p")[1])
        Rp = pwa._coil_R(p, args.cell)
        dd = np.load(os.path.join(ddir, f"delta_qh_patch_{tag}.npz"))
        added = float(dd["added_cells"]) * float(dd["delta_cm"])
        W_FD = -(Rp - R0dose) / added
        rows.append(dict(tag=tag, i=i, j=j, W_FD=W_FD,
                         A=float(A[i, j]), W=float(W[i, j]), Wsig=float(Wsig[i, j])))

    print(f"\n{'patch':7} {'W_FD':>11} {'A':>11} {'W(scalar)':>11} {'W_signed':>11}")
    for r in rows:
        print(f"{r['tag']:7} {r['W_FD']:+.3e} {r['A']:11.3e} {r['W']:11.3e} {r['Wsig']:+.3e}")

    wfd = np.array([r["W_FD"] for r in rows])
    def rho(name, key):
        v = np.array([r[key] for r in rows])
        r_, p_ = spearmanr(wfd, v)
        return f"  Spearman(W_FD, {name:9}) = {r_:+.3f}  (n={len(rows)})"
    print("\n=== does it predict the finite-difference worth? ===")
    print(rho("A", "A")); print(rho("W scalar", "W")); print(rho("W_signed", "Wsig"))
    print("\nExpected: scalar W anti-correlated (~-1); W_signed should recover a POSITIVE "
          "correlation if debiting the breeder band is the missing physics.")


if __name__ == "__main__":
    main()
