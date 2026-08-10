#!/usr/bin/env python
"""Step (2): does a FAST-group signed worth predict the n=12 finite-difference worth better than the
single-group (total-flux) signed worth?

Motivation (SIGNED_WORTH_ANALYSIS.md + the n=12 verdict in RESULTS_tier2.md): the coil response is
FAST flux, and the breeder's role in the fixed-envelope trade is to MODERATE fast neutrons -- so the
contributon C = phi*psi_dagger should use the FAST forward flux, not the total flux that over-weights
the thermalized breeder region. This driver recomputes the signed worth with the fast forward flux
(qh_fwd_meshflux_fast.npz, job 9354536) and compares its rank correlation with the finite-difference
worth against the total-flux baseline, across a sigma_shield/sigma_breeder sweep.

It is honest by construction: it reports BOTH fluxes' correlations and the sigma-sweep, so a null
result ("fast flux does not help") is as visible as a positive one. No tuning to a target.

Usage:
  python multigroup_signed_worth.py --adjoint A.npz --fwd-total T.npz --fwd-fast F.npz \
      --patch-glob 'cluster_pull/coil_patch_*.npz' --delta-dir cluster_pull \
      --baseline cluster_pull/coil_step1b_uniform.npz [--R0 1367 --a-minor 217]
"""
import argparse
import glob
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))     # shield_opt (plasma_geometry)
import patch_worth as pw            # noqa: E402
import patch_worth_analyze as pwa   # noqa: E402


def fd_worths(patch_glob, delta_dir, baseline, cell=20):
    """Finite-difference worth per dosed patch: W_FD = -(R_patch - R0) / added_material."""
    R0 = pwa._coil_R(baseline, cell)
    rows = []
    for p in sorted(glob.glob(patch_glob)):
        tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
        i = int(tag.split("t")[1].split("p")[0]); j = int(tag.split("p")[1])
        Rp = pwa._coil_R(p, cell)
        dd = np.load(os.path.join(delta_dir, f"delta_qh_patch_{tag}.npz"))
        added = float(dd["added_cells"]) * float(dd["delta_cm"])
        rows.append(dict(tag=tag, i=i, j=j, W_FD=-(Rp - R0) / added))
    return rows


def signed_grid(adjoint, fwd, R0, a_minor, te, pe, sig_sh, sig_br, nfp=4):
    Wsig, _ = pw.signed_worth_patches(adjoint, fwd, R0, te, pe, a_minor=a_minor,
                                      sigma_shield=sig_sh, sigma_breeder=sig_br, nfp=nfp)
    return Wsig


def correlate(rows, Wgrid):
    wfd = np.array([r["W_FD"] for r in rows])
    v = np.array([Wgrid[r["i"], r["j"]] for r in rows])
    rho, p = spearmanr(wfd, v)
    sign = float(np.mean(np.sign(wfd) == np.sign(v)))
    return rho, p, sign


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adjoint", required=True)
    ap.add_argument("--fwd-total", required=True)
    ap.add_argument("--fwd-fast", default=None, help="if omitted, only the total-flux baseline runs")
    ap.add_argument("--patch-glob", required=True)
    ap.add_argument("--delta-dir", required=True)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--R0", type=float, default=1367.0)
    ap.add_argument("--a-minor", type=float, default=217.0)
    ap.add_argument("--nfp", type=int, default=4)
    ap.add_argument("--cell", type=int, default=20)
    ap.add_argument("--sigma-shield", type=float, default=1.0 / 8.0)
    ap.add_argument("--sigma-breeder", type=float, default=1.0 / 17.0)
    args = ap.parse_args()

    te, pe = pw.define_patches(4, 4, nfp=args.nfp)
    rows = fd_worths(args.patch_glob, args.delta_dir, args.baseline, cell=args.cell)
    print(f"n={len(rows)} dosed patches\n")

    fluxes = [("total", args.fwd_total)]
    if args.fwd_fast and os.path.exists(args.fwd_fast):
        fluxes.append(("fast", args.fwd_fast))
    else:
        print("(fast forward flux not present yet -- only the total-flux baseline is shown)\n")

    # default-sigma comparison
    print(f"=== signed worth vs W_FD  (sigma_sh={args.sigma_shield:.3g}, sigma_br={args.sigma_breeder:.3g}) ===")
    print(f"{'flux':6}{'Spearman':>10}{'p':>9}{'sign-match':>12}")
    grids = {}
    for name, f in fluxes:
        W = signed_grid(args.adjoint, f, args.R0, args.a_minor, te, pe,
                        args.sigma_shield, args.sigma_breeder, nfp=args.nfp)
        grids[name] = W
        rho, p, sign = correlate(rows, W)
        print(f"{name:6}{rho:+10.3f}{p:9.2f}{sign:12.2f}")

    # sigma-ratio robustness sweep (does the conclusion depend on the assumed XS ratio?)
    print("\n=== sigma_sh/sigma_br robustness sweep (Spearman vs W_FD) ===")
    ratios = [0.5, 1.0, 2.13, 3.0, 5.0]     # 2.13 = WC/FLiBe lambda ratio
    hdr = "  ratio " + "".join(f"{n:>10}" for n, _ in fluxes)
    print(hdr)
    for r in ratios:
        sig_br = 1.0 / 17.0; sig_sh = r * sig_br
        line = f"  {r:5.2f} "
        for name, f in fluxes:
            W = signed_grid(args.adjoint, f, args.R0, args.a_minor, te, pe, sig_sh, sig_br, nfp=args.nfp)
            rho, _, _ = correlate(rows, W)
            line += f"{rho:+10.3f}"
        print(line)

    # per-patch table (total vs fast) for the eyeball check on the backfire patch
    if "fast" in grids:
        print("\n=== per-patch (sorted by W_FD) ===")
        print(f"{'patch':7}{'W_FD':>12}{'Wsig_total':>13}{'Wsig_fast':>12}")
        for r in sorted(rows, key=lambda r: -r["W_FD"]):
            wt = grids["total"][r["i"], r["j"]]; wf = grids["fast"][r["i"], r["j"]]
            print(f"{r['tag']:7}{r['W_FD']:+.3e}{wt:+13.3e}{wf:+12.3e}")


if __name__ == "__main__":
    main()
