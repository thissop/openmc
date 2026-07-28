#!/usr/bin/env python
"""Tier-2 analysis: finite-difference worth W_FD vs adjoint worth W vs attribution A.

Two modes:
  --pick-corners   read patch_worth_tier1.csv, choose ~6-8 patches spanning the
                   {high/low A} x {high/low W} corners, write patches.txt (the SLURM
                   array list). Deliberately includes disagreement corners (high-A/low-W
                   and low-A/high-W) so the experiment can DISTINGUISH attribution from worth.
  --finite-diff    read baseline coil flux + coil_patch_*.npz, compute
                   W_FD = -(R_patch - R_baseline)/added_cells, and correlate:
                     W_FD vs W  -> validates the adjoint scalar worth
                     W_FD vs A  -> tests attribution != actionability
                   reports Spearman + sign accuracy + relative-magnitude error, and the figure.

The decisive read: if W_FD tracks W but NOT A, "attribution is not a placement proxy;
worth is the lever" -- the novel result. If W_FD tracks both, attribution is a cheap proxy.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load_tier1(path):
    rows = list(csv.DictReader(open(path)))
    key = {}
    for r in rows:
        key[(int(r["i_tor"]), int(r["j_pol"]))] = (
            float(r["A_attribution"]),
            float(r["W_worth"]) if r["W_worth"] not in ("", "nan") else np.nan,
        )
    return key


def pick_corners(tier1_csv, out, n_per_corner=2):
    key = load_tier1(tier1_csv)
    ij = list(key)
    A = np.array([key[k][0] for k in ij])
    W = np.array([key[k][1] for k in ij])
    if np.all(np.isnan(W)):
        raise SystemExit("W is all-NaN: run fwd_meshflux + patch_worth.py --fwd first.")
    ra = A.argsort().argsort() / (len(A) - 1)     # normalized ranks in [0,1]
    rw = W.argsort().argsort() / (len(W) - 1)
    corners = {
        "hiA_hiW": (ra + rw),
        "hiA_loW": (ra + (1 - rw)),               # attribution says yes, worth says no
        "loA_hiW": ((1 - ra) + rw),               # worth says yes, attribution says no
        "loA_loW": ((1 - ra) + (1 - rw)),
    }
    chosen = []
    for name, score in corners.items():
        for idx in np.argsort(score)[-n_per_corner:]:
            if ij[idx] not in [c[0] for c in chosen]:
                chosen.append((ij[idx], name))
    with open(out, "w") as f:
        for (i, j), name in chosen:
            f.write(f"{i} {j}\n")
    print(f"wrote {out} with {len(chosen)} patches:")
    for (i, j), name in chosen:
        print(f"  ({i},{j})  {name}   A={key[(i,j)][0]:.3e}  W={key[(i,j)][1]:.3e}")


def finite_diff(tier1_csv, baseline_npz, patch_glob, cell=20, out_fig=None):
    key = load_tier1(tier1_csv)
    base = np.load(baseline_npz)
    cells = list(int(c) for c in base["cells"]); R0 = float(base["flux"][cells.index(cell)])
    rows = []
    for p in sorted(glob.glob(patch_glob)):
        tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]  # tIpJ
        i = int(tag.split("t")[1].split("p")[0]); j = int(tag.split("p")[1])
        d = np.load(p)
        cs = list(int(c) for c in d["cells"]); Rp = float(d["flux"][cs.index(cell)])
        dd = np.load(p.replace("coil_patch_", "delta_qh_patch_"))
        added = float(dd["added_cells"]) * float(dd["delta_cm"])
        W_FD = -(Rp - R0) / added                       # worth = response drop per unit added shield
        A, W = key[(i, j)]
        rows.append(dict(i=i, j=j, R=Rp, W_FD=W_FD, A=A, W=W))
        print(f"  patch({i},{j}): R={Rp:.4e}  dR={Rp-R0:+.3e}  W_FD={W_FD:+.3e}  "
              f"A={A:.3e}  W={W:.3e}")

    from scipy.stats import spearmanr
    wfd = np.array([r["W_FD"] for r in rows])
    A = np.array([r["A"] for r in rows]); W = np.array([r["W"] for r in rows])
    rWFD_W, pW = spearmanr(wfd, W)
    rWFD_A, pA = spearmanr(wfd, A)
    sign_W = float(np.mean(np.sign(wfd) == np.sign(W)))
    print("\n=== Tier-2 finite-difference validation ===")
    print(f"  baseline R(cell {cell}) = {R0:.4e}")
    print(f"  Spearman W_FD vs W (adjoint worth)   = {rWFD_W:+.3f}  (p={pW:.2e})  "
          f"sign-match={sign_W:.2f}")
    print(f"  Spearman W_FD vs A (attribution)     = {rWFD_A:+.3f}  (p={pA:.2e})")
    print("  VERDICT: W_FD~W & W_FD!~A -> attribution != actionability (novel);")
    print("           W_FD~W & W_FD~A  -> attribution is a validated cheap proxy;")
    print("           W_FD!~W          -> scalar worth inadequate (need angular / FD ranking).")

    if out_fig:
        import smplotlib  # noqa: F401
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        ax[0].scatter(W, wfd, s=40, c="#2a6f97", edgecolor="k", linewidth=0.3)
        ax[0].set_xlabel(r"Adjoint Worth $W=\phi\,\psi^\dagger$ (Shield)")
        ax[0].set_ylabel(r"Finite-Difference Worth $-\Delta R/\Delta M$")
        ax[0].text(0.05, 0.9, fr"$\rho={rWFD_W:+.2f}$", transform=ax[0].transAxes)
        ax[1].scatter(A, wfd, s=40, c="#e07a2f", edgecolor="k", linewidth=0.3)
        ax[1].set_xlabel(r"Attribution $A=S\,\psi^\dagger$ (Plasma)")
        ax[1].set_ylabel(r"Finite-Difference Worth $-\Delta R/\Delta M$")
        ax[1].text(0.05, 0.9, fr"$\rho={rWFD_A:+.2f}$", transform=ax[1].transAxes)
        fig.tight_layout(); fig.savefig(out_fig, dpi=150); plt.close(fig)
        print(f"wrote {out_fig}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier1", default=os.path.join(HERE, "patch_worth_tier1.csv"))
    ap.add_argument("--pick-corners", action="store_true")
    ap.add_argument("--patches-out", default=os.path.join(HERE, "patches.txt"))
    ap.add_argument("--finite-diff", action="store_true")
    ap.add_argument("--baseline", help="baseline coil npz (cells, flux)")
    ap.add_argument("--patch-glob", default=None, help="coil_patch_*.npz glob")
    ap.add_argument("--cell", type=int, default=20)
    ap.add_argument("--fig", default=os.path.join(HERE, "patch_worth_tier2.png"))
    args = ap.parse_args()

    if args.pick_corners:
        pick_corners(args.tier1, args.patches_out)
    if args.finite_diff:
        finite_diff(args.tier1, args.baseline, args.patch_glob, cell=args.cell,
                    out_fig=args.fig)


if __name__ == "__main__":
    main()
