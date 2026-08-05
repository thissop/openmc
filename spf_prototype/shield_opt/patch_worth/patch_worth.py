#!/usr/bin/env python
"""Tier 1 of the patch-worth experiment: attribution A(theta,phi) vs scalar shielding
worth W(theta,phi) on ~16 (theta,phi) shield patches, QH coil-20.

CLEAN IDEA (see README): the contributon density C(x) = phi(x) * psi_dagger(x) is
defined everywhere.
  A(patch) = contributon in the PLASMA  (S * psi_dagger, S = plasma source) projected
             to the shield surface = "where responsible neutrons are BORN" = our current
             `placed` attribution logic.
  W(patch) = contributon in the SHIELD  (phi_fwd * psi_dagger summed over the shield
             radial band per angular column) = "where blocking a neutron HELPS" =
             leading-order  -dR/dm  (times the added-material removal XS, which cancels
             in a ranking).
Both are (theta,phi) fields; they differ because transport carries the contributon from
plasma to shield along streaming paths. We rank-compare them. Tier 2 (build_patch.py +
slurm) validates W against finite-difference truth.

CAVEAT: phi*psi_dagger is the SCALAR leading-order worth; exact anisotropic worth needs
angular moments (Jack). This W is a heuristic we validate in Tier 2 -- not claimed exact.

Inputs (all things we have except the forward mesh flux, which slurm step 1 produces):
  --adjoint   adjoint_importance_flat_P0.npz   (key 'importance' = psi_dagger on mesh)
  --fluxmap   qh_freestream_fluxmap.npz        (plasma source S, uniform emissivity)
  --fwd       qh_fwd_meshflux.npz              (forward flux phi on the SAME mesh; opt.)
Outputs: patch_worth_tier1.csv  (+ prints Spearman rho and Chatterjee xi of A vs W).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # shield_opt/
from adjoint_placement import _voxel_centers, contributon, importance_angles  # noqa: E402
import plasma_geometry as pg  # noqa: E402


def _fold_tor(phi, nfp):
    """Fold toroidal angle into one field period [0, 2pi/nfp). The shield design
    delta(theta,phi) is stellarator-symmetric (per-period, repeated nfp times), so
    worth/attribution must be measured in the period frame to match the Tier-2 build.
    The coil-k contributon peaks in one period (psi_dagger is coil-localized), so
    folding concentrates that signal into the period-frame patch -- it does not mix
    coils."""
    return np.mod(phi, 2 * np.pi / nfp)


def define_patches(n_tor=4, n_pol=4, nfp=4):
    """~16 (theta,phi) patches: toroidal phi over ONE field period [0, 2pi/nfp) and
    poloidal theta in [-pi,pi). Returns (tor_edges, pol_edges) in radians."""
    return (np.linspace(0.0, 2 * np.pi / nfp, n_tor + 1),
            np.linspace(-np.pi, np.pi, n_pol + 1))


def _bin2d(phi, theta, w, tor_edges, pol_edges):
    """Sum weights w into the (phi,theta) patch grid. Returns (n_tor, n_pol)."""
    H, _, _ = np.histogram2d(phi, theta, bins=[tor_edges, pol_edges], weights=w)
    return H


def attribution_patches(adjoint, fluxmap, tor_edges, pol_edges, R0=None,
                        emissivity="uniform", nfp=4):
    """A(patch) = plasma contributon C = S*psi_dagger, binned by voxel (phi,theta),
    toroidal folded into one field period."""
    field = contributon(adjoint, fluxmap, emissivity=emissivity)
    phi, theta, w, R0 = importance_angles(field, R0=R0)  # nonzero-contributon voxels
    return _bin2d(_fold_tor(phi, nfp), theta, w, tor_edges, pol_edges), R0


def worth_patches(adjoint, fwd, R0, tor_edges, pol_edges,
                  a_minor=None, band_cm=None, nfp=4):
    """W(patch) = sum over SHIELD-band voxels of phi_fwd*psi_dagger, binned by (phi,theta).

    Shield band = minor-radius-like distance rho = hypot(R-R0, Z) in
    [a_minor + band_cm[0], a_minor + band_cm[1]]  (breeder_out..shield_out beyond LCFS).
    band_cm defaults to the ACTUAL radial-build offsets (plasma_geometry.shield_band_offsets:
    FW+mult+breeder+back_wall -> +shield), so the band tracks the real build. a_minor should be
    the geometry-derived LCFS minor radius (plasma_geometry.axis_minor_from_fluxmap); if None it
    falls back to a crude 5th-pct-rho voxel heuristic (kept only so the function runs standalone).
    Removal XS is a constant prefactor -> omitted (cancels in the A-vs-W ranking)."""
    if band_cm is None:
        band_cm = pg.shield_band_offsets()
    da = np.load(adjoint)
    psi = da["importance"].astype(float)
    ll, ur, dim = da["lower_left"], da["upper_right"], da["dimension"]
    df = np.load(fwd)
    phi_f = df["flux"].astype(float).reshape(psi.shape)

    xc, yc, zc = _voxel_centers(ll, ur, dim)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    R = np.hypot(X, Y)
    rho = np.hypot(R - R0, Z)                      # minor-radius-like distance from axis
    if a_minor is None:
        # crude LCFS estimate: where the forward flux collapses is inside the plasma;
        # use the radius enclosing the bulk of psi*phi below the band as a fallback.
        a_minor = float(np.percentile(rho[psi > 0], 5.0))  # replace w/ true a if known
    r_in, r_out = a_minor + band_cm[0], a_minor + band_cm[1]
    shield = (rho >= r_in) & (rho <= r_out)

    C = phi_f * psi                                # contributon density everywhere
    phi_ang = _fold_tor(np.arctan2(Y, X)[shield], nfp)
    theta_ang = np.arctan2(Z, R - R0)[shield]
    w = C[shield]
    W = _bin2d(phi_ang, theta_ang, w, tor_edges, pol_edges)
    return W, dict(a_minor=a_minor, r_in=r_in, r_out=r_out, n_shield_vox=int(shield.sum()))


def spearman(x, y):
    from scipy.stats import spearmanr
    r, p = spearmanr(x, y)
    return float(r), float(p)


def chatterjee_xi(x, y):
    """Chatterjee's xi(x->y) rank correlation (asymmetric, 0=indep, 1=y is f(x)).
    No-ties form: order y by ascending x, r_i = rank of y_(i); xi = 1 - 3*sum|dr|/(n^2-1)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    if n < 2:
        return float("nan")
    y_by_x = y[np.argsort(x, kind="mergesort")]           # y ordered by ascending x
    ry = np.argsort(np.argsort(y_by_x, kind="mergesort")) + 1  # rank of each y_(i), 1..n
    return float(1.0 - 3.0 * np.sum(np.abs(np.diff(ry))) / (n * n - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adjoint", required=True, help="adjoint_importance_flat_P0.npz")
    ap.add_argument("--fluxmap", required=True, help="plasma-source fluxmap npz")
    ap.add_argument("--fwd", default=None,
                    help="forward mesh-flux npz (same mesh). If absent, only A is computed.")
    ap.add_argument("--n-tor", type=int, default=4)
    ap.add_argument("--n-pol", type=int, default=4)
    ap.add_argument("--nfp", type=int, default=4, help="field periods (QH=4)")
    ap.add_argument("--R0", type=float, default=None,
                    help="magnetic-axis major radius (cm); default = geometry-derived from fluxmap")
    ap.add_argument("--a-minor", type=float, default=None,
                    help="LCFS minor radius (cm); default = geometry-derived from fluxmap")
    ap.add_argument("--band", type=float, nargs=2, default=None,
                    help="shield band beyond LCFS (cm); default = actual radial-build offsets")
    ap.add_argument("--emissivity", default="uniform")
    ap.add_argument("--out", default=os.path.join(HERE, "patch_worth_tier1.csv"))
    args = ap.parse_args()

    # Physically-grounded R0 / a_minor from the plasma boundary (replaces the voxel heuristics).
    # Falls back gracefully if the fluxmap lacks R/Z.
    R0_geo = a_geo = None
    try:
        R0_geo, a_geo, _amax = pg.axis_minor_from_fluxmap(args.fluxmap)
        print(f"plasma geometry (fluxmap): R0={R0_geo:.1f} cm  a_minor={a_geo:.1f} cm")
    except Exception as e:  # noqa: BLE001
        print(f"NOTE: could not derive R0/a_minor from fluxmap ({e}); using voxel heuristics.")
    R0_used = args.R0 if args.R0 is not None else R0_geo
    a_used = args.a_minor if args.a_minor is not None else a_geo

    tor_e, pol_e = define_patches(args.n_tor, args.n_pol, nfp=args.nfp)
    A, R0 = attribution_patches(args.adjoint, args.fluxmap, tor_e, pol_e,
                                R0=R0_used, emissivity=args.emissivity, nfp=args.nfp)
    print(f"R0 = {R0:.1f} cm   patches = {args.n_tor}x{args.n_pol} = {A.size}  "
          f"(toroidal folded to 1/{args.nfp} period)")

    if args.fwd and os.path.exists(args.fwd):
        band = tuple(args.band) if args.band is not None else None
        W, meta = worth_patches(args.adjoint, args.fwd, R0, tor_e, pol_e,
                                a_minor=a_used, band_cm=band, nfp=args.nfp)
        print(f"shield band rho in [{meta['r_in']:.0f}, {meta['r_out']:.0f}] cm "
              f"(a={meta['a_minor']:.0f}), {meta['n_shield_vox']} voxels")
        a, w = A.ravel(), W.ravel()
        keep = (a > 0) | (w > 0)
        a, w = a[keep], w[keep]
        rho_s, p_s = spearman(a, w)
        xi = chatterjee_xi(a, w)
        # top-k overlap (do attribution and worth agree on the hottest patches?)
        k = max(1, len(a) // 4)
        topA = set(np.argsort(a)[-k:]); topW = set(np.argsort(w)[-k:])
        overlap = len(topA & topW) / k
        print("\n=== Attribution A  vs  Shielding-worth W (patch ranking) ===")
        print(f"  Spearman rho   = {rho_s:+.3f}  (p={p_s:.2e})")
        print(f"  Chatterjee xi  = {xi:+.3f}")
        print(f"  top-{k} overlap = {overlap:.2f}")
        print("  interpretation: high rho/overlap -> attribution IS a worth proxy;")
        print("                  low  rho/overlap -> attribution != actionability (novel).")
    else:
        W = np.full_like(A, np.nan)
        print("NOTE: no --fwd forward mesh flux yet -> only attribution A computed. "
              "Run slurm_patch_worth.sh step 1 to produce qh_fwd_meshflux.npz, then rerun.")

    with open(args.out, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["i_tor", "j_pol", "phi_deg", "theta_deg", "A_attribution", "W_worth"])
        tc = np.degrees(0.5 * (tor_e[:-1] + tor_e[1:]))
        pc = np.degrees(0.5 * (pol_e[:-1] + pol_e[1:]))
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                wr.writerow([i, j, f"{tc[i]:.1f}", f"{pc[j]:.1f}",
                             f"{A[i, j]:.6e}", f"{W[i, j]:.6e}"])
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
