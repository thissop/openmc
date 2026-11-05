"""Prototype: spectral-truth-driven ADAPTIVE radial grid for the fusion source.

The representation study showed the gridded source's error is dominated by the
near-axis Jacobian sqrt(g), which barely resolves away under uniform refinement.
A uniform radial (s) grid spends nodes where sqrt(g) is smooth (large s) and
starves the near-axis region where it is not. This builds a NON-uniform radial
grid by classic error equidistribution against the spectral truth, and measures
whether it reaches a target sqrt(g)/B-hat accuracy at fewer nodes (less memory)
than a uniform grid -- i.e. an adaptive accuracy/memory frontier.

Method (equidistribution for piecewise-linear interpolation): the leading
interpolation error on a cell scales with |f''| * h^2, so the error is
equidistributed when node spacing h ~ |f''|^(-1/2), i.e. when the nodes equalise
the integral of the monitor M(s) = |f''(s)|^(1/2). We monitor f = sqrt(g)
averaged over (theta, zeta) (the birth-density weight, and the worst-behaved
field), floor M to keep some coverage everywhere, integrate, and invert.

This is a first, honest version: radial adaptivity only (theta, zeta stay
uniform), because the study localised the dominant error in s near the axis. A
full 3-D monitor is the next step.
"""
from __future__ import annotations
import numpy as np

from fluxrep import SpectralRep, GridRep
import metrics


def _sqrtg_profile(spectral: SpectralRep, ns=600, nang=512, nfp=1, seed=0):
    """Mean_|theta,zeta| |sqrt(g)|(s) on a fine s-grid, for the error monitor."""
    rng = np.random.default_rng(seed)
    th = 2 * np.pi * rng.random(nang)
    ze = (2 * np.pi / nfp) * rng.random(nang)
    s_fine = np.linspace(1e-4, 1.0, ns)
    prof = np.empty(ns)
    for i, s in enumerate(s_fine):
        _, _, g = spectral.eval(np.full(nang, s), th, ze, want_sqrtg=True)
        prof[i] = np.nanmean(g)
    return s_fine, prof


def adaptive_s_nodes(spectral: SpectralRep, nrho, nfp=1, floor_frac=0.15, ns=600):
    """Return nrho error-equidistributed s-nodes in [~0, 1] (endpoints included).

    Monitor M(s) = sqrt(|d2(sqrtg)/ds2|), floored at floor_frac of its max so the
    outer region keeps coverage. Nodes equalise the cumulative integral of M.
    """
    s_fine, f = _sqrtg_profile(spectral, ns=ns, nfp=nfp)
    d2 = np.gradient(np.gradient(f, s_fine), s_fine)
    M = np.sqrt(np.abs(d2))
    M = np.maximum(M, floor_frac * M.max())        # floor so nowhere is starved
    C = np.concatenate([[0.0], np.cumsum(0.5 * (M[1:] + M[:-1]) * np.diff(s_fine))])
    C /= C[-1]
    targets = np.linspace(0.0, 1.0, nrho)
    s_nodes = np.interp(targets, C, s_fine)
    s_nodes[0] = max(s_fine[0], 1e-4)
    s_nodes[-1] = 1.0
    return s_nodes


def _angular_monitor(spectral: SpectralRep, axis, nfp, nfine=720, nsamp=256, seed=1):
    """Periodic error monitor M(angle) on [0, L): mean over (s, other-angle) of the
    interpolation-relevant curvature |f''| of the geometry (R, Z) and sqrt(g).

    axis='theta' -> L=2*pi, vary theta; axis='zeta' -> L=2*pi/nfp, vary zeta.
    """
    L = 2 * np.pi if axis == "theta" else 2 * np.pi / nfp
    ang = np.linspace(0.0, L, nfine, endpoint=False)
    dA = L / nfine
    rng = np.random.default_rng(seed)
    s = rng.uniform(1e-3, 1.0, nsamp)
    other_L = 2 * np.pi / nfp if axis == "theta" else 2 * np.pi
    other = rng.uniform(0.0, other_L, nsamp)

    # Sample each field on the fine angular grid at the fixed (s, other-angle) set.
    # Include the field direction B-hat, whose angular variation (from the field-
    # line pitch) is what wrecks a coarse angular grid -- geometry curvature alone
    # is the wrong monitor for the polarized source.
    has_B = spectral.has_field
    Rr = np.empty((nfine, nsamp)); Zz = np.empty((nfine, nsamp)); Gg = np.empty((nfine, nsamp))
    Bx = np.empty((nfine, nsamp)); By = np.empty((nfine, nsamp)); Bz = np.empty((nfine, nsamp))
    for k, a in enumerate(ang):
        if axis == "theta":
            th, ze = np.full(nsamp, a), other
        else:
            th, ze = other, np.full(nsamp, a)
        Rr[k], Zz[k], Gg[k] = spectral.eval(s, th, ze, want_sqrtg=True)
        if has_B:
            b = spectral.eval_bhat(s, th, ze)
            Bx[k], By[k], Bz[k] = b[:, 0], b[:, 1], b[:, 2]

    def curv(arr):                       # periodic |second difference|, mean over the sample set
        d2 = (np.roll(arr, -1, 0) - 2 * arr + np.roll(arr, 1, 0)) / dA**2
        return np.nanmean(np.abs(d2), axis=1)

    def norm(a):                         # scale each field to unit mean so terms are comparable
        return a / (np.nanmean(a) + 1e-30)

    terms = norm(curv(Rr)) + norm(curv(Zz)) + norm(curv(Gg))
    if has_B:                            # B-hat direction curvature, the field-aware term
        terms = terms + norm(curv(Bx) + curv(By) + curv(Bz))
    # equidistribute sqrt(sum of normalized curvatures) -> node spacing ~ |f''|^(-1/2)
    M = np.sqrt(terms)
    return ang, M, L


def adaptive_angular_nodes(spectral: SpectralRep, n, axis, nfp, floor_frac=0.2):
    """n periodic error-equidistributed nodes on [0, L), endpoint-free."""
    ang, M, L = _angular_monitor(spectral, axis, nfp)
    M = np.maximum(M, floor_frac * M.max())
    # cumulative integral on the periodic grid (endpoint-free), then place nodes
    dA = L / len(ang)
    C = np.concatenate([[0.0], np.cumsum(M) * dA])   # C at ang[0..n], C[-1]=integral over full period
    ang_ext = np.concatenate([ang, [L]])
    targets = np.linspace(0.0, C[-1], n, endpoint=False)
    nodes = np.interp(targets, C, ang_ext)
    return np.mod(nodes, L)


def compare(spectral: SpectralRep, nfp, nrho_list=(6, 8, 10, 12, 16, 24),
            ntheta=96, nzeta=288, nquery=20000, full_3d=False):
    """Uniform vs adaptive grid: sqrt(g) RMS and B-hat p99 vs nrho (memory).

    full_3d=True also builds a grid adaptive in all three axes (radial + theta +
    zeta), at the SAME node counts, to isolate the angular contribution.
    """
    th_ad = adaptive_angular_nodes(spectral, ntheta, "theta", nfp) if full_3d else None
    ze_ad = adaptive_angular_nodes(spectral, nzeta, "zeta", nfp) if full_3d else None
    rows = []
    for nrho in nrho_list:
        uni = GridRep(spectral, nrho=nrho, ntheta=ntheta, nzeta=nzeta, nfp=nfp)
        s_ad = adaptive_s_nodes(spectral, nrho, nfp=nfp)
        ada = GridRep(spectral, nrho=nrho, ntheta=ntheta, nzeta=nzeta, nfp=nfp,
                      s_nodes=s_ad)
        au = metrics.accuracy(spectral, uni, n=nquery, nfp=nfp)
        aa = metrics.accuracy(spectral, ada, n=nquery, nfp=nfp)
        row = dict(
            nrho=nrho, MiB=float(uni.nbytes / 2**20),
            uni_sqrtg_rms=au["sqrtg"]["rms"],
            ada_sqrtg_rms=aa["sqrtg"]["rms"],
            uni_bhat_p99=au.get("bhat_degrees", {}).get("source_weighted_p99"),
            ada_bhat_p99=aa.get("bhat_degrees", {}).get("source_weighted_p99"),
        )
        if full_3d:
            a3 = GridRep(spectral, nrho=nrho, ntheta=ntheta, nzeta=nzeta, nfp=nfp,
                         s_nodes=s_ad, theta_nodes=th_ad, zeta_nodes=ze_ad)
            a3a = metrics.accuracy(spectral, a3, n=nquery, nfp=nfp)
            row["ada3d_sqrtg_rms"] = a3a["sqrtg"]["rms"]
            row["ada3d_bhat_p99"] = a3a.get("bhat_degrees", {}).get("source_weighted_p99")
        rows.append(row)
    return rows


def frontier_3d(spectral, nfp, configs, nquery=20000):
    """uniform vs radial-adaptive vs full-3D-adaptive at each (nrho,ntheta,nzeta)."""
    th_cache = {}; ze_cache = {}
    rows = []
    for (nr, nt, nz) in configs:
        uni = GridRep(spectral, nrho=nr, ntheta=nt, nzeta=nz, nfp=nfp)
        s_ad = adaptive_s_nodes(spectral, nr, nfp=nfp)
        rad = GridRep(spectral, nrho=nr, ntheta=nt, nzeta=nz, nfp=nfp, s_nodes=s_ad)
        th_ad = th_cache.setdefault(nt, adaptive_angular_nodes(spectral, nt, "theta", nfp))
        ze_ad = ze_cache.setdefault(nz, adaptive_angular_nodes(spectral, nz, "zeta", nfp))
        a3 = GridRep(spectral, nrho=nr, ntheta=nt, nzeta=nz, nfp=nfp,
                     s_nodes=s_ad, theta_nodes=th_ad, zeta_nodes=ze_ad)
        au = metrics.accuracy(spectral, uni, n=nquery, nfp=nfp)
        ar = metrics.accuracy(spectral, rad, n=nquery, nfp=nfp)
        a3a = metrics.accuracy(spectral, a3, n=nquery, nfp=nfp)
        bkey = ("bhat_degrees", "source_weighted_p99")
        rows.append(dict(
            config=[nr, nt, nz], MiB=float(uni.nbytes / 2**20),
            uni_sqrtg=au["sqrtg"]["rms"], rad_sqrtg=ar["sqrtg"]["rms"], a3_sqrtg=a3a["sqrtg"]["rms"],
            uni_bhat=au[bkey[0]][bkey[1]], rad_bhat=ar[bkey[0]][bkey[1]], a3_bhat=a3a[bkey[0]][bkey[1]],
        ))
    return rows


def main(argv=None):
    import argparse, glob, json
    from fluxrep import read_wout
    ap = argparse.ArgumentParser()
    ap.add_argument("wouts", nargs="+")
    ap.add_argument("--out", default="adaptive_frontier.json")
    args = ap.parse_args(argv)
    # representative configs: a small angular grid (where B-hat is angular-limited)
    # up to the production grid, at fixed radial node counts.
    configs = [(24, 48, 96), (24, 64, 144), (24, 96, 288), (16, 48, 96), (12, 48, 96)]
    paths = [p for pat in args.wouts for p in (glob.glob(pat) or [pat])]
    allres = []
    for p in paths:
        wout = read_wout(p); spec = SpectralRep(wout)
        rows = frontier_3d(spec, wout["nfp"], configs)
        name = p.split("quasr")[-1].split(".")[0]
        allres.append(dict(device=name, nfp=wout["nfp"], rows=rows))
        print(f"\n=== {name} (nfp={wout['nfp']}): uniform / radial-adaptive / FULL-3D-adaptive ===")
        print(f"{'config':>14} {'MiB':>6} | {'sqrtg u/r/3D':>26} | {'Bhat u/r/3D (deg)':>24}")
        for r in rows:
            print(f"{str(tuple(r['config'])):>14} {r['MiB']:>6.1f} | "
                  f"{r['uni_sqrtg']:.3e} {r['rad_sqrtg']:.3e} {r['a3_sqrtg']:.3e} | "
                  f"{r['uni_bhat']:>7.2f} {r['rad_bhat']:>7.2f} {r['a3_bhat']:>7.2f}")
    json.dump(allres, open(args.out, "w"), indent=2)
    print(f"\nwrote {args.out}")
    return allres


if __name__ == "__main__":
    main()
