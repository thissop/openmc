#!/usr/bin/env python
"""Exact free-streaming NWL on the REAL precise_QA via anarrima's quadrature.

Goal: the analytic free-streaming benchmark the OpenMC stellarator run must match,
on the SAME plasma OpenMC uses (precise_QA) -- and to show whether anarrima's
TRUNCATED SERIES applies there (Phase 1 found eps_eff~2.3 >> 1, so it should not).

Convention-safe: loop dR,dZ are built by FFT of the real DESC surface loops (all
harmonics, no boundary-basis guesswork); the field is sampled from the same field
map OpenMC uses. Exact = free_streaming_quadrature(384 nodes); hybrid = 192.
Run with the anarrima venv: $HOME/ana_venv/bin/python /tmp/ana_precise_qa.py
"""
import sys
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/python")
from fieldmap import FieldMapField
from anarrima.ripple import device as DV
from anarrima.ripple import trig as T

DATA = "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/data/equil_precise_qa"
K = 24                      # trig order (precise_QA has n up to 8 -> 8*Nfp=16 < K)
WALL_GAP = 0.12            # fraction of R0
N_PATCH = 24
N_EXACT, N_HYB = 384, 192


def loop_to_trig(vals, K):
    """Real periodic samples over [0,2pi) -> (mean, (re,im) trig poly of vals-mean)."""
    N = len(vals)
    X = np.fft.rfft(vals) / N
    re = np.zeros(2 * K + 1); im = np.zeros(2 * K + 1)
    # DC (mean) is returned as p and consumed by anarrima as `R = p + dR(phi)`,
    # so the ripple trig dR must NOT carry it (re[K] stays 0).
    for k in range(1, min(K, N // 2) + 1):
        a_k = 2 * X[k].real          # cos coeff
        b_k = -2 * X[k].imag         # sin coeff
        re[K + k] = a_k / 2; re[K - k] = a_k / 2
        im[K + k] = -b_k / 2; im[K - k] = b_k / 2
    return float(X[0].real), (jnp.asarray(re), jnp.asarray(im))


def main():
    fm = FieldMapField(DATA)
    d = np.load(DATA + "_surface.npz")
    R, Z = d["R"], d["Z"]            # (ntheta, nphi) cm, LCFS
    nth, nph = R.shape
    phg = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    R0 = R.mean()

    # ---- build per-loop trig geometry + reconstruction check ----
    loops = []
    recon_err = 0.0
    for it in range(nth):
        p, dR = loop_to_trig(R[it], K)
        zc, dZ = loop_to_trig(Z[it], K)
        # check: p + dR(phigrid) reproduces R[it]
        Rrec = p + np.asarray(T.evaluate(dR, jnp.asarray(phg)))
        recon_err = max(recon_err, np.max(np.abs(Rrec - R[it])))
        loops.append(dict(p=p, zc=zc, dR=dR, dZ=dZ, Rbar=R[it].mean(), Zbar=Z[it].mean()))
    print(f"[check] loop FFT->trig reconstruction max|err| = {recon_err:.2e} cm (want ~0)")

    # mean cross-section centroid (for outward normals / wall)
    Rc = np.mean([lp["Rbar"] for lp in loops]); Zc = np.mean([lp["Zbar"] for lp in loops])

    def Bhat_fn_for(p, zc, dR, dZ):
        def Bfn(phi):
            ph = np.asarray(phi)
            Rl = p + np.asarray(T.evaluate(dR, jnp.asarray(ph)))
            Zl = zc + np.asarray(T.evaluate(dZ, jnp.asarray(ph)))
            xs, ys, zs = Rl * np.cos(ph), Rl * np.sin(ph), Zl
            B = np.array([fm.bhat((xs[i], ys[i], zs[i])) for i in range(len(ph))])
            return (jnp.asarray(B[:, 0]), jnp.asarray(B[:, 1]), jnp.asarray(B[:, 2]))
        return Bfn

    def front_arc(p, z, r, psi, cap=2.5):
        phis = np.linspace(1e-3, np.pi, 400)
        N0 = np.cos(psi) * (p * np.cos(phis) - r) + z * np.sin(psi)
        if N0[0] <= 0:
            return 0.0
        sign = N0 > 0
        idx = np.argmax(~sign) if (~sign).any() else len(phis) - 1
        return float(min(phis[idx], cap))

    # ---- poloidal patch sweep at one toroidal cut; sum LCFS filaments ----
    thetas_w = np.linspace(0, 2 * np.pi, N_PATCH, endpoint=False)
    modes = ("iso", "A", "B")
    out = {m: {"exact": [], "hybrid": []} for m in modes}
    for itw, _ in enumerate(thetas_w):
        # wall patch = mean cross-section point at this poloidal index + gap outward
        Rb, Zb = loops[itw]["Rbar"], loops[itw]["Zbar"]
        nrm = np.hypot(Rb - Rc, Zb - Zc) + 1e-30
        nR, nZ = (Rb - Rc) / nrm, (Zb - Zc) / nrm
        r = Rb + WALL_GAP * R0 * nR
        z_w = Zb + WALL_GAP * R0 * nZ
        psi = np.arctan2(-nZ, -nR)          # inward normal
        acc = {m: {"exact": 0.0, "hybrid": 0.0} for m in modes}
        for lp in loops:
            z = lp["zc"] - z_w
            phiv = front_arc(lp["p"], z, r, psi)
            if phiv <= 0.05:
                continue
            Bfn = Bhat_fn_for(lp["p"], lp["zc"], lp["dR"], lp["dZ"])
            for m in modes:
                bf = None if m == "iso" else Bfn
                for tag, nn in (("exact", N_EXACT), ("hybrid", N_HYB)):
                    v = float(DV.free_streaming_quadrature(
                        p=lp["p"], z=z, r=r, psi=psi, dR=lp["dR"], dZ=lp["dZ"],
                        Bhat_fn=bf, phim=-phiv, phip=phiv, mode=m, n_nodes=nn))
                    acc[m][tag] += v
        for m in modes:
            out[m]["exact"].append(acc[m]["exact"])
            out[m]["hybrid"].append(acc[m]["hybrid"])
    # ---- report exact vs hybrid agreement (ground-truth self-consistency) ----
    print(f"\nprecise_QA free-streaming NWL (LCFS filament source, gap={WALL_GAP}*R0, "
          f"{N_PATCH} poloidal patches), per mode:")
    print(f"{'mode':>5} {'exact(384)':>12} {'hybrid(192)':>12} {'max rel diff':>13}")
    for m in modes:
        ex = np.array(out[m]["exact"]); hy = np.array(out[m]["hybrid"])
        rel = np.max(np.abs(ex - hy)) / (np.max(np.abs(ex)) + 1e-300)
        print(f"{m:>5} {ex.mean():12.4e} {hy.mean():12.4e} {rel:13.2e}")
    # directionality the OpenMC run must reproduce: A/iso and B/iso peak ratios
    exi, exA, exB = (np.array(out[m]["exact"]) for m in modes)
    print(f"\ndirectionality (exact): A/iso in [{(exA/exi).min():.3f},{(exA/exi).max():.3f}], "
          f"B/iso in [{(exB/exi).min():.3f},{(exB/exi).max():.3f}]")
    np.savez("/tmp/ana_precise_qa_nwl.npz", thetas_w=thetas_w,
             iso=exi, A=exA, B=exB)
    print("saved /tmp/ana_precise_qa_nwl.npz")


if __name__ == "__main__":
    main()
