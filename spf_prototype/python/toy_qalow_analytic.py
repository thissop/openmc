#!/usr/bin/env python
"""Analytic side of the toy QA-low free-streaming cross-check (run in ana_venv).

Uses anarrima.ripple.free_streaming_quadrature (exact, 384 nodes) for the toy
QA-low source (toy_qalow_config) onto the shared square-torus wall patches, per
mode. Saves per-patch NWL + directionality (A/iso, B/iso) to compare with OpenMC.

Run: $HOME/ana_venv/bin/python toy_qalow_analytic.py
"""
import sys
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/python")
import importlib
import os
C = importlib.import_module(os.environ.get("SPF_CONFIG", "toy_qalow_config"))
from anarrima.ripple import device as DV
from anarrima.ripple import trig as T

K = 12
N_EXACT = 384
N_PW = 40                         # wall patches per wall (fine -> smooth profile line)
PHI_GRID = np.linspace(0, 2 * np.pi, 128, endpoint=False)


def loop_to_trig(vals, K):
    """real periodic samples -> (mean, ripple trig (re,im)); DC stays in mean."""
    N = len(vals); X = np.fft.rfft(vals) / N
    re = np.zeros(2 * K + 1); im = np.zeros(2 * K + 1)
    for k in range(1, min(K, N // 2) + 1):
        a_k, b_k = 2 * X[k].real, -2 * X[k].imag
        re[K + k] = a_k / 2; re[K - k] = a_k / 2
        im[K + k] = -b_k / 2; im[K - k] = b_k / 2
    return float(X[0].real), (jnp.asarray(re), jnp.asarray(im))


def front_arc(p, z, r, psi, cap=2.5):
    phis = np.linspace(1e-3, np.pi, 400)
    N0 = np.cos(psi) * (p * np.cos(phis) - r) + z * np.sin(psi)
    if N0[0] <= 0:
        return 0.0
    sign = N0 > 0
    idx = np.argmax(~sign) if (~sign).any() else len(phis) - 1
    return float(min(phis[idx], cap))


def build_loops():
    K = getattr(C, "K_TRIG", 12)     # trig order; raise for rich (QUASR) spectra
    L = []
    recon = 0.0
    for rho, th, w in C.loops():
        R, Z = C.loop_RZ(th, PHI_GRID, rho)
        p, dR = loop_to_trig(R, K); zc, dZ = loop_to_trig(Z, K)
        recon = max(recon, np.max(np.abs(p + np.asarray(T.evaluate(dR, jnp.asarray(PHI_GRID))) - R)))
        L.append(dict(p=p, zc=zc, dR=dR, dZ=dZ, th=th, w=w))
    print(f"[check] loop reconstruction max|err| = {recon:.2e} (want ~0)")
    return L


def bhat_fn_for(theta_s):
    def Bfn(phi):
        ph = np.asarray(phi)
        B = C.field_bhat(theta_s, ph)
        return (jnp.asarray(B[:, 0]), jnp.asarray(B[:, 1]), jnp.asarray(B[:, 2]))
    return Bfn


def main():
    loops = build_loops()
    patches = C.wall_patches(n_per_wall=N_PW)
    modes = ("iso", "A", "B")
    nwl = {m: np.zeros(len(patches)) for m in modes}
    for ip, pt in enumerate(patches):
        r, psi, Zw = pt["R"], pt["psi"], pt["Z"]
        for lp in loops:
            z = lp["zc"] - Zw
            phiv = front_arc(lp["p"], z, r, psi)
            if phiv <= 0.02:
                continue
            Bfn = bhat_fn_for(lp["th"])
            for m in modes:
                bf = None if m == "iso" else Bfn
                # normalize=True puts the polarized kernels on the unit-emission
                # convention OpenMC samples (A *=3/2, B *=2 exactly; see device.py
                # _KERNEL_MEAN), so A/iso and B/iso are directly comparable.
                v = float(DV.free_streaming_quadrature(
                    p=lp["p"], z=z, r=r, psi=psi, dR=lp["dR"], dZ=lp["dZ"],
                    Bhat_fn=bf, phim=-phiv, phip=phiv, mode=m, n_nodes=N_EXACT,
                    normalize=(m != "iso")))
                nwl[m][ip] += lp["w"] * v
    iso, A, B = nwl["iso"], nwl["A"], nwl["B"]
    good = iso > 1e-9 * np.max(iso)
    A_iso = np.where(good, A / np.where(good, iso, 1), np.nan)
    B_iso = np.where(good, B / np.where(good, iso, 1), np.nan)
    walls = np.array([p["wall"] for p in patches])
    s = np.array([p["s"] for p in patches])
    np.savez(f"/tmp/{C.STEM}_analytic.npz", iso=iso, A=A, B=B,
             A_iso=A_iso, B_iso=B_iso, wall=walls, s=s)
    print(f"analytic A/iso range [{np.nanmin(A_iso):.3f},{np.nanmax(A_iso):.3f}], "
          f"B/iso range [{np.nanmin(B_iso):.3f},{np.nanmax(B_iso):.3f}]")
    for w in ("inboard", "outboard", "floor"):
        msk = walls == w
        print(f"  {w:8s}: A/iso {np.nanmean(A_iso[msk]):.3f}, B/iso {np.nanmean(B_iso[msk]):.3f}")
    print(f"saved /tmp/{C.STEM}_analytic.npz")


if __name__ == "__main__":
    main()
