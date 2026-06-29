#!/usr/bin/env python
"""Analytic side of the free-streaming cross-check (config-driven via SPF_CONFIG).

Computes the per-wall-patch neutron wall load by toroidal quadrature of the
filamentary-loop free-streaming kernel, then the directionality A/iso, B/iso.

Two design points that matter when comparing to MONTE CARLO (which traces real rays
through the non-convex square torus), beyond anarrima's own series-vs-quadrature
tests where the arc cancels:

  * VISIBILITY uses the TRUE shaped loop R(phi)=loop_RZ(...) (not the mean radius)
    AND the inner-cylinder occlusion (the central hole blocks the far-side source).
    The crude mean-p, N0>0, cap=2.5 `front_arc` biases the small INBOARD arc and the
    OUTBOARD far side; this finds the real N>0 sub-arcs intersected with "chord
    clears r_in".  Multiple sub-arcs are integrated separately.

  * The quadrature is a fast pure-numpy reimplementation of anarrima's
    free_streaming_quadrature, evaluating loop_RZ/field_bhat DIRECTLY at the
    Gauss-Legendre nodes (identical to what the OpenMC source samples -- no FFT/trig
    intermediary).  It is gated at startup against anarrima itself (asserts a match
    to <1e-4) so anarrima remains the reference; if anarrima is absent the check is
    skipped with a note.  Pure-numpy -> ~1000x faster than eager-JAX anarrima calls.

Run: SPF_CONFIG=quasr_config $HOME/ana_venv/bin/python toy_qalow_analytic.py
     (ana_venv enables the anarrima cross-check; any numpy env runs the computation)
"""
import importlib
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, "/Users/tkiker/Documents/GitHub/openmc/spf_prototype/python")
C = importlib.import_module(os.environ.get("SPF_CONFIG", "toy_qalow_config"))

N_EXACT = int(os.environ.get("SPF_N_EXACT", "384"))   # quadrature nodes (env override)
N_PW = int(os.environ.get("SPF_N_PW", "40"))          # wall patches/wall (env override)
N_SCAN = 1440                                          # visibility scan resolution
_DPH = 1e-5                                            # finite-diff for arc-length J
_GL = {}                                               # cache leggauss nodes/weights


def _nodes(n):
    if n not in _GL:
        _GL[n] = leggauss(n)
    return _GL[n]


def _loop(th, phi, rho):
    """R, Z, dR/dphi, dZ/dphi at angles phi (numpy, exact from the config)."""
    R, Z = C.loop_RZ(th, phi, rho)
    Rp, Zp = C.loop_RZ(th, phi + _DPH, rho)
    Rm, Zm = C.loop_RZ(th, phi - _DPH, rho)
    return R, Z, (Rp - Rm) / (2 * _DPH), (Zp - Zm) / (2 * _DPH)


def quad(th, rho, r, psi, zw, a, b, mode, n):
    """Free-streaming reduced intensity of one loop onto a wall point, over [a,b].
    Polarized modes use normalize=True kernels (A x3/2, B x2)."""
    x, w = _nodes(n)
    half = 0.5 * (b - a); phi = 0.5 * (b + a) + half * x
    R, Z, Rd, Zd = _loop(th, phi, rho)
    J = np.sqrt(R ** 2 + Rd ** 2 + Zd ** 2)
    Dx = R * np.cos(phi) - r; Dy = R * np.sin(phi); Dz = Z - zw
    D = Dx ** 2 + Dy ** 2 + Dz ** 2
    N = np.cos(psi) * Dx + np.sin(psi) * Dz
    base = J * N / D ** 1.5
    if mode == "iso":
        integ = base
    else:
        B = C.field_bhat(th, phi)
        cos2 = (B[:, 0] * Dx + B[:, 1] * Dy + B[:, 2] * Dz) ** 2 / D
        integ = base * (1.0 - cos2) * 1.5 if mode == "A" else base * (0.25 + 0.75 * cos2) * 2.0
    return half * np.sum(w * integ)


def visible_arcs(th, rho, r, psi, zw, r_in):
    """Sub-arcs of phi where the patch sees the source: N>0 on the TRUE shaped loop
    AND the chord (patch->source) clears the inner cylinder r_in (central-hole
    occlusion). Returns a list of (phim, phip); handles the +-pi wrap."""
    phi = np.linspace(-np.pi, np.pi, N_SCAN, endpoint=False)
    R, Z = C.loop_RZ(th, phi, rho)
    Dx = R * np.cos(phi) - r; Dy = R * np.sin(phi); Dz = Z - zw
    N = np.cos(psi) * Dx + np.sin(psi) * Dz
    L2 = Dx * Dx + Dy * Dy
    t = np.clip(-(r * Dx) / np.maximum(L2, 1e-30), 0.0, 1.0)
    mind = np.hypot(r + t * Dx, t * Dy)
    vis = (N > 1e-12) & (mind >= r_in - 1e-9)
    if not vis.any():
        return []
    if vis.all():
        return [(-np.pi, np.pi)]
    dphi = phi[1] - phi[0]
    arcs = []; i = 0
    while i < N_SCAN:
        if vis[i]:
            j = i
            while j < N_SCAN and vis[j]:
                j += 1
            arcs.append([phi[i] - 0.5 * dphi, phi[j - 1] + 0.5 * dphi])
            i = j
        else:
            i += 1
    if len(arcs) >= 2 and vis[0] and vis[-1]:      # merge across the wrap
        first = arcs.pop(0); last = arcs.pop(-1)
        arcs.append([last[0], first[1] + 2 * np.pi])
    return [tuple(a) for a in arcs]


def validate_against_anarrima():
    """Gate: numpy quad must equal anarrima free_streaming_quadrature (same arc)."""
    try:
        import jax
        jax.config.update("jax_enable_x64", True)
        import jax.numpy as jnp
        from anarrima.ripple import device as DV
        from anarrima.ripple import trig as T
    except Exception as e:
        print(f"[validate] anarrima unavailable ({str(e)[:40]}); skipping cross-check")
        return
    K = getattr(C, "K_TRIG", 12)
    pg = np.linspace(0, 2 * np.pi, 128, endpoint=False)

    def to_trig(v):
        Nn = len(v); X = np.fft.rfft(v) / Nn
        re = np.zeros(2 * K + 1); im = np.zeros(2 * K + 1)
        for k in range(1, min(K, Nn // 2) + 1):
            re[K + k] = re[K - k] = X[k].real
            im[K + k] = X[k].imag; im[K - k] = -X[k].imag
        return float(X[0].real), (jnp.asarray(re), jnp.asarray(im))

    maxrel = 0.0
    r, psi, zw = C.R_OUT, np.pi, 0.05
    for rho, th, _ in C.loops()[:4]:
        R, Z = C.loop_RZ(th, pg, rho)
        p, dR = to_trig(R); zc, dZ = to_trig(Z)
        Bfn = (lambda ph: tuple(jnp.asarray(c) for c in C.field_bhat(th, np.asarray(ph)).T))
        for mode in ("iso", "A", "B"):
            mine = quad(th, rho, r, psi, zw, -1.0, 1.0, mode, 256)
            ana = float(DV.free_streaming_quadrature(
                p=p, z=zc - zw, r=r, psi=psi, dR=dR, dZ=dZ,
                Bhat_fn=(None if mode == "iso" else Bfn),
                phim=-1.0, phip=1.0, mode=mode, n_nodes=256, normalize=(mode != "iso")))
            maxrel = max(maxrel, abs(mine - ana) / (abs(ana) + 1e-30))
    print(f"[validate] numpy quad vs anarrima max rel = {maxrel:.2e}")
    assert maxrel < 1e-4, "numpy free-streaming quad disagrees with anarrima!"


def main():
    validate_against_anarrima()
    patches = C.wall_patches(n_per_wall=N_PW)
    modes = ("iso", "A", "B")
    loops = C.loops()
    nwl = {m: np.zeros(len(patches)) for m in modes}
    for ip, pt in enumerate(patches):
        r, psi, zw = pt["R"], pt["psi"], pt["Z"]
        for rho, th, w in loops:
            arcs = visible_arcs(th, rho, r, psi, zw, C.R_IN)
            if not arcs:
                continue
            for m in modes:
                tot = 0.0
                for a, b in arcs:
                    if b - a > 1e-3:
                        tot += quad(th, rho, r, psi, zw, a, b, m, N_EXACT)
                nwl[m][ip] += w * tot
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
    for w in ("inboard", "outboard", "floor", "ceiling"):
        msk = walls == w
        print(f"  {w:8s}: A/iso {np.nanmean(A_iso[msk]):.3f}, B/iso {np.nanmean(B_iso[msk]):.3f}")
    print(f"saved /tmp/{C.STEM}_analytic.npz")


if __name__ == "__main__":
    main()
