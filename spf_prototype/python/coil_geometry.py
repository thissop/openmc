#!/usr/bin/env python
"""Build homogenized magnet (winding-pack) solids by sweeping a cross-section along the
REAL QUASR coil filament centerlines -- replacing the uniform conformal 'coil' shell with
geometry that follows the actual coils, so we can tally magnet peaking factors.

Ethan (2026-07-13): homogenized HTS bulk + steel casing, circular OR square cross-section,
fidelity need not be high; the one real pitfall is self-intersection ('pinching') when a
finite-build reactor-scale winding pack follows a tightly curved filament.

Frame twist: a Frenet frame flips/twists at inflection points, which would make the swept
tube spiral. We instead use the double-reflection ROTATION-MINIMIZING frame (Wang, Juttler,
Zheng & Liu, ACM TOG 2008) -- it carries the cross-section along with the least possible
twist -- then close the loop by distributing the residual holonomy angle linearly so the
last ring matches the first. Diagnostics report:
  holonomy_deg : unavoidable twist per loop (small => clean sweep)
  max_pinch    : max curvature*radius; > 1 means the inner wall folds through the axis
                 (local self-intersection) -> shrink the winding pack or the reactor size.
Reuses stellarator_geometry's binary-STL writer + watertight/volume checks.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                      # python/ (stellarator_geometry, quasr_geom)
sys.path.insert(0, str(HERE.parent / "sweep"))     # sweep/ (quasr_loader)
import stellarator_geometry as sg                   # noqa: E402
import quasr_loader as ql                           # noqa: E402


def _unit(v, axis=-1):
    return v / (np.linalg.norm(v, axis=axis, keepdims=True) + 1e-30)


def tangents(P):
    """Central-difference unit tangents of a closed polyline P (N,3)."""
    return _unit(np.roll(P, -1, 0) - np.roll(P, 1, 0))


def curvature(P):
    """Discrete curvature ~|dT/ds| at each vertex of a closed polyline (N,)."""
    t = tangents(P)
    dt = np.roll(t, -1, 0) - np.roll(t, 1, 0)
    ds = np.linalg.norm(np.roll(P, -1, 0) - np.roll(P, 1, 0), axis=1) + 1e-30
    return np.linalg.norm(dt, axis=1) / ds


def _reflect_step(P0, P1, t0, t1, r0):
    """One double-reflection RMF step: transport r0 from (P0,t0) to (P1,t1)."""
    v1 = P1 - P0
    c1 = v1 @ v1 + 1e-30
    rL = r0 - (2.0 / c1) * (v1 @ r0) * v1
    tL = t0 - (2.0 / c1) * (v1 @ t0) * v1
    v2 = t1 - tL
    c2 = v2 @ v2 + 1e-30
    return _unit(rL - (2.0 / c2) * (v2 @ rL) * v2)


def rmf_frames(P):
    """Rotation-minimizing frame along a CLOSED polyline, holonomy distributed so
    frame[N] == frame[0]. Returns (r, s, t) each (N,3) and the holonomy angle (deg)."""
    N = P.shape[0]
    t = tangents(P)
    r = np.zeros((N, 3))
    seed = np.array([0.0, 0.0, 1.0])
    if abs(t[0] @ seed) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    r[0] = _unit(seed - (seed @ t[0]) * t[0])
    for i in range(N - 1):
        r[i + 1] = _reflect_step(P[i], P[i + 1], t[i], t[i + 1], r[i])
    # transport once more around the seam to measure the closing defect
    r_close = _reflect_step(P[N - 1], P[0], t[N - 1], t[0], r[N - 1])
    ang = np.arctan2(r_close @ np.cross(t[0], r[0]), r_close @ r[0])  # holonomy
    # distribute -ang linearly about the tangent so the loop closes smoothly
    for i in range(N):
        th = -ang * (i / N)
        r[i] = _unit(np.cos(th) * r[i] + np.sin(th) * np.cross(t[i], r[i]))
    s = _unit(np.cross(t, r))
    return r, s, t, float(np.degrees(ang))


def cross_section(radius, section, m):
    """m unit offsets (m,2) in local (r,s) coordinates tracing the section outline."""
    if section == "circle":
        a = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
        return radius * np.stack([np.cos(a), np.sin(a)], 1)
    if section == "square":
        q = max(1, m // 4)
        e = np.linspace(-1.0, 1.0, q, endpoint=False)
        top = np.stack([e, np.ones_like(e)], 1)
        right = np.stack([np.ones_like(e), -e], 1)
        bot = np.stack([-e, -np.ones_like(e)], 1)
        left = np.stack([-np.ones_like(e), e], 1)
        return radius * np.concatenate([top, right, bot, left], 0)
    raise ValueError(f"section must be circle|square, got {section}")


def build_tube(P, r, s, radius, section="circle", m=16):
    """Closed tube mesh swept along P with frame (r,s). Returns (V, T) watertight."""
    N = P.shape[0]
    off = cross_section(radius, section, m)                       # (M,2)
    M = off.shape[0]
    rings = (P[:, None, :]
             + off[None, :, 0, None] * r[:, None, :]
             + off[None, :, 1, None] * s[:, None, :])             # (N,M,3)
    V = rings.reshape(N * M, 3)
    T = []
    for i in range(N):
        i2 = (i + 1) % N
        for j in range(M):
            j2 = (j + 1) % M
            a, b = i * M + j, i * M + j2
            c, d = i2 * M + j, i2 * M + j2
            T.append([a, b, d]); T.append([a, d, c])
    T = np.asarray(T, int)
    if sg.signed_volume(V, T) < 0:                               # outward normals
        T = T[:, ::-1]
    return V, T


def sweep_device(ID, radius=None, radius_frac=0.06, scale=1.0,
                 section="circle", m=16, n_samples=256):
    """Sweep every coil of a QUASR device. Returns (V, T, radius, diagnostics, meta)."""
    d = ql.load_device(ID, n_samples=n_samples, coils=True)
    coils = [poly * scale for poly, _ in d.coils]
    if radius is None:  # default winding-pack radius = frac of the median coil "size"
        sizes = [np.median(np.linalg.norm(c - c.mean(0), axis=1)) for c in coils]
        radius = radius_frac * float(np.median(sizes))
    allV, allT, voff, diag = [], [], 0, []
    for k, P in enumerate(coils):
        r, s, t, holo = rmf_frames(P)
        pinch = curvature(P) * radius
        V, T = build_tube(P, r, s, radius, section, m)
        allV.append(V); allT.append(T + voff); voff += V.shape[0]
        diag.append(dict(coil=k, holonomy_deg=holo, max_pinch=float(pinch.max()),
                         n_pinch=int((pinch > 1.0).sum()),
                         manifold=bool(sg.is_edge_manifold(T)),
                         volume=float(sg.signed_volume(V, T))))
    return np.vstack(allV), np.vstack(allT), radius, diag, d.meta


def preview_png(coils_V, path, title=""):
    """Quick 3-D wireframe-ish preview so twist/pinch is visually obvious."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    for V in coils_V:
        ax.plot_trisurf  # noqa (kept import-light: scatter the surface points)
        ax.scatter(V[:, 0], V[:, 1], V[:, 2], s=0.3, c="k", alpha=0.3)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    if title:
        ax.set_title(title)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ID", type=int)
    ap.add_argument("--section", choices=["circle", "square"], default="circle")
    ap.add_argument("--radius", type=float, default=None,
                    help="winding-pack outer radius in OUTPUT units; default = frac of coil size")
    ap.add_argument("--radius_frac", type=float, default=0.06)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="native->output scale (e.g. 1000 for cm at a 10 m reactor)")
    ap.add_argument("--m", type=int, default=16, help="cross-section vertices")
    ap.add_argument("--out", default=None, help="STL path (default data/quasr/coil_geom/<ID>.stl)")
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()

    V, T, radius, diag, meta = sweep_device(
        args.ID, radius=args.radius, radius_frac=args.radius_frac,
        scale=args.scale, section=args.section, m=args.m)

    outdir = HERE.parent / "data" / "quasr" / "coil_geom"
    outdir.mkdir(parents=True, exist_ok=True)
    out = Path(args.out) if args.out else outdir / f"{args.ID}_{args.section}.stl"
    sg._write_binary_stl(V, T, out)

    n_coils = len(diag)
    holo = np.array([d["holonomy_deg"] for d in diag])
    pinch = np.array([d["max_pinch"] for d in diag])
    allman = all(d["manifold"] for d in diag)
    print(f"device {args.ID}: class {meta['symmetry_class']} nfp {meta['nfp']} "
          f"aspect {meta['aspect']}  coils {n_coils}")
    print(f"  section {args.section}  radius(output units) {radius:.4f}  scale {args.scale}")
    print(f"  holonomy/twist per loop [deg]: min {holo.min():.1f} med {np.median(holo):.1f} "
          f"max {abs(holo).max():.1f}")
    print(f"  max pinch (curvature*radius): {pinch.max():.3f} "
          f"({'OK <1, no fold' if pinch.max() < 1 else 'PINCH >=1: shrink radius/scale up'})")
    print(f"  watertight (all coils edge-manifold): {allman}")
    print(f"  vertices {V.shape[0]}  triangles {T.shape[0]}  -> {out}")
    if args.preview:
        png = out.with_suffix(".png")
        # split V back per coil for the preview
        per = V.shape[0] // n_coils
        preview_png([V[i * per:(i + 1) * per] for i in range(n_coils)], png,
                    title=f"{args.ID} {args.section} coils")
        print(f"  preview -> {png}")


if __name__ == "__main__":
    main()
