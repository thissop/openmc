#!/usr/bin/env python
"""Convert Gil et al. 2026 (PRL 137, 065101; Zenodo 18497939) augmented-Lagrangian coil sets to
the MAKEGRID coils file ParaStell's construct_magnets_from_filaments reads (same format as our
coils_qh). Tested end-to-end on the cluster: a QA BiotSavart JSON (20 coils) -> valid MAKEGRID.

The archive stores each coil set as a simsopt BiotSavart JSON (biot_savart_optimized_auglag_*.json).
Load it with simsopt (in ~/vmec_venv), take each coil's curve.gamma() + current, and write MAKEGRID.

NEUTRONICS note: ParaStell only needs the coil GEOMETRY (filament paths + width/thickness), not the
current -- so `--scale` (target_major_radius / archived_major_radius) rescales the ~1 m optimization
coils to reactor size (e.g. ARIES-CS). Current is written but irrelevant to the DAGMC build.

Usage (in the simsopt env):
  python gil_to_makegrid.py IN.json OUT_coils [--nfp N] [--scale S | --target-major-radius R0_m]
  --target-major-radius reads the set's own R0 (mean cylindrical R of coil centroids) and picks
  the multiplier for you, so "put this set at ARIES-CS scale" is one flag, not a hand computation.
"""
import argparse
import numpy as np
from simsopt import load


def coil_set_major_radius(coils):
    """Major radius of a filament coil set = mean cylindrical R of the coil CENTROIDS (meters).

    Each Gil coil is a closed filament sitting at some toroidal location; its centroid's distance
    from the z-axis is the local major radius, and averaging over coils gives the set's R0. This is
    the quantity a reactor rescale is defined against (target_R0 / measured_R0)."""
    Rs = []
    for c in coils:
        ctr = c.curve.gamma().mean(axis=0)                # (x,y,z) centroid of this filament
        Rs.append(float(np.hypot(ctr[0], ctr[1])))
    return float(np.mean(Rs))


def resolve_scale(coils, scale=None, target_major_radius=None):
    """Turn a user request into a single multiplier. Exactly one of `scale` (raw ratio) or
    `target_major_radius` (meters -> auto ratio) may be given; default is 1.0 (no rescale)."""
    if scale is not None and target_major_radius is not None:
        raise ValueError("give either --scale or --target-major-radius, not both")
    if target_major_radius is not None:
        r0 = coil_set_major_radius(coils)
        if r0 <= 0:
            raise ValueError("measured major radius is non-positive; cannot auto-scale")
        return float(target_major_radius) / r0
    return 1.0 if scale is None else float(scale)


def gil_to_makegrid(json_path, out_path, nfp=1, scale=None, target_major_radius=None):
    bs = load(json_path)
    coils = bs.coils
    scale = resolve_scale(coils, scale=scale, target_major_radius=target_major_radius)
    with open(out_path, "w") as f:
        f.write(f"periods {nfp}\n"); f.write("begin filament\n"); f.write("mirror NIL\n")
        for k, c in enumerate(coils):
            xyz = c.curve.gamma() * float(scale)          # meters (optionally rescaled to reactor)
            I = float(c.current.get_value())
            for x, y, z in xyz:
                f.write(f"  {x: .15E}  {y: .15E}  {z: .15E}  {I: .15E}\n")
            x, y, z = xyz[0]                              # close the loop: current 0 + group label
            f.write(f"  {x: .15E}  {y: .15E}  {z: .15E}  {0.0: .15E} {k+1} Coil{k+1}\n")
        f.write("end\n")
    return len(coils), scale


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("json"); ap.add_argument("out")
    ap.add_argument("--nfp", type=int, default=1)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--scale", type=float, default=None,
                   help="raw multiplier = target_major_radius / archived_major_radius")
    g.add_argument("--target-major-radius", type=float, default=None,
                   help="target R0 in METERS; the multiplier is computed from the set's own R0")
    a = ap.parse_args()
    n, s = gil_to_makegrid(a.json, a.out, nfp=a.nfp, scale=a.scale,
                           target_major_radius=a.target_major_radius)
    print(f"wrote {a.out}: {n} coils (nfp={a.nfp}, scale={s:.6g})")
