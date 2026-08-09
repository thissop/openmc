#!/usr/bin/env python
"""Convert Gil et al. 2026 (PRL 137, 065101; Zenodo 18497939) augmented-Lagrangian coil sets to
the MAKEGRID coils file ParaStell's construct_magnets_from_filaments reads (same format as our
coils_qh). Tested end-to-end on the cluster: a QA BiotSavart JSON (20 coils) -> valid MAKEGRID.

The archive stores each coil set as a simsopt BiotSavart JSON (biot_savart_optimized_auglag_*.json).
Load it with simsopt (in ~/vmec_venv), take each coil's curve.gamma() + current, and write MAKEGRID.

NEUTRONICS note: ParaStell only needs the coil GEOMETRY (filament paths + width/thickness), not the
current -- so `--scale` (target_major_radius / archived_major_radius) rescales the ~1 m optimization
coils to reactor size (e.g. ARIES-CS). Current is written but irrelevant to the DAGMC build.

Usage (in the simsopt env):  python gil_to_makegrid.py IN.json OUT_coils [--nfp N] [--scale S]
"""
import argparse
import numpy as np
from simsopt import load


def gil_to_makegrid(json_path, out_path, nfp=1, scale=1.0):
    bs = load(json_path)
    coils = bs.coils
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
    return len(coils)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("json"); ap.add_argument("out")
    ap.add_argument("--nfp", type=int, default=1)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="target_major_radius / archived_major_radius (~1 m coils -> reactor)")
    a = ap.parse_args()
    n = gil_to_makegrid(a.json, a.out, nfp=a.nfp, scale=a.scale)
    print(f"wrote {a.out}: {n} coils (nfp={a.nfp}, scale={a.scale})")
