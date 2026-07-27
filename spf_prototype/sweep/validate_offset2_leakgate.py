#!/usr/bin/env python
"""Definitive before/after using the REAL 3D metric leak_gate.worst_shell_crossing
(trimesh.contains) -- the task's ground truth. Builds each device with the live normal
offset (OLD) and the phase-aligned morphological offset (NEW), scores both, prints a
table. Runs synchronously; leak_gate is slow (~1-2 min/device) but it is the oracle.
"""
from __future__ import annotations
import sys, time, tempfile
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
import stellarator_geometry as sg
import stellarator_geometry_offset2 as sg2
import validate_offset2 as V
import leak_gate

LEAKY = ["1505944", "2190977", "1880092"]
GOOD = ["59509", "115887_ld", "134412_ld", "157482_ld"]
POL_RES = int(sys.argv[sys.argv.index("--pol-res") + 1]) if "--pol-res" in sys.argv else 96


def score(stem, scale, builder, **kw):
    with tempfile.TemporaryDirectory() as td:
        builder(stem, layers=V.BASELINE, scale=scale, outdir=td, **kw)
        try:
            return leak_gate.worst_shell_crossing(td)
        except Exception as e:
            return float("nan")


def main():
    which = [a for a in sys.argv[1:] if not a.startswith("--")
             and a != str(POL_RES)] or (LEAKY + GOOD)
    print(f"pol_res={POL_RES}")
    print(f"{'device':11s} {'class':6s} {'OLD xcross':>10s} {'NEW xcross':>10s}  verdict")
    print("-" * 60)
    for cid in which:
        stem = f"quasr{cid}"
        klass = "leaky" if cid.split('_')[0] in [x.split('_')[0] for x in LEAKY] else "good"
        try:
            scale = V.scale_for(stem)
        except FileNotFoundError:
            print(f"{cid:11s} (no surface npz)"); continue
        t = time.time()
        old = score(stem, scale, sg.build_layers)
        new = score(stem, scale, sg2.build_layers_sdf, pol_res=POL_RES)
        if klass == "leaky":
            verdict = "RECOVERED" if new < 1e-3 else "still leaky"
        else:
            verdict = "OK" if new < 1e-3 else "REGRESSION"
        print(f"{cid:11s} {klass:6s} {old:10.4f} {new:10.4f}  {verdict}  "
              f"({time.time()-t:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
