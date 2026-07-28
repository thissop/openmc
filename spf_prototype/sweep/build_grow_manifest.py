#!/usr/bin/env python
"""Build the worker manifest from the native pre-audit PASSERS only. Production
fidelity (200k x 20, free streams) so the deltas are law-grade (matches the
committed Ginsburg 64). Sorted by C so a partial run still spans the predictor.

In : sweep_out_grow/preaudit/pa_*.json   (passed==True)
Out: configs/grow_passers.json
"""
import glob
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PADIR = HERE / "sweep_out_grow" / "preaudit"
OUT = HERE / "configs" / "grow_passers.json"

DEFAULTS = {
    "particles": 200000, "batches": 20, "streams": "free", "blanket": "baseline",
    "nR": 48, "nphi": 96, "nZ": 48, "scale_retries": 3,
    "lost_particle_max": 50, "vr": False,
}


def main():
    passers = []
    for f in sorted(glob.glob(str(PADIR / "pa_*.json"))):
        d = json.load(open(f))
        if d.get("ok") and d.get("passed"):
            passers.append(d)
    passers.sort(key=lambda d: d.get("C", 0.0))
    configs = [dict(id=d["ID"], symmetry_class=d.get("symmetry_class", ""),
                    C_preaudit=round(d.get("C", 0.0), 4),
                    S_phi_preaudit=round(d.get("S_phi", 0.0), 4))
               for d in passers]
    man = dict(_comment="grow-the-law manifest: native pre-audit passers, production "
               "fidelity. C_preaudit/S_phi_preaudit are the cheap coils-only predictor; "
               "the sweep re-measures C and eta_source under transport.",
               _n=len(configs), defaults=DEFAULTS, configs=configs)
    OUT.write_text(json.dumps(man, indent=2))
    cls = [c["symmetry_class"] for c in configs]
    cs = [c["C_preaudit"] for c in configs]
    print(f"grow manifest: {len(configs)} passers "
          f"(QA {cls.count('QA')}, QH {cls.count('QH')})  "
          f"C[{min(cs):.3f},{max(cs):.3f}]" if configs else "no passers")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
