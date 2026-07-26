"""Build the sweep manifest directly from selected_devices.csv (ID, class, nfp, S_phi, C
already computed by the flat-S_phi down-select). Uses the SAME defaults as
select_configs.build_manifest. We do NOT re-audit here: sweep.py runs field_audit as a
HARD GATE per device, so audit-failing devices are reported as audit_failed records
(a legitimate outcome) rather than silently pre-filtered -- and this avoids a redundant
~15 min Biot-Savart audit pass over all selected devices.

A small blanket-factorization subset (a few spread devices at thin/thick, streams=both)
is appended so analyze_sweep can trace the blanket factor A(tau).
"""
import csv, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SEL = HERE / "sweep_out" / "selected_devices.csv"
OUT = HERE / "configs" / "production_scaleup.json"

DEFAULTS = dict(particles=200000, batches=20, streams="free", blanket="baseline",
                nR=48, nphi=96, nZ=48, scale_retries=3, lost_particle_max=50,
                vr=False, field_source="coil")


def main():
    rows = [r for r in csv.DictReader(open(SEL))]
    for r in rows:
        r["S_phi"] = float(r["S_phi"]); r["C"] = float(r["C"])
    configs = [dict(id=int(r["ID"]), symmetry_class=r["class"],
                    C_predicted=round(r["C"], 4), S_phi_predicted=round(r["S_phi"], 4))
               for r in rows]
    # blanket subset: 3 devices spread across S_phi, each at thin+thick, streams=both
    srt = sorted(rows, key=lambda r: r["S_phi"])
    picks = [srt[int(t)] for t in np.linspace(0, len(srt) - 1, 3)]
    for r in picks:
        for blk in ("thin", "thick"):
            configs.append(dict(id=int(r["ID"]), symmetry_class=r["class"],
                                blanket=blk, streams="both",
                                C_predicted=round(r["C"], 4),
                                S_phi_predicted=round(r["S_phi"], 4),
                                _note="blanket-factorization subset"))
    Cs = np.array([r["C"] for r in rows]); Ss = np.array([r["S_phi"] for r in rows])
    cov = dict(n_scan=len(rows), n_total_configs=len(configs),
               C_range=[float(Cs.min()), float(Cs.max())],
               S_phi_range=[float(Ss.min()), float(Ss.max())],
               QA=sum(r["class"] == "QA" for r in rows),
               QH=sum(r["class"] == "QH" for r in rows),
               by_nfp={int(n): int((np.array([int(r["nfp"]) for r in rows]) == n).sum())
                       for n in sorted({int(r["nfp"]) for r in rows})})
    man = dict(_comment="Scale-up manifest from stratified flat-S_phi selection. "
               "C_predicted/S_phi_predicted are cheap field-only predictions; the sweep "
               "measures eta_source (free-stream coil-flux directional change) to test them.",
               defaults=DEFAULTS, _coverage=cov, configs=configs)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(man, indent=2))
    print(f"wrote {OUT}")
    print(f"  scan configs: {len(rows)}  (+{len(configs)-len(rows)} blanket subset) "
          f"= {len(configs)} total array tasks")
    print(f"  class: QA {cov['QA']} / QH {cov['QH']}   by nfp: {cov['by_nfp']}")
    print(f"  C span [{cov['C_range'][0]:.3f},{cov['C_range'][1]:.3f}]  "
          f"S_phi span [{cov['S_phi_range'][0]:.3f},{cov['S_phi_range'][1]:.3f}]")
    print(f"  submit: sbatch --array=0-{len(configs)-1}%25 sweep.sbatch")


if __name__ == "__main__":
    main()
