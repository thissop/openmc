#!/usr/bin/env python
"""Native (Mac, no VM, no transport) PRE-AUDIT of the fresh grow pool -- the efficiency
step. For each candidate run ONLY the cheap stages:
   quasr_loader.load_device -> biotsavart_field.CoilField
   -> field_audit.audit_coilfield  (the REAL thresholds, HARD gate)
   -> coherence_metrics.coherence_metrics  (C, S_phi)
Devices that PASS the field audit are the only ones worth transporting on the worker
(this skips the ~19% coil-fit audit-failers before they waste worker time). Resumable:
one JSON per device under sweep_out_grow/preaudit/.

Usage: python preaudit_grow.py [nproc]
Out:   sweep_out_grow/preaudit/pa_<ID>.json  and  sweep_out_grow/preaudit_summary.csv
"""
import csv
import json
import sys
import time
from pathlib import Path
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))

POOL = HERE / "sweep_out_grow" / "grow_candidate_pool.csv"
PADIR = HERE / "sweep_out_grow" / "preaudit"
SUMMARY = HERE / "sweep_out_grow" / "preaudit_summary.csv"


def one(rec):
    ID, cls, nfp = rec
    outp = PADIR / f"pa_{ID}.json"
    if outp.exists():
        try:
            return json.loads(outp.read_text())
        except Exception:
            pass
    import numpy as np
    import quasr_loader as ql
    import field_audit as fa
    import biotsavart_field as bsf
    import coherence_metrics as cm
    t = time.time()
    out = dict(ID=ID, cls=cls, nfp=nfp, ok=0, passed=False)
    try:
        dev = ql.load_device(ID)
        fld = bsf.CoilField(dev.coils, meta=dev.meta)
        passed, audit = fa.audit_coilfield(fld, dev)
        checks = {c["name"]: dict(value=float(c["value"]), passed=bool(c["passed"]),
                                  **({"max": float(c["max"])} if "max" in c else {}))
                  for c in audit["checks"]}
        out.update(ok=1, passed=bool(passed), checks=checks,
                   aspect=dev.meta.get("aspect"),
                   symmetry_class=dev.meta.get("symmetry_class"),
                   total_coil_length=dev.meta.get("total_coil_length"))
        if passed:
            xyz, w, rho = dev.source_sample()
            bhat = fld.bhat(xyz)
            m = cm.coherence_metrics(bhat, pos=xyz, weights=w, frame="cylindrical")
            out.update(C=float(m["C"]), S_phi=float(m["S_phi"]),
                       reversal_frac=float(m["reversal_frac"]))
    except Exception as e:  # noqa
        out.update(ok=0, err=f"{type(e).__name__}: {e}")
    out["secs"] = round(time.time() - t, 1)
    PADIR.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(out))
    tag = "PASS" if out.get("passed") else ("FAIL" if out.get("ok") else "ERR")
    print(f"  {ID} {cls} nfp{nfp}: {tag} "
          f"C={out.get('C','-')} bnd={out.get('checks',{}).get('boundary_flux',{}).get('value','-')} "
          f"({out['secs']}s)", flush=True)
    return out


def main(nproc=4):
    recs = [(int(r["ID"]), r["class"], int(r["nfp"]))
            for r in csv.DictReader(open(POOL))]
    print(f"pre-auditing {len(recs)} candidates on {nproc} procs ...", flush=True)
    PADIR.mkdir(parents=True, exist_ok=True)
    with Pool(nproc) as p:
        out = p.map(one, recs)
    ok = [o for o in out if o["ok"]]
    passed = [o for o in ok if o["passed"]]
    with open(SUMMARY, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "ok", "passed", "C", "S_phi",
                    "boundary_flux_rms", "boundary_flux_max", "aspect"])
        for o in out:
            bf = o.get("checks", {}).get("boundary_flux", {})
            w.writerow([o["ID"], o["cls"], o["nfp"], o["ok"], int(bool(o.get("passed"))),
                        o.get("C", ""), o.get("S_phi", ""),
                        bf.get("value", ""), bf.get("max", ""), o.get("aspect", "")])
    print(f"\nPRE-AUDIT: {len(out)} candidates | ok(loaded)={len(ok)} | "
          f"PASS={len(passed)} | FAIL={len(ok)-len(passed)} | ERR={len(out)-len(ok)}")
    if ok:
        rate = 100.0 * len(passed) / len(ok)
        print(f"pre-audit PASS rate = {len(passed)}/{len(ok)} = {rate:.1f}% (of loadable)")
    print("wrote", SUMMARY)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
