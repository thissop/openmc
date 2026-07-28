#!/usr/bin/env python
"""BATCH-2 grow: select ~200 MORE fresh engineering-relevant QUASR devices, pre-audit
them NATIVELY (field_audit + coherence C), keep only field-audit passers, and write
sweep/configs/grow_passers2.json for the ARIES-CS sweep (sweep_out_aries).

Single process, resumable, threaded serial-prefetch (I/O) + sequential audit (CPU,
OMP<=2). NO multiprocessing.Pool -> a timeout-kill leaves no orphan workers.

Exclusions (genuinely fresh): production_scaleup 178 + eta_C_zoo 64 + batch-1 grow
candidate pool 220 + every config_*.json id already present in sweep_out / sweep_out_grow
/ sweep_out_aries.

Run repeatedly (resumes) until it prints 'BATCH2 PREAUDIT COMPLETE'; then the manifest
configs/grow_passers2.json is ready.
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import csv, gzip, json, glob, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import collections
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "python"))
CAT = HERE.parent / "data" / "quasr" / "catalogue.csv.gz"
PROD = HERE / "configs" / "production_scaleup.json"
ZOO = HERE.parent / "shield_opt" / "data" / "eta_C_zoo.csv"
B1POOL = HERE / "sweep_out_grow" / "grow_candidate_pool.csv"
POOL2 = HERE / "sweep_out_grow2" / "grow_candidate_pool2.csv"
PADIR = HERE / "sweep_out_grow2" / "preaudit"
MANI = HERE / "configs" / "grow_passers2.json"
TARGET = 200
N_PER_STRATUM = 13
RNG = np.random.default_rng(20260727)

DEFAULTS = {"particles": 200000, "batches": 20, "streams": "free", "blanket": "baseline",
            "nR": 48, "nphi": 96, "nZ": 48, "scale_retries": 3,
            "lost_particle_max": 50, "vr": False}


def _f(r, k):
    try: return float(r[k])
    except Exception: return None


def excluded_ids():
    ids = set()
    for c in json.loads(PROD.read_text())["configs"]:
        ids.add(int(c["id"]))
    for r in csv.DictReader(open(ZOO)):
        ids.add(int(r["id"]))
    if B1POOL.exists():
        for r in csv.DictReader(open(B1POOL)):
            ids.add(int(r["ID"]))
    for d in ("sweep_out", "sweep_out_grow", "sweep_out_aries"):
        for f in glob.glob(str(HERE / d / "config_*_baseline.json")) + \
                 glob.glob(str(HERE / d / "config_*.json")):
            b = Path(f).stem  # config_<id>_baseline or config_<id>_...
            try:
                ids.add(int(b.split("_")[1]))
            except Exception:
                pass
    return ids


def select():
    if POOL2.exists():
        recs = [(int(r["ID"]), r["class"], int(r["nfp"])) for r in csv.DictReader(open(POOL2))]
        return recs
    excl = excluded_ids()
    rows = list(csv.DictReader(gzip.open(CAT, "rt")))
    for r in rows:
        try:
            r["_id"] = int(float(r["ID"]))
            r["_cls"] = "QA" if int(float(r["helicity"])) == 0 else "QH"
            r["_nfp"] = int(float(r["nfp"]))
        except Exception:
            r["_id"] = None
        r["_qs"] = _f(r, "qs_error"); r["_asp"] = _f(r, "aspect_ratio")
        r["_nchp"] = _f(r, "nc_per_hp"); r["_c2s"] = _f(r, "min_coil2surface_dist")
    rows = [r for r in rows if r["_id"] is not None and r["_qs"] is not None
            and r["_asp"] is not None and r["_id"] not in excl]
    rows = [r for r in rows if 2.5 <= r["_asp"] <= 12.0 and 1 <= r["_nfp"] <= 6
            and (r["_nchp"] or 0) >= 1 and (r["_c2s"] is None or r["_c2s"] > 0.05)]
    print(f"[select] catalogue after cuts + exclusion: {len(rows)} devices", flush=True)
    qs = np.array([r["_qs"] for r in rows])
    edges = np.quantile(qs, np.linspace(0, 1, 6)); edges[-1] += 1e-9
    qbin = lambda q: int(np.clip(np.searchsorted(edges, q, side="right") - 1, 0, 4))
    strata = {}
    for r in rows:
        strata.setdefault((r["_cls"], r["_nfp"], qbin(r["_qs"])), []).append(r)
    picked = []
    for key, g in sorted(strata.items()):
        idx = RNG.choice(len(g), size=min(N_PER_STRATUM, len(g)), replace=False)
        picked.extend(g[i] for i in idx)
    if len(picked) > TARGET:
        keep = sorted(RNG.choice(len(picked), size=TARGET, replace=False))
        picked = [picked[i] for i in keep]
    POOL2.parent.mkdir(parents=True, exist_ok=True)
    with open(POOL2, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ID", "class", "nfp", "qs_error", "aspect_ratio"])
        for r in picked:
            w.writerow([r["_id"], r["_cls"], r["_nfp"], f"{r['_qs']:.4f}", f"{r['_asp']:.3f}"])
    cls = [r["_cls"] for r in picked]
    print(f"[select] batch-2 pool: {len(picked)} (QA {cls.count('QA')}, QH {cls.count('QH')}) "
          f"nfp {dict(sorted(collections.Counter(r['_nfp'] for r in picked).items()))}", flush=True)
    return [(r["_id"], r["_cls"], r["_nfp"]) for r in picked]


def prefetch(recs):
    import quasr_loader as ql
    todo = [ID for ID, _, _ in recs]
    def _fetch(ID):
        try: ql.fetch_serial(ID); return (ID, True)
        except Exception as e: return (ID, f"{type(e).__name__}:{e}")
    ok = 0
    with ThreadPoolExecutor(max_workers=16) as ex:
        for fut in as_completed([ex.submit(_fetch, ID) for ID in todo]):
            ID, r = fut.result(); ok += 1 if r is True else 0
    print(f"[prefetch] serials cached: {ok}/{len(todo)}", flush=True)


def audit_one(ID, cls, nfp):
    outp = PADIR / f"pa_{ID}.json"
    if outp.exists():
        try: return json.loads(outp.read_text())
        except Exception: pass
    import quasr_loader as ql, field_audit as fa, biotsavart_field as bsf, coherence_metrics as cm
    out = dict(ID=ID, cls=cls, nfp=nfp, ok=0, passed=False)
    try:
        dev = ql.load_device(ID)
        fld = bsf.CoilField(dev.coils, meta=dev.meta)
        passed, audit = fa.audit_coilfield(fld, dev)
        out.update(ok=1, passed=bool(passed),
                   symmetry_class=dev.meta.get("symmetry_class"), aspect=dev.meta.get("aspect"),
                   boundary_flux=float([c for c in audit["checks"] if c["name"] == "boundary_flux"][0]["value"]))
        if passed:
            xyz, w, rho = dev.source_sample(); bhat = fld.bhat(xyz)
            m = cm.coherence_metrics(bhat, pos=xyz, weights=w, frame="cylindrical")
            out.update(C=float(m["C"]), S_phi=float(m["S_phi"]))
    except Exception as e:  # noqa
        out.update(ok=0, err=f"{type(e).__name__}: {e}")
    PADIR.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(out))
    return out


def build_manifest():
    passers = []
    for f in sorted(glob.glob(str(PADIR / "pa_*.json"))):
        d = json.load(open(f))
        if d.get("ok") and d.get("passed"):
            passers.append(d)
    passers.sort(key=lambda d: d.get("C", 0.0))
    configs = [dict(id=d["ID"], symmetry_class=d.get("symmetry_class", ""),
                    C_preaudit=round(d.get("C", 0.0), 4)) for d in passers]
    MANI.write_text(json.dumps(dict(_comment="batch-2 grow passers, ARIES-CS scale, production fidelity",
                                    _n=len(configs), defaults=DEFAULTS, configs=configs), indent=2))
    return configs


def main():
    recs = select()
    PADIR.mkdir(parents=True, exist_ok=True)
    remaining = [(ID, c, n) for (ID, c, n) in recs if not (PADIR / f"pa_{ID}.json").exists()]
    print(f"[audit] {len(recs)} candidates, {len(remaining)} remaining to audit", flush=True)
    if remaining:
        prefetch(remaining)
    t0 = time.time()
    for i, (ID, cls, nfp) in enumerate(remaining):
        o = audit_one(ID, cls, nfp)
        tag = "PASS" if o.get("passed") else ("FAIL" if o.get("ok") else "ERR")
        if i % 10 == 0 or tag != "PASS":
            print(f"[audit] {i+1}/{len(remaining)} {ID} {cls} {tag} "
                  f"C={o.get('C','-')} bnd={round(o.get('boundary_flux',0),4) if o.get('ok') else '-'}", flush=True)
    # summary
    recs_done = [json.load(open(f)) for f in glob.glob(str(PADIR / "pa_*.json"))]
    ok = [r for r in recs_done if r["ok"]]; p = [r for r in ok if r["passed"]]
    configs = build_manifest()
    if len(recs_done) >= len(recs):
        Cs = np.array([r["C"] for r in p if "C" in r])
        print(f"BATCH2 PREAUDIT COMPLETE: candidates={len(recs)} ok={len(ok)} "
              f"PASS={len(p)} FAIL={len(ok)-len(p)} ERR={len(recs_done)-len(ok)} "
              f"rate={100*len(p)/max(len(ok),1):.1f}% | manifest passers={len(configs)} "
              f"C[{Cs.min():.3f},{Cs.max():.3f}]" if len(p) else "BATCH2 PREAUDIT COMPLETE (0 pass)", flush=True)
    else:
        print(f"[audit] partial: {len(recs_done)}/{len(recs)} audited "
              f"(PASS={len(p)} so far) -- rerun to resume ({time.time()-t0:.0f}s this call)", flush=True)


if __name__ == "__main__":
    main()
