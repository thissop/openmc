#!/usr/bin/env python
"""PER-DEVICE MAGNET-SIDE NEUTRONICS DRIVER (orchestrator).

One command takes a QUASR device ID all the way through the debugged, hand-run QH/QA pipeline,
parameterized and made idempotent/resumable. Stages (each writes an artifact; a stage is
SKIPPED if its artifact exists, so re-running resumes):

  fetch      make_coils.py            QUASR coils -> reactor-scaled MAKEGRID + centroids
  equil      equil_device.py          DESC solve -> fluxmap(.npz/.bin) + VMEC wout (DESC venv)
  dagmc      build_device_dagmc.py    ParaStell 8-layer -> watertight DAGMC h5m   (pstl env)
  cells      step1/step1_cells.py     pymoab coil-cell discovery (cell_id=GLOBAL_ID+1)
  adjoint    fwcadis/adjoint_importance.py  random-ray ADJOINT, --response flat AND kerma
  finalize   (inline, numpy only)     Bosch-Hale contributon C=S*psi_dagger, concentration
                                      (Gini/PR/peaking), peak-coil adjoint, placement gate,
                                      per-device JSON record.
  [shield]   step1/{gen,build,coil_run}  adaptive placed-vs-uniform shield + coil dose (gated
                                      on the kill-shot verdict; --with-shield to enable)

Heavy stages (equil/dagmc/adjoint/shield) run as SLURM jobs, dependency-chained
(afterok). fetch + finalize run inline. Each stage's env is set in its generated sbatch.

Because the stages need DIFFERENT conda envs and SLURM, this driver is meant to run on the
Ginsburg LOGIN NODE. Paths are resolved from --suite (this dir on Ginsburg) and --work
(scratch for per-device artifacts). One JSON per device: <work>/dev<ID>/record.json.

Usage (Ginsburg login node):
  python run_device.py <ID> [--work DIR] [--suite DIR] [--reactor-a 1.704]
                        [--stages fetch,equil,dagmc,cells,adjoint,finalize]
                        [--with-shield] [--dry-run] [--submit]
Without --submit it prints the SLURM plan (dependency graph) and does the inline stages only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHIELD = HERE.parent                       # spf_prototype/shield_opt
CSV = SHIELD / "data" / "engineering_relevant_devices.csv"

# Ginsburg locations (defaults; override with flags)
DEF_SUITE = os.environ.get("SUITE", str(HERE))
DEF_WORK = os.environ.get("MAGNET_SUITE_WORK",
                          os.path.expanduser("~/pstl_test/magnet_suite"))
# equil uses VMEC (vmecpp) -- QUASR-native, converges where DESC continuation fails.
VMEC_PY = os.path.expanduser("~/vmec_venv/bin/python")
DESC_PY = os.path.expanduser("~/desc_venv/bin/python")
PSTL_ENV = "/ginsburg/astro/users/tjk2147/spf_work/envs/pstl"
ADJ_SRC = os.path.expanduser("~/src/GitHub/openmc-spf")
XS = os.path.expanduser("~/openmc_data/endfb-viii.0-hdf5/cross_sections.xml")
CONDA_SH = "/burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh"


# --------------------------------------------------------------------------- #
# device descriptors
# --------------------------------------------------------------------------- #
def device_row(ID):
    with open(CSV) as f:
        for row in csv.DictReader(f):
            if int(float(row["ID"])) == int(ID):
                return row
    return None


# --------------------------------------------------------------------------- #
# sbatch generation + submission
# --------------------------------------------------------------------------- #
SB_HEADER = """#!/bin/bash
#SBATCH -A astro
#SBATCH -J {job}
#SBATCH -p {part}
#SBATCH -t {time}
#SBATCH -N 1
#SBATCH -c {cpus}
#SBATCH --mem-per-cpu={mem}
#SBATCH -o {work}/%x_%j.out
#SBATCH -e {work}/%x_%j.err
set -uo pipefail
source {conda_sh}
"""


def write_sbatch(path, job, body, work, part="short", time="0-02:00", cpus=16, mem="4G"):
    txt = SB_HEADER.format(job=job, part=part, time=time, cpus=cpus, mem=mem,
                           work=work, conda_sh=CONDA_SH) + body + "\n"
    Path(path).write_text(txt)
    return path


def submit(sbatch, dep=None, dry=False):
    cmd = ["sbatch", "--parsable"]
    if dep:
        cmd += [f"--dependency=afterok:{dep}"]
    cmd += [str(sbatch)]
    if dry:
        print("  [dry] " + " ".join(cmd))
        return "DRYJOBID"
    jid = subprocess.check_output(cmd, text=True).strip().split(";")[0]
    print(f"  submitted {Path(sbatch).name} -> job {jid}"
          + (f" (afterok:{dep})" if dep else ""))
    return jid


# --------------------------------------------------------------------------- #
# inline finalize: contributon + concentration + record (numpy only, no openmc)
# --------------------------------------------------------------------------- #
def finalize(ID, work, suite, row):
    sys.path.insert(0, str(SHIELD))
    sys.path.insert(0, str(HERE))
    import numpy as np
    from adjoint_placement import contributon, placement_priority, importance_angles, coil_angles
    import concentration as conc

    d = Path(work) / f"dev{ID}"
    label = f"quasr{ID}"
    fluxmap = d / f"{label}_fluxmap.npz"
    rec = dict(ID=int(ID), label=label,
               nfp=int(float(row["nfp"])) if row else None,
               aspect=float(row["aspect"]) if row else None,
               qs_class=row["class"] if row else None,
               d_min_over_a=float(row["d_min_over_a"]) if row else None,
               gap_m=float(row["gap_m_ARIESCS"]) if row else None,
               blanket_fit_1p29=int(row["blanket_fit_1.29m"]) if row else None,
               gate_status={}, stages={})

    # ---- collect the pre-flight GATE statuses (GO / NO-GO-with-reason per stage) ----
    eq_status = d / f"{label}_status.json"
    if eq_status.exists():
        rec["gate_status"]["equil"] = json.loads(eq_status.read_text()).get("status")
    wt_status = d / f"dagmc_{label}_watertight.json"
    if wt_status.exists():
        wt = json.loads(wt_status.read_text())
        rec["gate_status"]["dagmc_watertight"] = wt.get("status")
        rec["dagmc_lost_frac"] = wt.get("lost_frac")

    # per-response adjoint maps produced by the adjoint stage
    cen_npz = d / f"{label}_coil_centroids.npz"
    centroids = np.load(cen_npz)["centroids"] if cen_npz.exists() else None
    for response in ("flat", "kerma"):
        amap = d / "adjoint" / f"adjoint_importance_{response}_P0.npz"
        if not amap.exists():
            rec["stages"][f"adjoint_{response}"] = "missing"
            continue
        entry = {}
        # Bosch-Hale contributon C = S(r) * psi_dagger(r) (reactivity-consistent w/ transport)
        if fluxmap.exists():
            C = contributon(str(amap), str(fluxmap), scale=100.0, emissivity="bosch_hale")
            imp = C["importance"]
            entry["emissivity"] = "bosch_hale"
        else:
            imp = np.load(amap)["importance"]
            entry["emissivity"] = "none(bare_adjoint)"
        vals = imp[imp > 0]
        pr = conc.participation_ratio(vals)
        pr, pr_norm = (pr if isinstance(pr, tuple) else (pr, None))
        entry["gini"] = float(conc.gini(vals))
        entry["participation_ratio"] = float(pr)
        entry["participation_ratio_norm"] = (float(pr_norm) if pr_norm is not None else None)
        entry["peaking"] = float(conc.peaking(vals))
        entry["nnz_frac"] = float((imp > 0).sum() / imp.size)
        # median rel-err from the adjoint statepoint export
        a = np.load(amap)
        if "err" in a and "flux" in a:
            fl, er = a["flux"], a["err"]
            rel = er[fl > 0] / fl[fl > 0]
            entry["median_relerr"] = float(np.median(rel)) if rel.size else None
        # placement priority + peak-coil match (structure gate proxy)
        field = C if fluxmap.exists() else str(amap)
        tor = np.linspace(0, 360, 48, endpoint=False)
        pol = np.linspace(0, 360, 48, endpoint=False)
        res = placement_priority(field, tor, pol)
        entry["placement_phi0_deg"] = res["phi0_deg"]
        entry["placement_theta0_deg"] = res["theta0_deg"]
        entry["R0_cm"] = res["R0"]
        if centroids is not None:
            # peak-coil = coil whose (phi) is nearest the placement peak
            dphis = []
            for c in centroids:
                cphi, _ = coil_angles(c, res["R0"])
                dphis.append(abs(((res["phi0_deg"] - cphi + 180) % 360) - 180))
            entry["peak_coil_index"] = int(np.argmin(dphis))
            entry["peak_coil_dphi_deg"] = float(min(dphis))
            entry["placement_points_at_coil"] = bool(min(dphis) < 45.0)
        # structure gate (localized + structured + converged)
        entry["gate_localized"] = entry["nnz_frac"] < 0.98
        entry["gate_structured"] = entry["peaking"] > 1.5
        entry["gate_relerr_ok"] = (entry.get("median_relerr") or 1.0) < 0.25
        entry["gate_pass"] = bool(entry["gate_localized"] and entry["gate_structured"]
                                  and entry["gate_relerr_ok"])
        rec["stages"][f"adjoint_{response}"] = entry

    # emissivity provenance for transport (item 3): the BH emission profile the source used
    try:
        from bosch_hale_emission import bosch_hale_emission
        if fluxmap.exists():
            rho, dens = bosch_hale_emission(str(fluxmap), "bosch_hale")
            rec["transport_emissivity"] = dict(kind="bosch_hale", mean_rho_birth=None,
                                               n_rho=int(len(rho)))
    except Exception as e:
        rec["transport_emissivity_error"] = str(e)

    # ---- OVERALL classification: GO only if every gate passed and >=1 adjoint gate passed ----
    reasons = []
    for stg, st in rec["gate_status"].items():
        if st and st != "go":
            reasons.append(f"{stg}:{st}")
    adj_entries = [v for k, v in rec["stages"].items()
                   if k.startswith("adjoint_") and isinstance(v, dict)]
    if adj_entries and not any(v.get("gate_pass") for v in adj_entries):
        reasons.append("adjoint:structure_gate_fail")
    if not adj_entries:
        reasons.append("adjoint:missing")
    rec["overall_status"] = "go" if not reasons else "no_go:" + ";".join(reasons)

    out = d / "record.json"
    out.write_text(json.dumps(rec, indent=2))
    print(f"[finalize] wrote {out}  OVERALL={rec['overall_status']}")
    print(json.dumps({k: rec[k] for k in ("ID", "nfp", "aspect", "qs_class")}, indent=0))
    for k, v in rec["stages"].items():
        if isinstance(v, dict):
            print(f"  {k}: gini={v.get('gini'):.3f} PR={v.get('participation_ratio'):.3f} "
                  f"peak={v.get('peaking'):.2f} relerr={v.get('median_relerr')} "
                  f"gate={'PASS' if v.get('gate_pass') else 'FAIL'} "
                  f"@coil{v.get('peak_coil_index')}(dphi {v.get('peak_coil_dphi_deg')})")
        else:
            print(f"  {k}: {v}")
    return rec


# --------------------------------------------------------------------------- #
# stage bodies (sbatch)
# --------------------------------------------------------------------------- #
def plan_and_run(args):
    ID = args.ID
    row = device_row(ID)
    if row is None:
        print(f"WARNING: device {ID} not in {CSV.name}; descriptors will be null")
    work = Path(args.work); d = work / f"dev{ID}"; d.mkdir(parents=True, exist_ok=True)
    (d / "adjoint").mkdir(exist_ok=True)
    suite = args.suite
    label = f"quasr{ID}"
    stages = args.stages.split(",")
    dry = args.dry_run or not args.submit
    print(f"=== device {ID} ({row['class'] if row else '?'} nfp"
          f"{row['nfp'] if row else '?'} aspect {row['aspect'] if row else '?'}) ===")
    print(f"work={d}  suite={suite}  stages={stages}  submit={args.submit}")

    scale_file = d / f"{label}_reactor_scale.txt"
    coils_file = d / f"{label}.coils"
    wout = d / f"wout_{label}.nc"
    dagmc = d / f"dagmc_{label}.h5m"
    cells = d / f"{label}_cells.npz"
    cen = d / f"{label}_coil_centroids.npz"

    jid = None  # dependency handle

    # ---- fetch (coils) FIRST: inline on the login node (needs QUASR internet); writes the
    #      shared reactor-scale file that equil then reads, so both share ONE scale factor ----
    if "fetch" in stages and not coils_file.exists():
        cmd = [sys.executable, str(HERE / "make_coils.py"), str(ID), str(d),
               "--reactor-a", str(args.reactor_a)]
        if dry:
            print("  [inline] fetch: " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"  fetch FAILED: {e}")
    else:
        print("  fetch: artifact present or not requested")

    # equil (SLURM) -------------------------------------------------------
    if "equil" in stages and not (wout.exists() and (d / f"{label}_fluxmap.meta").exists()):
        body = (f"{VMEC_PY} {suite}/equil_device.py {ID} {d} "
                f"--reactor-a {args.reactor_a} --solver vmec")
        sb = write_sbatch(d / "sb_equil.sh", f"eq{ID}", body, str(d),
                          time="0-01:00", cpus=8, mem="8G")
        jid = submit(sb, dep=None, dry=dry)
    else:
        print("  equil: artifact present or not requested")

    def _dep(j):   # None if resuming (prior artifact present) or dry
        return j if (j and j != "DRYJOBID") else None

    watertight_json = d / f"dagmc_{label}_watertight.json"

    # dagmc ---------------------------------------------------------------
    jid_d = None
    if "dagmc" in stages and not dagmc.exists():
        body = (f"conda activate {PSTL_ENV}\nexport PATH=\"$CONDA_PREFIX/bin:$PATH\"\n"
                f"python {suite}/build_device_dagmc.py --wout {wout} --coils {coils_file} "
                f"--nfp {row['nfp'] if row else 0} --outname {label} --export-dir {d}")
        sb = write_sbatch(d / "sb_dagmc.sh", f"dg{ID}", body, str(d),
                          time="0-04:00", cpus=16, mem="4G")
        jid_d = submit(sb, dep=_dep(jid), dry=dry)
    else:
        print("  dagmc: artifact present or not requested")

    # cells ---------------------------------------------------------------
    jid_c = None
    if "cells" in stages and not cells.exists():
        body = (f"conda activate spf-stellarator\n"
                f"python {SHIELD}/step1/step1_cells.py {dagmc} {cells}")
        sb = write_sbatch(d / "sb_cells.sh", f"cl{ID}", body, str(d),
                          time="0-00:30", cpus=4, mem="4G")
        jid_c = submit(sb, dep=_dep(jid_d), dry=dry)
    else:
        print("  cells: artifact present or not requested")

    # WATERTIGHTNESS PRE-FLIGHT GATE (cheap; BEFORE the expensive adjoint) ----
    # A leaky DAGMC exits nonzero here, so the adjoint's afterok dependency auto-cancels it --
    # the device is rejected at the right gate with a reason (dagmc_leaky) instead of burning
    # cluster time in transport and MPI_ABORTing on "Maximum number of lost particles reached".
    jid_w = None
    if "watertight" in stages and not watertight_json.exists():
        body = (
            f"conda activate spf-stellarator\n"
            f"export PYTHONPATH={ADJ_SRC}:$PYTHONPATH\n"
            f"export PATH={ADJ_SRC}/build/bin:$PATH\n"
            f"export LD_PRELOAD={ADJ_SRC}/build/lib/libopenmc.so\n"
            f"export OPENMC_CROSS_SECTIONS={XS}\n"
            f"python {suite}/dagmc_watertight_check.py {dagmc} "
            f"--particles 10000 --max-lost 20 --rel-max-lost 1e-3 --out {watertight_json}")
        sb = write_sbatch(d / "sb_watertight.sh", f"wt{ID}", body, str(d),
                          time="0-00:30", cpus=8, mem="4G")
        jid_w = submit(sb, dep=_dep(jid_d), dry=dry)
    elif watertight_json.exists():
        st = json.loads(watertight_json.read_text()).get("status")
        print(f"  watertight: artifact present (status={st})")
        if st and st != "go":
            print(f"  -> NO-GO: skipping adjoint (dagmc gate = {st})")
            return
    else:
        print("  watertight: not requested")

    # adjoint (flat + kerma) -- depends on BOTH cells AND watertight passing ----
    if "adjoint" in stages:
        adj_dep = ":".join(str(x) for x in (jid_c, jid_w) if x and x != "DRYJOBID") or None
        for response in ("flat", "kerma"):
            outmap = d / "adjoint" / f"adjoint_importance_{response}_P0.npz"
            if outmap.exists():
                print(f"  adjoint {response}: artifact present"); continue
            body = (
                f"conda activate spf-stellarator\n"
                # belt-and-suspenders: refuse to run transport if the watertight gate is not GO
                # (afterok already guards, but re-check the status file so a stale/leaky DAGMC is
                # never transported even on a manual/resumed submit).
                f"WT=$(python -c \"import json;print(json.load(open('{watertight_json}'))"
                f"['status'])\" 2>/dev/null || echo missing)\n"
                f"if [ \"$WT\" != go ]; then echo \"NO-GO: watertight=$WT; skipping adjoint\"; "
                f"exit 0; fi\n"
                f"export PYTHONPATH={ADJ_SRC}:$PYTHONPATH\n"
                f"export PATH={ADJ_SRC}/build/bin:$PATH\n"
                f"export LD_PRELOAD={ADJ_SRC}/build/lib/libopenmc.so\n"
                f"export OPENMC_CROSS_SECTIONS={XS}\n"
                f"CELLS=$(python -c \"import numpy as np;print(' '.join(map(str,"
                f"np.load('{cells}')['magnet_cells'])))\")\n"
                f"python {SHIELD}/fwcadis/adjoint_importance.py --dagmc {dagmc} "
                f"--workdir {d}/adjoint --coil-cells $CELLS --coils-file {coils_file} "
                f"--coil-centroids {cen} --response {response} --scatter-order 0")
            sb = write_sbatch(d / f"sb_adjoint_{response}.sh", f"aj{ID}{response[0]}",
                              body, str(d), time="0-06:00", cpus=16, mem="4G")
            submit(sb, dep=adj_dep, dry=dry)

    # finalize (inline) -- wrapped so an unexpected error becomes a clean status, never a crash
    if "finalize" in stages:
        any_map = any((d / "adjoint" / f"adjoint_importance_{r}_P0.npz").exists()
                      for r in ("flat", "kerma"))
        if any_map:
            try:
                finalize(ID, work, suite, row)
            except Exception as e:
                (d / "record.json").write_text(json.dumps(
                    dict(ID=int(ID), overall_status=f"error:finalize:{type(e).__name__}",
                         message=str(e)[:200]), indent=2))
                print(f"  finalize ERROR (recorded, no crash): {type(e).__name__}: {e}")
        else:
            print("  finalize: no adjoint map yet (run after the adjoint stage completes)")

    print("=== plan complete ===")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ID", type=int)
    ap.add_argument("--work", default=DEF_WORK)
    ap.add_argument("--suite", default=DEF_SUITE)
    ap.add_argument("--reactor-a", type=float, default=1.704)
    ap.add_argument("--stages",
                    default="fetch,equil,dagmc,cells,watertight,adjoint,finalize")
    ap.add_argument("--with-shield", action="store_true")
    ap.add_argument("--submit", action="store_true",
                    help="actually submit SLURM jobs (else dry plan + inline stages)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    plan_and_run(args)


if __name__ == "__main__":
    main()
