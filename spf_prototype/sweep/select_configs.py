#!/usr/bin/env python
"""Config-family selection: compute the CHEAP coherence C for candidate QUASR
devices (field + audit only, no neutronics), then tile the C axis roughly uniformly
and emit a sweep manifest. Compute C FIRST, then choose, so the family spans a real
range of C instead of clustering (per docs/EXPERIMENTAL_DESIGN.md).

Also appends the small blanket-factorization subset (2-3 configs re-run at the
thin/thick blankets with streams='both') so analyze_sweep can trace A(tau).

Usage:
  python select_configs.py --ids 59509 932746 ... --out configs/production.json
  python select_configs.py --auto 20 --out configs/production.json   # pick diverse candidates
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))

import quasr_loader as ql
import coherence_metrics as cm
import field_audit as fa
from biotsavart_field import CoilField

CAT = HERE.parent / "data" / "quasr" / "catalogue.csv.gz"


def catalogue_rows():
    with gzip.open(CAT, "rt") as f:
        return list(csv.DictReader(f))


def auto_candidates(n):
    """Pick a diverse candidate set from the catalogue: spread over nfp and qs_error
    (a proxy for shaping), so the resulting C axis is populated rather than clustered."""
    rows = [r for r in catalogue_rows() if r.get("ID")]
    # bucket by nfp, then within each spread over |qs_error|
    by_nfp = {}
    for r in rows:
        by_nfp.setdefault(int(float(r["nfp"])), []).append(r)
    picks = []
    for nfp, rs in sorted(by_nfp.items()):
        rs = sorted(rs, key=lambda r: abs(float(r["qs_error"])))
        idx = np.linspace(0, len(rs) - 1, max(1, n // len(by_nfp))).astype(int)
        picks += [int(float(rs[i]["ID"])) for i in np.unique(idx)]
    return picks[:n]


def compute_metrics(ID, class_map, include_failed=False, field_source="coil"):
    """Load, audit, compute the cheap predictors (C, S_phi, lambda_phi, reversal_frac).
    Returns a dict (with 'audit_passed') or None if the device can't be loaded, or if
    it failed the audit and include_failed is False. include_failed=True keeps
    audit-failed devices for the COVERAGE report only -- e.g. reversal devices like
    1190023 that the coil-field B.n audit filters (docs/SCOPE_AND_CAVEATS.md).
    field_source='equilibrium' (G3) reads the pre-generated data/quasr<ID>_equil map
    (a flux function, B.n=0 by construction) so those devices pass and enter the family."""
    try:
        dev = ql.load_device(ID)
    except Exception as e:  # noqa
        print(f"  {ID}: load failed ({e})"); return None
    if field_source == "equilibrium":
        import equilibrium_field as ef
        field = ef.EquilibriumField(str(HERE.parent / "data" / f"quasr{ID}_equil"),
                                    meta=dev.meta)
    else:
        field = CoilField(dev.coils, meta=dev.meta)
    passed, audit = fa.audit_coilfield(field, dev)
    bflux = [c for c in audit["checks"] if c["name"] == "boundary_flux"][0]
    if not passed and not include_failed:
        print(f"  {ID}: AUDIT FAILED (B.n rms {bflux['value']:.2e}) -> excluded")
        return None
    xyz, w, _ = dev.source_sample()
    m = cm.coherence_metrics(field.bhat(xyz), pos=xyz, weights=w, frame="cylindrical")
    cls = class_map.get(str(ID), class_map.get(ID, dev.meta["symmetry_class"]))
    flag = "" if passed else "  [AUDIT FAILED -- coverage only]"
    print(f"  {ID}: C={m['C']:.4f} S_phi={m['S_phi']:.4f} rev={m['reversal_frac']*100:.1f}% "
          f"nfp={dev.meta['nfp']} class={cls}  (B.n rms {bflux['value']:.1e}){flag}")
    return dict(id=ID, C=float(m["C"]), S_phi=float(m["S_phi"]),
                lambda_phi=float(m["lambda_phi"]), reversal_frac=float(m["reversal_frac"]),
                nfp=dev.meta["nfp"], iota=dev.meta["iota"], aspect=dev.meta["aspect"],
                symmetry_class=cls, anisotropy=float(m["anisotropy"]),
                audit_passed=bool(passed), bn_rms=float(bflux["value"]))


compute_C = compute_metrics   # backward-compatible alias


def coverage_report(cands, nbins=5):
    """Class x C-bin coverage. Exposes whether symmetry class is ALIASED with C (each
    class in a disjoint C range) -- which confounds a measured class 'split' with mere
    curvature in eta(C) -- and lists reversal_frac>0 devices where C decouples from the
    parity-correct S_phi. Returns a dict and prints a table."""
    if not cands:
        return {}
    Cs = [c["C"] for c in cands]
    lo, hi = min(Cs), max(Cs)
    span = (hi - lo) or 1.0
    edges = np.linspace(lo, hi, nbins + 1)
    classes = sorted({c["symmetry_class"] for c in cands})
    grid = {cl: [0] * nbins for cl in classes}
    for c in cands:
        b = min(nbins - 1, int((c["C"] - lo) / span * nbins))
        grid[c["symmetry_class"]][b] += 1
    ranges = {cl: [min(c["C"] for c in cands if c["symmetry_class"] == cl),
                   max(c["C"] for c in cands if c["symmetry_class"] == cl)]
              for cl in classes}
    cl_list = list(ranges)
    aliased = any(ranges[cl_list[i]][1] < ranges[cl_list[j]][0]
                  or ranges[cl_list[j]][1] < ranges[cl_list[i]][0]
                  for i in range(len(cl_list)) for j in range(i + 1, len(cl_list)))
    reversal_devs = [c["id"] for c in cands if c.get("reversal_frac", 0.0) > 0.01]
    print("\n[coverage] class x C-bin (C in [%.3f, %.3f], %d bins):" % (lo, hi, nbins))
    for cl in classes:
        print(f"    {cl:>3}  " + "  ".join(f"{n:3d}" for n in grid[cl]))
    print(f"[coverage] per-class C ranges: {ranges}")
    print(f"[coverage] class ALIASED with C (disjoint ranges): {aliased}"
          + ("  -> add QA-at-low-C / QH-at-high-C or a QI arm to de-alias"
             if aliased else ""))
    if reversal_devs:
        print(f"[coverage] reversal devices (C decouples from S_phi): {reversal_devs}"
              "  -> need equilibrium field to admit (G3)")
    return dict(c_range=[lo, hi], class_ranges=ranges, class_aliased_with_C=aliased,
                reversal_devices=reversal_devs, grid=grid,
                bin_edges=[float(e) for e in edges])


def tile_C(cands, n_tiles, key="C"):
    """Pick n_tiles configs spread roughly uniformly along the `key` axis (C or S_phi)."""
    cands = sorted(cands, key=lambda c: c[key])
    if len(cands) <= n_tiles:
        return cands
    targets = np.linspace(cands[0][key], cands[-1][key], n_tiles)
    chosen, used = [], set()
    for t in targets:
        i = min(range(len(cands)),
                key=lambda j: (abs(cands[j][key] - t), j in used))
        if i not in used:
            used.add(i); chosen.append(cands[i])
    return chosen


def build_manifest(chosen, n_blanket_subset=2, coverage=None, tile_on="C",
                   field_source="coil"):
    # field_source='coil' (Biot-Savart) is the default; 'equilibrium' (G3) uses a
    # VMEC/Boozer field map so reversal devices pass the audit.
    defaults = dict(particles=200000, batches=20, streams="free", blanket="baseline",
                    nR=48, nphi=96, nZ=48, scale_retries=3, lost_particle_max=50,
                    vr=False, field_source=field_source)
    configs = []
    for c in chosen:
        configs.append(dict(id=c["id"], symmetry_class=c["symmetry_class"],
                            C_predicted=round(c["C"], 4),
                            S_phi_predicted=round(c["S_phi"], 4)))
    # blanket-factorization subset: take a few spread configs, run BOTH streams at
    # baseline/thin/thick to trace A(tau) and test eta_source(C) invariance
    sub = tile_C(chosen, min(n_blanket_subset, len(chosen)), key=tile_on)
    for c in sub:
        for blk in ("baseline", "thin", "thick"):
            configs.append(dict(id=c["id"], symmetry_class=c["symmetry_class"],
                                blanket=blk, streams="both",
                                C_predicted=round(c["C"], 4),
                                S_phi_predicted=round(c["S_phi"], 4),
                                _note="blanket-factorization subset"))
    return dict(
        _comment="Auto-generated by select_configs.py. C_predicted/S_phi_predicted are "
                 "the cheap field-only predictions; the sweep measures eta to test them. "
                 "S_phi (nematic order) is the parity-correct predictor -- see docs/THEORY.md.",
        defaults=defaults, _coverage=coverage or {}, configs=configs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", type=int, nargs="*", default=None)
    ap.add_argument("--auto", type=int, default=None, help="auto-pick N diverse candidates")
    ap.add_argument("--tiles", type=int, default=12, help="configs in the scan")
    ap.add_argument("--tile-on", choices=["C", "S_phi"], default="C",
                    help="axis to tile the family uniformly along (S_phi is parity-correct)")
    ap.add_argument("--include-failed", action="store_true",
                    help="also compute metrics for audit-failed devices (coverage only, "
                         "not transported) -- surfaces reversal devices like 1190023")
    ap.add_argument("--field-source", choices=["coil", "equilibrium"], default="coil",
                    help="b_hat source; 'equilibrium' reads data/quasr<ID>_equil maps (G3, "
                         "flux function) so reversal / poor-coil-fit devices pass the audit")
    ap.add_argument("--class-map", default=None, help="JSON {id: 'QA'|'QH'|'QI'}")
    ap.add_argument("--out", default="configs/production.json")
    args = ap.parse_args()

    class_map = json.loads(Path(args.class_map).read_text()) if args.class_map else {}
    ids = args.ids or (auto_candidates(args.auto) if args.auto else [])
    if not ids:
        ap.error("give --ids or --auto N")
    print(f"[select] computing C/S_phi for {len(ids)} candidates ...")
    cands_all = [c for c in (compute_metrics(i, class_map, include_failed=args.include_failed,
                                             field_source=args.field_source)
                             for i in ids) if c]
    cov = coverage_report(cands_all)                      # coverage over ALL (incl failed)
    passed = [c for c in cands_all if c.get("audit_passed", True)]
    if not passed:
        ap.error("no candidates passed the field audit")
    chosen = tile_C(passed, args.tiles, key=args.tile_on)
    man = build_manifest(chosen, coverage=cov, tile_on=args.tile_on,
                         field_source=args.field_source)
    outp = Path(args.out); outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(man, indent=2))
    xs = sorted(c[args.tile_on] for c in chosen)
    print(f"[select] wrote {outp}: {len(chosen)} scan configs "
          f"({args.tile_on} in [{xs[0]:.3f},{xs[-1]:.3f}]) + blanket subset. "
          f"Set sbatch --array=0-{len(man['configs'])-1}.")


if __name__ == "__main__":
    main()
