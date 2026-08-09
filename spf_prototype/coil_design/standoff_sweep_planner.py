#!/usr/bin/env python
"""Pick a CONTROLLED coil-surface-distance (standoff) sweep from a Gil batch manifest.

The Gil archive varies many optimization knobs at once. For a clean neutronics standoff sweep we
want the coil-surface distance `cs` (= our fixed blanket envelope) to be the ONLY thing changing:
same symmetry (QA/QH), same coil count, everything else as constant as the archive allows. This
reads manifest.csv (from batch_gil_to_makegrid), finds the (symmetry, ncoils) family with the widest
usable cs span, and selects `n_points` sets whose cs values are as evenly spread as the archive
offers (nearest-available to evenly-spaced targets across [cs_min, cs_max], deduplicated).

Output: a sweep plan (list of chosen sets with cs, R0, coils file) + a one-line rationale, as CSV
and/or printed. This is a SELECTION tool -- it does not run neutronics; it names which coil sets to
feed the ParaStell->OpenMC pipeline for the standoff campaign.

Usage:  python standoff_sweep_planner.py manifest.csv [--symmetry QA] [--n-points 5] [--out plan.csv]
"""
import argparse
import csv


def _to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load_manifest(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return rows


def _usable(rows):
    """Keep only successfully-converted sets that carry a numeric cs."""
    out = []
    for r in rows:
        if r.get("status") != "ok":
            continue
        cs = _to_float(r.get("cs"))
        if cs is None:
            continue
        rr = dict(r); rr["_cs"] = cs
        out.append(rr)
    return out


def _family_key(r):
    return (r.get("symmetry", ""), str(r.get("ncoils", "")))


def choose_family(rows, symmetry=None):
    """Return (key, members) for the family with the widest cs span (optionally forced symmetry).
    A family = same (symmetry, ncoils). Ties broken by member count then key for determinism."""
    fams = {}
    for r in _usable(rows):
        if symmetry and r.get("symmetry", "").upper() != symmetry.upper():
            continue
        fams.setdefault(_family_key(r), []).append(r)
    if not fams:
        raise ValueError("no usable sets with a numeric cs" + (f" for symmetry={symmetry}" if symmetry else ""))

    def span(members):
        cs = [m["_cs"] for m in members]
        return max(cs) - min(cs)

    best = max(fams.items(), key=lambda kv: (span(kv[1]), len(kv[1]), kv[0]))
    return best[0], best[1]


def select_spread(members, n_points):
    """Pick up to n_points members whose cs values best cover [cs_min, cs_max]: for each evenly
    spaced target, take the nearest-available cs, then dedupe preserving ascending cs order."""
    ms = sorted(members, key=lambda m: m["_cs"])
    cs = [m["_cs"] for m in ms]
    lo, hi = cs[0], cs[-1]
    if n_points >= len(ms) or hi == lo:
        return ms if hi != lo else ms[:1]
    targets = [lo + (hi - lo) * i / (n_points - 1) for i in range(n_points)]
    chosen_idx = []
    for t in targets:
        j = min(range(len(cs)), key=lambda k: abs(cs[k] - t))
        if j not in chosen_idx:
            chosen_idx.append(j)
    return [ms[j] for j in sorted(chosen_idx)]


def plan_sweep(manifest_path, symmetry=None, n_points=5):
    rows = load_manifest(manifest_path)
    key, members = choose_family(rows, symmetry=symmetry)
    chosen = select_spread(members, n_points)
    cs = [m["_cs"] for m in members]
    rationale = (f"family symmetry={key[0]} ncoils={key[1]}: {len(members)} sets, "
                 f"cs in [{min(cs):.4g}, {max(cs):.4g}]; picked {len(chosen)} spanning the range")
    plan = [dict(order=i, symmetry=m.get("symmetry", ""), ncoils=m.get("ncoils", ""),
                 cs=m["_cs"], measured_R0_m=m.get("measured_R0_m", ""),
                 out_coils=m.get("out_coils", ""), source_json=m.get("source_json", ""))
            for i, m in enumerate(chosen)]
    return plan, rationale


def write_plan(plan, out_path):
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(plan[0].keys()))
        w.writeheader(); w.writerows(plan)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--symmetry", default=None, help="force QA/QH/QI (else widest-span family wins)")
    ap.add_argument("--n-points", type=int, default=5)
    ap.add_argument("--out", default=None, help="write the plan to this CSV")
    a = ap.parse_args()
    plan, rationale = plan_sweep(a.manifest, symmetry=a.symmetry, n_points=a.n_points)
    print(rationale)
    for p in plan:
        print(f"  [{p['order']}] {p['symmetry']} ncoils={p['ncoils']}  cs={p['cs']:.4g}  "
              f"R0={p['measured_R0_m']}  coils={p['out_coils']}")
    if a.out:
        write_plan(plan, a.out); print("wrote", a.out)
