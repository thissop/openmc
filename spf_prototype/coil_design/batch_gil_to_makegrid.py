#!/usr/bin/env python
"""Batch-convert the whole Gil et al. 2026 auglag archive to ParaStell-ready MAKEGRID files, and
emit a manifest so every set is neutronics-ready the moment the cluster is back.

The archive (Zenodo 18497939, extracted on the cluster at ~/gil_coils/) stores 86 coil sets, each a
simsopt BiotSavart JSON named biot_savart_optimized_auglag_*.json, filed under directories whose names
encode the optimization config (QA/QH, ncoils, and penalty knobs curvature/msc/force/length/cc/cs --
cc = coil-coil distance, cs = coil-surface distance = our fixed blanket standoff). We:
  * recursively find every biot_savart_optimized_auglag_*.json,
  * parse the config tokens out of the path,
  * convert each to a MAKEGRID coils file (optionally rescaled to a target reactor major radius),
  * measure each set's own R0 and coil count,
  * write a manifest CSV (one row per set) for downstream neutronics selection.

Config-token parsing is heuristic (the archive's own directory scheme) and NEVER blocks a conversion:
an unrecognised path still converts; its unknown tokens land in the manifest 'raw_tokens' column.

Run in the simsopt env (cluster ~/vmec_venv):
  python batch_gil_to_makegrid.py ARCHIVE_ROOT OUT_DIR [--target-major-radius R0_m] [--nfp-map QA=2,QH=4]
"""
import argparse
import csv
import glob
import os
import re

from gil_to_makegrid import gil_to_makegrid, coil_set_major_radius
from simsopt import load


# ---- config-token parsing from the archive directory names ------------------------------------
# tokens the archive uses; values are numbers (possibly in 1p5 / 1e-3 style). We keep it liberal.
_KNOWN_TOKENS = ("curvature", "msc", "force", "length", "cc", "cs", "ncoils", "n")
_SYMMETRY_RE = re.compile(r"\b(QA|QH|QI)\b", re.IGNORECASE)
_NCOILS_RE = re.compile(r"(?:ncoils|ncpp|n)[_=]?(\d+)", re.IGNORECASE)


def parse_config_from_path(path):
    """Best-effort: pull symmetry (QA/QH/QI), ncoils, and any token=value pairs out of a path.
    Returns a dict; unrecognised chunks are preserved under 'raw_tokens' so nothing is silently lost."""
    cfg = {"symmetry": "", "ncoils_hint": "", "raw_tokens": ""}
    sym = _SYMMETRY_RE.search(path)
    if sym:
        cfg["symmetry"] = sym.group(1).upper()
    nc = _NCOILS_RE.search(path)
    if nc:
        cfg["ncoils_hint"] = nc.group(1)
    # token=value or token_value pairs anywhere in the path
    tokens = {}
    for m in re.finditer(r"([A-Za-z]+)[=_]([0-9pPeE.+\-]+)", path):
        key, val = m.group(1).lower(), m.group(2)
        if key in _KNOWN_TOKENS:
            tokens[key] = val.replace("p", ".")            # archive writes 1p5 for 1.5
    cfg.update({k: tokens.get(k, "") for k in _KNOWN_TOKENS})
    # anything left that looks like a config token but isn't known -> raw_tokens (audit trail)
    leftover = [f"{m.group(1)}={m.group(2)}" for m in re.finditer(r"([A-Za-z]+)[=_]([0-9pPeE.+\-]+)", path)
                if m.group(1).lower() not in _KNOWN_TOKENS]
    cfg["raw_tokens"] = ";".join(leftover)
    return cfg


def _nfp_for(symmetry, nfp_map, default=1):
    return nfp_map.get(symmetry.upper(), default)


def batch_convert(archive_root, out_dir, target_major_radius=None, nfp_map=None):
    nfp_map = nfp_map or {}
    os.makedirs(out_dir, exist_ok=True)
    jsons = sorted(glob.glob(os.path.join(archive_root, "**", "biot_savart_optimized_auglag_*.json"),
                             recursive=True))
    rows = []
    for jp in jsons:
        cfg = parse_config_from_path(os.path.relpath(jp, archive_root))
        nfp = _nfp_for(cfg["symmetry"], nfp_map)
        # unique, filesystem-safe output name derived from the relative path
        rel = os.path.relpath(jp, archive_root)
        stem = re.sub(r"[^A-Za-z0-9]+", "_", os.path.splitext(rel)[0]).strip("_")
        out = os.path.join(out_dir, f"coils_{stem}")
        try:
            bs = load(jp)
            r0 = coil_set_major_radius(bs.coils)
            n, scale = gil_to_makegrid(jp, out, nfp=nfp, scale=None,
                                       target_major_radius=target_major_radius)
            status = "ok"
        except Exception as e:                              # never let one bad set kill the batch
            n, scale, r0, status = 0, float("nan"), float("nan"), f"ERROR:{type(e).__name__}:{e}"
            out = ""
        rows.append(dict(source_json=rel, out_coils=os.path.basename(out) if out else "",
                         symmetry=cfg["symmetry"], ncoils=n, ncoils_hint=cfg["ncoils_hint"],
                         nfp=nfp, measured_R0_m=r0, applied_scale=scale,
                         cs=cfg["cs"], cc=cfg["cc"], length=cfg["length"], curvature=cfg["curvature"],
                         msc=cfg["msc"], force=cfg["force"], raw_tokens=cfg["raw_tokens"], status=status))
    manifest = os.path.join(out_dir, "manifest.csv")
    if rows:
        with open(manifest, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    return rows, manifest


def _parse_nfp_map(s):
    if not s:
        return {}
    out = {}
    for pair in s.split(","):
        k, v = pair.split("=")
        out[k.strip().upper()] = int(v)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("archive_root"); ap.add_argument("out_dir")
    ap.add_argument("--target-major-radius", type=float, default=None,
                    help="rescale every set to this R0 in METERS (else keep archive scale)")
    ap.add_argument("--nfp-map", default="QA=2,QH=4,QI=1",
                    help="symmetry->nfp, e.g. 'QA=2,QH=4' (default QA=2,QH=4,QI=1)")
    a = ap.parse_args()
    rows, manifest = batch_convert(a.archive_root, a.out_dir,
                                   target_major_radius=a.target_major_radius,
                                   nfp_map=_parse_nfp_map(a.nfp_map))
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"converted {ok}/{len(rows)} sets -> {a.out_dir}")
    for r in rows:
        if r["status"] != "ok":
            print("  FAILED:", r["source_json"], r["status"])
    print("manifest:", manifest)
