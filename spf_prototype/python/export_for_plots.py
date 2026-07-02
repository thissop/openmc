#!/usr/bin/env python
"""export_for_plots.py -- bundle everything needed to make plots + slides LOCALLY.

Run on Ginsburg in the spf-stellarator env. Produces a self-contained directory
`spf_export/` (and a tarball) containing ONLY plain data -- CSV, JSON, NPZ -- that
can be read on your Mac with numpy/pandas/json alone. No OpenMC, DAGMC, MOAB, or
cluster needed to read any of it.

Contents:
  results/tallies_long.csv     every tally value: stream, mode, score, bin, mean, std
  results/tally_meta.json      each tally's scores/filters/shape (to reshape/plot)
  results/RESULTS_tier8_conformal.md   copy of the headline results file
  results/summary.json         parsed headline metrics (TBR, coil flux, heating, eta)
  geometry/boundaries_RZ.npz   (theta,phi) R,Z grids for LCFS + every layer boundary
  geometry/cross_sections.csv  poloidal (R,Z) contours per layer at several phi (2D)
  geometry/lcfs_wireframe.csv  LCFS as (x,y,z) points for a 3D plasma wireframe
  geometry/shells_xyz.npz      subsampled (x,y,z) meshes of each shell (3D nested view)
  geometry/layers.json         layer names, offsets, thickness, material, color, volume
  meta.json                    scale, nfp, particle count, provenance, file index

Usage:
  python export_for_plots.py                 # auto-detect repo paths
  python export_for_plots.py --results DIR --geom DIR --stem NAME --out DIR
"""
import argparse
import glob
import itertools
import json
import shutil
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


# ----------------------------------------------------------------------------
def _log(msg):
    print(f"[export] {msg}", flush=True)


def _binval(b):
    """Readable scalar for a filter bin (cell/material id, energy tuple, mesh idx)."""
    a = np.atleast_1d(b).ravel()
    if a.size == 1:
        return a[0].item() if hasattr(a[0], "item") else a[0]
    return tuple(x.item() if hasattr(x, "item") else x for x in a)


CELL_TO_MAT = {1: "W", 2: "steel", 3: "Be", 4: "FLiBe", 5: "shield", 6: "coil"}
LAYER_COLORS = {"W": "#6e6e6e", "steel": "#9aa7b0", "Be": "#7fbf7f",
                "FLiBe": "#4f9dde", "shield": "#c9a24b", "coil": "#b5651d",
                "sol": "#dddddd", "LCFS": "#e05a8a"}


# ============================================================================
# 1. TALLIES  ->  long CSV + metadata
# ============================================================================
def export_tallies(results_dir, out):
    try:
        import openmc
    except Exception as e:  # noqa
        _log(f"SKIP tallies: openmc import failed ({e}). Activate spf-stellarator.")
        return {}
    sps = sorted(glob.glob(str(results_dir / "**/statepoint*.h5"), recursive=True))
    if not sps:
        _log(f"SKIP tallies: no statepoints under {results_dir}")
        return {}
    _log(f"found {len(sps)} statepoints")
    rows = []
    meta = {}
    for sp in sps:
        tag = Path(sp).parent.name                     # e.g. free_parallel
        stream, _, mode = tag.partition("_")
        try:
            with openmc.StatePoint(sp) as s:
                nb = s.n_batches
                for tid, t in s.tallies.items():
                    scores = [str(x) for x in t.scores]
                    nuclides = [str(x) for x in t.nuclides]
                    filters = [(type(f).__name__, [_binval(b) for b in f.bins])
                               for f in t.filters]
                    nfb = int(np.prod([len(fb) for _, fb in filters])) if filters else 1
                    m = t.mean.reshape(nfb, len(nuclides), len(scores))
                    sd = t.std_dev.reshape(nfb, len(nuclides), len(scores))
                    combos = (list(itertools.product(*[fb for _, fb in filters]))
                              if filters else [()])
                    for fi in range(nfb):
                        binstr = ";".join(str(x) for x in combos[fi])
                        # material name if the single filter bin is a known cell id
                        matname = ""
                        if len(combos[fi]) == 1 and combos[fi][0] in CELL_TO_MAT:
                            matname = CELL_TO_MAT[combos[fi][0]]
                        for ni in range(len(nuclides)):
                            for si in range(len(scores)):
                                mv = float(m[fi, ni, si])
                                ev = float(sd[fi, ni, si])
                                rows.append((stream, mode, t.name, int(tid), scores[si],
                                             nuclides[ni], binstr, matname, mv, ev,
                                             (ev / mv if mv else float("nan")), nb))
                    meta[t.name] = dict(
                        scores=scores, nuclides=nuclides,
                        filters=[(fn, [str(x) for x in fb]) for fn, fb in filters],
                        mean_shape=list(t.mean.shape))
        except Exception as e:  # noqa
            _log(f"  WARN failed to read {sp}: {e}")
    if not rows:
        return {}
    hdr = ("stream,mode,tally,tally_id,score,nuclide,filter_bins,material,"
           "mean,std_dev,rel_err,n_batches\n")
    lines = [hdr]
    for r in rows:
        lines.append(",".join(str(x) for x in r) + "\n")
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "results/tallies_long.csv").write_text("".join(lines))
    (out / "results/tally_meta.json").write_text(json.dumps(meta, indent=2))
    _log(f"wrote results/tallies_long.csv ({len(rows)} rows), tally_meta.json")
    return meta


# ============================================================================
# 2. HEADLINE SUMMARY  (parse straight from tallies, no pandas)
# ============================================================================
def export_summary(results_dir, out):
    try:
        import openmc
    except Exception:  # noqa
        return
    def score(sp, tally, cid, subs):
        with openmc.StatePoint(sp) as s:
            try:
                t = s.get_tally(name=tally)
            except Exception:  # noqa
                return float("nan"), float("nan")
            scs = list(t.scores)
            si = next((i for i, sc in enumerate(scs)
                       if any(ss in str(sc) for ss in subs)), None)
            if si is None:
                return float("nan"), float("nan")
            flat, mult, found = 0, 1, False
            for f in reversed(t.filters):
                bins = list(f.bins); idx = 0
                for j, b in enumerate(bins):
                    try:
                        if int(np.atleast_1d(b).ravel()[0]) == int(cid):
                            idx, found = j, True; break
                    except (TypeError, ValueError):
                        pass
                flat += idx * mult; mult *= len(bins)
            if t.filters and not found:
                return float("nan"), float("nan")
            m = t.mean.reshape(mult, -1, len(scs))
            sd = t.std_dev.reshape(mult, -1, len(scs))
            return float(m[flat, :, si].sum()), float(np.sqrt((sd[flat, :, si]**2).sum()))

    modes = ["unpolarized", "perpendicular", "parallel"]
    streams = ["free", "scatter"]
    S = {"tbr": {}, "coil_fast": {}, "heat": {}}
    for st in streams:
        for md in modes:
            sp = results_dir / f"{st}_{md}" / "statepoint.2.h5"
            if not sp.exists():
                cand = sorted(glob.glob(str(results_dir / f"{st}_{md}" / "statepoint*.h5")))
                if not cand:
                    continue
                sp = Path(cand[-1])
            S["tbr"][f"{st}_{md}"] = score(sp, "cell_response", 4, ("(n,Xt)", "H3", "(n,t)"))
            S["coil_fast"][f"{st}_{md}"] = score(sp, "coil_fast", 6, ("flux",))
            S["heat"][f"{st}_{md}"] = {n: score(sp, "cell_response", cid, ("heating",))[0]
                                       for cid, n in CELL_TO_MAT.items()}
    # directional efficiency eta on coil fast flux
    try:
        f_par = S["coil_fast"]["free_parallel"][0]; f_unp = S["coil_fast"]["free_unpolarized"][0]
        s_par = S["coil_fast"]["scatter_parallel"][0]; s_unp = S["coil_fast"]["scatter_unpolarized"][0]
        d_free = (f_par - f_unp) / f_unp
        d_scat = (s_par - s_unp) / s_unp
        S["eta"] = dict(delta_free=d_free, delta_scatter=d_scat,
                        eta=(d_scat / d_free if d_free else float("nan")))
    except Exception:  # noqa
        S["eta"] = {}
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "results/summary.json").write_text(json.dumps(S, indent=2))
    _log("wrote results/summary.json (TBR, coil flux, heating, eta)")


# ============================================================================
# 3. GEOMETRY  ->  boundaries, cross-sections, wireframe, layers
# ============================================================================
def export_geometry(repo, geom_dir, stem, out):
    sys.path.insert(0, str(repo / "spf_prototype/python"))
    try:
        import stellarator_geometry as sg
    except Exception as e:  # noqa
        _log(f"SKIP geometry: stellarator_geometry import failed ({e})")
        return None
    try:
        R, Z, nfp = sg.load_surface(stem)
    except Exception as e:  # noqa
        _log(f"SKIP geometry: load_surface failed ({e})")
        return None
    nt, nph = R.shape
    theta = 2 * np.pi * np.arange(nt) / nt
    phi = 2 * np.pi * np.arange(nph) / nph
    nR, nZ = sg.poloidal_outward_normals(R, Z)

    man = {}
    mp = geom_dir / "manifest.json"
    if mp.exists():
        man = json.loads(mp.read_text())
    layers = man.get("layers", [])
    # boundaries: LCFS (0) + each layer's OUTER offset (cumulative), in radial order
    bnames = ["LCFS"]
    boffs = [0.0]
    for L in layers:
        if L["name"] == "sol":
            # sol outer is the first material (W) inner boundary -> keep it as a boundary
            bnames.append("plasma_edge")
            boffs.append(L["outer_offset_cm"])
        else:
            bnames.append(L["name"] + "_outer")
            boffs.append(L["outer_offset_cm"])

    (out / "geometry").mkdir(parents=True, exist_ok=True)

    # 3a. R,Z grids per boundary -> NPZ
    npz = dict(theta=theta, phi=phi, nfp=np.array(nfp))
    RZ = {}
    for name, delta in zip(bnames, boffs):
        Ro, Zo = sg.offset_surface(R, Z, nR, nZ, float(delta))
        RZ[name] = (Ro, Zo)
        npz[f"{name}_R"] = Ro
        npz[f"{name}_Z"] = Zo
    np.savez(out / "geometry/boundaries_RZ.npz", **npz)
    _log(f"wrote geometry/boundaries_RZ.npz ({len(bnames)} boundaries, grid {nt}x{nph})")

    # 3b. poloidal cross-sections at several phi -> long CSV (easy 2D plot)
    phi_pick = sorted(set([0, nph // 8, nph // 4, 3 * nph // 8]))
    lines = ["phi_deg,phi_idx,layer,theta_idx,R_cm,Z_cm\n"]
    for pidx in phi_pick:
        deg = np.degrees(phi[pidx])
        for name in bnames:
            Ro, Zo = RZ[name]
            for ti in range(nt):
                lines.append(f"{deg:.3f},{pidx},{name},{ti},{Ro[ti, pidx]:.4f},{Zo[ti, pidx]:.4f}\n")
    (out / "geometry/cross_sections.csv").write_text("".join(lines))
    _log(f"wrote geometry/cross_sections.csv (phi indices {phi_pick})")

    # 3c. LCFS 3D wireframe -> CSV (x,y,z) and closed in both directions
    Rl, Zl = RZ["LCFS"]
    X = Rl * np.cos(phi)[None, :]
    Y = Rl * np.sin(phi)[None, :]
    wf = ["theta_idx,phi_idx,x_cm,y_cm,z_cm\n"]
    for ti in range(nt):
        for pj in range(nph):
            wf.append(f"{ti},{pj},{X[ti, pj]:.4f},{Y[ti, pj]:.4f},{Zl[ti, pj]:.4f}\n")
    (out / "geometry/lcfs_wireframe.csv").write_text("".join(wf))
    _log("wrote geometry/lcfs_wireframe.csv (3D plasma wireframe)")

    # 3d. subsampled 3D shell meshes -> NPZ (nested translucent shells figure)
    sub_t = max(1, nt // 40)
    sub_p = max(1, nph // 80)
    shells = {}
    for name in bnames:
        Ro, Zo = RZ[name]
        Rs, Zs = Ro[::sub_t, ::sub_p], Zo[::sub_t, ::sub_p]
        ph = phi[::sub_p]
        Xs = Rs * np.cos(ph)[None, :]
        Ys = Rs * np.sin(ph)[None, :]
        shells[f"{name}_x"] = Xs
        shells[f"{name}_y"] = Ys
        shells[f"{name}_z"] = Zs
    np.savez(out / "geometry/shells_xyz.npz", **shells)
    _log("wrote geometry/shells_xyz.npz (subsampled nested shells)")

    # 3e. layers.json (materials, offsets, thickness, color)
    lj = []
    for L in layers:
        lj.append(dict(name=L["name"], material=L["name"],
                       inner_offset_cm=L.get("inner_offset_cm"),
                       outer_offset_cm=L.get("outer_offset_cm"),
                       thickness_cm=(L.get("outer_offset_cm", 0) - L.get("inner_offset_cm", 0)),
                       enclosed_volume_cm3=L.get("enclosed_volume_cm3"),
                       color=LAYER_COLORS.get(L["name"], "#888888"),
                       is_vacuum=(L["name"] == "sol")))
    (out / "geometry/layers.json").write_text(json.dumps(
        dict(stem=stem, nfp=int(nfp), scale=man.get("scale"),
             boundary_order=bnames, layers=lj), indent=2))
    _log("wrote geometry/layers.json")
    return dict(nfp=int(nfp), scale=man.get("scale"), nt=nt, nph=nph)


# ============================================================================
def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve()
    repo_guess = here.parents[2] if len(here.parents) >= 3 else Path.cwd()
    ap.add_argument("--repo", default=str(repo_guess))
    ap.add_argument("--results", default=None)
    ap.add_argument("--geom", default=None)
    ap.add_argument("--stem", default="equil_precise_qa")
    ap.add_argument("--out", default="spf_export")
    a = ap.parse_args()

    repo = Path(a.repo)
    results_dir = Path(a.results) if a.results else repo / "results_conformal"
    geom_dir = Path(a.geom) if a.geom else repo / f"spf_prototype/data/{a.stem}_geom"
    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    _log(f"repo    = {repo}")
    _log(f"results = {results_dir}")
    _log(f"geom    = {geom_dir}")

    meta = export_tallies(results_dir, out)
    export_summary(results_dir, out)
    geo = export_geometry(repo, geom_dir, stem := a.stem, out)

    # copy the human-readable results file if present
    rm = results_dir / "RESULTS_tier8_conformal.md"
    if rm.exists():
        (out / "results").mkdir(parents=True, exist_ok=True)
        shutil.copy(rm, out / "results/RESULTS_tier8_conformal.md")
        _log("copied results/RESULTS_tier8_conformal.md")

    # top-level meta + file index
    index = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
    (out / "meta.json").write_text(json.dumps(dict(
        stem=stem, exported_utc=datetime.now(timezone.utc).isoformat(),
        geometry=geo, tally_names=list(meta.keys()),
        note="Plain-data export for local plotting. Read with numpy/pandas/json; "
             "no OpenMC/DAGMC needed.", files=index), indent=2))
    _log("wrote meta.json")

    # tarball for easy scp
    tarpath = Path(str(out) + ".tar.gz")
    with tarfile.open(tarpath, "w:gz") as tf:
        tf.add(out, arcname=out.name)
    _log(f"DONE -> {out}/  and  {tarpath}")
    print("\nFiles in bundle:")
    for f in index:
        print("  ", f)
    print(f"\nscp it home with:\n  scp {tarpath} <you>@<mac>:  "
          "# or pull from your Mac side")


if __name__ == "__main__":
    main()