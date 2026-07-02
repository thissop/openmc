#!/usr/bin/env python
"""SINGLE OFFLINE ENTRYPOINT for a SLURM batch run on Ginsburg (no internet in the
job). It bootstraps everything offline, runs the conformal stellarator SPF
transport for unpolarized/perpendicular/parallel x free-streaming/scattering, and
writes a RESULTS markdown with the headline metrics. Driven by ginsburg_job.sh.

NETWORK: this script touches the network for NOTHING. All network steps (conda env,
pip/anarrima/stl_to_h5m install, cross-section download, generating field maps for
NEW equilibria via DESC) are done once on the LOGIN node -- see SETUP_GINSBURG.md.
For the committed precise_QA baseline, the field map + conformal geometry are
already in data/, so the job needs only the conda env + cross sections.

OFFLINE BOOTSTRAP (all local toolchain, no internet):
  1. build the compiled source .so if missing  (cmake + compiler; honors CONDA_PREFIX)
  2. build the conformal geometry STLs if missing  (numpy)
  3. build the DAGMC .h5m if missing  (stl_to_h5m)
  4. rescale the field map to the device size  (file I/O)

PERSISTENCE: statepoints + RESULTS go to --results-dir (default ./results_conformal),
NOT /tmp -- SLURM /tmp is node-local and is purged when the job ends.

Usage (inside the SLURM job, conda env active, OPENMC_CROSS_SECTIONS set):
  python run_ginsburg.py --stem equil_precise_qa --scale 10 \
      --particles 2000000 --batches 20 --results-dir $SLURM_SUBMIT_DIR/results_qa
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
DATADIR = REPO / "spf_prototype" / "data"

import build_and_run as br  # noqa: E402
import stellarator_geometry as sg  # noqa: E402
import run_conformal as rc  # noqa: E402

MODES = rc.MODES  # unpolarized / perpendicular / parallel -> (a,b,c)


# --------------------------------------------------------------------------
# Offline bootstrap
# --------------------------------------------------------------------------
def bootstrap(stem, scale, h5m_path):
    """Build .so, conformal geometry, .h5m, and the scaled field map -- all offline."""
    import json
    # .so: reuse the (login-node) prebuild if present. Building unconditionally would
    # (a) race when array tasks share one REPO/src/build, and (b) choke on a stale
    # CMakeCache from an rsync'd dev tree. Skip-if-exists avoids both.
    if br.SO.exists():
        print(f"[bootstrap] reusing existing .so {br.SO}", flush=True)
    else:
        print("[bootstrap] building compiled source .so ...", flush=True)
        br.build_so()  # honors OPENMC_PREFIX/CONDA_PREFIX

    # geometry: rebuild unless a cached manifest exists AT THE REQUESTED SCALE
    # (a cached scale-10 wall around a --scale!=10 source would silently misalign).
    geom_dir = DATADIR / f"{stem}_geom"
    mfp = geom_dir / "manifest.json"
    stale = True
    if mfp.exists():
        cached = abs(float(json.loads(mfp.read_text()).get("scale", -1)) - scale) < 1e-9
        stale = not cached
        if stale:
            print(f"[bootstrap] cached geometry scale != {scale}; rebuilding "
                  "(and the .h5m) ...", flush=True)
    if stale:
        m = sg.build_layers(stem, scale=scale, outdir=geom_dir)
        bad = [L["name"] for L in m["layers"]
               if not (L["watertight"] and L["simple_cross_section"])]
        if bad:
            raise SystemExit(f"geometry invalid (self-intersecting): {bad}; "
                             "increase --scale or thin the build (DEFERRED.md).")
        Path(h5m_path).unlink(missing_ok=True)  # force .h5m rebuild from new geometry

    if not Path(h5m_path).exists():
        print(f"[bootstrap] building DAGMC {h5m_path} ...", flush=True)
        import build_dagmc
        build_dagmc.build(geom_dir, h5m_path)

    fmap = f"{stem}_s{int(scale)}"
    if not (DATADIR / f"{fmap}.meta").exists():
        rc.scale_fieldmap(stem, scale, fmap)
    return fmap


# --------------------------------------------------------------------------
# Transport (persistent statepoints)
# --------------------------------------------------------------------------
def run_all(stem, scale, h5m_path, fmap, results_dir, particles, batches):
    import openmc
    R0_cm, a_cm = 103.0 * scale, 17.2 * scale
    sp = {}
    for stream, dscale in [("free", 1e-4), ("scatter", 1.0)]:
        for mode, abc in MODES.items():
            tag = f"{stream}_{mode}"
            cwd = results_dir / "statepoints" / tag
            cwd.mkdir(parents=True, exist_ok=True)
            model, mats = rc.build_model(abc, h5m_path, fmap, R0_cm, a_cm,
                                         particles=particles, batches=batches,
                                         density_scale=dscale)
            print(f"[run] {tag}: {particles}x{batches} ...", flush=True)
            sp[(stream, mode)] = (model.run(cwd=str(cwd), output=False), mats)
    return sp


# --------------------------------------------------------------------------
# Metrics (robust MaterialFilter extraction; statepoints already saved)
# --------------------------------------------------------------------------
def _mat_score(sp_path, tally_name, mat_id, score_substrs):
    """Mean+std of a MaterialFilter tally for one material id and a score whose
    label contains any of score_substrs (OpenMC renames H3-production -> (n,Xt))."""
    import openmc
    with openmc.StatePoint(sp_path) as s:
        df = s.get_tally(name=tally_name).get_pandas_dataframe()
    col_mat = "material" if "material" in df.columns else ("material", "")
    mid = df[col_mat].to_numpy()
    sc = df["score"].to_numpy() if "score" in df.columns else df[("score", "")].to_numpy()
    mean = df["mean"].to_numpy() if "mean" in df.columns else df[("mean", "")].to_numpy()
    sd = (df["std. dev."].to_numpy() if "std. dev." in df.columns
          else df[("std. dev.", "")].to_numpy())
    for i in range(len(mid)):
        if mid[i] == mat_id and any(ss in str(sc[i]) for ss in score_substrs):
            return float(mean[i]), float(sd[i])
    return float("nan"), float("nan")


def extract_metrics(sp, results_dir):
    """Headline metrics from the MaterialFilter tallies (robust path). The
    phi-resolved wall mesh + the per-patch analytic comparison are left as a
    documented postprocessing hook (needs the analytic stellarator NWL).

    Material ids are taken PER CONFIG from that config's own `mats` (the 6
    build_model calls assign their own ids), not cached from the first run."""
    M = {}
    for (stream, mode), (sp_path, mats) in sp.items():
        tbr, tbr_sd = _mat_score(sp_path, "cell_response", mats["FLiBe"].id,
                                 ("(n,Xt)", "H3", "(n,t)"))
        cf, cf_sd = _mat_score(sp_path, "coil_fast", mats["coil"].id, ("flux",))
        heat = {n: _mat_score(sp_path, "cell_response", mats[n].id, ("heating",))[0]
                for n in ("W", "steel", "Be", "FLiBe", "shield", "coil")}
        M[(stream, mode)] = dict(tbr=tbr, tbr_sd=tbr_sd, coil_fast=cf,
                                 coil_fast_sd=cf_sd, heat=heat)
    return M


def write_results(M, stem, scale, particles, batches, results_dir):
    def g(stream, mode, k):
        return M.get((stream, mode), {}).get(k, float("nan"))

    L = []; w = L.append
    w("# RESULTS — Tier 8 conformal stellarator SPF (Ginsburg)\n")
    w(f"Equilibrium `{stem}`, device scale {scale} (R0~{103*scale:.0f}cm, "
      f"a~{17.2*scale:.0f}cm), conformal DAGMC wall, {particles*batches:,} "
      "histories/config. Modes: unpolarized=iso, perpendicular=A(sin^2), "
      "parallel=B/C(1+3cos^2). Per source neutron, rate held fixed.\n")
    w("> Machinery validation on a public equilibrium (precise_QA), NOT a Helios "
      "physics claim. Coil tallies are deep behind the blanket+shield; without "
      "variance reduction expect large MC error -- read the +/- and prefer the "
      "trend (see DEFERRED.md).\n")

    w("## Global TBR per mode (FLiBe H3-production)\n")
    w("| mode | free-streaming | scattering |")
    w("|---|---|---|")
    for mode in MODES:
        w(f"| {mode} | {g('free',mode,'tbr'):.3f} ± {g('free',mode,'tbr_sd'):.3f} "
          f"| {g('scatter',mode,'tbr'):.3f} ± {g('scatter',mode,'tbr_sd'):.3f} |")

    w("\n## Coil fast flux (>0.1 MeV) per mode -- magnet-lifetime proxy\n")
    w("| mode | free-streaming [/cm²/src] | scattering [/cm²/src] |")
    w("|---|---|---|")
    for mode in MODES:
        fr, frs = g('free', mode, 'coil_fast'), g('free', mode, 'coil_fast_sd')
        sc, scs = g('scatter', mode, 'coil_fast'), g('scatter', mode, 'coil_fast_sd')
        w(f"| {mode} | {fr:.3e} (±{100*frs/fr if fr else float('nan'):.0f}%) "
          f"| {sc:.3e} (±{100*scs/sc if sc else float('nan'):.0f}%) |")

    # directional efficiency eta: parallel coil-flux reduction, scatter vs free
    def red(stream):
        u = g(stream, "unpolarized", "coil_fast"); p = g(stream, "parallel", "coil_fast")
        return (1 - p / u) if u else float("nan")
    rf, rs = red("free"), red("scatter")
    eta = (rs / rf) if rf else float("nan")
    w(f"\n**Directional efficiency η (coil fast flux):** parallel reduces coil flux "
      f"by {100*rf:.1f}% free-streaming → {100*rs:.1f}% with scattering ⇒ "
      f"**η = {eta:.2f}** (fraction of the steering benefit that survives "
      "non-axisymmetric conformal walls + scattering). The headline answer to "
      "'does SPF steering survive'. (Large coil MC error without variance reduction "
      "-- treat as indicative; rerun with weight windows for a precise η.)\n")

    w("## Per-material heating (scattering, per source neutron) [eV]\n")
    w("| layer | unpolarized | perpendicular | parallel |")
    w("|---|---|---|---|")
    for n in ("W", "steel", "Be", "FLiBe", "shield", "coil"):
        w(f"| {n} | " + " | ".join(
            f"{g('scatter',m,'heat').get(n, float('nan')):.3e}" for m in MODES) + " |")

    w("\n## Not in this table (postprocessing hooks)\n")
    w("- **φ-resolved first-wall load** (toroidal/poloidal map): the `wall_current_phi` "
      "mesh tally is saved in each statepoint; extract the poloidal map there.\n")
    w("- **Joint analytic↔MC free-streaming cross-check**: compare the free-streaming "
      "per-patch wall load to the analytic stellarator NWL (companion "
      "field-period-perturbation code / anarrima angled kernels) on the same "
      "precise_QA patches. See JOINT_VALIDATION.md.\n")
    (results_dir / "RESULTS_tier8_conformal.md").write_text("\n".join(L) + "\n")
    print(f"[results] wrote {results_dir / 'RESULTS_tier8_conformal.md'}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default="equil_precise_qa")
    ap.add_argument("--scale", type=float, default=10.0)
    ap.add_argument("--h5m", default=None, help="DAGMC file (default <stem>.h5m in results dir)")
    ap.add_argument("--particles", type=int, default=2_000_000)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--results-dir", default="results_conformal")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    h5m = args.h5m or str(results_dir / f"{args.stem}.h5m")

    print(f"[run_ginsburg] stem={args.stem} scale={args.scale} "
          f"particles={args.particles}x{args.batches} -> {results_dir}", flush=True)
    fmap = bootstrap(args.stem, args.scale, h5m)
    sp = run_all(args.stem, args.scale, h5m, fmap, results_dir, args.particles, args.batches)
    print("[run_ginsburg] transport done; statepoints saved.", flush=True)
    rfile = results_dir / "RESULTS_tier8_conformal.md"
    # statepoints are already on disk (run_all), so a metrics failure never loses
    # the expensive transport -- but it MUST exit nonzero so an sbatch-only user
    # sees FAILED in sacct, not a silently-garbage COMPLETED.
    try:
        M = extract_metrics(sp, results_dir)
    except Exception:
        print("[run_ginsburg] ERROR: metrics extraction raised; statepoints are "
              "preserved in results-dir/statepoints.", flush=True)
        traceback.print_exc()
        rfile.write_text("# RUN FAILED: metrics extraction raised; statepoints "
                         "preserved in statepoints/. See job stderr.\n")
        sys.exit(1)
    import math
    finite = all(math.isfinite(M.get((s, m), {}).get("tbr", float("nan")))
                 and math.isfinite(M.get((s, m), {}).get("coil_fast", float("nan")))
                 for s in ("free", "scatter") for m in MODES)
    if not finite:
        print("[run_ginsburg] ERROR: non-finite metrics (likely a material-id "
              "mismatch); statepoints preserved.", flush=True)
        rfile.write_text("# RUN FAILED: non-finite TBR/coil metrics (likely "
                         "material-id mismatch); statepoints preserved in "
                         "statepoints/.\n")
        sys.exit(1)
    write_results(M, args.stem, args.scale, args.particles, args.batches, results_dir)


if __name__ == "__main__":
    main()
