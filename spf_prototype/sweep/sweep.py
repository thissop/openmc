#!/usr/bin/env python
"""Per-config SPF sweep driver: QUASR coils -> field -> HARD-GATE audit -> coherence
C -> conformal geometry -> DAGMC -> transport -> two eta observables -> a per-config
JSON record. Idempotent and resumable (one JSON per config; a crash loses at most
one config). The SAME code runs the local smoke test and the Ginsburg production run,
parameterized by a config manifest + argparse.

Pipeline per config (reusing the debugged pipeline, not duplicating it):
  quasr_loader.load_device      coils + boundary + metadata (QUASR)
  biotsavart_field.CoilField    exact b_hat from the coils
  field_audit.audit_coilfield   |b|=1, div B~0, B.n~0 on LCFS   [HARD GATE: skip on fail]
  coherence_metrics             C, direction tensor, angular std   (the cheap predictor)
  biotsavart_field.write_fieldmap  -> spf_fieldmap_v1 the .so reads
  stellarator_geometry.build_layers -> conformal shell STLs   [self-intersection gate]
  dagmc_writer.build_from_stls  -> watertight DAGMC .h5m (pymoab-free)
  run_conformal.build_model + model.run  -> statepoints (Ginsburg only)
  run_ginsburg._mat_score       -> tallies (array API, no pandas)

Two eta observables (see docs/EXPERIMENTAL_DESIGN.md, the single eta definition):
  eta_source : FREE-stream (pre-blanket, near-void) coil-flux directional change.
               Fast-converging. PRIMARY predictor for eta_source(C).
  eta_coil   : SCATTER-stream (full blanket) coil-flux directional change. Deep,
               noisy, VR-pending. A = eta_coil / eta_source is the blanket factor.
The record stores RAW fractional changes (delta_free/scatter for each polarized
mode); analyze_sweep.py normalizes to the C~1 anchor to form eta.

Units: QUASR is metres; the geometry/source pipeline is cm. UNIT_CM converts. b_hat
is scale-invariant, so the audit and C run in native units; the field map and STLs
are written in cm.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))

import quasr_loader as ql
import coherence_metrics as cm
import field_audit as fa
import biotsavart_field as bsf
import stellarator_geometry as sg
import dagmc_writer

DATADIR = HERE.parent / "data"
UNIT_CM = 100.0                      # QUASR metres -> cm
TARGET_A_CM = 250.0                  # built minor radius target (precise_QA scale-15 ~ 258)
POL_MODES = ("perpendicular", "parallel")   # vs unpolarized

# blanket recipes (thickness knob; tau proxy from macroscopic scatter power * thickness).
# baseline == stellarator_geometry.DEFAULT_LAYERS. See docs/FACTORIZATION.md.
BLANKETS = {
    "baseline": [("sol", 5.0), ("W", 0.2), ("steel", 3.8), ("Be", 2.0),
                 ("FLiBe", 50.0), ("shield", 40.0), ("coil", 15.0)],
    "thin":     [("sol", 5.0), ("W", 0.2), ("steel", 2.0), ("Be", 1.0),
                 ("FLiBe", 30.0), ("shield", 25.0), ("coil", 15.0)],
    "thick":    [("sol", 5.0), ("W", 0.2), ("steel", 5.0), ("Be", 3.0),
                 ("FLiBe", 70.0), ("shield", 60.0), ("coil", 15.0)],
}
# crude macroscopic scattering power [1/cm] per material (fast-neutron, order-of-mag);
# tau ~ sum(Sigma_s * thickness). Proxy only (documented). FLiBe/steel/shield dominate.
SIGMA_S = {"sol": 0.0, "W": 0.10, "steel": 0.12, "Be": 0.30, "FLiBe": 0.13,
           "shield": 0.11, "coil": 0.10}


def tau_proxy(layers):
    """Scattering optical-depth proxy: sum of macroscopic-scatter-power * thickness
    over the material shells (vacuum 'sol' excluded)."""
    return float(sum(SIGMA_S.get(n, 0.0) * t for n, t in layers if n != "sol"))


def _make_field(dev, cfg, defaults, stem):
    """The b_hat source for a config: CoilField (Biot-Savart, DEFAULT) or, when
    field_source='equilibrium' (G3), EquilibriumField reading the pre-generated
    data/<stem>_equil map -- a flux function that admits the reversal devices the coil
    B.n audit filters. Both expose the same .B/.bhat interface downstream."""
    src = cfg.get("field_source", defaults.get("field_source", "coil"))
    if src == "equilibrium":
        import equilibrium_field as ef
        return ef.EquilibriumField(str(DATADIR / f"{stem}_equil"), meta=dev.meta)
    return bsf.CoilField(dev.coils, meta=dev.meta)


# --------------------------------------------------------------------------- #
# Geometry staging: write the surface npz + field map a config needs
# --------------------------------------------------------------------------- #
def stage_surface_and_field(dev, field_native, stem, nR=48, nphi=96, nZ=48,
                            ntheta=64, nph=96, margin=0.06):
    """Write data/<stem>_surface.npz (cm) for build_layers, and the biotsavart
    field map data/<stem>.meta/.bin (cm) for the source. Returns (a_native_cm)."""
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z = dev.device.RZ(TH, PH, 1.0)                 # native
    Rc, Zc = R * UNIT_CM, Z * UNIT_CM                  # cm
    np.savez(DATADIR / f"{stem}_surface.npz", R=Rc, Z=Zc, nfp=np.int64(dev.meta["nfp"]),
             theta=th, phi_axis=ph, phi=np.broadcast_to(ph[None, :], Rc.shape).copy())
    a_native = 0.5 * (R.max() - R.min())

    # field map in cm: bounds = LCFS bbox (cm) + margin; field from cm-scaled coils
    coils_cm = [(p * UNIT_CM, c) for p, c in dev.coils]
    field_cm = bsf.CoilField(coils_cm)
    dR = (Rc.max() - Rc.min()) * margin
    dZ = (Zc.max() - Zc.min()) * margin
    bounds = dict(R_min=float(Rc.min() - dR), R_max=float(Rc.max() + dR),
                  Z_min=float(Zc.min() - dZ), Z_max=float(Zc.max() + dZ))
    bsf.write_fieldmap(field_cm, bounds, dict(nR=nR, nphi=nphi, nZ=nZ),
                       str(DATADIR / stem), nfp=dev.meta["nfp"],
                       source=f"biotsavart:quasr{dev.ID}")
    return float(a_native)


def build_geometry(stem, layers, scale, geom_dir):
    """build_layers + hard self-intersection gate. Returns manifest or raises."""
    m = sg.build_layers(stem, layers=layers, scale=scale, outdir=geom_dir)
    bad = [L["name"] for L in m["layers"]
           if not (L["watertight"] and L["simple_cross_section"])]
    if bad:
        raise GeometryFold(bad)
    return m


class GeometryFold(Exception):
    def __init__(self, bad):
        super().__init__(f"conformal offset folds inboard: {bad}")
        self.bad = bad


# --------------------------------------------------------------------------- #
# Transport (Ginsburg only; guarded so the local smoke test skips it cleanly)
# --------------------------------------------------------------------------- #
def have_openmc():
    try:
        import openmc  # noqa
        return True
    except Exception:
        return False


def run_transport(stem, scale, h5m, layers, a_native_cm, results_dir, particles,
                  batches, streams, vr=False):
    """Run the requested streams and return raw tally deltas + lost-particle count.
    streams: 'free' (3 streams) or 'both' (6). Reuses run_conformal.build_model and
    run_ginsburg._mat_score."""
    import run_conformal as rc
    import run_ginsburg as rg
    import build_and_run as br

    if not br.SO.exists():
        raise SystemExit(
            f"compiled source .so missing at {br.SO}. Prebuild it ONCE on a compute "
            "node before the array job (array tasks must not race-build it): "
            "cd spf_prototype/src && cmake -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX . "
            "&& cmake --build build. See docs/README.md runbook.")

    fmap = f"{stem}_s{int(round(scale))}"
    rc.scale_fieldmap(stem, scale, fmap)
    R0_cm = float(a_native_cm) * UNIT_CM * 0.0  # placeholder; real R0 below
    # R0 from the surface npz mean major radius (cm) * scale
    surf = np.load(DATADIR / f"{stem}_surface.npz")
    R0_cm = float(surf["R"].mean()) * scale
    a_cm = float(a_native_cm) * UNIT_CM * scale

    stream_list = [("free", 1e-4)] if streams == "free" else [("free", 1e-4), ("scatter", 1.0)]
    vals = {}          # (stream, mode) -> dict(coil_fast=(m,sd), wall_dir=(m,sd))
    lost = 0
    for stream, dscale in stream_list:
        for mode, abc in rc.MODES.items():
            tag = f"{stream}_{mode}"
            cwd = Path(results_dir) / "statepoints" / tag
            cwd.mkdir(parents=True, exist_ok=True)
            model, mats = rc.build_model(abc, h5m, fmap, R0_cm, a_cm,
                                         particles=particles, batches=batches,
                                         density_scale=dscale)
            if vr and stream == "scatter":
                _attach_weight_windows(model, mats)
            sp = model.run(cwd=str(cwd), output=False)
            lost = max(lost, _lost_particles(sp))
            # BOTH observables per run (statepoints persist for re-analysis): the deep
            # coil_fast (eta_coil) and the near-source wall directional contrast (eta_source)
            vals[(stream, mode)] = dict(
                coil_fast=rg._mat_score(sp, "coil_fast", mats["coil"].id, ("flux",)),
                wall_dir=rg._wall_directional(sp, R0_cm))
    return vals, lost


def _lost_particles(sp_path):
    """Read the lost-particle count from the statepoint if available (0 otherwise)."""
    try:
        import openmc
        with openmc.StatePoint(sp_path) as s:
            return int(getattr(s, "n_lost_particles", 0) or 0)
    except Exception:
        return 0


def _attach_weight_windows(model, mats):
    """VR HOOK (MAGIC weight windows targeting the coil). Wired but conservative;
    tune before trusting a precise eta_coil. Enabled per config via --vr."""
    try:
        import openmc
        ww_mesh = openmc.RegularMesh()
        bb = model.geometry.bounding_box
        ww_mesh.lower_left = bb[0]
        ww_mesh.upper_right = bb[1]
        ww_mesh.dimension = (10, 10, 10)
        wwg = openmc.WeightWindowGenerator(
            method="magic", mesh=ww_mesh, max_realizations=model.settings.batches)
        model.settings.weight_window_generators = [wwg]
    except Exception as e:  # noqa
        print(f"[vr] weight-window hook unavailable: {e}", flush=True)


def deltas_from_vals(vals):
    """Directional steering deltas, polarized vs unpolarized (see docs/EXPERIMENTAL_DESIGN.md):
      delta_free_<mode>    := near-source WALL directional-contrast DIFFERENCE (free stream)
                              -> the PRIMARY, low-noise eta_source (parity-correct S_phi test).
      delta_scatter_<mode> := deep COIL_FAST FRACTIONAL change (scatter stream) -> eta_coil.
    Every raw (observable, stream, mode) value is also stored under 'observables_raw' so no
    output is lost -- results are pulled back to the Mac and re-analyzed there."""
    out = {}

    # PRIMARY eta_source: the near-source wall CONTRAST is already self-normalizing
    # (an in/out asymmetry in [-1,1]), so its steering signal is the DIFFERENCE from the
    # unpolarized load, not a fractional change.
    uw = vals.get(("free", "unpolarized"), {}).get("wall_dir")
    if uw and np.isfinite(uw[0]):
        for mode in POL_MODES:
            p = vals.get(("free", mode), {}).get("wall_dir")
            if p and np.isfinite(p[0]):
                out[f"delta_free_{mode}"] = float(p[0] - uw[0])
                out[f"delta_free_{mode}_sd"] = float(np.hypot(p[1], uw[1]))

    # eta_coil: deep coil_fast FRACTIONAL change (scatter stream)
    uc = vals.get(("scatter", "unpolarized"), {}).get("coil_fast")
    if uc and np.isfinite(uc[0]) and uc[0] != 0:
        for mode in POL_MODES:
            p = vals.get(("scatter", mode), {}).get("coil_fast")
            if p and np.isfinite(p[0]):
                out[f"delta_scatter_{mode}"] = float((p[0] - uc[0]) / uc[0])
                out[f"delta_scatter_{mode}_sd"] = float(
                    abs(p[0] / uc[0]) * np.hypot(p[1] / p[0] if p[0] else 0.0,
                                                 uc[1] / uc[0]))

    # save EVERY observable for EVERY (stream, mode) -- lossless for local re-analysis
    raw = {}
    for (stream, mode), d in vals.items():
        for obs, ms in d.items():
            if ms is not None:
                raw[f"{obs}__{stream}__{mode}"] = [float(ms[0]), float(ms[1])]
    out["observables_raw"] = raw
    return out


# --------------------------------------------------------------------------- #
# One config
# --------------------------------------------------------------------------- #
def run_config(cfg, defaults, out_dir, force=False):
    """Process one config dict; write out_dir/config_<id>.json. Returns the record."""
    cid = cfg["id"]
    blanket = cfg.get("blanket", defaults.get("blanket", "baseline"))
    stem = f"quasr{cid}"                       # field/surface are blanket-independent
    tag = f"{cid}_{blanket}"                   # record + geometry are blanket-specific
    rec_path = Path(out_dir) / f"config_{tag}.json"
    if rec_path.exists() and not force:
        print(f"[sweep] {tag}: record exists, skipping (use --force)", flush=True)
        return json.loads(rec_path.read_text())

    rec = dict(id=cid, stem=stem, blanket=blanket, status="started", warnings=[])
    try:
        # ---- field + audit (HARD GATE) --------------------------------------
        dev = ql.load_device(cid, n_samples=cfg.get("coil_samples", 256))
        rec.update(nfp=dev.meta["nfp"], iota=dev.meta["iota"],
                   aspect=dev.meta["aspect"], symmetry_class=cfg.get(
                       "symmetry_class", dev.meta["symmetry_class"]),
                   total_coil_length=dev.meta["total_coil_length"])
        field_native = _make_field(dev, cfg, defaults, stem)
        passed, audit = fa.audit_coilfield(field_native, dev)
        rec["field_audit"] = _audit_summary(audit)
        if not passed:
            rec["status"] = "audit_failed"
            _write(rec_path, rec)
            print(f"[sweep] {cid}: AUDIT FAILED -> skipped (not transported)", flush=True)
            return rec

        # ---- coherence metrics (the cheap predictor) ------------------------
        xyz, w, rho = dev.source_sample(model=cfg.get("source_model", "one_minus_rho2"))
        bhat = field_native.bhat(xyz)
        m = cm.coherence_metrics(bhat, pos=xyz, weights=w, frame="cylindrical")
        rec["C"] = m["C"]
        rec["tensor_evals"] = m["tensor_evals"].tolist()
        rec["anisotropy"] = m["anisotropy"]              # biaxiality (diagnostic only)
        rec["lambda_phi"] = m["lambda_phi"]             # <b_phi^2>_s (toroidal axis)
        rec["S_phi"] = m["S_phi"]                       # nematic order = HEADLINE predictor
        rec["reversal_frac"] = m["reversal_frac"]       # >0 => C under-predicts (use S_phi)
        rec["angular_std_deg"] = float(np.degrees(m["angular_std"]))
        rec["C_lab"] = cm.coherence_C(bhat, frame="lab")   # diagnostic (should be ~0)

        # ---- blanket recipe + tau -------------------------------------------
        layers = BLANKETS[blanket]
        rec["tau"] = tau_proxy(layers)

        # ---- stage geometry inputs + auto-scale -----------------------------
        a_native = stage_surface_and_field(
            dev, field_native, stem,
            nR=defaults.get("nR", 48), nphi=defaults.get("nphi", 96),
            nZ=defaults.get("nZ", 48))
        scale = cfg.get("scale") or (TARGET_A_CM / (a_native * UNIT_CM))
        geom_dir = DATADIR / f"{stem}_{blanket}_geom"    # blanket-specific build

        # ---- conformal geometry with self-intersection scale-retry ----------
        m_geo = None
        for attempt in range(defaults.get("scale_retries", 3)):
            try:
                m_geo = build_geometry(stem, layers, scale, geom_dir)
                break
            except GeometryFold as gf:
                rec["warnings"].append(f"geometry fold at scale={scale:.2f}: {gf.bad}")
                scale *= 1.4     # thinner build relative to minor radius
        if m_geo is None:
            rec["status"] = "geometry_folded"
            _write(rec_path, rec)
            print(f"[sweep] {cid}: geometry folds at all scales -> skipped", flush=True)
            return rec
        rec["scale"] = float(scale)
        rec["built_minor_cm"] = float(a_native * UNIT_CM * scale)

        # ---- DAGMC (pymoab-free) --------------------------------------------
        h5m = str(geom_dir / f"{stem}.h5m")
        info = dagmc_writer.build_from_stls(str(geom_dir), h5m)
        errs = dagmc_writer.verify_structure(h5m)
        if errs:
            rec["status"] = "dagmc_invalid"
            rec["warnings"] += errs
            _write(rec_path, rec)
            print(f"[sweep] {cid}: DAGMC structure errors -> skipped", flush=True)
            return rec
        rec["dagmc"] = dict(volumes=info["volumes"], tris=info["tris"])

        # ---- transport (Ginsburg only) --------------------------------------
        if not have_openmc():
            rec["status"] = "staged_no_transport"
            _write(rec_path, rec)
            print(f"[sweep] {cid}: staged (C={rec['C']:.3f}); transport needs "
                  "openmc+DAGMC (Ginsburg).", flush=True)
            return rec

        results_dir = Path(out_dir) / f"transport_{tag}"
        vals, lost = run_transport(
            stem, scale, h5m, layers, a_native,
            results_dir, cfg.get("particles", defaults.get("particles", 20000)),
            cfg.get("batches", defaults.get("batches", 5)),
            cfg.get("streams", defaults.get("streams", "free")),
            vr=cfg.get("vr", defaults.get("vr", False)))
        rec["lost_particles"] = lost
        if lost > defaults.get("lost_particle_max", 10):
            rec["warnings"].append(f"{lost} lost particles")
        rec.update(deltas_from_vals(vals))
        rec["status"] = "done"
        _write(rec_path, rec)
        print(f"[sweep] {cid}: DONE  C={rec['C']:.3f}  "
              f"delta_free_parallel={rec.get('delta_free_parallel')}", flush=True)
        return rec

    except Exception as e:  # noqa
        rec["status"] = "error"
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc()
        _write(rec_path, rec)
        print(f"[sweep] {cid}: ERROR {e}", flush=True)
        return rec


def _audit_summary(audit):
    return {c["name"]: dict(value=c["value"], passed=c["passed"],
                            **({"max": c["max"]} if "max" in c else {}))
            for c in audit["checks"]} | {"passed": audit["passed"]}


def _write(path, rec):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(rec, indent=2))


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="config manifest JSON")
    ap.add_argument("--out", default="sweep_out", help="per-config JSON output dir")
    ap.add_argument("--only", type=int, default=None, help="run a single config id")
    ap.add_argument("--index", type=int, default=None,
                    help="run the i-th config (for SLURM array $SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--force", action="store_true", help="recompute existing records")
    args = ap.parse_args()

    man = json.loads(Path(args.manifest).read_text())
    defaults = man.get("defaults", {})
    configs = man["configs"]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.only is not None:
        configs = [c for c in configs if c["id"] == args.only]
    elif args.index is not None:
        configs = [configs[args.index]]

    for cfg in configs:
        run_config(cfg, defaults, out_dir, force=args.force)


if __name__ == "__main__":
    main()
