#!/usr/bin/env python
"""Aggregate per-config sweep records and TEST the two hypotheses:
  (1) universality : eta_source collapses onto ONE curve eta_source(C) across the
      symmetry classes QA/QH/QI (else the direction-tensor structure, not scalar C,
      is the predictor).
  (2) separability : eta ~ eta_source(C) * A(tau) -- eta_source(C) invariant across
      blanket recipes while only A moves along an A(tau) curve.

THE eta DEFINITION (single source of truth; mirrored in docs/EXPERIMENTAL_DESIGN.md):
  For observable X (coil fast flux), the raw directional change of a config is
      delta_stream = (X_parallel - X_unpolarized) / X_unpolarized
  eta_source = delta_free(config)    / delta_free(anchor)      # pre-blanket, PRIMARY
  eta_coil   = delta_scatter(config) / delta_scatter(anchor)   # full blanket, VR-pending
  A          = eta_coil / eta_source                           # blanket attenuation
  where the anchor is the high-coherence (C ~ 1) config, so eta_*(anchor) = 1.

Outputs: sweep.csv, a fit report (stdout + report.json), and figures in --figs.
Runs on the Mac against the pulled-back per-config JSON (no OpenMC needed).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------- #
# Load + assemble
# --------------------------------------------------------------------------- #
def load_records(rec_dir):
    recs = []
    for p in sorted(Path(rec_dir).glob("config_*.json")):
        recs.append(json.loads(p.read_text()))
    return recs


def pick_anchor(recs, anchor_id=None):
    """The C~1 normalization anchor: an explicitly flagged config, else the highest-C
    record that has free-stream deltas."""
    have = [r for r in recs if r.get("delta_free_parallel") is not None]
    if anchor_id is not None:
        for r in have:
            if r["id"] == anchor_id:
                return r
    if not have:
        return None
    return max(have, key=lambda r: r.get("C", 0.0))


def assemble(recs, anchor):
    """One row per (config, blanket). eta_source/eta_coil/A normalized to anchor."""
    d_free_a = anchor.get("delta_free_parallel") if anchor else None
    d_scat_a = anchor.get("delta_scatter_parallel") if anchor else None
    rows = []
    for r in recs:
        row = dict(id=r["id"], C=r.get("C"), S_phi=r.get("S_phi"),
                   lambda_phi=r.get("lambda_phi"), reversal_frac=r.get("reversal_frac"),
                   nfp=r.get("nfp"), iota=r.get("iota"),
                   aspect=r.get("aspect"), symmetry_class=r.get("symmetry_class"),
                   blanket=r.get("blanket"), tau=r.get("tau"),
                   anisotropy=r.get("anisotropy"), angular_std_deg=r.get("angular_std_deg"),
                   status=r.get("status"), lost=r.get("lost_particles"),
                   audit_passed=(r.get("field_audit") or {}).get("passed"))
        df = r.get("delta_free_parallel"); ds = r.get("delta_scatter_parallel")
        row["delta_free"] = df
        row["delta_scatter"] = ds
        row["eta_source"] = (df / d_free_a) if (df is not None and d_free_a) else None
        row["eta_coil"] = (ds / d_scat_a) if (ds is not None and d_scat_a) else None
        row["A"] = (row["eta_coil"] / row["eta_source"]) \
            if (row["eta_coil"] is not None and row["eta_source"]) else None
        rows.append(row)
    return rows


def write_csv(rows, path):
    cols = ["id", "symmetry_class", "nfp", "iota", "aspect", "C", "S_phi",
            "lambda_phi", "reversal_frac", "anisotropy",
            "angular_std_deg", "blanket", "tau", "delta_free", "delta_scatter",
            "eta_source", "eta_coil", "A", "status", "lost", "audit_passed"]
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join("" if r.get(c) is None else str(r.get(c)) for c in cols))
    Path(path).write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# Fits + hypothesis tests
# --------------------------------------------------------------------------- #
def _finite(xs, ys):
    x = np.asarray(xs, float); y = np.asarray(ys, float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def fit_eta_source_of_C(rows):
    """Fit eta_source(C). Test whether the CHEAP C predicts the EXPENSIVE eta_source.
    Models: linear, and a 1-parameter power eta = C^p (anchored so C=1 -> eta=1).
    Report R^2 and residual scatter."""
    pts = [(r["C"], r["eta_source"]) for r in rows
           if r["C"] is not None and r["eta_source"] is not None
           and r["blanket"] in (None, "baseline")]
    C, eta = _finite(*zip(*pts)) if pts else (np.array([]), np.array([]))
    out = dict(n=len(C))
    if len(C) < 3:
        out["note"] = "insufficient transported configs for a fit"
        return out, (C, eta)
    # linear
    A = np.vstack([C, np.ones_like(C)]).T
    (a, b), *_ = np.linalg.lstsq(A, eta, rcond=None)
    pred = A @ [a, b]
    ss = 1 - np.sum((eta - pred) ** 2) / np.sum((eta - eta.mean()) ** 2)
    out["linear"] = dict(slope=float(a), intercept=float(b), R2=float(ss),
                         rms_resid=float(np.sqrt(np.mean((eta - pred) ** 2))))
    # power eta = C^p (fit p by least squares in log-log, C>0, eta>0)
    ok = (C > 0) & (eta > 0)
    if ok.sum() >= 3:
        p = float(np.sum(np.log(C[ok]) * np.log(eta[ok])) / np.sum(np.log(C[ok]) ** 2))
        predp = C ** p
        ssp = 1 - np.sum((eta - predp) ** 2) / np.sum((eta - eta.mean()) ** 2)
        out["power"] = dict(exponent=p, R2=float(ssp))
    return out, (C, eta)


def fit_eta_source_of_S(rows):
    """Fit eta_source(S_phi) -- the PARITY-CORRECT predictor (nematic order about the
    toroidal axis). The kernel is even in b_hat, so eta is a functional of the SECOND
    moment; S_phi, not the first-moment C, is the causal variable (docs/THEORY.md).
    Reports the linear fit AND the PARAMETER-FREE prediction eta_source ~ S_phi/S_phi
    (anchor): leading-order steering ~ S_phi, so this residual is a zero-free-parameter
    test of the physics."""
    tri = [(r["S_phi"], r["eta_source"], r["C"]) for r in rows
           if r.get("S_phi") is not None and r.get("eta_source") is not None
           and r.get("C") is not None and r["blanket"] in (None, "baseline")]
    out = dict(n=len(tri))
    if len(tri) < 3:
        out["note"] = "insufficient transported configs for a fit"
        return out, (np.array([]), np.array([]))
    S = np.array([t[0] for t in tri]); eta = np.array([t[1] for t in tri])
    Cs = np.array([t[2] for t in tri])
    A = np.vstack([S, np.ones_like(S)]).T
    (a, b), *_ = np.linalg.lstsq(A, eta, rcond=None)
    pred = A @ [a, b]
    ss = 1 - np.sum((eta - pred) ** 2) / np.sum((eta - eta.mean()) ** 2)
    out["linear"] = dict(slope=float(a), intercept=float(b), R2=float(ss),
                         rms_resid=float(np.sqrt(np.mean((eta - pred) ** 2))))
    S_anchor = float(S[np.argmax(Cs)])                 # anchor = highest-C config
    if S_anchor:
        predpf = S / S_anchor
        out["parameter_free"] = dict(
            model="eta_source = S_phi / S_phi(anchor)", S_phi_anchor=S_anchor,
            rms_resid=float(np.sqrt(np.mean((eta - predpf) ** 2))))
    return out, (S, eta)


def compare_predictors(rows):
    """Head-to-head: does the first-moment C or the second-moment S_phi predict
    eta_source better (tighter collapse = higher R^2)? The parity argument predicts
    S_phi wins, or ties C on the cap-shaped family where lambda_phi ~ C^2."""
    fitC, _ = fit_eta_source_of_C(rows)
    fitS, _ = fit_eta_source_of_S(rows)
    r2C = fitC.get("linear", {}).get("R2")
    r2Cp = fitC.get("power", {}).get("R2")
    r2S = fitS.get("linear", {}).get("R2")
    cand = [("C_linear", r2C), ("C_power", r2Cp), ("S_phi_linear", r2S)]
    best = max(cand, key=lambda kv: (kv[1] if kv[1] is not None else -np.inf))
    return dict(R2_C_linear=r2C, R2_C_power=r2Cp, R2_S_phi_linear=r2S,
                best_predictor=best[0], parameter_free_S_phi=fitS.get("parameter_free"),
                note=("S_phi is the parity-correct predictor; C ties it only on a "
                      "co-toroidal cap (lambda_phi ~ C^2). A C-vs-S_phi R^2 gap, or any "
                      "reversal_frac>0 device, favors S_phi -- see docs/THEORY.md."))


def multivariate(rows):
    """Does C alone predict eta_source, or do nfp/iota add information? Standardized
    linear regression eta_source ~ C + nfp + iota; compare R^2 to C-only."""
    data = [(r["eta_source"], r["C"], r["nfp"], r["iota"]) for r in rows
            if None not in (r.get("eta_source"), r.get("C"), r.get("nfp"), r.get("iota"))
            and r["blanket"] in (None, "baseline")]
    if len(data) < 5:
        return dict(note="insufficient data")
    Y = np.array([d[0] for d in data], float)
    X = np.array([[d[1], d[2], d[3]] for d in data], float)
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-12)

    def r2(cols):
        M = np.hstack([Xs[:, cols], np.ones((len(Y), 1))])
        beta, *_ = np.linalg.lstsq(M, Y, rcond=None)
        pred = M @ beta
        return 1 - np.sum((Y - pred) ** 2) / np.sum((Y - Y.mean()) ** 2), beta
    r2_C, _ = r2([0])
    r2_all, beta = r2([0, 1, 2])
    return dict(n=len(Y), R2_C_only=float(r2_C), R2_C_nfp_iota=float(r2_all),
                delta_R2=float(r2_all - r2_C),
                std_coeffs=dict(C=float(beta[0]), nfp=float(beta[1]), iota=float(beta[2])),
                verdict=("C suffices" if (r2_all - r2_C) < 0.05 else
                         "nfp/iota add signal -> scalar C incomplete"))


def _class_residuals(base, xkey):
    """Per-class mean residual about the pooled linear eta_source(xkey) fit."""
    pts = [(r[xkey], r["eta_source"], r.get("symmetry_class") or "unknown")
           for r in base if r.get(xkey) is not None]
    if len(pts) < 4:
        return None
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    A = np.vstack([x, np.ones_like(x)]).T
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    by = {}
    for (_, _, cls), e in zip(pts, resid):
        by.setdefault(cls, []).append(e)
    per = {k: dict(n=len(v), mean_resid=float(np.mean(v)), std_resid=float(np.std(v)))
           for k, v in by.items()}
    spread = (max(c["mean_resid"] for c in per.values())
              - min(c["mean_resid"] for c in per.values()))
    return dict(per_class=per, class_mean_spread=float(spread))


def universality(rows):
    """Do QA/QH/QI collapse onto ONE curve? Compare per-class residuals about the
    pooled linear fit for BOTH predictors (C and S_phi). A class with a coherent
    nonzero mean residual SEPARATES under that predictor; the predictor with the
    smaller class spread is the better one (and if neither collapses, the full tensor
    structure -- not any scalar -- is the predictor). NOTE: on this family class is
    aliased with C (QA high-C, QH low-C, disjoint), so a 'split' can be curvature in
    eta(x); de-alias (G2) before over-reading a separation."""
    base = [r for r in rows if r.get("eta_source") is not None
            and r["blanket"] in (None, "baseline")]
    if len(base) < 4:
        return dict(note="insufficient data")
    rC = _class_residuals(base, "C"); rS = _class_residuals(base, "S_phi")
    out = dict(by_C=rC, by_S_phi=rS)
    # verdict from whichever predictors are available (older records lack S_phi); the
    # tighter class collapse wins.
    avail = {k: v["class_mean_spread"] for k, v in (("C", rC), ("S_phi", rS)) if v}
    if avail:
        better = min(avail, key=avail.get)
        out["verdict"] = (f"classes collapse under {better} (universal)"
                          if avail[better] < 0.1 else
                          f"classes SEPARATE (min spread {avail[better]:.3f} under "
                          f"{better}) -> tensor structure, not a scalar, is the predictor")
        if rC and rS:
            out["class_alias_warning"] = ("class is aliased with C on this family; "
                                          "confirm with de-aliased configs (G2) before over-reading")
    return out


def separability(rows):
    """A(tau) from the blanket subset + the three failure-mode diagnostics.
    A config that appears at >1 blanket lets us test eta_source(C) invariance."""
    have_A = [r for r in rows if r.get("A") is not None and r.get("tau") is not None]
    out = dict(n_A=len(have_A))
    if len(have_A) >= 3:
        tau, A = _finite([r["tau"] for r in have_A], [r["A"] for r in have_A])
        # A(tau): expect decreasing; fit A = exp(-k tau) (log-linear)
        ok = A > 0
        if ok.sum() >= 3:
            k = -float(np.polyfit(tau[ok], np.log(A[ok]), 1)[0])
            out["A_of_tau"] = dict(model="A=exp(-k*tau)", k=k)
    # failure mode 3: eta_source(C) invariance across blankets (per config id)
    byid = {}
    for r in rows:
        if r.get("eta_source") is not None:
            byid.setdefault(r["id"], []).append(r)
    var = [np.std([x["eta_source"] for x in v]) for v in byid.values() if len(v) > 1]
    if var:
        out["eta_source_blanket_spread"] = float(np.mean(var))
        out["failure_near_source_contamination"] = bool(np.mean(var) > 0.1)
    # failure mode 2: A depends on C at fixed tau (geometric non-separability)
    if len(have_A) >= 5:
        C = np.array([r["C"] for r in have_A]); Av = np.array([r["A"] for r in have_A])
        tauv = np.array([r["tau"] for r in have_A])
        # partial correlation of A with C controlling tau (residualize both on tau)
        def resid_on(z, x):
            b = np.polyfit(x, z, 1); return z - np.polyval(b, x)
        rA = resid_on(Av, tauv); rC = resid_on(C, tauv)
        if np.std(rA) > 0 and np.std(rC) > 0:
            pcorr = float(np.corrcoef(rA, rC)[0, 1])
            out["A_vs_C_partial_corr_at_fixed_tau"] = pcorr
            out["failure_geometric_nonseparability"] = bool(abs(pcorr) > 0.5)
    out["note_energy_coupling"] = ("energy-coupling failure needs per-config spectra; "
                                   "not stored in these records (diagnostic pending)")
    return out


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(rows, fit, figs_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figs = Path(figs_dir); figs.mkdir(parents=True, exist_ok=True)
    CLASS_COLOR = {"QA": "#2a78d6", "QH": "#eb6834", "QI": "#1baf7a", "unknown": "#8a8884"}

    def cls_color(r):
        return CLASS_COLOR.get(r.get("symmetry_class"), "#8a8884")

    base = [r for r in rows if r.get("C") is not None and r["blanket"] in (None, "baseline")]

    # 1. central validation: cheap C predicts expensive eta_source
    tr = [r for r in base if r.get("eta_source") is not None]
    if tr:
        fig, ax = plt.subplots(figsize=(6.4, 5))
        for r in tr:
            ax.scatter(r["C"], r["eta_source"], color=cls_color(r), s=45, zorder=3,
                       edgecolor="white")
        if "linear" in fit:
            xs = np.linspace(min(r["C"] for r in tr), 1.0, 50)
            ax.plot(xs, fit["linear"]["slope"] * xs + fit["linear"]["intercept"],
                    color="#52514e", lw=1.5,
                    label=f"linear R^2={fit['linear']['R2']:.2f}")
        ax.set_xlabel("field-direction coherence  C  (cheap, coils only)")
        ax.set_ylabel("eta_source  (expensive, neutronics)")
        ax.set_title("Does C predict the SPF benefit without neutronics?")
        ax.legend(frameon=False)
        _classes_legend(ax, CLASS_COLOR, tr)
        fig.tight_layout(); fig.savefig(figs / "fig_C_predicts_eta.png", dpi=200)
        plt.close(fig)

    # 1b. the PARITY-CORRECT predictor: S_phi vs eta_source, with the parameter-free
    #     physics line eta = S_phi/S_phi(anchor) overlaid (no fit).
    trS = [r for r in base if r.get("eta_source") is not None and r.get("S_phi") is not None]
    if trS:
        fig, ax = plt.subplots(figsize=(6.4, 5))
        for r in trS:
            ax.scatter(r["S_phi"], r["eta_source"], color=cls_color(r), s=45, zorder=3,
                       edgecolor="white")
        Sa = max(trS, key=lambda r: r["C"])["S_phi"]          # anchor S_phi (highest C)
        if Sa:
            xs = np.linspace(min(r["S_phi"] for r in trS), max(r["S_phi"] for r in trS), 50)
            ax.plot(xs, xs / Sa, color="#52514e", lw=1.5, ls="--",
                    label="parameter-free  eta = S_phi / S_phi(anchor)")
        ax.set_xlabel("nematic order  S_phi = (3*lambda_phi - 1)/2  (parity-correct)")
        ax.set_ylabel("eta_source  (expensive, neutronics)")
        ax.set_title("Does the SECOND-moment S_phi predict the SPF benefit?")
        ax.legend(frameon=False)
        _classes_legend(ax, CLASS_COLOR, trS)
        fig.tight_layout(); fig.savefig(figs / "fig_S_predicts_eta.png", dpi=200)
        plt.close(fig)

    # 2. eta_source vs nfp
    tr_nfp = [r for r in base if r.get("eta_source") is not None and r.get("nfp")]
    if tr_nfp:
        fig, ax = plt.subplots(figsize=(6, 4.4))
        for r in tr_nfp:
            ax.scatter(r["nfp"], r["eta_source"], color=cls_color(r), s=40, edgecolor="white")
        ax.set_xlabel("field periods nfp"); ax.set_ylabel("eta_source")
        ax.set_title("eta_source vs nfp"); fig.tight_layout()
        fig.savefig(figs / "fig_eta_vs_nfp.png", dpi=200); plt.close(fig)

    # 3. A(tau) blanket-factorization curve
    hasA = [r for r in rows if r.get("A") is not None and r.get("tau") is not None]
    if hasA:
        fig, ax = plt.subplots(figsize=(6, 4.4))
        for r in hasA:
            ax.scatter(r["tau"], r["A"], color=cls_color(r), s=40, edgecolor="white")
        ax.set_xlabel("scattering optical-depth proxy  tau")
        ax.set_ylabel("A = eta_coil / eta_source")
        ax.set_title("Blanket attenuation factor A(tau)")
        fig.tight_layout(); fig.savefig(figs / "fig_A_of_tau.png", dpi=200); plt.close(fig)

    # 4. direction-spread visualization: C vs angular_std, sized by anisotropy
    ds = [r for r in base if r.get("angular_std_deg") is not None]
    if ds:
        fig, ax = plt.subplots(figsize=(6, 4.4))
        for r in ds:
            ax.scatter(r["C"], r["angular_std_deg"], color=cls_color(r),
                       s=30 + 120 * (r.get("anisotropy") or 0), edgecolor="white", alpha=0.8)
        ax.set_xlabel("C"); ax.set_ylabel("source-weighted angular std of b_hat [deg]")
        ax.set_title("Field-direction spread (size ~ tensor anisotropy)")
        fig.tight_layout(); fig.savefig(figs / "fig_direction_spread.png", dpi=200)
        plt.close(fig)
    return sorted(str(p.name) for p in figs.glob("*.png"))


def _classes_legend(ax, colors, rows):
    import matplotlib.patches as mp
    present = sorted(set(r.get("symmetry_class") or "unknown" for r in rows))
    ax.legend(handles=ax.get_legend_handles_labels()[0] +
              [mp.Patch(color=colors.get(c, "#8a8884"), label=c) for c in present],
              frameon=False, fontsize=9)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True, help="dir of config_*.json")
    ap.add_argument("--out", default="analysis", help="output dir (csv, report, figs)")
    ap.add_argument("--anchor", type=int, default=None, help="anchor config id (C~1)")
    ap.add_argument("--no-figs", action="store_true")
    args = ap.parse_args()

    recs = load_records(args.records)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    anchor = pick_anchor(recs, args.anchor)
    rows = assemble(recs, anchor)
    write_csv(rows, out / "sweep.csv")

    fit, _ = fit_eta_source_of_C(rows)
    fit_S, _ = fit_eta_source_of_S(rows)
    comparison = compare_predictors(rows)
    report = dict(
        n_configs=len(recs),
        n_transported=sum(1 for r in rows if r.get("eta_source") is not None),
        anchor_id=(anchor["id"] if anchor else None),
        eta_source_of_C=fit,
        eta_source_of_S_phi=fit_S,
        predictor_comparison=comparison,
        multivariate=multivariate(rows),
        universality=universality(rows),
        separability=separability(rows),
    )
    (out / "report.json").write_text(json.dumps(report, indent=2, default=float))
    figs = [] if args.no_figs else make_figures(rows, fit, out / "figs")

    print(f"[analyze] {report['n_configs']} configs, "
          f"{report['n_transported']} transported, anchor={report['anchor_id']}")
    print(f"[analyze] eta_source(C):     {fit.get('linear', fit.get('note'))}")
    print(f"[analyze] eta_source(S_phi): {fit_S.get('linear', fit_S.get('note'))}")
    print(f"[analyze] predictor winner: {comparison.get('best_predictor')} "
          f"(R2 C={comparison.get('R2_C_linear')}, S_phi={comparison.get('R2_S_phi_linear')})")
    print(f"[analyze] universality: {report['universality'].get('verdict', report['universality'].get('note'))}")
    print(f"[analyze] separability: A(tau)={report['separability'].get('A_of_tau')}")
    print(f"[analyze] wrote {out/'sweep.csv'}, {out/'report.json'}"
          + (f", figs: {figs}" if figs else ""))


if __name__ == "__main__":
    main()
