"""Validate analyze_sweep on SYNTHETIC records with a KNOWN injected law
(no transport needed): eta_source(C) = C, A(tau) = exp(-k tau), class-independent.
The analyzer must recover the linear law (slope~1, intercept~0), declare the
classes universal, recover k, and find eta_source invariant across blankets."""
import json
from pathlib import Path

import numpy as np
import pytest

import analyze_sweep as az


def _make_records(tmp):
    rng = np.random.default_rng(1)
    d = Path(tmp); d.mkdir(parents=True, exist_ok=True)
    k1 = 0.135                      # anchor free-stream delta at C=1
    ktau = 0.05                     # A = exp(-ktau * tau)
    classes = ["QA", "QH", "QI"]
    # main C-scan (baseline blanket): eta_source = C exactly (+ tiny noise)
    Cs = np.linspace(0.35, 1.0, 12)
    for i, C in enumerate(Cs):
        df = k1 * C                 # delta_free ~ C  -> eta_source = C
        tau = 13.0
        ds = df * np.exp(-ktau * tau)
        rec = dict(id=1000 + i, stem=f"q{i}", nfp=int(2 + i % 4),
                   iota=0.1 + 0.02 * (i % 3), aspect=6.0,
                   symmetry_class=classes[i % 3], blanket="baseline", tau=tau,
                   C=float(C + rng.normal(0, 0.002)), anisotropy=0.9,
                   S_phi=float((3 * C * C - 1) / 2), lambda_phi=float(C * C),
                   reversal_frac=0.0,
                   angular_std_deg=float(30 * (1 - C)),
                   field_audit={"passed": True}, status="done",
                   delta_free_parallel=float(df), delta_scatter_parallel=float(ds))
        (d / f"config_{rec['id']}_baseline.json").write_text(json.dumps(rec))
    # blanket subset: one device (C=0.8) at 3 blankets/taus -> A varies, eta_source fixed
    for blk, tau in [("baseline", 13.0), ("thin", 8.0), ("thick", 20.0)]:
        C = 0.8; df = k1 * C
        rec = dict(id=2000, stem="qsub", nfp=3, iota=0.12, aspect=6.0,
                   symmetry_class="QA", blanket=blk, tau=tau, C=C, anisotropy=0.9,
                   S_phi=float((3 * C * C - 1) / 2), lambda_phi=float(C * C),
                   reversal_frac=0.0,
                   angular_std_deg=6.0, field_audit={"passed": True}, status="done",
                   delta_free_parallel=float(df),
                   delta_scatter_parallel=float(df * np.exp(-ktau * tau)))
        (d / f"config_2000_{blk}.json").write_text(json.dumps(rec))
    return d


def test_analyze_recovers_injected_law(tmp_path):
    d = _make_records(tmp_path / "recs")
    recs = az.load_records(d)
    anchor = az.pick_anchor(recs)
    assert anchor["C"] == pytest.approx(1.0, abs=0.01)      # highest-C is the anchor
    rows = az.assemble(recs, anchor)

    fit, _ = az.fit_eta_source_of_C(rows)
    assert fit["linear"]["slope"] == pytest.approx(1.0, abs=0.05)     # eta_source = C
    assert abs(fit["linear"]["intercept"]) < 0.05
    assert fit["linear"]["R2"] > 0.98

    uni = az.universality(rows)
    assert "universal" in uni["verdict"]                    # class-independent law

    sep = az.separability(rows)
    assert sep["A_of_tau"]["k"] == pytest.approx(0.05, abs=0.01)      # recover k
    assert sep["eta_source_blanket_spread"] < 1e-6          # eta_source invariant
    assert sep["failure_near_source_contamination"] is False


def test_figures_generate(tmp_path):
    d = _make_records(tmp_path / "recs")
    recs = az.load_records(d)
    rows = az.assemble(recs, az.pick_anchor(recs))
    fit, _ = az.fit_eta_source_of_C(rows)
    names = az.make_figures(rows, fit, tmp_path / "figs")
    assert "fig_C_predicts_eta.png" in names
    assert "fig_A_of_tau.png" in names
    assert "fig_S_predicts_eta.png" in names          # G1: the parity-correct predictor


# --------------------------------------------------------------------------- #
# G1: when eta_source ~ S_phi (parity-correct), analyze must PREFER S_phi over C
# --------------------------------------------------------------------------- #
def _make_S_records(tmp):
    """eta_source proportional to S_phi = (3 C^2 - 1)/2 (nonlinear in C), across the
    real-family-shaped C/class values -- so a linear S_phi fit beats any C fit."""
    d = Path(tmp); d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    devs = [(0.999, "QA"), (0.97, "QA"), (0.90, "QA"),
            (0.78, "QH"), (0.72, "QH"), (0.65, "QH")]
    for i, (C, cls) in enumerate(devs):
        lam = C * C; S = (3 * lam - 1) / 2
        df = 0.30 * S + rng.normal(0, 0.001)
        rec = dict(id=3000 + i, stem=f"s{i}", nfp=3, iota=1.0, aspect=6.0,
                   symmetry_class=cls, blanket="baseline", tau=13.0,
                   C=float(C), S_phi=float(S), lambda_phi=float(lam), reversal_frac=0.0,
                   anisotropy=0.03, angular_std_deg=float(np.degrees(np.arccos(C))),
                   field_audit={"passed": True}, status="done",
                   delta_free_parallel=float(df), delta_scatter_parallel=float(0.6 * df))
        (d / f"config_{3000 + i}_baseline.json").write_text(json.dumps(rec))
    return d


def test_S_phi_is_preferred_predictor(tmp_path):
    d = _make_S_records(tmp_path / "recs")
    rows = az.assemble(az.load_records(d), az.pick_anchor(az.load_records(d)))
    fitS, _ = az.fit_eta_source_of_S(rows)
    assert fitS["linear"]["R2"] > 0.99                        # tight linear collapse
    assert fitS["parameter_free"]["rms_resid"] < 0.02         # zero-free-param physics
    cmp = az.compare_predictors(rows)
    assert cmp["R2_S_phi_linear"] >= cmp["R2_C_linear"]       # S_phi at least as good
    assert cmp["best_predictor"] == "S_phi_linear"            # ... and wins outright
    uni = az.universality(rows)
    assert uni["by_S_phi"] is not None and "verdict" in uni
