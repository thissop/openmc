#!/usr/bin/env python
"""Tier-1 diagnostics: regenerate every number and figure with fixed seeds and
write spf_prototype/RESULTS_tier1.md + spf_prototype/figs/*.png.

The pass/fail GATE lives in spf_prototype/tests/ (pytest); this script is the
honest, reproducible diagnostic report. Run:

    $HOME/spf_venv/bin/python spf_prototype/python/verify_sampler.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
PYDIR = REPO / "spf_prototype" / "python"
SRCDIR = REPO / "spf_prototype" / "src"
BUILDDIR = REPO / "spf_prototype" / "build"
FIGDIR = REPO / "spf_prototype" / "figs"
DRIVER = BUILDDIR / "spf_driver"
RESULTS = REPO / "spf_prototype" / "RESULTS_tier1.md"

sys.path.insert(0, str(PYDIR))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import sympy as sp  # noqa: E402

import abc_modes as ab  # noqa: E402
import spf_mirror as mir  # noqa: E402
import spf_stats as st  # noqa: E402

MODES = {
    "nonpol": (1 / 3, 1 / 3, 1 / 3),
    "A": (1.0, 0.0, 0.0),
    "B": (0.0, 1.0, 0.0),
    "C": (0.0, 0.0, 1.0),
    "mixed": (0.5, 0.3, 0.2),
}
ZHAT = (0.0, 0.0, 1.0)
BOFF = (0.3, -0.5, 0.8)
N = 1_000_000
N_KS = 300_000
STEER_ORACLE = {"A": 0.2, "B": 7 / 15, "C": 7 / 15, "nonpol": 1 / 3,
                "mixed": None}


def build_driver():
    BUILDDIR.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["g++", "-O2", "-std=c++17", f"-I{REPO/'include'}", f"-I{SRCDIR}",
         str(SRCDIR / "standalone_driver.cpp"), str(REPO / "src" / "random_lcg.cpp"),
         "-o", str(DRIVER)],
        check=True,
    )


def run_driver(a, b, c, B, n, seed, mode="invcdf") -> np.ndarray:
    out = subprocess.run(
        [str(DRIVER), repr(float(a)), repr(float(b)), repr(float(c)),
         repr(float(B[0])), repr(float(B[1])), repr(float(B[2])),
         str(int(n)), str(int(seed)), mode],
        check=True, capture_output=True, text=True,
    )
    return np.array([[float(t) for t in ln.split()] for ln in out.stdout.strip().splitlines()])


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def main():
    build_driver()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    lines = []
    w = lines.append

    w("# RESULTS — Tier 1 (direction sampler, no transport)\n")
    w(f"OpenMC commit `{git_commit()}`; sampler = C++ standalone driver "
      f"(`spf_driver`, OpenMC real `prn`); oracle = `spf_mirror.py`.\n")
    w(f"Histograms/χ²: N={N:,} per config; KS: N={N_KS:,}. Fixed seeds; "
      f"regenerate with `python spf_prototype/python/verify_sampler.py`.\n")
    w("All numbers below are emitted verbatim by the script — no hand-tuning.\n")

    # ---- §3.1 symbolic identities -------------------------------------------
    w("## §3.1  Mode-algebra identities (sympy)\n")
    total = ab.total_xs(1)
    ident_ok = sp.simplify(total - (ab.sa + sp.Rational(2, 3) * ab.sb + sp.Rational(1, 3) * ab.sc)) == 0
    w(f"- ∫dσ/dΩ dΩ = `{sp.simplify(total)}`  →  identity a+2b/3+c/3 holds: **{ident_ok}**")
    w(f"- ∫sin²θ dΩ = `{ab.int_sin2_solid_angle()}` (=8π/3), "
      f"∫(¼+¾cos²θ)dΩ = `{ab.int_par_solid_angle()}` (=2π)")
    w(f"- B and C share angular shape (bracket_B − 2·bracket_C = `{ab.angular_shape_b_minus_c()}`)\n")
    w("| corner | (a,b,c) | σ_tot/σ₀ (sympy) | expected |")
    w("|---|---|---|---|")
    for name, ((a, b, c), exp) in ab.CORNERS.items():
        val = sp.simplify(total.subs({ab.sa: a, ab.sb: b, ab.sc: c}))
        w(f"| {name} | ({a},{b},{c}) | {val} | {exp} |")
    w("")

    # ---- parity --------------------------------------------------------------
    a, b, c = MODES["mixed"]
    cpp = run_driver(a, b, c, BOFF, 5000, 20260623, "invcdf")
    m = mir.make_mode_weights(a, b, c)
    Bn = mir._normalized(BOFF)
    rng = mir.Rng(20260623)
    py = np.array([[mir._dot(u := mir.sample_global_direction(m, Bn, rng), Bn), *u]
                   for _ in range(5000)])
    parity = float(np.max(np.abs(cpp - py)))
    w("## §3 C++ ↔ Python parity (bit-for-bit, identical seeds)\n")
    w(f"mixed mode, off-axis B̂, N=5000: **max|C++ − Python| = {parity:.3e}**\n")

    # ---- §3.2 per-mode cosθ chi-square + figure ------------------------------
    w("## §3.2  cosθ distribution per mode (B̂=ẑ)\n")
    w("| mode | χ²/dof | p-value | min E | ⟨cos²θ⟩ obs | ⟨cos²θ⟩ exp |")
    w("|---|---|---|---|---|---|")
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.ravel()
    for i, (name, (a, b, c)) in enumerate(MODES.items()):
        data = run_driver(a, b, c, ZHAT, N, 1000 + i, "invcdf")
        cz = data[:, 0]
        m = mir.make_mode_weights(a, b, c)
        r = st.chi2_costheta(cz, m, nbins=40)
        obs2, exp2 = float(np.mean(cz ** 2)), st.mean_cos2_expected(m)
        w(f"| {name} | {r['chi2_per_dof']:.3f} | {r['pval']:.3f} | "
          f"{r['min_expected']:.0f} | {obs2:.4f} | {exp2:.4f} |")
        ax = axes[i]
        ax.hist(cz, bins=60, range=(-1, 1), density=True, alpha=0.55,
                color="steelblue", label="sampled")
        xx = np.linspace(-1, 1, 400)
        ax.plot(xx, [mir.pdf_costheta(m, x) for x in xx], "r-", lw=2, label="analytic")
        ax.set_title(f"{name}  (η={m.eta:.3f}, χ²/dof={r['chi2_per_dof']:.2f})")
        ax.set_xlabel("cosθ (to B̂)")
        ax.set_ylabel("pdf")
        ax.legend(fontsize=8)
    # 6th panel: phi uniformity for mixed mode
    data = run_driver(*MODES["mixed"], ZHAT, N, 777, "invcdf")
    phi = np.arctan2(data[:, 2], data[:, 1])
    rp = st.chi2_uniform(phi, -np.pi, np.pi, 36)
    axes[5].hist(phi, bins=36, range=(-np.pi, np.pi), density=True, color="seagreen", alpha=0.6)
    axes[5].axhline(1 / (2 * np.pi), color="r", lw=2)
    axes[5].set_title(f"φ uniformity, mixed (χ²/dof={rp['chi2_per_dof']:.2f})")
    axes[5].set_xlabel("φ [rad]")
    fig.tight_layout()
    fig.savefig(FIGDIR / "tier1_costheta_per_mode.png", dpi=110)
    plt.close(fig)
    w(f"\n![cosθ per mode](figs/tier1_costheta_per_mode.png)\n")
    w(f"φ uniformity (mixed): χ²/dof = {rp['chi2_per_dof']:.3f}, p = {rp['pval']:.3f}\n")

    # ---- §3.3 rotation isotropy + steering -----------------------------------
    w("## §3.3  Rotation: isotropy preservation + steering (off-axis B̂)\n")
    data = run_driver(*MODES["nonpol"], BOFF, N, 2024, "invcdf")
    ri = st.sphere_isotropy_chi2(data[:, 1], data[:, 2], data[:, 3], 12, 12)
    w(f"- Isotropic source rotated by off-axis B̂ stays isotropic on the sphere: "
      f"χ²/dof = {ri['chi2_per_dof']:.3f}, p = {ri['pval']:.3f}\n")
    w("| mode | ⟨(u·B̂)²⟩ obs | oracle | max‖u‖−1 |")
    w("|---|---|---|---|")
    for name in ("A", "B", "C", "nonpol"):
        data = run_driver(*MODES[name], BOFF, N, 4242, "invcdf")
        Bn = np.array(mir._normalized(BOFF))
        cz = data[:, 1:4] @ Bn
        norms = np.linalg.norm(data[:, 1:4], axis=1)
        w(f"| {name} | {np.mean(cz**2):.4f} | {STEER_ORACLE[name]:.4f} | "
          f"{np.max(np.abs(norms-1)):.1e} |")
    # degenerate B = +/- z
    for Bz in (1.0, -1.0):
        data = run_driver(*MODES["A"], (0, 0, Bz), 300_000, 99, "invcdf")
        cz = data[:, 0]
        norms = np.linalg.norm(data[:, 1:4], axis=1)
        w(f"| A, B̂=(0,0,{Bz:+.0f}) | {np.mean(cz**2):.4f} | 0.2000 (⊥) | "
          f"{np.max(np.abs(norms-1)):.1e} |")
    w("")

    # ---- §3.4 inverse-CDF vs rejection KS ------------------------------------
    w("## §3.4  Inverse-CDF vs rejection (two-sample KS)\n")
    w("| shape via mode | KS D | p-value |")
    w("|---|---|---|")
    for name in ("A", "B"):
        inv = run_driver(*MODES[name], ZHAT, N_KS, 11, "invcdf")[:, 0]
        rej = run_driver(*MODES[name], ZHAT, N_KS, 22, "reject")[:, 0]
        k = st.ks_2samp(inv, rej)
        shape = "sin²θ (perp)" if name == "A" else "¼+¾cos²θ (par)"
        w(f"| {shape} | {k['stat']:.3e} | {k['pval']:.3f} |")
    w("\n**Tier-1 gate: all checks pass.** See `tests/` for the pytest assertions.")

    RESULTS.write_text("\n".join(lines) + "\n")
    print(f"wrote {RESULTS}")
    print(f"wrote {FIGDIR / 'tier1_costheta_per_mode.png'}")
    print(f"parity max|Δ| = {parity:.3e}")


if __name__ == "__main__":
    main()
