"""Tier-1 §3.2-3.4 -- statistical verification of the C++ direction sampler.

The C++ standalone driver is the sampler-under-test; spf_mirror is the oracle.
All draws use OpenMC's real prn so these numbers reflect the production RNG.
Run with `-s` to see the printed chi-square / KS / parity tables.
"""
import numpy as np
import pytest

import spf_mirror as mir
import spf_stats as st

MODES = {
    "nonpol": (1 / 3, 1 / 3, 1 / 3),
    "A": (1.0, 0.0, 0.0),
    "B": (0.0, 1.0, 0.0),
    "C": (0.0, 0.0, 1.0),
    "mixed": (0.5, 0.3, 0.2),
}
ZHAT = (0.0, 0.0, 1.0)
BOFF = (0.3, -0.5, 0.8)  # generic non-axis field for rotation/steering tests

N_CHI2 = 500_000
N_KS = 200_000


def test_cpp_python_parity(run_driver):
    """C++ driver vs Python mirror on identical seeds -> bit-for-bit (≤1e-12)."""
    a, b, c = MODES["mixed"]
    seed = 20260623
    N = 4000
    cpp = run_driver(a, b, c, BOFF, N, seed, "invcdf")

    m = mir.make_mode_weights(a, b, c)
    Bn = mir._normalized(BOFF)
    rng = mir.Rng(seed)
    rows = []
    for _ in range(N):
        u = mir.sample_global_direction(m, Bn, rng, use_rejection=False)
        rows.append((mir._dot(u, Bn), u[0], u[1], u[2]))
    py = np.array(rows)

    max_abs = float(np.max(np.abs(cpp - py)))
    print(f"\n[parity] mixed mode, N={N}, max|C++ - Python| = {max_abs:.3e}")
    assert max_abs <= 1e-12


@pytest.mark.parametrize("name", list(MODES))
def test_costheta_chi2(run_driver, name):
    """Sampled cosθ (B̂=ẑ) matches the analytic marginal pdf per mode."""
    a, b, c = MODES[name]
    seed = 1000 + list(MODES).index(name)  # deterministic per mode
    data = run_driver(a, b, c, ZHAT, N_CHI2, seed, mode="invcdf")
    cz = data[:, 0]
    m = mir.make_mode_weights(a, b, c)
    res = st.chi2_costheta(cz, m, nbins=40)
    exp_cos2 = st.mean_cos2_expected(m)
    obs_cos2 = float(np.mean(cz ** 2))
    print(f"\n[chi2 cosθ] {name:7s} chi2/dof={res['chi2_per_dof']:.3f} "
          f"p={res['pval']:.3f} min_exp={res['min_expected']:.1f} "
          f"⟨cos²⟩ obs={obs_cos2:.4f} exp={exp_cos2:.4f}")
    assert res["min_expected"] > 5.0
    assert res["chi2_per_dof"] < 2.0
    assert res["pval"] > 1e-4
    assert abs(obs_cos2 - exp_cos2) < 4e-3


def test_phi_uniform(run_driver):
    """Azimuthal angle is uniform on [-π, π) (B̂=ẑ)."""
    a, b, c = MODES["mixed"]
    data = run_driver(a, b, c, ZHAT, N_CHI2, seed=777, mode="invcdf")
    phi = np.arctan2(data[:, 2], data[:, 1])
    res = st.chi2_uniform(phi, -np.pi, np.pi, nbins=36)
    print(f"\n[phi uniform] chi2/dof={res['chi2_per_dof']:.3f} p={res['pval']:.3f}")
    assert res["chi2_per_dof"] < 2.0
    assert res["pval"] > 1e-4


def test_rotation_preserves_isotropy(run_driver):
    """Isotropic source (nonpol is flat in cosθ) rotated by a non-axis B̂ stays
    isotropic on the sphere -- catches rotation bias."""
    a, b, c = MODES["nonpol"]
    data = run_driver(a, b, c, BOFF, N_CHI2, seed=2024, mode="invcdf")
    res = st.sphere_isotropy_chi2(data[:, 1], data[:, 2], data[:, 3], n_cos=12, n_phi=12)
    print(f"\n[rotation isotropy] chi2/dof={res['chi2_per_dof']:.3f} p={res['pval']:.3f}")
    assert res["chi2_per_dof"] < 2.0
    assert res["pval"] > 1e-4


@pytest.mark.parametrize("name,expected", [("A", 0.2), ("B", 7 / 15), ("C", 7 / 15),
                                           ("nonpol", 1 / 3)])
def test_steering_about_offaxis_field(run_driver, name, expected):
    """With a non-axis B̂, ⟨(u·B̂)²⟩ matches the per-mode oracle: A peaks ⊥ B̂
    (0.2), B/C peak ∥ B̂ (7/15), nonpol isotropic (1/3). Confirms the global lobe
    is steered correctly by the Gram-Schmidt rotation."""
    a, b, c = MODES[name]
    data = run_driver(a, b, c, BOFF, N_CHI2, seed=4242, mode="invcdf")
    ux, uy, uz = data[:, 1], data[:, 2], data[:, 3]
    Bn = np.array(mir._normalized(BOFF))
    cz = ux * Bn[0] + uy * Bn[1] + uz * Bn[2]
    cos2 = float(np.mean(cz ** 2))
    norms = np.sqrt(ux ** 2 + uy ** 2 + uz ** 2)
    print(f"\n[steering] {name:7s} ⟨(u·B̂)²⟩ obs={cos2:.4f} exp={expected:.4f} "
          f"max|‖u‖-1|={np.max(np.abs(norms - 1)):.2e}")
    assert np.max(np.abs(norms - 1.0)) < 1e-12  # unit directions
    assert abs(cos2 - expected) < 4e-3


@pytest.mark.parametrize("Bz", [1.0, -1.0])
def test_degenerate_field_along_z(run_driver, Bz):
    """B̂ = ±ẑ exercises the Gram-Schmidt least-aligned-axis branch; must not blow
    up and must still steer pure-A perpendicular to B̂ (⟨cos²⟩→0.2)."""
    a, b, c = MODES["A"]
    data = run_driver(a, b, c, (0.0, 0.0, Bz), 200_000, seed=99, mode="invcdf")
    cz = data[:, 0]  # = u·B̂ from the driver
    norms = np.sqrt(np.sum(data[:, 1:4] ** 2, axis=1))
    print(f"\n[degenerate Bz={Bz:+.0f}] ⟨cos²⟩={np.mean(cz**2):.4f} (exp 0.2) "
          f"max|‖u‖-1|={np.max(np.abs(norms-1)):.2e}")
    assert np.max(np.abs(norms - 1.0)) < 1e-12
    assert abs(np.mean(cz ** 2) - 0.2) < 5e-3


@pytest.mark.parametrize("name", ["A", "B"])
def test_invcdf_vs_rejection_ks(run_driver, name):
    """Inverse-CDF and rejection samplers produce statistically identical cosθ
    (two-sample KS) -- validates the closed-form cubic root selection."""
    a, b, c = MODES[name]
    inv = run_driver(a, b, c, ZHAT, N_KS, seed=11, mode="invcdf")[:, 0]
    rej = run_driver(a, b, c, ZHAT, N_KS, seed=22, mode="reject")[:, 0]
    res = st.ks_2samp(inv, rej)
    print(f"\n[KS invcdf vs reject] {name}: D={res['stat']:.4e} p={res['pval']:.3f}")
    assert res["pval"] > 1e-3
