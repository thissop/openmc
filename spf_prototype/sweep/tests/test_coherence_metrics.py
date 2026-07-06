"""Analytic unit tests for coherence_metrics -- every target is hand-derived.

Gates (from the brief):
 (a) uniform field  -> C = 1 exactly
 (b) field rotating uniformly by total angle Phi -> C = sin(Phi/2)/(Phi/2) (sinc)
 (c) source weighting actually changes C when s is nonuniform
plus: the cylindrical-frame calibration (toroidal field -> C=1; uniform lab field
-> C<1 in cyl frame), direction-tensor eigenvalues, angular std, source_weights.
"""
import numpy as np
import pytest

import coherence_metrics as cm


# --------------------------------------------------------------------------- #
# (a) uniform field -> C = 1 exactly
# --------------------------------------------------------------------------- #
def test_uniform_lab_field_C_is_one():
    N = 500
    b = np.tile([0.0, 0.0, 1.0], (N, 1))
    assert cm.coherence_C(b, frame="lab") == pytest.approx(1.0, abs=1e-12)


def test_purely_toroidal_field_C_is_one_in_cyl_frame():
    # tokamak limit: b = e_phi at every point -> (0,1,0) in local cylindrical basis
    phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    pos = np.stack([np.cos(phi), np.sin(phi), 0.2 * np.sin(3 * phi)], axis=1)
    b = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=1)
    m = cm.coherence_metrics(b, pos=pos, frame="cylindrical")
    assert m["C"] == pytest.approx(1.0, abs=1e-12)
    assert m["tensor_evals"] == pytest.approx([1.0, 0.0, 0.0], abs=1e-9)
    assert m["angular_std"] == pytest.approx(0.0, abs=1e-6)  # arccos-near-1 roundoff


def test_uniform_lab_field_is_incoherent_in_cyl_frame():
    # a uniform lab-x field winds through the cylindrical basis -> C ~ 0 over 2pi
    phi = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
    pos = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=1)
    b = np.tile([1.0, 0.0, 0.0], (len(phi), 1))
    assert cm.coherence_C(b, pos=pos, frame="cylindrical") < 1e-2
    # ... but the SAME field is perfectly coherent in the lab frame
    assert cm.coherence_C(b, frame="lab") == pytest.approx(1.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# (b) uniformly rotating field -> sinc law
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("Phi", [0.5, 1.0, np.pi, 2.0])
def test_rotating_field_sinc_law(Phi):
    N = 20000
    a = Phi * (np.arange(N) + 0.5) / N            # midpoint rule on [0, Phi]
    b = np.stack([np.cos(a), np.sin(a), np.zeros_like(a)], axis=1)
    C = cm.coherence_C(b, frame="lab")
    C_exact = np.sin(Phi / 2) / (Phi / 2)         # |<e^{i a}>| over uniform [0,Phi]
    assert C == pytest.approx(C_exact, abs=2e-3)


def test_rotation_limits():
    # Phi -> 0 gives C -> 1; Phi = pi gives C = 2/pi
    small = np.linspace(0, 1e-3, 100)
    b0 = np.stack([np.cos(small), np.sin(small), np.zeros_like(small)], axis=1)
    assert cm.coherence_C(b0, frame="lab") == pytest.approx(1.0, abs=1e-6)
    a = np.pi * (np.arange(20000) + 0.5) / 20000
    bpi = np.stack([np.cos(a), np.sin(a), np.zeros_like(a)], axis=1)
    assert cm.coherence_C(bpi, frame="lab") == pytest.approx(2 / np.pi, abs=2e-3)


# --------------------------------------------------------------------------- #
# (c) source weighting changes C
# --------------------------------------------------------------------------- #
def test_source_weighting_changes_C():
    # two clusters: A=(1,0,0), B=(0,1,0)
    b = np.array([[1.0, 0, 0], [0, 1.0, 0]])
    C_equal = cm.coherence_C(b, weights=[1.0, 1.0], frame="lab")
    C_skew = cm.coherence_C(b, weights=[3.0, 1.0], frame="lab")
    assert C_equal == pytest.approx(np.sqrt(0.5), abs=1e-12)      # |(.5,.5)|
    assert C_skew == pytest.approx(np.sqrt(0.75 ** 2 + 0.25 ** 2), abs=1e-12)
    assert abs(C_skew - C_equal) > 0.05                          # weighting matters


def test_direction_tensor_two_clusters():
    b = np.array([[1.0, 0, 0], [0, 1.0, 0]])
    T, evals, _ = cm.direction_tensor(b, weights=[1.0, 1.0], frame="lab")
    assert np.trace(T) == pytest.approx(1.0, abs=1e-12)
    assert evals == pytest.approx([0.5, 0.5, 0.0], abs=1e-12)


def test_anisotropy_planar_fan_vs_isotropic():
    # planar fan (all in-plane) has a large l1-l2 gap only if concentrated; an
    # isotropic 3D smear has evals ~ (1/3,1/3,1/3). Different structure, and the
    # tensor tells them apart even when C is similar.
    rng = np.random.default_rng(0)
    iso = rng.normal(size=(4000, 3)); iso /= np.linalg.norm(iso, axis=1, keepdims=True)
    m_iso = cm.coherence_metrics(iso, frame="lab")
    assert m_iso["tensor_evals"] == pytest.approx([1 / 3, 1 / 3, 1 / 3], abs=3e-2)
    assert m_iso["C"] < 0.05


# --------------------------------------------------------------------------- #
# source weights + angular std
# --------------------------------------------------------------------------- #
def test_source_weights_profile():
    assert cm.source_weights(0.0) == pytest.approx(1.0)
    assert cm.source_weights(1.0) == pytest.approx(0.0)
    assert cm.source_weights(np.array([0.5]))[0] == pytest.approx(0.75)
    assert cm.source_weights(0.5, model="n2sigmav") == pytest.approx(0.5625)


def test_angular_std_zero_for_uniform():
    b = np.tile([0.3, 0.4, np.sqrt(1 - 0.25)], (100, 1))
    m = cm.coherence_metrics(b, frame="lab")
    assert m["angular_std"] == pytest.approx(0.0, abs=1e-9)
    assert m["circular_std"] == pytest.approx(0.0, abs=1e-6)


def test_invalid_inputs_raise():
    with pytest.raises(ValueError):
        cm.coherence_C(np.zeros((3, 3)), frame="cylindrical")   # needs pos
    with pytest.raises(ValueError):
        cm.coherence_C(np.array([[0.0, 0, 0]]), frame="lab")    # zero vector
    with pytest.raises(ValueError):
        cm.coherence_C(np.array([[1.0, 0, 0]]), weights=[-1.0], frame="lab")


# --------------------------------------------------------------------------- #
# G1: nematic order S_phi / lambda_phi and the PARITY property (headless kernel)
# --------------------------------------------------------------------------- #
def _cap_about_ephi(half_angle, n=40000, seed=0):
    """n unit vectors in a cap of `half_angle` about the toroidal slot (column 1,
    e_phi), azimuthally uniform, drawn uniform in cos(theta) over [cos(HA), 1].
    Returned as raw (N,3) so coherence_metrics(..., frame='lab') treats column 1 as
    the toroidal axis directly (lambda_phi = T[1,1] then measures the cap width)."""
    rng = np.random.default_rng(seed)
    cos_t = 1.0 - rng.random(n) * (1.0 - np.cos(half_angle))   # cos in [cos(HA), 1]
    sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, 1.0))
    az = rng.random(n) * 2.0 * np.pi
    return np.stack([sin_t * np.cos(az), cos_t, sin_t * np.sin(az)], axis=1)


def test_Sphi_limits_aligned_and_isotropic():
    aligned = np.tile([0.0, 1.0, 0.0], (500, 1))               # all e_phi
    m = cm.coherence_metrics(aligned, frame="lab")
    assert m["lambda_phi"] == pytest.approx(1.0, abs=1e-12)
    assert m["S_phi"] == pytest.approx(1.0, abs=1e-12)         # aligned -> full order
    assert m["reversal_frac"] == pytest.approx(0.0, abs=1e-12)
    rng = np.random.default_rng(1)
    iso = rng.normal(size=(8000, 3)); iso /= np.linalg.norm(iso, axis=1, keepdims=True)
    mi = cm.coherence_metrics(iso, frame="lab")
    assert mi["lambda_phi"] == pytest.approx(1 / 3, abs=3e-2)  # isotropic -> 1/3
    assert mi["S_phi"] == pytest.approx(0.0, abs=5e-2)         # -> 0 (no steering)


def test_cap_sandwich_C2_le_lamphi_le_C():
    # symmetric cap about e_phi: perp means ~0 so C = <cos_t>, lambda_phi = <cos_t^2>,
    # giving C^2 <= lambda_phi <= C (Jensen, and cos<=1 on the cap). No reversal.
    m = cm.coherence_metrics(_cap_about_ephi(np.radians(50)), frame="lab")
    assert m["C"] ** 2 <= m["lambda_phi"] + 1e-6
    assert m["lambda_phi"] <= m["C"] + 1e-6
    assert m["reversal_frac"] == pytest.approx(0.0, abs=1e-12)


def test_parity_Sphi_invariant_but_C_collapses_under_bhat_flip():
    """THE parity argument, executable: the SPF kernel is EVEN in b_hat, so eta is a
    functional of the SECOND moment only. Flipping b_hat -> -b_hat on half the source
    must leave lambda_phi/S_phi/tensor_evals UNCHANGED while the first moment C
    collapses -- proving C is not the causal variable and S_phi is."""
    cap = _cap_about_ephi(np.radians(30), n=40000, seed=2)
    m0 = cm.coherence_metrics(cap, frame="lab")
    flipped = cap.copy(); flipped[::2] *= -1.0                 # flip exactly half
    m1 = cm.coherence_metrics(flipped, frame="lab")
    # 2nd-moment predictors are invariant under the sign flip ...
    assert m1["lambda_phi"] == pytest.approx(m0["lambda_phi"], abs=1e-9)
    assert m1["S_phi"] == pytest.approx(m0["S_phi"], abs=1e-9)
    assert m1["tensor_evals"] == pytest.approx(m0["tensor_evals"], abs=1e-9)
    # ... but the first moment C collapses (half the vectors now oppose) ...
    assert m0["C"] > 0.9 and m1["C"] < 0.1
    # ... and half the source now reads as reversed toroidal sense.
    assert m1["reversal_frac"] == pytest.approx(0.5, abs=1e-9)
