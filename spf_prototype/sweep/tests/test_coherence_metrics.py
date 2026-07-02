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
