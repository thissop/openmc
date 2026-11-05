"""Correctness gates for the representation harness. Run: pytest -q test_fluxrep.py

Uses the local QH wout if present; otherwise skips (the harness logic is still exercised by
the synthetic-torus tests, which need no data file)."""
import os
import numpy as np
import pytest

from fluxrep import read_wout, SpectralRep, GridRep
from validate_vmec_bhat import validate_wout
import metrics

# candidate local wout locations (repo layout may vary)
_WOUT_CANDIDATES = [
    os.path.expanduser("~/Documents/GitHub/openmc/spf_prototype/talk/geom/wout_qh.nc"),
    "spf_prototype/talk/geom/wout_qh.nc",
]
WOUT = next((p for p in _WOUT_CANDIDATES if os.path.exists(p)), None)


def _analytic_torus_wout(R0=6.0, a=1.5, nfp=1, ns=64):
    """A circular axisymmetric torus as a minimal VMEC-like wout dict: R = R0 + a s^.5 cos th,
    Z = a s^.5 sin th. Two modes (m=0,n=0)->R0 and (m=1,n=0). Exact volume = 2 pi^2 R0 a^2."""
    s = np.linspace(0.0, 1.0, ns)
    xm = np.array([0.0, 1.0]); xn = np.array([0.0, 0.0])
    rmnc = np.zeros((ns, 2)); zmns = np.zeros((ns, 2))
    rmnc[:, 0] = R0
    rmnc[:, 1] = a * np.sqrt(s)
    zmns[:, 1] = a * np.sqrt(s)
    return dict(xm=xm, xn=xn, nfp=nfp, ns=ns, s=s, rmnc=rmnc, zmns=zmns,
               rmns=None, zmnc=None, iotaf=None, lmns=None), R0, a


def test_sqrtg_integrates_to_volume():
    """|sqrt(g)| integrated over (s,theta,zeta) must recover the analytic torus volume."""
    wout, R0, a = _analytic_torus_wout()
    spec = SpectralRep(wout)
    ns, nt, nz = 200, 256, 4
    s = np.linspace(1e-6, 1, ns); th = np.linspace(0, 2*np.pi, nt, endpoint=False)
    ze = np.linspace(0, 2*np.pi, nz, endpoint=False)
    S, TH, ZE = np.meshgrid(s, th, ze, indexing="ij")
    _, _, G = spec.eval(S.ravel(), TH.ravel(), ZE.ravel())
    dV = (1.0/ (ns-1)) * (2*np.pi/nt) * (2*np.pi/nz)
    V = np.sum(G) * dV
    V_true = 2 * np.pi**2 * R0 * a**2
    assert abs(V - V_true) / V_true < 2e-3, (V, V_true)


def test_spectral_exact_at_nodes():
    """Grid filled from spectral must equal spectral exactly AT its own nodes (no interp error)."""
    wout, _, _ = _analytic_torus_wout()
    spec = SpectralRep(wout)
    grid = GridRep(spec, nrho=10, ntheta=32, nzeta=8, nfp=1, s_min=1e-6)
    s_nodes, th_nodes, ze_nodes = grid._grid
    s = s_nodes[3]; th = th_nodes[5]; ze = ze_nodes[2]
    Rs, Zs, Gs = spec.eval([s], [th], [ze])
    Rg, Zg, Gg = grid.eval([s], [th], [ze])
    assert abs(Rg[0]-Rs[0]) < 1e-9 and abs(Zg[0]-Zs[0]) < 1e-9 and abs(Gg[0]-Gs[0]) < 1e-6


def test_grid_accuracy_improves_with_resolution():
    """Coarse->fine grid must monotonically reduce the R interpolation error vs spectral."""
    wout, _, _ = _analytic_torus_wout()
    spec = SpectralRep(wout)
    coarse = GridRep(spec, 6, 16, 4, nfp=1)
    fine = GridRep(spec, 24, 96, 4, nfp=1)
    a_c = metrics.accuracy(spec, coarse, n=5000, nfp=1)["R"]["rms"]
    a_f = metrics.accuracy(spec, fine, n=5000, nfp=1)["R"]["rms"]
    assert a_f < a_c


def test_source_pdf_and_periodic_seams():
    """Source-weight metrics are bounded and interpolation wraps at both seams."""
    wout, _, _ = _analytic_torus_wout()
    spec = SpectralRep(wout)
    grid = GridRep(spec, 12, 32, 8, nfp=1)
    result = metrics.accuracy(spec, grid, n=2000, nfp=1)
    assert 0.0 <= result["source_pdf"]["total_variation"] <= 1.0
    assert 0.0 <= result["source_pdf"]["hellinger"] <= 1.0

    s = np.array([0.5, 0.5])
    theta = np.array([1e-12, 2 * np.pi - 1e-12])
    zeta = np.array([1e-12, 2 * np.pi - 1e-12])
    R, Z, G = grid.eval(s, theta, zeta)
    np.testing.assert_allclose(R[0], R[1], atol=1e-9)
    np.testing.assert_allclose(Z[0], Z[1], atol=1e-9)
    np.testing.assert_allclose(G[0], G[1], atol=1e-7)


@pytest.mark.skipif(WOUT is None, reason="local QH wout not found")
def test_real_wout_loads_and_evaluates():
    wout = read_wout(WOUT)
    assert wout["iotaf"] is not None and wout["lmns"] is not None  # iota + lambda present for B-hat backend
    spec = SpectralRep(wout)
    R, Z, G = spec.eval([0.5, 0.9], [0.1, 1.0], [0.0, 0.3])
    assert np.all(np.isfinite(R)) and np.all(G > 0)


@pytest.mark.skipif(WOUT is None, reason="local QH wout not found")
def test_vmec_bhat_matches_exported_components():
    result = validate_wout(WOUT, n=1000)
    assert result["p99_degrees"] < 0.05, result
