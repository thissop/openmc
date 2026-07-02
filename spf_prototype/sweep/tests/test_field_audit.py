"""field_audit hard gates on synthetic good/bad fields (network-free) + a
network-optional end-to-end check on a cached real QUASR device."""
import numpy as np
import pytest

import field_audit as fa


# --- synthetic fields ---------------------------------------------------------
class Toroidal:
    """B = e_phi: divergence-free and tangent to any axisymmetric torus -> GOOD."""
    def B(self, pts):
        x, y = pts[:, 0], pts[:, 1]
        r = np.hypot(x, y) + 1e-30
        return np.stack([-y / r, x / r, np.zeros_like(x)], axis=1)
    def bhat(self, pts):
        return self.B(pts)


class Radial:
    """B = e_R: div != 0 AND large normal component on a torus -> BAD."""
    def B(self, pts):
        x, y = pts[:, 0], pts[:, 1]
        r = np.hypot(x, y) + 1e-30
        return np.stack([x / r, y / r, np.zeros_like(x)], axis=1)
    def bhat(self, pts):
        return self.B(pts)


def _torus_grid(R0=1.0, a=0.3, nt=48, nph=64):
    th = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R = R0 + a * np.cos(TH); Z = a * np.sin(TH)
    xyz = np.stack([R * np.cos(PH), R * np.sin(PH), Z], axis=-1)
    return dict(xyz=xyz)


def _interior_pts(R0=1.0, a=0.15, nt=10, nph=12):
    th = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R = R0 + a * np.cos(TH); Z = a * np.sin(TH)
    return np.stack([R * np.cos(PH), R * np.sin(PH), Z], axis=-1).reshape(-1, 3)


# --- good field passes every check -------------------------------------------
def test_good_field_passes():
    f = Toroidal()
    pts = _interior_pts()
    assert fa.check_unit(f, pts)["passed"]
    assert fa.check_divergence(f, pts, h=1e-4, L=0.3)["passed"]
    assert fa.check_boundary_flux(f, _torus_grid())["passed"]


# --- deliberately-bad field fails the gate -----------------------------------
def test_bad_field_fails_boundary_and_divergence():
    f = Radial()
    bnd = fa.check_boundary_flux(f, _torus_grid())
    div = fa.check_divergence(f, _interior_pts(), h=1e-4, L=0.3)
    assert not bnd["passed"]           # e_R has O(1) normal component on the torus
    assert not div["passed"]           # div e_R = 1/r != 0
    assert bnd["value"] > 0.5          # and it fails by a wide margin, not a hair


def test_nonunit_field_flagged():
    class Weak:
        def B(self, pts):
            return 0.5 * Toroidal().B(pts)
        def bhat(self, pts):
            return 0.5 * Toroidal().B(pts)   # deliberately NOT unit
    assert not fa.check_unit(Weak(), _interior_pts())["passed"]


# --- network-optional end-to-end on a real cached device ---------------------
def test_real_device_audit_if_cached():
    from pathlib import Path
    import quasr_loader as ql
    serial = ql.SERIAL_DIR / "serial0059509.json"
    if not serial.exists():
        pytest.skip("QUASR serial not cached; run quasr_loader.py 59509 first")
    from biotsavart_field import CoilField
    dev = ql.load_device(59509)
    fld = CoilField(dev.coils, meta=dev.meta)
    passed, report = fa.audit_coilfield(fld, dev)
    assert passed
    bflux = [c for c in report["checks"] if c["name"] == "boundary_flux"][0]
    assert bflux["value"] < 1e-2       # coils reproduce the flux surface to ~1%
