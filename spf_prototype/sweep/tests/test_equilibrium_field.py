"""G3: EquilibriumField reads the portable spf_fieldmap_v1 map and plugs into the audit,
with the div-B gate skipped by construction. A synthetic TOROIDAL map has known answers
(C = 1 in the cylindrical frame; B.n = 0 on any axisymmetric torus), so no DESC needed."""
from pathlib import Path

import numpy as np
import pytest

import equilibrium_field as ef
import coherence_metrics as cm
import field_audit as fa


def _write_toroidal_map(stem, nR=8, nphi=48, nZ=8,
                        Rmin=50.0, Rmax=150.0, Zmin=-50.0, Zmax=50.0):
    """A purely toroidal field b_hat = (-sin phi, cos phi, 0) on the (R,phi,Z) grid."""
    phi = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
    _, Pg, _ = np.meshgrid(np.linspace(Rmin, Rmax, nR), phi,
                           np.linspace(Zmin, Zmax, nZ), indexing="ij")
    bx = (-np.sin(Pg)).astype("<f8"); by = np.cos(Pg).astype("<f8")
    bz = np.zeros_like(bx)
    with open(str(stem) + ".bin", "wb") as f:
        f.write(bx.tobytes()); f.write(by.tobytes()); f.write(bz.tobytes())
    Path(str(stem) + ".meta").write_text(
        "magic spf_fieldmap_v1\nndim 3\n"
        f"nR {nR}\nnphi {nphi}\nnZ {nZ}\n"
        f"R_min {float(Rmin)!r}\nR_max {float(Rmax)!r}\n"
        f"phi_min 0.0\nperiod {2 * np.pi!r}\n"
        f"Z_min {float(Zmin)!r}\nZ_max {float(Zmax)!r}\n"
        "nfp 1\nendian little\nsource test:toroidal\n")


class _FakeDevice:
    """Minimal QuasrDevice stand-in: an axisymmetric circular torus (R0=100, a=20 cm)."""
    ID = 999
    meta = {"minor_radius": 20.0, "R0": 100.0}

    class _Dev:
        def RZ(self, TH, PH, rho):
            return 100.0 + 20.0 * rho * np.cos(TH), 20.0 * rho * np.sin(TH)
    device = _Dev()

    def boundary_grid(self, ntheta=64, nphi=128, rho=1.0):
        th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
        ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        R, Z = self.device.RZ(TH, PH, rho)
        return dict(theta=th, phi=ph, R=R, Z=Z,
                    xyz=np.stack([R * np.cos(PH), R * np.sin(PH), Z], axis=-1))


def test_equilibrium_field_toroidal_C_is_one(tmp_path):
    stem = tmp_path / "equil_test"
    _write_toroidal_map(stem)
    fld = ef.EquilibriumField(stem)
    ph = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    pts = np.stack([100 * np.cos(ph), 100 * np.sin(ph), np.zeros_like(ph)], axis=1)
    bh = fld.bhat(pts)
    assert np.allclose(np.linalg.norm(bh, axis=1), 1.0, atol=1e-9)      # unit
    assert cm.coherence_C(bh, pos=pts, frame="cylindrical") == pytest.approx(1.0, abs=1e-6)


def test_audit_skips_div_and_passes_for_equilibrium(tmp_path):
    stem = tmp_path / "equil_test"
    _write_toroidal_map(stem)
    fld = ef.EquilibriumField(stem)
    passed, report = fa.audit_coilfield(fld, _FakeDevice(), ntheta=32, nphi=64)
    assert passed                                                       # B.n = 0 on the torus
    div = [c for c in report["checks"] if c["name"] == "div_free"][0]
    assert div["passed"] and str(div.get("note", "")).startswith("skipped")


def test_missing_map_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ef.EquilibriumField(tmp_path / "does_not_exist")
