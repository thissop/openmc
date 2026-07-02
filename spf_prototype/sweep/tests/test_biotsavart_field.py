"""Biot-Savart kernel vs the analytic circular loop, and the field-map round-trip."""
import numpy as np
import pytest

import biotsavart_field as bs


def test_loop_onaxis_direction():
    coils = bs.circular_loop(radius=1.0, current=1.0, n=800)
    z = np.array([0.0, 0.2, 0.5, 1.0, 2.0])
    pts = np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1)
    B = bs.biot_savart_B(coils, pts)
    bh = B / np.linalg.norm(B, axis=1, keepdims=True)
    assert np.allclose(bh, np.tile([0, 0, 1.0], (len(z), 1)), atol=1e-6)  # +z on axis
    assert np.allclose(B[:, :2], 0.0, atol=1e-6)                          # no transverse


def test_loop_onaxis_magnitude_matches_analytic():
    coils = bs.circular_loop(radius=1.0, current=1.0, n=1200)
    z = np.array([0.0, 0.3, 0.7, 1.5])
    pts = np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1)
    Bz = bs.biot_savart_B(coils, pts)[:, 2]
    assert np.allclose(Bz, bs.loop_Bz_onaxis(z, 1.0), rtol=1e-3)


def test_loop_current_and_orientation_signs():
    # reversing current flips the field; reversing radius sign is nonsense but
    # a loop about x-axis points its on-axis field along +x
    z0 = np.array([[0.0, 0.0, 0.5]])
    Bpos = bs.biot_savart_B(bs.circular_loop(1.0, +1.0, 600), z0)[0, 2]
    Bneg = bs.biot_savart_B(bs.circular_loop(1.0, -1.0, 600), z0)[0, 2]
    assert Bpos == pytest.approx(-Bneg, rel=1e-9)
    xax = bs.circular_loop(1.0, 1.0, 600, axis="x")
    Bx = bs.biot_savart_B(xax, np.array([[0.5, 0.0, 0.0]]))[0]
    assert Bx[0] > 0 and abs(Bx[1]) < 1e-6 and abs(Bx[2]) < 1e-6


def test_fieldmap_roundtrip():
    # write a loop field to spf_fieldmap_v1, read it back with the .so's Python
    # mirror, and confirm b_hat agrees at a grid point (interop contract).
    import tempfile
    from pathlib import Path
    import fieldmap  # ../python, via conftest
    coils = bs.circular_loop(radius=1.0, current=1.0, n=800)
    field = bs.CoilField(coils)
    with tempfile.TemporaryDirectory() as d:
        stem = str(Path(d) / "loopmap")
        Rmin, Rmax, Zmin, Zmax = 0.2, 2.0, -1.0, 1.0
        nR, nphi, nZ = 24, 32, 24
        bs.write_fieldmap(field, dict(R_min=Rmin, R_max=Rmax, Z_min=Zmin, Z_max=Zmax),
                          dict(nR=nR, nphi=nphi, nZ=nZ), stem, nfp=1)
        fm = fieldmap.FieldMapField(stem)
        # evaluate at an EXACT grid node (no interpolation) -> validates the byte
        # layout / index order to machine precision
        R = Rmin + 12 * (Rmax - Rmin) / (nR - 1)
        Z = Zmin + 15 * (Zmax - Zmin) / (nZ - 1)
        p = (R, 0.0, Z)                      # phi node 0
        got = np.array(fm.bhat(p))
        exp = field.bhat(np.array([p]))[0]
        assert np.allclose(got, exp, atol=1e-9)
