"""Local smoke test of gil_to_makegrid WITHOUT the 481 MB Gil archive or the cluster.

We can't ship the Gil coils here, but the converter's only real dependency is that the archive is a
simsopt BiotSavart JSON whose .coils each expose curve.gamma() + current.get_value(). So we SYNTHESIZE
exactly that object (planar circular filaments on a major-radius ring, like a set of TF-style coils),
save it with simsopt.save (the same serialization Gil used), run the real converter, and parse the
MAKEGRID output back -- asserting geometry survives, the loop closes with I=0, the group labels are
present, and --scale multiplies coordinates as documented.

Run:  python -m pytest test_gil_to_makegrid.py -q      (needs simsopt; local python3 has 1.10.6)
"""
import os
import numpy as np
import pytest

simsopt = pytest.importorskip("simsopt")
from simsopt import save                                  # noqa: E402
from simsopt.geo import CurveXYZFourier                   # noqa: E402
from simsopt.field import Current, Coil, BiotSavart       # noqa: E402

from gil_to_makegrid import gil_to_makegrid, coil_set_major_radius, resolve_scale   # noqa: E402


def _ring_coil(R, phi0, npts_dofs=1):
    """A planar circular filament of radius R centred on the origin ring at toroidal angle phi0.
    Built as a CurveXYZFourier so gamma() traces a real circle (matches how Gil's curves serialize)."""
    c = CurveXYZFourier(64, npts_dofs)
    dofs = dict(zip(c.local_dof_names, c.x))
    # circle of radius R in a plane containing z_hat and the radial direction at phi0:
    #   x = (Rmaj + R cos t) cos phi0 ... but for a pure test we just want a closed 3-D loop.
    #   Use a tilted planar circle: center at (Rmaj cos phi0, Rmaj sin phi0, 0).
    Rmaj = 3.0
    dofs["xc(0)"] = Rmaj * np.cos(phi0); dofs["yc(0)"] = Rmaj * np.sin(phi0); dofs["zc(0)"] = 0.0
    # radial+vertical circle: e_r = (cos phi0, sin phi0, 0), e_z = (0,0,1)
    dofs["xc(1)"] = R * np.cos(phi0); dofs["yc(1)"] = R * np.sin(phi0); dofs["zs(1)"] = R
    c.x = np.array([dofs[n] for n in c.local_dof_names])
    return c


def _synth_biotsavart(ncoils=5, R=1.0, I=1.2e6):
    coils = []
    for k in range(ncoils):
        phi0 = 2 * np.pi * k / ncoils
        coils.append(Coil(_ring_coil(R, phi0), Current(I * (1 + 0.01 * k))))
    return BiotSavart(coils)


def _parse_makegrid(path):
    """Return (periods, list-of-coils) where each coil is (Nx4 array of x y z I, group_label)."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    assert lines[0].startswith("periods")
    periods = int(lines[0].split()[1])
    assert lines[1].strip() == "begin filament"
    assert lines[2].strip() == "mirror NIL"
    assert lines[-1].strip() == "end"
    coils, cur = [], []
    label = None
    for ln in lines[3:-1]:
        tok = ln.split()
        xyzI = [float(t) for t in tok[:4]]
        cur.append(xyzI)
        if len(tok) > 4:                                   # loop-closing row carries the group label
            label = " ".join(tok[4:])
            coils.append((np.array(cur), label))
            cur = []
    return periods, coils


def test_roundtrip_geometry_and_format(tmp_path):
    bs = _synth_biotsavart(ncoils=5, R=1.0)
    jpath = str(tmp_path / "biot_savart_optimized_auglag_synth.json")
    save(bs, jpath)

    out = str(tmp_path / "coils_synth")
    n, s = gil_to_makegrid(jpath, out, nfp=1, scale=1.0)
    assert n == 5 and s == pytest.approx(1.0)

    periods, coils = _parse_makegrid(out)
    assert periods == 1
    assert len(coils) == 5
    for k, (arr, label) in enumerate(coils):
        # last row of each coil closes the loop: I=0 and matches the first vertex
        assert arr[-1, 3] == 0.0
        np.testing.assert_allclose(arr[-1, :3], arr[0, :3], atol=1e-9)
        assert label == f"{k+1} Coil{k+1}"
        # geometry: every vertex sits at |point - ringcenter| == R (planar circle), radius ~1.0 m
        ctr = arr[:-1, :3].mean(axis=0)
        radii = np.linalg.norm(arr[:-1, :3] - ctr, axis=1)
        assert radii.max() == pytest.approx(1.0, rel=0.02)


def test_scale_multiplies_coordinates(tmp_path):
    bs = _synth_biotsavart(ncoils=3, R=1.0)
    jpath = str(tmp_path / "biot_savart_optimized_auglag_synth.json")
    save(bs, jpath)

    o1 = str(tmp_path / "coils_1x"); o5 = str(tmp_path / "coils_5x")
    gil_to_makegrid(jpath, o1, nfp=1, scale=1.0)
    gil_to_makegrid(jpath, o5, nfp=1, scale=5.0)
    _, c1 = _parse_makegrid(o1)
    _, c5 = _parse_makegrid(o5)
    for (a1, _), (a5, _) in zip(c1, c5):
        np.testing.assert_allclose(a5[:, :3], a1[:, :3] * 5.0, rtol=1e-10)   # coords scale
        np.testing.assert_allclose(a5[:, 3], a1[:, 3], rtol=1e-10)           # current unchanged


def test_nfp_written_to_periods(tmp_path):
    bs = _synth_biotsavart(ncoils=2, R=1.0)
    jpath = str(tmp_path / "bs.json"); save(bs, jpath)
    out = str(tmp_path / "coils_nfp4")
    gil_to_makegrid(jpath, out, nfp=4, scale=1.0)
    periods, _ = _parse_makegrid(out)
    assert periods == 4


def test_major_radius_measured_from_centroids():
    # _ring_coil centres each filament on a ring of major radius Rmaj=3.0 (see helper)
    bs = _synth_biotsavart(ncoils=6, R=1.0)
    assert coil_set_major_radius(bs.coils) == pytest.approx(3.0, rel=1e-6)


def test_target_major_radius_autoscales(tmp_path):
    # measured R0 = 3.0 m; ask for ARIES-CS-like R0 = 7.75 m -> multiplier 7.75/3.0
    bs = _synth_biotsavart(ncoils=4, R=1.0)
    jpath = str(tmp_path / "bs.json"); save(bs, jpath)
    out = str(tmp_path / "coils_aries")
    n, s = gil_to_makegrid(jpath, out, target_major_radius=7.75)
    assert s == pytest.approx(7.75 / 3.0, rel=1e-6)
    # and the written geometry indeed sits at the requested major radius
    _, coils = _parse_makegrid(out)
    r0_out = np.mean([np.hypot(*a[:-1, :2].mean(axis=0)) for a, _ in coils])
    assert r0_out == pytest.approx(7.75, rel=1e-6)


def test_scale_and_target_are_mutually_exclusive(tmp_path):
    bs = _synth_biotsavart(ncoils=2, R=1.0)
    with pytest.raises(ValueError):
        resolve_scale(bs.coils, scale=2.0, target_major_radius=7.0)
