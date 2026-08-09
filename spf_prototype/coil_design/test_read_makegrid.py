"""Round-trip test: gil_to_makegrid writes -> read_makegrid recovers the coil layout.

Uses a deliberately NON-uniform set of toroidal angles so centroid recovery and the phi-sort are
actually exercised (evenly-spaced coils would hide an ordering bug).
"""
import numpy as np
import pytest

pytest.importorskip("simsopt")
from simsopt import save                                  # noqa: E402
from simsopt.geo import CurveXYZFourier                   # noqa: E402
from simsopt.field import Current, Coil, BiotSavart       # noqa: E402

from gil_to_makegrid import gil_to_makegrid               # noqa: E402
from read_makegrid import read_makegrid, coil_centroids   # noqa: E402


def _ring_coil(R, phi0, Rmaj, Zc=0.0):
    c = CurveXYZFourier(64, 1)
    d = dict(zip(c.local_dof_names, c.x))
    d["xc(0)"] = Rmaj * np.cos(phi0); d["yc(0)"] = Rmaj * np.sin(phi0); d["zc(0)"] = Zc
    d["xc(1)"] = R * np.cos(phi0); d["yc(1)"] = R * np.sin(phi0); d["zs(1)"] = R
    c.x = np.array([d[n] for n in c.local_dof_names])
    return c


def _biotsavart(phis, Rmaj=3.0):
    return BiotSavart([Coil(_ring_coil(1.0, p, Rmaj), Current(1e6)) for p in phis])


def test_read_recovers_layout(tmp_path):
    phis = np.array([0.3, 2.9, 1.1, 5.0])                  # unsorted, non-uniform
    bs = _biotsavart(phis, Rmaj=3.0)
    jp = str(tmp_path / "bs.json"); save(bs, jp)
    out = str(tmp_path / "coils_layout")
    n, _ = gil_to_makegrid(jp, out, nfp=4)
    assert n == 4

    periods, coils = read_makegrid(out)
    assert periods == 4 and len(coils) == 4

    lay = coil_centroids(out)
    assert lay["n_coils"] == 4 and lay["periods"] == 4
    # recovered angles equal the input set, sorted ascending
    np.testing.assert_allclose(lay["phi"], np.sort(phis % (2 * np.pi)), atol=1e-6)
    np.testing.assert_allclose(lay["R"], 3.0, rtol=1e-6)   # all coils on the same major radius


def test_read_recovers_z_offsets(tmp_path):
    # coils at distinct Z to confirm centroid Z is real, not zeroed
    c0 = Coil(_ring_coil(1.0, 0.0, 3.0, Zc=0.5), Current(1e6))
    c1 = Coil(_ring_coil(1.0, np.pi, 3.0, Zc=-0.7), Current(1e6))
    jp = str(tmp_path / "bs.json"); save(BiotSavart([c0, c1]), jp)
    out = str(tmp_path / "coils_z")
    gil_to_makegrid(jp, out, nfp=2)
    lay = coil_centroids(out)
    # sorted by phi: coil at phi=0 first (Z=0.5), then phi=pi (Z=-0.7)
    np.testing.assert_allclose(lay["Z"], [0.5, -0.7], atol=1e-6)


def test_bad_file_raises(tmp_path):
    p = tmp_path / "junk.txt"; p.write_text("not a coils file\n")
    with pytest.raises(ValueError):
        read_makegrid(str(p))
