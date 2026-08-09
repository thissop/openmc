"""Local test of batch_gil_to_makegrid against a synthetic Gil-shaped archive (no 481 MB download).

Builds a nested directory tree mimicking the archive's config-in-path scheme, drops real
simsopt BiotSavart JSONs into it plus one corrupt file, and checks: every good set converts,
the manifest captures symmetry/R0/nfp/config tokens, target-major-radius rescales all sets,
and a single bad file is isolated (status=ERROR) without killing the batch.
"""
import csv
import os

import numpy as np
import pytest

pytest.importorskip("simsopt")
from simsopt import save                                  # noqa: E402
from simsopt.geo import CurveXYZFourier                   # noqa: E402
from simsopt.field import Current, Coil, BiotSavart       # noqa: E402

from batch_gil_to_makegrid import batch_convert, parse_config_from_path   # noqa: E402


def _ring_coil(R, phi0, Rmaj):
    c = CurveXYZFourier(64, 1)
    d = dict(zip(c.local_dof_names, c.x))
    d["xc(0)"] = Rmaj * np.cos(phi0); d["yc(0)"] = Rmaj * np.sin(phi0)
    d["xc(1)"] = R * np.cos(phi0); d["yc(1)"] = R * np.sin(phi0); d["zs(1)"] = R
    c.x = np.array([d[n] for n in c.local_dof_names])
    return c


def _save_set(path, ncoils, Rmaj):
    coils = [Coil(_ring_coil(1.0, 2 * np.pi * k / ncoils, Rmaj), Current(1e6)) for k in range(ncoils)]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    save(BiotSavart(coils), path)


def test_parse_config_from_path():
    cfg = parse_config_from_path("QA/ncoils_3/cs_1p5_cc_0p2/biot_savart_optimized_auglag_x.json")
    assert cfg["symmetry"] == "QA"
    assert cfg["cs"] == "1.5" and cfg["cc"] == "0.2"
    # QI + unknown token preserved in raw_tokens
    cfg2 = parse_config_from_path("QI/weirdknob_9/biot_savart_optimized_auglag_y.json")
    assert cfg2["symmetry"] == "QI"
    assert "weirdknob=9" in cfg2["raw_tokens"]


def test_batch_converts_all_and_isolates_bad(tmp_path):
    root = tmp_path / "gil_coils"
    _save_set(str(root / "QA" / "ncoils_3" / "cs_1p5" / "biot_savart_optimized_auglag_a.json"), 3, 3.0)
    _save_set(str(root / "QH" / "ncoils_4" / "cs_2p0" / "biot_savart_optimized_auglag_b.json"), 4, 5.0)
    # a corrupt "set" that load() will choke on
    bad = root / "QA" / "ncoils_5" / "biot_savart_optimized_auglag_bad.json"
    os.makedirs(bad.parent, exist_ok=True); bad.write_text("{not valid simsopt json")

    out = tmp_path / "makegrid"
    rows, manifest = batch_convert(str(root), str(out), nfp_map={"QA": 2, "QH": 4})

    assert len(rows) == 3
    ok = [r for r in rows if r["status"] == "ok"]
    bad_rows = [r for r in rows if r["status"] != "ok"]
    assert len(ok) == 2 and len(bad_rows) == 1
    assert bad_rows[0]["status"].startswith("ERROR")       # isolated, not fatal

    by_sym = {r["symmetry"]: r for r in ok}
    assert by_sym["QA"]["ncoils"] == 3 and by_sym["QA"]["nfp"] == 2
    assert by_sym["QH"]["ncoils"] == 4 and by_sym["QH"]["nfp"] == 4
    assert by_sym["QA"]["measured_R0_m"] == pytest.approx(3.0, rel=1e-6)
    assert by_sym["QH"]["measured_R0_m"] == pytest.approx(5.0, rel=1e-6)

    # manifest on disk has a header + one row per set, and the output coils files exist
    assert os.path.exists(manifest)
    with open(manifest) as f:
        man = list(csv.DictReader(f))
    assert len(man) == 3
    for r in ok:
        assert os.path.exists(out / r["out_coils"])


def test_target_major_radius_rescales_every_set(tmp_path):
    root = tmp_path / "gil_coils"
    _save_set(str(root / "QA" / "biot_savart_optimized_auglag_a.json"), 3, 3.0)
    _save_set(str(root / "QH" / "biot_savart_optimized_auglag_b.json"), 4, 5.0)
    out = tmp_path / "makegrid"
    rows, _ = batch_convert(str(root), str(out), target_major_radius=7.75, nfp_map={"QA": 2, "QH": 4})
    ok = {r["symmetry"]: r for r in rows if r["status"] == "ok"}
    assert ok["QA"]["applied_scale"] == pytest.approx(7.75 / 3.0, rel=1e-6)
    assert ok["QH"]["applied_scale"] == pytest.approx(7.75 / 5.0, rel=1e-6)
