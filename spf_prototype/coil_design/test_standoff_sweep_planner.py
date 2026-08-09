"""Tests for the standoff-sweep planner: family selection, cs-spread coverage, edge cases."""
import csv

import pytest

from standoff_sweep_planner import (choose_family, select_spread, plan_sweep, load_manifest)


def _row(sym, ncoils, cs, status="ok", r0="6.0"):
    return dict(symmetry=sym, ncoils=str(ncoils), cs=("" if cs is None else str(cs)),
                status=status, measured_R0_m=r0, out_coils=f"coils_{sym}_{ncoils}_{cs}",
                source_json=f"{sym}/{ncoils}/{cs}.json")


def _write_manifest(tmp_path, rows):
    p = tmp_path / "manifest.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return str(p)


def test_choose_family_picks_widest_cs_span():
    rows = [
        _row("QA", 3, 1.0), _row("QA", 3, 1.2),                    # QA/3 span 0.2
        _row("QH", 4, 1.0), _row("QH", 4, 2.5), _row("QH", 4, 1.7),  # QH/4 span 1.5 (widest)
    ]
    key, members = choose_family(rows)
    assert key == ("QH", "4") and len(members) == 3


def test_choose_family_respects_forced_symmetry():
    rows = [_row("QA", 3, 1.0), _row("QA", 3, 2.0), _row("QH", 4, 1.0), _row("QH", 4, 9.0)]
    key, members = choose_family(rows, symmetry="QA")
    assert key[0] == "QA" and len(members) == 2


def test_usable_filters_bad_and_missing_cs():
    rows = [_row("QA", 3, 1.0), _row("QA", 3, None), _row("QA", 3, 2.0, status="ERROR:x")]
    key, members = choose_family(rows)
    assert len(members) == 1                                        # only the ok + numeric-cs row


def test_select_spread_covers_range():
    members = [dict(_cs=c, symmetry="QH", ncoils="4", out_coils="", measured_R0_m="6",
                    source_json="") for c in [1.0, 1.1, 1.2, 2.0, 3.0, 5.0]]
    chosen = select_spread(members, 3)
    cs = [m["_cs"] for m in chosen]
    assert cs[0] == 1.0 and cs[-1] == 5.0                           # endpoints kept
    assert cs == sorted(cs)                                         # ascending
    assert len(chosen) == 3


def test_select_spread_dedupes_when_points_exceed_members():
    members = [dict(_cs=c, symmetry="QH", ncoils="4", out_coils="", measured_R0_m="6",
                    source_json="") for c in [1.0, 2.0]]
    chosen = select_spread(members, 5)
    assert len(chosen) == 2                                         # no duplicates invented


def test_plan_sweep_end_to_end(tmp_path):
    rows = [_row("QH", 4, c) for c in [1.0, 1.4, 1.9, 2.6, 3.3]] + [_row("QA", 3, 1.0), _row("QA", 3, 1.1)]
    mpath = _write_manifest(tmp_path, rows)
    plan, rationale = plan_sweep(mpath, n_points=3)
    assert "QH" in rationale and len(plan) == 3
    assert plan[0]["cs"] == 1.0 and plan[-1]["cs"] == 3.3
    assert [p["order"] for p in plan] == [0, 1, 2]


def test_plan_sweep_raises_when_no_usable(tmp_path):
    rows = [_row("QA", 3, None), _row("QH", 4, 1.0, status="ERROR:boom")]
    mpath = _write_manifest(tmp_path, rows)
    with pytest.raises(ValueError):
        plan_sweep(mpath)
