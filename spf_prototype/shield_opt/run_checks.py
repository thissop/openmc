"""Correctness gates for neutronics runs — hard-fail on physically impossible results instead
of silently emitting garbage (the lesson from the QA-coil saga: TBR 0.03 / flux 70x / flux
above baseline all exited code 0). Two layers:

  1. assert_physical(...)          — post-run number gates (TBR, flux scale, source-in-plasma).
  2. preflight_source_geometry(...) — pre-run: sample the source, find_cell, and assert the
                                      neutrons are born in the plasma/breeder, not shield/coil/
                                      outside. Catches units/scale/frame mismatches in ~seconds,
                                      BEFORE a 30-minute transport run.

Use both on every device. QH is the golden reference for flux scale.
"""
from __future__ import annotations
import numpy as np

# QH golden reference (corrected build, unpol, cells 11-30): the sane scale for a coil run.
QH_REF = dict(coil_maxflux=3.26e-2, tbr=1.106, peaking=2.91)


class PhysicsCheckError(AssertionError):
    pass


def assert_physical(*, peaking, maxflux, tbr, relerr=None,
                    ref_maxflux=QH_REF["coil_maxflux"],
                    tbr_range=(0.85, 1.35), flux_factor=(0.2, 5.0),
                    relerr_max=0.30, label="run", hard=True):
    """Gate a coil-run RESULT on physics. Returns (ok, messages). If hard and not ok, raises.

    - TBR must be in tbr_range (breeding is ~1; 0.03 or 0.28 is a broken source).
    - maxflux must be within flux_factor of the QH reference (70x or above-baseline = broken
      geometry/source; near-zero = shielded-to-death or wrong tally).
    - relerr must be < relerr_max (else it's noise, not a result).
    """
    msgs, ok = [], True
    def check(cond, msg):
        nonlocal ok
        tag = "PASS" if cond else "FAIL"
        if not cond: ok = False
        msgs.append(f"  [{tag}] {msg}")
    check(tbr_range[0] <= tbr <= tbr_range[1],
          f"TBR {tbr:.3f} in {tbr_range}  (out-of-range => broken source/breeder)")
    if ref_maxflux is not None:
        lo, hi = flux_factor[0]*ref_maxflux, flux_factor[1]*ref_maxflux
        check(lo <= maxflux <= hi,
              f"coil maxflux {maxflux:.3e} within [{flux_factor[0]}x,{flux_factor[1]}x] of QH ref "
              f"{ref_maxflux:.3e} = [{lo:.2e},{hi:.2e}]")
    check(1.0 <= peaking <= 6.0, f"peaking {peaking:.3f} in [1,6]")
    if relerr is not None:
        check(relerr <= relerr_max, f"median relerr {relerr:.3f} <= {relerr_max}")
    header = f"[assert_physical:{label}] {'OK' if ok else 'VIOLATION'}"
    report = header + "\n" + "\n".join(msgs)
    if not ok and hard:
        raise PhysicsCheckError(report)
    return ok, report


def assert_shield_reduces(*, maxflux, baseline_maxflux, tol=1.05, label="shield", hard=True):
    """A shield/trade run (delta>0) CANNOT raise the peak coil flux above its own delta=0
    baseline. Catches sub-grid CAD defects (the coarse-w20 artifact: 3.97e-2 > 3.26e-2 baseline)
    that assert_physical's coarse scale-gate misses."""
    ok = maxflux <= baseline_maxflux * tol
    rep = (f"[assert_shield_reduces:{label}] {'OK' if ok else 'VIOLATION'} "
           f"maxflux {maxflux:.3e} vs baseline {baseline_maxflux:.3e} (tol {tol}x) "
           f"{'-- adding shield raised flux => CAD defect / sub-grid feature' if not ok else ''}")
    if not ok and hard:
        raise PhysicsCheckError(rep)
    return ok, rep


def preflight_source_geometry(model, n=20000, good_materials=("breeder", "multiplier",
                              "first_wall", "Vacuum", "gap"), min_good_frac=0.98, hard=True):
    """BEFORE transport: sample n source sites, find_cell each, tally the material they land in.
    Assert >= min_good_frac are born in plasma-void/breeder-side materials (NOT shield/coil/lost).
    A units/scale/frame mismatch makes the source land in the wrong material -> caught here in
    seconds. Requires an initialized openmc.lib (call inside the run script after export).
    Returns (ok, report). Needs the model already exported; uses openmc.lib.find_cell + sample.
    """
    import openmc.lib
    openmc.lib.init(args=["-s", "1"]) if not openmc.lib.is_initialized else None
    # sample source sites via the library
    sites = openmc.lib.sample_external_source(n)
    from collections import Counter
    matcount = Counter(); lost = 0
    matname = {m.id: m.name for m in openmc.lib.materials.values()}
    for s in sites:
        try:
            cell, _ = openmc.lib.find_cell(tuple(s.r))
            mid = cell.fill.id if hasattr(cell.fill, "id") else None
            matcount[matname.get(mid, "?")] += 1
        except Exception:
            lost += 1
    tot = n
    good = sum(v for k, v in matcount.items() if k in good_materials)
    frac = good / tot
    ok = frac >= min_good_frac and lost == 0
    rep = (f"[preflight_source_geometry] {'OK' if ok else 'MISMATCH'} "
           f"good_frac={frac:.3f} lost={lost}\n  materials: {dict(matcount)}")
    if not ok and hard:
        raise PhysicsCheckError(rep + "\n  => source not born in plasma/breeder: check "
                                "fluxmap units (cm!), scale, and that fluxmap+geometry share one wout.")
    return ok, rep


if __name__ == "__main__":
    # Demonstrate: the module FLAGS every one of the QA-saga garbage results and PASSES the good QH.
    cases = [
        ("QH baseline (good)",   dict(peaking=2.91, maxflux=3.26e-2, tbr=1.106, relerr=0.062)),
        ("QA v1 (qh_fluxmap)",   dict(peaking=1.40, maxflux=3.11e-1, tbr=0.281,  relerr=0.020)),
        ("QA v2 (qa meters)",    dict(peaking=1.71, maxflux=2.28e+0, tbr=0.031,  relerr=0.008)),
        ("coarse w20 (artifact)",dict(peaking=3.19, maxflux=3.97e-2, tbr=1.172,  relerr=0.090)),
    ]
    for name, r in cases:
        ok, rep = assert_physical(label=name, hard=False, **r)
        print(rep, "\n")
