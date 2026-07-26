#!/usr/bin/env python
"""PRE-FLIGHT GATE (openmc+DAGMC env): cheap watertightness probe on a DAGMC h5m BEFORE any
expensive adjoint/forward transport. Catches the exact failure the conformal sweep hit only at
transport time ("Maximum number of lost particles reached / No intersection found with DAGMC
cell N") -- i.e. a non-watertight build -- and classifies the device GO / dagmc_leaky CLEANLY.

How: build a minimal model (low-density filler in every DAGMC material tag so neutrons stream
far and probe surfaces), an isotropic 14 MeV point source at the geometry center, and run ONE
small fixed-source batch with a low max_lost_particles. If OpenMC exceeds the lost-particle cap
(hard leak) the run aborts -- but THIS SCRIPT is meant to be launched as a SUBPROCESS by
run_device.py, so that abort is contained: we catch it, and either way we write
<h5m_dir>/<stem>_watertight.json with status and lost-particle fraction and exit with a code
run_device interprets (0=GO, 7=leaky, 8=error). Nothing crashes the orchestrator.

Usage: python dagmc_watertight_check.py <dagmc.h5m> [--particles 10000] [--max-lost 20]
                                         [--rel-max-lost 1e-3] [--out STATUS.json]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def build_probe_model(h5m, particles, max_lost, rel_max_lost):
    import numpy as np
    import openmc

    dag = openmc.DAGMCUniverse(str(h5m), auto_geom_ids=True)
    # a low-density filler for EVERY material tag the DAGMC references, so the probe never fails
    # for a missing-material reason (that is a different, also-clean error) and neutrons stream.
    try:
        names = list(dag.material_names)
    except Exception:
        names = []
    mats = []
    for nm in names:
        if nm.lower() in ("void", "vacuum", "graveyard"):
            continue
        m = openmc.Material(name=nm)
        m.add_nuclide("H1", 1.0)
        m.set_density("g/cm3", 1e-4)      # very thin -> long mean free path -> probes surfaces
        mats.append(m)
    materials = openmc.Materials(mats)

    # bound the DAGMC universe in a VACUUM sphere + root cell, exactly as adjoint_importance.py
    # (the DAGMC h5m has no boundary condition of its own; without this OpenMC MPI_ABORTs on
    # "No boundary conditions were applied", which is a setup error, NOT a watertightness leak).
    bb = dag.bounding_box
    ll, ur = np.asarray(bb.lower_left), np.asarray(bb.upper_right)
    ext = float(np.max(np.abs(np.concatenate([ll, ur]))))
    sph = openmc.Sphere(r=ext * 1.5 + 10.0, boundary_type="vacuum")
    root = openmc.Cell(region=-sph, fill=dag)
    geometry = openmc.Geometry([root])
    center = tuple(0.5 * (ll + ur))

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = 1
    settings.particles = int(particles)
    settings.max_lost_particles = int(max_lost)
    settings.rel_max_lost_particles = float(rel_max_lost)
    settings.output = {"summary": False, "tallies": False}
    src = openmc.IndependentSource()
    src.space = openmc.stats.Point(center)
    src.angle = openmc.stats.Isotropic()
    src.energy = openmc.stats.Discrete([14.06e6], [1.0])
    settings.source = src
    return openmc.Model(geometry, materials, settings), center


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5m")
    ap.add_argument("--particles", type=int, default=10000)
    ap.add_argument("--max-lost", type=int, default=20)
    ap.add_argument("--rel-max-lost", type=float, default=1e-3)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    h5m = Path(a.h5m)
    out = Path(a.out) if a.out else h5m.with_name(h5m.stem + "_watertight.json")

    def write(status, **extra):
        rec = dict(stage="dagmc_watertight", h5m=str(h5m), status=status, **extra)
        out.write_text(json.dumps(rec, indent=2))
        print(f"[watertight] STATUS={status} " + " ".join(f"{k}={v}" for k, v in extra.items()),
              flush=True)
        return rec

    if not h5m.exists():
        write("error:dagmc_watertight:missing_h5m")
        return 8

    try:
        import openmc
        model, center = build_probe_model(h5m, a.particles, a.max_lost, a.rel_max_lost)
    except Exception as e:
        write(f"error:dagmc_watertight:setup:{type(e).__name__}", message=str(e)[:200])
        return 8

    try:
        sp = model.run(output=True, threads=None)
        # completed without hitting the lost-particle cap -> watertight enough to transport.
        import openmc
        n_lost = None
        with openmc.StatePoint(sp) as s:
            n_lost = int(getattr(s, "n_lost_particles", 0) or 0)
        frac = (n_lost / a.particles) if a.particles else 0.0
        leaky = frac > a.rel_max_lost
        write("dagmc_leaky" if leaky else "go",
              n_lost=n_lost, particles=a.particles, lost_frac=round(frac, 6),
              rel_max_lost=a.rel_max_lost, center=[round(c, 2) for c in center])
        return 7 if leaky else 0
    except Exception as e:
        # OpenMC aborted (typically: max lost particles reached / no DAGMC intersection).
        # THIS is the leak we are pre-detecting; classify cleanly (contained in this subprocess).
        msg = str(e)[:300]
        is_leak = any(k in msg.lower() for k in
                      ("lost particle", "no intersection", "dagmc", "max_lost"))
        write("dagmc_leaky" if is_leak else f"error:dagmc_watertight:{type(e).__name__}",
              particles=a.particles, message=msg)
        return 7 if is_leak else 8


if __name__ == "__main__":
    rc = main()
    # os._exit bypasses Python/C++ atexit cleanup -- libopenmc's static destructors trigger a
    # BENIGN exit-time double-free (SIGABRT -> 134) AFTER our verdict is already written. Exiting
    # with our own code keeps SLURM afterok honest (GO=0 must not look like a failure).
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(rc)
