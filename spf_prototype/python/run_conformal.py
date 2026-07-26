#!/usr/bin/env python
"""STEP (x86/Ginsburg, conda): run OpenMC+DAGMC on the CONFORMAL stellarator wall
with the spin-polarized compiled source driven by the real-equilibrium field map.

This is the deliverable that answers the science question: do SPF directional-
steering benefits SURVIVE non-axisymmetric (conformal) wall smearing, and how much
does scattering change the answer from the free-streaming limit? It is built to run
UNCHANGED on Ginsburg; only the conda env + OPENMC_CROSS_SECTIONS resolve at run
time. Requires OpenMC built WITH DAGMC (environment.yml); imports are guarded so
the aarch64 sandbox gives a clear message instead of a traceback.

Pipeline (each step has its own script; see RUN_ON_GINSBURG.md):
  desc_to_fieldmap.py  -> field map + LCFS surface grid   (isolated desc env)
  stellarator_geometry.py -> per-layer conformal STLs     (any env, numpy)
  build_dagmc.py       -> stellarator.h5m                  (conda, DAGMC)
  run_conformal.py     -> transport + metrics             (conda, OpenMC+DAGMC)

Configs: unpolarized / perpendicular(A) / parallel(B/C), rate held fixed, per
source neutron (Bae 2025 convention). Anchors (3D has no analytic ground truth in
the collided case): (A) free-streaming (near-void) reproduces the analytic
free-streaming NWL of the companion paper, per patch per mode -- the joint
analytic<->MC validation; (B) sampled-direction moments about the local field;
then scattering ON quantifies the departure and the surviving steering.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
DATADIR = REPO / "spf_prototype" / "data"

import build_and_run as br  # noqa: E402  (compiled .so + XS path)
import reactor_model as rm  # noqa: E402  (W/steel/Be/FLiBe materials)
import stellarator_model as sm  # noqa: E402  (shield/coil materials)

# Bae naming -> Schwartz collision modes: unpol=iso, perpendicular=A, parallel=B/C
MODES = {"unpolarized": (1 / 3, 1 / 3, 1 / 3),
         "perpendicular": (1.0, 0.0, 0.0),
         "parallel": (0.0, 1.0, 0.0)}


def scale_fieldmap(stem, scale, out_stem):
    """Rescale a field map's GRID BOUNDS by `scale` (B-hat DIRECTION is invariant
    under uniform geometric scaling, so the .bin is copied verbatim and only the
    .meta R/Z bounds change). Lets one DESC map drive any device size."""
    src_meta = (DATADIR / f"{stem}.meta").read_text().splitlines()
    out = []
    for line in src_meta:
        k = line.split()
        if k and k[0] in ("R_min", "R_max", "Z_min"):
            out.append(f"{k[0]} {float(k[1]) * scale!r}")
        elif k and k[0] == "Z_max":
            out.append(f"{k[0]} {float(k[1]) * scale!r}")
        else:
            out.append(line)
    (DATADIR / f"{out_stem}.meta").write_text("\n".join(out) + "\n")
    # copy the .bin unchanged
    (DATADIR / f"{out_stem}.bin").write_bytes((DATADIR / f"{stem}.bin").read_bytes())
    return DATADIR / out_stem


def make_geometry(h5m_path):
    """Load the DAGMC conformal wall and assign materials by group tag.

    DAGMC matches each volume's 'mat:NAME' group to an openmc.Material whose .name
    == NAME (NOT the python dict key). build_dagmc.py tags volumes with the LAYER
    names (W/steel/Be/FLiBe/shield/coil), so we force each Material.name to equal
    its layer key here (the verifier's F2/D3 fix; otherwise 5 of 6 mismatch and the
    run aborts at geometry load)."""
    try:
        import openmc
    except ImportError as e:
        raise SystemExit(f"openmc unavailable (run on Ginsburg conda env): {e}")
    mats = sm.make_materials()  # keys W/steel/Be/FLiBe/shield/coil
    for tag, m in mats.items():
        m.name = tag            # Material.name must equal the DAGMC 'mat:TAG'
    # implicit complement (SOL + plasma interior + exterior) -> vacuum, else OpenMC
    # errors on an unassigned complement (verifier D10).
    vac = openmc.Material(name="vacuum"); vac.add_nuclide("H1", 1.0)
    vac.set_density("g/cm3", 1e-12)
    dag = openmc.DAGMCUniverse(h5m_path)
    dag.material_names  # noqa: B018 (touch to surface a clear error if tags missing)
    bounded = dag.bounded_universe()  # adds a vacuum-boundary bounding cell
    geom = openmc.Geometry(bounded)
    return geom, openmc.Materials(list(mats.values()) + [vac]), mats


def build_model(abc, h5m_path, fieldmap_stem, R0_cm, a_cm,
                particles=200_000, batches=10, density_scale=1.0):
    import openmc
    openmc.reset_auto_ids()  # stable material/cell ids per build (matches the
    # other models; metrics read ids per-config so this is belt-and-suspenders)
    geom, materials, mats = make_geometry(h5m_path)
    if density_scale != 1.0:  # free-streaming anchor: near-void materials
        for m in materials:
            m.set_density("g/cm3", m.density * density_scale)

    a, b, c = abc
    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.particles = int(particles); s.batches = int(batches); s.inactive = 0
    # Lost-particle tolerance. OpenMC defaults to max_lost_particles=10, so a
    # conformal shell with a few facet-tangle traps on the concave inboard side
    # hard-ABORTS (RuntimeError) even when the lost fraction is negligible -- the
    # sweep's own lost_particle_max=50 warn-threshold never runs because the abort
    # fires first. Env override lets the leaky-recovery pass raise the cap so the
    # run COMPLETES and the sweep can classify honestly on the recorded lost count.
    import os as _os
    s.max_lost_particles = int(_os.environ.get("SPF_MAX_LOST", 10))
    s.rel_max_lost_particles = float(_os.environ.get("SPF_REL_MAX_LOST", 1e-6))
    # Cap per-particle restart-file dumps so a leaky build cannot flood the run dir
    # with thousands of particle_*.h5 files while we probe recoverability.
    try:
        s.max_write_lost_particles = int(_os.environ.get("SPF_MAX_WRITE_LOST", 5))
    except Exception:
        pass
    # source: births on (1-rho^2) flux surfaces (circular MVP, sized inside the FW),
    # B-hat from the real-equilibrium field map. INJECT(helios): a true LCFS
    # flux-surface source is the upgrade (DEFERRED.md).
    a_src = 0.85 * a_cm
    s.source = openmc.CompiledSource(
        str(br.SO),
        parameters=(f"a={a},b={b},c={c},bmode=fieldmap,"
                    f"fieldmap={DATADIR / fieldmap_stem},"
                    f"shape=plasma,R0={R0_cm},aminor={a_src}"))

    tallies = []
    # phi-RESOLVED first-wall current (the 3D signature; toroidal x poloidal mesh)
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0.3 * R0_cm, R0_cm + 1.5 * a_cm, 30),
        phi_grid=np.linspace(0.0, 2 * np.pi, 33),   # 32 toroidal bins -> phi-resolved
        z_grid=np.linspace(-1.5 * a_cm, 1.5 * a_cm, 30))
    tcur = openmc.Tally(name="wall_current_phi")
    tcur.filters = [openmc.MeshSurfaceFilter(mesh)]; tcur.scores = ["current"]
    tallies.append(tcur)
    # per-material response: heating / damage / tritium
    tcell = openmc.Tally(name="cell_response")
    tcell.filters = [openmc.MaterialFilter([mats["W"], mats["steel"], mats["Be"],
                                            mats["FLiBe"], mats["shield"], mats["coil"]])]
    tcell.scores = ["heating", "damage-energy", "H3-production"]
    tallies.append(tcell)
    # coil fast flux (>0.1 MeV) -- magnet-lifetime proxy (deep; the eta_coil observable)
    tcf = openmc.Tally(name="coil_fast")
    tcf.filters = [openmc.MaterialFilter([mats["coil"]]),
                   openmc.EnergyFilter([0.1e6, 20.0e6])]
    tcf.scores = ["flux"]
    tallies.append(tcf)
    # NEAR-SOURCE directional observable (G4): poloidally-resolved first-wall-region
    # fast flux on a cylindrical mesh, reduced to the inboard/outboard-midplane steering
    # contrast in run_ginsburg._wall_directional. Close-in + large solid angle => far
    # lower variance than coil_fast, so this is the low-noise eta_source (coil_fast stays
    # as the deep eta_coil). See docs/EXPERIMENTAL_DESIGN.md.
    # radial band SYMMETRIC about R0 so the inboard (r<R0) and outboard (r>R0) cells are
    # geometrically fair and bracket the wall; the OLD [0.3R0, R0+1.5a] span made the
    # inboard half ~3x wider, sweeping up plasma-hole flux -> a spurious +0.59 inboard
    # contrast baseline (unphysical vs Schwartz's slight OUTBOARD isotropic bias). G4 fix.
    r_half = 1.5 * a_cm
    dmesh = openmc.CylindricalMesh(
        r_grid=np.linspace(max(0.05 * R0_cm, R0_cm - r_half), R0_cm + r_half, 24),
        phi_grid=np.linspace(0.0, 2 * np.pi, 17),
        z_grid=np.linspace(-1.5 * a_cm, 1.5 * a_cm, 24))
    tdir = openmc.Tally(name="firstwall_dir")
    tdir.filters = [openmc.MeshFilter(dmesh), openmc.EnergyFilter([0.1e6, 20.0e6])]
    tdir.scores = ["flux"]
    tallies.append(tdir)

    return openmc.Model(geometry=geom, settings=s, materials=materials,
                        tallies=openmc.Tallies(tallies)), mats


def main():
    h5m = sys.argv[1] if len(sys.argv) > 1 else "stellarator.h5m"
    stem = sys.argv[2] if len(sys.argv) > 2 else "equil_precise_qa"
    scale = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
    # native precise_QA ~ R0 103cm, a 17cm -> scaled device
    R0_cm, a_cm = 103.0 * scale, 17.2 * scale
    fmap = f"{stem}_s{int(scale)}"
    scale_fieldmap(stem, scale, fmap)
    br.build_so()

    res = {}
    for stream, dscale in [("free", 1e-4), ("scatter", 1.0)]:
        for mode, abc in MODES.items():
            model, mats = build_model(abc, h5m, fmap, R0_cm, a_cm, density_scale=dscale)
            sp = model.run(cwd=f"/tmp/spf_conformal_{stream}_{mode}", output=False)
            res[(stream, mode)] = sp
            print(f"ran {stream}/{mode}")

    # metrics + RESULTS_tier8_conformal.md are assembled here on Ginsburg
    # (inboard/outboard steering per mode, free vs scatter; coil fast flux; TBR;
    #  eta = retained steering; joint cross-check vs the analytic free-streaming NWL).
    print("transport complete; assemble metrics into RESULTS_tier8_conformal.md")
    return res


if __name__ == "__main__":
    main()
