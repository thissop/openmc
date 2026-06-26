#!/usr/bin/env python
"""Tier 8 / Part B -- Helios-class QA-stellarator radial-build test rig.

PUBLIC placeholder values only (no proprietary Helios/Thea data); every
device-specific slot is tagged INJECT(helios). Public anchors: the Helios design
overview arXiv:2512.08027 (2 field periods, aspect ~4.5, R0~8 m, B0~6 T,
plasma-coil standoff >=1.2 m) and the ARIES-CS radial build (El-Guebaly,
UWFDM-1336). See DATA_NEEDED.md and DEFERRED.md.

IMPORTANT honesty note: the GEOMETRY here is an *axisymmetric* layered box-torus
-- a standoff / radial-build test rig. The 3-D / quasi-axisymmetric physics enters
only through the SOURCE: the spin-polarized birth direction follows the 3-D field
map B-hat(x) (2 field periods, from python/make_standin_field.py), and birth
positions are flux-surface-like ((1-rho^2)-weighted circular surfaces). True 3-D
plasma-boundary shaping (D-shape, field-period geometry, DAGMC) is DEFERRED; do
NOT read field-period/QA *geometry* effects from this model. Units: cm.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
DATADIR = REPO / "spf_prototype" / "data"

import openmc  # noqa: E402
import build_and_run as br  # noqa: E402  (compiled .so + XS path + build_so)
import reactor_model as rm  # noqa: E402  (reuse W/steel/Be/FLiBe materials)

openmc.config["cross_sections"] = str(br.XS)

# --- Helios-class geometry (cm). INJECT(helios): real equilibrium + radial build ---
R0, AMINOR = 800.0, 180.0          # Helios ~8 m / 1.8 m (arXiv:2512.08027)
SOL = 5.0                           # scrape-off gap (ARIES-CS)
U = R0 - AMINOR - SOL              # cavity inboard radius
W = R0 + AMINOR + SOL              # cavity outboard radius
ZC = AMINOR + SOL                  # cavity half-height
# Radial build, cavity outward (cm). ~117 cm => the >=1.2 m standoff.
# INJECT(helios): real per-layer thicknesses, coolant/structure fractions, WC+B4C shield.
LAYERS = [("W", 0.1), ("steel", 1.9), ("Be", 5.0),
          ("FLiBe", 50.0), ("shield", 40.0), ("coil", 15.0)]


def make_materials(li6_enrich=65.0, density_scale=1.0):
    """Reuse the validated W/steel/Be/FLiBe stack; add a WC shield and a coil
    block. INJECT(helios): Li-6 enrichment (Helios ~65%), real shield stack,
    HTS (ReBCO) winding-pack composition for the coil."""
    m = rm.make_materials(li6_enrich=li6_enrich, density_scale=density_scale)
    wc = openmc.Material(name="WC_shield")     # INJECT(helios): WC + B4C + steel + ...
    wc.add_element("W", 1.0); wc.add_element("C", 1.0)
    wc.set_density("g/cm3", 15.6 * density_scale)
    m["shield"] = wc
    coil = openmc.Material(name="coil_block")  # INJECT(helios): ReBCO/Cu/steel winding pack
    coil.add_element("Fe", 0.70, percent_type="wo")
    coil.add_element("Cr", 0.18, percent_type="wo")
    coil.add_element("Ni", 0.12, percent_type="wo")
    coil.set_density("g/cm3", 8.0 * density_scale)
    m["coil"] = coil
    return m


def _box(off, vac=False):
    bc = "vacuum" if vac else "transmission"
    rin = openmc.ZCylinder(r=U - off, boundary_type=bc)
    rout = openmc.ZCylinder(r=W + off, boundary_type=bc)
    zt = openmc.ZPlane(z0=ZC + off, boundary_type=bc)
    zb = openmc.ZPlane(z0=-(ZC + off), boundary_type=bc)
    return +rin & -rout & +zb & -zt


def build_model(abc, fieldmap_stem, li6_enrich=65.0, density_scale=1.0,
                particles=400_000, batches=10, nr=40, nz=40):
    """abc = (a,b,c) polarization; fieldmap_stem = path stem to a B-hat map.
    Source: shape=plasma ((1-rho^2) circular flux surfaces) + bmode=fieldmap."""
    openmc.reset_auto_ids()
    mats = make_materials(li6_enrich, density_scale)
    names = [n for n, _ in LAYERS]; thick = [t for _, t in LAYERS]
    offs = np.cumsum([0.0] + thick)
    boxes = [_box(o, vac=(i == len(offs) - 1)) for i, o in enumerate(offs)]
    cavity = openmc.Cell(name="cavity", region=boxes[0], fill=None)
    cells = {}; layer_cells = []
    for i, name in enumerate(names):
        c = openmc.Cell(name=name, region=boxes[i + 1] & ~boxes[i], fill=mats[name])
        cells[name] = c; layer_cells.append(c)
    geom = openmc.Geometry([cavity] + layer_cells)

    a, b, c = abc
    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.particles = int(particles); s.batches = int(batches); s.inactive = 0
    s.source = openmc.CompiledSource(
        str(br.SO),
        parameters=(f"a={a},b={b},c={c},bmode=fieldmap,fieldmap={fieldmap_stem},"
                    f"shape=plasma,R0={R0},aminor={AMINOR}"))

    tallies = []
    # 1) first-wall incident current, poloidal (first-interaction / NWL map)
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(U, W, nr + 1), phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-ZC, ZC, nz + 1))
    tcur = openmc.Tally(name="wall_current")
    tcur.filters = [openmc.MeshSurfaceFilter(mesh)]; tcur.scores = ["current"]
    tallies.append(tcur)
    # 2) per-cell heating / damage / tritium production
    tcell = openmc.Tally(name="cell_response")
    tcell.filters = [openmc.CellFilter(layer_cells)]
    tcell.scores = ["heating", "damage-energy", "H3-production"]
    tallies.append(tcell)
    # 3) coil fast flux (>0.1 MeV) -- magnet-lifetime proxy
    tcf = openmc.Tally(name="coil_fast")
    tcf.filters = [openmc.CellFilter([cells["coil"]]),
                   openmc.EnergyFilter([0.1e6, 20.0e6])]
    tcf.scores = ["flux"]
    tallies.append(tcf)
    # 4) spatial (poloidal) TBR: z-binned H3 in the FLiBe band
    fl_in = U - offs[names.index("FLiBe")]
    fl_out = W + offs[names.index("FLiBe") + 1]
    smesh = openmc.CylindricalMesh(
        r_grid=np.array([fl_in, fl_out]), phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-(ZC + 50.0), (ZC + 50.0), 11))
    ttbr = openmc.Tally(name="tbr_poloidal")
    ttbr.filters = [openmc.MeshFilter(smesh)]
    ttbr.scores = ["H3-production"]
    tallies.append(ttbr)

    model = openmc.Model(geometry=geom, settings=s,
                         materials=openmc.Materials(list(mats.values())),
                         tallies=openmc.Tallies(tallies))
    return model, cells, (nr, nz)


if __name__ == "__main__":
    # smoke: build + a tiny run, print TBR + coil-fast relative error
    br.build_so()
    if not (DATADIR / "standin_qa.meta").exists():
        import make_standin_field; make_standin_field.main()
    model, cells, _ = build_model((1 / 3, 1 / 3, 1 / 3),
                                  str(DATADIR / "standin_qa"),
                                  particles=40_000, batches=5)
    sp = model.run(cwd="/tmp/spf_t8_3d_smoke", output=False)
    tbr, tbr_sd = rm.tbr(sp, cells["FLiBe"].id)
    print(f"smoke TBR = {tbr:.3f} +/- {tbr_sd:.3f}")
