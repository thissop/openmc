#!/usr/bin/env python
"""Tier 5/6 reactor model: a layered box-torus with the spin-polarized compiled
source, for SCATTERING (Tier 5: W first wall, heating/damage/spectrum) and a
FLiBe breeder blanket (Tier 6: TBR per polarization mode).

Geometry keeps the validated square cross-section (so the thin-wall limit recovers
the free-streaming Tier-2b result), scaled to ARC-class size (minor radius ~1 m).
Layers lining the cavity, outward (paper-grounded; see BLANKET_PARAMS.md):
    cavity(void) -> W 0.1cm -> Fe-9Cr steel 1cm -> Be 1cm -> FLiBe 40cm -> vacuum

All parameters traceable to Sorbom 2015 (ARC), Segantin 2020 (TBR opt), and
Peterson 2022 (LIBRA). Units: cm.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
import openmc  # noqa: E402
import build_and_run as br  # noqa: E402  (reuses the compiled .so + XS path + build_so)

openmc.config["cross_sections"] = str(br.XS)

# --- geometry (cm): same proportions as Tier 2b (u/R0=0.4, w/R0=1.6, z/R0=0.6) ---
R0, AMINOR = 200.0, 100.0          # plasma (ARC-class minor radius ~1 m)
U, W, ZC = 80.0, 320.0, 120.0      # cavity walls (20 cm SOL gap)
# layer thicknesses (cm), cavity outward
T_W, T_STEEL, T_BE, T_FLIBE = 0.1, 1.0, 1.0, 40.0
LAYERS = ["W", "steel", "Be", "FLiBe"]
THICK = [T_W, T_STEEL, T_BE, T_FLIBE]


def make_materials(li6_enrich=30.0, density_scale=1.0):
    """W / Fe-9Cr steel / Be / FLiBe(Li2BeF4, Li-6 enriched). density_scale<<1
    gives the free-streaming 'thin-wall' anchor."""
    w = openmc.Material(name="tungsten")
    w.add_element("W", 1.0)
    w.set_density("g/cm3", 19.3 * density_scale)

    steel = openmc.Material(name="FeCr_steel")     # Eurofer-like, data-safe (Fe,Cr)
    steel.add_element("Fe", 0.91, percent_type="wo")
    steel.add_element("Cr", 0.09, percent_type="wo")
    steel.set_density("g/cm3", 7.8 * density_scale)

    be = openmc.Material(name="beryllium")
    be.add_element("Be", 1.0)
    be.set_density("g/cm3", 1.85 * density_scale)

    flibe = openmc.Material(name="FLiBe")          # Li2BeF4, atoms Li:Be:F = 2:1:4
    flibe.add_element("Li", 2.0, percent_type="ao",
                      enrichment=li6_enrich, enrichment_target="Li6", enrichment_type="ao")
    flibe.add_element("Be", 1.0, percent_type="ao")
    flibe.add_element("F", 4.0, percent_type="ao")
    flibe.set_density("g/cm3", 1.94 * density_scale)
    return {"W": w, "steel": steel, "Be": be, "FLiBe": flibe}


def _box(off, vac=False):
    """Square-cross-section toroidal box region at offset `off` from the cavity:
    u-off < R < w+off, |z| < ZC+off. Outermost gets vacuum boundaries."""
    bc = "vacuum" if vac else "transmission"
    rin = openmc.ZCylinder(r=U - off, boundary_type=bc)
    rout = openmc.ZCylinder(r=W + off, boundary_type=bc)
    zt = openmc.ZPlane(z0=ZC + off, boundary_type=bc)
    zb = openmc.ZPlane(z0=-(ZC + off), boundary_type=bc)
    return +rin & -rout & +zb & -zt


def build_model(abc, li6_enrich=30.0, density_scale=1.0, particles=100_000,
                batches=10, nr=40, nz=40, spectrum_bins=None):
    openmc.reset_auto_ids()
    mats = make_materials(li6_enrich, density_scale)

    offs = np.cumsum([0.0] + THICK)          # [0, tW, tW+tS, ..., total]
    boxes = [_box(o, vac=(i == len(offs) - 1)) for i, o in enumerate(offs)]
    cavity = openmc.Cell(name="cavity", region=boxes[0], fill=None)
    cells = {}
    layer_cells = []
    for i, name in enumerate(LAYERS):
        c = openmc.Cell(name=name, region=boxes[i + 1] & ~boxes[i], fill=mats[name])
        cells[name] = c
        layer_cells.append(c)
    geom = openmc.Geometry([cavity] + layer_cells)

    a, b, c = abc
    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.particles = int(particles)
    s.batches = int(batches)
    s.inactive = 0
    s.source = openmc.CompiledSource(
        str(br.SO),
        parameters=f"a={a},b={b},c={c},bmode=toroidal,shape=plasma,R0={R0},aminor={AMINOR}")

    tallies = []
    # 1) cavity-boundary surface current (NWL / free-streaming anchor)
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(U, W, nr + 1),
        phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-ZC, ZC, nz + 1))
    tcur = openmc.Tally(name="wall_current")
    tcur.filters = [openmc.MeshSurfaceFilter(mesh)]
    tcur.scores = ["current"]
    tallies.append(tcur)
    # 2) per-cell heating / damage / tritium production
    tcell = openmc.Tally(name="cell_response")
    tcell.filters = [openmc.CellFilter(layer_cells)]
    tcell.scores = ["heating", "damage-energy", "H3-production"]
    tallies.append(tcell)
    # 3) FLiBe flux spectrum (down-scattered tail = scattering signature)
    if spectrum_bins is not None:
        tspec = openmc.Tally(name="flibe_spectrum")
        tspec.filters = [openmc.CellFilter([cells["FLiBe"]]), openmc.EnergyFilter(spectrum_bins)]
        tspec.scores = ["flux"]
        tallies.append(tspec)
    # 4) tritium production split by Li-6 / Li-7 (internal verification of TBR)
    tnuc = openmc.Tally(name="tbr_nuclide")
    tnuc.filters = [openmc.CellFilter([cells["FLiBe"]])]
    tnuc.nuclides = ["Li6", "Li7"]
    tnuc.scores = ["H3-production"]
    tallies.append(tnuc)

    model = openmc.Model(geometry=geom, settings=s,
                         materials=openmc.Materials(list(mats.values())),
                         tallies=openmc.Tallies(tallies))
    return model, cells, (nr, nz)


# --------------------------------------------------------------------------
# Extraction helpers
# --------------------------------------------------------------------------
# OpenMC labels the 'H3-production' score as '(n,Xt)' in the pandas dataframe.
_TRIT_SCORES = ("H3-production", "(n,Xt)", "(n,t)")


def _flatcol(df, name):
    """Column access tolerant of flat or MultiIndex frames."""
    if name in df.columns:
        return df[name]
    return df[(name, "")]


def wall_ratio_inboard_outboard(sp_path, nr, nz):
    """Incident-current intensity outboard/inboard at the midplane (free-stream anchor)."""
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="wall_current").get_pandas_dataframe()
    xi = df[("mesh 1", "x")].to_numpy(); zi = df[("mesh 1", "z")].to_numpy()
    surf = df[("mesh 1", "surf")].to_numpy(); mean = df[("mean", "")].to_numpy()
    midbins = (nz // 2, nz // 2 + 1)  # two central z bins straddling the midplane
    inb = mean[(surf == "x-min out") & (xi == 1) & np.isin(zi, midbins)].sum()
    outb = mean[(surf == "x-max out") & (xi == nr) & np.isin(zi, midbins)].sum()
    dz = 2 * ZC / nz
    return (outb / (2 * np.pi * W * dz)) / (inb / (2 * np.pi * U * dz))


def midplane_intensities(sp_path, nr, nz, nwin=4):
    """Incident-current intensity (current/area) at the inboard and outboard walls,
    averaged over the central 2*nwin z-bins (the midplane band). Mode/iso ratios of
    these = directionality with scattering on. nwin>1 reduces statistical noise."""
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="wall_current").get_pandas_dataframe()
    xi = df[("mesh 1", "x")].to_numpy(); zi = df[("mesh 1", "z")].to_numpy()
    surf = df[("mesh 1", "surf")].to_numpy(); mean = df[("mean", "")].to_numpy()
    midbins = list(range(nz // 2 - nwin + 1, nz // 2 + nwin + 1))  # 2*nwin central bins
    dz = 2 * ZC / nz
    nb = len(midbins)
    inb = mean[(surf == "x-min out") & (xi == 1) & np.isin(zi, midbins)].sum() / (2 * np.pi * U * dz * nb)
    outb = mean[(surf == "x-max out") & (xi == nr) & np.isin(zi, midbins)].sum() / (2 * np.pi * W * dz * nb)
    return inb, outb


def cell_table(sp_path, id_to_name):
    """Per-cell {heating[eV/src], damage[eV/src], tbr[trit/src]} (+ std for tbr)."""
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="cell_response").get_pandas_dataframe()
    cid = _flatcol(df, "cell").to_numpy()
    score = _flatcol(df, "score").to_numpy()
    mean = _flatcol(df, "mean").to_numpy()
    sd = _flatcol(df, "std. dev.").to_numpy()
    out = {}
    for the_id, name in id_to_name.items():
        rows = cid == the_id
        rec = {}
        for sc, m, s in zip(score[rows], mean[rows], sd[rows]):
            if sc == "heating":
                rec["heating"] = float(m)
            elif sc == "damage-energy":
                rec["damage"] = float(m)
            elif sc in _TRIT_SCORES:
                rec["tbr"] = float(m); rec["tbr_sd"] = float(s)
        out[name] = rec
    return out


def tbr(sp_path, flibe_cell_id):
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="cell_response").get_pandas_dataframe()
    cid = _flatcol(df, "cell").to_numpy()
    score = _flatcol(df, "score").to_numpy()
    m = (cid == flibe_cell_id) & np.isin(score, _TRIT_SCORES)
    return float(_flatcol(df, "mean").to_numpy()[m][0]), float(_flatcol(df, "std. dev.").to_numpy()[m][0])


def _selftest():
    br.build_so()
    print("=== self-test: iso, real density, 40 cm FLiBe @30% Li-6 ===")
    model, cells, (nr, nz) = build_model((1/3, 1/3, 1/3), li6_enrich=30.0,
                                         density_scale=1.0, particles=40000, batches=10)
    sp = model.run(cwd="/tmp/spf_t56_real", output=False)
    t, ts = tbr(sp, cells["FLiBe"].id)
    print(f"  TBR (H3-production in FLiBe) = {t:.3f} +/- {ts:.3f}  (expect ~1.1-1.4)")
    # heating split
    with openmc.StatePoint(sp) as spo:
        df = spo.get_tally(name="cell_response").get_pandas_dataframe()
    for name in LAYERS:
        cid = cells[name].id
        h = df[(df[("cell", "id")] == cid) & (df[("score", "")] == "heating")][("mean", "")]
        if len(h):
            print(f"  heating in {name:6s} = {float(h.iloc[0]):.3e} eV/src")
    print(f"  outboard/inboard incident ratio (real) = {wall_ratio_inboard_outboard(sp, nr, nz):.3f}")

    print("=== self-test: iso, THIN walls (x1e-4) -> must recover free-streaming ~1.13 ===")
    model2, cells2, (nr2, nz2) = build_model((1/3, 1/3, 1/3), li6_enrich=30.0,
                                             density_scale=1e-4, particles=40000, batches=10)
    sp2 = model2.run(cwd="/tmp/spf_t56_thin", output=False)
    print(f"  outboard/inboard incident ratio (thin) = {wall_ratio_inboard_outboard(sp2, nr2, nz2):.3f}"
          f"  (free-streaming analytic ~1.146)")


if __name__ == "__main__":
    _selftest()
