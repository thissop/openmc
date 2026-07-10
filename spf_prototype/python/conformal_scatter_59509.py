#!/usr/bin/env python
"""Volumetric flux tally on the conformal blanket, QUASR 59509. Native StellaratorSource
(spf_fluxmap_v1, real sqrt(g) births + local Bhat -- same source as the Sec 9 free-streaming
maps) transported WITH SCATTERING through the multi-layer DAGMC blanket
(W/steel/Be/FLiBe/shield/coil), for unpolarized / perpendicular(A) / parallel(B-C). Saves the
full volumetric fast-flux mesh (poloidal map) and the per-layer heating. Device scaled x10 so
the ~110 cm ARC blanket fits (Bhat is scale-invariant; the fluxmap geometry is rescaled)."""
import sys
from pathlib import Path
import numpy as np

PP = Path("/ginsburg/astro/users/tjk2147/spf_work/spf_pp")
DATA = PP / "data"
sys.path.insert(0, str(PP / "python"))
import stellarator_geometry as sg; sg.DATADIR = DATA
import stellarator_model as sm
from conformal_wallmap import write_fluxmap_bin
import openmc

SCALE = 10.0
CM = 100.0
H5M = str(DATA / "quasr59509_blanket.h5m")
MODES = {"unpolarized": (1/3, 1/3, 1/3), "perpendicular": (1.0, 0.0, 0.0),
         "parallel": (0.0, 1.0, 0.0)}
LAYERS = ["W", "steel", "Be", "FLiBe", "shield", "coil"]


def scaled_fluxmap():
    f = dict(np.load(DATA / "quasr59509_vmec_fluxmap.npz"))
    stem = str(DATA / "quasr59509_cm_fluxmap_s10")
    write_fluxmap_bin(stem, f["rho"], f["theta"], f["zeta"], f["sqrtg"],
                      f["R"] * CM * SCALE, f["Z"] * CM * SCALE, f["phi"],
                      f["BR"], f["Bphi"], f["BZ"], int(f["nfp"]), float(f["sign_sqrtg"]))
    Rl = f["R"][-1] * CM * SCALE
    return stem, 0.5 * (Rl.max() + Rl.min()), 0.5 * (Rl.max() - Rl.min())


def build(abc, fluxstem, R0, a, particles=300_000):
    openmc.reset_auto_ids()
    mats = sm.make_materials()
    for tag, m in mats.items():
        m.name = tag
    vac = openmc.Material(name="vacuum"); vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-12)
    geom = openmc.Geometry(openmc.DAGMCUniverse(H5M).bounded_universe())
    src = openmc.StellaratorSource(fluxmap=fluxstem, polarization=abc, field_model="fluxmap",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.particles = int(particles)
    s.batches = 10; s.inactive = 0; s.photon_transport = False; s.source = [src]
    rh = 1.5 * a
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(max(0.05 * R0, R0 - rh), R0 + rh, 40),
        phi_grid=np.linspace(0.0, 2 * np.pi, 33),
        z_grid=np.linspace(-rh, rh, 40))
    tflux = openmc.Tally(name="flux_mesh")
    tflux.filters = [openmc.MeshFilter(mesh), openmc.EnergyFilter([0.1e6, 20.0e6])]
    tflux.scores = ["flux"]
    theat = openmc.Tally(name="heat")
    theat.filters = [openmc.MaterialFilter([mats[k] for k in LAYERS])]
    theat.scores = ["heating"]
    model = openmc.Model(geom, openmc.Materials(list(mats.values()) + [vac]),
                         s, openmc.Tallies([tflux, theat]))
    return model, mesh


def main():
    fluxstem, R0, a = scaled_fluxmap()
    print(f"[cfg] 59509 x{SCALE:.0f}: R0={R0:.0f}cm a={a:.0f}cm", flush=True)
    out = {"R0": R0, "a": a, "layers": np.array(LAYERS)}
    for mode, abc in MODES.items():
        model, mesh = build(abc, fluxstem, R0, a)
        sp = model.run(cwd=f"/tmp/cs59509_scatter_{mode}", output=False)
        with openmc.StatePoint(sp) as st:
            flux = st.get_tally(name="flux_mesh").mean.ravel()
            heat = st.get_tally(name="heat").mean.ravel()
        nr, nph, nz = len(mesh.r_grid) - 1, len(mesh.phi_grid) - 1, len(mesh.z_grid) - 1
        out[f"flux_{mode}"] = flux.reshape(nr, nph, nz)
        out[f"heat_{mode}"] = heat
        if "r_grid" not in out:
            out["r_grid"] = np.asarray(mesh.r_grid); out["z_grid"] = np.asarray(mesh.z_grid)
            out["phi_grid"] = np.asarray(mesh.phi_grid)
        print(f"[scatter/{mode}] flux sum={flux.sum():.3e}  heat/src={heat.sum():.3e} eV"
              f"  (FLiBe {100*heat[LAYERS.index('FLiBe')]/heat.sum():.0f}% of heating)", flush=True)
    np.savez(DATA / "quasr59509_conformal_flux.npz", **out)
    print("SAVED quasr59509_conformal_flux.npz", flush=True)


if __name__ == "__main__":
    main()
