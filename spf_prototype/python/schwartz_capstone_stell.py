#!/usr/bin/env python
"""StellaratorSource transport capstone (README build-gate 2): the same square-torus Schwartz NWL
test as schwartz_capstone.py, but driven by the native C++ StellaratorSource on a CIRCULAR-torus
fluxmap with field_model='toroidal'. Confirms the StellaratorSource reproduces the anarrima oracle
(+40.6% inboard / -21.2% outboard) -- symmetric with the native TokamakSource result.

Cluster run env identical to schwartz_capstone.py (PR openmc + built lib via LD_PRELOAD).
"""
import sys
import json
from pathlib import Path
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stellarator_source import circular_torus_fluxmap
from quasr_fluxmap import write_fluxmap_bin

SCALE = 100.0
R0, A = 1.0 * SCALE, 0.5 * SCALE
RIN, ROUT, ZTB = 0.4 * SCALE, 1.6 * SCALE, 0.6 * SCALE
DZ_MID, NZ = 0.05 * SCALE, 40
MODES = {"unpol": (1/3, 1/3, 1/3), "perp": (1.0, 0.0, 0.0), "par": (0.0, 1.0, 0.0)}


def make_fluxmap(stem):
    fm, _ = circular_torus_fluxmap(R0=R0, a=A, nr=16, nt=64, nz=96, nfp=1)
    write_fluxmap_bin(stem, fm["rho"], fm["theta"], fm["zeta"], fm["sqrtg"], fm["R"], fm["Z"],
                      fm["phi"], fm["BR"], fm["Bphi"], fm["BZ"], fm["nfp"], 1.0)


def build(pol, stem):
    void = openmc.Material(name="void"); void.add_nuclide("H1", 1.0)
    void.set_density("g/cm3", 1e-10)
    cin = openmc.ZCylinder(r=RIN, boundary_type="vacuum")
    cout = openmc.ZCylinder(r=ROUT, boundary_type="vacuum")
    zt = openmc.ZPlane(z0=ZTB, boundary_type="vacuum")
    zb = openmc.ZPlane(z0=-ZTB, boundary_type="vacuum")
    cell = openmc.Cell(fill=void, region=+cin & -cout & +zb & -zt)
    src = openmc.StellaratorSource(fluxmap=stem, polarization=pol, field_model="toroidal",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.source = [src]
    s.particles = 200_000; s.batches = 20; s.photon_transport = False
    mesh = openmc.CylindricalMesh(r_grid=[RIN, ROUT], phi_grid=[0.0, 2 * np.pi],
                                  z_grid=list(np.linspace(-ZTB, ZTB, NZ + 1)))
    t = openmc.Tally(name="wall"); t.filters = [openmc.MeshSurfaceFilter(mesh)]
    t.scores = ["current"]
    return openmc.Model(openmc.Geometry([cell]), openmc.Materials([void]), s, openmc.Tallies([t]))


def midplane_persurf(sp):
    with openmc.StatePoint(sp) as s:
        mean = np.asarray(s.get_tally(name="wall").mean).ravel()
    nsurf = mean.size // NZ
    arr = mean.reshape(NZ, nsurf)
    zc = np.linspace(-ZTB, ZTB, NZ + 1); zmid = 0.5 * (zc[:-1] + zc[1:])
    return arr[np.abs(zmid) < DZ_MID].sum(axis=0)


def main():
    stem = str(HERE.parent / "data" / "circular_fluxmap")
    make_fluxmap(stem)
    res = {}
    for name, pol in MODES.items():
        sp = build(pol, stem).run(cwd=f"capstell_{name}")
        res[name] = midplane_persurf(sp).tolist()
    json.dump(res, open("schwartz_capstone_stell.json", "w"), indent=2)
    u, p, q = (np.array(res[k]) for k in ("unpol", "perp", "par"))
    print("=== StellaratorSource capstone (circular fluxmap, toroidal b_hat) ===")
    for si, label in [(0, "inboard "), (2, "outboard")]:
        print(f"{label} (surf {si}): unpol {u[si]:.4g} | perp {100*(p[si]/u[si]-1):+.1f}% "
              f"| par {100*(q[si]/u[si]-1):+.1f}%")
    print("oracle : perp +40.6% inboard / -21.2% outboard ; par the mirror")


if __name__ == "__main__":
    main()
