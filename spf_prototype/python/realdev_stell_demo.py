#!/usr/bin/env python
"""Real-equilibrium demonstration: run the native C++ StellaratorSource in transport on the REAL
precise_QA DESC equilibrium (real sqrt(g) births + real b_hat from the fluxmap), free-streaming, and
tally the 3-D neutron flux map for unpolarized / perpendicular / parallel. Shows the StellaratorSource
working end-to-end on an actual stellarator equilibrium and the polarization redistributing the load
in 3-D. (Not a conformal-wall NWL -- a bounded-cylinder flux map, illustrative.)

Cluster run env = schwartz_capstone.py's (PR openmc + built lib via LD_PRELOAD).
"""
import sys
from pathlib import Path
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from quasr_fluxmap import write_fluxmap_bin

CM = 100.0
MODES = {"unpol": (1/3, 1/3, 1/3), "perp": (1.0, 0.0, 0.0), "par": (0.0, 1.0, 0.0)}


def make_cm_fluxmap(stem):
    """Rescale the meters-unit precise_QA fluxmap to cm (positions); b_hat/sqrt(g) are scale-free."""
    f = dict(np.load(HERE.parent / "data" / "precise_QA_fluxmap.npz"))
    write_fluxmap_bin(stem, f["rho"], f["theta"], f["zeta"], f["sqrtg"],
                      f["R"] * CM, f["Z"] * CM, f["phi"],
                      f["BR"], f["Bphi"], f["BZ"], int(f["nfp"]), float(f["sign_sqrtg"]))
    R = f["R"] * CM; Z = f["Z"] * CM
    return float(R.min()), float(R.max()), float(np.abs(Z).max())


def build(pol, stem, Rlo, Rhi, Zmax):
    void = openmc.Material(name="void"); void.add_nuclide("H1", 1.0)
    void.set_density("g/cm3", 1e-10)
    Rwall = Rhi + 20.0; Zwall = Zmax + 20.0
    outer = openmc.ZCylinder(r=Rwall, boundary_type="vacuum")
    zt = openmc.ZPlane(z0=Zwall, boundary_type="vacuum")
    zb = openmc.ZPlane(z0=-Zwall, boundary_type="vacuum")
    cell = openmc.Cell(fill=void, region=-outer & +zb & -zt)
    src = openmc.StellaratorSource(fluxmap=stem, polarization=pol, field_model="fluxmap",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.source = [src]
    s.particles = 400_000; s.batches = 20; s.photon_transport = False
    mesh = openmc.CylindricalMesh(
        r_grid=list(np.linspace(max(1.0, Rlo - 20), Rwall, 40)),
        phi_grid=list(np.linspace(0, 2 * np.pi, 49)),
        z_grid=list(np.linspace(-Zwall, Zwall, 40)))
    t = openmc.Tally(name="flux"); t.filters = [openmc.MeshFilter(mesh)]; t.scores = ["flux"]
    return openmc.Model(openmc.Geometry([cell]), openmc.Materials([void]),
                        s, openmc.Tallies([t]))


def main():
    stem = str(HERE.parent / "data" / "precise_QA_cm_fluxmap")
    Rlo, Rhi, Zmax = make_cm_fluxmap(stem)
    print(f"precise_QA (cm): R[{Rlo:.0f},{Rhi:.0f}] Zmax {Zmax:.0f}; running 3 modes ...", flush=True)
    import json
    dims = None
    for name, pol in MODES.items():
        m = build(pol, stem, Rlo, Rhi, Zmax)
        sp = m.run(cwd=f"rd_{name}")
        with openmc.StatePoint(sp) as s:
            t = s.get_tally(name="flux")
            mesh = t.find_filter(openmc.MeshFilter).mesh
            dims = tuple(int(x) for x in mesh.dimension)
            arr = np.asarray(t.mean).reshape(dims)   # (nR, nphi, nZ)
        np.save(HERE.parent / "data" / f"precise_QA_flux_{name}.npy", arr)
        print(f"[{name}] flux map {dims}  total {arr.sum():.3e}", flush=True)
    print("wrote data/precise_QA_flux_{unpol,perp,par}.npy  dims", dims)
    print("DONE_REALDEV")


if __name__ == "__main__":
    main()
