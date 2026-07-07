#!/usr/bin/env python
"""CAPSTONE: does the NATIVE SPF TokamakSource reproduce Schwartz's inboard/outboard NWL steering?

Square-cross-section torus, parabolic (1-r^2) plasma, toroidal field -- the exact anarrima
oracle geometry (anarrima/examples/square_torus.py), scaled x100 to cm. Free-streaming (near-void).
Runs unpolarized / perpendicular(A) / parallel(B) and compares the inboard & outboard MIDPLANE NWL
to the analytic oracle from python/anarrima_oracle.py (perp/A: +40.6% inboard, -21.2% outboard;
par/BC the mirror).

Run on the cluster (built OpenMC PR + native SPF TokamakSource):
    export PYTHONPATH=$SRC:$PYTHONPATH LD_PRELOAD=$SRC/build/lib/libopenmc.so PATH=$SRC/build/bin:$PATH
    export OPENMC_CROSS_SECTIONS=.../cross_sections.xml
    python schwartz_capstone.py
"""
import numpy as np
import openmc

SCALE = 100.0                      # anarrima units -> cm
R0, A = 1.0 * SCALE, 0.5 * SCALE
RIN, ROUT, ZTB = 0.4 * SCALE, 1.6 * SCALE, 0.6 * SCALE
DZ_MID = 0.05 * SCALE              # +/- midplane band for the "midplane" NWL
NZ = 40                            # poloidal (z) resolution of the wall tally
MODES = {"unpol": (1/3, 1/3, 1/3), "perp": (1.0, 0.0, 0.0), "par": (0.0, 1.0, 0.0)}


def build(polarization):
    void = openmc.Material(name="void"); void.add_nuclide("H1", 1.0)
    void.set_density("g/cm3", 1e-10)                      # ~free-streaming
    mats = openmc.Materials([void])

    cin = openmc.ZCylinder(r=RIN, boundary_type="vacuum")
    cout = openmc.ZCylinder(r=ROUT, boundary_type="vacuum")
    ztop = openmc.ZPlane(z0=ZTB, boundary_type="vacuum")
    zbot = openmc.ZPlane(z0=-ZTB, boundary_type="vacuum")
    cell = openmc.Cell(fill=void, region=+cin & -cout & +zbot & -ztop)
    geom = openmc.Geometry([cell])

    roa = np.linspace(0.0, 1.0, 40)
    src = openmc.TokamakSource(
        major_radius=R0, minor_radius=A, elongation=1.0, triangularity=0.0,
        shafranov_shift=0.0, r_over_a=roa, emission_density=1.0 - roa ** 2,
        energy=openmc.stats.Discrete([14.1e6], [1.0]),
        polarization=polarization, field_model="toroidal")

    s = openmc.Settings()
    s.run_mode = "fixed source"; s.source = [src]
    s.particles = 200_000; s.batches = 20
    s.photon_transport = False

    # NWL: net current on the inboard (r=RIN) and outboard (r=ROUT) walls, poloidally (z) resolved.
    mesh = openmc.CylindricalMesh(
        r_grid=[RIN, ROUT], phi_grid=[0.0, 2 * np.pi],
        z_grid=list(np.linspace(-ZTB, ZTB, NZ + 1)))
    t = openmc.Tally(name="wall"); t.filters = [openmc.MeshSurfaceFilter(mesh)]
    t.scores = ["current"]
    return openmc.Model(geom, mats, s, openmc.Tallies([t]))


def wall_midplane(sp_path, mesh_nz=NZ):
    """Return (inboard, outboard) midplane NWL from the MeshSurfaceFilter 'current' tally.
    CylindricalMesh (1 r-cell, 1 phi, NZ z-cells): reduce the surface bins to the r-min face
    (inboard) and r-max face (outboard) in the midplane z-band. Prints the surface structure once."""
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="wall")
        mean = np.asarray(t.mean).ravel()
    nz = mesh_nz
    nsurf = mean.size // nz                        # surfaces per (r,phi,z) cell
    arr = mean.reshape(1, 1, nz, nsurf)            # (nr=1, nphi=1, nz, nsurf)
    zc = np.linspace(-ZTB, ZTB, nz + 1); zmid = 0.5 * (zc[:-1] + zc[1:])
    midband = np.abs(zmid) < DZ_MID
    # per-surface midplane totals (diagnostic); the physical faces are the r-min/r-max ones.
    persurf = np.array([arr[0, 0, midband, si].sum() for si in range(nsurf)])
    return persurf, nsurf


def main():
    import json
    # oracle (from python/anarrima_oracle.py, run locally): perp/A +40.6% inboard / -21.2% outboard
    results = {}
    for name, pol in MODES.items():
        m = build(pol)
        sp = m.run(cwd=f"cap_{name}")
        persurf, nsurf = wall_midplane(sp)
        results[name] = persurf.tolist()
        print(f"[{name}] midplane per-surface currents ({nsurf} surfaces): "
              + "  ".join(f"{v:+.3e}" for v in persurf))
    json.dump(results, open("schwartz_capstone_persurf.json", "w"), indent=2)
    print("\nORACLE (target):  perp/A +40.6% inboard / -21.2% outboard ; par/BC the mirror")
    print("Wrote schwartz_capstone_persurf.json -- pick the r-min (inboard) & r-max (outboard) "
          "face indices from the unpol run, then form perp/unpol and par/unpol ratios.")


if __name__ == "__main__":
    main()
