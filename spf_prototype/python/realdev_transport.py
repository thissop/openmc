#!/usr/bin/env python
"""Real-device StellaratorSource transport with a WALL-CURRENT tally (steering is visible, unlike a
volume-flux map). Births from the real 803097 VMEC equilibrium (real sqrt(g) + real b_hat), free-
streaming to a cylindrical wall; the MeshSurfaceFilter 'current' on the inboard (r-min) and outboard
(r-max) faces gives the NWL redistribution under unpolarized / perpendicular / parallel emission.

Cluster env = schwartz_capstone.py's (PR openmc + built lib via LD_PRELOAD).
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
ID = 803097


def make_cm_fluxmap(stem):
    f = dict(np.load(HERE.parent / "data" / f"quasr{ID}_vmec_fluxmap.npz"))
    write_fluxmap_bin(stem, f["rho"], f["theta"], f["zeta"], f["sqrtg"],
                      f["R"] * CM, f["Z"] * CM, f["phi"], f["BR"], f["Bphi"], f["BZ"],
                      int(f["nfp"]), float(f["sign_sqrtg"]))
    R = f["R"] * CM; Z = f["Z"] * CM
    return float(R.min()), float(R.max()), float(np.abs(Z).max()), int(f["nfp"])


NZ = 60


def build(pol, stem, Rlo, Rhi, Zmax):
    gap = 0.25 * (Rhi - Rlo) + 3.0
    Rin, Rout, Zw = Rlo - gap, Rhi + gap, Zmax + gap
    void = openmc.Material(name="void"); void.add_nuclide("H1", 1.0); void.set_density("g/cm3", 1e-10)
    cin = openmc.ZCylinder(r=Rin, boundary_type="vacuum")
    cout = openmc.ZCylinder(r=Rout, boundary_type="vacuum")
    zt = openmc.ZPlane(z0=Zw, boundary_type="vacuum"); zb = openmc.ZPlane(z0=-Zw, boundary_type="vacuum")
    cell = openmc.Cell(fill=void, region=+cin & -cout & +zb & -zt)
    src = openmc.StellaratorSource(fluxmap=stem, polarization=pol, field_model="fluxmap",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.source = [src]
    s.particles = 800_000; s.batches = 25; s.photon_transport = False
    mesh = openmc.CylindricalMesh(r_grid=[Rin, Rout], phi_grid=[0.0, 2 * np.pi],
                                  z_grid=list(np.linspace(-Zw, Zw, NZ + 1)))
    t = openmc.Tally(name="wall"); t.filters = [openmc.MeshSurfaceFilter(mesh)]; t.scores = ["current"]
    return openmc.Model(openmc.Geometry([cell]), openmc.Materials([void]), s,
                        openmc.Tallies([t])), (Rin, Rout, Zw)


def main():
    stem = str(HERE.parent / "data" / f"quasr{ID}_cm_fluxmap")
    Rlo, Rhi, Zmax, nfp = make_cm_fluxmap(stem)
    print(f"{ID} (cm): R[{Rlo:.0f},{Rhi:.0f}] Zmax {Zmax:.0f} nfp {nfp}", flush=True)
    prof = {}
    geom = None
    for name, pol in MODES.items():
        m, geom = build(pol, stem, Rlo, Rhi, Zmax)
        sp = m.run(cwd=f"rt_{ID}_{name}")   # per-device cwd: array tasks must not share a run dir
        with openmc.StatePoint(sp) as s:
            mean = np.asarray(s.get_tally(name="wall").mean).ravel()
        nsurf = mean.size // NZ
        arr = mean.reshape(NZ, nsurf)
        prof[name] = dict(inb=arr[:, 0], out=arr[:, 2])   # surface 0 = inboard (r-min), 2 = outboard
        print(f"[{name}] nsurf={nsurf}", flush=True)
    Rin, Rout, Zw = geom
    zc = 0.5 * (np.linspace(-Zw, Zw, NZ + 1)[:-1] + np.linspace(-Zw, Zw, NZ + 1)[1:])
    np.savez(HERE.parent / "data" / f"quasr{ID}_wallcurrent.npz", zc=zc, geom=np.array([Rin, Rout, Zw]),
             **{f"{k}_{s}": prof[k][s] for k in MODES for s in ("inb", "out")})
    mid = int(np.argmin(np.abs(zc)))
    print(f"--- 803097 midplane steering (vs unpolarized) ---")
    for s, lab in [("inb", "inboard "), ("out", "outboard")]:
        u = prof["unpol"][s][mid]
        for name in ("perp", "par"):
            print(f"  {lab} {name}: {100*(prof[name][s][mid]-u)/u:+.1f}%", flush=True)
    print("wrote data/quasr%d_wallcurrent.npz" % ID)
    print("DONE_REALDEV_TRANSPORT")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        ID = int(sys.argv[1])
    main()
