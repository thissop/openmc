#!/usr/bin/env python
"""REAL neutronics (θ,φ) first-wall NWL maps: native StellaratorSource on a device's VMEC equilibrium
(real √g births + real b̂), free-streaming through near-void to a TORUS first wall, with every neutron's
wall crossing captured by surface_source_write and histogrammed in (poloidal θ, toroidal φ). No analytic
free-streaming kernel -- this is OpenMC transport of the compiled source. Three modes (unpol/perp/par)
-> the maps the a₂ optimizer combines (NWL is exactly linear in a₂; verified downstream).

Run per device on Ginsburg (capstone env). Saves data/quasr<ID>_wallmap.npz (small histograms only).
"""
import sys
from pathlib import Path
import numpy as np
import h5py
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from quasr_fluxmap import write_fluxmap_bin

CM = 100.0
MODES = {"unpol": (1/3, 1/3, 1/3), "perp": (1.0, 0.0, 0.0), "par": (0.0, 1.0, 0.0)}
NTH, NPH = 48, 96
WALL_GAP = 1.2           # torus minor radius = WALL_GAP * plasma cross-section half-size
PARTICLES = 3_000_000


def make_cm_fluxmap(ID, stem):
    f = dict(np.load(HERE.parent / "data" / f"quasr{ID}_vmec_fluxmap.npz"))
    write_fluxmap_bin(stem, f["rho"], f["theta"], f["zeta"], f["sqrtg"],
                      f["R"] * CM, f["Z"] * CM, f["phi"], f["BR"], f["Bphi"], f["BZ"],
                      int(f["nfp"]), float(f["sign_sqrtg"]))
    R = f["R"] * CM; Z = f["Z"] * CM
    R0 = 0.5 * (R.max() + R.min())
    aR = 0.5 * (R.max() - R.min()); aZ = 0.5 * (Z.max() - Z.min())
    return R0, max(aR, aZ)   # plasma cross-section half-size (encloses both R and Z extent)


def run_mode(pol, stem, R0, bwall, cwd):
    void = openmc.Material(name="void"); void.add_nuclide("H1", 1.0); void.set_density("g/cm3", 1e-10)
    torus = openmc.ZTorus(x0=0, y0=0, z0=0, a=R0, b=bwall, c=bwall, boundary_type="vacuum")
    cell = openmc.Cell(fill=void, region=-torus)
    src = openmc.StellaratorSource(fluxmap=stem, polarization=pol, field_model="fluxmap",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.source = [src]
    s.particles = PARTICLES; s.batches = 1; s.photon_transport = False
    s.surf_source_write = {"surface_ids": [torus.id], "max_particles": PARTICLES}  # 0.15.3 name
    m = openmc.Model(openmc.Geometry([cell]), openmc.Materials([void]), s)
    m.run(cwd=cwd, output=False)
    with h5py.File(Path(cwd) / "surface_source.h5", "r") as fh:
        sb = fh["source_bank"][()]
    r = sb["r"]
    x, y, z = r["x"], r["y"], r["z"]
    phi = np.mod(np.arctan2(y, x), 2 * np.pi)
    R = np.hypot(x, y)
    th = np.mod(np.arctan2(z, R - R0), 2 * np.pi)
    H, _, _ = np.histogram2d(th, phi, bins=[NTH, NPH], range=[[0, 2*np.pi], [0, 2*np.pi]])
    return H, x.size


def main():
    ID = int(sys.argv[1]) if len(sys.argv) > 1 else 803097
    stem = str(HERE.parent / "data" / f"quasr{ID}_cm_fluxmap")
    R0, awall = make_cm_fluxmap(ID, stem)
    bwall = min(WALL_GAP * awall, 0.80 * R0)   # snug torus first wall, kept non-degenerate (bwall<R0)
    assert bwall > 1.02 * awall, f"plasma a={awall:.0f} does not fit non-degenerate torus b={bwall:.0f}"
    thc = (np.arange(NTH) + 0.5) * (2 * np.pi / NTH)
    dth, dph = 2 * np.pi / NTH, 2 * np.pi / NPH
    dA = (R0 + bwall * np.cos(thc)) * bwall * dth * dph   # (NTH,) torus area element
    print(f"{ID}: R0={R0:.0f}cm torus b={bwall:.0f}cm; 3 modes x {PARTICLES:.0e}", flush=True)
    out = dict(theta=thc, phi=(np.arange(NPH) + 0.5) * dph, dA=dA, R0=R0, bwall=bwall)
    for name, pol in MODES.items():
        H, n = run_mode(pol, stem, R0, bwall, f"/tmp/wm_{ID}_{name}")
        nwl = H / dA[:, None]                              # NWL(θ,φ) = crossings / wall area
        out[name] = nwl
        print(f"[{name}] {n} wall crossings, peak/mean={nwl.max()/nwl.mean():.3f}", flush=True)
    np.savez(HERE.parent / "data" / f"quasr{ID}_wallmap.npz", **out)
    # linearity check (the physics, not an assumption): perp + par == 2*unpol
    lin = np.abs(out["perp"] + out["par"] - 2 * out["unpol"])
    scale = out["unpol"].mean()
    print(f"linearity |perp+par-2unpol|/<unpol>: max {lin.max()/scale:.3f} mean {lin.mean()/scale:.4f}")
    print(f"wrote data/quasr{ID}_wallmap.npz")
    print("DONE_WALLMAP")


if __name__ == "__main__":
    main()
