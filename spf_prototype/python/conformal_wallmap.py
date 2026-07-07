#!/usr/bin/env python
"""FREE-STREAMING NWL on the DAGMC CONFORMAL wall (real transport, real conformal geometry).

Per device: build a conformal first wall = the VMEC LCFS offset outward along its poloidal normal
(stellarator_geometry.build_layers), convert the shell STL to a watertight DAGMC .h5m
(dagmc_writer.build_from_stls, pymoab-free), then run the native StellaratorSource (VMEC √g births +
real b̂) free-streaming through near-void inside it. Every crossing of the plasma-facing wall surface
is captured with surf_source_write and mapped (KD-tree on the wall grid) to (poloidal θ, toroidal φ)
-> NWL(θ,φ). Three modes (unpol/perp/par) -> the maps the a₂ optimizer combines.

Run on Ginsburg (PR openmc built WITH DAGMC + the native StellaratorSource).
  python conformal_wallmap.py <ID>
"""
import sys
from pathlib import Path
import numpy as np
import h5py
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stellarator_geometry as sg
import dagmc_writer as dw
from quasr_fluxmap import write_fluxmap_bin

sg.DATADIR = HERE.parent / "data"        # anchor geometry I/O to THIS checkout's data dir

CM = 100.0
GAP_FRAC = 0.30                 # first wall = LCFS + 0.30 * a_minor along poloidal normal
FW_THICK = 1.0                  # cm shell (near-void; only its inner surface is the wall)
PARTICLES = 3_000_000
MODES = {"unpol": (1/3, 1/3, 1/3), "perp": (1.0, 0.0, 0.0), "par": (0.0, 1.0, 0.0)}


def build_wall(ID):
    """Conformal wall .h5m + the (θ,φ) wall grid (for crossing->(θ,φ) mapping + areas)."""
    f = dict(np.load(HERE.parent / "data" / f"quasr{ID}_vmec_fluxmap.npz"))
    Rl = f["R"][-1] * CM; Zl = f["Z"][-1] * CM          # LCFS (n_theta, n_zeta) in cm
    nfp = int(f["nfp"])
    nR, nZ = sg.poloidal_outward_normals(Rl, Zl)
    a = 0.5 * (Rl.max() - Rl.min())
    Rw, Zw = sg.offset_surface(Rl, Zl, nR, nZ, GAP_FRAC * a)   # conformal first wall
    stem = f"quasr{ID}_wall"
    np.savez(sg.DATADIR / f"{stem}_surface.npz", R=Rw, Z=Zw, nfp=nfp)
    sg.build_layers(stem, layers=[("firstwall", FW_THICK)], scale=1.0)
    h5m = sg.DATADIR / f"{stem}.h5m"
    dw.build_from_stls(sg.DATADIR / f"{stem}_geom", str(h5m))
    # wall grid points (θ,φ) -> Cartesian, + per-cell area
    nt, nz = Rw.shape
    phi = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    P = sg._verts(Rw, Zw, phi)                            # (nt, nz, 3)
    dth, dph = 2 * np.pi / nt, 2 * np.pi / nz
    Pt = (np.roll(P, -1, 0) - np.roll(P, 1, 0)) / (2 * dth)
    Pp = (np.roll(P, -1, 1) - np.roll(P, 1, 1)) / (2 * dph)
    dA = np.linalg.norm(np.cross(Pt, Pp), axis=-1) * dth * dph   # (nt, nz)
    return str(h5m), P.reshape(-1, 3), dA, (nt, nz)


def wall_surface_id(dag):
    """DAGMC surface id of the plasma-facing (inner) wall surface S0."""
    ids = sorted(int(s) for s in dag.surface_ids) if hasattr(dag, "surface_ids") else []
    return ids[0] if ids else 1


def run_mode(pol, fluxstem, h5m, cwd, sid):
    import openmc
    fw = openmc.Material(name="firstwall"); fw.add_nuclide("H1", 1.0); fw.set_density("g/cm3", 1e-10)
    vac = openmc.Material(name="vacuum"); vac.add_nuclide("H1", 1.0); vac.set_density("g/cm3", 1e-12)
    dag = openmc.DAGMCUniverse(h5m)
    geom = openmc.Geometry(dag.bounded_universe())
    src = openmc.StellaratorSource(fluxmap=fluxstem, polarization=pol, field_model="fluxmap",
                                   energy=openmc.stats.Discrete([14.1e6], [1.0]))
    s = openmc.Settings(); s.run_mode = "fixed source"; s.source = [src]
    s.particles = PARTICLES; s.batches = 1; s.photon_transport = False
    s.surf_source_write = {"surface_ids": [sid], "max_particles": PARTICLES}
    m = openmc.Model(geom, openmc.Materials([fw, vac]), s)
    m.run(cwd=cwd, output=False)
    with h5py.File(Path(cwd) / "surface_source.h5", "r") as fh:
        r = fh["source_bank"][()]["r"]
    return np.stack([r["x"], r["y"], r["z"]], axis=1)


def main():
    ID = int(sys.argv[1]) if len(sys.argv) > 1 else 803097
    fluxstem = str(HERE.parent / "data" / f"quasr{ID}_cm_fluxmap")
    f = dict(np.load(HERE.parent / "data" / f"quasr{ID}_vmec_fluxmap.npz"))
    write_fluxmap_bin(fluxstem, f["rho"], f["theta"], f["zeta"], f["sqrtg"],
                      f["R"] * CM, f["Z"] * CM, f["phi"], f["BR"], f["Bphi"], f["BZ"],
                      int(f["nfp"]), float(f["sign_sqrtg"]))
    h5m, Pw, dA, (nt, nz) = build_wall(ID)
    import openmc
    sid = wall_surface_id(openmc.DAGMCUniverse(h5m))
    tree = cKDTree(Pw)
    print(f"{ID}: conformal wall {nt}x{nz}, DAGMC surface id {sid}; 3 modes x {PARTICLES:.0e}", flush=True)
    out = dict(dA=dA, shape=np.array([nt, nz]))
    for name, pol in MODES.items():
        xyz = run_mode(pol, fluxstem, h5m, f"/tmp/cw_{ID}_{name}", sid)
        _, idx = tree.query(xyz)                          # nearest wall grid cell
        H = np.bincount(idx, minlength=nt * nz).reshape(nt, nz).astype(float)
        nwl = H / dA
        out[name] = nwl
        print(f"[{name}] {xyz.shape[0]} crossings, PF={nwl.max()/nwl.mean():.3f}", flush=True)
    np.savez(HERE.parent / "data" / f"quasr{ID}_conformalmap.npz", **out)
    lin = np.abs(out["perp"] + out["par"] - 2 * out["unpol"]).mean() / out["unpol"].mean()
    print(f"MC linearity |perp+par-2unpol|/<unpol>: {lin:.4f}")
    print(f"wrote data/quasr{ID}_conformalmap.npz\nDONE_CONFORMAL")


if __name__ == "__main__":
    main()
