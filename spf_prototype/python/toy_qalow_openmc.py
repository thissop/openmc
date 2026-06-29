#!/usr/bin/env python
"""OpenMC side of the toy QA-low free-streaming cross-check (run in spf_venv).

Pre-samples the toy QA-low source (positions on the flux-surface loops, weighted by
loop weight x arc-length; directions drawn about the LOCAL toy field by the
bit-parity Python mirror of the compiled SPF sampler) as an OpenMC source bank, one
per mode, free-streams them into the SAME square-torus wall (void interior, CSG --
no DAGMC), tallies first-flight current per wall, and compares the directionality
(A/iso, B/iso) to the analytic result (toy_qalow_analytic.npz).

Run: PATH=$HOME/spf_venv/bin:$PATH OMP_NUM_THREADS=4 \
     $HOME/spf_venv/bin/python toy_qalow_openmc.py
"""
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
import importlib  # noqa: E402
import openmc  # noqa: E402
C = importlib.import_module(os.environ.get("SPF_CONFIG", "toy_qalow_config"))  # noqa: E402
import spf_mirror as mir  # noqa: E402

N_PER_MODE = 1_000_000
NR = NZ = 16                  # match analytic n_per_wall


def presample(mode_abc, seed_pos, seed_dir):
    """Birth bank for one mode: same positions (seed_pos) for all modes; directions
    (seed_dir) about the local field. wgt = w_loop * J(phi) (arc-length).

    Geometry+field are evaluated VECTORIZED (per rho group) -- the configs accept an
    array of theta -- so large banks are affordable; only the bit-parity direction
    draw stays per particle, preserving the mirror RNG sequence."""
    m = mir.make_mode_weights(*mode_abc)
    loops = C.loops()
    rng_pos = np.random.default_rng(seed_pos)
    rng_dir = mir.Rng(seed_dir)
    li = rng_pos.integers(0, len(loops), N_PER_MODE)
    phi = rng_pos.uniform(0, 2 * np.pi, N_PER_MODE)
    rho_arr = np.array([loops[i][0] for i in li])
    th_arr = np.array([loops[i][1] for i in li])
    w_arr = np.array([loops[i][2] for i in li])
    R = np.empty(N_PER_MODE); Z = np.empty(N_PER_MODE); J = np.empty(N_PER_MODE)
    Bh = np.empty((N_PER_MODE, 3))
    dph = 1e-4
    for rho in np.unique(rho_arr):
        sel = rho_arr == rho
        th = th_arr[sel]; ph = phi[sel]
        Rs, Zs = C.loop_RZ(th, ph, rho)
        Rp, Zp = C.loop_RZ(th, ph + dph, rho)
        Rm, Zm = C.loop_RZ(th, ph - dph, rho)
        dR = (Rp - Rm) / (2 * dph); dZ = (Zp - Zm) / (2 * dph)
        R[sel] = Rs; Z[sel] = Zs
        J[sel] = np.sqrt(Rs ** 2 + dR ** 2 + dZ ** 2)
        Bh[sel] = C.field_bhat(th, ph)
    x = R * np.cos(phi); y = R * np.sin(phi)
    parts = []
    for k in range(N_PER_MODE):
        u = mir.sample_global_direction(m, (Bh[k, 0], Bh[k, 1], Bh[k, 2]), rng_dir)
        parts.append(openmc.SourceParticle(r=(x[k], y[k], Z[k]), u=u, E=14.1e6,
                                           wgt=w_arr[k] * J[k],
                                           particle=openmc.ParticleType.NEUTRON))
    return parts


def make_model(src_path):
    openmc.reset_auto_ids()
    inner = openmc.ZCylinder(r=C.R_IN, boundary_type="vacuum")
    outer = openmc.ZCylinder(r=C.R_OUT, boundary_type="vacuum")
    bot = openmc.ZPlane(z0=-C.Z_W, boundary_type="vacuum")
    top = openmc.ZPlane(z0=C.Z_W, boundary_type="vacuum")
    vessel = openmc.Cell(region=+inner & -outer & +bot & -top, fill=None)
    geom = openmc.Geometry([vessel])
    s = openmc.Settings()
    s.run_mode = "fixed source"; s.particles = N_PER_MODE; s.batches = 5; s.inactive = 0
    s.source = openmc.FileSource(src_path)
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(C.R_IN, C.R_OUT, NR + 1),
        phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-C.Z_W, C.Z_W, NZ + 1))
    t = openmc.Tally(name="wall"); t.filters = [openmc.MeshSurfaceFilter(mesh)]
    t.scores = ["current"]
    return openmc.Model(geometry=geom, settings=s, tallies=openmc.Tallies([t]))


_ZC = 0.5 * (np.linspace(-C.Z_W, C.Z_W, NZ + 1)[:-1] + np.linspace(-C.Z_W, C.Z_W, NZ + 1)[1:])
_RC = 0.5 * (np.linspace(C.R_IN, C.R_OUT, NR + 1)[:-1] + np.linspace(C.R_IN, C.R_OUT, NR + 1)[1:])


def extract(sp_path):
    """Per wall: (s, per-bin current profile). Bin centers s = z (verticals) or r
    (horizontals); the wall total is profile.sum()."""
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="wall").get_pandas_dataframe()
    xi = df[("mesh 1", "x")].to_numpy(); zi = df[("mesh 1", "z")].to_numpy()
    surf = df[("mesh 1", "surf")].to_numpy(); mean = df[("mean", "")].to_numpy()

    def prof(mask, idx, n, s):
        v = np.zeros(n)
        v[idx[mask] - 1] = mean[mask]
        return s, v

    return {
        "inboard": prof((surf == "x-min out") & (xi == 1), zi, NZ, _ZC),
        "outboard": prof((surf == "x-max out") & (xi == NR), zi, NZ, _ZC),
        "floor": prof((surf == "z-min out") & (zi == 1), xi, NR, _RC),
        "ceiling": prof((surf == "z-max out") & (zi == NZ), xi, NR, _RC),
    }


def main():
    openmc.config["cross_sections"] = os.environ.get(
        "OPENMC_CROSS_SECTIONS", os.path.expanduser("~/nndc_hdf5/cross_sections.xml"))
    walls = ("inboard", "outboard", "floor", "ceiling")
    cur = {}        # cur[mode][wall] = (s, per-bin profile)
    for mode, abc in C.MODES.items():
        parts = presample(abc, seed_pos=12345, seed_dir=999 + hash(mode) % 1000)
        src = f"/tmp/{C.STEM}_src_{mode}.h5"
        openmc.write_source_file(parts, src)
        cwd = f"/tmp/{C.STEM}_mc_{mode}"; os.makedirs(cwd, exist_ok=True)
        sp = make_model(src).run(cwd=cwd, output=False)
        cur[mode] = extract(sp)
        print(f"ran {mode}: " + ", ".join(f"{w}={cur[mode][w][1].sum():.3e}" for w in walls))

    # save per-bin profiles (for the paper figure) -> s and per-mode current per wall
    save = {}
    for w in walls:
        save[f"s_{w}"] = cur["iso"][w][0]
        for m in C.MODES:
            save[f"{m}_{w}"] = cur[m][w][1]
    np.savez(f"/tmp/{C.STEM}_openmc.npz", **save)
    print(f"saved /tmp/{C.STEM}_openmc.npz (per-bin profiles for plotting)")

    apath = f"/tmp/{C.STEM}_analytic.npz"
    if not os.path.exists(apath):
        print(f"\n[skip comparison] run the analytic side first: {apath} not found")
        return
    an = np.load(apath)
    an_wall = an["wall"]
    print(f"\n{'wall':>9} {'A/iso MC':>10} {'A/iso ana':>10} {'relerr':>7} | "
          f"{'B/iso MC':>10} {'B/iso ana':>10} {'relerr':>7}")
    worst = 0.0
    for w in walls:
        iso_w = cur["iso"][w][1].sum()
        a_mc = cur["A"][w][1].sum() / iso_w; b_mc = cur["B"][w][1].sum() / iso_w
        # apples-to-apples: ratio-of-sums on BOTH sides (total A load / total iso
        # load over the wall), not MC ratio-of-sums vs analytic mean-of-ratios.
        msk = an_wall == w
        a_an = an["A"][msk].sum() / an["iso"][msk].sum()
        b_an = an["B"][msk].sum() / an["iso"][msk].sum()
        ea = abs(a_mc - a_an) / abs(a_an); eb = abs(b_mc - b_an) / abs(b_an)
        worst = max(worst, ea, eb)
        print(f"{w:>9} {a_mc:10.3f} {a_an:10.3f} {ea:7.1%} | "
              f"{b_mc:10.3f} {b_an:10.3f} {eb:7.1%}")
    print(f"\nworst-wall directionality discrepancy: {worst:.1%} "
          f"(OpenMC free-streaming vs analytic, config={C.STEM})")


if __name__ == "__main__":
    main()
