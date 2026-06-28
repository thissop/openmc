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
import openmc  # noqa: E402
import toy_qalow_config as C  # noqa: E402
import spf_mirror as mir  # noqa: E402

N_PER_MODE = 300_000
NR = NZ = 16                  # match analytic n_per_wall
PHI_J = np.linspace(0, 2 * np.pi, 256, endpoint=False)


def presample(mode_abc, seed_pos, seed_dir):
    """Birth bank for one mode: same positions (seed_pos) for all modes; directions
    (seed_dir) about the local toy field. wgt = w_loop * J(phi) (arc-length)."""
    m = mir.make_mode_weights(*mode_abc)
    loops = C.loops()
    rng_pos = np.random.default_rng(seed_pos)
    rng_dir = mir.Rng(seed_dir)
    # choose loops uniformly; weight carries w_loop. phi uniform; weight carries J.
    li = rng_pos.integers(0, len(loops), N_PER_MODE)
    phi = rng_pos.uniform(0, 2 * np.pi, N_PER_MODE)
    parts = []
    for k in range(N_PER_MODE):
        rho, th, w = loops[li[k]]
        ph = np.array([phi[k]])
        R, Z = C.loop_RZ(th, ph, rho)
        # arc-length Jacobian J = sqrt(R^2 + R'^2 + Z'^2) via finite diff in phi
        dph = 1e-4
        Rp, Zp = C.loop_RZ(th, ph + dph, rho)
        Rm, Zm = C.loop_RZ(th, ph - dph, rho)
        dR = (Rp - Rm) / (2 * dph); dZ = (Zp - Zm) / (2 * dph)
        J = float(np.sqrt(R[0] ** 2 + dR[0] ** 2 + dZ[0] ** 2))
        x, y, z = R[0] * np.cos(phi[k]), R[0] * np.sin(phi[k]), Z[0]
        Bhat = tuple(C.field_bhat(th, ph)[0])
        u = mir.sample_global_direction(m, Bhat, rng_dir)
        parts.append(openmc.SourceParticle(r=(x, y, z), u=u, E=14.1e6,
                                           wgt=w * J, particle=openmc.ParticleType.NEUTRON))
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
    s.run_mode = "fixed source"; s.particles = N_PER_MODE; s.batches = 10; s.inactive = 0
    s.source = openmc.FileSource(src_path)
    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(C.R_IN, C.R_OUT, NR + 1),
        phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-C.Z_W, C.Z_W, NZ + 1))
    t = openmc.Tally(name="wall"); t.filters = [openmc.MeshSurfaceFilter(mesh)]
    t.scores = ["current"]
    return openmc.Model(geometry=geom, settings=s, tallies=openmc.Tallies([t]))


def extract(sp_path):
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="wall").get_pandas_dataframe()
    xi = df[("mesh 1", "x")].to_numpy(); zi = df[("mesh 1", "z")].to_numpy()
    surf = df[("mesh 1", "surf")].to_numpy(); mean = df[("mean", "")].to_numpy()
    out = {}
    out["inboard"] = mean[(surf == "x-min out") & (xi == 1)].sum()
    out["outboard"] = mean[(surf == "x-max out") & (xi == NR)].sum()
    out["floor"] = mean[(surf == "z-min out") & (zi == 1)].sum()
    out["ceiling"] = mean[(surf == "z-max out") & (zi == NZ)].sum()
    return out


def main():
    openmc.config["cross_sections"] = os.environ.get(
        "OPENMC_CROSS_SECTIONS", os.path.expanduser("~/nndc_hdf5/cross_sections.xml"))
    cur = {}
    for mode, abc in C.MODES.items():
        parts = presample(abc, seed_pos=12345, seed_dir=999 + hash(mode) % 1000)
        src = f"/tmp/toy_qalow_src_{mode}.h5"
        openmc.write_source_file(parts, src)
        sp = make_model(src).run(cwd=f"/tmp/toy_qalow_mc_{mode}", output=False)
        cur[mode] = extract(sp)
        print(f"ran {mode}: " + ", ".join(f"{w}={cur[mode][w]:.3e}" for w in cur[mode]))

    an = np.load("/tmp/toy_qalow_analytic.npz")
    an_wall = an["wall"]
    print(f"\n{'wall':>9} {'A/iso MC':>10} {'A/iso ana':>10} {'relerr':>7} | "
          f"{'B/iso MC':>10} {'B/iso ana':>10} {'relerr':>7}")
    worst = 0.0
    for w in ("inboard", "outboard", "floor", "ceiling"):
        a_mc = cur["A"][w] / cur["iso"][w]; b_mc = cur["B"][w] / cur["iso"][w]
        msk = an_wall == w
        a_an = float(np.nanmean(an["A_iso"][msk])); b_an = float(np.nanmean(an["B_iso"][msk]))
        ea = abs(a_mc - a_an) / abs(a_an); eb = abs(b_mc - b_an) / abs(b_an)
        worst = max(worst, ea, eb)
        print(f"{w:>9} {a_mc:10.3f} {a_an:10.3f} {ea:7.1%} | "
              f"{b_mc:10.3f} {b_an:10.3f} {eb:7.1%}")
    print(f"\nworst-wall directionality discrepancy: {worst:.1%} "
          f"(OpenMC free-streaming vs analytic, toy QA-low, eps_eff~0.08)")


if __name__ == "__main__":
    main()
