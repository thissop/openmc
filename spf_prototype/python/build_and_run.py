#!/usr/bin/env python
"""Tier-2b driver: build the polarized-fusion CompiledSource .so, build the
square-cross-section torus OpenMC model (near-void free-streaming), run each
polarization config, and extract per-wall poloidal NWL patterns.

Geometry matches analytic_nwl exactly (R0=1, a=0.5, u=0.4, w=1.6, z=±0.6); units
are arbitrary because the interior is void (infinite mean free path). Self-test
(`python build_and_run.py`) checks current conservation and the isotropic
inboard/outboard ratio against the analytic.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
PYDIR = REPO / "spf_prototype" / "python"
SRCDIR = REPO / "spf_prototype" / "src"
SO = SRCDIR / "build" / "libpolarized_fusion_source.so"
VENV = Path(os.path.expanduser("~/spf_venv"))
XS = Path(os.path.expanduser("~/nndc_hdf5/cross_sections.xml"))

sys.path.insert(0, str(PYDIR))
import openmc  # noqa: E402
import analytic_nwl as an  # noqa: E402

# Geometry (normalized, == analytic_nwl module constants)
U, W, ZW = an.R_IN, an.R_OUT, an.Z_WALL   # 0.4, 1.6, 0.6
R0, AM = an.R0, an.A_MINOR                 # 1.0, 0.5


def build_so():
    SO.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cmake", "-B", str(SO.parent), f"-DCMAKE_PREFIX_PATH={VENV}", str(SRCDIR)],
                   check=True, capture_output=True)
    subprocess.run(["cmake", "--build", str(SO.parent)], check=True, capture_output=True)
    return SO


def make_model(params, particles, batches=20, nr=60, nz=60):
    openmc.reset_auto_ids()
    inner = openmc.ZCylinder(r=U, boundary_type="vacuum")
    outer = openmc.ZCylinder(r=W, boundary_type="vacuum")
    bot = openmc.ZPlane(z0=-ZW, boundary_type="vacuum")
    top = openmc.ZPlane(z0=ZW, boundary_type="vacuum")
    vessel = openmc.Cell(region=+inner & -outer & +bot & -top, fill=None)
    geom = openmc.Geometry([vessel])

    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.particles = int(particles)
    s.batches = int(batches)
    s.inactive = 0
    s.source = openmc.CompiledSource(str(SO), parameters=params)

    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(U, W, nr + 1),
        phi_grid=np.array([0.0, 2 * np.pi]),
        z_grid=np.linspace(-ZW, ZW, nz + 1),
    )
    t = openmc.Tally(name="wall")
    t.filters = [openmc.MeshSurfaceFilter(mesh)]
    t.scores = ["current"]
    model = openmc.Model(geometry=geom, settings=s, tallies=openmc.Tallies([t]))
    return model, nr, nz


def run(params, particles, cwd, **kw):
    model, nr, nz = make_model(params, particles, **kw)
    sp_path = model.run(cwd=cwd, output=False)
    return sp_path, nr, nz


def extract_walls(sp_path, nr, nz):
    """Return dict of per-wall (centers, current, std) using the leaking ('out')
    partial current on the four boundary faces of the cylindrical mesh."""
    with openmc.StatePoint(sp_path) as sp:
        df = sp.get_tally(name="wall").get_pandas_dataframe()
    xi = df[("mesh 1", "x")].to_numpy()
    zi = df[("mesh 1", "z")].to_numpy()
    surf = df[("mesh 1", "surf")].to_numpy()
    mean = df[("mean", "")].to_numpy()
    sd = df[("std. dev.", "")].to_numpy()

    def pick(mask, order_idx, n):
        m = mean[mask]; s = sd[mask]; o = order_idx[mask]
        out_m = np.zeros(n); out_s = np.zeros(n)
        out_m[o - 1] = m; out_s[o - 1] = s
        return out_m, out_s

    z_centers = 0.5 * (np.linspace(-ZW, ZW, nz + 1)[:-1] + np.linspace(-ZW, ZW, nz + 1)[1:])
    r_centers = 0.5 * (np.linspace(U, W, nr + 1)[:-1] + np.linspace(U, W, nr + 1)[1:])
    dz = (2 * ZW) / nz
    redges = np.linspace(U, W, nr + 1)

    res = {}
    m = (surf == "x-min out") & (xi == 1)
    res["inboard"] = dict(centers=z_centers, area=2 * np.pi * U * dz,
                          **dict(zip(("cur", "std"), pick(m, zi, nz))))
    m = (surf == "x-max out") & (xi == nr)
    res["outboard"] = dict(centers=z_centers, area=2 * np.pi * W * dz,
                           **dict(zip(("cur", "std"), pick(m, zi, nz))))
    m = (surf == "z-min out") & (zi == 1)
    res["floor"] = dict(centers=r_centers, area=np.pi * (redges[1:] ** 2 - redges[:-1] ** 2),
                        **dict(zip(("cur", "std"), pick(m, xi, nr))))
    m = (surf == "z-max out") & (zi == nz)
    res["ceiling"] = dict(centers=r_centers, area=np.pi * (redges[1:] ** 2 - redges[:-1] ** 2),
                         **dict(zip(("cur", "std"), pick(m, xi, nr))))
    return res


def _selftest():
    openmc.config["cross_sections"] = str(XS)
    build_so()
    print("built", SO.name)
    params = f"a=0.333333,b=0.333333,c=0.333333,bmode=toroidal,shape=plasma,R0={R0},aminor={AM}"
    sp, nr, nz = run(params, 200000, "/tmp/spf_2b_selftest", batches=10, nr=50, nz=50)
    walls = extract_walls(sp, nr, nz)
    total = sum(walls[w]["cur"].sum() for w in walls)
    print(f"conservation: sum of boundary 'out' currents = {total:.4f} (expect ~1.0 per source particle)")
    # isotropic intensity at midplane inboard vs outboard
    mid = nz // 2
    ii = walls["inboard"]["cur"][mid] / walls["inboard"]["area"]
    io = walls["outboard"]["cur"][mid] / walls["outboard"]["area"]
    print(f"iso intensity outboard/inboard (MC, midplane) = {io/ii:.4f}  (analytic ~1.146)")


if __name__ == "__main__":
    _selftest()
