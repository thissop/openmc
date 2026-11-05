"""Validate analytic VMEC B-hat reconstruction against wout field components.

The production ``StellaratorSource`` reconstructs the field-line direction from
``iotas`` and ``lmns``.  This script compares that route with the independently
exported VMEC contravariant components ``bsupumnc`` and ``bsupvmnc`` at native
half-mesh surfaces.  Both are mapped through the same physical coordinate basis,
so the reported angular error isolates the field reconstruction rather than the
choice of radial geometry interpolant.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


def _read(path):
    names = (
        "ns", "xm", "xn", "xm_nyq", "xn_nyq", "rmnc", "zmns",
        "rmns", "zmnc", "lmns", "lmnc", "iotas", "bsupumnc",
        "bsupvmnc", "bsupumns", "bsupvmns",
    )
    with netcdf_file(str(path), "r", mmap=False) as data:
        return {
            name: np.asarray(data.variables[name].data).copy()
            for name in names if name in data.variables
        }


def _geometry_basis(data, rho, theta, zeta):
    """Evaluate e_theta and e_zeta with production's linear-in-rho geometry."""
    ns = int(data["ns"])
    rho_grid = np.sqrt(np.linspace(0.0, 1.0, ns))
    upper = np.searchsorted(rho_grid, rho, side="right")
    lower = np.clip(upper - 1, 0, ns - 2)
    fraction = ((rho - rho_grid[lower])
                / (rho_grid[lower + 1] - rho_grid[lower]))
    angle = theta[:, None] * data["xm"] - zeta[:, None] * data["xn"]
    cosine, sine = np.cos(angle), np.sin(angle)
    m, n = data["xm"], data["xn"]

    def interp(name):
        values = data[name]
        return values[lower] + fraction[:, None] * (
            values[lower + 1] - values[lower])

    rc, zs = interp("rmnc"), interp("zmns")
    radius = np.einsum("ij,ij->i", rc, cosine)
    r_theta = np.einsum("ij,ij->i", -m * rc, sine)
    r_zeta = np.einsum("ij,ij->i", n * rc, sine)
    z_theta = np.einsum("ij,ij->i", m * zs, cosine)
    z_zeta = np.einsum("ij,ij->i", -n * zs, cosine)
    if "rmns" in data and "zmnc" in data:
        rs, zc = interp("rmns"), interp("zmnc")
        radius += np.einsum("ij,ij->i", rs, sine)
        r_theta += np.einsum("ij,ij->i", m * rs, cosine)
        r_zeta += np.einsum("ij,ij->i", -n * rs, cosine)
        z_theta += np.einsum("ij,ij->i", -m * zc, sine)
        z_zeta += np.einsum("ij,ij->i", n * zc, sine)

    cos_zeta, sin_zeta = np.cos(zeta), np.sin(zeta)
    e_theta = np.column_stack((
        r_theta * cos_zeta,
        r_theta * sin_zeta,
        z_theta,
    ))
    e_zeta = np.column_stack((
        r_zeta * cos_zeta - radius * sin_zeta,
        r_zeta * sin_zeta + radius * cos_zeta,
        z_zeta,
    ))
    return e_theta, e_zeta


def validate_wout(path, n=10_000, seed=20260831):
    """Return sign-invariant angular errors for one VMEC equilibrium."""
    data = _read(path)
    required = ("lmns", "iotas", "bsupumnc", "bsupvmnc")
    missing = [name for name in required if name not in data]
    if missing:
        raise ValueError(f"{path} lacks required variables: {missing}")

    rng = np.random.default_rng(seed)
    ns = int(data["ns"])
    # Row zero of classic half-mesh quantities is padding.
    surface = rng.integers(1, ns, size=n)
    rho = np.sqrt((surface - 0.5) / (ns - 1))
    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    zeta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    e_theta, e_zeta = _geometry_basis(data, rho, theta, zeta)

    angle = theta[:, None] * data["xm"] - zeta[:, None] * data["xn"]
    cosine, sine = np.cos(angle), np.sin(angle)
    lmns = data["lmns"][surface]
    lambda_theta = np.einsum(
        "ij,ij->i", data["xm"] * lmns, cosine)
    lambda_zeta = np.einsum(
        "ij,ij->i", -data["xn"] * lmns, cosine)
    if "lmnc" in data:
        lmnc = data["lmnc"][surface]
        lambda_theta += np.einsum(
            "ij,ij->i", -data["xm"] * lmnc, sine)
        lambda_zeta += np.einsum(
            "ij,ij->i", data["xn"] * lmnc, sine)

    analytic = ((data["iotas"][surface] - lambda_zeta)[:, None] * e_theta
                + (1.0 + lambda_theta)[:, None] * e_zeta)

    nyq_angle = (theta[:, None] * data["xm_nyq"]
                 - zeta[:, None] * data["xn_nyq"])
    nyq_cosine, nyq_sine = np.cos(nyq_angle), np.sin(nyq_angle)
    bsupu = np.einsum(
        "ij,ij->i", data["bsupumnc"][surface], nyq_cosine)
    bsupv = np.einsum(
        "ij,ij->i", data["bsupvmnc"][surface], nyq_cosine)
    if "bsupumns" in data:
        bsupu += np.einsum(
            "ij,ij->i", data["bsupumns"][surface], nyq_sine)
        bsupv += np.einsum(
            "ij,ij->i", data["bsupvmns"][surface], nyq_sine)
    oracle = bsupu[:, None] * e_theta + bsupv[:, None] * e_zeta

    analytic /= np.linalg.norm(analytic, axis=1)[:, None]
    oracle /= np.linalg.norm(oracle, axis=1)[:, None]
    cosine_error = np.clip(
        np.abs(np.einsum("ij,ij->i", analytic, oracle)), 0.0, 1.0)
    degrees = np.degrees(np.arccos(cosine_error))
    return {
        "path": str(path),
        "samples": int(n),
        "median_degrees": float(np.median(degrees)),
        "p95_degrees": float(np.percentile(degrees, 95.0)),
        "p99_degrees": float(np.percentile(degrees, 99.0)),
        "max_degrees": float(np.max(degrees)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("wouts", nargs="+", type=Path)
    parser.add_argument("--samples", type=int, default=10_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = [validate_wout(path, args.samples) for path in args.wouts]
    text = json.dumps(results, indent=2)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
