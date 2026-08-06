"""Physically-grounded plasma major/minor radius from the boundary geometry.

Motivation (patch-worth correctness): the worth metric bins shield voxels by a minor-radius-like
coordinate rho = hypot(R - R0, Z) and selects the shield band rho in [a_minor + d_in, a_minor +
d_out]. Both R0 and a_minor were previously HEURISTICS:
  * R0    = importance-weighted mean cylindrical radius of the nonzero adjoint voxels
            (biased toward wherever the source/importance happened to live: gave 1298 vs 1337 cm
            for the coil-20 contributon vs the multi-coil map -- inconsistent);
  * a_minor = 5th percentile of rho over psi>0 voxels (a rho-of-the-source quantile, not a plasma
            dimension at all).
Neither is a property of the plasma. The magnetic-axis major radius R0 and the LCFS minor radius
a_minor are geometric quantities of the boundary surface (R(theta,phi), Z(theta,phi)) -- so we
compute them there, once, and feed both to the worth binning. The band OFFSETS (d_in, d_out) are
the shield layer's radial position from the actual radial build (FW+mult+breeder+back_wall to
+shield), so the band tracks the real geometry rather than a percentile.

All functions are unit-agnostic; pass scale=100.0 to convert VMEC metres -> cm to match the mesh.
"""
from __future__ import annotations

import numpy as np

# Corrected 8-layer QH radial build (cm) -- MUST match build_step1.py / build_patch.py.
# The shield layer sits beyond the LCFS by [FW+mult+breeder+back_wall, +shield].
RADIAL_BUILD_CM = dict(first_wall=3.2, multiplier=2.0, breeder=50.0, back_wall=4.0,
                       shield=40.0, gap=2.0, vacuum_vessel=25.0, thermal_shield=3.0)


def shield_band_offsets(build=RADIAL_BUILD_CM):
    """(d_in, d_out) cm: radial distance from the LCFS to the shield layer's inner/outer face."""
    d_in = build["first_wall"] + build["multiplier"] + build["breeder"] + build["back_wall"]
    d_out = d_in + build["shield"]
    return float(d_in), float(d_out)


def breeder_band_offsets(build=RADIAL_BUILD_CM):
    """(d_in, d_out) cm: radial distance from the LCFS to the BREEDER layer's inner/outer face.
    The fixed-envelope trade thins the breeder as it thickens the shield, so the signed worth
    must debit the breeder band (shallower -> higher flux) against the shield band's benefit."""
    d_in = build["first_wall"] + build["multiplier"]
    d_out = d_in + build["breeder"]
    return float(d_in), float(d_out)


def _boundary(arr):
    """Return the outermost flux surface as a 2-D (n_theta, n_phi) array. Accepts a 2-D boundary
    already, or a 3-D (n_rho, n_theta, n_phi) fluxmap volume (LCFS = last rho index)."""
    a = np.asarray(arr, float)
    if a.ndim == 3:
        return a[-1]
    if a.ndim == 2:
        return a
    raise ValueError(f"expected 2-D boundary or 3-D volume, got shape {a.shape}")


def axis_minor(R, Z, scale=1.0):
    """Major/minor radius from a boundary surface R(theta,phi), Z(theta,phi).

    The per-phi geometric axis is the poloidal centroid (mean over theta); R0 is that axis' major
    radius averaged over phi. The minor radius is the RMS distance of the boundary from its per-phi
    axis (exact 'a' for a circular cross-section; a fair central value for a shaped one). a_max is
    the largest boundary excursion, reported for the band's outer reach.

    Returns dict: R0, a_minor (RMS), a_max, R0_phi (per-phi axis major radius), Z0_phi. All *scale.
    theta is taken to be the FIRST axis of R/Z (surface npz: R (n_theta, n_phi); fluxmap LCFS the
    same after _boundary())."""
    R = _boundary(R) * scale
    Z = _boundary(Z) * scale
    R0_phi = R.mean(axis=0)                       # per-phi geometric axis (centroid over theta)
    Z0_phi = Z.mean(axis=0)
    dR = R - R0_phi[None, :]
    dZ = Z - Z0_phi[None, :]
    r2 = dR * dR + dZ * dZ
    a_rms = float(np.sqrt(r2.mean()))
    a_max = float(np.sqrt(r2.max()))
    return dict(R0=float(R0_phi.mean()), a_minor=a_rms, a_max=a_max,
                R0_phi=R0_phi, Z0_phi=Z0_phi)


def axis_minor_from_surface(npz, scale=100.0):
    """(R0, a_minor, a_max) from a *_surface.npz (keys R, Z in metres). scale=100 -> cm."""
    d = np.load(npz) if isinstance(npz, str) else npz
    g = axis_minor(d["R"], d["Z"], scale=scale)
    return g["R0"], g["a_minor"], g["a_max"]


def axis_minor_from_fluxmap(fluxmap, scale=100.0):
    """(R0, a_minor, a_max) from a VMEC fluxmap npz (keys R, Z; boundary = outer flux surface).
    Same geometry the worth attribution already loads -> R0/a_minor become self-consistent with
    the source rasterization instead of a separate voxel heuristic."""
    d = np.load(fluxmap) if isinstance(fluxmap, str) else fluxmap
    g = axis_minor(d["R"], d["Z"], scale=scale)
    return g["R0"], g["a_minor"], g["a_max"]
