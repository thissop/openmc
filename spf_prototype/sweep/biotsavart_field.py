#!/usr/bin/env python
"""Exact b_hat(x) anywhere in the plasma volume from filament coils via Biot-Savart.

This REPLACES the fitted / flux-surface-tangent field with the field the real QUASR
coils actually produce, removing the field-fidelity risk. The plasma volume is far
from the coil filaments, so the near-filament singularity is benign (documented in
field_audit). Two backends share one interface (.B(points) -> (Q,3) Tesla-up-to-a
-constant; .bhat(points) -> (Q,3) unit):

  CoilField        : Biot-Savart from a list of (polyline, current) filaments
                     (the QUASR coils, loaded by quasr_loader). PRIMARY.
  FluxTangentField : the existing quasr_geom flux-surface-tangent model, wrapped as
                     an audited FALLBACK so the sweep still runs if coils are absent.

Both write the `.so`-readable field map (spf_fieldmap_v1) via write_fieldmap, so the
trusted compiled source consumes either unchanged.

Units: everything is in the coils' native length unit (QUASR is metres). b_hat is
dimensionless and scale-invariant; the field-map GRID BOUNDS carry the unit and are
rescaled downstream by run_conformal.scale_fieldmap. The overall mu0*I/4pi constant
is dropped (it cancels in b_hat and in relative-current weighting).
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------- #
# Biot-Savart kernel (midpoint rule on closed polylines)
# --------------------------------------------------------------------------- #
def biot_savart_B(coils, points, eps=1e-12):
    """B (up to the common mu0/4pi) at `points` from filament `coils`.

    coils   : list of (poly (M,3), current). poly is a CLOSED filament sampled at M
              vertices (the wrap segment M-1 -> 0 is added here).
    points  : (Q,3) query points.
    Returns : (Q,3). Midpoint rule: for each segment (P_k -> P_{k+1}),
              dB = I * (dl x (r - r_mid)) / |r - r_mid|^3, summed over segments/coils.
    Second-order accurate; exact-in-the-limit for a straight segment far from it.
    """
    pts = np.asarray(points, float)
    B = np.zeros_like(pts)
    for poly, cur in coils:
        P = np.asarray(poly, float)
        Pn = np.roll(P, -1, axis=0)          # next vertex (closed)
        dl = Pn - P                           # (M,3)
        mid = 0.5 * (P + Pn)                  # (M,3)
        # r - r_mid : (Q,M,3)
        rvec = pts[:, None, :] - mid[None, :, :]
        r3 = (np.linalg.norm(rvec, axis=2) ** 3 + eps)[..., None]
        cross = np.cross(dl[None, :, :], rvec)     # (Q,M,3)  dl x (r-r_mid)
        B += cur * np.sum(cross / r3, axis=1)
    return B


def bhat_from_coils(coils, points, eps=1e-12):
    B = biot_savart_B(coils, points, eps=eps)
    n = np.linalg.norm(B, axis=1, keepdims=True)
    if np.any(n == 0):
        raise ValueError("Biot-Savart field vanishes at a query point")
    return B / n


# --------------------------------------------------------------------------- #
# Analytic reference (for the unit test): a single circular current loop
# --------------------------------------------------------------------------- #
def circular_loop(radius=1.0, current=1.0, n=400, center=(0, 0, 0), axis="z"):
    """Closed polyline for a circular loop of `radius` in the plane normal to `axis`.
    On-axis analytic field: |B|(z) = mu0 I R^2 / (2 (R^2+z^2)^{3/2}); direction +axis
    (right-hand rule). Returns [(poly, current)]."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = radius * np.cos(t); y = radius * np.sin(t); z = np.zeros_like(t)
    if axis == "z":
        poly = np.stack([x, y, z], axis=1)
    elif axis == "x":
        poly = np.stack([z, x, y], axis=1)
    elif axis == "y":
        poly = np.stack([y, z, x], axis=1)
    else:
        raise ValueError("axis must be x/y/z")
    return [(poly + np.asarray(center, float), current)]


def loop_Bz_onaxis(z, radius=1.0):
    """Analytic on-axis |B| profile SHAPE (drops mu0 I/4pi -> the 2pi*R^2 form
    matching biot_savart_B's dropped constant). Returns the value biot_savart_B
    gives on axis, so tests compare kernel-to-analytic directly:
        B_z(z) = 2*pi*R^2 / (R^2 + z^2)^{3/2}
    (mu0 I/4pi dropped in the kernel => the mu0 I/2 * R^2/(...)^{3/2} analytic
    becomes 2*pi*R^2/(...)^{3/2}.)"""
    return 2 * np.pi * radius ** 2 / (radius ** 2 + np.asarray(z, float) ** 2) ** 1.5


# --------------------------------------------------------------------------- #
# Backends
# --------------------------------------------------------------------------- #
class CoilField:
    """Biot-Savart field from QUASR coil filaments (the exact field)."""

    def __init__(self, coils, meta=None):
        self.coils = coils          # list of (poly (M,3), current)
        self.meta = meta or {}

    def B(self, points):
        return biot_savart_B(self.coils, points)

    def bhat(self, points):
        return bhat_from_coils(self.coils, points)


class FluxTangentField:
    """Fallback: quasr_geom flux-surface-tangent b_hat. Only bhat() is meaningful
    (it is a DIRECTION model, not a magnitude), so B() raises."""

    def __init__(self, device):
        self.device = device

    def bhat_flux(self, theta, phi, rho=1.0):
        return self.device.bhat(theta, phi, rho)

    def B(self, points):
        raise NotImplementedError("FluxTangentField is direction-only; use bhat_flux")


# --------------------------------------------------------------------------- #
# spf_fieldmap_v1 writer (drives the compiled source unchanged)
# --------------------------------------------------------------------------- #
def write_fieldmap(field, bounds, dims, out_stem, nfp, source="biotsavart:coils"):
    """Sample b_hat on a cylindrical grid and write <out_stem>.meta + .bin in the
    spf_fieldmap_v1 format that fieldmap.FieldMapField (and the C++ source) read.

    field  : object with .bhat(points (Q,3)) -> (Q,3)
    bounds : dict(R_min,R_max,Z_min,Z_max) in the field's length unit; phi is full 2pi.
    dims   : dict(nR,nphi,nZ). period is 2*pi (lab Cartesian components are NOT
             nfp-periodic, so the map must span the full torus).
    """
    nR, nphi, nZ = dims["nR"], dims["nphi"], dims["nZ"]
    R = np.linspace(bounds["R_min"], bounds["R_max"], nR)
    phi = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
    Z = np.linspace(bounds["Z_min"], bounds["Z_max"], nZ)
    # grid in the reader's index order (nR, nphi, nZ), C-order flatten
    Rg, Pg, Zg = np.meshgrid(R, phi, Z, indexing="ij")
    X = (Rg * np.cos(Pg)).ravel()
    Y = (Rg * np.sin(Pg)).ravel()
    pts = np.stack([X, Y, Zg.ravel()], axis=1)
    bh = field.bhat(pts)                              # (Q,3), unit
    bx = bh[:, 0].astype("<f8"); by = bh[:, 1].astype("<f8"); bz = bh[:, 2].astype("<f8")
    out_stem = str(out_stem)
    with open(out_stem + ".bin", "wb") as f:
        f.write(bx.tobytes()); f.write(by.tobytes()); f.write(bz.tobytes())
    meta = [
        "magic spf_fieldmap_v1", "ndim 3",
        f"nR {nR}", f"nphi {nphi}", f"nZ {nZ}",
        f"R_min {bounds['R_min']!r}", f"R_max {bounds['R_max']!r}",
        "phi_min 0.0", f"period {2 * np.pi!r}",
        f"Z_min {bounds['Z_min']!r}", f"Z_max {bounds['Z_max']!r}",
        f"nfp {nfp}", "endian little", f"source {source}",
    ]
    Path(out_stem + ".meta").write_text("\n".join(meta) + "\n")
    return out_stem
