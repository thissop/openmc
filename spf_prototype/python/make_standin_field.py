#!/usr/bin/env python
"""Generate the checked-in STAND-IN magnetic field map(s) for the 3D tier, so the
field-map path runs deterministically WITHOUT DESC. INJECT(helios): a real DESC
(or VMEC) QA equilibrium would write the same format; everything downstream
(C++ FieldMapField, python fieldmap, the plugin) is producer-agnostic.

Field-map format (portable, hand-read identically by C++ and Python):
  <name>.meta : text key=value (magic, grid sizes/bounds, nfp, period, endian)
  <name>.bin  : little-endian float64, Bx then By then Bz, each shape
                (nR, nphi, nZ) in C-order, UNIT-normalized at write time.

Grids:
  R[iR]   = R_min + iR*(R_max-R_min)/(nR-1)      (uniform, inclusive; clamp interp)
  Z[iZ]   = Z_min + iZ*(Z_max-Z_min)/(nZ-1)      (uniform, inclusive; clamp interp)
  phi[ip] = phi_min + ip*(2pi/nphi)              (FULL torus; phi=2pi wraps to ip=0)

NOTE: stored over the FULL torus (period=2pi), not one field period -- Cartesian B
components are not invariant under a field-period rotation, so a one-period
Cartesian map would give a wrong (rotated) B-hat outside period 1. (Same reasoning
as desc_to_fieldmap.py.)

The field is a pure closed-form function of position (no RNG, no timestamp), so
the .bin is reproducible byte-for-byte (tests re-run this and diff).

Analytic QA-like field (cylindrical (R,phi,Z); center R0c, scale a):
  rho = sqrt((R-R0c)^2 + Z^2);  theta = atan2(Z, R-R0c)
  B_phi = B0 * R0c / R                                   (1/R toroidal)
  B_pol = iota * (rho/a) * B0 * (1 + delta*cos(nfp*phi - theta))   (rot. transform,
                                                       2-field-period modulation)
  B_R = -B_pol*sin(theta);  B_Z = +B_pol*cos(theta)
  -> Cartesian, normalized. delta=iota=0 gives a purely toroidal field (the
     reduce-to-axisymmetric check).

Run:  python spf_prototype/python/make_standin_field.py   (writes data/standin_*.{meta,bin})
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DATADIR = REPO / "spf_prototype" / "data"

# Helios-class-ish domain (cm), aligned with reactor_model Part-B scale (R0~800).
# INJECT(helios): real plasma extent / equilibrium domain.
DEFAULTS = dict(R0c=800.0, a=300.0, B0=1.0,
                R_min=500.0, R_max=1100.0, Z_min=-300.0, Z_max=300.0,
                nR=25, nphi=32, nZ=25, nfp=2)  # nphi over FULL torus


def _field_cyl(R, phi, Z, R0c, a, B0, iota, delta, nfp):
    """Unit B-hat Cartesian components at cylindrical (R,phi,Z)."""
    rho = np.sqrt((R - R0c) ** 2 + Z ** 2)
    theta = np.arctan2(Z, R - R0c)
    Bphi = B0 * R0c / R
    Bpol = iota * (rho / a) * B0 * (1.0 + delta * np.cos(nfp * phi - theta))
    BR = -Bpol * np.sin(theta)
    BZ = Bpol * np.cos(theta)
    Bx = BR * np.cos(phi) - Bphi * np.sin(phi)
    By = BR * np.sin(phi) + Bphi * np.cos(phi)
    n = np.sqrt(Bx * Bx + By * By + BZ * BZ)
    return Bx / n, By / n, BZ / n


def write_map(name, iota, delta, **kw):
    p = dict(DEFAULTS, **kw)
    nR, nphi, nZ, nfp = p["nR"], p["nphi"], p["nZ"], p["nfp"]
    period = 2.0 * np.pi                                  # FULL torus (see header)
    R = np.linspace(p["R_min"], p["R_max"], nR)
    Z = np.linspace(p["Z_min"], p["Z_max"], nZ)
    phi = 0.0 + np.arange(nphi) * (period / nphi)        # phi_min = 0, full torus
    # C-order grids (nR, nphi, nZ)
    RR = R[:, None, None]; PP = phi[None, :, None]; ZZ = Z[None, None, :]
    Bx, By, Bz = _field_cyl(RR, PP, ZZ, p["R0c"], p["a"], p["B0"], iota, delta, nfp)
    Bx, By, Bz = np.broadcast_arrays(Bx, By, Bz)

    DATADIR.mkdir(parents=True, exist_ok=True)
    meta = (DATADIR / f"{name}.meta"); binf = (DATADIR / f"{name}.bin")
    meta.write_text(
        "magic spf_fieldmap_v1\n"
        f"ndim 3\nnR {nR}\nnphi {nphi}\nnZ {nZ}\n"
        f"R_min {p['R_min']!r}\nR_max {p['R_max']!r}\n"
        f"phi_min 0.0\nperiod {period!r}\n"
        f"Z_min {p['Z_min']!r}\nZ_max {p['Z_max']!r}\n"
        f"nfp {nfp}\nendian little\n")
    # little-endian float64, Bx then By then Bz, C-order
    with open(binf, "wb") as fh:
        for arr in (Bx, By, Bz):
            fh.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())
    # sanity: every stored vector is unit
    norms = np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    assert np.allclose(norms, 1.0, atol=1e-12)
    print(f"wrote {meta.name} + {binf.name}  ({nR}x{nphi}x{nZ}, nfp={nfp}, "
          f"iota={iota}, delta={delta}, {binf.stat().st_size} bytes)")
    return binf


def main():
    # 3D QA-like stand-in (rotational transform + 2-field-period modulation)
    write_map("standin_qa", iota=0.3, delta=0.15)
    # purely-toroidal stand-in (reduce-to-axisymmetric anchor: must match bmode=toroidal)
    write_map("standin_toroidal", iota=0.0, delta=0.0)


if __name__ == "__main__":
    sys.exit(main())
