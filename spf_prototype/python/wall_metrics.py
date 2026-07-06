#!/usr/bin/env python
"""Pure-numpy reduction of the first-wall directional mesh tally to a steering scalar.

Kept OpenMC-free ON PURPOSE so the reduction MATH is unit-testable off-cluster (the
aarch64/osx dev box has no OpenMC). The OpenMC-side reader that pulls the tally out of
a statepoint and calls this is run_ginsburg._wall_directional.

Why a near-source observable: the sweep's PRIMARY eta_source should be measured close
to the plasma with large solid angle (low MC variance), NOT at the deep coil (huge
variance without variance reduction). The inboard/outboard-midplane fast-flux contrast
below is polarization-sensitive -- mode A (perp) raises the inboard-midplane load, mode
B/C (parallel) lowers it -- and converges fast. See docs/EXPERIMENTAL_DESIGN.md.
"""
from __future__ import annotations

import numpy as np


def inboard_outboard_contrast(mean, sd, r_grid, z_grid, R0_cm, midplane_frac=0.34):
    """Near-source directional steering scalar from a cylindrical-mesh fast-flux tally.

    mean, sd : (nR, nphi, nZ) flux mean and 1-sigma (energy already collapsed to 1 bin).
    r_grid   : mesh R edges (len nR+1). z_grid: mesh Z edges (len nZ+1).
    R0_cm    : device major radius; R < R0 is 'inboard', R > R0 'outboard'.
    Returns (value, sd) with
        value = (I_inboard_mid - I_outboard_mid) / (I_inboard_mid + I_outboard_mid),
    the inboard/outboard imbalance of the first-wall fast flux over the midplane band
    |z| < midplane_frac * max|z_mid|. In [-1, 1]; 0 for an up/down- and in/out-symmetric
    (unpolarized-limit) load. The '(a-b)/(a+b)' form is self-normalizing, so the source
    strength and the total-rate factor cancel -- it measures DIRECTION only.
    """
    mean = np.asarray(mean, float)
    sd = np.asarray(sd, float)
    r_grid = np.asarray(r_grid, float)
    z_grid = np.asarray(z_grid, float)
    if mean.ndim != 3:
        raise ValueError(f"mean must be (nR,nphi,nZ); got {mean.shape}")
    r_mid = 0.5 * (r_grid[:-1] + r_grid[1:])
    z_mid = 0.5 * (z_grid[:-1] + z_grid[1:])
    inb = r_mid < R0_cm
    zmax = float(np.abs(z_mid).max()) or 1.0
    midz = np.abs(z_mid) < midplane_frac * zmax
    f = mean.sum(axis=1)                       # sum over phi -> (nR, nZ)
    var = (sd ** 2).sum(axis=1)                # variances add
    Iin = float(f[np.ix_(inb, midz)].sum())
    Iout = float(f[np.ix_(~inb, midz)].sum())
    Vin = float(var[np.ix_(inb, midz)].sum())
    Vout = float(var[np.ix_(~inb, midz)].sum())
    denom = Iin + Iout
    if denom <= 0:
        return float("nan"), float("nan")
    val = (Iin - Iout) / denom
    # error propagation for (a-b)/(a+b): d/da = 2b/(a+b)^2, d/db = -2a/(a+b)^2
    sd_val = float(np.sqrt((2 * Iout / denom ** 2) ** 2 * Vin
                           + (2 * Iin / denom ** 2) ** 2 * Vout))
    return float(val), sd_val
