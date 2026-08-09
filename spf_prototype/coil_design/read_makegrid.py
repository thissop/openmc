#!/usr/bin/env python
"""Read a MAKEGRID coils file (the format gil_to_makegrid writes / ParaStell consumes) back into
per-coil filament arrays and centroid layout. The inverse of gil_to_makegrid, so the placement
pipeline can run on a REAL coil arrangement (count + toroidal positions) instead of a synthetic ring.

MAKEGRID layout:
  periods N
  begin filament
  mirror NIL
    x y z I           <- filament vertices, current I repeated
    ...
    x y z 0  g Label   <- loop-closing row: I=0 and a trailing group id + label
  ... (repeat per coil)
  end
"""
import numpy as np


def read_makegrid(path):
    """Return (periods, coils) where coils is a list of (Nx4 float array [x y z I], label)."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    if not lines or not lines[0].startswith("periods"):
        raise ValueError(f"{path}: not a MAKEGRID file (missing 'periods' header)")
    periods = int(lines[0].split()[1])
    # body is between 'mirror NIL' (or 'begin filament') and 'end'
    try:
        start = next(i for i, ln in enumerate(lines) if ln.strip() == "mirror NIL") + 1
    except StopIteration:
        start = 3
    end = next(i for i, ln in enumerate(lines) if ln.strip() == "end")
    coils, cur, label = [], [], None
    for ln in lines[start:end]:
        tok = ln.split()
        if len(tok) < 4:
            continue
        cur.append([float(t) for t in tok[:4]])
        if len(tok) > 4:                                   # loop-closing row -> end of this coil
            label = " ".join(tok[4:])
            coils.append((np.array(cur), label))
            cur, label = [], None
    if not coils:
        raise ValueError(f"{path}: no coils parsed")
    return periods, coils


def coil_centroids(path):
    """Per-coil geometric centroid in cylindrical coords. Returns dict of arrays (length = n_coils):
    R [same units as file], phi [rad in [0,2pi)], Z, plus periods and n_coils. The filament's
    closing vertex duplicates the first, so it is dropped before averaging."""
    periods, coils = read_makegrid(path)
    R, PHI, Z = [], [], []
    for arr, _ in coils:
        xyz = arr[:-1, :3]                                 # drop duplicated closing vertex
        cx, cy, cz = xyz.mean(axis=0)
        R.append(float(np.hypot(cx, cy)))
        PHI.append(float(np.arctan2(cy, cx) % (2 * np.pi)))
        Z.append(float(cz))
    order = np.argsort(PHI)                                 # sort coils by toroidal angle
    return dict(R=np.array(R)[order], phi=np.array(PHI)[order], Z=np.array(Z)[order],
                periods=periods, n_coils=len(coils))


if __name__ == "__main__":
    import sys
    lay = coil_centroids(sys.argv[1])
    print(f"{lay['n_coils']} coils, periods={lay['periods']}")
    for k in range(lay["n_coils"]):
        print(f"  coil {k:2d}: R={lay['R'][k]:8.3f}  phi={np.degrees(lay['phi'][k]):7.2f} deg  "
              f"Z={lay['Z'][k]:8.3f}")
