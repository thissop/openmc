"""Python mirror of src/spf_fieldmap.hpp (FieldMapField). Reads the same
<stem>.meta + <stem>.bin and interpolates IDENTICALLY (clamped R/Z, periodic phi,
fixed term order R->phi->Z, normalize as (1/|B|)*B). Parity vs C++ is checked in
tests/test_field_parity.py. Kept out of spf_mirror.py so the tier-1 sampler tests
stay numpy-free; loaded lazily by spf_mirror.make_field(bmode='fieldmap')."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


class FieldMapField:
    def __init__(self, stem: str):
        meta = {}
        for line in Path(stem + ".meta").read_text().splitlines():
            if line.strip():
                k, v = line.split()
                meta[k] = v
        if meta.get("magic") != "spf_fieldmap_v1":
            raise ValueError(f"FieldMapField: bad/missing magic in {stem}.meta")
        self.nR = int(meta["nR"]); self.nphi = int(meta["nphi"]); self.nZ = int(meta["nZ"])
        self.Rmin = float(meta["R_min"]); self.Zmin = float(meta["Z_min"])
        self.phimin = float(meta["phi_min"]); self.period = float(meta["period"])
        self.dR = (float(meta["R_max"]) - self.Rmin) / (self.nR - 1)
        self.dZ = (float(meta["Z_max"]) - self.Zmin) / (self.nZ - 1)
        self.dphi = self.period / self.nphi
        n = self.nR * self.nphi * self.nZ
        raw = np.fromfile(stem + ".bin", dtype="<f8")
        if raw.size != 3 * n:
            raise ValueError(f"FieldMapField: {stem}.bin size {raw.size} != 3*{n}")
        shape = (self.nR, self.nphi, self.nZ)
        self.bx = raw[0:n].reshape(shape)
        self.by = raw[n:2 * n].reshape(shape)
        self.bz = raw[2 * n:3 * n].reshape(shape)

    @staticmethod
    def _clamp(i, n):
        if i < 0:
            return 0
        if i > n - 2:
            return n - 2
        return i

    @staticmethod
    def _tri(a, iR, iR1, ip, ip1, iZ, iZ1, fR, fp, fZ):
        c00 = float(a[iR, ip, iZ]) * (1 - fZ) + float(a[iR, ip, iZ1]) * fZ
        c01 = float(a[iR, ip1, iZ]) * (1 - fZ) + float(a[iR, ip1, iZ1]) * fZ
        c10 = float(a[iR1, ip, iZ]) * (1 - fZ) + float(a[iR1, ip, iZ1]) * fZ
        c11 = float(a[iR1, ip1, iZ]) * (1 - fZ) + float(a[iR1, ip1, iZ1]) * fZ
        c0 = c00 * (1 - fp) + c01 * fp
        c1 = c10 * (1 - fp) + c11 * fp
        return c0 * (1 - fR) + c1 * fR

    def bhat(self, r):
        x, y, z = r
        R = math.sqrt(x * x + y * y)
        phi = math.atan2(y, x)
        Z = z
        tR = (R - self.Rmin) / self.dR
        iR = self._clamp(int(math.floor(tR)), self.nR); fR = tR - iR
        tZ = (Z - self.Zmin) / self.dZ
        iZ = self._clamp(int(math.floor(tZ)), self.nZ); fZ = tZ - iZ
        pp = math.fmod(phi - self.phimin, self.period)
        if pp < 0.0:
            pp += self.period
        tphi = pp / self.dphi
        ip = int(math.floor(tphi)); fphi = tphi - ip
        ip = ip % self.nphi
        if ip < 0:
            ip += self.nphi
        ip1 = (ip + 1) % self.nphi
        iR1 = iR + 1; iZ1 = iZ + 1
        bx = self._tri(self.bx, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ)
        by = self._tri(self.by, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ)
        bz = self._tri(self.bz, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ)
        n = math.sqrt(bx * bx + by * by + bz * bz)
        inv = 1.0 / n  # match spf::normalized = (1/|B|)*B
        return (bx * inv, by * inv, bz * inv)
