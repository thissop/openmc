"""Pure-numpy geometry + field for a QUASR device, from its VMEC input namelist.

We avoid simsopt's compiled core entirely (its wheel SIGILLs on this CPU): the QUASR
per-device VMEC namelist  nml/<zpad7(ID)[:4]>/input.<zpad7(ID)>  carries the boundary
Fourier spectrum RBC(n,m)/ZBS(n,m) and NFP, which we evaluate directly:

    R(theta,phi) = sum RBC(n,m) cos(m*theta - n*NFP*phi)
    Z(theta,phi) = sum ZBS(n,m) sin(m*theta - n*NFP*phi)        (stellarator symmetry)

Interior flux surfaces are taken self-similar: the m>=1 (shaping) harmonics scale with
rho, the m=0 major-radius/axis terms stay put -- the same rho-scaling used for the toy
QA-low source loops.

Field model: on a flux surface the magnetic field is tangent and winds with the
rotational transform iota, so a field line is theta = theta0 + iota*phi and its 3-D
tangent is  dr/dphi|_theta + iota * dr/dtheta|_phi .  We use that direction as B-hat.
This is the equilibrium-implied field direction from the real boundary shape + the
device's iota (NOT a toy field); it uses the VMEC poloidal angle rather than the
straight-field-line angle, an O(shaping) approximation that is small for the gentle
(low-eps_eff) configs we select. Both the analytic and OpenMC sides call THIS function,
so the cross-check is valid regardless; B-hat only sets how realistic the pattern is.

Only depends on numpy + urllib.
"""
from __future__ import annotations

import re
import urllib.request
from pathlib import Path

import numpy as np

BASE = "https://quasr.flatironinstitute.org/nml"
NML_DIR = Path(__file__).resolve().parents[1] / "data" / "quasr" / "nml"
_RE = re.compile(r"(R[BC]+|Z[BS]+)\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*=\s*"
                 r"([-+0-9.eEdD]+)")


def fetch_nml(ID):
    z = f"{int(ID):07d}"
    NML_DIR.mkdir(parents=True, exist_ok=True)
    fn = NML_DIR / f"input.{z}"
    if not fn.exists():
        urllib.request.urlretrieve(f"{BASE}/{z[:4]}/input.{z}", fn)
    return fn


def parse_nml(fn):
    """-> dict(nfp:int, n:array, m:array, rbc:array, zbs:array)."""
    txt = Path(fn).read_text()
    mnfp = re.search(r"NFP\s*=\s*(\d+)", txt, re.I)
    nfp = int(mnfp.group(1)) if mnfp else 1
    rbc, zbs = {}, {}
    for kind, n, m, val in _RE.findall(txt):
        v = float(val.replace("D", "E").replace("d", "e"))
        n, m = int(n), int(m)
        if kind.upper().startswith("RB"):
            rbc[(n, m)] = v
        else:
            zbs[(n, m)] = v
    keys = sorted(set(rbc) | set(zbs))
    n = np.array([k[0] for k in keys]); m = np.array([k[1] for k in keys])
    R = np.array([rbc.get(k, 0.0) for k in keys])
    Z = np.array([zbs.get(k, 0.0) for k in keys])
    return dict(nfp=nfp, n=n, m=m, rbc=R, zbs=Z)


class Device:
    """Boundary geometry + flux-surface field for one QUASR equilibrium."""

    def __init__(self, spec, iota):
        self.nfp = spec["nfp"]
        self.n = spec["n"]; self.m = spec["m"]
        self.rbc = spec["rbc"]; self.zbs = spec["zbs"]
        self.iota = float(iota)
        self.R0 = float(self.rbc[(self.m == 0) & (self.n == 0)].sum())

    def _scale(self, rho):
        # m>=1 (and any n!=0) shaping scales with rho; the m=0,n=0 axis term fixed.
        s = np.where((self.m == 0) & (self.n == 0), 1.0, rho)
        return s

    def RZ(self, theta, phi, rho=1.0):
        """R,Z at (theta,phi) arrays (broadcast); flux surface rho."""
        th = np.asarray(theta)[..., None]; ph = np.asarray(phi)[..., None]
        a = self.m * th - self.n * self.nfp * ph
        s = self._scale(rho)
        R = np.sum(self.rbc * s * np.cos(a), axis=-1)
        Z = np.sum(self.zbs * s * np.sin(a), axis=-1)
        return R, Z

    def bhat(self, theta, phi, rho=1.0):
        """Unit flux-surface-tangent field at (theta,phi); shape (...,3) Cartesian."""
        th = np.asarray(theta)[..., None]; ph = np.asarray(phi)[..., None]
        a = self.m * th - self.n * self.nfp * ph
        s = self._scale(rho)
        nfp = self.nfp
        R = np.sum(self.rbc * s * np.cos(a), axis=-1)
        Rt = np.sum(self.rbc * s * (-self.m) * np.sin(a), axis=-1)
        Rp = np.sum(self.rbc * s * (self.n * nfp) * np.sin(a), axis=-1)
        Zt = np.sum(self.zbs * s * (self.m) * np.cos(a), axis=-1)
        Zp = np.sum(self.zbs * s * (-self.n * nfp) * np.cos(a), axis=-1)
        phi = np.asarray(phi); cp, sp = np.cos(phi), np.sin(phi)
        # dr/dtheta and dr/dphi in Cartesian
        xt = Rt * cp; yt = Rt * sp; zt = Zt
        xp = Rp * cp - R * sp; yp = Rp * sp + R * cp; zp = Zp
        tx = xp + self.iota * xt; ty = yp + self.iota * yt; tz = zp + self.iota * zt
        t = np.stack([tx, ty, tz], axis=-1)
        return t / np.linalg.norm(t, axis=-1, keepdims=True)


def load_device(ID, iota):
    return Device(parse_nml(fetch_nml(ID)), iota)


def eps_eff(dev, nphi=96, ntheta=48, gap=0.25, rho=1.0):
    """Geometric gentleness: max helical excursion / minor radius (e_geom) and a
    square-wall eps_eff (wall gap*a outside the bbox)."""
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z = dev.RZ(TH, PH, rho)          # (ntheta, nphi)
    a = 0.5 * (R.max() - R.min())
    rin, rout = R.min() - gap * a, R.max() + gap * a
    ztop = np.abs(Z).max() + gap * a
    exc = np.zeros(ntheta); eps = np.zeros(ntheta)
    for it in range(ntheta):
        Rl, Zl = R[it], Z[it]
        p, zc = Rl.mean(), Zl.mean()
        exc[it] = np.max(np.hypot(Rl - p, Zl - zc))
        so = min(p - rin, rout - p, ztop - zc, ztop + zc)
        eps[it] = exc[it] / max(so, 1e-9)
    return float(exc.max() / a), float(eps.max()), float(a)
