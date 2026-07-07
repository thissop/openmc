#!/usr/bin/env python
"""Conformal StellaratorSource position sampler (Python reference for the C++ core class).

Birth positions are drawn from  p(rho, theta, zeta) proportional to  S(rho) * sqrt(g)(rho,theta,zeta)
via a REJECTION-FREE marginal -> conditional CDF cascade  p(rho) p(zeta|rho) p(theta|rho,zeta),
so every birth weight is exactly 1. Consumes the spf_fluxmap_v1 .npz written by quasr_fluxmap.py:
    rho[nr], theta[nt], zeta[nz], sqrtg[nr,nt,nz], R/Z/phi[nr,nt,nz], BR/Bphi/BZ[nr,nt,nz], nfp.

The field b_hat is read from the SAME grid, so position and field are self-consistent on one
equilibrium (the unification PLAN_B calls for). Direction sampling (the SPF P2 mixture about b_hat)
is layered on top by the caller and is validated separately.
"""
from __future__ import annotations

import numpy as np


def _wrap_index(x, grid, period):
    """Return (i0, i1, frac) for periodic-linear interpolation of value x on a uniform `grid`."""
    n = grid.size
    d = period / n
    t = (x - grid[0]) / d
    i0 = np.floor(t).astype(int) % n
    frac = t - np.floor(t)
    i1 = (i0 + 1) % n
    return i0, i1, frac


class ConformalStellaratorSampler:
    def __init__(self, fluxmap, S_of_rho=None):
        f = fluxmap
        self.rho = np.asarray(f["rho"], float)
        self.theta = np.asarray(f["theta"], float)
        self.zeta = np.asarray(f["zeta"], float)
        self.sqrtg = np.asarray(f["sqrtg"], float)                 # (nr,nt,nz), >= 0
        self.R = np.asarray(f["R"], float); self.Z = np.asarray(f["Z"], float)
        self.phi = np.asarray(f["phi"], float)
        self.BR = np.asarray(f["BR"], float); self.Bphi = np.asarray(f["Bphi"], float)
        self.BZ = np.asarray(f["BZ"], float)
        self.nfp = int(f["nfp"])
        self.nr, self.nt, self.nz = self.sqrtg.shape
        assert self.sqrtg.min() >= -1e-12, "sqrt(g) must be non-negative (store |sqrt(g)|)"

        if S_of_rho is None:
            S = np.clip(1.0 - self.rho ** 2, 0.0, None)            # parabolic flux function
        else:
            S = np.asarray([float(S_of_rho(r)) for r in self.rho], float)
        self.S = S

        # per-cell emission weight (uniform (rho,theta,zeta) spacing -> constants cancel in CDFs)
        w = self.S[:, None, None] * self.sqrtg                     # (nr,nt,nz)
        self._w = w
        m_rho = w.sum(axis=(1, 2))                                 # p(rho)
        self._cdf_rho = np.concatenate([[0.0], np.cumsum(m_rho)]); self._cdf_rho /= self._cdf_rho[-1]
        cz = w.sum(axis=1)                                         # p(zeta|rho): (nr,nz)
        self._cdf_zeta = np.concatenate([np.zeros((self.nr, 1)), np.cumsum(cz, axis=1)], axis=1)
        dz = self._cdf_zeta[:, -1:]                                # 0 for rho slices with S(rho)=0
        self._cdf_zeta = self._cdf_zeta / np.where(dz > 0, dz, 1.0)
        ct = np.concatenate([np.zeros((self.nr, 1, self.nz)), np.cumsum(w, axis=1)], axis=1)
        dt = ct[:, -1:, :]                                         # p(theta|rho,zeta): (nr,nt+1,nz)
        self._cdf_theta = ct / np.where(dt > 0, dt, 1.0)

    # -- discrete-cell sampling of (i,j,k) proportional to the cell weight, via the cascade --
    def _sample_cells(self, n, rng):
        u = rng.random(n)
        i = np.searchsorted(self._cdf_rho, u, side="right") - 1
        i = np.clip(i, 0, self.nr - 1)
        uz = rng.random(n)
        k = np.array([np.searchsorted(self._cdf_zeta[ii], uu, side="right") - 1
                      for ii, uu in zip(i, uz)])
        k = np.clip(k, 0, self.nz - 1)
        ut = rng.random(n)
        j = np.array([np.searchsorted(self._cdf_theta[ii, :, kk], uu, side="right") - 1
                      for ii, kk, uu in zip(i, k, ut)])
        j = np.clip(j, 0, self.nt - 1)
        return i, j, k

    def sample(self, n, rng=None, jitter=True):
        """Return dict with cartesian positions (x,y,z) [same units as R,Z] and unit b_hat.
        Birth weight is exactly 1 for every particle."""
        rng = rng or np.random.default_rng()
        i, j, k = self._sample_cells(n, rng)
        # continuous (rho,theta,zeta): jitter uniformly inside the selected cell (piecewise-const
        # density per cell -> unbiased to the grid discretization, which -> 0 under refinement)
        drho = self.rho[1] - self.rho[0]
        dth = self.theta[1] - self.theta[0]
        dze = self.zeta[1] - self.zeta[0]
        if jitter:
            rr = self.rho[i] + (rng.random(n) - 0.5) * drho
            tt = self.theta[j] + (rng.random(n) - 0.5) * dth
            zz = self.zeta[k] + (rng.random(n) - 0.5) * dze
        else:
            rr, tt, zz = self.rho[i], self.theta[j], self.zeta[k]
        # trilinear interp of geometry+field at (rr,tt,zz): rho non-periodic, theta/zeta periodic
        ir = np.clip(np.searchsorted(self.rho, rr) - 1, 0, self.nr - 2)
        fr = (rr - self.rho[ir]) / drho
        jt0, jt1, ft = _wrap_index(tt, self.theta, 2 * np.pi)
        kz0, kz1, fz = _wrap_index(zz, self.zeta, 2 * np.pi)

        def interp(A):
            c00 = A[ir, jt0, kz0] * (1 - ft) + A[ir, jt1, kz0] * ft
            c01 = A[ir, jt0, kz1] * (1 - ft) + A[ir, jt1, kz1] * ft
            c0 = c00 * (1 - fz) + c01 * fz
            d00 = A[ir + 1, jt0, kz0] * (1 - ft) + A[ir + 1, jt1, kz0] * ft
            d01 = A[ir + 1, jt0, kz1] * (1 - ft) + A[ir + 1, jt1, kz1] * ft
            d0 = d00 * (1 - fz) + d01 * fz
            return c0 * (1 - fr) + d0 * fr

        R = interp(self.R); Z = interp(self.Z); PH = interp(self.phi)
        x, y = R * np.cos(PH), R * np.sin(PH)
        br, bp, bz = interp(self.BR), interp(self.Bphi), interp(self.BZ)
        bx = br * np.cos(PH) - bp * np.sin(PH)
        by = br * np.sin(PH) + bp * np.cos(PH)
        bn = np.sqrt(bx ** 2 + by ** 2 + bz ** 2) + 1e-300
        return dict(x=x, y=y, z=Z, R=R, phi=PH, rho=rr, theta=tt, zeta=zz,
                    bhat=np.stack([bx / bn, by / bn, bz / bn], axis=-1),
                    weight=np.ones(n))


def circular_torus_fluxmap(R0=6.0, a=2.0, nr=16, nt=64, nz=96, nfp=1):
    """Analytic circular-torus fluxmap with a KNOWN sqrt(g) for validating the sampler.
    R = R0 + a*rho*cos(theta), Z = a*rho*sin(theta), phi = zeta.  Poloidal Jacobian = a^2*rho,
    so sqrt(g) = R * a^2 * rho  (exact). Toroidal B_phi ~ 1/R (toroidal field)."""
    rho = np.linspace(1.0 / nr, 1.0, nr)
    theta = np.linspace(0.0, 2 * np.pi, nt, endpoint=False)
    zeta = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    RHO, TH, ZE = np.meshgrid(rho, theta, zeta, indexing="ij")
    R = R0 + a * RHO * np.cos(TH); Z = a * RHO * np.sin(TH); PHI = ZE
    sqrtg = R * (a ** 2) * RHO
    Bphi = 1.0 / R
    return dict(rho=rho, theta=theta, zeta=zeta, sqrtg=sqrtg, R=R, Z=Z, phi=PHI,
                BR=np.zeros_like(R), Bphi=Bphi, BZ=np.zeros_like(R), nfp=nfp), (R0, a)
