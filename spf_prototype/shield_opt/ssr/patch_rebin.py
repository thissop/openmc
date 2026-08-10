#!/usr/bin/env python
"""Offline re-binning of a frozen surface source into an arbitrary (toroidal, poloidal) patch mesh.

This is the "play with resolution for free" layer of the decoupled pipeline: once the master run has
written the first-wall (or shield-exit) surface source, re-bin it into any patch resolution without
re-running transport. Handles the two things Peterson flagged:
  * patch AREA normalization -- equal-angular patches cover unequal physical area on a twisted wall,
    so a load map must divide by true patch area (from the VMEC surface metric), else the map is a
    geometry artifact. area_normalize=True does this.
  * MAX-in-patch reduction -- for the optimizer cost, take the peak load inside a patch, not the mean,
    so a localized hotspot inside a coarse patch is not averaged away (reduce="max").

Reads OpenMC's surface_source.h5 ('source_bank' compound dataset: r, u, E, wgt).
"""
import numpy as np
import h5py


def load_surface_source(path):
    """Return dict r (N,3 cm), u (N,3), E (N eV), wgt (N)."""
    with h5py.File(path, "r") as f:
        sb = f["source_bank"][...]
    return dict(r=np.asarray(sb["r"], float), u=np.asarray(sb["u"], float),
                E=np.asarray(sb["E"], float), wgt=np.asarray(sb["wgt"], float))


def source_to_angles(r, R0, nfp=1):
    """(phi_folded, theta) for wall points r. phi in [0, 2pi/nfp); theta = atan2(Z, R-R0) in [-pi,pi)."""
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    phi = np.mod(np.arctan2(y, x), 2 * np.pi / nfp)
    R = np.hypot(x, y)
    theta = np.arctan2(z, R - R0)
    return phi, theta


def patch_edges(ntor, npol, nfp=1):
    return np.linspace(0.0, 2 * np.pi / nfp, ntor + 1), np.linspace(-np.pi, np.pi, npol + 1)


def bin_patches(phi, theta, val, tor_edges, pol_edges, reduce="sum"):
    """Bin per-particle `val` into (ntor, npol) patches. reduce in {sum, max, mean, count}."""
    it = np.clip(np.digitize(phi, tor_edges) - 1, 0, len(tor_edges) - 2)
    jp = np.clip(np.digitize(theta, pol_edges) - 1, 0, len(pol_edges) - 2)
    ntor, npol = len(tor_edges) - 1, len(pol_edges) - 1
    out = np.zeros((ntor, npol)); cnt = np.zeros((ntor, npol))
    if reduce == "max":
        out[:] = -np.inf
    for i, j, v in zip(it, jp, val):
        cnt[i, j] += 1
        if reduce == "max":
            out[i, j] = max(out[i, j], v)
        else:
            out[i, j] += v
    if reduce == "max":
        out[~np.isfinite(out)] = 0.0
    elif reduce == "mean":
        out = np.where(cnt > 0, out / np.maximum(cnt, 1), 0.0)
    elif reduce == "count":
        out = cnt
    return out


def _vmec_boundary(wout_path):
    import netCDF4
    d = netCDF4.Dataset(wout_path)
    xm = d.variables["xm"][:]; xn = d.variables["xn"][:]
    rmnc = d.variables["rmnc"][-1, :]; zmns = d.variables["zmns"][-1, :]
    nfp = int(d.variables["nfp"][:])
    R0 = float(rmnc[np.argmax((xm == 0) & (xn == 0))])

    def surf(theta, phi):
        ang = np.multiply.outer(xm, theta)[..., None] - np.multiply.outer(xn, phi)[:, None, :]
        R = np.tensordot(rmnc, np.cos(ang), axes=(0, 0))
        Z = np.tensordot(zmns, np.sin(ang), axes=(0, 0))
        return R, Z
    return surf, R0, nfp


def patch_areas(tor_edges, pol_edges, wout_path, sub=8):
    """Physical area (same length units as VMEC, i.e. m^2) of each (tor,pol) patch on the wall,
    from the VMEC boundary surface metric |dr/dtheta x dr/dphi|. Addresses the equal-angle != equal-area
    confounder. sub = fine sub-samples per patch for the area quadrature."""
    surf, R0, _ = _vmec_boundary(wout_path)
    ntor, npol = len(tor_edges) - 1, len(pol_edges) - 1
    A = np.zeros((ntor, npol))
    for i in range(ntor):
        ph = np.linspace(tor_edges[i], tor_edges[i + 1], sub + 1)
        for j in range(npol):
            th = np.linspace(pol_edges[j], pol_edges[j + 1], sub + 1)
            R, Z = surf(th, ph)                          # (sub+1, sub+1)
            X = R * np.cos(ph)[None, :]; Y = R * np.sin(ph)[None, :]
            # triangulate the sub-grid, sum triangle areas
            P = np.stack([X, Y, Z], axis=-1)             # (nth, nph, 3)
            a = 0.0
            for ii in range(sub):
                for jj in range(sub):
                    p00, p10, p01, p11 = P[ii, jj], P[ii + 1, jj], P[ii, jj + 1], P[ii + 1, jj + 1]
                    a += 0.5 * np.linalg.norm(np.cross(p10 - p00, p01 - p00))
                    a += 0.5 * np.linalg.norm(np.cross(p10 - p11, p01 - p11))
            A[i, j] = a
    return A


def rebin_load(source, ntor, npol, nfp=4, R0=None, weight="energy", reduce="sum",
               wout_path=None, area_normalize=False):
    """Frozen surface source -> (ntor,npol) load map.
    weight: 'energy' (wgt*E, ~ wall loading), 'flux' (wgt), 'fast' (wgt where E>0.1 MeV).
    reduce: 'sum' (patch total) or 'max' (peak-in-patch for the optimizer guardrail).
    area_normalize: divide by VMEC patch area (needs wout_path)."""
    r = source["r"]
    if R0 is None:
        if wout_path is None:
            R0 = float(np.median(np.hypot(r[:, 0], r[:, 1])))
        else:
            R0 = _vmec_boundary(wout_path)[1] * 100.0     # VMEC in m -> cm to match source (cm)
    phi, theta = source_to_angles(r, R0, nfp=nfp)
    w = source["wgt"]
    if weight == "energy":
        val = w * source["E"]
    elif weight == "fast":
        val = w * (source["E"] > 0.1e6)
    else:
        val = w
    m = bin_patches(phi, theta, val, *patch_edges(ntor, npol, nfp), reduce=reduce)
    if area_normalize:
        assert wout_path is not None, "area_normalize needs wout_path"
        A = patch_areas(*patch_edges(ntor, npol, nfp), wout_path)
        m = np.where(A > 0, m / A, 0.0)
    return m, R0
