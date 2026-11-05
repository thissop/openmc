"""Flux-surface representations under head-to-head test.

Two ways to answer the same physical question a fusion source must answer at every
birth: given internal coordinates (s, theta, zeta), where in real space (R, Z, phi) is
this point, and what is the birth weight sqrt(g) there?

  * SpectralRep  -- Ethan's `StellaratorSource` style: keep the VMEC Fourier coefficients
    (rmnc, zmns[, rmns, zmnc]) and evaluate R, Z, and the metric Jacobian sqrt(g)
    ANALYTICALLY on the fly. Exact (to VMEC's own truncation); compact in memory; costs a
    Fourier sum per evaluation.
  * GridRep      -- Thaddaeus' `StellaratorSource` style: precompute R, Z, sqrt(g) on a
    (nrho x ntheta x nzeta) grid once and TRILINEARLY interpolate at query time. O(1) per
    evaluation and identical to the fluxmap the compiled source consumes, at the cost of
    interpolation error and a larger stored footprint.

This module is pure NumPy/SciPy (no compiled OpenMC needed) so the representation-level
accuracy/speed/memory comparison runs anywhere. The sampler-level comparison (actual birth
draws from each compiled source, KS-tested across the QUASR set) lives in compare_samplers.py
and does need both builds.

VMEC conventions used here: angle kernel is (xm*theta - xn*zeta), with xn already carrying the
nfp factor (as stored in the wout); zeta is the cylindrical toroidal angle phi; s = normalized
toroidal flux in [0, 1] on ns surfaces s_j = j/(ns-1). sqrt(g) is returned as |sqrt(g)| because
the source samples birth position proportional to the volume element |sqrt(g)| ds dtheta dzeta.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator


# --------------------------------------------------------------------------- wout I/O
def read_wout(path):
    """Read the fields we need from a VMEC wout: a wout_*.nc (classic NetCDF3 or NetCDF4) OR a
    wout_arrays_*.npz produced by gen_wout_arrays.py (canonical (ns, nmodes) coeff orientation).

    Returns a dict with xm, xn (mode numbers), nfp, s (radial grid), rmnc, zmns, and the
    optional non-stellarator-symmetric rmns, zmnc (None if absent), plus iotaf and lmns
    (the rotational transform and lambda stream function -- the two free functions from which
    B can be reconstructed on the fly; carried here so the B-hat backend can use them later).
    """
    if str(path).endswith(".npz"):
        d = np.load(str(path), allow_pickle=False)
        ns = int(d["ns"])
        return {
            "xm": d["xm"].astype(float), "xn": d["xn"].astype(float),
            "nfp": int(d["nfp"]), "ns": ns, "s": np.linspace(0.0, 1.0, ns),
            "rmnc": d["rmnc"], "zmns": d["zmns"],
            "rmns": d["rmns"] if "rmns" in d.files else None,
            "zmnc": d["zmnc"] if "zmnc" in d.files else None,
            "iotaf": d["iotaf"] if "iotaf" in d.files else None,
            "iotas": d["iotas"] if "iotas" in d.files else None,
            "lmns": d["lmns"] if "lmns" in d.files else None,
            "lmnc": d["lmnc"] if "lmnc" in d.files else None,
        }
    try:
        from scipy.io import netcdf_file
        f = netcdf_file(str(path), "r", mmap=False)
        get = lambda k: np.array(f.variables[k][()]) if k in f.variables else None
    except (ValueError, TypeError):
        import netCDF4
        f = netCDF4.Dataset(str(path), "r")
        get = lambda k: np.array(f.variables[k][:]) if k in f.variables else None

    ns = int(get("ns"))
    out = {
        "xm": get("xm").astype(float),
        "xn": get("xn").astype(float),
        "nfp": int(get("nfp")),
        "ns": ns,
        "s": np.linspace(0.0, 1.0, ns),
        "rmnc": get("rmnc"),            # (ns, nmodes)
        "zmns": get("zmns"),
        "rmns": get("rmns"),            # None for stellarator-symmetric equilibria
        "zmnc": get("zmnc"),
        "iotaf": get("iotaf"),          # (ns,) rotational transform, full grid
        "iotas": get("iotas"),          # (ns,) half grid with row-0 padding
        "lmns": get("lmns"),            # (ns,nmodes), half grid, row-0 padding
        "lmnc": get("lmnc"),            # asymmetric lambda cosine coefficients
    }
    return out


# --------------------------------------------------------------------------- spectral (exact)
class SpectralRep:
    """On-the-fly Fourier evaluation of R, Z, sqrt(g) -- the compact/exact representation."""

    def __init__(self, wout):
        self.xm = wout["xm"]
        self.xn = wout["xn"]
        self.s = wout["s"]
        self.asym = wout["rmns"] is not None and wout["zmnc"] is not None
        # cubic splines in s for every mode so we can evaluate at arbitrary s AND get d/ds
        self._rmnc = CubicSpline(self.s, wout["rmnc"], axis=0)
        self._zmns = CubicSpline(self.s, wout["zmns"], axis=0)
        self._rmns = CubicSpline(self.s, wout["rmns"], axis=0) if self.asym else None
        self._zmnc = CubicSpline(self.s, wout["zmnc"], axis=0) if self.asym else None
        self.has_field = wout.get("lmns") is not None
        self._lmns = self._lmnc = self._iota = None
        field_nbytes = 0
        if self.has_field:
            lmns = np.asarray(wout["lmns"])
            # Classic VMEC stores half-mesh quantities with row zero reserved
            # as padding. Keep the actual native half-mesh coordinates.
            if lmns.shape[0] == len(self.s):
                field_s = ((np.arange(1, len(self.s)) - 0.5)
                           / (len(self.s) - 1))
                lmns = lmns[1:]
                lmnc = (np.asarray(wout["lmnc"])[1:]
                        if wout.get("lmnc") is not None else None)
                if wout.get("iotas") is not None:
                    iota = np.asarray(wout["iotas"])[1:]
                else:
                    iotaf = np.asarray(wout["iotaf"])
                    iota = 0.5 * (iotaf[:-1] + iotaf[1:])
            else:
                field_s = np.linspace(0.0, 1.0, lmns.shape[0])
                lmnc = (np.asarray(wout["lmnc"])
                        if wout.get("lmnc") is not None else None)
                iota = np.asarray(wout["iotaf"])
            self._lmns = CubicSpline(field_s, lmns, axis=0)
            self._lmnc = (CubicSpline(field_s, lmnc, axis=0)
                          if lmnc is not None else None)
            self._iota = CubicSpline(field_s, iota)
            field_nbytes = lmns.nbytes + iota.nbytes
            if lmnc is not None:
                field_nbytes += lmnc.nbytes
        # memory footprint of what THIS representation must store (the coefficient arrays)
        self.nbytes = int(wout["rmnc"].nbytes + wout["zmns"].nbytes
                          + (wout["rmns"].nbytes + wout["zmnc"].nbytes if self.asym else 0)
                          + field_nbytes)

    def _kernels(self, theta, zeta):
        ang = np.multiply.outer(theta, self.xm) - np.multiply.outer(zeta, self.xn)  # (N, nmodes)
        return np.cos(ang), np.sin(ang)

    def eval(self, s, theta, zeta, want_sqrtg=True):
        """Evaluate at N points (arrays s, theta, zeta of equal length).

        Returns (R, Z) and, if want_sqrtg, |sqrt(g)| via the analytic (s, theta) Jacobian
        sqrt(g) = R * (R_s Z_theta - R_theta Z_s).
        """
        s = np.atleast_1d(np.asarray(s, float))
        theta = np.atleast_1d(np.asarray(theta, float))
        zeta = np.atleast_1d(np.asarray(zeta, float))
        cos, sin = self._kernels(theta, zeta)
        rc = self._rmnc(s); zs = self._zmns(s)                       # (N, nmodes)
        R = np.einsum("ij,ij->i", rc, cos)
        Z = np.einsum("ij,ij->i", zs, sin)
        if self.asym:
            R += np.einsum("ij,ij->i", self._rmns(s), sin)
            Z += np.einsum("ij,ij->i", self._zmnc(s), cos)
        if not want_sqrtg:
            return R, Z, None
        # analytic derivatives
        R_th = np.einsum("ij,ij->i", rc, -self.xm * sin)
        Z_th = np.einsum("ij,ij->i", zs,  self.xm * cos)
        rc_s = self._rmnc(s, 1); zs_s = self._zmns(s, 1)             # d/ds via spline
        R_s = np.einsum("ij,ij->i", rc_s, cos)
        Z_s = np.einsum("ij,ij->i", zs_s, sin)
        if self.asym:
            R_th += np.einsum("ij,ij->i", self._rmns(s),  self.xm * cos)
            Z_th += np.einsum("ij,ij->i", self._zmnc(s), -self.xm * sin)
            R_s += np.einsum("ij,ij->i", self._rmns(s, 1), sin)
            Z_s += np.einsum("ij,ij->i", self._zmnc(s, 1), cos)
        sqrtg = np.abs(R * (R_s * Z_th - R_th * Z_s))
        return R, Z, sqrtg

    def eval_bhat(self, s, theta, zeta):
        """Evaluate the analytic VMEC field-line unit direction in Cartesian coordinates."""
        if not self.has_field:
            raise ValueError("wout does not contain iota/lambda field data")
        s = np.atleast_1d(np.asarray(s, float))
        theta = np.atleast_1d(np.asarray(theta, float))
        zeta = np.atleast_1d(np.asarray(zeta, float))
        cosine, sine = self._kernels(theta, zeta)
        rc, zs = self._rmnc(s), self._zmns(s)
        R = np.einsum("ij,ij->i", rc, cosine)
        R_th = np.einsum("ij,ij->i", -self.xm * rc, sine)
        R_ze = np.einsum("ij,ij->i", self.xn * rc, sine)
        Z_th = np.einsum("ij,ij->i", self.xm * zs, cosine)
        Z_ze = np.einsum("ij,ij->i", -self.xn * zs, cosine)
        if self.asym:
            rs, zc = self._rmns(s), self._zmnc(s)
            R += np.einsum("ij,ij->i", rs, sine)
            R_th += np.einsum("ij,ij->i", self.xm * rs, cosine)
            R_ze += np.einsum("ij,ij->i", -self.xn * rs, cosine)
            Z_th += np.einsum("ij,ij->i", -self.xm * zc, sine)
            Z_ze += np.einsum("ij,ij->i", self.xn * zc, sine)

        ls = self._lmns(s)
        lambda_th = np.einsum("ij,ij->i", self.xm * ls, cosine)
        lambda_ze = np.einsum("ij,ij->i", -self.xn * ls, cosine)
        if self._lmnc is not None:
            lc = self._lmnc(s)
            lambda_th += np.einsum("ij,ij->i", -self.xm * lc, sine)
            lambda_ze += np.einsum("ij,ij->i", self.xn * lc, sine)
        b_th = self._iota(s) - lambda_ze
        b_ze = 1.0 + lambda_th

        cos_ze, sin_ze = np.cos(zeta), np.sin(zeta)
        radial = b_th * R_th + b_ze * R_ze
        toroidal = b_ze * R
        vertical = b_th * Z_th + b_ze * Z_ze
        field = np.column_stack((
            radial * cos_ze - toroidal * sin_ze,
            radial * sin_ze + toroidal * cos_ze,
            vertical,
        ))
        field /= np.linalg.norm(field, axis=1)[:, None]
        return field


# --------------------------------------------------------------------------- grid (precomputed)
class GridRep:
    """Precompute R, Z, sqrt(g) on a (nrho x ntheta x nzeta) grid, then trilinearly interpolate.

    This mirrors the fluxmap the compiled `StellaratorSource` consumes. `sqrtg_from` selects the
    truth used to fill the grid (a SpectralRep); accuracy is then measured as the interpolation
    error of this grid against that same SpectralRep at off-node query points.
    """

    def __init__(self, spectral: SpectralRep, nrho=12, ntheta=48, nzeta=144, nfp=1,
                 s_min=1e-4, s_nodes=None, theta_nodes=None, zeta_nodes=None):
        self.nrho, self.ntheta, self.nzeta = nrho, ntheta, nzeta
        # s_nodes / theta_nodes / zeta_nodes let a caller supply NON-uniform grids
        # (e.g. error-equidistributed ones); default is the uniform grid on each
        # axis. Angular nodes must be endpoint-free (no duplicate at the 2*pi/nfp
        # seam); the periodic padding below closes the wrap.
        if s_nodes is not None:
            s = np.asarray(s_nodes, float)
            self.nrho = nrho = len(s)
        else:
            s = np.linspace(s_min, 1.0, nrho)
        if theta_nodes is not None:
            th = np.asarray(theta_nodes, float)
            self.ntheta = ntheta = len(th)
        else:
            th = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
        if zeta_nodes is not None:
            ze = np.asarray(zeta_nodes, float)
            self.nzeta = nzeta = len(ze)
        else:
            ze = np.linspace(0.0, 2 * np.pi / nfp, nzeta, endpoint=False)
        S, TH, ZE = np.meshgrid(s, th, ze, indexing="ij")
        flat_s, flat_th, flat_ze = S.ravel(), TH.ravel(), ZE.ravel()
        n_points = flat_s.size
        R = np.empty(n_points); Z = np.empty(n_points); G = np.empty(n_points)
        B = np.empty((n_points, 3)) if spectral.has_field else None
        # High-mode QA equilibria otherwise materialize several multi-GiB
        # (point x mode) temporaries at production grid resolutions.
        for start in range(0, n_points, 10_000):
            stop = min(start + 10_000, n_points)
            slc = slice(start, stop)
            R[slc], Z[slc], G[slc] = spectral.eval(
                flat_s[slc], flat_th[slc], flat_ze[slc], want_sqrtg=True)
            if B is not None:
                B[slc] = spectral.eval_bhat(
                    flat_s[slc], flat_th[slc], flat_ze[slc])
        shp = (nrho, ntheta, nzeta)
        self._R = R.reshape(shp); self._Z = Z.reshape(shp); self._G = G.reshape(shp)
        self.has_field = spectral.has_field
        self._B = None
        if self.has_field:
            self._B = B.reshape((*shp, 3))

        # Periodically pad theta and zeta. Endpoint-free arrays passed directly
        # to RegularGridInterpolator extrapolate at the seam and bias the
        # comparison; the appended first slice makes the wrap explicit.
        th_pad = np.append(th, 2 * np.pi)
        ze_pad = np.append(ze, 2 * np.pi / nfp)
        self._grid = (s, th_pad, ze_pad)
        self._nfp = nfp

        def periodic_pad(values):
            values = np.concatenate((values, values[:, :1]), axis=1)
            return np.concatenate((values, values[:, :, :1]), axis=2)

        self._iR = RegularGridInterpolator(
            self._grid, periodic_pad(self._R), bounds_error=False,
            fill_value=None)
        self._iZ = RegularGridInterpolator(
            self._grid, periodic_pad(self._Z), bounds_error=False,
            fill_value=None)
        self._iG = RegularGridInterpolator(
            self._grid, periodic_pad(self._G), bounds_error=False,
            fill_value=None)
        self._iB = None
        if self.has_field:
            self._iB = [
                RegularGridInterpolator(
                    self._grid, periodic_pad(self._B[..., component]),
                    bounds_error=False, fill_value=None)
                for component in range(3)
            ]
        field_bytes = self._B.nbytes if self._B is not None else 0
        self.nbytes = int(
            self._R.nbytes + self._Z.nbytes + self._G.nbytes + field_bytes)

    def eval(self, s, theta, zeta, want_sqrtg=True):
        s = np.atleast_1d(np.asarray(s, float))
        theta = np.mod(np.atleast_1d(np.asarray(theta, float)), 2 * np.pi)
        zeta = np.mod(np.atleast_1d(np.asarray(zeta, float)), 2 * np.pi / self._nfp)
        pts = np.column_stack([np.clip(s, self._grid[0][0], self._grid[0][-1]), theta, zeta])
        R = self._iR(pts); Z = self._iZ(pts)
        G = self._iG(pts) if want_sqrtg else None
        return R, Z, G

    def eval_bhat(self, s, theta, zeta):
        """Trilinearly interpolate and renormalize the gridded field direction."""
        if not self.has_field:
            raise ValueError("grid was built without magnetic-field data")
        s = np.atleast_1d(np.asarray(s, float))
        theta = np.mod(np.atleast_1d(np.asarray(theta, float)), 2 * np.pi)
        zeta = np.mod(np.atleast_1d(np.asarray(zeta, float)),
                      2 * np.pi / self._nfp)
        pts = np.column_stack([
            np.clip(s, self._grid[0][0], self._grid[0][-1]), theta, zeta])
        field = np.column_stack([interp(pts) for interp in self._iB])
        field /= np.linalg.norm(field, axis=1)[:, None]
        return field
