#!/usr/bin/env python
"""Load a QUASR device record: coils (filament geometry + currents), plasma
boundary, and metadata (nfp, iota, aspect, symmetry class).

QUASR serves each device as a simsopt-serialized JSON at
  https://quasr.flatironinstitute.org/simsopt_serials/<ID[:4]>/serial<ID>.json
which we parse WITHOUT importing simsopt (its compiled core SIGILLs on some CPUs;
the pure-JSON parse also runs anywhere). The serial is a graph of:
  CurveXYZFourier (base modular coils, order N, 3*(2N+1) dofs, simsopt layout
                   per coord [c0, s1, c1, s2, c2, ... sN, cN]),
  RotatedCurve{curve, phi, flip}  (field-period rotation + stellarator-symmetry),
  Current / ScaledCurrent{current_to_scale, scale}, and Coil{curve, current}.
The full coil set is the list of Coil objects (base * nfp * stellsym).

We VALIDATE the parse two ways: (1) summed filament length vs the catalogue
`total_coil_length`; (2) downstream, field_audit checks B.n ~ 0 on the LCFS (the
coils make the boundary a flux surface) -- that catches any rotation/current error.

Boundary + metadata reuse the existing pure-numpy quasr_geom loader.
"""
from __future__ import annotations

import csv
import gzip
import json
import urllib.request
from pathlib import Path

import numpy as np

import quasr_geom as qg  # existing: VMEC-boundary loader + flux-tangent bhat

HERE = Path(__file__).resolve().parent
DATADIR = HERE.parent / "data" / "quasr"
SERIAL_DIR = DATADIR / "serials"
BASE = "https://quasr.flatironinstitute.org/simsopt_serials"


# --------------------------------------------------------------------------- #
# Fetch + low-level graph helpers
# --------------------------------------------------------------------------- #
def fetch_serial(ID):
    z = f"{int(ID):07d}"
    SERIAL_DIR.mkdir(parents=True, exist_ok=True)
    fn = SERIAL_DIR / f"serial{z}.json"
    if not fn.exists():
        url = f"{BASE}/{z[:4]}/serial{z}.json"
        # explicit timeout: urlretrieve has none by default, so a hung socket
        # would block a worker forever. Fetch to a tmp path then atomic-rename
        # so a partial download can never masquerade as a cached serial.
        tmp = fn.with_suffix(".json.part")
        with urllib.request.urlopen(url, timeout=20) as resp, open(tmp, "wb") as out:
            out.write(resp.read())
        tmp.replace(fn)
    return fn


def _objs(serial_json):
    return json.loads(Path(serial_json).read_text())["simsopt_objs"]


def _is_ref(v):
    return isinstance(v, dict) and v.get("$type") == "ref"


def _arr(v):
    """Decode a simsopt/GSON numpy-array dict or a plain list -> np.ndarray."""
    if isinstance(v, dict) and v.get("@class") == "array":
        return np.asarray(v["data"], float)
    return np.asarray(v, float)


# --------------------------------------------------------------------------- #
# Curve / current reconstruction
# --------------------------------------------------------------------------- #
def _curve_xyz_fourier_points(obj, objs, n_samples=256):
    """Reconstruct a CurveXYZFourier as a closed polyline (n_samples,3).
    simsopt dof layout per coordinate (x,y,z): [c0, s1, c1, s2, c2, ..., sN, cN],
    X(t) = c0 + sum_{m=1}^N c_m cos(m*2pi t) + s_m sin(m*2pi t)."""
    order = int(obj["order"])
    dofs = objs[obj["dofs"]["value"]]
    x = _arr(dofs["x"])
    per = 2 * order + 1
    assert x.size == 3 * per, f"CurveXYZFourier dof size {x.size} != 3*{per}"
    t = np.linspace(0.0, 1.0, n_samples, endpoint=False)
    th = 2 * np.pi * t
    out = []
    for c in range(3):
        coeff = x[c * per:(c + 1) * per]
        val = np.full_like(t, coeff[0])
        for m in range(1, order + 1):
            s_m = coeff[2 * m - 1]; c_m = coeff[2 * m]
            val = val + c_m * np.cos(m * th) + s_m * np.sin(m * th)
        out.append(val)
    return np.stack(out, axis=1)


def _rotmat(phi, flip):
    """simsopt RotatedCurve rotation: Rz(phi), then (if flip) the stellarator-
    symmetry reflection diag(1,-1,-1)."""
    c, s = np.cos(phi), np.sin(phi)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    if flip:
        R = R @ np.diag([1.0, -1.0, -1.0])
    return R


def _resolve_curve(ref_id, objs, n_samples=256):
    obj = objs[ref_id]
    cls = obj["@class"]
    if cls == "CurveXYZFourier":
        return _curve_xyz_fourier_points(obj, objs, n_samples)
    if cls == "RotatedCurve":
        base = _resolve_curve(obj["curve"]["value"], objs, n_samples)
        R = _rotmat(float(obj["phi"]), bool(obj["flip"]))
        return base @ R.T
    raise ValueError(f"unsupported curve class {cls}")


def _resolve_current(ref_id, objs):
    obj = objs[ref_id]
    cls = obj["@class"]
    if cls == "Current":
        return float(obj["current"])
    if cls == "ScaledCurrent":
        return float(obj["scale"]) * _resolve_current(obj["current_to_scale"]["value"], objs)
    raise ValueError(f"unsupported current class {cls}")


def load_coils(ID, n_samples=256):
    """Return the full coil set as a list of (poly (n_samples,3), current)."""
    objs = _objs(fetch_serial(ID))
    coils = []
    for v in objs.values():
        if isinstance(v, dict) and v.get("@class") == "Coil":
            poly = _resolve_curve(v["curve"]["value"], objs, n_samples)
            cur = _resolve_current(v["current"]["value"], objs)
            coils.append((poly, cur))
    if not coils:
        raise ValueError(f"no Coil objects in serial for device {ID}")
    return coils


def total_filament_length(coils):
    L = 0.0
    for poly, _ in coils:
        seg = np.roll(poly, -1, axis=0) - poly
        L += float(np.sum(np.linalg.norm(seg, axis=1)))
    return L


# --------------------------------------------------------------------------- #
# Catalogue metadata
# --------------------------------------------------------------------------- #
def _catalogue_row(ID):
    with gzip.open(DATADIR / "catalogue.csv.gz", "rt") as f:
        for row in csv.DictReader(f):
            if int(float(row["ID"])) == int(ID):
                return row
    return None


def _symmetry_class(row):
    """Map QUASR's `helicity` column to the QS class. QUASR helicity encodes the
    quasisymmetry type of |B| in Boozer coordinates: 0 = QA (quasi-axisymmetric,
    N=0), nonzero = QH (quasi-helical, N=nfp). QUASR is a QS database, so QI is not
    represented (a QI device would need a separate source). This is the QS TYPE, a
    distinct axis from both the real-space gentleness (eps_eff) and the field-
    direction coherence C -- see docs/THEORY.md. A config may still override the
    class in the manifest."""
    if row and row.get("helicity") not in (None, ""):
        try:
            return "QA" if int(float(row["helicity"])) == 0 else "QH"
        except (TypeError, ValueError):
            return str(row["helicity"])
    return "unknown"


# --------------------------------------------------------------------------- #
# Public loader
# --------------------------------------------------------------------------- #
class QuasrDevice:
    def __init__(self, ID, coils, device, meta):
        self.ID = int(ID)
        self.coils = coils            # [(poly,current)]
        self.device = device          # quasr_geom.Device (boundary + flux-tangent bhat)
        self.meta = meta              # dict

    def source_sample(self, n_rho=16, n_theta=32, n_phi=48, rho_max=0.85,
                      model="one_minus_rho2"):
        """Points (xyz) and weights sampling the FUSION SOURCE VOLUME, mirroring the
        compiled source (births on (1-rho^2) flux surfaces up to rho_max*a). Weight
        per point = s(rho) * rho  (the rho*drho volume Jacobian for self-similar
        surfaces whose cross-section area ~ rho^2), theta/phi on a uniform grid.
        This is the s(x) used for the coherence metric C; a stage-4 consistency
        check compares it to the .so's actual emitted birth sites. Returns
        (xyz (N,3), weights (N,), rho (N,))."""
        import coherence_metrics as _cm
        rho_mid = (np.arange(n_rho) + 0.5) / n_rho * rho_max
        th = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
        ph = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        xyz, w, rr = [], [], []
        for rho in rho_mid:
            R, Z = self.device.RZ(TH, PH, float(rho))
            pts = np.stack([R * np.cos(PH), R * np.sin(PH), Z], axis=-1).reshape(-1, 3)
            xyz.append(pts)
            w.append(np.full(pts.shape[0], _cm.source_weights(rho, model) * rho))
            rr.append(np.full(pts.shape[0], rho))
        return np.vstack(xyz), np.concatenate(w), np.concatenate(rr)

    def boundary_grid(self, ntheta=64, nphi=128, rho=1.0):
        """LCFS grid: R,Z (ntheta,nphi) and Cartesian points + outward normals for
        the field audit (B.n ~ 0 on a flux surface)."""
        th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
        ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        R, Z = self.device.RZ(TH, PH, rho)
        X = R * np.cos(PH); Y = R * np.sin(PH)
        return dict(theta=th, phi=ph, R=R, Z=Z,
                    xyz=np.stack([X, Y, Z], axis=-1))


def load_device(ID, n_samples=256, coils=True):
    """Load coils + boundary + metadata for one QUASR device."""
    row = _catalogue_row(ID)
    nfp = int(float(row["nfp"])) if row else None
    iota = float(row["mean_iota"]) if row and row.get("mean_iota") else 0.0
    dev = qg.load_device(ID, iota)          # boundary + flux-tangent bhat
    if nfp is None:
        nfp = dev.nfp
    meta = dict(
        ID=int(ID), nfp=nfp, iota=iota,
        aspect=float(row["aspect_ratio"]) if row else None,
        minor_radius=float(row["minor_radius"]) if row else None,
        qs_error=float(row["qs_error"]) if row else None,
        nc_per_hp=int(float(row["nc_per_hp"])) if row and row.get("nc_per_hp") else None,
        total_coil_length=float(row["total_coil_length"]) if row and row.get("total_coil_length") else None,
        symmetry_class=_symmetry_class(row),
        R0=float(dev.R0),
    )
    cset = load_coils(ID, n_samples) if coils else None
    return QuasrDevice(ID, cset, dev, meta)


if __name__ == "__main__":
    import sys
    ID = int(sys.argv[1]) if len(sys.argv) > 1 else 59509
    d = load_device(ID)
    L = total_filament_length(d.coils)
    print(f"device {ID}: nfp {d.meta['nfp']} iota {d.meta['iota']} "
          f"aspect {d.meta['aspect']} class {d.meta['symmetry_class']}")
    print(f"coils: {len(d.coils)}  summed filament length {L:.4f} "
          f"(catalogue total_coil_length {d.meta['total_coil_length']})")
