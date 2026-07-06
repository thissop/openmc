#!/usr/bin/env python
"""HARD-GATE field audit. A config that fails is logged and SKIPPED, never
transported (no warn-and-continue). Two tiers:

CODE-INVARIANT gates (a real magnetostatic field / a correct parse MUST satisfy
these; a failure means a BUG, and they clear by orders of magnitude):
  1. UNIT     : |b_hat| = 1 everywhere sampled.
  2. DIV-FREE : div B ~ 0 on a finite-difference stencil.

PHYSICS coil-fit-QUALITY gate (measures a real property of the device, not a bug):
  3. BOUNDARY : B . n_hat / |B| on the LCFS. The QUASR coils are optimized so the
                nominal boundary is a flux surface, so the coil field should be
                nearly tangent to it. A WRONG parse (bad rotation/stellsym/current)
                gives an O(1) normal field; a correct parse gives the device's true
                residual. That residual is grid-converged and genuinely a few percent
                for aggressive (high-nfp/high-iota QH) configs, and discretization-
                limited (true value ~0) for gentle QA configs -- so we gate on the
                FINE-grid value (closest to truth) and report grid convergence. The
                threshold means 'the coil field is too far from this device's nominal
                equilibrium to trust b_hat for C', so the config is skipped. C itself
                is an integral over the source and is robust to a few-percent local
                field error, so passing configs at a few percent is defensible.
"""
from __future__ import annotations

import numpy as np

import coherence_metrics as cm

# code-invariant tolerances: a correct field clears these by orders of magnitude.
TOL_UNIT = 1e-9          # |b_hat|-1
TOL_DIV = 5e-3           # |div B| * L / <|B|>   (finite-difference-limited)
# coil-fit-quality tolerances on the LCFS B.n/|B|. A wrong parse fails by ~O(1);
# gentle QA clears ~1e-3; aggressive QH sits at a few e-2 (genuine, grid-converged).
TOL_BND_RMS = 0.08       # rms(B.n / |B|) -- skip devices whose coils fit worse than ~8%
TOL_BND_MAX = 0.30       # max(|B.n| / |B|)


def check_unit(field, pts):
    bh = field.bhat(pts)
    err = float(np.max(np.abs(np.linalg.norm(bh, axis=1) - 1.0)))
    return dict(name="unit", value=err, tol=TOL_UNIT, passed=err <= TOL_UNIT)


def check_divergence(field, pts, h, L):
    """Central-difference div B, normalized by <|B|>/L (so a smooth field gives a
    small dimensionless residual set by the stencil size h)."""
    B0 = field.B(pts)
    Bmag = float(np.mean(np.linalg.norm(B0, axis=1)))
    div = np.zeros(len(pts))
    for k, e in enumerate(np.eye(3)):
        Bp = field.B(pts + h * e)
        Bm = field.B(pts - h * e)
        div += (Bp[:, k] - Bm[:, k]) / (2 * h)
    resid = float(np.max(np.abs(div))) * L / (Bmag + 1e-300)
    return dict(name="div_free", value=resid, tol=TOL_DIV, passed=resid <= TOL_DIV)


def _surface_normals(xyz):
    """Outward unit normals of a closed (periodic theta,phi) surface grid
    xyz (ntheta,nphi,3), from dr/dtheta x dr/dphi. The (theta,phi) parameterization
    gives a globally CONSISTENT normal direction (up to one overall sign), so we set
    that single sign by the divergence theorem -- integral(r . n dA) = 3*Volume > 0
    iff outward. This is robust for NON-CONVEX (bean) cross-sections, where a
    per-point 'away from the centroid' test can flip individual normals."""
    dtheta = 0.5 * (np.roll(xyz, -1, axis=0) - np.roll(xyz, 1, axis=0))
    dphi = 0.5 * (np.roll(xyz, -1, axis=1) - np.roll(xyz, 1, axis=1))
    n = np.cross(dtheta, dphi)
    n /= (np.linalg.norm(n, axis=2, keepdims=True) + 1e-300)
    if np.sum(xyz * n) < 0.0:          # one global flip -> outward
        n = -n
    return n


def check_boundary_flux(field, boundary):
    """B . n_hat on the LCFS, normalized by |B|. `boundary` from
    QuasrDevice.boundary_grid (has 'xyz' (ntheta,nphi,3))."""
    xyz = boundary["xyz"]
    n = _surface_normals(xyz)
    pts = xyz.reshape(-1, 3)
    B = field.B(pts)
    Bmag = np.linalg.norm(B, axis=1)
    bdotn = np.sum(B * n.reshape(-1, 3), axis=1) / (Bmag + 1e-300)
    rms = float(np.sqrt(np.mean(bdotn ** 2)))
    mx = float(np.max(np.abs(bdotn)))
    passed = (rms <= TOL_BND_RMS) and (mx <= TOL_BND_MAX)
    return dict(name="boundary_flux", value=rms, max=mx,
                tol=TOL_BND_RMS, tol_max=TOL_BND_MAX, passed=passed)


def audit_coilfield(field, device, ntheta=64, nphi=128, verbose=False):
    """Full audit of a CoilField against its QuasrDevice boundary. Gates on the FINE
    grid; also evaluates a COARSE grid to report whether B.n is grid-converged
    (genuine coil-fit residual) or discretization-limited (true value even smaller).
    Returns (passed, report)."""
    bnd = device.boundary_grid(ntheta=ntheta, nphi=nphi, rho=1.0)
    bnd_coarse = device.boundary_grid(ntheta=max(24, ntheta // 2),
                                      nphi=max(48, nphi // 2), rho=1.0)
    # interior sample points on a few flux surfaces for unit + divergence
    thi = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    phi = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    TH, PH = np.meshgrid(thi, phi, indexing="ij")
    ipts = []
    for rho in (0.3, 0.6, 0.9):
        R, Z = device.device.RZ(TH, PH, rho)
        ipts.append(np.stack([R * np.cos(PH), R * np.sin(PH), Z], axis=-1).reshape(-1, 3))
    ipts = np.vstack(ipts)
    minor = float(device.meta.get("minor_radius") or 0.1)
    h = 1e-4 * max(device.meta.get("R0", 1.0), 1.0)

    bflux = check_boundary_flux(field, bnd)
    bflux_coarse = check_boundary_flux(field, bnd_coarse)
    # converged if fine ~ coarse; if fine << coarse the residual is discretization-
    # limited and the TRUE B.n is <= the fine value (so gating on fine is safe).
    bflux["coarse_value"] = bflux_coarse["value"]
    bflux["grid_converged"] = bool(
        bflux["value"] > 0 and abs(bflux_coarse["value"] - bflux["value"])
        / max(bflux["value"], 1e-30) < 0.25)
    checks = [check_unit(field, ipts)]
    # A solved equilibrium field (EquilibriumField, G3) is divergence-free BY
    # CONSTRUCTION and its portable map stores only the unit direction, so a
    # finite-difference div-B on it is meaningless -- skip the gate (mark it passed +
    # noted). The boundary B.n gate below still applies and is the meaningful check.
    if getattr(field, "div_free_by_construction", False):
        checks.append(dict(name="div_free", value=0.0, tol=TOL_DIV, passed=True,
                           note="skipped: divergence-free by construction "
                                "(equilibrium map stores unit b_hat only)"))
    else:
        checks.append(check_divergence(field, ipts, h=h, L=minor))
    checks.append(bflux)
    passed = all(c["passed"] for c in checks)
    report = dict(ID=device.ID, passed=passed, checks=checks)
    if verbose:
        print(f"[field_audit] device {device.ID}: PASS={passed}")
        for c in checks:
            extra = f" max={c.get('max'):.3e}" if "max" in c else ""
            print(f"   {c['name']:14s} value={c['value']:.3e} tol={c['tol']:.1e}"
                  f"{extra}  {'ok' if c['passed'] else 'FAIL'}")
    return passed, report


if __name__ == "__main__":
    import sys
    import quasr_loader as ql
    from biotsavart_field import CoilField
    ID = int(sys.argv[1]) if len(sys.argv) > 1 else 59509
    dev = ql.load_device(ID)
    fld = CoilField(dev.coils, meta=dev.meta)
    audit_coilfield(fld, dev, verbose=True)
