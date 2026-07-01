#!/usr/bin/env python
"""Targeted demonstration that the inboard MC<->analytic residual is the limit of the
FILAMENTARY POINT-PATCH-FLUX analytic on a strongly-curved wall under strong shaping --
shown by direct construction, not only by elimination.

Minimal, controlled setting: ONE helical source filament with a tunable shaping
amplitude `s`, emitting the spin-polarized kernel about a toroidal field, onto three
bare walls at the same standoff:
  * inner cylinder  (curves toward the source -- the 'inboard' case),
  * outer cylinder  (curves away from the source -- the 'outboard' case),
  * flat plane       (zero curvature -- the control).
Two evaluations of the same physics:
  * ANALYTIC  = the point-patch free-streaming flux integral (the method under test;
    kernel evaluated at the patch-centre line of sight), identical in form to anarrima.
  * RAYTRACE  = an exact analytic ray/surface intersection of the SAME sampled births
    and directions (the bit-parity SPF sampler), i.e. what OpenMC transport does, with
    no point-patch approximation. This is the geometry-exact reference.

Sweeping `s` shows the analytic<->raytrace directionality discrepancy stay ~0 on the
flat wall at all shaping, stay small on the outer cylinder, and grow on the inner
cylinder -- reproducing the toy(gentle)->QUASR(shaped) inboard trend in isolation and
pinning the cause to curvature x shaping in the point-patch model.

Run: $HOME/spf_venv/bin/python spf_prototype/python/bare_cylinder_demo.py
"""
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
import spf_mirror as mir  # noqa: E402

NFP = 3
R0 = 1.0
AR, AZ = 0.20, 0.20          # shaping harmonic amplitudes (x s)
R_IN, R_OUT = 0.74, 1.30     # inner / outer cylinder radii
X_PLANE = 0.74               # flat wall (tangent to inner cylinder at phi=0)
ZW = 0.34                    # half-height of the wall band we score
NZ = 14                      # z patches/bins
NB = 600_000                 # raytrace births per mode
_GL = leggauss(256)


def loop(phi, s):
    """R, Z of the single shaped filament; s scales the non-axisymmetric excursion."""
    R = R0 + s * AR * np.cos(NFP * phi)
    Z = s * AZ * np.sin(NFP * phi)
    return R, Z


def loop_d(phi, s):
    R, Z = loop(phi, s)
    Rd = -s * AR * NFP * np.sin(NFP * phi)
    Zd = s * AZ * NFP * np.cos(NFP * phi)
    return R, Z, Rd, Zd


def bhat(phi):
    """Toroidal field direction phi-hat (unit)."""
    return np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=-1)


# --------------------------------------------------------------------------- #
# ANALYTIC: point-patch free-streaming flux integral (kernel at patch centre)
# --------------------------------------------------------------------------- #
def _vis_arc(s, r, psi, zw, wall):
    """Visible phi-arc (N>0) for the patch; + inner-cylinder occlusion for inner wall."""
    phi = np.linspace(-np.pi, np.pi, 1440, endpoint=False)
    R, Z = loop(phi, s)
    Dx = R * np.cos(phi) - r; Dy = R * np.sin(phi); Dz = Z - zw
    N = np.cos(psi) * Dx + np.sin(psi) * Dz
    vis = N > 1e-12
    if wall == "inner":
        L2 = Dx * Dx + Dy * Dy
        t = np.clip(-(r * Dx) / np.maximum(L2, 1e-30), 0.0, 1.0)
        vis &= np.hypot(r + t * Dx, t * Dy) >= r - 1e-9
    if not vis.any():
        return []
    d = phi[1] - phi[0]; arcs = []; i = 0
    while i < len(phi):
        if vis[i]:
            j = i
            while j < len(phi) and vis[j]:
                j += 1
            arcs.append((phi[i] - 0.5 * d, phi[j - 1] + 0.5 * d)); i = j
        else:
            i += 1
    return arcs


def _quad(s, r, psi, zw, a, b, mode):
    x, w = _GL
    half = 0.5 * (b - a); phi = 0.5 * (b + a) + half * x
    R, Z, Rd, Zd = loop_d(phi, s)
    J = np.sqrt(R ** 2 + Rd ** 2 + Zd ** 2)
    Dx = R * np.cos(phi) - r; Dy = R * np.sin(phi); Dz = Z - zw
    D = Dx ** 2 + Dy ** 2 + Dz ** 2
    N = np.cos(psi) * Dx + np.sin(psi) * Dz
    base = J * N / D ** 1.5
    if mode == "iso":
        integ = base
    else:
        B = bhat(phi)
        cos2 = (B[:, 0] * Dx + B[:, 1] * Dy + B[:, 2] * Dz) ** 2 / D
        integ = base * (1 - cos2) * 1.5 if mode == "A" else base * (0.25 + 0.75 * cos2) * 2.0
    return half * np.sum(w * integ)


def analytic_profile(wall, s):
    """A/iso per z-patch on a wall (inner/outer cylinder or plane)."""
    if wall == "inner":
        r, psi = R_IN, 0.0           # inward normal +R
    elif wall == "outer":
        r, psi = R_OUT, np.pi        # inward normal -R
    else:
        r, psi = X_PLANE, 0.0        # plane x=const, inward normal +x == +R at phi=0
    zs = np.linspace(-ZW, ZW, NZ)
    out = {}
    for m in ("iso", "A"):
        v = np.zeros(NZ)
        for i, zp in enumerate(zs):
            for a, b in _vis_arc(s, r, psi, zp, wall):
                if b - a > 1e-3:
                    v[i] += _quad(s, r, psi, zp, a, b, m)
        out[m] = v
    return zs, out["A"] / np.maximum(out["iso"], 1e-30)


# --------------------------------------------------------------------------- #
# RAYTRACE: exact ray/surface intersection of the sampled births (the MC truth)
# --------------------------------------------------------------------------- #
def _hit_z(x0, y0, z0, u, wall):
    """z of first intersection with the wall, or nan if it misses."""
    ux, uy, uz = u[:, 0], u[:, 1], u[:, 2]
    if wall == "plane":
        t = np.where(np.abs(ux) > 1e-12, (X_PLANE - x0) / ux, -1.0)
        return np.where(t > 0, z0 + t * uz, np.nan)
    r = R_IN if wall == "inner" else R_OUT
    a = ux ** 2 + uy ** 2
    b = 2 * (x0 * ux + y0 * uy)
    c = x0 ** 2 + y0 ** 2 - r ** 2
    disc = b ** 2 - 4 * a * c
    sq = np.sqrt(np.maximum(disc, 0.0))
    t1 = (-b - sq) / (2 * a); t2 = (-b + sq) / (2 * a)
    t = np.where(t1 > 1e-9, t1, t2)          # first forward crossing
    ok = (disc > 0) & (t > 1e-9)
    return np.where(ok, z0 + t * uz, np.nan)


def _sample_dirs(s, mode, seed):
    """Births (shared across walls) and SPF-sampled directions for one mode."""
    rng_pos = np.random.default_rng(1)                     # SAME births for all modes
    phi = rng_pos.uniform(0, 2 * np.pi, NB)
    R, Z, Rd, Zd = loop_d(phi, s)
    J = np.sqrt(R ** 2 + Rd ** 2 + Zd ** 2)
    x0, y0, z0 = R * np.cos(phi), R * np.sin(phi), Z
    B = bhat(phi)
    mw = mir.make_mode_weights(*({"iso": (1 / 3, 1 / 3, 1 / 3), "A": (1.0, 0.0, 0.0)}[mode]))
    rng = mir.Rng(seed)
    u = np.array([mir.sample_global_direction(mw, tuple(B[k]), rng) for k in range(NB)])
    return x0, y0, z0, J, u


def raytrace_band_load(s):
    """Total in-band (|z|<0.7 ZW) A and iso loads on inner and outer cylinders, from a
    single shared set of ray-traced births (directions sampled once per mode)."""
    zc = 0.7 * ZW
    loads = {"inner": {}, "outer": {}}
    counts = {"inner": 0, "outer": 0}
    for mode, seed in (("iso", 100), ("A", 200)):
        x0, y0, z0, J, u = _sample_dirs(s, mode, seed)
        for wall in ("inner", "outer"):
            zh = _hit_z(x0, y0, z0, u, wall)
            good = np.isfinite(zh) & (np.abs(zh) <= zc)
            loads[wall][mode] = J[good].sum()
            if mode == "iso":
                counts[wall] = int(good.sum())
    aiso = {w: loads[w]["A"] / loads[w]["iso"] for w in loads}
    return aiso, counts


def analytic_band_load(wall, s):
    """Total in-band A/iso from the point-patch flux integral (fine z patches)."""
    if wall == "inner":
        r, psi = R_IN, 0.0
    else:
        r, psi = R_OUT, np.pi
    tot = {"iso": 0.0, "A": 0.0}
    for zp in np.linspace(-0.7 * ZW, 0.7 * ZW, 40):
        for m in ("iso", "A"):
            for a, b in _vis_arc(s, r, psi, zp, wall):
                if b - a > 1e-3:
                    tot[m] += _quad(s, r, psi, zp, a, b, m)
    return tot["A"] / tot["iso"]


if __name__ == "__main__":
    print("In-band A/iso: analytic point-patch flux vs exact ray-trace, per wall, "
          "vs shaping s.")
    print(f"{'s':>5} | {'INNER (convex->src)':>28} | {'OUTER (concave->src)':>28}")
    print(f"{'':>5} | {'ana':>8} {'ray':>8} {'|Δ|':>6} {'hits':>6} | "
          f"{'ana':>8} {'ray':>8} {'|Δ|':>6} {'hits':>6}")
    rows = []
    for s in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        ray, cnt = raytrace_band_load(s)
        ai = analytic_band_load("inner", s); ao = analytic_band_load("outer", s)
        di = abs(ai - ray["inner"]); do = abs(ao - ray["outer"])
        rows.append((s, ai, ray["inner"], di, ao, ray["outer"], do))
        print(f"{s:5.1f} | {ai:8.3f} {ray['inner']:8.3f} {di:6.3f} {cnt['inner']:6d} | "
              f"{ao:8.3f} {ray['outer']:8.3f} {do:6.3f} {cnt['outer']:6d}")
    np.savez("/tmp/bare_cylinder_demo.npz", rows=np.array(rows))
    print("saved /tmp/bare_cylinder_demo.npz")
