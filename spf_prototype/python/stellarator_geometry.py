#!/usr/bin/env python
"""Conformal stellarator wall builder (PORTABLE, architecture-independent).

Takes a convention-free LCFS surface grid R(theta,phi), Z(theta,phi) (written by
desc_to_fieldmap.py) and builds nested CONFORMAL radial-build surfaces by offsetting
the boundary outward along the poloidal cross-section normal. Each layer is emitted
as a watertight closed triangle-mesh SHELL SOLID in an STL file; build_dagmc.py
(Ginsburg, conda) converts the STLs to a DAGMC .h5m, and run_conformal.py runs
OpenMC+DAGMC with our SPF compiled source.

This module is pure numpy (no DAGMC/OpenMC), so the geometry math is built and
UNIT-TESTED in the aarch64 sandbox; only the STL->.h5m and transport steps need
the compiled stack on x86. Conformal (not axisymmetric) walls are the point: they
let us test whether SPF directional-steering benefits survive non-axisymmetric
"smearing" of the wall.

Offset convention: at fixed phi the cross-section is a closed curve (R(theta),
Z(theta)); we offset along its outward unit normal in the (R,Z) plane. The small
toroidal (d/dphi) component of the true 3-D surface normal is neglected -- standard
for a first conformal build and a documented approximation (DEFERRED.md).

INJECT(helios): swap the surface grid + scale + layer thicknesses for real Helios
values; everything downstream is identical.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DATADIR = REPO / "spf_prototype" / "data"

# Default reactor-scale conformal radial build (cm), plasma -> out. INJECT(helios).
# Reuses the material identities of reactor_model (W/steel/Be/FLiBe) + shield+coil.
DEFAULT_LAYERS = [
    ("sol", 5.0),       # scrape-off gap (vacuum)
    ("W", 0.2),         # tungsten first-wall armor
    ("steel", 3.8),     # RAFM structural first wall
    ("Be", 2.0),        # multiplier
    ("FLiBe", 50.0),    # breeder
    ("shield", 40.0),   # WC/steel shield
    ("coil", 15.0),     # coil block (damage tallies)
]


def load_surface(stem):
    d = np.load(DATADIR / f"{stem}_surface.npz")
    return np.asarray(d["R"]), np.asarray(d["Z"]), int(d["nfp"])


def poloidal_outward_normals(R, Z):
    """Outward unit normals in the (R,Z) cross-section plane at each (theta,phi).
    Tangent t = d/dtheta (R,Z); normal n = (t_Z, -t_R) or its negation, chosen to
    point away from the per-phi cross-section centroid. R,Z shape (ntheta,nphi)."""
    dR = np.gradient(R, axis=0, edge_order=2)   # periodic-ish; theta is periodic
    dZ = np.gradient(Z, axis=0, edge_order=2)
    # use a true periodic difference for the theta derivative
    dR = 0.5 * (np.roll(R, -1, axis=0) - np.roll(R, 1, axis=0))
    dZ = 0.5 * (np.roll(Z, -1, axis=0) - np.roll(Z, 1, axis=0))
    nR, nZ = dZ.copy(), -dR.copy()
    mag = np.sqrt(nR ** 2 + nZ ** 2) + 1e-30
    nR, nZ = nR / mag, nZ / mag
    # orient outward: away from the centroid of each phi cross-section
    Rc = R.mean(axis=0, keepdims=True); Zc = Z.mean(axis=0, keepdims=True)
    out = (R - Rc) * nR + (Z - Zc) * nZ
    sign = np.where(out >= 0, 1.0, -1.0)
    return nR * sign, nZ * sign


def offset_surface(R, Z, nR, nZ, delta):
    return R + delta * nR, Z + delta * nZ


def _verts(R, Z, phi_axis):
    """(ntheta,nphi) R,Z -> Cartesian vertex array, shape (ntheta,nphi,3)."""
    ph = phi_axis[None, :]
    X = R * np.cos(ph); Y = R * np.sin(ph)
    return np.stack([X, Y, np.broadcast_to(Z, X.shape)], axis=-1)


def _torus_triangles(ntheta, nphi, flip=False):
    """Triangle connectivity for a closed (periodic theta & phi) structured grid.
    Returns (Ntri,3) int index array into a flattened (ntheta*nphi,3) vertex list."""
    tris = []
    for i in range(ntheta):
        i1 = (i + 1) % ntheta
        for j in range(nphi):
            j1 = (j + 1) % nphi
            a = i * nphi + j; b = i1 * nphi + j
            c = i1 * nphi + j1; d = i * nphi + j1
            if not flip:
                tris.append((a, b, c)); tris.append((a, c, d))
            else:
                tris.append((a, c, b)); tris.append((a, d, c))
    return np.asarray(tris, dtype=np.int64)


def write_stl_shell(inner_surf, outer_surf, phi_axis, path):
    """Write a watertight shell-solid STL = region between inner and outer toroidal
    surfaces (outer normals out, inner normals in). Each surf is (R,Z) tuple."""
    Ri, Zi = inner_surf; Ro, Zo = outer_surf
    nt, npH = Ri.shape
    Vi = _verts(Ri, Zi, phi_axis).reshape(-1, 3)
    Vo = _verts(Ro, Zo, phi_axis).reshape(-1, 3)
    Ti = _torus_triangles(nt, npH, flip=True)    # inner: inward-facing
    To = _torus_triangles(nt, npH, flip=False)   # outer: outward-facing
    V = np.vstack([Vo, Vi]); T = np.vstack([To, Ti + len(Vo)])
    _write_binary_stl(V, T, path)
    return V, T


def _write_binary_stl(V, T, path):
    tri = V[T]  # (Ntri,3,3)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-30)
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(T)))
        for k in range(len(T)):
            f.write(struct.pack("<3f", *n[k]))
            for v in tri[k]:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


# ---- correctness checks (run in the sandbox; no compiled stack needed) ----
def is_edge_manifold(T):
    """Watertight closed mesh <=> every edge shared by exactly 2 triangles."""
    from collections import Counter
    e = Counter()
    for a, b, c in T:
        for u, v in ((a, b), (b, c), (c, a)):
            e[(min(u, v), max(u, v))] += 1
    return all(cnt == 2 for cnt in e.values())


def signed_volume(V, T):
    tri = V[T]
    return float(np.sum(np.einsum("ij,ij->i",
                 tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))) / 6.0)


def _poly_is_simple(rr, zz):
    """True if the closed polygon (rr,zz) has no non-adjacent edge crossings.
    Works for NON-CONVEX polygons (stellarator cross-sections are bean-shaped),
    unlike a turn-sign test. This is the real detector of inboard fold-over from
    over-offsetting. O(n^2) per cross-section."""
    n = len(rr)
    P = np.column_stack([rr, zz])
    def x2(o, a, b):  # 2D cross (b-o) x (a-o)
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    def seg_cross(a, b, c, d):
        return ((x2(c, d, a) > 0) != (x2(c, d, b) > 0)) and \
               ((x2(a, b, c) > 0) != (x2(a, b, d) > 0))
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        for j in range(i + 1, n):
            if j == i or (i + 1) % n == j or (j + 1) % n == i:
                continue
            if seg_cross(a, b, P[j], P[(j + 1) % n]):
                return False
    return True


def cross_sections_simple(R, Z):
    """No-self-intersection guard: every per-phi cross-section polygon must stay
    SIMPLE. Naive normal-offset folds on the concave inboard side once the offset
    exceeds the local radius of curvature; that shows up as a self-crossing
    cross-section. Returns (all_simple, n_bad_phi)."""
    npH = R.shape[1]
    bad = sum(0 if _poly_is_simple(R[:, j], Z[:, j]) else 1 for j in range(npH))
    return bad == 0, bad


def build_layers(stem, layers=DEFAULT_LAYERS, scale=1.0, outdir=None):
    """Build per-layer conformal shell STLs from the LCFS surface grid.
    scale multiplies the device size (field B-hat is scale-invariant; rescale the
    field-map grid bounds by the same factor on Ginsburg). Returns a manifest."""
    R, Z, nfp = load_surface(stem)
    R = R * scale; Z = Z * scale
    nt, npH = R.shape
    phi_axis = np.linspace(0.0, 2 * np.pi, npH, endpoint=False)
    nR, nZ = poloidal_outward_normals(R, Z)
    outdir = Path(outdir) if outdir else (DATADIR / f"{stem}_geom")
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {"stem": stem, "nfp": nfp, "scale": scale, "layers": []}
    cum = 0.0
    Rin, Zin = R.copy(), Z.copy()
    prev_vol = signed_volume(_verts(Rin, Zin, phi_axis).reshape(-1, 3),
                             _torus_triangles(nt, npH))
    for name, thick in layers:
        cum_out = cum + thick
        Ro, Zo = offset_surface(R, Z, nR, nZ, cum_out)
        Ri, Zi = offset_surface(R, Z, nR, nZ, cum)
        stl = outdir / f"{name}.stl"
        V, T = write_stl_shell((Ri, Zi), (Ro, Zo), phi_axis, stl)
        wt = is_edge_manifold(T)
        simple, n_bad = cross_sections_simple(Ro, Zo)
        vol = abs(signed_volume(_verts(Ro, Zo, phi_axis).reshape(-1, 3),
                                _torus_triangles(nt, npH)))
        manifest["layers"].append(dict(
            name=name, inner_offset_cm=cum, outer_offset_cm=cum_out,
            stl=str(stl.name), watertight=bool(wt), simple_cross_section=simple,
            n_self_intersecting_phi=n_bad, enclosed_volume_cm3=vol))
        cum = cum_out
    import json
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    stem = sys.argv[1] if len(sys.argv) > 1 else "equil_precise_qa"
    # default scale -> Helios-like minor radius (precise_QA a~17cm *10 = ~170cm ~ Helios
    # a=180cm); build 116cm < a so the conformal offset is self-intersection-free.
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    m = build_layers(stem, scale=scale)
    print(f"built {len(m['layers'])} conformal layers (scale={scale}) -> "
          f"data/{stem}_geom/")
    allok = True
    for L in m["layers"]:
        ok = L["watertight"] and L["simple_cross_section"]
        allok &= ok
        print(f"  {L['name']:7s} [{L['inner_offset_cm']:.1f},{L['outer_offset_cm']:.1f}]cm "
              f"watertight={L['watertight']} simple={L['simple_cross_section']} "
              f"(bad_phi={L['n_self_intersecting_phi']}) V={L['enclosed_volume_cm3']:.3e}cm3")
    print(f"ALL LAYERS VALID: {allok}" + ("" if allok else
          "  <- thick build self-intersects on the concave inboard side; increase "
          "scale or thin the build (naive normal-offset limit; see DEFERRED.md)"))
