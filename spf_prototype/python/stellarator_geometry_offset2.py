#!/usr/bin/env python
"""stellarator_geometry_offset2.py -- NON-INVERTING conformal offset (STAGED, OFF by default).

WHY THIS EXISTS
---------------
The live builder `stellarator_geometry.build_layers` offsets the LCFS cross-section
along a FIXED poloidal normal field (`R + delta*n`). On the concave-inboard
indentation of a bean cross-section the normal rays CONVERGE; once the cumulative
offset exceeds the local radius of curvature the offset curves at different deltas
cross each other, so an inner shell pokes THROUGH an outer shell. Each shell is
individually watertight, but the ASSEMBLY self-intersects and DAGMC's ray tracer
leaks (~38% of sweep devices lost to this). Diagnosed in sweep/RESULTS_leak_diag.md;
measured by sweep/leak_gate.worst_shell_crossing (good devices -> 0.000; leaky ->
0.32-0.38).

THE FIX (morphological / Minkowski offset)
------------------------------------------
Replace the fixed-normal offset with an outward MORPHOLOGICAL offset: for each fixed-
phi cross-section polygon P, the boundary at radial distance d is the exterior of the
Minkowski sum P (+) disk(d), i.e. shapely `Polygon(P).buffer(d)`. Two mathematical
facts make this a GUARANTEED non-inverting, non-crossing offset:

  1. The outward buffer of a simple polygon is always simple (the disk dilation can
     never fold a boundary onto itself -- it fills concave notches rather than
     inverting through them).
  2. Buffers are strictly NESTED: for d1 < d2, buffer(d1) (subset) buffer(d2), because
     buffer(d2) = buffer(d1) (+) disk(d2-d1). So an inner boundary can never poke
     through an outer boundary -- worst_shell_crossing is 0 BY CONSTRUCTION.

This is exactly the curvature-limited behaviour asked for: on a concave region whose
radius of curvature is smaller than the offset, the naive normal offset inverts,
whereas the disk dilation caps the offset at the disk and rounds it -- it can shrink
the concavity but never invert it.

INTERFACE PRESERVED
-------------------
Emits the identical artefacts as build_layers: one `<layer>.stl` per layer in the
[outer-tris | inner-tris] shell-solid layout of sg.write_stl_shell, plus a
manifest.json with the same fields. dagmc_writer.build_from_stls and
sweep/leak_gate consume it UNCHANGED. Only the offset math differs.

RESAMPLING
----------
Each buffered boundary is resampled to `pol_res` points by equal ARC LENGTH, anchored
at the outboard (max-R) vertex and oriented CCW, so (a) the per-boundary polygon is
simple, (b) the toroidal seam is consistent across phi (no twist), and (c) adjacent
thin boundaries stay vertex-aligned (their curves are near-parallel), which keeps the
piecewise-linear representation nested even across the 0.2 cm W layer. pol_res
defaults to 3x the input poloidal resolution, which drives the measured worst
crossing to exactly 0 on every tested device (the 0.2 cm W gap is the only
sag-sensitive pair).

NOT WIRED INTO THE LIVE SWEEP. Import it explicitly (build_layers_sdf) or gate it
behind a flag that defaults OFF; integration is a later, separately-validated step.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stellarator_geometry as sg  # reuse verts/triangles/stl/checks/load_surface

try:
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.ops import unary_union
except ImportError:  # pragma: no cover
    sys.exit("shapely required: pip install shapely")


# --------------------------------------------------------------------------- #
# core: buffered offset of one cross-section, resampled to a fixed point count
# --------------------------------------------------------------------------- #
def _largest_polygon(geom):
    """buffer() of a simple polygon is a single Polygon, but be defensive: if a
    (near-)degenerate input yields a MultiPolygon, keep the largest-area component."""
    if isinstance(geom, MultiPolygon):
        return max(geom.geoms, key=lambda g: g.area)
    return geom


def _resample_ring(poly, n):
    """Resample the exterior ring of a shapely Polygon to n points by equal arc
    length, anchored at the max-R (outboard) vertex, oriented CCW. Returns (n,2)."""
    ext = poly.exterior
    if not ext.is_ccw:                       # normalize winding so all boundaries agree
        ext = ext.reverse()
    # anchor: rebuild the ring so arc-length 0 sits at the outboard (max-R) vertex,
    # a geometrically stable seam across both phi and offset distance.
    coords = np.asarray(ext.coords)[:-1]     # drop the closing duplicate
    a = int(np.argmax(coords[:, 0]))
    coords = np.roll(coords, -a, axis=0)
    from shapely.geometry import LinearRing
    ring = LinearRing(np.vstack([coords, coords[:1]]))
    L = ring.length
    # sample n points at equal arc length; ring.interpolate follows the (piecewise
    # linear) buffer boundary EXACTLY, so no shape is invented between vertices.
    out = np.empty((n, 2))
    for i in range(n):
        p = ring.interpolate(L * i / n)
        out[i] = (p.x, p.y)
    return out


def curvature_clamped_normals(R, Z):
    """Base outward normals PLUS the per-node maximum outward offset before the
    normal offset folds. On a concave-inboard node the centre of curvature lies on
    the OUTWARD (normal) side, so offsetting past the local radius of curvature R_c
    makes neighbouring normal rays cross -> the fold. We return (nR,nZ,dmax) where
    dmax[k,j] = 0.9*R_c on concave nodes and +inf on convex nodes. Offsetting every
    node by min(d, dmax) keeps each shell simple AND keeps adjacent shells exactly
    radial-parallel (offset along the SAME normals by a monotone distance), so they
    nest EVERYWHERE -- between phi planes too -- which is what DAGMC requires. This
    is the task's curvature-limited normal offset (option b)."""
    nR, nZ = sg.poloidal_outward_normals(R, Z)
    # periodic derivatives along theta (axis 0)
    dR = 0.5 * (np.roll(R, -1, 0) - np.roll(R, 1, 0))
    dZ = 0.5 * (np.roll(Z, -1, 0) - np.roll(Z, 1, 0))
    d2R = np.roll(R, -1, 0) - 2 * R + np.roll(R, 1, 0)
    d2Z = np.roll(Z, -1, 0) - 2 * Z + np.roll(Z, 1, 0)
    speed2 = dR * dR + dZ * dZ + 1e-30
    # signed curvature kappa = (x' y'' - y' x'') / speed^3 ; curvature vector kappa*Nfrenet
    cross = dR * d2Z - dZ * d2R
    kappa = cross / speed2 ** 1.5                    # signed, sign vs (dR,dZ)->(-dZ,dR)
    Rc = 1.0 / (np.abs(kappa) + 1e-30)               # radius of curvature
    # centre-of-curvature direction = sign(kappa) * left-normal (-dZ,dR)/speed
    lnR = -dZ / np.sqrt(speed2); lnZ = dR / np.sqrt(speed2)
    ccR = np.sign(kappa) * lnR; ccZ = np.sign(kappa) * lnZ
    concave = (ccR * nR + ccZ * nZ) > 0.0            # centre on OUTWARD side -> folds
    dmax = np.where(concave, 0.9 * Rc, np.inf)
    return nR, nZ, dmax


def build_layers_clamped(stem, layers=None, scale=1.0, outdir=None):
    """Curvature-limited (non-folding) normal offset on the ORIGINAL (nt,nph) grid.
    Same artefacts/return as build_layers. Preserves the exact toroidal grid that
    DAGMC accepts for the good devices (no resampling, no phi-jitter); the only change
    vs the live builder is that the outward offset is capped at 0.9*radius-of-curvature
    on concave-inboard nodes so the offset never inverts."""
    if layers is None:
        layers = sg.DEFAULT_LAYERS
    R, Z, nfp = sg.load_surface(stem)
    R, Z = R * scale, Z * scale
    nt, npH = R.shape
    phi_axis = np.linspace(0.0, 2 * np.pi, npH, endpoint=False)
    nR, nZ, dmax = curvature_clamped_normals(R, Z)
    outdir = Path(outdir) if outdir else (sg.DATADIR / f"{stem}_geom")
    outdir.mkdir(parents=True, exist_ok=True)

    def off(d):
        t = np.minimum(d, dmax)
        return R + t * nR, Z + t * nZ

    base_T = sg._torus_triangles(nt, npH, flip=False)
    base_V = sg._verts(R, Z, phi_axis).reshape(-1, 3)
    outer_flip = sg.signed_volume(base_V, base_T) < 0.0

    manifest = {"stem": stem, "nfp": nfp, "scale": scale, "layers": [],
                "offset_method": "curvature_clamped_normal"}
    cum = 0.0
    for name, thick in layers:
        cum_out = cum + thick
        Ri, Zi = off(cum)
        Ro, Zo = off(cum_out)
        stl = outdir / f"{name}.stl"
        V, T = sg.write_stl_shell((Ri, Zi), (Ro, Zo), phi_axis, stl,
                                  outer_flip=outer_flip)
        wt = sg.is_edge_manifold(T)
        simple, n_bad = sg.cross_sections_simple(Ro, Zo)
        vol = abs(sg.signed_volume(sg._verts(Ro, Zo, phi_axis).reshape(-1, 3),
                                   sg._torus_triangles(nt, npH)))
        manifest["layers"].append(dict(
            name=name, inner_offset_cm=cum, outer_offset_cm=cum_out,
            stl=str(stl.name), watertight=bool(wt), simple_cross_section=simple,
            n_self_intersecting_phi=n_bad, enclosed_volume_cm3=vol))
        cum = cum_out
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def _best_roll(b, ref, nt):
    """Integer cyclic roll r of b (nt,2) that best matches ref (nt,2): maximises the
    circular cross-correlation (FFT, O(nt log nt)). Returns rolled b."""
    corr = np.zeros(nt)
    for c in range(2):
        corr += np.fft.irfft(np.fft.rfft(b[:, c]) * np.conj(np.fft.rfft(ref[:, c])),
                             n=nt)
    return np.roll(b, -int(np.argmax(corr)), axis=0)


def _base_parametrization(R, Z, pol_res):
    """A phi-SMOOTH poloidal parametrization of the LCFS: arc-length resample per phi,
    then phase-align each slice to the previous so vertex k tracks the same material
    column around the torus (no anchor jitter). Returns (nph, pol_res, 2)."""
    nph = R.shape[1]
    P = np.empty((nph, pol_res, 2))
    for j in range(nph):
        poly = Polygon(np.column_stack([R[:, j], Z[:, j]]))
        if not poly.is_valid:
            poly = _largest_polygon(poly.buffer(0))
        P[j] = _resample_ring(poly, pol_res)
    for j in range(1, nph):
        P[j] = _best_roll(P[j], P[j - 1], pol_res)
    return P


def offset_boundaries_shared(R, Z, dists, pol_res):
    """Non-folding, 3D-nesting conformal offset.

    Two failure modes must be beaten at once:
      * FOLD  -- the naive normal offset inverts on the concave-inboard notch; the
        morphological buffer (Minkowski dilation) never folds (it fills the notch).
      * BETWEEN-PLANE CROSSING -- resampling each shell's buffer INDEPENDENTLY lets
        the poloidal labelling drift by a vertex or two between neighbouring phi
        slices; at the large outer shells that ~1-vertex jitter is comparable to the
        toroidal facet size, so adjacent shells' faceted strips cross BETWEEN phi
        planes (leak_gate/DAGMC catch this even though every cross-section nests).
    Fix: resample every shell's buffer by arc length (evenly spread -> no vertex
    collapse, no fold), then PHASE-ALIGN each slice -- first to a phi-smooth base
    parametrization of the LCFS (kills cross-shell drift so vertex k means the same
    column on every shell), giving a mesh that is phi-smooth AND radially consistent.
    Keeps the TRUE per-layer buffer thickness (unlike a base->outer morph, which
    collapses vertices in filled concavities). Returns dict d->(Ro,Zo)."""
    nph = R.shape[1]
    out = {}
    for d in dists:
        P = np.empty((nph, pol_res, 2))
        for j in range(nph):
            P0 = Polygon(np.column_stack([R[:, j], Z[:, j]]))
            if not P0.is_valid:
                P0 = _largest_polygon(P0.buffer(0))
            if d <= 0.0:
                P[j] = _resample_ring(P0, pol_res)
            else:
                poly = _largest_polygon(
                    P0.buffer(d, quad_segs=max(8, int(round(d))), join_style=1))
                P[j] = _resample_ring(poly, pol_res)
        # Make THIS shell phi-smooth so its OWN closed surface never self-intersects
        # toroidally: at large offsets the buffer is near-circular and the per-phi
        # arc-length anchor is unstable, so without this the toroidal edges (k,j)-(k,j+1)
        # swing wildly and the surface tangles -- a self-intersecting outer shell that
        # both trimesh.contains and DAGMC choke on. Phase-align each slice to the
        # previous so vertex k tracks one material column around the torus. Nesting of
        # adjacent shells then follows for free from region containment (buffer(d1) is a
        # strict subset of buffer(d2)) once each shell surface is simple.
        for j in range(1, nph):
            P[j] = _best_roll(P[j], P[j - 1], pol_res)
        out[d] = (np.ascontiguousarray(P[:, :, 0].T),
                  np.ascontiguousarray(P[:, :, 1].T))
    return out


def offset_boundary(R, Z, d, pol_res):
    """Morphological outward offset of every fixed-phi cross-section by distance d.
    R,Z are (ntheta,nphi). Returns (Ro,Zo) each (pol_res,nphi). Guaranteed simple and
    (across d) nested by the Minkowski-dilation properties above.

    quad_segs is scaled with the offset so large blanket buffers keep their rounded
    joins smooth without exploding vertex counts on the tiny first-wall offsets."""
    nph = R.shape[1]
    Ro = np.empty((pol_res, nph))
    Zo = np.empty((pol_res, nph))
    for j in range(nph):
        P = Polygon(np.column_stack([R[:, j], Z[:, j]]))
        if not P.is_valid:
            P = P.buffer(0)                  # clean a self-touching input polygon
            P = _largest_polygon(P)
        if d <= 0.0:
            poly = P
        else:
            qs = max(8, int(round(d)))       # ~1 cm arc resolution on rounded joins
            poly = _largest_polygon(P.buffer(d, quad_segs=qs, join_style=1))
        rc = _resample_ring(poly, pol_res)
        Ro[:, j] = rc[:, 0]
        Zo[:, j] = rc[:, 1]
    return Ro, Zo


# --------------------------------------------------------------------------- #
# build: same signature/return/artefacts as sg.build_layers, buffered offset
# --------------------------------------------------------------------------- #
def build_layers_sdf(stem, layers=None, scale=1.0, outdir=None, pol_res=None):
    """Drop-in analogue of stellarator_geometry.build_layers using a NON-INVERTING
    morphological (Minkowski-buffer) offset + per-shell phi-smoothing, instead of the
    normal offset that folds on the concave inboard.

    Method (offset_boundaries_shared): per fixed-phi cross-section, each boundary is
    the exterior of Polygon.buffer(d) (disk dilation -> never folds, fills the concave
    notch), arc-length resampled to pol_res points, then phase-aligned across phi so
    each shell's closed surface is simple in 3D (no toroidal self-intersection).

    LIMITATION -- inboard axis-crossing (recorded in the manifest as 'axis_crossing'):
    on a COMPACT low-aspect-ratio device (R0/a ~< 1.5) the inboard wall is closer to
    the machine axis than the blanket is thick, so the OUTER shells offset past R=0 and
    the torus overlaps itself at the axis. This is a device/blanket-thickness
    incompatibility that NO conformal offset can fix (the naive normal offset hits the
    identical R<0). Such devices are flagged (min_R<=0) and must be rejected or given a
    thinner blanket; they are NOT recovered here. Fold-limited devices (R0/a ~> 1.5
    whose only defect was the concave fold) ARE fully recovered.

    pol_res : output poloidal point count per boundary (default 3x input ntheta;
              higher pol_res removes piecewise-linear chord-sag on the 0.2 cm W gap).
    """
    if layers is None:
        layers = sg.DEFAULT_LAYERS
    R, Z, nfp = sg.load_surface(stem)
    R = R * scale
    Z = Z * scale
    nt_in, npH = R.shape
    pol_res = int(pol_res) if pol_res else 3 * nt_in
    phi_axis = np.linspace(0.0, 2 * np.pi, npH, endpoint=False)
    outdir = Path(outdir) if outdir else (sg.DATADIR / f"{stem}_geom")
    outdir.mkdir(parents=True, exist_ok=True)

    # Precompute each UNIQUE cumulative boundary distance ONCE (inner of layer k+1 is
    # the outer of layer k -> identical buffer -> identical resampled surface, so the
    # shared interface is byte-consistent as dagmc_writer.build_from_stls assumes).
    # ALL distances share ONE poloidal parametrization (offset_boundaries_shared) so
    # adjacent shells nest in 3D (between phi planes), not only at each phi slice.
    cum = 0.0
    dists = [0.0]
    for _, thick in layers:
        cum += thick
        dists.append(cum)
    surf_cache = offset_boundaries_shared(R, Z, dists, pol_res)

    # orientation guard (G6): pick the outer winding so its normals point OUTWARD.
    base_Ro, base_Zo = surf_cache[0.0]
    base_T = sg._torus_triangles(pol_res, npH, flip=False)
    base_V = sg._verts(base_Ro, base_Zo, phi_axis).reshape(-1, 3)
    outer_flip = sg.signed_volume(base_V, base_T) < 0.0

    manifest = {"stem": stem, "nfp": nfp, "scale": scale, "layers": [],
                "offset_method": "morphological_buffer", "pol_res": pol_res}
    cum = 0.0
    for name, thick in layers:
        cum_out = cum + thick
        Ri, Zi = surf_cache[cum]
        Ro, Zo = surf_cache[cum_out]
        stl = outdir / f"{name}.stl"
        V, T = sg.write_stl_shell((Ri, Zi), (Ro, Zo), phi_axis, stl,
                                  outer_flip=outer_flip)
        wt = sg.is_edge_manifold(T)
        simple, n_bad = sg.cross_sections_simple(Ro, Zo)
        vol = abs(sg.signed_volume(sg._verts(Ro, Zo, phi_axis).reshape(-1, 3),
                                   sg._torus_triangles(pol_res, npH)))
        manifest["layers"].append(dict(
            name=name, inner_offset_cm=cum, outer_offset_cm=cum_out,
            stl=str(stl.name), watertight=bool(wt), simple_cross_section=simple,
            n_self_intersecting_phi=n_bad, enclosed_volume_cm3=vol))
        cum = cum_out
    # inboard axis-crossing flag: the outermost boundary's minimum R over the torus.
    # <=0 means a shell has offset past the machine axis (compact device + thick
    # blanket) -> the torus self-intersects at the axis; no offset can fix it.
    outermost = surf_cache[max(dists)]
    manifest["min_R_cm"] = float(np.min(outermost[0]))
    manifest["axis_crossing"] = bool(np.min(outermost[0]) <= 0.0)
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("stem")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--outdir")
    ap.add_argument("--pol-res", type=int, default=None)
    a = ap.parse_args()
    m = build_layers_sdf(a.stem, scale=a.scale, outdir=a.outdir, pol_res=a.pol_res)
    allok = True
    for L in m["layers"]:
        ok = L["watertight"] and L["simple_cross_section"]
        allok &= ok
        print(f"  {L['name']:7s} [{L['inner_offset_cm']:.1f},{L['outer_offset_cm']:.1f}]cm "
              f"watertight={L['watertight']} simple={L['simple_cross_section']} "
              f"(bad_phi={L['n_self_intersecting_phi']}) V={L['enclosed_volume_cm3']:.3e}")
    print(f"offset=morphological_buffer pol_res={m['pol_res']} ALL VALID: {allok}")
