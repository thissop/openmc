"""Turn an adjoint magnet-importance map into a shield-PLACEMENT priority field.

This is the through-shield replacement for the free-streaming culprit map, wired to
the shield optimizer's control geometry. The adjoint solve gives psi_dagger(r) over
plasma voxels = the importance of a source neutron born at r to a target coil's
heating THROUGH the shield+blanket (adjoint_importance.py). The optimizer's control
is a smooth thickness field delta(theta, phi) on the (poloidal, toroidal) shield
surface (thickness_field.ThicknessField). This module projects the plasma-volume
importance onto that (theta, phi) grid, so the optimizer knows WHERE to trade
breeder for shield to protect the coil.

Mapping: each plasma voxel at (x,y,z) has toroidal angle phi = atan2(y, x) and
poloidal angle theta = atan2(z, R - R0), R = hypot(x, y), about the magnetic axis
major radius R0. Its importance is accumulated (Gaussian-splatted for smoothness)
onto the (theta, phi) grid -> priority(theta, phi). The peak of that field is the
data-driven target (theta0, phi0) for ThicknessField.gaussian_bump, and the whole
field can weight the optimizer's objective.

Drop-in vs the free-streaming culprit map: culprit_59509.py attributes the FIRST-WALL
load to plasma voxels by free streaming (ignores the shield); this attributes the
COIL load through the shield via the adjoint transport solve. Same idea (which plasma
region drives the hotspot), correct physics for magnets.
"""
from __future__ import annotations

import numpy as np

_WARNED = {"uniform": False}   # one-time guard so uniform emissivity is never used silently


def _voxel_centers(ll, ur, dim):
    return [np.linspace(ll[i], ur[i], dim[i] + 1)[:-1] + 0.5 * (ur[i] - ll[i]) / dim[i]
            for i in range(3)]


def plasma_source_on_mesh(fluxmap, ll, ur, dim, scale=100.0, emissivity="uniform"):
    """Rasterize the plasma neutron SOURCE S(r) onto the importance mesh.

    Source density = |sqrt(g)| (plasma volume element) times a reactivity weight r(rho):
      emissivity='uniform'    -> r=1 (the S_v the free-streaming culprit map uses)
      emissivity='bosch_hale' -> core-peaked DT reaction rate n^2<sigma v>(T) on the
                                 Miralles-Dolz profiles (shield_opt/emissivity.py).
    NOTE the ADJOINT importance psi_dagger(r) is reactivity-INDEPENDENT (pure geometry +
    coil response); the reactivity profile enters ONLY here, in S(r) -- so switching
    profiles needs no new transport run, just a re-weight. The VMEC fluxmap gives
    (R, Z, phi, rho, sqrt(g)); we scale m->cm and histogram onto the mesh. Returns S.
    """
    fm = np.load(fluxmap) if isinstance(fluxmap, str) else fluxmap
    R = fm["R"] * scale
    Z = fm["Z"] * scale
    phi = fm["phi"]
    w = np.abs(fm["sqrtg"])
    if emissivity == "uniform" and not _WARNED["uniform"]:
        import sys
        print("NOTE: plasma source S(r) uses UNIFORM emissivity (non-physical -- "
              "over-weights the shaped plasma edge; ~3x optimistic vs core-peaked "
              "Bosch-Hale on QH). Legitimate only for validation/geometry-isolation "
              "or matching a uniform-source forward run.", file=sys.stderr, flush=True)
        _WARNED["uniform"] = True
    if emissivity != "uniform":
        import os
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from emissivity import emissivity_weight
        rho = np.asarray(fm["rho"], dtype=float)            # (nrho,) per flux surface
        r_of_rho = emissivity_weight(rho, kind=emissivity)  # normalized reactivity r(rho)
        w = w * r_of_rho[:, None, None]                     # broadcast over (rho, theta, zeta)
    xs = (R * np.cos(phi)).ravel()
    ys = (R * np.sin(phi)).ravel()
    zs = Z.ravel()
    edges = [np.linspace(ll[i], ur[i], dim[i] + 1) for i in range(3)]
    S, _ = np.histogramdd(np.column_stack([xs, ys, zs]), bins=edges, weights=w.ravel())
    return S


def contributon(map_npz, fluxmap, scale=100.0, emissivity="uniform"):
    """Coil-attribution CONTRIBUTON C(r) = S(r) * psi_dagger(r).

    The bare adjoint flux psi_dagger peaks in the optically-thin central vacuum
    (physically irrelevant); weighting by the plasma source S(r) restricts it to
    where coil-bound neutrons are actually BORN -- the true "which plasma region
    drives this coil" map, the through-shield analog of the free-streaming culprit
    map. Returns a dict mirroring the map npz but with 'importance' -> the contributon,
    so the (theta, phi) placement pipeline runs on it unchanged.
    """
    d = np.load(map_npz) if isinstance(map_npz, str) else map_npz
    imp = d["importance"].astype(float)
    ll, ur, dim = d["lower_left"], d["upper_right"], d["dimension"]
    S = plasma_source_on_mesh(fluxmap, ll, ur, dim, scale=scale, emissivity=emissivity)
    C = S * imp
    out = {k: d[k] for k in d.files} if hasattr(d, "files") else dict(d)
    out["importance"] = C
    out["bare_importance"] = imp
    out["plasma_source"] = S
    return out


def importance_angles(npz, R0=None):
    """Return (phi, theta, weight) for every nonzero importance voxel (radians).

    phi   = toroidal angle atan2(y, x) in [-pi, pi)
    theta = poloidal angle atan2(z, R - R0) about the magnetic axis
    R0    = magnetic-axis major radius (cm); if None, use the importance-weighted
            mean cylindrical radius of the nonzero voxels (a robust auto-estimate).
    """
    d = np.load(npz) if isinstance(npz, str) else npz
    imp = d["importance"].astype(float)
    ll, ur, dim = d["lower_left"], d["upper_right"], d["dimension"]
    xc, yc, zc = _voxel_centers(ll, ur, dim)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    m = imp > 0
    x, y, z, w = X[m], Y[m], Z[m], imp[m]
    R = np.hypot(x, y)
    if R0 is None:
        R0 = float(np.sum(R * w) / np.sum(w))     # importance-weighted mean major radius
    phi = np.arctan2(y, x)
    theta = np.arctan2(z, R - R0)
    return phi, theta, w, R0


def placement_priority(npz, toroidal_angles_deg, poloidal_angles_deg,
                       R0=None, width_deg=15.0):
    """Splat the voxel importance onto the optimizer's (theta, phi) control grid.

    toroidal_angles_deg, poloidal_angles_deg: the SAME grids the ThicknessField uses.
    width_deg: Gaussian splat width (smoothness of the priority field).
    Returns dict with the priority field P[n_tor, n_pol] (matching ThicknessField.shape),
    the target (theta0, phi0) at the peak (deg), and R0.
    """
    phi, theta, w, R0 = importance_angles(npz, R0=R0)
    tor = np.radians(np.asarray(toroidal_angles_deg, float))
    pol = np.radians(np.asarray(poloidal_angles_deg, float))
    TOR, POL = np.meshgrid(tor, pol, indexing="ij")           # matches ThicknessField
    sig = np.radians(width_deg)

    def dperiodic(a, b):
        return (a - b + np.pi) % (2 * np.pi) - np.pi

    P = np.zeros(TOR.shape)
    # Gaussian splat each voxel onto the grid (periodic in both angles).
    for pj, tj, wj in zip(phi, theta, w):
        dph = dperiodic(TOR, pj)
        dth = dperiodic(POL, tj)
        P += wj * np.exp(-(dph ** 2 + dth ** 2) / (2 * sig ** 2))
    P /= P.max() + 1e-300

    k = np.unravel_index(np.argmax(P), P.shape)
    phi0 = float(np.degrees(tor[k[0]]))
    theta0 = float(np.degrees(pol[k[1]]))
    return dict(priority=P, phi0_deg=phi0, theta0_deg=theta0, R0=R0,
                TOR=TOR, POL=POL)


def coil_angles(coil_centroid, R0):
    """(phi, theta) of a coil centroid (cm) in the same convention, degrees."""
    x, y, z = coil_centroid
    R = np.hypot(x, y)
    return float(np.degrees(np.arctan2(y, x))), float(np.degrees(np.arctan2(z, R - R0)))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Adjoint importance -> shield placement priority")
    ap.add_argument("npz", help="adjoint_importance_*.npz map")
    ap.add_argument("--fluxmap", default=None,
                    help="VMEC fluxmap npz for the plasma source S(r); forms the "
                         "contributon C=S*psi_dagger (STRONGLY recommended -- bare flux "
                         "peaks in the central vacuum and mis-attributes).")
    ap.add_argument("--fluxmap-scale", type=float, default=100.0, help="m->cm = 100")
    ap.add_argument("--coil-centroid", type=float, nargs=3, default=None)
    ap.add_argument("--R0", type=float, default=None)
    ap.add_argument("--n-tor", type=int, default=48)
    ap.add_argument("--n-pol", type=int, default=48)
    args = ap.parse_args()

    if args.fluxmap is not None:
        field = contributon(args.npz, args.fluxmap, scale=args.fluxmap_scale)
        print("  using CONTRIBUTON C = S(r) * psi_dagger(r)")
    else:
        field = args.npz
        print("  WARNING: using BARE adjoint flux (no --fluxmap); peaks in central "
              "vacuum and mis-attributes. Pass --fluxmap for the physical contributon.")

    tor = np.linspace(0, 360, args.n_tor, endpoint=False)
    pol = np.linspace(0, 360, args.n_pol, endpoint=False)
    res = placement_priority(field, tor, pol, R0=args.R0)
    print(f"  R0 (auto) = {res['R0']:.1f} cm")
    print(f"  placement target (peak priority): phi0={res['phi0_deg']:.1f} deg  "
          f"theta0={res['theta0_deg']:.1f} deg")
    if args.coil_centroid is not None:
        cphi, cth = coil_angles(args.coil_centroid, res["R0"])
        dphi = ((res["phi0_deg"] - cphi + 180) % 360) - 180
        print(f"  coil centroid angles:            phi={cphi:.1f} deg  theta={cth:.1f} deg")
        print(f"  |phi0 - phi_coil| = {abs(dphi):.1f} deg  -> "
              f"{'MATCHES coil toroidal angle (placement points at coil)' if abs(dphi) < 45 else 'MISMATCH (!)'}")
