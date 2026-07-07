#!/usr/bin/env python
"""Exploratory, FULLY LOCAL study: is the SPF wall-load PEAKING FACTOR predictable
from cheap field/geometry statistics across QUASR quasisymmetric devices?

No cluster, no OpenMC transport. Everything here is an ANALYTIC free-streaming
integral plus the cheap coherence/direction-tensor statistics that already drive the
sweep. The point is to see whether the SPF-induced change in wall peaking factor
(perpendicular / parallel emission vs unpolarized) tracks any cheap predictor
(symmetry class, S_phi, nfp, aspect ratio, direction-tensor spectrum).

Physics recap (headless P2 kernel, rate normalized out):
  A neutron born at source point x with local field direction b_hat(x) is emitted
  with angular density  w(theta_B) = 1 + a2 * P2(cos theta_B),  P2(c)=1.5 c^2-0.5,
  theta_B = angle(emission dir, b_hat).  Because integral of P2 over the sphere is 0,
  every mode emits the SAME number of neutrons -- so comparing peaking factors across
  modes isolates DIRECTION, not rate.  a2 = 0 (unpolarized), a2 = -1 (perpendicular,
  the sin^2 A-mode shape), a2 = +1 (parallel, the 1/4+3/4cos^2 B/C-mode shape).

Free-streaming wall load at wall point r_w (no scattering, no occlusion -- convex
free-streaming limit, matching Schwartz's caveats):
  Q(r_w) = sum_x  s(x) * incid(x,r_w) * [1 + a2 P2(n.b_hat(x))] / |r_w - x|^2
  n = (r_w - x)/|r_w - x|                 (emission direction, source -> wall)
  incid = max(0, n . w_out(r_w))          (cosine of incidence on the wall; the
                                           genuine NWL current factor, forward faces)
PEAK FACTOR = max(Q) / area_weighted_mean(Q)  over the wall mesh.

Cheap stats reuse coherence_metrics (C, S_phi, lambda_phi, direction tensor evals)
with b_hat from the REAL coils (Biot-Savart, CoilField) sampled on the (1-rho^2)-
weighted fusion source volume -- exactly the source the compiled neutron source uses.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))

import coherence_metrics as cm  # noqa: E402
import quasr_loader as ql  # noqa: E402
from biotsavart_field import CoilField  # noqa: E402

# Device set (QA helicity=0, QH helicity!=0). Skip any that fail to load.
# The first five QA + five QH are the requested set; 166515 (nfp=1, near-axisymmetric
# anchor) and 806854 (nfp=3) are the two extra local serials, both QA.
QA_IDS = [803097, 886079, 923602, 932746, 59509, 166515, 806854]
QH_IDS = [1960314, 1190023, 2588494, 2630694, 2632243]
ALL_IDS = QA_IDS + QH_IDS

MODES = {"unpol": 0.0, "perp": -1.0, "par": +1.0}  # a2 per mode

WALL_GAP_FRAC = 0.3  # conformal wall = LCFS pushed out by 0.3*a_minor along its
                     # OUTWARD NORMAL (uniform gap => wall minor radius ~1.3*a).
                     # A fixed normal offset (not self-similar rho-scaling) is used
                     # because self-similar expansion of a strongly-shaped QH boundary
                     # folds the surface to within ~0.004*a of the source, producing a
                     # pure 1/r^2 discretization spike (one source sample sets 98% of
                     # the peak). The normal offset keeps a uniform source<->wall gap.


# --------------------------------------------------------------------------- #
# Source sampling + cheap field/geometry stats
# --------------------------------------------------------------------------- #
def sample_source(dev, n_rho=8, n_theta=16, n_phi=64, rho_max=0.85):
    """Source points (xyz), (1-rho^2) volume weights, and b_hat from the real coils.

    b_hat via Biot-Savart from the loaded QUASR filament coils (CoilField), i.e. the
    exact coil field direction on the fusion source volume."""
    xyz, w, rho = dev.source_sample(n_rho=n_rho, n_theta=n_theta, n_phi=n_phi,
                                    rho_max=rho_max)
    cf = CoilField(dev.coils, dev.meta)
    bhat = cf.bhat(xyz)  # (N,3) unit, Cartesian
    return xyz, w, rho, bhat


def cheap_stats(xyz, w, bhat, meta, dev):
    """All cheap predictors: coherence C, S_phi, lambda_phi, direction-tensor evals,
    plus nfp / aspect / class / boundary-derived R0,a."""
    m = cm.coherence_metrics(bhat, pos=xyz, weights=w, frame="cylindrical")
    R0b, a_b = boundary_R0_a(dev)
    return dict(
        C=m["C"], S_phi=m["S_phi"], lambda_phi=m["lambda_phi"],
        evals=[float(x) for x in m["tensor_evals"]],
        anisotropy=m["anisotropy"], reversal_frac=m["reversal_frac"],
        circular_std=m["circular_std"],
        nfp=int(meta["nfp"]), qs_class=meta["symmetry_class"],
        aspect_cat=meta.get("aspect"),
        R0=R0b, a_minor=a_b, aspect_boundary=(R0b / a_b if a_b else None),
    )


def boundary_R0_a(dev, ntheta=128, nphi=128):
    """R0/a from the LCFS: R0 = 0.5(Rmax+Rmin), a = 0.5(Rmax-Rmin) at rho=1."""
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, _ = dev.device.RZ(TH, PH, 1.0)
    Rmax, Rmin = float(R.max()), float(R.min())
    return 0.5 * (Rmax + Rmin), 0.5 * (Rmax - Rmin)


# --------------------------------------------------------------------------- #
# Conformal wall mesh + outward normals + area element
# --------------------------------------------------------------------------- #
def _surface_normals_area(P, dth, dph, R, Z, ph):
    """Outward unit normals + area element for a periodic (theta,phi) surface grid P.
    Normal oriented outward relative to the per-phi cross-section centroid."""
    nphi = P.shape[1]
    dPt = (np.roll(P, -1, axis=0) - np.roll(P, 1, axis=0)) / (2 * dth)
    dPp = (np.roll(P, -1, axis=1) - np.roll(P, 1, axis=1)) / (2 * dph)
    nrm = np.cross(dPt, dPp)                              # area-scaled normal
    dA = np.linalg.norm(nrm, axis=-1) * dth * dph
    nhat = nrm / np.linalg.norm(nrm, axis=-1, keepdims=True)
    Rc = R.mean(axis=0, keepdims=True); Zc = Z.mean(axis=0, keepdims=True)
    cen = np.stack([Rc * np.cos(ph[None, :]), Rc * np.sin(ph[None, :]),
                    Zc * np.ones((1, nphi))], axis=-1)
    outward = np.sign(np.sum(nhat * (P - cen), axis=-1, keepdims=True))
    outward[outward == 0] = 1.0
    return nhat * outward, dA


def wall_mesh(dev, ntheta=72, nphi=120, gap_frac=WALL_GAP_FRAC):
    """Conformal wall = LCFS offset outward by gap_frac * a_minor along its normal.
    Returns wall points (M,3), outward unit normals (M,3), area weights (M,), shape.
    Guarantees a uniform source<->wall gap (no self-intersection / near-field spike)."""
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")          # (nt, np)
    dth = th[1] - th[0]; dph = ph[1] - ph[0]

    # LCFS points + outward normals
    R0s, Z0s = dev.device.RZ(TH, PH, 1.0)
    P0 = np.stack([R0s * np.cos(PH), R0s * np.sin(PH), Z0s], axis=-1)
    n0, _ = _surface_normals_area(P0, dth, dph, R0s, Z0s, ph)

    a = boundary_R0_a(dev)[1]
    gap = gap_frac * a
    P = P0 + gap * n0                                     # offset wall points
    # recompute R,Z and outward normals/area ON the offset surface
    Rw = np.hypot(P[..., 0], P[..., 1]); Zw = P[..., 2]
    nhat, dA = _surface_normals_area(P, dth, dph, Rw, Zw, ph)

    return (P.reshape(-1, 3), nhat.reshape(-1, 3), dA.reshape(-1),
            (ntheta, nphi))


# --------------------------------------------------------------------------- #
# Analytic free-streaming wall load (all three modes share the geometry)
# --------------------------------------------------------------------------- #
def wall_load(src_xyz, src_w, src_bhat, wall_pts, wall_nhat, block=256, eps=1e-30):
    """Q_mode(r_w) for unpol/perp/par. Returns dict mode -> (M,) wall-load array.

    Vectorized over source points; chunked over wall points (block at a time) to
    bound memory. incid = max(0, n.w_out) is the wall-incidence current factor."""
    N = src_xyz.shape[0]; M = wall_pts.shape[0]
    out = {k: np.empty(M) for k in MODES}
    for i0 in range(0, M, block):
        i1 = min(i0 + block, M)
        rw = wall_pts[i0:i1]                              # (m,3)
        nw = wall_nhat[i0:i1]                             # (m,3)
        d = rw[:, None, :] - src_xyz[None, :, :]          # (m,N,3) source->wall
        r2 = np.sum(d * d, axis=2)                        # (m,N)
        r = np.sqrt(r2)
        n = d / (r[..., None] + eps)                      # emission unit dir
        incid = np.clip(np.sum(n * nw[:, None, :], axis=2), 0.0, None)  # (m,N)
        geom = src_w[None, :] * incid / (r2 + eps)        # (m,N) mode-independent
        cB = np.sum(n * src_bhat[None, :, :], axis=2)     # (m,N) = n.b_hat
        P2 = 1.5 * cB * cB - 0.5
        for mode, a2 in MODES.items():
            out[mode][i0:i1] = np.sum(geom * (1.0 + a2 * P2), axis=1)
    return out


def peak_factor(Q, dA):
    """max(Q) / area-weighted mean(Q)."""
    mean = float(np.sum(Q * dA) / np.sum(dA))
    return float(Q.max() / mean), mean


# --------------------------------------------------------------------------- #
# Per-device driver
# --------------------------------------------------------------------------- #
# Production resolution (see convergence table in the findings note). QA is fully
# converged here; QH PF-ratios are within ~7% of the finest grid and their sign/
# ordering is resolution-invariant.
SRC_RES = dict(n_rho=10, n_theta=24, n_phi=96)
WALL_RES = dict(ntheta=96, nphi=160)


def analyze_device(ID, src_res=None, wall_res=None, block=128, keep_maps=False):
    src_res = src_res or SRC_RES
    wall_res = wall_res or WALL_RES
    dev = ql.load_device(ID)
    xyz, w, rho, bhat = sample_source(dev, **src_res)
    stats = cheap_stats(xyz, w, bhat, dev.meta, dev)

    wpts, wnhat, wdA, wshape = wall_mesh(dev, **wall_res)
    Q = wall_load(xyz, w, bhat, wpts, wnhat, block=block)
    pf = {}; mean = {}
    for mode in MODES:
        pf[mode], mean[mode] = peak_factor(Q[mode], wdA)

    rec = dict(ID=int(ID), **stats,
               PF_unpol=pf["unpol"], PF_perp=pf["perp"], PF_par=pf["par"],
               PF_perp_over_unpol=pf["perp"] / pf["unpol"],
               PF_par_over_unpol=pf["par"] / pf["unpol"],
               n_src=int(xyz.shape[0]), n_wall=int(wpts.shape[0]))

    # --- directional (field-intrinsic) observables: NOT the local peak. These measure how the
    # emission bias redistributes the load, which should track the field (S_phi) far better than
    # the geometry-dominated peak extremum. ---
    Rw = np.hypot(wpts[:, 0], wpts[:, 1])
    inb = Rw < stats["R0"]
    tot = float(np.sum(wdA))
    Qu_mean = float(np.sum(Q["unpol"] * wdA) / tot)

    def io_asym(Qm):                 # inboard/outboard asymmetry (Schwartz's +-43% is this shift)
        qin = float(np.sum(Qm[inb] * wdA[inb])); qout = float(np.sum(Qm[~inb] * wdA[~inb]))
        return (qin - qout) / (qin + qout + 1e-30)

    def rms_redist(Qm):              # global RMS fractional redistribution vs unpol
        return float(np.sqrt(np.sum((Qm - Q["unpol"]) ** 2 * wdA) / tot) / (Qu_mean + 1e-30))

    au, ap, aq = io_asym(Q["unpol"]), io_asym(Q["perp"]), io_asym(Q["par"])
    rec.update(io_asym_unpol=au, io_shift_perp=ap - au, io_shift_par=aq - au,
               redist_perp=rms_redist(Q["perp"]), redist_par=rms_redist(Q["par"]))
    maps = None
    if keep_maps:
        maps = dict(shape=wshape,
                    Q_unpol=Q["unpol"], Q_perp=Q["perp"], Q_par=Q["par"],
                    dA=wdA, R0=stats["R0"])
    return rec, maps


def run_all(ids=ALL_IDS, keep_map_ids=()):
    recs = []; maps = {}
    for ID in ids:
        try:
            rec, mp = analyze_device(ID, keep_maps=(ID in keep_map_ids))
        except Exception as e:  # noqa: BLE001 -- skip unloadable devices, report
            print(f"[skip] device {ID}: {type(e).__name__}: {e}")
            continue
        recs.append(rec)
        if mp is not None:
            maps[ID] = mp
        print(f"[ok] {ID} {rec['qs_class']} nfp={rec['nfp']} "
              f"S_phi={rec['S_phi']:.3f} C={rec['C']:.3f} "
              f"PF_unpol={rec['PF_unpol']:.3f} "
              f"perp/unpol={rec['PF_perp_over_unpol']:.4f} "
              f"par/unpol={rec['PF_par_over_unpol']:.4f}")
    return recs, maps


if __name__ == "__main__":
    out_dir = HERE / "sweep_out"
    out_dir.mkdir(exist_ok=True)
    recs, _ = run_all()
    (out_dir / "geometry_peaking.json").write_text(json.dumps(recs, indent=2))
    print(f"\nwrote {out_dir / 'geometry_peaking.json'}  ({len(recs)} devices)")
