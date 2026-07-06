#!/usr/bin/env python
"""Field-direction coherence metrics -- the CHEAP predictors of the SPF benefit.

These use field direction + source weight ONLY (no neutronics), and are the
hypothesized predictors of the directional efficiency eta. The central quantity:

    C = | integral s(x) bhat(x) dV | / integral s(x) dV        in [0,1]

is the source-weighted mean RESULTANT LENGTH (R-bar) of the field-direction
distribution on the unit sphere -- the magnitude of its first moment
(Mardia-style directional statistics). C = 1 => perfectly aligned field direction
(tokamak limit, maximal coherent SPF steering); C -> 0 => dispersed (steering from
different source points cancels at a fixed detector).

FRAME (decisive, documented in docs/THEORY.md): bhat is expressed in the LOCAL
CYLINDRICAL basis (e_R, e_phi, e_Z) before the moment is taken. A purely toroidal
field is (0,1,0) everywhere in that basis, so the tokamak limit calibrates to
C = 1; a stellarator's helical/pitch variation across the source drives C < 1.
(In the fixed lab xyz frame a toroidal field winds around the torus and integrates
to ~0, so lab-frame C could not be the discriminator -- see THEORY.md.)

Also computed:
  * the source-weighted DIRECTION TENSOR  T = <bhat bhat^T>_s  (3x3, symmetric,
    trace 1); its eigenvalues capture the STRUCTURE of the spread (a planar fan vs
    an isotropic smear give the same |first moment| but different eigenvalue
    spectra), a stronger predictor than the scalar C when symmetry classes split.
  * the source-weighted angular standard deviation of bhat about its mean
    direction, and the directional-statistics circular std sqrt(-2 ln C).

Use the SAME source weights s(x) the neutronics source uses (see source_weights)
so C and eta are consistent by construction.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# Frame handling
# --------------------------------------------------------------------------- #
def to_cylindrical(bhat_xyz, pos_xyz):
    """Express Cartesian unit vectors in the LOCAL cylindrical basis at each point.
    bhat_xyz, pos_xyz: (N,3). Returns (N,3) = (b.e_R, b.e_phi, b.e_Z).
    A rotation of basis preserves the norm, so unit vectors stay unit."""
    bhat_xyz = np.asarray(bhat_xyz, float)
    pos_xyz = np.asarray(pos_xyz, float)
    phi = np.arctan2(pos_xyz[:, 1], pos_xyz[:, 0])
    c, s = np.cos(phi), np.sin(phi)
    bx, by, bz = bhat_xyz[:, 0], bhat_xyz[:, 1], bhat_xyz[:, 2]
    bR = bx * c + by * s
    bP = -bx * s + by * c
    return np.stack([bR, bP, bz], axis=1)


def _prep(bhat, pos, weights, frame):
    b = np.asarray(bhat, float)
    if b.ndim != 2 or b.shape[1] != 3:
        raise ValueError("bhat must be (N,3)")
    nrm = np.linalg.norm(b, axis=1, keepdims=True)
    if np.any(nrm == 0):
        raise ValueError("zero-length bhat sample")
    b = b / nrm  # hard-normalize (audit enforces |b|=1 upstream; belt-and-suspenders)
    if frame == "cylindrical":
        if pos is None:
            raise ValueError("frame='cylindrical' requires pos (N,3)")
        b = to_cylindrical(b, pos)
    elif frame != "lab":
        raise ValueError(f"unknown frame {frame!r} (use 'cylindrical' or 'lab')")
    if weights is None:
        w = np.ones(len(b))
    else:
        w = np.asarray(weights, float)
        if w.shape != (len(b),):
            raise ValueError("weights must be (N,)")
        if np.any(w < 0):
            raise ValueError("negative source weight")
    tot = w.sum()
    if tot <= 0:
        raise ValueError("source weights sum to <= 0")
    return b, w / tot


# --------------------------------------------------------------------------- #
# The metrics
# --------------------------------------------------------------------------- #
def coherence_C(bhat, pos=None, weights=None, frame="cylindrical"):
    """Scalar order parameter C = |<bhat>_s| (mean resultant length) in [0,1]."""
    b, w = _prep(bhat, pos, weights, frame)
    return float(np.linalg.norm((w[:, None] * b).sum(axis=0)))


def direction_tensor(bhat, pos=None, weights=None, frame="cylindrical"):
    """Source-weighted direction tensor T = <bhat bhat^T>_s and its eigenvalues.
    Returns (T (3,3), evals descending, evecs columns). trace(T) = 1 exactly."""
    b, w = _prep(bhat, pos, weights, frame)
    T = np.einsum("i,ij,ik->jk", w, b, b)
    T = 0.5 * (T + T.T)  # symmetrize against roundoff
    evals, evecs = np.linalg.eigh(T)
    order = np.argsort(evals)[::-1]
    return T, evals[order], evecs[:, order]


def coherence_metrics(bhat, pos=None, weights=None, frame="cylindrical"):
    """All cheap predictors in one call. Returns a dict:
        C            : mean resultant length |<b>| in [0,1]  (FIRST moment / polar)
        mean_dir     : unit mean field direction (in `frame`)
        tensor_evals : (l1>=l2>=l3), sum = 1  (structure of the spread)
        anisotropy   : l1 - l2   (BIAXIALITY: planar-fan vs isotropic-smear; this is
                       NOT the coherence-loss scalar -- use S_phi for that, see below)
        lambda_phi   : <b_phi^2>_s = T[1,1], the direction-tensor weight along the
                       toroidal axis e_phi (the frame's 2nd axis).  SECOND moment.
        S_phi        : (3*lambda_phi - 1)/2, the NEMATIC (P2 / director) order about
                       e_phi and the PHYSICALLY-MATCHED predictor of SPF steering. The
                       emission kernel w ~ 1 + a2 P2(cos theta_B) is EVEN in b_hat
                       (headless / quadrupolar), so every wall observable is a linear
                       functional of the SECOND moment <b b^T> ONLY and is INDEPENDENT
                       of the first moment C. C tracks eta only via its correlation
                       with S_phi, which holds for a co-toroidal cap where
                       C^2 <~ lambda_phi <~ C. See docs/THEORY.md (parity argument).
        reversal_frac: source weight with b_phi < 0 (reversed toroidal sense). 0 for a
                       cap (C and S_phi monotone-locked); >0 only where the field
                       reverses, and THERE C under-predicts while S_phi still tracks eta.
        angular_std  : source-weighted RMS angle [rad] of bhat about mean_dir
        circular_std : directional-statistics sqrt(-2 ln C) [rad]

    lambda_phi/S_phi/reversal_frac reference the frame's 2nd axis (the toroidal e_phi
    in the cylindrical frame -- the tokamak C=1 calibration axis), so they are only
    physically meaningful with frame='cylindrical'.
    """
    b, w = _prep(bhat, pos, weights, frame)
    mean_vec = (w[:, None] * b).sum(axis=0)
    C = float(np.linalg.norm(mean_vec))
    mean_dir = mean_vec / (C if C > 0 else 1.0)
    T = np.einsum("i,ij,ik->jk", w, b, b)
    T = 0.5 * (T + T.T)
    evals = np.sort(np.linalg.eigvalsh(T))[::-1]
    lam_phi = float(T[1, 1])                           # <b_phi^2>_s (toroidal axis)
    S_phi = float((3.0 * lam_phi - 1.0) / 2.0)         # nematic order about e_phi
    reversal_frac = float(w[b[:, 1] < 0.0].sum())      # source weight with b_phi < 0
    cosang = np.clip(b @ mean_dir, -1.0, 1.0)
    ang = np.arccos(cosang)
    angular_std = float(np.sqrt((w * ang ** 2).sum()))
    circular_std = float(np.sqrt(max(-2.0 * np.log(max(C, 1e-300)), 0.0)))
    return dict(C=C, mean_dir=mean_dir, tensor_evals=evals,
                anisotropy=float(evals[0] - evals[1]),
                lambda_phi=lam_phi, S_phi=S_phi, reversal_frac=reversal_frac,
                angular_std=angular_std, circular_std=circular_std)


# --------------------------------------------------------------------------- #
# Source weighting -- MUST match what the neutronics CompiledSource uses
# --------------------------------------------------------------------------- #
def source_weights(rho, model="one_minus_rho2"):
    """Fusion-source density proxy s(x) as a function of the flux label rho in
    [0,1]. The compiled source births on (1 - rho^2) flux surfaces (run_conformal
    'shape=plasma'), so the DEFAULT here is the same (1 - rho^2) profile -- keep C
    and the neutronics source consistent BY CONSTRUCTION. `n2sigmav` is a hotter
    core proxy for n^2<sigma v> ~ (1 - rho^2)^2 for later use.
    """
    rho = np.asarray(rho, float)
    if model == "one_minus_rho2":
        w = 1.0 - rho ** 2
    elif model == "n2sigmav":
        w = (1.0 - rho ** 2) ** 2
    elif model == "uniform":
        w = np.ones_like(rho)
    else:
        raise ValueError(f"unknown source model {model!r}")
    return np.clip(w, 0.0, None)
