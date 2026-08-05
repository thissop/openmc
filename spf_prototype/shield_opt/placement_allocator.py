"""Multi-coil-aware greedy shield-placement allocator.

Turns a (theta, phi) PRIORITY field into a fixed-envelope shield-thickness field delta(theta,phi),
ready for ThicknessField.trade / .radial_build_dict -> DAGMC build. The priority field is whatever
ranking we trust -- a single coil's contributon, the scalar shielding-worth W, or (the point of
this module) a MULTI-COIL importance combining every coil, so the placement protects the whole
coil set rather than over-protecting one coil (the kill-shot's peak-migration failure: crushing
coil 20 pushed the peak to the unprotected coil 25).

Two allocation policies:
  * greedy   -- spend a FIXED material budget where it buys the most protection: fill the
                highest-priority cells to delta_max first, then partially fill the marginal cell so
                the added-shield volume hits the budget exactly. Bang-bang, material-optimal for a
                linear coverage objective. This is the design lever.
  * proportional -- delta = delta_max * (P/P_max)^gamma. Smooth, budget-free; a soft baseline.

All output delta is clipped to [0, tf.delta_max] (= t_breeder0 - t_breeder_min), so the breeder
floor and the fixed envelope are respected BY CONSTRUCTION (ThicknessField.trade then just holds).
The allocator is coil-agnostic; 'multi-coil-aware' means you feed it combine_priorities(...).
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------------------------
# Multi-coil combination (at the (theta,phi) grid level -- mirrors multicoil_importance but for
# already-projected priority fields).
# ----------------------------------------------------------------------------------------------
def combine_priorities(fields, mode="sum", weights=None):
    """Combine a stack of per-coil (theta,phi) priority grids into one.

    fields : (K, nTor, nPol) or list of K (nTor,nPol) arrays -- one per coil.
    mode   : 'sum'    -> protect total coil load (every coil equally),
             'max'    -> protect the worst coil at each (theta,phi) (envelope),
             'weighted' (needs weights, len K) -> e.g. flux-weight so hot coils dominate -> peak.
    Returns a single (nTor,nPol) field, normalized to max 1 (priority is a ranking, scale-free)."""
    F = np.asarray(fields, float)
    if F.ndim != 3:
        raise ValueError(f"fields must stack to (K,nTor,nPol); got {F.shape}")
    if mode == "sum":
        P = F.sum(axis=0)
    elif mode == "max":
        P = F.max(axis=0)
    elif mode == "weighted":
        if weights is None or len(weights) != F.shape[0]:
            raise ValueError("mode='weighted' needs weights of length K")
        w = np.asarray(weights, float)[:, None, None]
        P = (w * F).sum(axis=0)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return P / (P.max() + 1e-300)


# ----------------------------------------------------------------------------------------------
# Allocation policies
# ----------------------------------------------------------------------------------------------
def _cell_area(priority, area):
    if area is None:
        return np.ones_like(priority, float)   # ThicknessField's uniform-cell convention
    a = np.asarray(area, float)
    if a.shape != priority.shape:
        raise ValueError(f"area shape {a.shape} != priority shape {priority.shape}")
    return a


def greedy_allocate(priority, tf, budget_frac=0.5, area=None):
    """Fixed-budget greedy allocation.

    Spend budget_frac of the maximum addable shield volume (sum over cells of delta_max*area) by
    filling the highest-priority cells to delta_max first; the marginal cell is partially filled so
    the total added volume equals the budget exactly. Returns delta (tf.shape) in [0, delta_max].

    budget_frac in (0,1]; budget_frac=1 -> full trade everywhere (== uniform max shield)."""
    P = np.asarray(priority, float)
    if P.shape != tf.shape:
        raise ValueError(f"priority shape {P.shape} != tf.shape {tf.shape}")
    dmax = float(tf.delta_max)
    A = _cell_area(P, area)
    budget = float(np.clip(budget_frac, 0.0, 1.0)) * dmax * float(A.sum())

    order = np.argsort(P, axis=None)[::-1]        # highest priority first (ties: arbitrary stable)
    delta = np.zeros(P.size, float)
    Af = A.ravel()
    spent = 0.0
    for idx in order:
        cell_full = dmax * Af[idx]
        if spent + cell_full <= budget + 1e-12:
            delta[idx] = dmax
            spent += cell_full
        else:                                     # partial fill to hit the budget exactly
            remaining = budget - spent
            if remaining > 0:
                delta[idx] = remaining / Af[idx]
                spent = budget
            break
    return np.clip(delta.reshape(P.shape), 0.0, dmax)


def proportional_allocate(priority, tf, gamma=1.0):
    """Smooth budget-free allocation: delta = delta_max * (P/P_max)^gamma, clipped to [0,delta_max].
    gamma>1 concentrates on the peak; gamma<1 spreads. A soft baseline vs the greedy lever."""
    P = np.asarray(priority, float)
    Pn = P / (P.max() + 1e-300)
    return np.clip(float(tf.delta_max) * Pn ** float(gamma), 0.0, float(tf.delta_max))


def smooth(delta, tf, passes=1):
    """Manufacturability smoothing: periodic 3x3 box blur (both angles periodic), re-clipped to
    [0, delta_max]. Reduces ThicknessField.smoothness (sharp steps) at a small coverage cost."""
    d = np.asarray(delta, float)
    for _ in range(int(passes)):
        d = (d
             + np.roll(d, 1, 0) + np.roll(d, -1, 0)
             + np.roll(d, 1, 1) + np.roll(d, -1, 1)) / 5.0
    return np.clip(d, 0.0, float(tf.delta_max))


# ----------------------------------------------------------------------------------------------
# Figure of merit
# ----------------------------------------------------------------------------------------------
def coverage(priority, delta, tf, area=None):
    """Priority-weighted coverage in [0,1]: how much of the priority mass got shielded, relative to
    filling every cell to delta_max. FoM for comparing allocations at equal budget -- a
    priority-aware allocation beats a uniform one here."""
    P = np.asarray(priority, float)
    A = _cell_area(P, area)
    num = float(np.sum(P * delta * A))
    den = float(np.sum(P * float(tf.delta_max) * A)) + 1e-300
    return num / den


def allocate(priority, tf, policy="greedy", **kw):
    """Dispatch: policy in {'greedy','proportional'}. Convenience wrapper."""
    if policy == "greedy":
        return greedy_allocate(priority, tf, **kw)
    if policy == "proportional":
        return proportional_allocate(priority, tf, **kw)
    raise ValueError(f"unknown policy {policy!r}")
