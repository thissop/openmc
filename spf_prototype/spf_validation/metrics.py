"""Accuracy / speed / memory metrics for the representation head-to-head.

Pure functions, unit-tested in test_fluxrep.py. Accuracy compares the grid representation
against the spectral one at random off-node query points (the spectral eval is the truth,
since the grid is filled FROM it); speed is evaluation throughput; memory is the stored
footprint each representation carries.
"""
from __future__ import annotations
import time
import numpy as np


def random_query(n, nfp, seed=0, s_min=1e-3):
    """n random interior points (s, theta, zeta), theta in [0,2pi), zeta in one field period."""
    rng = np.random.default_rng(seed)
    s = rng.uniform(s_min, 1.0, n)
    theta = rng.uniform(0.0, 2 * np.pi, n)
    zeta = rng.uniform(0.0, 2 * np.pi / nfp, n)
    return s, theta, zeta


def accuracy(spectral, grid, n=20000, nfp=1, seed=0):
    """RMS and max RELATIVE deviation of the grid vs spectral for R, Z, sqrt(g).

    R, Z are normalized by the minor-radius scale (max|R-R0|) so the relative error is
    meaningful near the axis; sqrt(g) by its own mean.
    """
    s, th, ze = random_query(n, nfp, seed)
    Rs, Zs, Gs = spectral.eval(s, th, ze, want_sqrtg=True)
    Rg, Zg, Gg = grid.eval(s, th, ze, want_sqrtg=True)
    R0 = np.mean(Rs)
    scaleR = np.max(np.abs(Rs - R0)) + 1e-30
    scaleZ = np.max(np.abs(Zs)) + 1e-30
    scaleG = np.mean(np.abs(Gs)) + 1e-30

    def stat(a, b, scale):
        d = (a - b) / scale
        return {"rms": float(np.sqrt(np.mean(d**2))), "max": float(np.max(np.abs(d)))}

    exact_pdf = np.maximum(Gs, 0.0)
    grid_pdf = np.maximum(Gg, 0.0)
    exact_norm = np.sum(exact_pdf)
    grid_norm = np.sum(grid_pdf)
    exact_pdf /= exact_norm
    grid_pdf /= grid_norm
    source_pdf = {
        "normalization_bias": float(grid_norm / exact_norm - 1.0),
        "total_variation": float(0.5 * np.sum(np.abs(grid_pdf - exact_pdf))),
        "hellinger": float(np.sqrt(
            0.5 * np.sum((np.sqrt(grid_pdf) - np.sqrt(exact_pdf))**2))),
    }

    return {
        "n_query": int(n),
        "R": stat(Rg, Rs, scaleR),
        "Z": stat(Zg, Zs, scaleZ),
        "sqrtg": stat(Gg, Gs, scaleG),
        "source_pdf": source_pdf,
        **(_bhat_accuracy(spectral, grid, s, th, ze, exact_pdf)
           if getattr(spectral, "has_field", False) else {}),
    }


def _weighted_quantile(values, weights, quantile):
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    return float(np.interp(quantile, cumulative, values))


def _bhat_accuracy(spectral, grid, s, theta, zeta, source_weights):
    exact = spectral.eval_bhat(s, theta, zeta)
    approximate = grid.eval_bhat(s, theta, zeta)
    # SPF angular laws are invariant under B-hat -> -B-hat, so score the
    # physically relevant unoriented-axis error.
    cosine = np.clip(np.abs(np.einsum("ij,ij->i", exact, approximate)),
                     0.0, 1.0)
    degrees = np.degrees(np.arccos(cosine))
    return {"bhat_degrees": {
        "median": float(np.median(degrees)),
        "p95": float(np.percentile(degrees, 95.0)),
        "p99": float(np.percentile(degrees, 99.0)),
        "max": float(np.max(degrees)),
        "source_weighted_rms": float(np.sqrt(
            np.sum(source_weights * degrees**2))),
        "source_weighted_p95": _weighted_quantile(
            degrees, source_weights, 0.95),
        "source_weighted_p99": _weighted_quantile(
            degrees, source_weights, 0.99),
    }}


def throughput(rep, n=100000, nfp=1, seed=1, repeats=3, want_sqrtg=True):
    """Evaluation throughput (points/second): best of `repeats` timed runs of n points."""
    s, th, ze = random_query(n, nfp, seed)
    best = np.inf
    for _ in range(repeats):
        t0 = time.perf_counter()
        rep.eval(s, th, ze, want_sqrtg=want_sqrtg)
        best = min(best, time.perf_counter() - t0)
    return {"points": int(n), "seconds": float(best), "points_per_sec": float(n / best)}


def bhat_throughput(rep, n=100000, nfp=1, seed=1, repeats=3):
    """Magnetic-field direction evaluation throughput (points/second)."""
    s, th, ze = random_query(n, nfp, seed)
    best = np.inf
    for _ in range(repeats):
        t0 = time.perf_counter()
        rep.eval_bhat(s, th, ze)
        best = min(best, time.perf_counter() - t0)
    return {"points": int(n), "seconds": float(best),
            "points_per_sec": float(n / best)}


def footprint(rep):
    """Stored memory footprint of the representation, in bytes and MiB."""
    b = int(getattr(rep, "nbytes", 0))
    return {"bytes": b, "MiB": round(b / 2**20, 3)}
