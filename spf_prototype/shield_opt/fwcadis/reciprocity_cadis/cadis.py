"""Construct manual CADIS weight-window bounds from a reciprocal flux map."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CadisBounds:
    lower: np.ndarray
    upper: np.ndarray
    relative_error: np.ndarray
    reliable: np.ndarray


def build_cadis_bounds(
    mean: np.ndarray,
    std_dev: np.ndarray,
    cell_volume: float | np.ndarray,
    relative_error_threshold: float = 0.5,
    upper_bound_ratio: float = 5.0,
) -> CadisBounds:
    """Apply OpenMC's CADIS algebra, rejecting noisy bins before normalization.

    OpenMC uses a volume-normalized importance, takes its reciprocal, and scales
    the largest reciprocal to 1/2.  This implementation first removes bins that
    will be ignored so one barely sampled bin cannot collapse every useful
    weight window.  Ignored bins use OpenMC's negative-bound convention.
    """
    mean = np.asarray(mean, dtype=float)
    std_dev = np.asarray(std_dev, dtype=float)
    if mean.shape != std_dev.shape:
        raise ValueError("mean and std_dev must have the same shape")
    if relative_error_threshold <= 0.0:
        raise ValueError("relative_error_threshold must be positive")
    if upper_bound_ratio <= 1.0:
        raise ValueError("upper_bound_ratio must be greater than one")

    volume = np.asarray(cell_volume, dtype=float)
    if np.any(~np.isfinite(volume)) or np.any(volume <= 0.0):
        raise ValueError("cell_volume must be finite and positive")

    relative_error = np.full(mean.shape, np.inf, dtype=float)
    sampled = np.isfinite(mean) & np.isfinite(std_dev) & (mean > 0.0) & (std_dev >= 0.0)
    relative_error[sampled] = std_dev[sampled] / mean[sampled]
    reliable = sampled & (relative_error <= relative_error_threshold)
    if not np.any(reliable):
        raise ValueError("No positive mesh bins satisfy the relative-error threshold")

    flux_density = mean / volume
    inverse_importance = np.zeros(mean.shape, dtype=float)
    inverse_importance[reliable] = 1.0 / flux_density[reliable]
    normalizer = 2.0 * np.max(inverse_importance[reliable])

    lower = np.full(mean.shape, -1.0, dtype=float)
    lower[reliable] = inverse_importance[reliable] / normalizer
    upper = upper_bound_ratio * lower
    return CadisBounds(lower, upper, relative_error, reliable)
