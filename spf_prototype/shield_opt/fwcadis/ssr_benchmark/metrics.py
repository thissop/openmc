"""Small, dependency-free metrics for SSR variance-reduction comparisons."""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable, Sequence


def _same_length(first: Sequence[float], second: Sequence[float]) -> None:
    if len(first) != len(second):
        raise ValueError("metric inputs must have the same length")


def relative_errors(
    means: Sequence[float], standard_deviations: Sequence[float]
) -> list[float]:
    """Return standard deviation divided by absolute mean.

    Zero-valued bins are unresolved and receive positive infinity rather than
    an apparently favorable zero uncertainty.
    """

    _same_length(means, standard_deviations)
    return [
        abs(float(sigma) / float(mean)) if mean != 0.0 else math.inf
        for mean, sigma in zip(means, standard_deviations)
    ]


def figures_of_merit(
    errors: Sequence[float], simulation_seconds: float
) -> list[float]:
    """Return the conventional Monte Carlo FOM, ``1 / (R**2 * T)``."""

    if simulation_seconds <= 0.0:
        raise ValueError("simulation_seconds must be positive")
    return [
        1.0 / (error * error * simulation_seconds)
        if math.isfinite(error) and error > 0.0
        else 0.0
        for error in errors
    ]


def normalized_pattern(values: Sequence[float]) -> list[float]:
    """Normalize a non-negative response vector to unit sum."""

    if any(value < 0.0 for value in values):
        raise ValueError("response values must be non-negative")
    total = math.fsum(float(value) for value in values)
    if total <= 0.0:
        return [math.nan for _ in values]
    return [float(value) / total for value in values]


def rms_fractional_bias(
    estimate: Sequence[float], reference: Sequence[float]
) -> float:
    """Root-mean-square fractional difference from a positive reference."""

    _same_length(estimate, reference)
    if not estimate:
        return math.nan
    if any(value <= 0.0 for value in reference):
        raise ValueError("reference values must be positive")
    terms = [
        ((float(value) - float(target)) / float(target)) ** 2
        for value, target in zip(estimate, reference)
    ]
    return math.sqrt(math.fsum(terms) / len(terms))


def median_positive_finite(values: Iterable[float | None]) -> float:
    """Median of positive finite values, or zero when none are available."""

    usable = [
        float(value)
        for value in values
        if value is not None and value > 0.0 and math.isfinite(value)
    ]
    return statistics.median(usable) if usable else 0.0
