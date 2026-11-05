"""Shared one-period geometry helpers for FW-CADIS experiments."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read_makegrid(path: str | Path) -> list[np.ndarray]:
    """Read MAKEGRID filament coils, including the terminating point of each coil."""
    lines = Path(path).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "mirror NIL") + 1
    end = next(i for i, line in enumerate(lines[start:], start) if line.strip() == "end")

    coils: list[np.ndarray] = []
    current: list[list[float]] = []
    for line in lines[start:end]:
        fields = line.split()
        if len(fields) < 4:
            continue
        current.append([float(value) for value in fields[:4]])
        if len(fields) > 4:
            coils.append(np.asarray(current, dtype=float))
            current = []
    if current:
        raise ValueError(f"Unterminated coil in {path}")
    if not coils:
        raise ValueError(f"No coils found in {path}")
    return coils


def sample_coil_points(
    coils: list[np.ndarray], points_per_coil: int, scale: float = 100.0
) -> np.ndarray:
    """Sample equal-arc-length interior points from each filament polyline."""
    if points_per_coil < 1:
        raise ValueError("points_per_coil must be positive")

    sampled = []
    for coil in coils:
        xyz = np.asarray(coil, dtype=float)[:, :3] * scale
        segment_vectors = np.diff(xyz, axis=0)
        segment_lengths = np.linalg.norm(segment_vectors, axis=1)
        keep = segment_lengths > 0.0
        segment_vectors = segment_vectors[keep]
        segment_lengths = segment_lengths[keep]
        starts = xyz[:-1][keep]
        if not segment_lengths.size:
            raise ValueError("Coil contains no nonzero-length segments")

        cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
        targets = (np.arange(points_per_coil) + 0.5) * cumulative[-1] / points_per_coil
        indices = np.searchsorted(cumulative, targets, side="right") - 1
        indices = np.minimum(indices, len(segment_lengths) - 1)
        fractions = (targets - cumulative[indices]) / segment_lengths[indices]
        sampled.append(starts[indices] + fractions[:, None] * segment_vectors[indices])

    return np.concatenate(sampled, axis=0)


def stellarator_source_point_cloud(
    stem: str | Path, max_points: int = 20_000
) -> tuple[np.ndarray, np.ndarray]:
    """Approximate ``StellaratorSource`` births as an IndependentSource cloud.

    Positions use the fluxmap grid and the native source weight
    ``max(1-rho**2, 0) * sqrt(g)``. If needed, deterministic systematic
    resampling bounds the XML size. Toroidal positions are folded into the
    first field period for the one-period DAGMC model.
    """
    if max_points < 1:
        raise ValueError("max_points must be positive")
    stem = Path(stem)
    if stem.suffix in {".meta", ".bin"}:
        stem = stem.with_suffix("")

    metadata = {}
    for line in stem.with_suffix(".meta").read_text().splitlines():
        if line.strip():
            key, value = line.split()
            metadata[key] = value
    if metadata.get("magic") != "spf_fluxmap_v1":
        raise ValueError(f"Invalid fluxmap metadata in {stem}.meta")
    nr = int(metadata["nR"])
    ntheta = int(metadata["ntheta"])
    nzeta = int(metadata["nzeta"])
    nfp = int(metadata.get("nfp", "1"))

    raw = np.fromfile(stem.with_suffix(".bin"), dtype="<f8")
    offset = 0

    def take(size: int) -> np.ndarray:
        nonlocal offset
        values = raw[offset:offset + size]
        offset += size
        if values.size != size:
            raise ValueError(f"Short fluxmap payload in {stem}.bin")
        return values

    rho = take(nr)
    take(ntheta)  # theta grid
    take(nzeta)   # zeta grid
    n3 = nr * ntheta * nzeta

    def take3() -> np.ndarray:
        return take(n3).reshape(nr, ntheta, nzeta)

    sqrtg = take3()
    major_radius = take3()
    vertical = take3()
    phi = take3()
    take3()  # BR
    take3()  # Bphi
    take3()  # BZ
    if offset != raw.size:
        raise ValueError(f"Trailing data in {stem}.bin")

    period = 2.0 * np.pi / nfp
    folded_phi = np.mod(phi, period)
    angular_epsilon = 1.0e-10
    folded_phi = np.clip(folded_phi, angular_epsilon, period - angular_epsilon)
    positions = np.column_stack([
        (major_radius * np.cos(folded_phi)).ravel(),
        (major_radius * np.sin(folded_phi)).ravel(),
        vertical.ravel(),
    ])
    weights = (
        np.maximum(1.0 - rho * rho, 0.0)[:, None, None] * np.abs(sqrtg)
    ).ravel()
    valid = np.all(np.isfinite(positions), axis=1) & np.isfinite(weights) & (weights > 0.0)
    positions = positions[valid]
    weights = weights[valid]
    if not weights.size:
        raise ValueError("Fluxmap has no positive finite source weights")

    if weights.size > max_points:
        cumulative = np.cumsum(weights)
        targets = (np.arange(max_points) + 0.5) * cumulative[-1] / max_points
        indices = np.searchsorted(cumulative, targets, side="left")
        positions = positions[indices]
        weights = np.ones(max_points)
    return positions, weights


def build_materials(openmc, li6: float = 60.0):
    """Build the material names expected by ``dagmc_qh_op_uniform.h5m``."""
    tungsten = openmc.Material()
    tungsten.add_element("W", 1.0)
    tungsten.set_density("g/cm3", 19.3)

    rafm = openmc.Material()
    for element, fraction in [
        ("Fe", 0.885), ("Cr", 0.09), ("W", 0.011), ("Mn", 0.006),
        ("V", 0.002), ("C", 0.006),
    ]:
        rafm.add_element(element, fraction, "wo")
    rafm.set_density("g/cm3", 7.9)
    first_wall = openmc.Material.mix_materials(
        [tungsten, rafm], [0.2 / 3.2, 3.0 / 3.2], "vo"
    )
    first_wall.name = "first_wall"

    multiplier = openmc.Material(name="multiplier")
    multiplier.add_element("Be", 1.0)
    multiplier.set_density("g/cm3", 1.85)

    breeder = openmc.Material(name="breeder")
    breeder.add_element(
        "Li", 2.0, enrichment=li6, enrichment_target="Li6", enrichment_type="ao"
    )
    breeder.add_element("Be", 1.0)
    breeder.add_element("F", 4.0)
    breeder.set_density("g/cm3", 1.94)

    back_wall = openmc.Material(name="back_wall")
    for element, fraction in [
        ("Fe", 0.885), ("Cr", 0.09), ("W", 0.011), ("Mn", 0.006),
        ("V", 0.002), ("C", 0.006),
    ]:
        back_wall.add_element(element, fraction, "wo")
    back_wall.set_density("g/cm3", 7.9)

    shield = openmc.Material(name="shield")
    shield.add_element("W", 1.0)
    shield.add_element("C", 1.0)
    shield.set_density("g/cm3", 15.6)

    gap = openmc.Material(name="gap")
    gap.add_nuclide("H1", 1.0)
    gap.set_density("g/cm3", 1.0e-8)

    vessel = openmc.Material(name="vac_vessel")
    for element, fraction in [
        ("Fe", 0.62), ("Cr", 0.16), ("Ni", 0.11), ("Mo", 0.02),
        ("Mn", 0.015), ("B", 0.03), ("O", 0.03), ("H", 0.005),
    ]:
        vessel.add_element(element, fraction, "wo")
    vessel.set_density("g/cm3", 7.9)

    thermal_shield = openmc.Material(name="thermal_shield")
    thermal_shield.add_element("Fe", 0.89, "wo")
    thermal_shield.add_element("Cr", 0.09, "wo")
    thermal_shield.add_element("Mn", 0.02, "wo")
    thermal_shield.set_density("g/cm3", 7.9)

    magnets = openmc.Material(name="magnets")
    for element, fraction in [
        ("Cu", 0.58), ("Fe", 0.30), ("Cr", 0.04), ("Ni", 0.03),
        ("Y", 0.02), ("Ba", 0.02), ("O", 0.01),
    ]:
        magnets.add_element(element, fraction, "wo")
    magnets.set_density("g/cm3", 8.5)

    vacuum = openmc.Material(name="Vacuum")
    vacuum.add_nuclide("H1", 1.0)
    vacuum.set_density("g/cm3", 1.0e-10)

    return openmc.Materials([
        first_wall, multiplier, breeder, back_wall, shield, gap, vessel,
        thermal_shield, magnets, vacuum,
    ])
