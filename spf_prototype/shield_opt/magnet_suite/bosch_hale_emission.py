#!/usr/bin/env python
"""Real Bosch-Hale DT emissivity for the TRANSPORT source, consistent with the adjoint.

Closes emissivity gap item (3): the native ``StellaratorSource`` births neutrons with
per-cell weight ``w = S(rho) * |sqrt(g)|`` (see native_spf/stellarator_source.cpp
``build_cdf_cascade``). ``S(rho)`` is the RADIAL emission profile: it defaults to the
crude parabolic ``1 - rho^2`` and is OVERRIDDEN by the ``emission=(rho_grid, density)``
API of ``openmc.StellaratorSource``. Setting

    S(rho) = r_BH(rho) = n(rho)^2 * <sigma v>_DT(T(rho))      (shield_opt/emissivity.py)

makes the transport birth density equal the REAL DT reaction-rate density
r_BH(rho)*|sqrt(g)| -- the same reactivity the adjoint contributon already folds in via
``adjoint_placement.plasma_source_on_mesh(..., emissivity='bosch_hale')`` (which weights
|sqrt(g)| by the identical r(rho)). So BOTH attribution and transport use ONE emissivity.

Primary API (do this in the transport driver):

    from bosch_hale_emission import bosch_hale_emission
    rho, dens = bosch_hale_emission(fluxmap_npz_or_rho, kind='bosch_hale')
    src = openmc.StellaratorSource(fluxmap=stem, emission=(rho, dens), ...)

We deliberately do NOT bake r(rho) into the fluxmap's ``sqrtg`` array (that array also
carries the geometric volume used by the volume rigor gates); the ``emission`` profile is
the designed, side-effect-free knob. A ``bake_weighted_fluxmap`` helper is provided only
for producers that read ``sqrtg`` directly and cannot pass an emission profile.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

# reuse the single source of truth for r(rho)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from emissivity import emissivity_weight  # noqa: E402


def _rho_of(fluxmap):
    """Accept a fluxmap npz path/handle or a bare rho array; return the 1-D rho grid."""
    if isinstance(fluxmap, (str, Path)):
        fluxmap = np.load(fluxmap)
    if hasattr(fluxmap, "files"):          # npz
        return np.asarray(fluxmap["rho"], float)
    return np.asarray(fluxmap, float)


def bosch_hale_emission(fluxmap, kind="bosch_hale"):
    """Return ``(rho_grid, density)`` to pass as ``StellaratorSource(emission=...)``.

    density = r(rho) normalized to max 1 (only the shape matters; the C++ re-normalizes
    the CDF).  ``kind='uniform'`` returns ones (the geometry-only S_v that matches a
    uniform-source forward/reciprocity run).
    """
    rho = _rho_of(fluxmap)
    dens = emissivity_weight(rho, kind=kind)
    return rho, np.asarray(dens, float)


def birth_pdf_rho(fluxmap, kind="bosch_hale"):
    """Marginal birth pdf p(rho) actually produced by the C++ cascade for emissivity ``kind``.

    p(rho_i) proportional to  S(rho_i) * sum_{theta,zeta} |sqrt(g)|_{i,.,.}  (build_cdf_cascade),
    normalized to sum 1.  Used for verification / core-peaking diagnostics.
    """
    fm = np.load(fluxmap) if isinstance(fluxmap, (str, Path)) else fluxmap
    rho = np.asarray(fm["rho"], float)
    sqrtg = np.abs(np.asarray(fm["sqrtg"], float))       # (nr, nt, nz)
    S = emissivity_weight(rho, kind=kind)
    m = S * sqrtg.sum(axis=(1, 2))                        # unnormalized marginal mass
    return rho, m / m.sum()


def bake_weighted_fluxmap(src_npz, dst_stem, kind="bosch_hale"):
    """Write a standalone fluxmap whose ``sqrtg`` is multiplied by r(rho) (for producers
    that cannot pass an emission profile). Emits ``<dst_stem>.npz`` + ``<dst_stem>.{meta,bin}``.
    Prefer the emission API instead; use this only when a tool reads sqrtg directly.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "python"))
    from quasr_fluxmap import write_fluxmap_bin  # noqa: E402
    d = np.load(src_npz)
    rho = np.asarray(d["rho"], float)
    r = emissivity_weight(rho, kind=kind)
    sg = np.abs(np.asarray(d["sqrtg"], float)) * r[:, None, None]
    out = {k: d[k] for k in d.files}
    out["sqrtg"] = sg
    out["emissivity_kind"] = kind
    np.savez(str(dst_stem) + ".npz", **out)
    write_fluxmap_bin(str(dst_stem), rho, d["theta"], d["zeta"], sg,
                      d["R"], d["Z"], d["phi"], d["BR"], d["Bphi"], d["BZ"],
                      int(d["nfp"]), float(d.get("sign_sqrtg", 1.0)))
    return str(dst_stem) + ".npz"


def _verify(fluxmap):
    """Print a core-peaking table proving BH shifts births inward vs uniform / parabolic."""
    print(f"# Bosch-Hale emissivity verification on {os.path.basename(str(fluxmap))}")
    print(f"{'kind':>12} | {'<rho>':>7} | {'frac rho<0.5':>12} | {'p(rho) peak at':>14}")
    print("-" * 56)
    ref = {}
    for kind in ("uniform", "one_minus_rho2", "bosch_hale"):
        if kind == "one_minus_rho2":
            fm = np.load(fluxmap) if isinstance(fluxmap, str) else fluxmap
            rho = np.asarray(fm["rho"], float)
            S = 1.0 - rho ** 2
            m = S * np.abs(fm["sqrtg"]).sum(axis=(1, 2))
            p = m / m.sum()
        else:
            rho, p = birth_pdf_rho(fluxmap, kind=kind)
        mean_rho = float((rho * p).sum())
        frac_core = float(p[rho < 0.5].sum())
        peak = float(rho[np.argmax(p)])
        ref[kind] = mean_rho
        print(f"{kind:>12} | {mean_rho:7.3f} | {frac_core:12.3f} | {peak:14.3f}")
    print("-" * 56)
    shift = ref["uniform"] - ref["bosch_hale"]
    print(f"BH pulls <rho> inward by {shift:+.3f} vs uniform "
          f"({'CORE-PEAKED OK' if shift > 0 else 'FAIL: not core-peaked'})")
    assert shift > 0, "Bosch-Hale must move births toward the core (rho=0)"
    print("PASS: Bosch-Hale emissivity is core-peaked relative to uniform.")


if __name__ == "__main__":
    fx = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "data", "quasr112734_vmec_fluxmap.npz")
    _verify(fx)
