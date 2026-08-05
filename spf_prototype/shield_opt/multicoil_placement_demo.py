#!/usr/bin/env python
"""End-to-end multi-coil placement: coil importance maps -> combine -> greedy allocate ->
ParaStell radial-build dict, with the peak-migration failure made explicit.

Pipeline (all offline; real per-coil adjoint maps drop in where the synthetic ones are built):

    per-coil priority P_c(theta,phi)   [adjoint_placement.placement_priority on each coil map]
        -> combine_priorities(...)      [sum / max / flux-weighted -> one field]
        -> greedy_allocate(P, tf, budget)   [fixed material where it buys the most protection]
        -> tf.radial_build_dict(delta, base)  [DAGMC-ready breeder<->shield trade]

The point: SINGLE-coil placement (allocate on one coil's priority) over-protects that coil and
leaves the others exposed -- the peak coil dose MIGRATES to an unprotected coil (the kill-shot
9193680 result: crushing coil 20 pushed the peak to coil 25). MULTI-coil placement (combined
priority, same material budget) protects the whole set and gives a lower peak. This script
quantifies that with the surrogate coil-dose model and emits the radial-build dict.
"""
from __future__ import annotations

import argparse

import numpy as np

import attenuation_surrogate as att
import placement_allocator as pa
from thickness_field import ThicknessField


# ----------------------------------------------------------------------------------------------
# Synthetic per-coil fields (stand-ins for real adjoint placement priorities). Each coil c sits at
# a distinct toroidal angle phi_c; its shield-load originates near phi_c (a poloidal inboard band
# x a toroidal lobe at phi_c), and its dose baseline shares that footprint. Real runs replace
# `coil_fields` with adjoint_placement.placement_priority(map_c) per coil.
# ----------------------------------------------------------------------------------------------
def coil_fields(tf, n_coils=5, hot=1.6, tor_width_deg=22.0, hotness=None):
    """Return (priorities, baselines, phi_c_deg): each (K, nTor, nPol). priorities are placement
    priorities per coil; baselines are per-coil surrogate dose footprints (same shape). `hotness`
    (len K, default all 1) scales each coil's baseline dose -> a non-uniform peak set."""
    TOR, POL = tf.TOR, tf.POL                                   # radians, (nTor, nPol)
    phi_c = np.linspace(0.0, 2 * np.pi, n_coils, endpoint=False)
    sig = np.radians(tor_width_deg)
    hotness = np.ones(n_coils) if hotness is None else np.asarray(hotness, float)

    def dperiodic(a, b):
        return (a - b + np.pi) % (2 * np.pi) - np.pi

    inboard = np.exp(hot * (np.cos(POL - np.pi) - 1.0))         # peaks at theta=pi (inboard)
    P, B = [], []
    for c, pc in enumerate(phi_c):
        lobe = np.exp(-0.5 * (dperiodic(TOR, pc) / sig) ** 2)   # toroidal lobe at coil c
        f = inboard * lobe
        P.append(f / (f.max() + 1e-300))
        B.append(hotness[c] * f)
    return np.array(P), np.array(B), np.degrees(phi_c)


def coil_doses(baselines, delta, k):
    """Per-coil integrated surrogate dose = sum_(theta,phi) baseline_c * exp(-k*delta)."""
    atten = np.exp(-k * np.asarray(delta))
    return np.array([float(np.sum(b * atten)) for b in baselines])


def peak_coil_dose(baselines, delta, k):
    return float(coil_doses(baselines, delta, k).max())


def run_pipeline(priorities, baselines, tf, budget_frac=0.3, combine_mode="sum",
                 weights=None, base_dict=None):
    """priorities -> combined field -> greedy delta -> radial-build dict + diagnostics."""
    P = pa.combine_priorities(priorities, mode=combine_mode, weights=weights)
    delta = pa.greedy_allocate(P, tf, budget_frac=budget_frac)
    k = tf.k if hasattr(tf, "k") else att.trade_k()
    out = dict(
        priority=P, delta=delta,
        peak_dose=peak_coil_dose(baselines, delta, k),
        coil_doses=coil_doses(baselines, delta, k),
        coverage=pa.coverage(P, delta, tf),
        added_volume=tf.added_shield_volume(delta),
    )
    if base_dict is not None:
        out["radial_build"] = tf.radial_build_dict(delta, base_dict)
    return out


def base_radial_build(tf):
    """A minimal ParaStell-style base dict the trade slots breeder/shield into."""
    ones = np.ones(tf.shape)
    return {
        "first_wall": {"thickness_matrix": ones * 3.2},
        "multiplier": {"thickness_matrix": ones * 2.0},
        "breeder":    {"thickness_matrix": ones * tf.t_b0},
        "back_wall":  {"thickness_matrix": ones * 4.0},
        "shield":     {"thickness_matrix": ones * tf.t_s0},
        "vacuum_vessel": {"thickness_matrix": ones * 25.0},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-coils", type=int, default=5)
    ap.add_argument("--nfp", type=int, default=4)
    ap.add_argument("--n-tor", type=int, default=48)
    ap.add_argument("--n-pol", type=int, default=36)
    ap.add_argument("--budget-frac", type=float, default=0.3,
                    help="material budget as a fraction of the max addable shield volume")
    ap.add_argument("--combine", default="sum", choices=["sum", "max", "weighted"])
    ap.add_argument("--fig", default=None)
    args = ap.parse_args()

    tor = np.linspace(0, 360, args.n_tor, endpoint=False)
    pol = np.linspace(0, 360, args.n_pol, endpoint=False)
    tf = ThicknessField(nfp=args.nfp, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                        t_breeder0=80.0, t_shield0=20.0, t_breeder_min=15.0)
    tf.k = att.trade_k()

    # make one coil genuinely hottest so a naive single-coil placement has an obvious target
    hotness = np.ones(args.n_coils); hotness[0] = 1.35
    P, B, phi_c = coil_fields(tf, n_coils=args.n_coils, hotness=hotness)
    base = base_radial_build(tf)

    no_shield_peak = peak_coil_dose(B, np.zeros(tf.shape), tf.k)
    print(f"coils at phi = {np.round(phi_c,0)} deg;  hottest = coil 0 (x{hotness[0]:.2f})")
    print(f"baseline peak coil dose (no shield) = {no_shield_peak:.2f}\n")

    # (1) SINGLE-coil placement: greedy on the hottest coil's priority alone.
    single = run_pipeline(P[[0]], B, tf, budget_frac=args.budget_frac, combine_mode="sum",
                          base_dict=base)
    # (2) MULTI-coil placement: same budget, combined priority over all coils.
    multi = run_pipeline(P, B, tf, budget_frac=args.budget_frac, combine_mode=args.combine,
                         base_dict=base)

    def report(name, r):
        cd = r["coil_doses"]
        print(f"[{name:22s}] peak coil dose {r['peak_dose']:.2f}  "
              f"({100*(1-r['peak_dose']/no_shield_peak):+.1f}% vs no-shield)  "
              f"hottest=coil {int(np.argmax(cd))}  coverage={r['coverage']:.2f}  "
              f"added_vol={r['added_volume']:.1f} cm")

    print("fixed material budget = "
          f"{args.budget_frac:.2f} x max addable shield  (both cases identical):")
    report("single-coil (naive)", single)
    report("multi-coil (combined)", multi)

    imp = 100 * (1 - multi["peak_dose"] / single["peak_dose"])
    migrated = int(np.argmax(single["coil_doses"])) != 0
    print(f"\nHEADLINE: at equal material, multi-coil placement lowers the peak coil dose "
          f"{imp:+.1f}% vs single-coil.")
    print(f"  single-coil peak {'MIGRATED to coil %d (over-protected coil 0)' % np.argmax(single['coil_doses']) if migrated else 'stayed on coil 0'}; "
          f"multi-coil peak on coil {int(np.argmax(multi['coil_doses']))}.")
    print(f"  radial-build dict emitted (breeder/shield traded, envelope conserved): "
          f"keys={list(multi['radial_build'].keys())}")

    if args.fig:
        _figure(tf, tor, pol, single, multi, phi_c, args.fig)
        print("  wrote", args.fig)
    return single, multi


def _figure(tf, tor, pol, single, multi, phi_c, out):
    import smplotlib  # noqa: F401
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, r, ttl in ((ax[0], single, "Single-Coil Placement (Peak Migrates)"),
                      (ax[1], multi, "Multi-Coil Placement (Whole Set Protected)")):
        im = a.pcolormesh(tor, pol, r["delta"].T, cmap="gray_r", shading="auto")
        a.set_xlabel("Toroidal Angle (deg)"); a.set_ylabel("Poloidal Angle (deg)")
        a.set_title(ttl)
        for pc in phi_c:
            a.axvline(pc, color="#b2182b", lw=0.8, ls="--")
        fig.colorbar(im, ax=a, fraction=0.046, pad=0.04, label="Shield Added (cm)")
    fig.tight_layout(); fig.savefig(out, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
