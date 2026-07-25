#!/usr/bin/env python
"""Closed loop: adjoint magnet-importance -> Bayesian shield optimizer.

The optimizer trades breeder for shield at each (theta, phi) to cut peak coil dose
under a fixed envelope + TBR floor. Its key physical input is the baseline coil-dose
field per (theta, phi) sightline -- i.e. WHERE the coil load originates. That is
exactly what the adjoint importance map provides (through the shield): the
contributon projected onto the (theta, phi) control grid.

This script closes the loop and demonstrates the payoff:
  * ADJOINT-INFORMED: baseline coil-dose field = the adjoint placement priority
    (peaked on the coil-facing sightlines). The optimizer localizes shield there.
  * UNIFORM (no attribution): baseline field = flat (what you get WITHOUT the
    adjoint map -- a uniform guess). The optimizer has no spatial signal.
Both run under the identical envelope, TBR floor, and regularizers. With the adjoint
attribution the optimizer kills the coil-dose PEAK for a small breeder sacrifice;
the uniform baseline cannot target and must spread shield thin -> worse peak.

Run:  python adjoint_closed_loop.py [--map ...] [--fluxmap ...] [--coil-centroid x y z]
"""
from __future__ import annotations

import argparse
import os

import numpy as np

import adjoint_placement as apl
import bayes_optimizer as bo
from thickness_field import ThicknessField

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_MAP = os.path.join(HERE, "..", "data", "adjoint", "coil15_fil_flat.npz")
DEF_FLUX = os.path.join(HERE, "..", "data", "qh_freestream_fluxmap.npz")


def build_problem(baseline_field, tf, basis, tbr0=1.15, tbr_floor=1.05):
    # Pure shield-PLACEMENT demo: no standoff lever, and light manufacturability/material
    # regularizers so the objective is dominated by the engineering goal (peak coil dose)
    # and the TBR floor -- not by cm^3-scaled volume terms. (The default weights are tuned
    # for real OpenMC dose magnitudes; here the dose field is normalized to peak 1.)
    return bo.OptProblem(baseline_coil_dose_field=baseline_field, tf=tf, basis=basis,
                         tbr0=tbr0, tbr_floor=tbr_floor, n_standoff=0,
                         lam_smooth=0.03, lam_material=1.0e-4, lam_standoff=0.0,
                         penalty=50.0)


def run_case(name, baseline_field, tf, basis, seed=0):
    p = build_problem(baseline_field, tf, basis)
    x_uniform = p.uniform_field_design(bias=-6.0)     # ~no shield: the starting build
    J0 = bo.objective(x_uniform, None, None, None, None, problem=p)
    peak0 = bo.peak_dose(x_uniform, p)
    res = bo.optimize(p, n_bootstrap=40, n_verify=10, n_init=10, seed=seed, verbose=False)
    peak1 = res["peak_dose"]
    sac = bo.breeder_sacrificed_fraction(res["best_design"], p)
    print(f"[{name:16s}] peak dose {peak0:.4f} -> {peak1:.4f}  "
          f"({100*(1-peak1/peak0):+.1f}%)   breeder sacrificed {100*sac:.1f}%   "
          f"TBR {res['tbr']:.3f}   J {J0:.3f}->{res['best_value']:.3f}")
    return dict(name=name, peak0=peak0, peak1=peak1, reduction=1 - peak1 / peak0,
                sacrificed=sac, tbr=res["tbr"], design=res["best_design"], problem=p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", default=DEF_MAP)
    ap.add_argument("--fluxmap", default=DEF_FLUX)
    ap.add_argument("--coil-centroid", type=float, nargs=3, default=[-891.95, -743.44, -74.47])
    ap.add_argument("--nfp", type=int, default=4)
    ap.add_argument("--n-tor", type=int, default=36)
    ap.add_argument("--n-pol", type=int, default=36)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--emissivity", default="bosch_hale", choices=["bosch_hale", "uniform"],
                    help="plasma reactivity profile for S(r) in the contributon. "
                         "bosch_hale = realistic core-peaked (PHYSICAL, the default); "
                         "uniform over-weights the shaped edge and OVER-STATES the benefit.")
    ap.add_argument("--fig", default=None, help="save the closed-loop figure to this path")
    args = ap.parse_args()

    tor = np.linspace(0, 360, args.n_tor, endpoint=False)
    pol = np.linspace(0, 360, args.n_pol, endpoint=False)

    # --- adjoint attribution -> placement priority on the optimizer's control grid ---
    # Use the PHYSICAL (Bosch-Hale core-peaked) reactivity by default: uniform emissivity
    # over-weights the strongly-shaped plasma edge and overstates the localizable benefit
    # (~55% vs ~19% realistic on QH coil-15). The adjoint psi_dagger is reactivity-
    # independent; only this S(r) weighting changes.
    print(f"reactivity profile (S(r) weighting): {args.emissivity}")
    field = apl.contributon(args.map, args.fluxmap, scale=100.0, emissivity=args.emissivity)
    res = apl.placement_priority(field, tor, pol)
    P = res["priority"]                                # (n_tor, n_pol) in [0,1]
    cphi, _ = apl.coil_angles(args.coil_centroid, res["R0"])
    dphi = ((res["phi0_deg"] - cphi + 180) % 360) - 180
    print(f"adjoint placement target phi0={res['phi0_deg']:.0f} deg "
          f"(coil phi={cphi:.0f}, |dphi|={abs(dphi):.0f}); R0={res['R0']:.0f} cm\n")

    # baseline coil-dose field: adjoint-informed (the priority) vs uniform (flat guess).
    # Scale both to the same PEAK so the comparison is about the SHAPE (localization),
    # not the overall dose level.
    adj_field = P / P.max()
    uni_field = np.full_like(P, 1.0)                   # no spatial information

    tf = ThicknessField(nfp=args.nfp, toroidal_angles_deg=tor, poloidal_angles_deg=pol,
                        t_breeder0=80.0, t_shield0=20.0, t_breeder_min=10.0)
    basis = tf.fourier_basis(M=3, N=2)

    print("closed-loop shield optimization (surrogate objective, fixed envelope + TBR floor):")
    a = run_case("adjoint-informed", adj_field, tf, basis, seed=args.seed)
    u = run_case("uniform (naive)", uni_field, tf, basis, seed=args.seed)

    print(f"\nHEADLINE: adjoint-informed cuts peak coil dose {100*a['reduction']:.1f}% "
          f"vs {100*u['reduction']:.1f}% for the uniform guess "
          f"(breeder sacrificed {100*a['sacrificed']:.1f}% vs {100*u['sacrificed']:.1f}%).")

    if args.fig:
        _figure(a, tf, basis, P, tor, pol, args.fig)
        print("  wrote", args.fig)
    return a, u


def _figure(a, tf, basis, P, tor, pol, out):
    import smplotlib  # noqa: F401  house style before pyplot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, coeffs = a["problem"].split(a["design"])
    delta = tf.delta(coeffs, basis)                    # optimized shield thickness (cm)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    im0 = ax[0].pcolormesh(tor, pol, P.T, cmap="magma", shading="auto")
    ax[0].set_title("Adjoint Placement Priority (Where Coil Load Originates)")
    ax[0].set_xlabel("Toroidal Angle (deg)"); ax[0].set_ylabel("Poloidal Angle (deg)")
    fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04, label="Priority")
    im1 = ax[1].pcolormesh(tor, pol, delta.T, cmap="viridis", shading="auto")
    ax[1].set_title("Optimized Shield Added $\\delta(\\theta,\\phi)$ (cm)")
    ax[1].set_xlabel("Toroidal Angle (deg)"); ax[1].set_ylabel("Poloidal Angle (deg)")
    fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04, label="Shield Added (cm)")
    fig.suptitle(f"Closed Loop: Adjoint-Informed Shield Cuts Peak Coil Dose "
                 f"{100*a['reduction']:.0f}% (Breeder $-${100*a['sacrificed']:.0f}%, "
                 f"TBR {a['tbr']:.2f})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
