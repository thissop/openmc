#!/usr/bin/env python
"""Generate the MULTI-COIL placement delta on the corrected (13,19) build grid -- the fix for the
kill-shot peak migration (single-coil coil-20 placement over-protected coil 20; the peak jumped to
the unprotected coil 25).

Same recipe as gen_step1_delta_corr.py (contributon -> priority -> delta), but:
  * the importance is the MULTI-COIL map (sum/max/flux-weighted over all 20 coils), not coil-20;
  * the contributon is FOLDED into one field period [0, 2pi/nfp) before splatting, because the
    multi-coil importance has mass in every period (coil-20's did not) and the per-period delta is
    stellarator-symmetric.

Two deltas, both at the SAME material budget as the single-coil 'placed' field (630.66 cm) so the
peak-coil-flux comparison isolates WHICH coils the placement protects, not how much material:
  * delta_multi_prop   = DELTA_MAX * P_multi              (proportional -- matches placed's shape)
  * delta_multi_greedy = greedy_allocate(P_multi, budget) (bang-bang -- validates placement_allocator)
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))     # shield_opt
from adjoint_placement import contributon, importance_angles  # noqa: E402
import placement_allocator as pa  # noqa: E402
from thickness_field import ThicknessField  # noqa: E402

TOR = list(np.linspace(0.0, 90.0, 13))
POL = list(np.linspace(0.0, 360.0, 19))
T_BR0, T_SH0, DELTA_MAX, T_BR_MIN = 50.0, 40.0, 35.0, 15.0
PLACED_BUDGET = 630.66                          # sum(delta_placed) of the single-coil field


def splat_folded(field, tor_deg, pol_deg, nfp, width_deg=15.0, R0=None):
    """Priority on the (tor,pol) grid from a contributon field, folding phi into one field period.
    Mirrors adjoint_placement.placement_priority but folds phi -> [0, 2pi/nfp)."""
    phi, theta, w, R0 = importance_angles(field, R0=R0)
    phi = np.mod(phi, 2 * np.pi / nfp)                         # FOLD into one period
    tor = np.radians(np.asarray(tor_deg)); pol = np.radians(np.asarray(pol_deg))
    TORg, POLg = np.meshgrid(tor, pol, indexing="ij")
    sig = np.radians(width_deg)

    def dper(a, b):
        return (a - b + np.pi) % (2 * np.pi) - np.pi

    P = np.zeros(TORg.shape)
    for pj, tj, wj in zip(phi, theta, w):
        P += wj * np.exp(-(dper(TORg, pj) ** 2 + dper(POLg, tj) ** 2) / (2 * sig ** 2))
    return P / (P.max() + 1e-300), float(R0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, help="multi-coil importance npz (importance = psi_multi)")
    ap.add_argument("--fluxmap", required=True)
    ap.add_argument("--nfp", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    field = contributon(args.map, args.fluxmap, scale=100.0, emissivity="uniform")
    P, R0 = splat_folded(field, TOR, POL, args.nfp)
    assert P.shape == (13, 19)

    tf = ThicknessField(nfp=args.nfp, toroidal_angles_deg=TOR, poloidal_angles_deg=POL,
                        t_breeder0=T_BR0, t_shield0=T_SH0, t_breeder_min=T_BR_MIN)

    # proportional, scaled to the placed budget (then clipped to delta_max)
    c = PLACED_BUDGET / float(P.sum())
    delta_prop = np.clip(c * P, 0.0, DELTA_MAX)
    # greedy at the same budget (fraction of max addable = budget / (delta_max * ncells))
    bf = PLACED_BUDGET / (DELTA_MAX * P.size)
    delta_greedy = pa.greedy_allocate(P, tf, budget_frac=bf)

    for name, d in (("prop", delta_prop), ("greedy", delta_greedy)):
        br = T_BR0 - d
        assert br.min() >= T_BR_MIN - 1e-9, f"{name}: breeder floor {br.min():.2f}"
        print(f"delta_multi_{name}: sum={d.sum():.2f} cm  max={d.max():.2f}  mean={d.mean():.3f}  "
              f"breeder_min={br.min():.2f}  n_active={(d>1e-6).sum()}")

    np.savez(args.out, delta_multicoil=delta_prop, delta_multi_greedy=delta_greedy,
             delta_uniform=np.full_like(delta_prop, float(delta_prop.mean())),
             priority=P, R0=R0, toroidal_angles=np.array(TOR), poloidal_angles=np.array(POL),
             t_breeder0=T_BR0, t_shield0=T_SH0, delta_max=DELTA_MAX, t_breeder_min=T_BR_MIN,
             importance="multicoil", build="corrected_w35f")
    print(f"WROTE {args.out}")


if __name__ == "__main__":
    main()
