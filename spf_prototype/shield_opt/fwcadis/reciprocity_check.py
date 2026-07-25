#!/usr/bin/env python
"""Forward<->adjoint reciprocity validation of the magnet-importance maps (#2).

The adjoint contributon integral for a coil, INT S(r) psi_dagger(r) dr, is the coil's
neutron load from the full plasma source -- the SAME physical quantity an independent
FORWARD transport calculation gives (the per-coil coil flux, from emitting the plasma
source and tallying each coil). Reciprocity therefore predicts:

    INT S * psi_dagger  (adjoint, per coil)   ~   forward per-coil coil flux

across coils. This script computes the adjoint prediction from the committed per-coil
adjoint maps (the kerma sweep) + the plasma source S(r), and correlates it against the
committed forward per-coil flux (percoil_*.npz). A strong rank correlation validates
the adjoint attribution's MAGNITUDE/ranking (the direction was validated separately by
adjoint_placement). This uses only committed data -- no new transport run.

Run (on Ginsburg, after the kerma sweep):
  python reciprocity_check.py --sweep-dir ~/adjoint_test --coils 15 25 30 20 \
      --fluxmap .../qh_freestream_fluxmap.npz --percoil .../percoil_unpol_corr_w35f.npz
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # shield_opt/
from adjoint_placement import plasma_source_on_mesh   # noqa: E402


def adjoint_coil_load(map_npz, fluxmap, scale=100.0, emissivity="uniform"):
    """INT S(r) psi_dagger(r) dr for one coil's adjoint map (the contributon integral).

    emissivity='uniform' to MATCH a uniform-source forward run (reciprocity validation);
    'bosch_hale' for the PHYSICAL per-coil load ranking (free -- psi_dagger is
    reactivity-independent, so this needs no new transport)."""
    d = np.load(map_npz)
    imp = d["importance"].astype(float)
    S = plasma_source_on_mesh(fluxmap, d["lower_left"], d["upper_right"], d["dimension"],
                              scale=scale, emissivity=emissivity)
    return float((S * imp).sum())


def _parse_coils(path, scale):
    """MAKEGRID coils file -> list of (Npts,3) filament polylines (scaled), same as
    adjoint_importance._parse_coils. Used to length-normalize the reciprocity integral."""
    coils, cur = [], []
    for line in open(path):
        p = line.split()
        if len(p) >= 4:
            try:
                x, y, z, c = map(float, p[:4])
            except ValueError:
                continue
            cur.append((x, y, z))
            if c == 0.0:
                coils.append(np.array(cur) * scale)
                cur = []
    return coils


def coil_filament_lengths(coils_file, coil_centroids, cells, scale=100.0):
    """Arc length (cm) of each requested coil's guide curve, matched by nearest centroid
    (the SAME matching adjoint_importance uses). Longer windings spread the fixed 40-point
    source over more plasma, so raw INT(S*psi) over-weights them vs the forward flux
    DENSITY -- dividing by length removes that geometric bias."""
    fils = _parse_coils(coils_file, scale)
    fcen = np.array([f.mean(0) for f in fils])
    cc = np.load(coil_centroids)
    cids = [int(x) for x in cc["cell_ids"]]
    dcen = cc["centroids"]
    out = {}
    for cid in cells:
        j = int(np.linalg.norm(fcen - dcen[cids.index(cid)], axis=1).argmin())
        f = fils[j]
        out[cid] = float(np.sum(np.linalg.norm(np.diff(f, axis=0), axis=1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-dir", required=True,
                    help="dir holding sweep_coil<ID>_kerma/adjoint_importance_kerma.npz")
    ap.add_argument("--coils", type=int, nargs="+", required=True)
    ap.add_argument("--fluxmap", required=True)
    ap.add_argument("--fluxmap-scale", type=float, default=100.0)
    ap.add_argument("--percoil", required=True, help="forward per-coil flux npz (cells, flux)")
    ap.add_argument("--response", default="flat",
                    help="response tag in the map filename (match the forward tally: 'flat' "
                         "= total/fast flux, 'kerma' = heating)")
    ap.add_argument("--map-pattern",
                    default="recip_w35f_coil{cid}_{response}/adjoint_importance_{response}_P0.npz",
                    help="per-coil map path under --sweep-dir; {cid},{response} substituted")
    ap.add_argument("--coils-file", default="/burg-archive/home/tjk2147/pstl_test/coils_qh")
    ap.add_argument("--coil-centroids",
                    default="/burg-archive/home/tjk2147/pstl_test/corrected/coil_centroids_corr.npz")
    ap.add_argument("--emissivity", default="uniform", choices=["uniform", "bosch_hale"],
                    help="uniform to MATCH the uniform-source forward run (validation); "
                         "bosch_hale for the physical per-coil ranking (free, reactivity-"
                         "independent adjoint).")
    args = ap.parse_args()

    pc = np.load(args.percoil)
    fcells = [int(c) for c in pc["cells"]]
    fflux = np.asarray(pc["flux"], float)

    rows = []
    for cid in args.coils:
        mp = os.path.join(args.sweep_dir,
                          args.map_pattern.format(cid=cid, response=args.response))
        if not os.path.exists(mp):
            print(f"  coil {cid}: MISSING {mp} (sweep not done?)")
            continue
        adj = adjoint_coil_load(mp, args.fluxmap, scale=args.fluxmap_scale,
                                emissivity=args.emissivity)
        rows.append((cid, adj, fflux[fcells.index(cid)]))

    if len(rows) < 3:
        print("need >=3 coils for a correlation; sweep incomplete.")
        return

    cells = [r[0] for r in rows]
    lengths = coil_filament_lengths(args.coils_file, args.coil_centroids, cells,
                                    scale=args.fluxmap_scale)
    adj = np.array([r[1] for r in rows])
    fwd = np.array([r[2] for r in rows])
    length = np.array([lengths[c] for c in cells])
    adj_density = adj / length                          # length-normalized (flux-density analog)

    print(f"  {'coil':>4}  {'INT(S*psi)':>12}  {'/length':>12}  {'fwd fast flux':>13}  {'len(m)':>7}")
    for c, a, ad, f, L in zip(cells, adj, adj_density, fwd, length):
        print(f"  {c:>4}  {a:>12.3e}  {ad:>12.3e}  {f:>13.3e}  {L/100:>7.1f}")

    def corr(x, y):
        pear = float(np.corrcoef(x, y)[0, 1])
        rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
        spear = float(np.corrcoef(rx, ry)[0, 1])
        return pear, spear

    pr, sr = corr(adj, fwd)
    pd, sd = corr(adj_density, fwd)
    print(f"\nRECIPROCITY across {len(rows)} coils (matched geometry + response):")
    print(f"  raw    INT(S*psi)  vs forward flux : Pearson {pr:+.3f}  Spearman {sr:+.3f}")
    print(f"  length-normalized  vs forward flux : Pearson {pd:+.3f}  Spearman {sd:+.3f}")
    best = max(sr, sd)
    print("  -> adjoint attribution reproduces the forward per-coil loading (reciprocity holds)"
          if best > 0.6 else
          "  -> correlation still weak: investigate further")


if __name__ == "__main__":
    main()
