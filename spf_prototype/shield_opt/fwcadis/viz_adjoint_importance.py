#!/usr/bin/env python
"""Visualize an adjoint magnet-importance map and check it points where physics says.

Input: an adjoint_importance_*.npz written by adjoint_importance.py (fields:
importance[nx,ny,nz], lower_left, upper_right, dimension, response, coil_cells).

The adjoint flux psi_dagger(r) is the importance of a source neutron born at r to
the chosen coil's response (heating), THROUGH the shield+blanket. Physically it must
be biased toward the plasma region facing the target coil and fall off with
shielding/distance. We (a) project it (top + side views, log color) with the coil
marked, and (b) quantify the bias: the importance-weighted centroid should sit on the
coil side of the device, i.e. its unit vector should have positive projection on the
coil-centroid direction.

House style: smplotlib imported before any pyplot use; Title Case labels.

Usage:
  python viz_adjoint_importance.py <map.npz> [--coil-centroid X Y Z] [--out fig.png]
"""
import argparse
from pathlib import Path

import numpy as np
import smplotlib  # noqa: F401  house serif style — before pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402


def voxel_centers(ll, ur, dim):
    return [np.linspace(ll[i], ur[i], dim[i] + 1)[:-1] + 0.5 * (ur[i] - ll[i]) / dim[i]
            for i in range(3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--fluxmap", default=None,
                    help="VMEC fluxmap npz -> plot the CONTRIBUTON C=S*psi_dagger "
                         "(physical attribution) instead of the bare adjoint flux.")
    ap.add_argument("--fluxmap-scale", type=float, default=100.0)
    ap.add_argument("--coil-centroid", type=float, nargs=3, default=None,
                    help="Coil centroid (cm) to overlay + test the bias against.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.fluxmap is not None:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
        from adjoint_placement import contributon
        d = contributon(args.npz, args.fluxmap, scale=args.fluxmap_scale)
    else:
        d = np.load(args.npz)
    imp = d["importance"].astype(float)
    ll, ur, dim = d["lower_left"], d["upper_right"], d["dimension"]
    resp = str(d["response"]) if "response" in d else "?"
    coils = d["coil_cells"].tolist() if "coil_cells" in d else "?"
    xc, yc, zc = voxel_centers(ll, ur, dim)

    # --- quantitative "points toward the coil" test -------------------------
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    w = imp
    tot = w.sum()
    cen = np.array([(X * w).sum(), (Y * w).sum(), (Z * w).sum()]) / tot   # importance centroid
    msg = [f"map: {Path(args.npz).name}",
           f"response={resp}  coil_cells={coils}",
           f"nonzero={(imp>0).sum()}/{imp.size} ({(imp>0).mean():.0%})  peak/mean(nz)={imp.max()/imp[imp>0].mean():.1f}",
           f"importance-weighted centroid (cm) = [{cen[0]:.0f}, {cen[1]:.0f}, {cen[2]:.0f}]"]
    bias = None
    if args.coil_centroid is not None:
        coil = np.array(args.coil_centroid)
        u_coil = coil / np.linalg.norm(coil)
        # project centroid offset (from device center) onto coil direction
        proj = float(cen @ u_coil)
        cos = float((cen @ coil) / (np.linalg.norm(cen) * np.linalg.norm(coil) + 1e-30))
        bias = (proj, cos)
        msg.append(f"coil centroid (cm) = [{coil[0]:.0f}, {coil[1]:.0f}, {coil[2]:.0f}]")
        msg.append(f"centroid.projection onto coil dir = {proj:+.1f} cm  (cos={cos:+.2f})  "
                   f"-> {'POINTS TOWARD coil (physical)' if proj > 0 else 'does NOT point toward coil (!)'}")
    print("\n".join("  " + m for m in msg))

    # --- projections --------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    vmax = imp.max()
    vmin = max(imp[imp > 0].min(), vmax * 1e-4)

    top = imp.max(axis=2).T                      # (ny,nx): view down z
    im0 = axes[0].pcolormesh(xc, yc, top, norm=LogNorm(vmin, vmax), cmap="magma", shading="auto")
    axes[0].set_xlabel("X (cm)"); axes[0].set_ylabel("Y (cm)")
    axes[0].set_title("Top View: Max Importance Over Z")
    axes[0].set_aspect("equal")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label="Adjoint Flux (Importance)")

    side = imp.max(axis=1).T                      # (nz,nx): view along y
    im1 = axes[1].pcolormesh(xc, zc, side, norm=LogNorm(vmin, vmax), cmap="magma", shading="auto")
    axes[1].set_xlabel("X (cm)"); axes[1].set_ylabel("Z (cm)")
    axes[1].set_title("Side View: Max Importance Over Y")
    axes[1].set_aspect("equal")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label="Adjoint Flux (Importance)")

    # overlays: coil (red star) + importance centroid (cyan +)
    if args.coil_centroid is not None:
        c = args.coil_centroid
        axes[0].plot(c[0], c[1], "*", color="red", ms=16, mec="white", label="Coil")
        axes[1].plot(c[0], c[2], "*", color="red", ms=16, mec="white", label="Coil")
    axes[0].plot(cen[0], cen[1], "+", color="cyan", ms=14, mew=2.5, label="Importance Centroid")
    axes[1].plot(cen[0], cen[2], "+", color="cyan", ms=14, mew=2.5)
    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.7)

    title = f"Adjoint Magnet-Importance Map --- Coil {coils}, {resp.capitalize()} Response"
    if bias is not None:
        title += f"   (Centroid $\\cdot$ Coil Dir = {bias[0]:+.0f} cm)"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out = args.out or str(Path(args.npz).with_suffix("")).replace("adjoint_importance_", "fig_adjoint_") + ".png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("  wrote", out)


if __name__ == "__main__":
    main()
