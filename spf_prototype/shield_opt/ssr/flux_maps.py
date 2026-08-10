#!/usr/bin/env python
"""Poloidal-toroidal load maps from frozen surface sources: the first-wall load (NWL), the shield-exit
load (leakage / transmitted wall loading), and their transmission ratio -- the three-panel figure
Peterson asked for, showing how the blanket smears the peaking before it reaches the magnets.

Both maps are energy-weighted and area-normalized (VMEC surface metric), so they are physical loads,
not angular-binning artifacts. Cross-talk (the peak shifting in (phi,theta) from wall to exit) reads
straight off the difference of the two panels.

Usage: flux_maps.py FIRSTWALL_SRC.h5 SHIELDEXIT_SRC.h5 OUT.png [--ntor 24 --npol 24 --wout PATH]
"""
import argparse
import numpy as np

from patch_rebin import load_surface_source, rebin_load, patch_edges


def maps(fw_src, exit_src, ntor, npol, nfp=4, wout=None):
    kw = dict(nfp=nfp, weight="energy", reduce="sum", wout_path=wout, area_normalize=wout is not None)
    nwl, R0 = rebin_load(load_surface_source(fw_src), ntor, npol, **kw)
    leak, _ = rebin_load(load_surface_source(exit_src), ntor, npol, **kw)
    with np.errstate(divide="ignore", invalid="ignore"):
        trans = np.where(nwl > 0, leak / nwl, 0.0)
    return nwl, leak, trans


def plot(nwl, leak, trans, out, nfp=4):
    import smplotlib  # noqa
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    te, pe = patch_edges(nwl.shape[0], nwl.shape[1], nfp)
    tor = np.degrees(0.5 * (te[:-1] + te[1:])); pol = np.degrees(0.5 * (pe[:-1] + pe[1:]))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, M, ttl, cm in ((ax[0], nwl, "First-Wall Load (NWL)", "inferno"),
                          (ax[1], leak, "Shield-Exit Load (Leakage)", "inferno"),
                          (ax[2], trans, "Transmission Ratio", "viridis")):
        im = a.pcolormesh(tor, pol, M.T, cmap=cm, shading="auto")
        a.set_xlabel("Toroidal Angle (deg, one period)"); a.set_ylabel("Poloidal Angle (deg)")
        a.set_title(ttl); fig.colorbar(im, ax=a, fraction=0.046, pad=0.04)
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("fw_src"); ap.add_argument("exit_src"); ap.add_argument("out")
    ap.add_argument("--ntor", type=int, default=24); ap.add_argument("--npol", type=int, default=24)
    ap.add_argument("--nfp", type=int, default=4); ap.add_argument("--wout", default=None)
    a = ap.parse_args()
    nwl, leak, trans = maps(a.fw_src, a.exit_src, a.ntor, a.npol, nfp=a.nfp, wout=a.wout)
    plot(nwl, leak, trans, a.out, nfp=a.nfp)
    print(f"wrote {a.out}  peak transmission={trans.max():.3f}  "
          f"NWL peak@{np.unravel_index(np.argmax(nwl), nwl.shape)}  "
          f"leak peak@{np.unravel_index(np.argmax(leak), leak.shape)} (cross-talk offset)")
