#!/usr/bin/env python
"""Results figure (smplotlib) for the real-device StellaratorSource transport on 803097: inboard and
outboard poloidal wall-load profiles for unpolarized / perpendicular / parallel, and the fractional
steering (mode - unpol)/unpol. Reads data/quasr803097_wallcurrent.npz.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
try:
    import smplotlib  # noqa: F401
except Exception:
    pass
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = HERE.parent / "data"
FIGD = HERE.parent / "figs" / "stellarator_source"; FIGD.mkdir(parents=True, exist_ok=True)
ID = 803097
COL = {"unpol": "k", "perp": "#2a6fdb", "par": "#e2711d"}


def main():
    f = dict(np.load(D / f"quasr{ID}_wallcurrent.npz"))
    zc = f["zc"]
    mid = int(np.argmin(np.abs(zc)))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

    for mode in ("unpol", "perp", "par"):
        ax[0].plot(f[f"{mode}_inb"], zc, color=COL[mode], label=mode, lw=1.6)
        ax[1].plot(f[f"{mode}_out"], zc, color=COL[mode], label=mode, lw=1.6)
    ax[0].set_title("inboard wall (r-min)"); ax[1].set_title("outboard wall (r-max)")
    for a in ax[:2]:
        a.set_xlabel("neutron current  [arb.]"); a.set_ylabel("Z [cm]"); a.legend(frameon=False)

    u_in, u_out = f["unpol_inb"], f["unpol_out"]
    for mode in ("perp", "par"):
        ax[2].plot(100 * (f[f"{mode}_inb"] - u_in) / np.where(u_in > 0, u_in, np.nan), zc,
                   color=COL[mode], ls="-", lw=1.6, label=f"{mode} inboard")
        ax[2].plot(100 * (f[f"{mode}_out"] - u_out) / np.where(u_out > 0, u_out, np.nan), zc,
                   color=COL[mode], ls="--", lw=1.6, label=f"{mode} outboard")
    ax[2].axvline(0, color="grey", ls=":", lw=0.8)
    ax[2].set_xlabel("steering  (mode - unpol)/unpol  [%]"); ax[2].set_ylabel("Z [cm]")
    ax[2].set_title("fractional redistribution"); ax[2].legend(frameon=False, fontsize=8)

    mp = {}
    for s, lab in [("inb", "inboard"), ("out", "outboard")]:
        u = f[f"unpol_{s}"][mid]
        mp[lab] = {m: 100 * (f[f"{m}_{s}"][mid] - u) / u for m in ("perp", "par")}
    fig.suptitle(f"SPF on the real 803097 VMEC equilibrium (native StellaratorSource, transport).  "
                 f"Midplane perp: {mp['inboard']['perp']:+.0f}% inboard, {mp['outboard']['perp']:+.0f}% outboard",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(FIGD / f"realdev_{ID}_wallsteering.png", dpi=175, bbox_inches="tight")
    print(f"midplane steering: {mp}")
    print(f"plot -> figs/stellarator_source/realdev_{ID}_wallsteering.png")


if __name__ == "__main__":
    main()
