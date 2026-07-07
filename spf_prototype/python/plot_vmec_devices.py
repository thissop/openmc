#!/usr/bin/env python
"""Per-device figures from the VMEC fluxmaps (smplotlib): a 3D wireframe (coil filaments + LCFS) and
three poloidal cross-sections (nested VMEC flux surfaces at three toroidal angles across a half field
period). One figure per device -> figs/vmec_devices/. Plus a montage of all wireframes.
"""
import sys
import glob
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
try:
    import smplotlib  # noqa: F401
except Exception:
    pass
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sweep"))
import quasr_loader as ql  # noqa: E402
DATA = HERE.parent / "data"
FIGD = HERE.parent / "figs" / "vmec_devices"; FIGD.mkdir(parents=True, exist_ok=True)
BLUE = "#2a6fdb"


def one(npz):
    ID = int(re.search(r"quasr(\d+)_vmec", npz.name).group(1))
    f = dict(np.load(npz))
    R, Z, phi = f["R"], f["Z"], f["phi"]; nfp = int(f["nfp"])
    nr, nt, nz = R.shape
    fig = plt.figure(figsize=(16, 4.2))

    # --- panel 1: 3D wireframe (coils + LCFS) ---
    ax = fig.add_subplot(1, 4, 1, projection="3d")
    try:
        dev = ql.load_device(ID)
        for poly, _ in dev.coils:
            p = np.vstack([poly, poly[:1]])
            ax.plot(p[:, 0], p[:, 1], p[:, 2], color=BLUE, lw=0.6, alpha=0.55)
    except Exception:
        pass
    Rl, Zl, Pl = R[-1], Z[-1], phi[-1]
    Xs = Rl * np.cos(Pl); Ys = Rl * np.sin(Pl)
    ax.plot_surface(Xs, Ys, Zl, color="#e0782a", alpha=0.85, rcount=24, ccount=48,
                    linewidth=0, antialiased=True, shade=True)
    ax.set_axis_off(); ax.set_box_aspect((1, 1, 0.55)); ax.view_init(elev=26, azim=40)
    ax.set_title(f"{ID}  (nfp {nfp})", fontsize=10)

    # --- panels 2-4: nested flux-surface cross-sections at 3 toroidal angles ---
    ks = [0, int(round(nz / (4 * nfp))), int(round(nz / (2 * nfp)))]
    zang = [f"$\\zeta=0$", f"$\\zeta=\\pi/2N_{{fp}}$", f"$\\zeta=\\pi/N_{{fp}}$"]
    cmap = plt.cm.viridis
    for j, k in enumerate(ks):
        ax = fig.add_subplot(1, 4, 2 + j)
        for i in range(nr):
            c = cmap(i / (nr - 1))
            ax.plot(np.append(R[i, :, k], R[i, 0, k]), np.append(Z[i, :, k], Z[i, 0, k]),
                    color=c, lw=0.7)
        ax.set_aspect("equal"); ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        ax.set_title(zang[j], fontsize=10)
    fig.suptitle(f"VMEC equilibrium, device {ID}: coils + plasma boundary and nested flux surfaces", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGD / f"device_{ID}.png", dpi=170, bbox_inches="tight"); plt.close(fig)
    return ID


def main():
    maps = sorted(DATA.glob("quasr*_vmec_fluxmap.npz"))
    print(f"{len(maps)} VMEC fluxmaps found")
    ids = [one(m) for m in maps]
    # montage of wireframes
    n = len(ids)
    if n:
        cols = min(5, n); rows = int(np.ceil(n / cols))
        fig = plt.figure(figsize=(3.2 * cols, 3.0 * rows))
        for idx, m in enumerate(maps):
            ID = int(re.search(r"quasr(\d+)_vmec", m.name).group(1))
            f = dict(np.load(m)); R, Z, phi = f["R"], f["Z"], f["phi"]; nfp = int(f["nfp"])
            ax = fig.add_subplot(rows, cols, idx + 1, projection="3d")
            try:
                dev = ql.load_device(ID)
                for poly, _ in dev.coils:
                    p = np.vstack([poly, poly[:1]]); ax.plot(p[:, 0], p[:, 1], p[:, 2], color=BLUE, lw=0.4, alpha=0.8)
            except Exception:
                pass
            Rl, Zl, Pl = R[-1], Z[-1], phi[-1]
            ax.plot_surface(Rl*np.cos(Pl), Rl*np.sin(Pl), Zl, color="#e0782a", alpha=0.85,
                            rcount=20, ccount=40, linewidth=0, antialiased=True)
            ax.set_axis_off(); ax.set_box_aspect((1, 1, 0.55)); ax.view_init(elev=26, azim=40)
            ax.set_title(f"{ID} (nfp{nfp})", fontsize=8)
        fig.suptitle("VMEC-solved QUASR devices: coils + plasma boundary", y=1.01)
        fig.tight_layout(); fig.savefig(FIGD / "montage_wireframes.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
    print(f"wrote {n} device figures + montage -> figs/vmec_devices/")


if __name__ == "__main__":
    main()
