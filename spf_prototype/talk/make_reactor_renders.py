#!/usr/bin/env python
"""Black-and-white wireframe reactor renders for the SPF Beamer deck.

House style (per T.K.): `import smplotlib` is the VERY FIRST plotting-related
import so the whole deck shares the serif Super-Mongo look. Works entirely from
committed local data (VMEC wout + MAKEGRID coils); no OpenMC, no DAGMC, no net.

Design: the plasma last-closed-flux-surface is drawn as a thin BLACK wireframe
on white (see-through, so the shaped cross-section reads); coils are RED tube
volumes wrapping outside the plasma. No in-figure titles or legends (those live
in the LaTeX frame). Art fills the frame with minimal whitespace.

Outputs (dpi>=200, tight bbox):
  talk/figs/render_qa_wire.png    QA plasma wireframe, no coils          (nfp=2)
  talk/figs/render_qa_coils.png   QA plasma wireframe + red coil volumes (nfp=2)
  talk/figs/render_qh_wire.png    QH plasma wireframe, no coils          (nfp=4)
  talk/figs/render_qh_coils.png   QH plasma wireframe + red coil volumes (nfp=4)

Run:  python talk/make_reactor_renders.py
"""
from __future__ import annotations

import smplotlib  # noqa: F401  house serif style — MUST precede pyplot use
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
PROTO = HERE.parent
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

PLASMA_COLOR = "black"      # see-through wireframe on white
COIL_COLOR = "#c0392b"      # red magnet volumes
CUTAWAY_DEG = 90.0          # angular width of the removed wedge -> 3/4 cutaway
THETA_STRIDE = 6            # draw every 6th poloidal line
PHI_STRIDE = 4             # draw every 4th toroidal line


# --------------------------------------------------------------------------- #
# VMEC boundary reconstruction
# --------------------------------------------------------------------------- #
def lcfs_surface(wout_path, n_theta=121, n_phi=241, cutaway_deg=CUTAWAY_DEG):
    """Reconstruct the last-flux-surface Cartesian grid from a VMEC wout.

    Returns (X, Y, Z, R, nfp) on a (n_theta, n_phi) grid. A toroidal wedge of
    width `cutaway_deg` is removed (phi runs continuously past the gap) so the
    shaped cross-section is visible — a clean 3/4 cutaway.
    """
    from scipy.io import netcdf_file

    f = netcdf_file(str(wout_path), mmap=False)
    xm = f.variables["xm"].data
    xn = f.variables["xn"].data  # already includes the nfp factor (xn = n*nfp)
    rmnc = f.variables["rmnc"].data
    zmns = f.variables["zmns"].data
    rmns = f.variables["rmns"].data if "rmns" in f.variables else None
    zmnc = f.variables["zmnc"].data if "zmnc" in f.variables else None
    nfp = int(f.variables["nfp"].data)
    f.close()

    ns = rmnc.shape[0]
    s = ns - 1  # last flux surface

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=True)
    gap = np.deg2rad(cutaway_deg)
    phi = np.linspace(gap, 2.0 * np.pi, n_phi, endpoint=True)

    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    ang = xm[None, None, :] * TH[..., None] - xn[None, None, :] * PH[..., None]
    R = np.sum(rmnc[s][None, None, :] * np.cos(ang), axis=-1)
    Z = np.sum(zmns[s][None, None, :] * np.sin(ang), axis=-1)
    if rmns is not None:
        R += np.sum(rmns[s][None, None, :] * np.sin(ang), axis=-1)
    if zmnc is not None:
        Z += np.sum(zmnc[s][None, None, :] * np.cos(ang), axis=-1)

    X = R * np.cos(PH)
    Y = R * np.sin(PH)
    return X, Y, Z, R, nfp


# --------------------------------------------------------------------------- #
# MAKEGRID coils parser + stellarator-symmetry completion
# --------------------------------------------------------------------------- #
def parse_coils(coils_path):
    """Parse a MAKEGRID/focus 'coils.' file into a list of (M,3) polylines.

    Rows are 'x y z current' (meters). A coil closes on the row whose current is
    0.0; 'end' terminates the file.
    """
    coils = []
    cur = []
    with open(coils_path) as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            if t[0] in ("periods", "begin", "mirror"):
                continue
            if t[0] == "end":
                break
            try:
                x, y, z, current = (float(t[0]), float(t[1]),
                                    float(t[2]), float(t[3]))
            except (ValueError, IndexError):
                continue
            cur.append((x, y, z))
            if current == 0.0:
                if len(cur) >= 2:
                    coils.append(np.asarray(cur))
                cur = []
    if len(cur) >= 2:
        coils.append(np.asarray(cur))
    return coils


def complete_coils_stellsym(coils):
    """Fill the torus using stellarator symmetry (R,phi,Z)->(R,-phi,-Z), i.e.
    Cartesian (x,y,z)->(x,-y,-z). Some `coils.` exports store only the
    symmetry-unique half of each field period; without the mirror partners the
    coils cover only part of the torus. Partners that coincide with an existing
    coil (already-complete sets like QA) are de-duplicated, so this is a no-op
    when the file is already full.
    """
    centroids = [c.mean(0) for c in coils]
    tol = 1.0  # metres
    completed = list(coils)
    added = 0
    for c in coils:
        mirror = c.copy()
        mirror[:, 1] *= -1.0
        mirror[:, 2] *= -1.0
        mc = mirror.mean(0)
        if min(np.linalg.norm(mc - g) for g in centroids) > tol:
            completed.append(mirror)
            centroids.append(mc)
            added += 1
    return completed, added


def coil_diagnostics(coils):
    """Cylindrical R0 and radial spread of the whole coil cloud."""
    allpts = np.vstack(coils)
    Rcyl = np.hypot(allpts[:, 0], allpts[:, 1])
    return float(Rcyl.mean()), float(Rcyl.std()), float(Rcyl.min()), float(Rcyl.max())


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def render(device, wout_path, coils_path, out_png, show_coils):
    tag = f"{device}/{'coils' if show_coils else 'wire'}"
    print(f"\n===== {tag} =====")
    X, Y, Z, R, nfp = lcfs_surface(wout_path)

    lcfs_R0 = 0.5 * (R.max() + R.min())
    lcfs_a = 0.5 * (R.max() - R.min())
    print(f"[{tag}] LCFS  R0 = {lcfs_R0:.3f} m   a = {lcfs_a:.3f} m   (nfp={nfp})")

    coils = None
    if show_coils:
        coils = parse_coils(coils_path)
        coils, n_added = complete_coils_stellsym(coils)
        if n_added:
            print(f"[{tag}] stellarator-symmetry completion: +{n_added} mirror "
                  f"coils (now {len(coils)}) to wrap the full torus.")
        else:
            print(f"[{tag}] coil set already covers the full torus.")
        coil_R0, coil_spread, coil_Rmin, coil_Rmax = coil_diagnostics(coils)
        ratio = coil_R0 / lcfs_R0
        print(f"[{tag}] Coils R0 = {coil_R0:.3f} m (ratio {ratio:.2f})  "
              f"n={len(coils)}  Rmin={coil_Rmin:.2f} Rmax={coil_Rmax:.2f}")
        # Same scale-consistency gate as before; verified in-tolerance already.
        if not (0.4 <= ratio <= 3.0):
            factor = lcfs_R0 / coil_R0
            coils = [c * factor for c in coils]
            coil_R0, coil_spread, _, _ = coil_diagnostics(coils)
            print(f"[{tag}] WARNING: ratio {ratio:.2f} outside 0.4..3x — RESCALED "
                  f"coils by {factor:.4f}. New coil R0 = {coil_R0:.3f} m.")
        else:
            print(f"[{tag}] scale OK (ratio in 0.4..3x), no rescale.")
        loop = float(np.median([np.linalg.norm(c - c.mean(0), axis=1).mean()
                                for c in coils]))
        if loop < lcfs_a:
            coils = None
            print(f"[{tag}] DROPPED coils: median loop radius {loop:.2f} m < "
                  f"plasma minor radius {lcfs_a:.2f} m.")

    # ---- figure: fill the frame, no chrome --------------------------------- #
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Black see-through wireframe plasma (subsampled poloidal/toroidal lines).
    ax.plot_wireframe(
        X, Y, Z,
        rstride=THETA_STRIDE, cstride=PHI_STRIDE,
        color=PLASMA_COLOR, linewidth=0.6, alpha=0.9,
    )

    if coils is not None:
        for c in coils:
            ax.plot(c[:, 0], c[:, 1], c[:, 2],
                    color=COIL_COLOR, linewidth=4.5, solid_capstyle="round",
                    solid_joinstyle="round", alpha=0.95)

    # ---- equal geometry, snug limits --------------------------------------- #
    if coils is not None:
        allx = np.concatenate([X.ravel()] + [c[:, 0] for c in coils])
        ally = np.concatenate([Y.ravel()] + [c[:, 1] for c in coils])
        allz = np.concatenate([Z.ravel()] + [c[:, 2] for c in coils])
    else:
        allx, ally, allz = X.ravel(), Y.ravel(), Z.ravel()
    pad = 1.01
    cx = 0.5 * (allx.max() + allx.min())
    cy = 0.5 * (ally.max() + ally.min())
    cz = 0.5 * (allz.max() + allz.min())
    xr = (allx.max() - allx.min()) * pad
    yr = (ally.max() - ally.min()) * pad
    zr = (allz.max() - allz.min()) * pad
    ax.set_box_aspect((xr, yr, zr))
    ax.set_xlim(cx - xr / 2, cx + xr / 2)
    ax.set_ylim(cy - yr / 2, cy + yr / 2)
    ax.set_zlim(cz - zr / 2, cz + zr / 2)

    ax.view_init(elev=24, azim=35)
    ax.set_axis_off()
    # Zoom in so the flat torus fills the square frame (mpl3d reserves a cube).
    try:
        ax.set_box_aspect((xr, yr, zr), zoom=1.6)
    except TypeError:
        pass  # older mpl without zoom kwarg
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    fig.savefig(out_png, dpi=220, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print(f"[{tag}] wrote {out_png}  (coils {'ON' if coils is not None else 'OFF'})")
    return out_png


if __name__ == "__main__":
    jobs = [
        ("QA", "data/wout_QA_reactor.nc",
         "shield_opt/data/processed/Wechsung_QA24.coils"),
        ("QH", "data/wout_qh.nc",
         "shield_opt/data/processed/Wiedman_LandremanPaulQH.coils"),
    ]
    for device, wout, coilf in jobs:
        wpath = PROTO / wout
        cpath = PROTO / coilf
        render(device, wpath, cpath,
               FIGS / f"render_{device.lower()}_wire.png", show_coils=False)
        render(device, wpath, cpath,
               FIGS / f"render_{device.lower()}_coils.png", show_coils=True)
    print("\nDone.")
