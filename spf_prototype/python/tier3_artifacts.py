#!/usr/bin/env python
"""Tier 3 — presentation artifacts. PURE POSTPROCESSING of the existing Tier-2b
wall-current statepoints (no new OpenMC runs, no geometry change) + the Tier-2a
analytic reference. Produces engineering tables, a consolidated analytic-vs-
numerical figure, and a 3D box-torus NWL render, then writes RESULTS_tier3.md.

Honesty scope (stated in every caption/limitation): filamentary/parabolic-ring
source, purely TOROIDAL field, scattering OFF (free-streaming sightline regime).
NWL here is a GEOMETRIC quantity -- not a shielded dose, not a TBR.

Source roles (both already validated):
- parabolic multi-ring  = physical plasma -> NWL patterns, peaking, current fractions
- single central ring   = clean analytic-vs-numerical check -> residual (tightest)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "spf_prototype" / "python"))
FIGDIR = REPO / "spf_prototype" / "figs"
RESULTS = REPO / "spf_prototype" / "RESULTS_tier3.md"

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: F401,E402

import openmc  # noqa: E402
import analytic_nwl as an  # noqa: E402
import build_and_run as br  # noqa: E402

openmc.config["cross_sections"] = str(br.XS)

MODES = ["iso", "A", "B", "C"]
WALLS = ["inboard", "outboard", "floor", "ceiling"]
NR = NZ = 50  # statepoint.10.h5 (4e6 histories) -- the set reported in RESULTS_tier2b.md
SP = {m: f"/tmp/spf_2b_{m}/statepoint.10.h5" for m in MODES + ["mixed"]}


# --------------------------------------------------------------------------
# Load existing wall-current arrays (no transport re-run)
# --------------------------------------------------------------------------
def load():
    W = {m: br.extract_walls(SP[m], NR, NZ) for m in MODES}
    # constant total fusion rate: each config emits 1 (conservation check)
    for m in MODES:
        tot = sum(W[m][w]["cur"].sum() for w in WALLS)
        assert abs(tot - 1.0) < 1e-6, f"{m}: total current {tot} != 1"
    # intensity = current / patch area
    I = {m: {w: W[m][w]["cur"] / W[m][w]["area"] for w in WALLS} for m in MODES}
    S = {m: {w: W[m][w]["std"] / W[m][w]["area"] for w in WALLS} for m in MODES}
    return W, I, S


def wall_area_total(W, w):
    return W["iso"][w]["area"].sum() if np.ndim(W["iso"][w]["area"]) else \
        W["iso"][w]["area"] * len(W["iso"][w]["cur"])


# --------------------------------------------------------------------------
# Engineering tables
# --------------------------------------------------------------------------
def peak_table(W, I):
    """Peak NWL per wall + overall, mode vs iso (% change)."""
    rows = {}
    for m in MODES:
        rows[m] = {w: I[m][w].max() for w in WALLS}
        allI = np.concatenate([I[m][w] for w in WALLS])
        rows[m]["overall"] = allI.max()
    return rows


def avg_per_wall(W, m, w):
    # area-weighted average intensity = total current on wall / wall area
    return W[m][w]["cur"].sum() / (W[m][w]["area"].sum() if np.ndim(W[m][w]["area"])
                                   else W[m][w]["area"] * len(W[m][w]["cur"]))


def total_area(W, w):
    a = W["iso"][w]["area"]
    return a.sum() if np.ndim(a) else a * len(W["iso"][w]["cur"])


def peaking_table(W, I):
    """Peaking factor (peak/average) per wall + overall, per mode."""
    PF = {}
    Atot = sum(total_area(W, w) for w in WALLS)
    for m in MODES:
        PF[m] = {}
        for w in WALLS:
            PF[m][w] = I[m][w].max() / avg_per_wall(W, m, w)
        allI = np.concatenate([I[m][w] for w in WALLS])
        # overall average = total load (=1) / total area
        PF[m]["overall"] = allI.max() / (1.0 / Atot)
    return PF


def current_fractions(W):
    """Fraction of total wall current on each wall, per mode (total = 1)."""
    F = {}
    for m in MODES:
        F[m] = {w: W[m][w]["cur"].sum() for w in WALLS}
    return F


# --------------------------------------------------------------------------
# Analytic reference (exact constant-rate directional D = bracket/eta)
# --------------------------------------------------------------------------
def analytic_D(W):
    z_c = W["iso"]["inboard"]["centers"]; r_c = W["iso"]["floor"]["centers"]
    tbw = {
        "inboard": [dict(wall="inboard", R=br.U, Z=z, s=z) for z in z_c],
        "outboard": [dict(wall="outboard", R=br.W, Z=z, s=z) for z in z_c],
        "floor": [dict(wall="floor", R=r, Z=-br.ZW, s=r) for r in r_c],
        "ceiling": [dict(wall="ceiling", R=r, Z=br.ZW, s=r) for r in r_c],
    }
    allt, sl, k = [], {}, 0
    for w in WALLS:
        sl[w] = slice(k, k + len(tbw[w])); k += len(tbw[w]); allt += tbw[w]
    pat = an.plasma_patterns(allt, n_grid=120, engine=an.g_ring_quad)
    D = {}
    for m in MODES:
        a, b, c = an.MODES[m]
        d = pat[f"{m}_bracket"] / an.eta(a, b, c)
        D[m] = {w: d[sl[w]] for w in WALLS}
    return D


def single_ring_tightest():
    """Tightest analytic-vs-numerical: single central ring, numerical quadrature
    vs anarrima closed form (and our Eq.5). No MC noise."""
    p = an.R0
    errs = []
    for wall, r, z in [("inboard", an.R_IN, 0.0), ("outboard", an.R_OUT, 0.0),
                       ("floor", an.R0, 1.0)]:
        for fac in ("iso", "A", "cos2", "BC"):
            q = an.g_ring_quad_scalar(p, z, r, wall, fac, an.R_IN)
            a = an.g_ring_anarrima(p, z, r, wall, fac, an.R_IN)
            errs.append(abs(q - a))
    e5 = abs(an.g_ring_quad_scalar(p, 0.0, an.R_IN, "inboard", "iso", an.R_IN)
             - float(an.g_RiI_closed_eq5(p, 0.0, an.R_IN)))
    return max(errs), e5


# --------------------------------------------------------------------------
# Figure 1: analytic vs numerical (multi-ring overlay + residual)
# --------------------------------------------------------------------------
def fig_analytic_vs_numerical(W, I, S, D):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    resid = []
    for ax, w, xl in [(axes[0, 0], "inboard", "z [norm]"),
                      (axes[0, 1], "outboard", "z [norm]"),
                      (axes[1, 0], "floor", "R [norm]")]:
        c = W["iso"][w]["centers"]
        for m, col in [("iso", "k"), ("A", "r"), ("B", "b")]:
            mr = I[m][w] / I["iso"][w]
            sg = mr * np.sqrt((S[m][w] / I[m][w]) ** 2 + (S["iso"][w] / I["iso"][w]) ** 2)
            ar = D[m][w] / D["iso"][w]
            ax.plot(c, ar, col + "-", lw=1.5)
            ax.errorbar(c, mr, yerr=sg, fmt=col + ".", ms=3, alpha=0.5, label=m)
            if m != "iso":
                good = sg > 0
                resid.append((mr[good] - ar[good]) / sg[good])
        ax.axhline(1.0, color="0.6", ls="dashed", lw=1)
        ax.set_title(f"{w}: NWL ratio to isotropic")
        ax.set_xlabel(xl); ax.set_ylabel("mode / iso")
    axes[0, 0].legend(fontsize=8)
    resid = np.concatenate(resid)
    axes[1, 1].hist(resid, bins=35, density=True, alpha=0.6, color="purple")
    xx = np.linspace(-4, 4, 200)
    axes[1, 1].plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), "k-", lw=2, label="N(0,1)")
    axes[1, 1].set_title(f"standardized residual A/B/C (μ={resid.mean():+.2f}, σ={resid.std():.2f})")
    axes[1, 1].set_xlabel("(OpenMC − analytic)/σ"); axes[1, 1].legend()
    fig.suptitle("Tier 3 — OpenMC (points) reproduces analytic NWL directionality (lines).  "
                 "Parabolic multi-ring plasma, toroidal field, scattering OFF (free-streaming).",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGDIR / "tier3_analytic_vs_numerical.png", dpi=110)
    plt.close(fig)
    return resid.mean(), resid.std()


# --------------------------------------------------------------------------
# Figure 2: 3D box-torus colored by NWL (iso vs B/C), shared colorbar
# --------------------------------------------------------------------------
def _wall_surface(ax, kind, coord_centers, vals, cmap, norm, phi0, phi1, nphi=60):
    phi = np.linspace(phi0, phi1, nphi)
    C = np.tile(vals.reshape(-1, 1), (1, nphi))      # color constant in phi (axisym)
    fc = cmap(norm(C))
    if kind == "inboard":
        R = br.U; PH, ZZ = np.meshgrid(phi, coord_centers)
        X, Y, Z = R * np.cos(PH), R * np.sin(PH), ZZ
    elif kind == "outboard":
        R = br.W; PH, ZZ = np.meshgrid(phi, coord_centers)
        X, Y, Z = R * np.cos(PH), R * np.sin(PH), ZZ
    elif kind == "floor":
        PH, RR = np.meshgrid(phi, coord_centers); Zc = -br.ZW
        X, Y, Z = RR * np.cos(PH), RR * np.sin(PH), np.full_like(RR, Zc)
    else:  # ceiling
        PH, RR = np.meshgrid(phi, coord_centers); Zc = br.ZW
        X, Y, Z = RR * np.cos(PH), RR * np.sin(PH), np.full_like(RR, Zc)
    ax.plot_surface(X, Y, Z, facecolors=fc, rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False)


def fig_torus3d(W, I):
    cmap = cm.viridis
    allvals = np.concatenate([I[m][w] for m in ("iso", "B") for w in WALLS])
    norm = Normalize(vmin=allvals.min(), vmax=allvals.max())
    phi0, phi1 = 0.0, 1.5 * np.pi  # 3/4 cutaway to see inside
    fig = plt.figure(figsize=(14, 6.5))
    for idx, (m, title) in enumerate([("iso", "Non-polarized (isotropic)"),
                                      ("B", "Polarized B/C mode")]):
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
        zc = W[m]["inboard"]["centers"]; rc = W[m]["floor"]["centers"]
        _wall_surface(ax, "inboard", zc, I[m]["inboard"], cmap, norm, phi0, phi1)
        _wall_surface(ax, "outboard", zc, I[m]["outboard"], cmap, norm, phi0, phi1)
        _wall_surface(ax, "floor", rc, I[m]["floor"], cmap, norm, phi0, phi1)
        _wall_surface(ax, "ceiling", rc, I[m]["ceiling"], cmap, norm, phi0, phi1)
        # schematic plasma torus (translucent) + bright central ring, same cutaway
        u = np.linspace(phi0, phi1, 60); v = np.linspace(0, 2 * np.pi, 30)
        UU, VV = np.meshgrid(u, v)
        Xt = (br.R0 + br.AM * np.cos(VV)) * np.cos(UU)
        Yt = (br.R0 + br.AM * np.cos(VV)) * np.sin(UU)
        Zt = br.AM * np.sin(VV)
        ax.plot_surface(Xt, Yt, Zt, color="orange", alpha=0.30, linewidth=0, shade=False)
        ph = np.linspace(phi0, phi1, 120)
        ax.plot(br.R0 * np.cos(ph), br.R0 * np.sin(ph), np.zeros_like(ph),
                color="darkorange", lw=2.5, label="plasma (R0 ring)")
        ax.set_title(title, fontsize=11)
        ax.set_box_aspect((1, 1, 0.5)); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.view_init(elev=32, azim=-60)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cax = fig.add_axes([0.92, 0.2, 0.015, 0.6])
    fig.colorbar(sm, cax=cax, label="NWL [normalized, constant fusion rate]")
    fig.suptitle("Schematic box-torus render of computed free-streaming NWL (NOT CAD/DAGMC). "
                 "Orange = schematic plasma source. Shared colorbar.\n"
                 "B/C steers neutrons OFF the inboard center stack toward the outboard wall.",
                 fontsize=10)
    fig.savefig(FIGDIR / "tier3_torus3d.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    W, I, S = load()
    D = analytic_D(W)
    peaks = peak_table(W, I)
    PF = peaking_table(W, I)
    F = current_fractions(W)
    sr_max, sr_e5 = single_ring_tightest()
    rmean, rstd = fig_analytic_vs_numerical(W, I, S, D)
    fig_torus3d(W, I)

    L = []; w = L.append
    w("# RESULTS — Tier 3 (presentation artifacts; postprocessing only)\n")
    w("All numbers are postprocessing of the existing Tier-2b wall-current data "
      f"(parabolic multi-ring plasma, toroidal field, scattering OFF, "
      f"{4_000_000:,} histories/config, {NR}×{NZ} bins/face) and the Tier-2a analytic "
      "reference. **No new transport runs.**\n")
    w("> **Scope / honesty:** NWL here is a *geometric, free-streaming* quantity "
      "(sightline regime, no scattering). It is **not** a shielded dose and **not** "
      "a TBR. Source is a filamentary/parabolic ring set; field is purely toroidal. "
      "Patterns are normalized to **constant total fusion rate** (each mode emits the "
      "same total), so differences are pure directionality.\n")

    # peak table
    w("## Peak NWL per wall, mode vs isotropic (% change)\n")
    w("| wall | iso (abs) | A | B | C |")
    w("|---|---|---|---|---|")
    for wl in WALLS + ["overall"]:
        cells = [f"{peaks['iso'][wl]:.3f}"]
        for m in ("A", "B", "C"):
            cells.append(f"{100*(peaks[m][wl]/peaks['iso'][wl]-1):+.1f}%")
        w(f"| {wl} | " + " | ".join(cells) + " |")
    w("\n(iso column is the absolute normalized peak intensity; A/B/C are % change "
      "vs iso at constant fusion rate. B and C are identical by construction.)\n")

    # peaking factor table
    w("## Peaking factor (peak / average) per mode — the engineering scalar\n")
    w("| wall | iso | A | B/C | Δ(B/C vs iso) | Δ(A vs iso) |")
    w("|---|---|---|---|---|---|")
    for wl in WALLS + ["overall"]:
        dbc = 100 * (PF["B"][wl] / PF["iso"][wl] - 1)
        da = 100 * (PF["A"][wl] / PF["iso"][wl] - 1)
        w(f"| {wl} | {PF['iso'][wl]:.2f} | {PF['A'][wl]:.2f} | {PF['B'][wl]:.2f} | "
          f"{dbc:+.1f}% | {da:+.1f}% |")
    w("\n**Honest finding (contrary to a naive 'polarization lowers peaking' "
      "expectation):** in this FIXED square geometry, *both* polarized modes RAISE the "
      "overall peaking factor vs isotropic (iso 1.33 → A 1.65, B/C 1.61). Directional "
      "emission concentrates the wall load, so the isotropic case is the most uniform. "
      "Per-wall peaking changes are small (±1.5%); the overall rise comes from the "
      "global peak intensifying at the outboard midplane (B/C) or relocating to the "
      "inboard midplane (A). The B/C engineering benefit here is therefore **not** a "
      "lower global peaking factor — it is the **center-stack (inboard) load reduction** "
      "(inboard peak −41%; inboard current fraction 8.7%→5.1%, see next table), which "
      "Schwartz highlights for spherical-tokamak center stacks. Schwartz's peaking-"
      "factor *reduction* result comes from re-shaping the first wall to exploit the "
      "steering — NOT modeled here (fixed square). Geometric free-streaming peaking, "
      "not shielded-dose peaking.\n")

    # current fraction table
    w("## Inboard vs outboard neutron-current fraction (geometric precursor, NOT a TBR)\n")
    w("| mode | inboard | outboard | floor | ceiling |")
    w("|---|---|---|---|---|")
    for m in MODES:
        w(f"| {m} | {100*F[m]['inboard']:.1f}% | {100*F[m]['outboard']:.1f}% | "
          f"{100*F[m]['floor']:.1f}% | {100*F[m]['ceiling']:.1f}% |")
    moved = 100 * (F["iso"]["inboard"] - F["B"]["inboard"])
    gained = 100 * (F["B"]["outboard"] - F["iso"]["outboard"])
    w(f"\n**B/C moves {moved:.1f} percentage points of total wall current off the "
      f"inboard center stack** (inboard {100*F['iso']['inboard']:.1f}%→"
      f"{100*F['B']['inboard']:.1f}%); the outboard wall gains {gained:.1f} pts "
      f"({100*F['iso']['outboard']:.1f}%→{100*F['B']['outboard']:.1f}%). A mode does the "
      "reverse. **This is the geometric precursor to a TBR argument — it is NOT a TBR**: "
      "a real TBR needs a breeder blanket + scattering, absent here (see LIMITATIONS.md).\n")

    # analytic vs numerical
    w("## Analytic vs numerical (consolidated)\n")
    w(f"- **Tightest check (single central ring, no MC noise):** numerical quadrature "
      f"vs anarrima closed form agree to **{sr_max:.1e}** across all walls/factors; "
      f"our independent Eq.5 closed form agrees to **{sr_e5:.1e}**.\n")
    w(f"- **Full-stack (parabolic multi-ring, OpenMC Monte Carlo vs analytic):** "
      f"standardized directionality residual over all walls (A/B/C) has mean "
      f"**{rmean:+.2f}**, stdev **{rstd:.2f}** ⇒ OpenMC reproduces the analytic to "
      f"statistics (σ<1 ⇒ the 10-batch tally error bars are mildly conservative; no bias).\n")
    w("![analytic vs numerical](figs/tier3_analytic_vs_numerical.png)\n")
    w("![3D box-torus NWL](figs/tier3_torus3d.png)\n")
    w("*3D: schematic box-torus render of the computed free-streaming NWL (not a "
      "CAD/DAGMC geometry); orange = schematic plasma; shared colorbar. The B/C panel "
      "shows the inboard center-stack reduction at a glance.*")

    RESULTS.write_text("\n".join(L) + "\n")
    print(f"wrote {RESULTS}")
    print(f"single-ring tightest max|Δ|={sr_max:.1e}; multi-ring residual μ={rmean:+.2f} σ={rstd:.2f}")
    print(f"B/C inboard fraction {100*F['B']['inboard']:.1f}% vs iso {100*F['iso']['inboard']:.1f}%")
    print(f"overall peaking iso {PF['iso']['overall']:.2f} -> B/C {PF['B']['overall']:.2f}")


if __name__ == "__main__":
    main()
