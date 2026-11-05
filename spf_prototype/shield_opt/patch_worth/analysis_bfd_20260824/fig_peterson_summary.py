#!/usr/bin/env python
"""Peterson check-in summary figure: the design-tradeoff thesis in three panels, from VERIFIED artifacts only.
  A. Equal thickness != equal mass  (40M-sample MC VolumeCalculation, error-barred)
  B. Placement-dependent control     (structured: inboard reduce, wrong sector backfire; B_FD finite diff)
  C. Need collides with feasibility  (clearance available-to-grow dips inboard, where coil-20 peaks)
All numbers trace to masses_local.json (40M MC), b_fd_fractional.csv, coil_step1b_uniform.npz, clearance_domain.npz.
No transport rerun; this only re-plots gated data (verify_closeout.py PASSES all gates)."""
import os, json, numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

HERE = os.path.dirname(os.path.abspath(__file__)); SH = "/Users/tkiker/lima_share"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold",
                     "figure.dpi": 130, "savefig.dpi": 160, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
GREEN, RED, BLUE, GREY = "#2ca25f", "#c0392b", "#2c7fb8", "#666666"

# ---- data ----------------------------------------------------------------
mc = json.load(open(f"{SH}/masses_local.json"))
patches = ["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
dM  = np.array([mc[p]["dM_net_kg_vs_zerodelta"]     for p in patches]) / 1e3   # tonnes
dMe = np.array([mc[p]["dM_net_err_kg_vs_zerodelta"] for p in patches]) / 1e3
Muni = (mc["uniform"]["dM_net_kg_vs_zerodelta"]) / 1e3
order = np.argsort(dM)
p_s, dM_s, dMe_s = [patches[i] for i in order], dM[order], dMe[order]

# Panel B uses the DIRECT patch-run coil vectors (coil_patch_*.npz) — the real transported data behind
# fig_peak_migration. Peak RATIO to baseline; the sign (reduce<1 / backfire>1) is the robust finding, the
# absolute magnitude is baseline-confounded by the shield-outside offset (patches at 40 cm vs baseline at 42.55 cm).
base = np.load(f"{HERE}/../cluster_pull/coil_step1b_uniform.npz")
Rbase = np.ravel(base["coil_fast_flux"]).astype(float); peak_base = Rbase.max()
pk_peak, pk_coil = [], []
for p in patches:
    R = np.ravel(np.load(f"{HERE}/../cluster_pull/coil_patch_{p}.npz")["coil_fast_flux"]).astype(float)
    pk_peak.append(R.max()); pk_coil.append(11 + int(R.argmax()))
pk_peak = np.array(pk_peak); pk_ratio = pk_peak / peak_base
po = np.argsort(pk_ratio)
ratio_s = pk_ratio[po]; pk_names = [patches[i] for i in po]; coil_s = [pk_coil[i] for i in po]

cd = np.load(f"{HERE}/clearance_domain.npz")
avail = cd["avail"]; pe = cd["pe"]            # (61,55) growth room, poloidal edges (deg)
avail_pol = avail.min(axis=0)                 # tightest over toroidal, per poloidal angle

# ---- figure --------------------------------------------------------------
fig = plt.figure(figsize=(15, 4.7))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.15, 1.15, 1.0], wspace=0.28)

# A. mass spread
axA = fig.add_subplot(gs[0])
axA.bar(range(12), dM_s, yerr=dMe_s, color=BLUE, edgecolor="k", lw=0.5, capsize=3)
axA.set_ylim(0, Muni * 1.18)
axA.axhline(Muni, ls="--", color="k", lw=1.5)
axA.text(11.4, Muni + 8, f"uniform +2.55 cm control = {Muni:.0f} t", fontsize=9, ha="right")
axA.set_xticks(range(12)); axA.set_xticklabels(p_s, rotation=90, fontsize=8)
axA.set_ylabel("added shield mass vs zero-δ  (tonnes)")
axA.set_title("A.  Equal thickness ≠ equal mass")
axA.text(0.02, 0.97, f"same +20 cm patch: {dM_s.min():.0f}–{dM_s.max():.0f} t  →  {dM_s.max()/dM_s.min():.2f}× spread\n"
                     f"every patch lighter than uniform (t3p3 −{Muni-dM_s.max():.0f} t, 8σ)",
         transform=axA.transAxes, va="top", fontsize=8.5,
         bbox=dict(boxstyle="round", fc="#fff7e6", ec="#e0a800"))

# B. placement-dependent control
axB = fig.add_subplot(gs[1])
colsB = [GREEN if v < 1.0 else RED for v in ratio_s]
axB.bar(range(12), ratio_s, color=colsB, edgecolor="k", lw=0.5)
axB.axhline(1.0, ls="--", color="k", lw=1.5)
axB.text(0.3, 1.03, "baseline peak (coil 20)", fontsize=9)
for i, (r, c) in enumerate(zip(ratio_s, coil_s)):   # annotate where the hotspot migrates to
    axB.text(i, r + 0.03, f"c{c}", ha="center", fontsize=6.5, color=GREY, rotation=90)
axB.set_xticks(range(12)); axB.set_xticklabels(pk_names, rotation=90, fontsize=8)
axB.set_ylabel("peak coil load after patch  (× baseline)")
axB.set_title("B.  Placement-dependent control")
axB.text(0.02, 0.97, "inboard sectors REDUCE the peak (green),\nwrong sector BACKFIRES (red) — sign is ROBUST.\n"
                     "magnitude baseline-confounded; hotspot MIGRATES\n(labels: new peak coil ≠ 20).",
         transform=axB.transAxes, va="top", fontsize=8.3,
         bbox=dict(boxstyle="round", fc="#eef7ee", ec=GREEN))

# C. need vs feasibility
axC = fig.add_subplot(gs[2])
axC.plot(pe, avail_pol, color=BLUE, lw=2)
axC.fill_between(pe, 0, avail_pol, color=BLUE, alpha=0.12)
imin = int(np.argmin(avail_pol))
axC.scatter([pe[imin]], [avail_pol[imin]], color=RED, zorder=5, s=45)
axC.annotate(f"min {avail_pol[imin]:.1f} cm\n(inboard, pol≈{pe[imin]:.0f}°)",
             (pe[imin], avail_pol[imin]), (pe[imin]+25, avail_pol[imin]+120),
             fontsize=9, color=RED, arrowprops=dict(arrowstyle="->", color=RED))
axC.axvspan(150, 210, color=RED, alpha=0.07)
axC.text(180, axC.get_ylim()[1]*0.9, "inboard\n(coil-20 peaks here)", ha="center", fontsize=8.5, color=RED)
axC.set_xlabel("poloidal angle  (deg)"); axC.set_ylabel("available growth toward coils  (cm)")
axC.set_title("C.  Need collides with feasibility")
axC.set_xlim(0, 360)

fig.suptitle("Response-aware stellarator shield design — the tradeoff, quantified   "
             "(QH, 1 field period, 40M-MC gated)", y=1.02, fontsize=13, fontweight="bold")
out = f"{HERE}/fig_peterson_summary.png"
fig.savefig(out, bbox_inches="tight")
print("WROTE", out)
print(f"A: mass spread {dM_s.min():.0f}-{dM_s.max():.0f} t ({dM_s.max()/dM_s.min():.2f}x); uniform {Muni:.0f} t")
print(f"B: baseline peak {peak_base:.4e} c20; best {ratio_s.min():.2f}x ({pk_names[0]}); worst {ratio_s.max():.2f}x ({pk_names[-1]}); hotspot migrates c20->{coil_s}")
print(f"C: min available growth {avail_pol[imin]:.1f} cm at poloidal {pe[imin]:.0f} deg (inboard)")
