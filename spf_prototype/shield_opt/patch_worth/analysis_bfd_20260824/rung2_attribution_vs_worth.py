#!/usr/bin/env python
"""RUNG 2 (charter 2026-08-24): compare the coil-20 ADJOINT IMPORTANCE psi_dagger(r) against the empirical
intervention effect (coarse B_FD) across MULTIPLE worth definitions -- not just coil-20.

adjoint_coil20_flat_P0.npz = coil-20 adjoint importance (flat response) on a 30x30x20 CARTESIAN mesh
(x,y,z in cm). This is an IMPORTANCE-OVERLAP proxy, NOT exact material worth (no scattering/GPT term).

Projection: each patch (i_tor,j_pol) adds +20cm shield over a coarse (theta_pol, phi_tor) sector, nfp-replicated.
We sample psi_dagger on the SHIELD surface over that angular sector and average (the mean coil-20 importance a
patch shields). APPROXIMATION (documented): circular-torus surface (R0, r_shield) -- NOT the true twisted VMEC
boundary -- so this is a first-order geometric projection, adequate for a rank/sign comparison, not a worth map.

Compares projected importance vs: coil-20 dR; global-peak change; inboard-cluster mean dR; #coils improved/worsened;
dR per added mass; TBR change. Reports Pearson & Spearman, sign accuracy, top-k overlap. Produces mechanism figure.
"""
import os, json, glob
import numpy as np
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__)); PW = os.path.dirname(HERE); D = os.path.join(PW, "cluster_pull")
R0_MAJOR = 1367.0        # major radius (cm), from patch-worth notes
R_SHIELD = 277.0         # ~shield minor radius (plasma a~217 + build ~60cm to shield mid)

# ---- adjoint importance ----
adj = np.load(f"{D}/adjoint_coil20_flat_P0.npz")
psi = adj["importance"].astype(float)            # (30,30,20)
ll = adj["lower_left"]; ur = adj["upper_right"]; dim = adj["dimension"]
def psi_at(x, y, z):
    ix = np.clip(((np.asarray(x)-ll[0])/(ur[0]-ll[0])*dim[0]).astype(int), 0, dim[0]-1)
    iy = np.clip(((np.asarray(y)-ll[1])/(ur[1]-ll[1])*dim[1]).astype(int), 0, dim[1]-1)
    iz = np.clip(((np.asarray(z)-ll[2])/(ur[2]-ll[2])*dim[2]).astype(int), 0, dim[2]-1)
    return psi[ix, iy, iz]

# ---- baseline + patches (coarse B_FD) ----
base = np.load(f"{D}/coil_step1b_uniform.npz")
Rb = base["coil_fast_flux"].astype(float); cells = base["coil_cells"].astype(int); tbr0 = float(base["tbr"])
pk = int(np.argmax(Rb)); peakcell = cells[pk]
inboard = np.argsort(Rb)[-4:]                    # the 4 hottest (inboard cluster)

rows = []
for p in sorted(glob.glob(f"{D}/coil_patch_t*p*.npz")):
    tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
    d = np.load(p); Rp = d["coil_fast_flux"].astype(float); tbr = float(d["tbr"])
    dd = np.load(f"{D}/delta_qh_patch_{tag}.npz")
    dl = dd["delta"]; tor_ang = dd["toroidal_angles"]; pol_ang = dd["poloidal_angles"]
    added = int(dd["added_cells"]) * float(dd["delta_cm"])
    # patch angular footprint (deg) from the nonzero delta grid cells
    ii, jj = np.where(dl > 0)
    th = np.radians(pol_ang[jj])                 # poloidal
    ph = np.radians(tor_ang[ii])                 # toroidal within one period (coil-20's own period)
    xr = (R0_MAJOR + R_SHIELD*np.cos(th)); x = xr*np.cos(ph); y = xr*np.sin(ph); z = R_SHIELD*np.sin(th)
    imp = float(np.mean(psi_at(x, y, z)))        # mean coil-20 importance over the patch shield sector
    E = Rp - Rb
    rows.append(dict(tag=tag, imp=imp, added=added,
                     dR20=E[pk], dpeak=(Rp.max()-Rb[pk]),
                     dcluster=float(np.mean(E[inboard])),
                     n_improved=int(np.sum(E < 0)), n_worsened=int(np.sum(E > 0)),
                     dR20_per_mass=E[pk]/added, dTBR=(tbr-tbr0)))

tags = [r["tag"] for r in rows]
imp = np.array([r["imp"] for r in rows])
metrics = {k: np.array([r[k] for r in rows]) for k in
           ["dR20", "dpeak", "dcluster", "n_improved", "dR20_per_mass", "dTBR"]}

print(f"adjoint importance projected onto {len(rows)} patch sectors (circular-torus approx).")
print(f"{'metric':14} {'Pearson':>9} {'Spearman':>9} {'sign_acc':>9}  interpretation")
report = {}
for k, v in metrics.items():
    # worth = -dR (reduction). Does HIGH importance predict LARGE reduction (negative dR)?
    pr = pearsonr(imp, -v)[0] if k not in ("n_improved",) else pearsonr(imp, v)[0]
    sp = spearmanr(imp, -v)[0] if k not in ("n_improved",) else spearmanr(imp, v)[0]
    sacc = float(np.mean((imp > np.median(imp)) == ((-v if k != "n_improved" else v) > np.median(-v if k != "n_improved" else v))))
    report[k] = {"pearson": round(float(pr), 3), "spearman": round(float(sp), 3), "sign_acc": round(sacc, 2)}
    print(f"{k:14} {pr:>9.2f} {sp:>9.2f} {sacc:>9.2f}")

# top-k overlap: do the top-3 highest-importance patches match the top-3 best coil-20 reducers?
top_imp = set(np.argsort(imp)[-3:]); top_reduce = set(np.argsort(metrics["dR20"])[:3])  # most negative dR20
overlap = len(top_imp & top_reduce)
print(f"\ntop-3 importance vs top-3 coil-20 reducers overlap: {overlap}/3")
print(f"  highest-importance patches: {[tags[i] for i in sorted(top_imp)]}")
print(f"  best coil-20 reducers:      {[tags[i] for i in sorted(top_reduce)]}")

# outcome
sp20 = report["dR20"]["spearman"]; sppk = report["dpeak"]["spearman"]
if abs(sp20) < 0.4 and abs(sppk) < 0.4:
    outcome = "IMPORTANCE PREDICTS NEITHER local nor global worth (in this coarse projection) -> attribution != placement, but SEE caveat"
elif abs(sp20) >= 0.5 and abs(sppk) < 0.4:
    outcome = "IMPORTANCE predicts LOCAL coil-20 worth but NOT global peak -> local-vs-system distinction"
elif abs(sp20) >= 0.5 and abs(sppk) >= 0.5:
    outcome = "IMPORTANCE predicts BOTH -> usable as a refinement guide"
else:
    outcome = "MIXED / ambiguous"

json.dump({"metrics_vs_importance": report, "top3_overlap": overlap, "outcome": outcome,
           "projection": "circular-torus approx (NOT true VMEC surface); first-order rank comparison only",
           "importance_quantity": "coil-20 adjoint importance psi_dagger (flat response), NOT contributon or exact worth"},
          open(f"{HERE}/rung2_attribution_vs_worth.json", "w"), indent=2)

# figure
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].scatter(imp, -metrics["dR20"], c="#c0392b");
    for i, t in enumerate(tags): ax[0].annotate(t, (imp[i], -metrics["dR20"][i]), fontsize=6)
    ax[0].set_xlabel("mean coil-20 adjoint importance over patch"); ax[0].set_ylabel("coil-20 reduction  -dR20")
    ax[0].set_title(f"local: Spearman {sp20:+.2f}")
    ax[1].scatter(imp, -metrics["dpeak"], c="#27ae60")
    for i, t in enumerate(tags): ax[1].annotate(t, (imp[i], -metrics["dpeak"][i]), fontsize=6)
    ax[1].set_xlabel("mean coil-20 adjoint importance over patch"); ax[1].set_ylabel("global peak reduction")
    ax[1].set_title(f"global: Spearman {sppk:+.2f}")
    fig.suptitle("Rung 2: coil-20 adjoint importance vs empirical intervention worth (coarse patches)")
    fig.tight_layout(); fig.savefig(f"{HERE}/fig_rung2_attribution_vs_worth.png", dpi=130); plt.close(fig)
    figmsg = "fig_rung2_attribution_vs_worth.png"
except Exception as e:
    figmsg = f"(fig skipped {e!r})"
print(f"\nOUTCOME: {outcome}\nwrote rung2_attribution_vs_worth.json, {figmsg}")
