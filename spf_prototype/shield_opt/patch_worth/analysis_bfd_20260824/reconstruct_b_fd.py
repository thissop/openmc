#!/usr/bin/env python
"""Reconstruct the full finite-difference response matrix from the existing patch study, per the
Execution Charter (2026-08-24). Audits data, recovers full per-coil response vectors, builds the
intervention-effect matrix E and (denominator-normalized) B_FD, analyzes structure / peak migration /
single-intervention controllability, and classifies the outcome.

Baseline caveats (documented, not hidden):
 - patch builds: shield 40 cm base, +20 cm in one (i_tor,j_pol) sector (breeder 50->30). added_cells,
   delta_cm=20 stored -> real denominators.
 - coil-run baseline = coil_step1b_uniform.npz (shield 42.55, delta_uniform~2.55) -> ~2.55 cm uniform
   confound in absolute dR; patches used FW-CADIS WW, baseline did NOT -> correlated-seed cancellation
   likely broken -> uncertainties treated as INDEPENDENT (sigma(dR)=sqrt(sp^2+sb^2)); this is an upper
   bound. Patch-to-patch STRUCTURE is robust to the common baseline offset; absolute reductions are not.
 - response = cell-INTEGRATED coil_fast_flux (track-length x cell; NOT volume-normalized density).
 - patches are SINGLE, mostly-disjoint sectors built one-at-a-time -> they are NOT an additive basis;
   arbitrary linear combinations are NOT geometrically validated (any combined design needs a new run).

Usage: python reconstruct_b_fd.py [--patch-dir DIR] [--baseline NPZ] [--out DIR]
"""
import os, sys, json, glob, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PW = os.path.dirname(HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patch-dir", default=os.path.join(PW, "cluster_pull"))
    ap.add_argument("--baseline", default=None, help="baseline coil npz (default: coil_step1b_uniform.npz in patch-dir)")
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--sig", type=float, default=2.0, help="significance threshold in sigma")
    a = ap.parse_args()
    D = a.patch_dir; OUT = a.out; os.makedirs(OUT, exist_ok=True)
    base_path = a.baseline or os.path.join(D, "coil_step1b_uniform.npz")

    # ---------- Step 0: audit ----------
    base = np.load(base_path)
    R0 = base["coil_fast_flux"].astype(float); E0 = base["coil_fast_flux_err"].astype(float)
    H0 = base["coil_heating"].astype(float); HE0 = base["coil_heating_err"].astype(float)
    cells = base["coil_cells"].astype(int); tbr0 = float(base["tbr"])
    NC = len(R0)
    patch_files = sorted(glob.glob(os.path.join(D, "coil_patch_t*p*.npz")))
    manifest = {"baseline": os.path.basename(base_path), "baseline_shield_cm": 42.55,
                "baseline_coil_cells": cells.tolist(), "baseline_tbr": tbr0,
                "response_metric": "coil_fast_flux (cell-integrated, NOT volume-normalized)",
                "n_coils": NC, "patches": [], "caveats": [
                    "patch shield-outside=40cm vs baseline 42.55cm -> ~2.55cm uniform confound in absolute dR",
                    "patches used FW-CADIS WW; baseline did not -> correlated seed likely broken -> independent-sigma bound",
                    "single disjoint-sector patches -> NOT an additive basis; combined designs need a new run",
                ]}
    names, E, SIG, EH, SIGH = [], [], [], [], []
    tbr, added, dcm, i_tor, j_pol = [], [], [], [], []
    for p in patch_files:
        tag = os.path.basename(p).split("coil_patch_")[1].split(".npz")[0]
        d = np.load(p); Rp = d["coil_fast_flux"].astype(float); Ep = d["coil_fast_flux_err"].astype(float)
        Hp = d["coil_heating"].astype(float); HEp = d["coil_heating_err"].astype(float)
        dd = np.load(os.path.join(D, f"delta_qh_patch_{tag}.npz"))
        names.append(tag)
        E.append(Rp - R0); SIG.append(np.sqrt(Ep**2 + E0**2))
        EH.append(Hp - H0); SIGH.append(np.sqrt(HEp**2 + HE0**2))
        tbr.append(float(d["tbr"])); added.append(int(dd["added_cells"])); dcm.append(float(dd["delta_cm"]))
        i_tor.append(int(dd["i_tor"])); j_pol.append(int(dd["j_pol"]))
        manifest["patches"].append({"tag": tag, "i_tor": int(dd["i_tor"]), "j_pol": int(dd["j_pol"]),
                                    "added_cells": int(dd["added_cells"]), "delta_cm": float(dd["delta_cm"]),
                                    "tbr": float(d["tbr"]), "peak_coil_flux": float(Rp.max()),
                                    "peak_cell": int(cells[np.argmax(Rp)])})
    E = np.array(E); SIG = np.array(SIG); EH = np.array(EH); SIGH = np.array(SIGH)
    tbr = np.array(tbr); added = np.array(added); dcm = np.array(dcm)
    NP = len(names); pk0 = int(np.argmax(R0))
    denom = added * dcm                      # added shield-cell-cm (per-patch, unequal)
    B = -E / denom[:, None]                  # response reduction per unit added shield (charter sign)
    frac = E / np.maximum(R0, 1e-30)[None, :]

    json.dump(manifest, open(os.path.join(OUT, "patch_manifest.json"), "w"), indent=2)
    # CSVs
    np.savetxt(os.path.join(OUT, "baseline_responses.csv"),
               np.column_stack([cells, R0, E0, H0, HE0]), delimiter=",",
               header="coil_cell,fast_flux,fast_flux_err,heating,heating_err", comments="")
    with open(os.path.join(OUT, "patch_response_vectors.csv"), "w") as f:
        f.write("patch,i_tor,j_pol,added_cells,delta_cm,tbr," + ",".join(f"dR_c{c}" for c in cells) + "\n")
        for k in range(NP):
            f.write(f"{names[k]},{i_tor[k]},{j_pol[k]},{added[k]},{dcm[k]},{tbr[k]:.4f}," +
                    ",".join(f"{E[k,j]:.4e}" for j in range(NC)) + "\n")
    np.savetxt(os.path.join(OUT, "b_fd_raw.csv"), B, delimiter=",",
               header=",".join(f"c{c}" for c in cells), comments="")
    np.savetxt(os.path.join(OUT, "b_fd_fractional.csv"), frac, delimiter=",",
               header=",".join(f"c{c}" for c in cells), comments="")

    # ---------- Step 1-2: significance, peak migration ----------
    sig_frac = float(np.mean(np.abs(E) > a.sig * SIG))
    R = R0[None, :] + E                       # per-patch full new response
    new_peak = R.max(axis=1); new_peak_cell = cells[np.argmax(R, axis=1)]
    peak_reduced = new_peak < R0[pk0]
    # tradeoff detector: patch lowers the (old) peak coil but raises >=1 other coil significantly
    lowers_peak = (E[:, pk0] < -a.sig * SIG[:, pk0])
    raises_other = np.array([np.any((E[k] > a.sig * SIG[k]) & (np.arange(NC) != pk0)) for k in range(NP)])

    # ---------- Step 3: SVD of uncertainty-whitened E ----------
    Ew = E / np.maximum(SIG, 1e-30)
    U, S, Vt = np.linalg.svd(Ew, full_matrices=False)
    p = S / S.sum(); reff = float(np.exp(-np.sum(p * np.log(p + 1e-30))))
    cum = np.cumsum(S**2) / np.sum(S**2)
    # is the top mode a COMMON reduction (all-same-sign response) or a TRADEOFF (mixed sign)?
    v1 = Vt[0]; common_mode = float(np.mean(np.sign(v1) == np.sign(np.median(v1))))

    # ---------- Step 4: single-intervention control authority (NO superposition) ----------
    best = int(np.argmin(new_peak))
    best_peak_ratio = float(R0[pk0] / new_peak[best])

    # ---------- outcome classification ----------
    n_sig_reduce_peak = int(np.sum(lowers_peak & ~raises_other & peak_reduced))
    n_tradeoff = int(np.sum(lowers_peak & raises_other))
    if sig_frac < 0.2:
        outcome = "4. EFFECTS UNRESOLVED (most dR within noise)"
    elif n_sig_reduce_peak >= 1 and best_peak_ratio > 1.3:
        outcome = "1. HIDDEN STRUCTURED CONTROL: single sector patches significantly reduce the peak (placement-dependent)"
    elif common_mode > 0.8 and cum[0] > 0.6:
        outcome = "2. COMMON/GLOBAL MODE dominates (bulk attenuation reduces many coils together)"
    elif n_tradeoff >= max(1, NP // 2):
        outcome = "3. TRADEOFF-ONLY: patches move load among coils"
    else:
        outcome = "1/2 MIXED (structured + common); see per-patch table"

    # ---------- results.md ----------
    with open(os.path.join(OUT, "results.md"), "w") as f:
        f.write("# B_FD reconstruction — results (2026-08-24)\n\n")
        f.write(f"**Baseline:** `{os.path.basename(base_path)}` (shield 42.55 cm). "
                f"**Response:** cell-integrated coil_fast_flux (NOT density). **Patches:** {NP} single +20 cm sectors.\n\n")
        f.write("## Caveats (charter §0.5, §Statistical)\n")
        for c in manifest["caveats"]: f.write(f"- {c}\n")
        f.write(f"\n## Significance\n- {100*sig_frac:.0f}% of dR entries exceed {a.sig}sigma (independent-error bound).\n")
        f.write(f"\n## Peak migration & control (per patch)\n\n")
        f.write(f"baseline peak: coil {cells[pk0]} = {R0[pk0]:.4e}, TBR {tbr0:.3f}\n\n")
        f.write(f"| patch | (i,j) | added·cm | dR_peak | new_peak | peakX | ->cell | dTBR% | class |\n|---|---|---|---|---|---|---|---|---|\n")
        for k in np.argsort(new_peak):
            cls = "reduce" if (lowers_peak[k] and not raises_other[k]) else ("tradeoff" if (lowers_peak[k] and raises_other[k]) else ("backfire" if E[k,pk0]>a.sig*SIG[k,pk0] else "weak"))
            f.write(f"| {names[k]} | ({i_tor[k]},{j_pol[k]}) | {denom[k]:.0f} | {E[k,pk0]:+.2e} | {new_peak[k]:.3e} | "
                    f"{R0[pk0]/new_peak[k]:.1f}x | {new_peak_cell[k]} | {100*(tbr[k]/tbr0-1):+.1f} | {cls} |\n")
        f.write(f"\n## SVD of uncertainty-whitened E ({NP}x{NC})\n")
        f.write(f"- singular values: {np.round(S,1).tolist()}\n- effective rank (entropy): {reff:.2f}\n")
        f.write(f"- cumulative energy: 1SV={cum[0]:.2f}, 2={cum[1]:.2f}, 3={cum[2]:.2f}\n")
        f.write(f"- top-mode common-sign fraction: {common_mode:.2f} (>0.8 => common reduction mode; ~0.5 => tradeoff)\n")
        f.write(f"\n## Single-intervention control authority (superposition NOT assumed)\n")
        f.write(f"- best single patch: **{names[best]}**  peak {R0[pk0]:.4e} -> {new_peak[best]:.4e} "
                f"(**{best_peak_ratio:.1f}x**, migrates to cell {new_peak_cell[best]}), dTBR {100*(tbr[best]/tbr0-1):+.1f}%\n")
        f.write(f"- patches that reduce the peak (no significant off-peak rise): {n_sig_reduce_peak}\n")
        f.write(f"- tradeoff patches (lower peak, raise another): {n_tradeoff}\n")
        f.write(f"\n## OUTCOME\n**{outcome}**\n\n")
        f.write("Note: the absolute peakX is baseline-confounded (~2.55 cm + WW-vs-noWW); the placement-dependent "
                "STRUCTURE (which sectors reduce vs backfire) is the robust result. Combined designs and clean "
                "derivatives require Gate A: a fresh, consistently-baselined OpenMC run.\n")

    # ---------- figures ----------
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        order = np.argsort(new_peak)
        # 1) fractional dR heatmap
        fig, ax = plt.subplots(figsize=(10, 5))
        im = ax.imshow(frac[order]*100, aspect="auto", cmap="RdBu_r", vmin=-100, vmax=100)
        ax.set_yticks(range(NP)); ax.set_yticklabels([names[k] for k in order], fontsize=7)
        ax.set_xticks(range(NC)); ax.set_xticklabels([f"c{c}" for c in cells], fontsize=6, rotation=90)
        ax.set_title("Fractional coil-response change per patch (%)  [blue=reduced]")
        fig.colorbar(im, label="% change vs baseline"); fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_frac_heatmap.png"), dpi=130); plt.close(fig)
        # 2) peak migration
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.axhline(R0[pk0], color="k", ls="--", label=f"baseline peak (coil {cells[pk0]})")
        ax.bar(range(NP), new_peak[order], color=["#27ae60" if new_peak[order][m]<R0[pk0] else "#c0392b" for m in range(NP)])
        ax.set_xticks(range(NP)); ax.set_xticklabels([names[k] for k in order], rotation=90, fontsize=7)
        ax.set_ylabel("new peak coil flux"); ax.set_title("Peak coil load after each single patch"); ax.legend()
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_peak_migration.png"), dpi=130); plt.close(fig)
        # 3) singular spectrum
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.semilogy(range(1, len(S)+1), S, "o-"); ax.set_xlabel("mode"); ax.set_ylabel("singular value (whitened)")
        ax.set_title(f"E singular spectrum (eff. rank {reff:.2f})"); fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_singular_spectrum.png"), dpi=130); plt.close(fig)
        figs = "fig_frac_heatmap.png, fig_peak_migration.png, fig_singular_spectrum.png"
    except Exception as e:
        figs = f"(figures skipped: {e!r})"

    print(f"OUTCOME: {outcome}")
    print(f"sig_frac={sig_frac:.2f} eff_rank={reff:.2f} cum1={cum[0]:.2f} common_mode={common_mode:.2f}")
    print(f"best single patch {names[best]}: {best_peak_ratio:.1f}x (baseline-confounded)")
    print(f"wrote: results.md, patch_manifest.json, *.csv, {figs}  -> {OUT}")


if __name__ == "__main__":
    main()
