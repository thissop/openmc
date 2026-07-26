#!/usr/bin/env python
"""Engineering-relevance (blanket-fit) filter + geometric SELECTION-EFFECT analysis.

The user's criterion for the 200-300 device magnet-side study: keep only devices where a
REAL radial build (~1 m+ blanket+shield) physically FITS between plasma and coils -- exclude
the tiny (~16 cm) plasma-coil-standoff devices -- AND quantify how this filter biases the
device population (the geometric selection effect), because that bias confounds any law we
derive on the filtered set.

KEY MODELING SUBTLETY (this IS the selection effect): "does a blanket fit" depends on the
reactor-scaling convention and on as-shipped vs min-achievable coil placement:
  (A) AS-SHIPPED: the QUASR coil solution's actual plasma-coil gap d_min, scaled to a reactor
      that fixes the minor radius to a_reactor (ARIES-CS a=1.704 m, Kappel convention).
      physical_gap = (d_min/a) * a_reactor.  Operative "is there room with THESE coils."
  (B) MIN-ACHIEVABLE: Kappel L_REGCOIL (the closest coils CAN sit) / the L_gradB proxy.
      The floor -- "could coils EVER be placed to leave blanket room." Strong nfp/class bias.
Both are reported; the gap between them is a real result (QUASR ships generous coils).

Radial build = the validated 8-layer corrected build (~1.29 m); also report a 1.0 m relaxed cut.

In : data/quasr/coil_standoff.csv (366 QA/QH, d_min/a etc.), data/kappel/kappel_configs.csv (45, L_REGCOIL/L_gradB, incl QI).
Out: shield_opt/data/engineering_relevant_devices.csv + a printed selection-effect report + fig.
"""
import csv
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
A_REACTOR = 1.704          # m, ARIES-CS minor radius (Kappel scaling)
RADIAL_BUILD = 1.29        # m, validated corrected 8-layer build (FW..thermal_shield)
RELAXED = 1.0              # m, relaxed cut


def load_standoff():
    rows = list(csv.DictReader(open(ROOT / "data/quasr/coil_standoff.csv")))
    for r in rows:
        r["nfp"] = int(float(r["nfp"]))
        r["aspect"] = float(r["aspect"])
        r["d_over_a"] = float(r["d_min_over_a"])
        r["gap_m"] = r["d_over_a"] * A_REACTOR          # (A) as-shipped physical gap
        r["qs_error"] = float(r["qs_error"])
        r["cls"] = r["symmetry_class"]
    return rows


def selection_effect(rows, thr):
    """How the blanket-fit cut (gap>=thr) shifts the population. Report kept vs dropped
    distributions of the geometry knobs -- the confound the user wants quantified."""
    g = np.array([r["gap_m"] for r in rows])
    keep = g >= thr
    def dist(mask, key):
        v = np.array([r[key] for r in rows])[mask]
        return v
    rep = {"threshold_m": thr, "n_total": len(rows), "n_keep": int(keep.sum()),
           "frac_keep": float(keep.mean())}
    # class + nfp composition shift
    for c in sorted(set(r["cls"] for r in rows)):
        m = np.array([r["cls"] == c for r in rows])
        rep[f"class_{c}_keep_frac"] = float((keep & m).sum() / max(m.sum(), 1))
    # geometry-knob shift kept vs dropped (the bias)
    for key in ("nfp", "aspect", "qs_error"):
        kv, dv = dist(keep, key), dist(~keep, key)
        rep[f"{key}_kept_med"] = float(np.median(kv)) if kv.size else float("nan")
        rep[f"{key}_dropped_med"] = float(np.median(dv)) if dv.size else float("nan")
    return keep, rep


def main():
    rows = load_standoff()
    print("=== ENGINEERING-RELEVANCE (blanket-fit) + selection effect ===")
    print(f"scaling: a_reactor={A_REACTOR} m (ARIES-CS); radial build={RADIAL_BUILD} m "
          f"(relaxed {RELAXED} m)\n")
    print("(A) AS-SHIPPED QUASR coils (d_min/a * a_reactor):")
    keep_ship = None
    for thr in (RELAXED, RADIAL_BUILD):
        keep, rep = selection_effect(rows, thr)
        if thr == RADIAL_BUILD:
            keep_ship = keep
        print(f"  gap>={thr} m: keep {rep['n_keep']}/{rep['n_total']} "
              f"({100*rep['frac_keep']:.0f}%) | "
              + " ".join(f"{c}={rep.get('class_'+c+'_keep_frac',0)*100:.0f}%"
                         for c in ("QA", "QH"))
              + f" | nfp med kept {rep['nfp_kept_med']:.0f} vs dropped {rep['nfp_dropped_med']:.0f}"
              + f" | aspect med kept {rep['aspect_kept_med']:.1f} vs dropped {rep['aspect_dropped_med']:.1f}")

    # (B) min-achievable view from Kappel (has QI + tokamak + stellarator)
    print("\n(B) MIN-ACHIEVABLE separation (Kappel L_REGCOIL, m; the floor incl QI):")
    kap = list(csv.DictReader(open(ROOT / "data/kappel/kappel_configs.csv")))
    L = np.array([float(r["L_REGCOIL_m"]) for r in kap])
    Lg = np.array([float(r["L_gradB_m"]) for r in kap])
    cls = np.array([r["qs_type"] for r in kap])
    nfp = np.array([int(r["nfp"]) for r in kap])
    for thr in (RELAXED, RADIAL_BUILD):
        ok = L >= thr
        print(f"  L_REGCOIL>={thr} m: {ok.sum()}/{len(kap)} ({100*ok.mean():.0f}%) | "
              + " ".join(f"{c}={int((ok&(cls==c)).sum())}/{int((cls==c).sum())}"
                         for c in sorted(set(cls))))
    print(f"  SELECTION EFFECT (min-achievable): corr(nfp, L_REGCOIL)={np.corrcoef(nfp,L)[0,1]:+.3f}, "
          f"corr(nfp, L_gradB)={np.corrcoef(nfp,Lg)[0,1]:+.3f}  -> small nfp => more blanket room "
          f"(Kappel's reactor-attractive trend). Blanket-fit filter PREFERENTIALLY KEEPS low-nfp/QI/QA.")

    # emit the engineering-relevant device list (as-shipped criterion, the operative one for QUASR)
    out = ROOT / "shield_opt/data/engineering_relevant_devices.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "class", "nfp", "aspect", "d_min_over_a", "gap_m_ARIESCS",
                    "blanket_fit_1.29m", "blanket_fit_1.0m"])
        for r, k in zip(rows, keep_ship):
            w.writerow([r["ID"], r["cls"], r["nfp"], f"{r['aspect']:.3f}",
                        f"{r['d_over_a']:.4f}", f"{r['gap_m']:.3f}",
                        int(r["gap_m"] >= RADIAL_BUILD), int(r["gap_m"] >= RELAXED)])
    ne = sum(1 for r in rows if r["gap_m"] >= RADIAL_BUILD)
    print(f"\nWROTE {out}  ({ne} blanket-fit @ {RADIAL_BUILD} m of {len(rows)} QA/QH)")
    print("NOTE: QUASR-366 is QA/QH only. QI blanket-fit devices come from Kappel-45 (median "
          "L_REGCOIL 3.58 m -> QI fits easily) + ConStellaration; the scaled magnet-side set "
          "(200-300, QI/QA/QH) needs d_min computed for MORE devices incl QI coils (agent task).")

    # figure: gap distribution + selection effect
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        g = np.array([r["gap_m"] for r in rows]); C = np.array([r["cls"] for r in rows])
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for c, col in (("QA", "C0"), ("QH", "C1")):
            ax[0].hist(g[C == c], bins=np.linspace(0, 12, 25), alpha=0.6, label=c, color=col)
        ax[0].axvline(RADIAL_BUILD, color="k", ls="--", label=f"{RADIAL_BUILD} m build")
        ax[0].set_xlabel("Plasma-Coil Gap at ARIES-CS Scale (m)")
        ax[0].set_ylabel("Devices"); ax[0].set_title("As-Shipped Blanket-Fit"); ax[0].legend()
        ax[1].scatter(nfp, L, c=["C2" if x == "QI" else "C0" if x == "QA" else "C1" if x == "QH" else "C7" for x in cls], s=28)
        ax[1].axhline(RADIAL_BUILD, color="k", ls="--")
        ax[1].set_xlabel("Field Periods nfp"); ax[1].set_ylabel("Min-Achievable Separation L_REGCOIL (m)")
        ax[1].set_title("Selection Effect: Small nfp -> More Room (Kappel)")
        fig.tight_layout()
        fp = ROOT / "figs/engineering_relevance.png"
        fp.parent.mkdir(exist_ok=True); fig.savefig(fp, dpi=110)
        print(f"WROTE {fp}")
    except Exception as e:
        print("fig skipped:", e)


if __name__ == "__main__":
    main()
