#!/usr/bin/env python
"""Multi-device SPF steering comparison (smplotlib): native StellaratorSource in transport on each
converged VMEC equilibrium. Left: midplane perpendicular steering (inboard +, outboard -) per device,
against the idealized axisymmetric capstone reference (+40.6% / -21.2%). Right: the inboard-perp
poloidal steering profile per device (Z normalized). Reads data/quasr<ID>_wallcurrent.npz + fluxmaps.
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
IDS = [803097, 886079, 932746, 59509, 1960314]
CAP_IN, CAP_OUT = 40.6, -21.2   # idealized axisymmetric capstone oracle


def main():
    recs = []
    for ID in IDS:
        wf = D / f"quasr{ID}_wallcurrent.npz"; ff = D / f"quasr{ID}_vmec_fluxmap.npz"
        if not wf.exists():
            print(f"skip {ID}: no wallcurrent"); continue
        w = dict(np.load(wf)); f = dict(np.load(ff))
        zc = w["zc"]; mid = int(np.argmin(np.abs(zc)))
        ip = 100 * (w["perp_inb"][mid] - w["unpol_inb"][mid]) / w["unpol_inb"][mid]
        op = 100 * (w["perp_out"][mid] - w["unpol_out"][mid]) / w["unpol_out"][mid]
        R = f["R"]; nfp = int(f["nfp"])
        R0 = 0.5 * (R.max() + R.min()); asp = R0 / (0.5 * (R.max() - R.min()))
        recs.append(dict(ID=ID, nfp=nfp, asp=asp, ip=ip, op=op, w=w, zc=zc))
    recs.sort(key=lambda r: r["asp"])
    print("device   nfp  aspect  inboard-perp  outboard-perp")
    for r in recs:
        print(f"{r['ID']:>8} {r['nfp']:>3} {r['asp']:>6.1f}   {r['ip']:+6.1f}%     {r['op']:+6.1f}%")

    fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.2))
    x = np.arange(len(recs))
    ax[0].bar(x - 0.2, [r["ip"] for r in recs], 0.4, color="#2a6fdb", label="inboard (perp)")
    ax[0].bar(x + 0.2, [r["op"] for r in recs], 0.4, color="#9ec3ef", label="outboard (perp)")
    ax[0].axhline(CAP_IN, ls=":", color="#c0392b", lw=1.2, label="ideal axisym +40.6%")
    ax[0].axhline(CAP_OUT, ls="--", color="#c0392b", lw=1.2, label="ideal axisym -21.2%")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xticks(x); ax[0].set_xticklabels([f"{r['ID']}\n nfp{r['nfp']}, A={r['asp']:.1f}" for r in recs], fontsize=8)
    ax[0].set_ylabel("midplane perp steering  [%]")
    ax[0].set_title("perp steering per device (sorted by aspect)")
    ax[0].legend(frameon=False, fontsize=8)

    for r in recs:
        u = r["w"]["unpol_inb"]
        frac = 100 * (r["w"]["perp_inb"] - u) / np.where(u > 0, u, np.nan)
        ax[1].plot(frac, r["zc"] / np.abs(r["zc"]).max(),
                   lw=1.6, label=f"{r['ID']} (nfp{r['nfp']}, A{r['asp']:.1f})")
    ax[1].axvline(0, color="grey", ls=":", lw=0.8)
    ax[1].set_xlabel("inboard perp steering  [%]"); ax[1].set_ylabel("Z / Z$_{max}$")
    ax[1].set_title("inboard steering profile per device"); ax[1].legend(frameon=False, fontsize=8)

    fig.suptitle("Multi-device SPF steering: native StellaratorSource on 5 real VMEC equilibria "
                 "(transport)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGD / "multidevice_steering.png", dpi=175, bbox_inches="tight")
    print("plot -> figs/stellarator_source/multidevice_steering.png")


if __name__ == "__main__":
    main()
