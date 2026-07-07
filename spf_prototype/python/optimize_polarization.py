#!/usr/bin/env python
"""Per-device polarization optimum for MINIMUM first-wall peaking, from the real-neutronics (θ,φ) NWL
maps (realdev_wallmap.py). NWL is exactly linear in the emission quadrupole a₂ (kernel 1+a₂P₂), so
NWL(a₂)=NWL_unpol + a₂(NWL_par-NWL_unpol) is built from the three TRANSPORT maps and the linearity is
verified in the MC data. Scans a₂∈[-1,+1] (physically achievable: -1 pure perpendicular / A mode,
+1 pure parallel / B-C mode) for min peaking factor PF=max/⟨·⟩. Emits grayscale (θ,φ) maps of
NWL/⟨NWL⟩ (unpol | optimal) per device and a markdown summary of the optimum.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
try:
    import smplotlib  # noqa: F401
except Exception:
    pass
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
D = HERE.parent / "data"
FIGD = HERE.parent / "figs" / "wallmaps"; FIGD.mkdir(parents=True, exist_ok=True)
DOCS = HERE.parent / "docs" / "notes"
IDS = [803097, 886079, 932746, 59509, 1960314]
NFP = {803097: 3, 886079: 2, 932746: 3, 59509: 3, 1960314: 5}


def pf(nwl, dA2):
    mean = np.sum(nwl * dA2) / np.sum(dA2)
    return float(nwl.max() / mean), float(mean)


def a2_to_abc(a2):
    """A representative achievable (a,b,c): mix pure-A(perp) or pure-B(par) with unpol so the emission
    quadrupole equals a2.  a2(w_perp,w_par) = (w_par/2 - 2 w_perp/3)/(w_par/2 + 2 w_perp/3)."""
    ts = np.linspace(0, 1, 2001)
    best, babc = 9, (1/3, 1/3, 1/3)
    for t in ts:
        if a2 <= 0:                       # blend unpol->pure A (a=1)
            a, b, c = (1 - t) / 3 + t, (1 - t) / 3, (1 - t) / 3
        else:                             # blend unpol->pure B (b=1)
            a, b, c = (1 - t) / 3, (1 - t) / 3 + t, (1 - t) / 3
        wp, wq = 0.75 * a, (2/3) * b + (1/3) * c
        av = (0.5 * wq - (2/3) * wp) / (0.5 * wq + (2/3) * wp)
        if abs(av - a2) < best:
            best, babc = abs(av - a2), (a, b, c)
    return babc


def main():
    rows = []
    a2s = np.linspace(-1, 1, 801)
    for ID in IDS:
        wm = D / f"quasr{ID}_wallmap.npz"
        if not wm.exists():
            print(f"skip {ID}: no wallmap"); continue
        f = dict(np.load(wm))
        u, p, q = f["unpol"], f["perp"], f["par"]
        dA2 = np.broadcast_to(f["dA"][:, None], u.shape)
        lin = np.abs(p + q - 2 * u).mean() / u.mean()
        PFs = [pf(u + a2 * (q - u), dA2)[0] for a2 in a2s]
        j = int(np.argmin(PFs))
        a2star = float(a2s[j]); pf_opt = PFs[j]
        pf_unpol = pf(u, dA2)[0]
        abc = a2_to_abc(a2star)
        rows.append(dict(ID=ID, nfp=NFP[ID], a2=a2star, pf_unpol=pf_unpol, pf_opt=pf_opt,
                         red=100 * (pf_unpol - pf_opt) / pf_unpol, abc=abc, lin=lin,
                         u=u, opt=u + a2star * (q - u), dA2=dA2))
        print(f"{ID}: lin={lin:.3f}  a2*={a2star:+.2f}  PF {pf_unpol:.2f}->{pf_opt:.2f} "
              f"({100*(pf_unpol-pf_opt)/pf_unpol:+.1f}%)")

    # grayscale (θ,φ) maps: unpol | optimal, shared scale
    for r in rows:
        um = r["u"] / (np.sum(r["u"] * r["dA2"]) / np.sum(r["dA2"]))
        om = r["opt"] / (np.sum(r["opt"] * r["dA2"]) / np.sum(r["dA2"]))
        vmax = max(um.max(), om.max())
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ext = [0, 360, 0, 360]
        for a, M, t in [(ax[0], um, "unpolarized"),
                        (ax[1], om, f"optimal  (a$_2$={r['a2']:+.2f})")]:
            im = a.imshow(M, origin="lower", extent=ext, aspect="auto", cmap="gray_r",
                          vmin=0, vmax=vmax)
            a.set_xlabel(r"toroidal $\phi$ [deg]"); a.set_ylabel(r"poloidal $\theta$ [deg]")
            a.set_title(f"{t}    PF={M.max():.2f}")
            fig.colorbar(im, ax=a, label="NWL / mean")
        fig.suptitle(f"Device {r['ID']} (nfp{r['nfp']}) first-wall NWL, transport (StellaratorSource).  "
                     f"Peaking {r['pf_unpol']:.2f} to {r['pf_opt']:.2f}", y=1.03)
        fig.tight_layout()
        fig.savefig(FIGD / f"wallmap_{r['ID']}.png", dpi=170, bbox_inches="tight")
        plt.close(fig)

    # markdown
    lines = ["# Optimal SPF polarization for minimum first-wall peaking",
             "",
             "Real OpenMC transport of the native `StellaratorSource` on each device's **VMEC** "
             "equilibrium (real √g births + real b̂), free-streaming to a torus first wall; every "
             "neutron's wall crossing captured via `surface_source_write` and histogrammed in "
             "(poloidal θ, toroidal φ). NWL is exactly linear in the emission quadrupole "
             "`a₂` (kernel `1+a₂P₂(cosθ_B)`), so `NWL(a₂)=NWL_unpol+a₂(NWL_par−NWL_unpol)` is built "
             "from the three transport maps; `a₂∈[−1,+1]` (−1 = pure perpendicular / A mode, "
             "+1 = pure parallel / B–C mode, 0 = unpolarized). Optimum = min peaking factor "
             "`PF = max(NWL)/⟨NWL⟩`.",
             "",
             "| device | nfp | optimal a₂ | achievable (a,b,c) | PF unpol | PF optimal | reduction | MC linearity |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        a, b, c = r["abc"]
        lean = "perpendicular" if r["a2"] < -0.05 else ("parallel" if r["a2"] > 0.05 else "≈unpolarized")
        lines.append(f"| {r['ID']} | {r['nfp']} | **{r['a2']:+.2f}** ({lean}) | "
                     f"({a:.2f}, {b:.2f}, {c:.2f}) | {r['pf_unpol']:.2f} | {r['pf_opt']:.2f} | "
                     f"**{r['red']:+.1f}%** | {r['lin']:.3f} |")
    lines += ["",
              "**MC linearity** = mean|NWL_perp+NWL_par−2·NWL_unpol| / ⟨NWL_unpol⟩ (→0 confirms the "
              "load is linear in a₂; residual is Monte-Carlo noise).",
              "",
              "Grayscale (θ,φ) NWL/⟨NWL⟩ maps (unpolarized | optimal) per device: "
              "`figs/wallmaps/wallmap_<ID>.png`."]
    (DOCS / "OPTIMAL_POLARIZATION.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {len(rows)} maps -> figs/wallmaps/ + docs/notes/OPTIMAL_POLARIZATION.md")


if __name__ == "__main__":
    main()
