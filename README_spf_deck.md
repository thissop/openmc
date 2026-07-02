# SPF conformal-stellarator deck — figures, slides, script, Q&A

Everything here is built **only** from the plain-data bundle in `./spf_export/`
(OpenMC Monte-Carlo transport of spin-polarized-fuel DT neutrons through a
conformal breeding blanket on the `precise_QA` stellarator, nfp = 2, 10,000
histories/config). No OpenMC, DAGMC, or cluster is needed — the scripts read
numpy/pandas/json only.

## Deliverables
| file | what |
|---|---|
| `figures.py` | regenerates every figure into `figs/` (200 dpi PNG) + `figs/derived_numbers.json` |
| `figs/fig1_lcfs_wireframe.png` | 3D plasma boundary (LCFS), nfp = 2 QA shape |
| `figs/fig2_cross_section.png` | poloidal blanket cross-sections at φ = 0 and 90° (non-axisymmetry) |
| `figs/fig3_nested_shells.png` | 3D translucent nested shells (plasma / FLiBe / coil), 300° cutaway |
| `figs/fig4_radial_schematic.png` | to-scale radial layer stack + inset zoom of the thin inner layers |
| `figs/fig5a_tbr.png` | tritium breeding ratio vs mode (scatter) |
| `figs/fig5b_coil_flux.png` | coil fast flux (>0.1 MeV) vs mode, free vs scatter |
| `figs/fig5c_heating.png` | per-layer nuclear heating (scatter), log-y, by mode |
| `figs/fig5d_eta.png` | directional signature + hero η = 0.63 ± 0.22 |
| `slides/talk.tex` → `slides/talk.pdf` | Beamer deck (default theme, 12 frames) |
| `talk_script.md` | first-person speaker notes, one block per slide |
| `qa.md` | 11 anticipated expert questions with honest answers |

## Regenerate everything
```bash
# 1. figures (writes figs/*.png and figs/derived_numbers.json)
cd /Users/tkiker/Documents/GitHub/openmc
python3 figures.py

# 2. compile the deck (run pdflatex twice for the nav bar; needs stock TeX Live)
cd slides
pdflatex -interaction=nonstopmode talk.tex
pdflatex -interaction=nonstopmode talk.tex
# -> slides/talk.pdf  (12 pages)
```
Requirements: Python 3 with numpy, pandas, matplotlib, scipy; a TeX Live
`pdflatex` with `graphicx` + `booktabs`. Figures are deterministic (no RNG).

## Headline numbers (all recomputed from the bundle; see `figs/derived_numbers.json`)
- **TBR** (tritium/source-neutron, FLiBe, scattering): unpolarized 0.680 ± 0.008,
  perpendicular 0.687 ± 0.007, parallel 0.693 ± 0.015. Free-streaming ≈ 0 by
  construction. Mode spread is within ~1σ — **polarization does not measurably
  change breeding here.**
- **Coil fast flux** (>0.1 MeV, /cm²/src): free 31.68 / 26.84 / 35.96,
  scatter 11.45 / 10.16 / 12.43 (unpol / perp / par).
- **Directional signature** (coil flux): parallel vs unpolarized +13.5% free →
  +8.5% scatter; perpendicular vs unpolarized −15.3% → −11.3%; parallel–perpendicular
  swing 34% free / 22% scatter.
- **η ≡ (scatter signal)/(free signal) for parallel = 8.5/13.5 = 0.63 ± 0.22** — the
  fraction of the ideal free-streaming steering that survives the conformal wall +
  scattering.
- **Coil heating** (scatter, eV/src): parallel 457,856 > unpolarized 414,510 >
  perpendicular 371,994 — same directional sign as the coil flux (independent
  corroboration).

## Honest framing (enforced across all deliverables)
- **Sign:** parallel (B/C, 1 + 3cos²θ) **raises** the coil load; perpendicular
  (A, sin²θ) **lowers** it. A magnet-protection benefit in this geometry would come
  from the *perpendicular* mode. (The bundle's `RESULTS_tier8_conformal.md` has a
  "parallel reduces coil flux" typo — it is wrong; we do not repeat it.)
- **Smoke run:** 10k histories; coil tallies sit behind 116 cm of blanket+shield+coil,
  so MC error is large. η = 0.63 ± 0.22 is **indicative** — needs variance reduction
  (weight windows / FW-CADIS) to converge. The robust result is the *trend* (a
  coherent directional ordering in both coil flux and coil heating, in both streams).
- **Machinery validation on a public equilibrium (`precise_QA`), not a reactor claim.**
- **Geometry caveat:** the wall is a naive constant-normal-offset ("conformal")
  construction. At φ = 0 the thick outer shells (shield, coil) develop a small
  self-fold at the inboard midplane where the offset overshoots the machine axis
  (coil_outer reaches R ≈ −2.4 cm). Enclosed areas stay monotonically nested and it
  is invisible at device scale, but a physics-grade run needs a proper
  non-self-intersecting offset. (Checked via poloidal-contour segment-crossing +
  area-monotonicity tests at every exported φ.)

## Notes on the bundle
- `summary.json` in the bundle is empty (`{}`); all metrics are recomputed from
  `results/tallies_long.csv` and cross-checked against
  `results/RESULTS_tier8_conformal.md` (they agree to the printed precision).
- `tallies_long.csv` is machine-written and not RFC-clean (unescaped commas in the
  mesh filter tuples, the `(n,Xt)` reaction name, and the coil energy filter);
  `figures.py` parses the well-formed physics rows tally-by-tally against the known
  layout. The 1.9M noisy `wall_current_phi` mesh rows are not used for these figures.
- Geometry is exported at the **base equilibrium scale** (LCFS R ≈ 63–133 cm),
  plotted in cm as-is; the transport run itself used device scale 15.
