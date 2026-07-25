# RESULTS_adjoint — frozen provenance for the §4 adjoint magnet-shielding results

*Citable ledger for every headline number in §4 (adjoint coil-kerma importance map + closed-loop
breeder-for-shield optimizer) of `PAPER_DIRECTION.md` / `METHODS_DRAFT.md`. One subsection per result.
Each records device, geometry (shield build), emissivity, response, **objective type (transport solve
vs surrogate)**, rel.err / n, and the honest caveat. Frozen 2026-07-25; supersedes ad-hoc quotes of
these numbers elsewhere. If a number in the manuscript disagrees with this file, this file wins (or
one of them is stale — reconcile before submission).*

All results are for the **Landreman–Paul quasi-helically symmetric (QH) reactor**, a watertight DAGMC
model built from the real QUASR equilibrium (conformal radial build: plasma-facing wall / steel / Be
multiplier / FLiBe breeder / WC shield / coil) with the Wiedman QH coil filaments. Transport is
OpenMC 0.15.x on Ginsburg (x86, DAGMC build). The adjoint importance map is the OpenMC **random-ray
adjoint** solve (localized adjoint source on the coil guide curve, mainline
`FlatSourceDomain::set_local_adjoint_sources` — no feature-branch rebuild). The **quasi-axisymmetric
(QA) device is in progress** and appears in no number below.

## Objective-type legend (read before quoting any benefit)

- **[transport]** — a full OpenMC transport (or random-ray adjoint transport) solve. The importance
  map, reciprocity, and the angular P0/P1 comparison are all transport.
- **[surrogate]** — the closed-loop optimizer's *inner objective*: an exponential-attenuation model
  with **uncalibrated literature removal MFPs** (λ_shield = 8 cm, λ_breeder = 17 cm) and a linearized
  TBR. The importance map that *seeds* it is transport; the dose-reduction number itself is not yet
  OpenMC-verified. Calibration to OpenMC uniform-thickness scans is Phase C of the execution plan.

## Cross-cutting property — ψ†(r) is reactivity-independent

The adjoint importance ψ†(r) depends only on geometry and the coil response function; it contains **no
plasma-source information**. Switching the plasma reactivity profile (uniform ↔ Bosch–Hale core-peaked,
or any user scenario) is therefore a **free re-weight** of the attribution — the contributon
C(r) = S(r)·ψ†(r) is recomputed from the *single* committed adjoint solve with no new transport run.
This underlies R4 (the per-coil ranking is stable under reactivity re-weighting) and makes the
emissivity sensitivity in R2 cheap. Cited in `METHODS_DRAFT.md §2.x.ii` and `PAPER_DIRECTION.md §4`.

## Summary table

| # | Result | Device / geometry | Emissivity | Response | Objective | Value | rel.err / n | Caveat (short) |
|---|---|---|---|---|---|---|---|---|
| R1 | Forward↔adjoint reciprocity | QH / `dagmc_corr_w35f` (35 cm shield) | uniform (matched) | flat fast band | **transport** | Spearman **+0.76** raw, **+0.70** length-norm. | n=20 coils; QH forward tight (~4–5%) | coils 15 & 28 per-coil outliers; ranking holds |
| R2 | Closed-loop peak coil-dose reduction | QH coil-15 / `dagmc_corr_w20` | Bosch–Hale (physical) | flat fast band | **surrogate** | **−19%** | single peak coil | surrogate objective, uncalibrated MFPs |
| R2b | " (uniform, non-physical) | QH coil-15 / `dagmc_corr_w20` | uniform | flat fast band | **surrogate** | **−55%** | — | over-optimistic; over-weights shaped edge; bound only |
| R2c | " (no-attribution baseline) | QH coil-15 / `dagmc_corr_w20` | (flat field) | flat fast band | **surrogate** | **−1%** | — | no spatial signal to target; control |
| R3 | Angular P1 vs P0 on real coil (M7) | QH / `openmc-m7` build | (importance map) | flat fast band | **transport** | median P1/P0 **~1.006**; φ₀ **210°** unchanged | — | plasma-side only; not deep-shield |
| R3b | Angular P1 vs P0 on slab (M6) | controlled slab | — | — | **transport** | P1/P0 → **~1.7×** at 5 mfp (scattering bulk); **~3%** in clean streaming gap | — | slab test; motivates R3 |
| R4 | Reactivity robustness of ranking | QH | uniform → Bosch–Hale | flat fast band | **transport** (re-weight) | hot-5 set + peak coil (cell 20) **unchanged** | n=20 | free re-weight; ψ† reactivity-independent |
| R5 | Peak-coil identity nuance | QH / `dagmc_corr_w35f` | uniform (matched) | flat fast band | **transport** | forward-hottest **cell 15**; adjoint-hottest **cell 20** | n=20 | the two 15/28 outliers drive the mismatch |

---

## R1 — Forward↔adjoint reciprocity (magnitude/ranking validation)

- **What.** Reciprocity predicts INT S(r)·ψ†(r) dr for a coil equals that coil's neutron load from
  the full plasma source computed by an independent forward transport run:
  `INT S·ψ† ~ forward per-coil coil flux`. Validates the *magnitude/ranking* of the adjoint
  attribution (the *direction* — that φ₀ points at the coil — is validated separately in R3/placement).
- **Device / geometry.** QH, DAGMC `dagmc_corr_w35f` (35 cm shield build), matched between the forward
  per-coil tally and the adjoint maps.
- **Emissivity.** **Uniform**, deliberately, to match the uniform-source forward run (this is a
  validation, not a physical-ranking, use of S(r)).
- **Response.** Flat fast band (matched to the forward tally band).
- **Objective.** Transport — committed per-coil adjoint maps (kerma sweep) correlated against the
  committed forward per-coil flux (`percoil_*_corr_w35f.npz`); no new transport run.
- **Value.** Spearman **+0.76** (raw INT S·ψ†) and **+0.70** (length-normalized). A longer winding
  spreads the fixed ~40-point filament source over more plasma, so the raw integral over-weights long
  coils vs the forward flux *density*; dividing by coil arc length removes that geometric bias (both
  correlations are reported so the bias is visible, not hidden).
- **rel.err / n.** n = 20 coils. QH forward per-coil flux is statistically tight (rel.err ~4–5%).
- **Caveat (honest).** Coils **15 and 28** are per-coil outliers, most likely from filament
  source-matching on highly non-planar windings (coil-15's winding sweeps R ≈ 470→1775 cm; the matched
  filament centroid offset is ~40 cm). The overall ranking holds; individual coils can be misranked.
  Reducing these outliers (segment-length-weighted sampling / per-coil alignment check) is Phase D of
  the execution plan; do not tune to a target.
- **Provenance.** `shield_opt/fwcadis/reciprocity_check.py` (`--emissivity uniform`, matched w35f maps
  + `--percoil percoil_unpol_corr_w35f.npz`).

## R2 — Closed-loop peak coil-dose reduction (the optimizer payoff)

- **What.** The closed-loop breeder-for-shield optimizer trades breeder for shield on the (θ,φ) control
  grid to cut peak coil dose under a fixed outer envelope + TBR floor, seeded by the adjoint placement
  priority (contributon projected to the grid).
- **Device / geometry.** QH, coil-15, DAGMC `dagmc_corr_w20` (20 cm shield build).
- **Emissivity.** **Bosch–Hale** core-peaked reactivity (the physical default) for the headline;
  uniform reported only to bound the sensitivity.
- **Response.** Flat fast band.
- **Objective.** **SURROGATE.** The optimizer's inner loop uses D(δ) = D₀·exp(−kδ),
  k = 1/λ_shield − 1/λ_breeder, with **uncalibrated literature MFPs λ_shield = 8 cm, λ_breeder = 17 cm**,
  plus a linearized TBR and light regularizers, minimized by Bayesian optimization (GP + EI). The
  adjoint importance map that provides the baseline dose *shape* D₀ is a transport solve; the
  dose-reduction number itself is **not yet OpenMC-verified**.
- **Value.**
  - **−19%** (Bosch–Hale) — **the honest headline**, adjoint-informed baseline shape.
  - **−55%** (R2b, uniform emissivity) — non-physical; uniform over-weights the strongly-shaped plasma
    edge where the contributon is most localizable, ~3× over-stating the localizable benefit. Reported
    only as an upper bound on the emissivity sensitivity.
  - **−1%** (R2c, no-attribution uniform baseline field) — the control: a flat D₀ has no spatial signal
    to target, spreads shield thin, and achieves essentially nothing. This is the contrast that shows
    the *adjoint attribution*, not the optimizer, is what buys the reduction.
- **rel.err / n.** Single peak coil (coil-15); the number is a surrogate-model peak reduction, not a
  transport-tallied dose, so a transport rel.err does not yet apply — that is exactly what Phase C
  calibration will attach.
- **Caveat (critical).** This is a **surrogate-objective** number. Calibration of λ_shield / λ_breeder
  to OpenMC uniform-thickness dose-vs-δ scans (Phase C) is required before −19% is a quantitatively
  trusted, OpenMC-verified dose cut; calibration may move it up or down. Report it as "surrogate" in
  every venue until Phase C closes.
- **Provenance.** `shield_opt/adjoint_closed_loop.py` (`--emissivity bosch_hale` headline;
  `--emissivity uniform` for R2b; uniform baseline field for R2c). Surrogate MFPs in
  `attenuation_surrogate.py`.

## R3 — Angular sensitivity, P1 vs P0 on the real coil (M7)

- **What.** Does direction-resolved (P1, first-order anisotropic scatter) adjoint change the coil
  importance map vs scalar (P0, transport-corrected)? Compares P0 and P1 adjoint importance maps for
  the same coil.
- **Device / geometry.** QH, real DAGMC, `openmc-m7` build (angular C++; `legendre_order = 1` forced
  via a targeted MGXS-library patch; both orders use an untransport-corrected total XS so the
  comparison isolates the explicit angular effect and avoids the near-void negative-XS reject).
- **Response.** Flat fast band.
- **Objective.** Transport (two random-ray adjoint solves compared).
- **Value.** Median **P1/P0 ≈ 1.006**; the shield-placement target **φ₀ = 210° is unchanged** between
  P0 and P1. Conclusion: **a scalar adjoint suffices for coil placement** — the plasma-side attribution
  the optimizer actually uses lives in the optically thin region where angular effects are small.
- **R3b (controlled slab, M6).** On a deep-penetration slab the P1/P0 importance ratio grows with depth
  into scattering shield — up to **~1.7×** at 5 mfp in the scattering bulk — but stays within **~3%**
  across a clean streaming gap. This is why R3's near-unity result on the coil is a *finding* (the
  attribution lives in the thin region), not a null measurement.
- **Caveat (honest).** The check is **P1-only**; higher moments are untested. "Scalar suffices" is
  specific to coil *placement* (plasma-side attribution) and does **not** extend to deep in-shield
  fields, where M6 shows real angular growth.
- **Provenance.** `shield_opt/fwcadis/compare_p0_p1.py` (M7); the slab test is M6.

## R4 — Reactivity robustness of the per-coil ranking

- **What.** Does realistic reactivity change *which coil to protect*? Because ψ†(r) is
  reactivity-independent, re-weighting S(r) from uniform to Bosch–Hale is free.
- **Device / geometry.** QH; re-weight of the committed adjoint maps.
- **Value.** The per-coil **hot-5 set and the peak coil (cell 20 by adjoint) are unchanged**
  uniform → Bosch–Hale. Realistic reactivity re-scales magnitudes but does not reorder which coil is
  hottest.
- **Objective.** Transport (single adjoint solve) + free re-weight; no new transport run.
- **Caveat.** Follows directly from the reactivity-independence property; it is a consistency result,
  not an independent measurement.
- **Provenance.** `shield_opt/adjoint_placement.py` (`contributon(..., emissivity=...)`),
  `reciprocity_check.py --emissivity {uniform,bosch_hale}`.

## R5 — Peak-coil identity nuance (state it, don't paper over it)

- **What.** The forward-hottest and adjoint-hottest coils are **not the same cell**: forward per-coil
  flux peaks at **cell 15**, while the adjoint attribution peaks at **cell 20**.
- **Why.** Cell 15 is one of the two per-coil reciprocity outliers (R1). The mismatch is a direct
  consequence of the filament source-matching error on the highly non-planar cell-15 winding, not a
  contradiction of reciprocity — the overall rank correlation is still Spearman +0.76.
- **Caveat.** Quote the peak coil carefully: "adjoint-hottest cell 20 / forward-hottest cell 15,
  overall Spearman +0.76." The Phase-D filament fix targets exactly this.
- **Provenance.** `reciprocity_check.py` (per-coil table), `dagmc_corr_w35f`, n=20.

---

## Provenance & reproduction notes

- **Code.** `shield_opt/fwcadis/adjoint_importance.py` (the random-ray adjoint importance solve +
  structure gate), `shield_opt/adjoint_placement.py` (contributon + (θ,φ) placement priority),
  `shield_opt/adjoint_closed_loop.py` (closed-loop optimizer), `shield_opt/fwcadis/reciprocity_check.py`
  (R1/R4/R5), `shield_opt/fwcadis/compare_p0_p1.py` (R3).
- **Structure gate.** Every saved adjoint map passes a hard gate before downstream use: non-empty
  (`total > 0`), localized (nonzero-voxel fraction < 0.98), structured (peak/mean > 1.5), and converged
  (median rel.err < 0.5; QH runs comfortably inside this). See `adjoint_importance.structure_gate`.
- **Group structure.** 13-group fusion structure with fast detail above 0.1 MeV
  (`GROUP_EDGES` in `adjoint_importance.py`), shared between the adjoint and forward decks so the maps
  are directly comparable.
- **Raw artifacts** (per-coil adjoint maps, forward `percoil_*_corr_w35f.npz`, coil centroids/filaments)
  live on Ginsburg; committed drop-ins used by the figure/optimizer scripts are under
  `shield_opt/data/adjoint/` (e.g. `coil15_fil_flat.npz`) and `shield_opt/data/qh_freestream_fluxmap.npz`.
- **Pending (do not report as done).** Phase C surrogate-MFP calibration (converts R2 to a transport
  number); Phase B QA bring-up (a second device); Phase D filament source fix (the R1/R5 outliers).
