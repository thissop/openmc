# PAPER_DIRECTION — the single capstone paper

*Authoritative plan for the one paper that closes out this summer's stellarator-neutronics work.
Rewritten 2026-07-23 after the Peterson advising meeting (see "What changed" below). Supersedes both
the 2026-07-08 SPF-centric version of this file and the exploratory
`docs/notes/planning/PAPER_DIRECTION.md`.*

## What changed and why (2026-07-23 Peterson meeting)

Ethan (advisor, MIT PSFC) reframed the whole paper. The old thesis put the **spin-polarized source**
at the center ("where and why MC is necessary for polarized stellarator wall loading"). Two problems:
Bae et al. 2025 already published a polarized source in OpenMC, so polarization is not our
methodological novelty; and the SPF story, while real, is a narrow lever, not a spine a paper can
stand on. Ethan flagged **broad geometric stellarator neutronics** — how device geometry sets both
first-wall **and magnet (coil)** neutron loading — as an underdeveloped area with room for an
original contribution, and pointed at **adjoint importance mapping for magnet protection** (an adjoint
source = the coil kerma response, solved with OpenMC's random-ray adjoint solver, whose adjoint flux
everywhere *is* the importance map) feeding a **closed-loop breeder-for-shield optimizer** as the
genuinely novel methodological center. A prior-art check supports that the closed-loop
breeder↔shield magnet optimizer plus adjoint coil→plasma attribution on a 3-D stellarator is
unpublished.

So: **geometry is now the spine; the adjoint importance map + optimizer is the methodological
novelty; SPF is demoted to one supplementary physics lever layered on top.**

## Thesis (one sentence)

> Device **geometry** — elongation, rotational transform, coil standoff — sets the *structure* of both
> first-wall and magnet neutron loading on a stellarator, and an **adjoint coil-kerma importance map**
> (random-ray adjoint in OpenMC) turns that structure into an actionable, closed-loop
> breeder-for-shield magnet-protection design lever; spin polarization is a modest, device-dependent
> actuator layered on top, not the driver.

This is a **methods + honest-findings** paper (same genre as Peterson's TokamakSource paper), not a
discovery paper. It delivers (i) a geometry→load-structure decomposition across a device zoo, (ii) an
adjoint importance-mapping method for coil protection, and (iii) SPF as a characterized secondary
lever.

## Abstract sketch (one paragraph)

> First-wall and magnet neutron loading on stellarators are usually reported device-by-device. We
> instead decompose the load into what **geometry** fixes versus what actuators can move. Across an
> n=83 free-streaming device zoo we find elongation drives first-wall load *concentration* (Chatterjee
> ξ = 0.96) and rotational transform drives the *peak*, so geometry sets the load structure before any
> plasma or shielding choice. We then introduce an **adjoint importance map for magnet protection**:
> an adjoint source equal to the coil-kerma response, localized in the coil cell and solved with
> OpenMC's random-ray adjoint solver, yields an importance field over the plasma that attributes coil
> heating to its phase-space origins and drives a closed-loop breeder-for-shield optimizer. Finally we
> layer a native spin-polarized DT source on top and show polarization is a ~10–15% device-dependent
> lever (parallel B/C emission flattens QA magnet peaking ~17%), distinct from Bae et al. 2025's
> precomputed tokamak source. Realistic core-peaked (Bosch–Hale) emissivity compresses the device
> spread and collapses the uniform-emissivity QA-vs-QH gap, and we flag that the raw QA-vs-QH magnet
> comparison is confounded by coil standoff, not quasisymmetry.

## Novelty framing (what to claim, what NOT to claim)

| Claim | Status | Prior art to cite & differentiate |
|---|---|---|
| **Geometry→load-structure decomposition** (geometry sets structure; polarization/emissivity are bounded modulators) across a device zoo | our contribution | ParaStell (Davis et al.) = manual parametric radial-build sweeps, per-device, no structure decomposition; Lyytinen 2024 (Serpent2, parametric geometry, isotropic) |
| **Adjoint coil-kerma importance map** (adjoint source = magnet response, random-ray adjoint, flux = importance field) for a **3-D stellarator** | novel (methodological center); **built + validated on real QH DAGMC** (reciprocity +0.76) | FW-CADIS / CADIS use adjoint for variance reduction, not as a coil→plasma *attribution* map; Miralles-Dolz / Thea use random-ray FW-CADIS for WW, not closed-loop shield design |
| **Closed-loop breeder-for-shield magnet-protection optimizer** driven by that importance map | novel; **built** (surrogate benefit ~19% Bosch–Hale; OpenMC calibration pending) | no published closed-loop breeder↔shield magnet optimizer on a 3-D stellarator |
| First-wall load **concentration law** (elongation→ξ, iota→peak) on a real device zoo | our result | — |
| **Native, on-the-fly** polarized fusion source in OpenMC (supplementary) | defensible increment | Bae 2025 = precomputed static point cloud, one fixed tokamak |

**Do NOT write:** "first polarized source in OpenMC" (Bae); "first MC stellarator source"
(Lyytinen/ParaStell); "QA magnets are intrinsically better protected than QH" (that gap is a
coil-standoff confound — see Gaps).
**DO write:** geometry sets load structure; adjoint importance mapping + closed-loop breeder-for-shield
is the novel method; SPF is a native, general, on-the-fly *increment* and a bounded device-dependent
lever.

## The spine (state this explicitly)

Three linked claims, in order:

1. **Geometry sets the structure of the load.** Elongation → first-wall *concentration* (Chatterjee
   ξ = 0.96, shape-anisotropy ξ = 0.83, n=83 free-streaming zoo); rotational transform (iota) → *peak*.
   This is prior-to-actuator: it holds before plasma profile or shielding choices.
2. **The load is attributable, and therefore optimizable, via adjoint importance mapping.** An adjoint
   source = coil-kerma response in the coil cell → random-ray adjoint solve → adjoint flux = the
   importance map (which plasma/phase-space regions drive coil heating) → closed-loop breeder-for-shield
   optimizer. *(Methodological centerpiece; **built + validated** on the real QH DAGMC — reciprocity
   Spearman +0.76, scalar-adjoint-suffices angular check; closed-loop benefit still surrogate pending
   OpenMC calibration. See `RESULTS_adjoint.md`.)*
3. **Actuators are bounded modulators on top of the geometric structure.** A
   geometry × polarization × emissivity decomposition: geometry sets structure; polarization is a
   ~10–15% device-dependent lever (parallel B/C flattens QA magnet peaking ~17%); realistic
   (Bosch–Hale core-peaked) emissivity is a ~6% modulator that compresses the device spread and
   collapses the uniform-emissivity QA-vs-QH gap.

## Structure (sections)

1. **Introduction** — geometric neutronics motivation (first wall AND magnet); prior art up front
   (ParaStell manual sweeps, Lyytinen, Bae, Peterson TokamakSource, Schwartz analytic); the gap =
   no structure decomposition, no adjoint coil-attribution / closed-loop shield optimizer.
2. **Geometry sets the load structure** *(spine, part 1)* — the n=83 free-streaming zoo; elongation→ξ
   concentration (Chatterjee ξ=0.96), iota→peak; concentration/Chatterjee metrics defined.
3. **Native sources and verification** *(the vehicle, condensed)* — the compiled StellaratorSource
   (VMEC √g/geometry/b̂, TokamakSource-pattern extension, thread-safe sampler), Tier 1–2 sampler +
   NWL-vs-Schwartz verification, free-streaming analytic↔MC cross-check on real QUASR 59509. Kept
   tight — it is infrastructure, not the headline.
4. **Adjoint importance mapping for magnet protection** *(methodological centerpiece)* — adjoint
   source = coil-kerma response localized in the coil cell; random-ray adjoint solve; adjoint flux =
   importance map; coil→plasma attribution; the closed-loop breeder-for-shield optimizer built on it.
5. **Actuators as bounded modulators** *(spine, part 3)* — geometry × polarization × emissivity
   decomposition (24-device factorial); the complete QA/QH × polarization magnet matrix (B/C flattens
   QA peaking ~17%, device-dependent); Bosch–Hale emissivity compression.
6. **Spin polarization as a supplementary lever** — native (a,b,c) sampler physics, angular/rate
   separation, a₂-linearity identity, *carefully* differentiated from Bae 2025 (precomputed static
   tokamak point cloud vs native/general/on-the-fly/stellarator). Explicitly framed as one actuator
   inside §5, not a headline.
7. **Discussion / Conclusion** — geometry sets structure; adjoint importance mapping makes coil
   heating attributable and optimizable; SPF is a bounded device-dependent actuator; the standoff
   confound and emissivity compression bound the QA-vs-QH story.

**Appendices:** concentration/Chatterjee metric definitions; the n=83 zoo table; the coil-standoff
attribution note; sampler stratified-sampling + a₂-linearity derivation; QUASR mining tooling;
Bae (a,b,c) symbol-collision footnote.

## IN / OUT

- **IN:** §2 geometry decomposition + ξ concentration law; §4 adjoint importance map + closed-loop
  breeder-for-shield optimizer; §5 geometry×polarization×emissivity decomposition + the QA/QH magnet
  matrix; §6 native SPF as a supplementary lever framed against Bae; the sampler verification (§3).
- **OUT (agreed / demoted):** the old SPF-centric thesis ("where MC is necessary for polarized
  wall loading") — demoted to one point inside §3/§6; the anarrima analytic-companion arm (keep as a
  Discussion outlook / cross-check only, not a section); TBR arithmetic beyond what the breeder-for-
  shield loop needs; the box-torus 3× scattering-dilution as a *headline* (it stays as a supporting
  SPF-lever detail, not the paper's spine).

## Committed results → sections (what we already have)

- **§2 (geometry):** n=83 free-streaming zoo; Chatterjee ξ=0.96 (elongation→concentration), shape
  anisotropy ξ=0.83, iota→peak. Committed (`FINDINGS_geometry_peaking.md`, zoo JSON).
- **§5 (actuators):** 24-device geometry×polarization×emissivity factorial (`zoo_factorial.json`);
  complete QA/QH × {unpol, A, B/C} magnet matrix (peaking/maxflux table, `OVERNIGHT_2026-07-22.md`);
  B==C confirmed on both devices (sampler check); Bosch–Hale emissivity compression.
- **§3 (verification):** Tier 1–2 sampler + NWL-vs-Schwartz; free-streaming analytic↔MC on QUASR 59509.
  Committed.
- **§6 (SPF):** native (a,b,c) StellaratorSource, a₂-linearity identity. Committed.
- **§4 (adjoint) — BUILT + VALIDATED.** The novel centerpiece now runs. The adjoint coil-kerma
  importance map is solved with OpenMC's **random-ray adjoint solver** (localized adjoint source =
  coil response on the coil guide curve, mainline `set_local_adjoint_sources`, no feature-branch
  rebuild) on the **real QH DAGMC** geometry, and the full chain — importance map → contributon
  C=S·ψ† → placement → closed-loop breeder-for-shield optimizer — is committed and validated.
  Validation: forward↔adjoint **reciprocity Spearman +0.76** (n=20 coils, matched geometry `w35f` +
  response band, uniform-matched emissivity); the **angular P1-vs-P0 check (M7)** shows a scalar
  adjoint suffices for coil placement (median P1/P0 ~1, placement φ₀ unchanged); the closed-loop
  optimizer cuts peak coil dose ~19% under realistic Bosch–Hale emissivity. **Key methodological
  property:** the adjoint importance ψ†(r) is **reactivity-independent** (geometry + coil response
  only; no plasma-source information), so switching the reactivity profile (uniform ↔ Bosch–Hale, or
  any scenario) is a *free re-weight* of the attribution needing no new transport solve — which is why
  the emissivity study and per-coil ranking are cheap, and why the hot-coil set and adjoint-peak coil
  (cell 20) are unchanged uniform→Bosch–Hale. All headline numbers,
  their run parameters, and their honest caveats are frozen in `paper/RESULTS_adjoint.md`. The
  earlier WW-reciprocity attempt (MAGIC weight windows: exploded, then solved-but-slow) was the
  precursor that redirected us to this random-ray adjoint (TRRM) path. **Remaining caveat:** the
  closed-loop *benefit* is a surrogate-objective number (uncalibrated literature MFPs); OpenMC
  transport calibration is pending (see Honest gaps / Open questions).

## Honest gaps (state in Limitations)

- **§4 closed-loop benefit is a surrogate number, not yet OpenMC-calibrated.** The adjoint importance
  map and reciprocity/angular validations are full transport solves and are done. But the optimizer's
  *inner loop* uses an exponential-attenuation surrogate with **uncalibrated literature removal MFPs**
  (λ_shield=8 cm, λ_breeder=17 cm); the ~19% (Bosch–Hale) / ~55% (uniform) peak-dose reductions are
  **surrogate-objective** figures given the adjoint-informed vs uniform baseline *shape*, not
  OpenMC-verified dose cuts. Calibration to OpenMC uniform-thickness scans is the pending step;
  label the benefit as surrogate until then.
- **§4 is demonstrated on a single device (QH).** The full adjoint→optimizer→validation chain runs on
  the Landreman–Paul QH reactor only; QA is in progress. Do not generalize the method's numbers beyond
  QH.
- **§4 per-coil reciprocity has two outliers.** Reciprocity holds in rank (Spearman +0.76 over 20
  coils), but coils 15 and 28 are per-coil outliers — most likely filament source-matching on highly
  non-planar windings; the forward-hottest coil is cell 15 while the adjoint-hottest is cell 20. The
  overall ranking holds; individual coil magnitudes carry that caveat.
- **§4 angular check is P1-only.** "Scalar suffices for coil placement" is specific to the plasma-side
  attribution (optically thin, median P1/P0 ~1); it does not extend to deep in-shield fields, where a
  controlled slab (M6) shows P1/P0 growth to ~1.7× at 5 mfp in the scattering bulk.
- **QA-vs-QH magnet gap is a coil-standoff confound.** QA coils sit ~1.9× further (global-min 3.11 m
  vs 1.63 m; ~148 cm extra standoff); 1/r² plus the longer attenuated path explains the ~25–40× lower
  QA coil flux. This is **not** a quasisymmetry advantage. The *within-device* polarization effects
  (B/C flattens QA peaking ~17% vs QH ~3%) survive because each device is compared to its own unpol.
- **QA statistics are loose** (relerr ~23–25%): trends solid, exact factors ±15%. Tighten any QA
  number that enters the paper. QH is tight (relerr ~4–5%).
- Free-streaming zoo is uncollided/single-energy (14.1 MeV); the scattering correction is
  characterized only on box-torus/conformal-59509 geometry, not the full zoo.
- Emissivity compression measured on the factorial subset, not all 83 devices.

## Open questions (updated 2026-07-23)

1. **Calibrate the closed-loop benefit to OpenMC transport (new critical path).** The adjoint
   importance map is built and reciprocity/angular-validated; what remains is turning the optimizer's
   *surrogate* benefit into an OpenMC-grounded number. Run full-transport coil dose vs uniform
   breeder→shield thickness, fit λ_shield/λ_breeder, replace the literature-scale defaults (8/17 cm),
   and re-run the closed loop. The ~19% must be reported as surrogate until this closes; the
   calibration may move it up or down.
2. **Closed-loop breeder-for-shield optimizer scope** — how many design iterations / which shield
   parameters for a first demonstration? Single coil cell or the full coil set?
3. **Fair QA-vs-QH comparison** — normalize by standoff or use same-standoff coil sets, so §5 has a
   clean configuration comparison rather than a confounded one. In or out of first submission?
4. **Device-set scope** — n=83 free-streaming for §2; how many transport/adjoint devices for §4/§5?
   The n=25 VMEC/predictor batch (Ginsburg 8885457/58) feeds §5.
5. **SPF framing paragraph vs Bae 2025** — still the reviewer's first objection; now easier because
   SPF is explicitly supplementary. Draft it early, in §6.
6. **Target venue** — Nuclear Fusion still fits (same journal as Bae, Lyytinen, Lion, Parisi), and the
   geometry+adjoint framing arguably fits it better than the SPF-only version did.

## Data/figure discipline

`paper/data/` holds **frozen** result snapshots (CSV/JSON extracted from the RESULTS/overnight docs
and the conformalmap npz) with `PROVENANCE.md` naming the source of every number. `paper/src/` figure
scripts read only from `paper/data/` and write to `paper/figures/`. smplotlib, minimalist. One script
per figure, fixed inputs, regenerable. New figures needed for the pivot: the ξ concentration-vs-
elongation zoo scatter (§2), the adjoint importance map / contributon + closed-loop placement (§4,
built — regenerate from the committed adjoint maps), the geometry×pol×emissivity
decomposition bars (§5).

---

*Duplicate note: `docs/notes/planning/PAPER_DIRECTION.md` is the **stale** 2026-07-07 two-paper
reconstruction (pre-pivot, and even pre the 2026-07-08 single-paper consolidation). It is now doubly
out of date. **This file (`paper/PAPER_DIRECTION.md`) is canonical.** Recommend deleting or adding a
one-line "superseded" banner to the notes copy.*
