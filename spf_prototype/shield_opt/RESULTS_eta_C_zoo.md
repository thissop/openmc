# Zoo-wide polarization law: cheap coherence C predicts the transport steering effect η (n=64)

**Result.** Across **64 QA/QH reactor devices** run through the full conformal free-streaming
transport sweep on Ginsburg, the cheap **field-coherence C** (computed from the coils alone, no
transport) predicts the **polarization steering effect η** (the fractional change in free-streaming
coil flux when the DT source emits perpendicular vs parallel to B̂):

- **|η| vs C: Spearman ρ = +0.464 (p = 1.15×10⁻⁴), Pearson r = +0.470 (p = 9×10⁻⁵), n = 64.**
- Effect magnitude: perpendicular emission changes free-streaming coil flux by up to **−16%**,
  parallel by up to **+16%** (symmetric, |η| ∈ [0, 0.161]).
- Coherence-stratified: **low-C devices |η| = 0.072 ± 0.022; high-C = 0.095 ± 0.044** — higher field
  coherence ⇒ more polarization leverage, as the mechanism predicts (a coherent B̂ field lets the
  anisotropic emission steer coil flux; a decoherent field averages it out).

**Why it matters.** This validates the geometry-first "cheap predictor" thesis at zoo scale: C is a
coils-only scalar (seconds to compute), η is an expensive per-device transport quantity (conformal
DAGMC + Monte-Carlo). C→η with ρ≈0.46, p~1e-4 means we can screen the polarization payoff of a device
without running transport. It is the polarization-side companion to the first-wall elongation→ξ law and
the magnet-side adjoint attribution.

## Provenance / honesty
- Device set: the 178-device flat-coherence conformal-sweep manifest (`production_scaleup.json`), run as
  SLURM array 9193039 on Ginsburg. Of 179 records: **64 done (36% usable)**, 69 error, 34 audit_failed,
  12 geometry_folded.
- The **69 "error" are DAGMC-leaky-geometry aborts** (lost particles) — the watertightness pre-flight
  gate (`magnet_suite/dagmc_watertight_check.py`) would convert these to clean early rejects and, if the
  leaks are fixable (higher mesh resolution), recover much of that yield → the honest usable-n could rise
  well above 64. Not yet wired into the conformal sweep (the sweep was running live).
- The 34 audit_failed are legitimate coil-fit-quality rejections (B·n on the LCFS too large) — see
  `RESULTS_zoo_coil_law.md`; not tuned away.
- ρ is Spearman (monotone, robust to the η distribution). Data: `data/eta_C_zoo.csv` (64 rows: id, class,
  nfp, aspect, C, S_phi, Δφ_perp/par ± sd). Figure: `figs/eta_C_zoo.png`. Regenerate:
  `python shield_opt/eta_C_zoo_analyze.py`.

## Caveats
- QA/QH only (QUASR shipped coils); QI enters once own-coil design lands. Effect is *free-streaming*
  (pre-blanket) coil flux — the through-shield magnet attribution is the separate adjoint pipeline.
- n=64 with C ∈ [0.66, 1.0]; the low-coherence tail (C<0.66) is thin (audit/geometry attrition there).
