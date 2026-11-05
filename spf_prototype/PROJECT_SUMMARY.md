# Project Summary — Physics-Informed Fusion Sources & Response-Aware Neutronics in OpenMC

**Abstract.** This project extends OpenMC with physics-informed fusion neutron sources —
a native stellarator plasma-fluxmap source and a spin-polarized-fuel (SPF) birth-direction
kernel — and uses them to ask whether magnet/coil neutron loading in a stellarator can be
reduced by *where* shielding is placed, not just how much is added. Three arcs: (1) a
compiled `StellaratorSource` enabling device-scale QUASR stellarator neutronics directly in
OpenMC; (2) a polarization-aware DT birth-direction sampler (Schwartz 2025 physics),
source-agnostic and validated against an analytic tokamak test problem; (3) an automated,
attribution- and mass-aware shield-placement study on real QUASR/ParaStell/DAGMC geometry,
whose honest headline is a *negative/structured* result: naive attribution does not rank
placements, and the mass accounting needed to even ask "equal-mass" questions required
building an independent volume-calculation pipeline. All three arcs share one OpenMC
transport stack and a hard-gated, non-cherry-picked verification discipline.

---

## 1. Stellarator neutron source for OpenMC (`StellaratorSource`)

**What was built.** A compiled `CompiledSource`-style birth-position/direction sampler that
draws fusion-neutron births proportional to `S(ρ)|√g|` over the true VMEC plasma equilibrium
(R,Z fluxmap; default `S(ρ)=1−ρ²`) for QUASR QH and QA configurations, with local B̂ derived from
`(B_R, B_φ, B_Z)`, and supports the polarized (a,b,c) parameterization. It runs as a native
source inside a source-built OpenMC transport tree, enabling first-wall loading and
coil/magnet loading calculations directly on DAGMC stellarator geometry (no external
coupling to a private codebase).

**Method / integration reality.**
- Lives in exactly one tree: `spf_work/openmc_src` on Ginsburg (`openmc.StellaratorSource`).
  A second, superficially identical clone (`~/src/GitHub/openmc-spf`, from a bundle) never
  received this code and is used only for adjoint/importance runs — the two trees must not
  be confused (cost a real forensic detour once).
- New-device standup is fragile and was debugged by hand for QA: hardcoded QH assumptions
  (chdir, coil-cell ranges, fluxmap choice) and a units bug (fluxmap in meters vs OpenMC's
  cm convention, giving TBR≈0 and ~70x coil flux) all had to be fixed per device.
- Correctness gates now exist and are run on every device: `assert_physical` (TBR range,
  coil-flux range vs a QH reference, relative error), `assert_shield_reduces`, and a
  preflight that samples the source + `find_cell` before committing to a full transport run.

**Status & caveats.**
- Working and is the workhorse source for all Arc 3 transport (coil dose, per-coil flux,
  FW-CADIS studies).
- `openmc.lib.Cell.fill` is not implemented for DAGMC in this build; coil/cell identification
  goes through pymoab or `summary.h5`, not `.fill`.
- Build/runtime is x86 (Ginsburg)-specific and version-sensitive (surface-source attribute
  name, DAGMC build flags — see Arc 3 caveats); this integration work is itself a real,
  hard-won deliverable, not a side detail.

---

## 2. Spin-polarized fusion (SPF) source kernel

**What was built.** A polarization-aware DT neutron **birth-direction** sampler implementing
the Schwartz (2025) differential cross section, designed to be **source-agnostic**: a kernel
that composes with any spatial source (analytic filamentary ring, or the Arc-1 stellarator
fluxmap source).

**Physics (authoritative form used).**
- Deuteron/triton spin fractions → three collision-mode fractions `a, b, c` (`a+b+c=1`).
- `dσ/dΩ = (σ₀/2π)·[(3/4)a·sin²θ + ((2/3)b + (1/3)c)·(1/4 + (3/4)cos²θ)]`, θ measured from
  the local B̂.
- Total-rate factor decoupled from angular shape: `η = a + (2/3)b + (1/3)c` (in units of σ₀).
  Corner checks (the sampler's oracle): non-polarized `η=2/3` (**not** σ₀); pure-A `η=1`
  (peaks ⊥ B̂); pure-B `η=2/3`; pure-C `η=1/3` (**half** the non-polarized rate).
- **Key physics correction relative to informal notes:** B and C modes share the *same*
  angular shape `1/4 + (3/4)cos²θ` — differing only in total rate. C is **not isotropic**,
  a common secondary-source mistake.

**Sampler design.** Stratified/inverse-CDF sampling of the two angular shapes in a
B̂-aligned local frame; robust Gram-Schmidt (least-aligned-axis) rotation to global frame to
avoid the degenerate-alignment pitfall; thread-safe by construction (no mutable state,
`openmc::prn` only). Validated against an analytic tokamak filamentary-ring test problem
(the Schwartz square-cross-section-torus NWL comparison) as the tiered verification target
(direction-sampler-alone tier, then full-stack NWL-vs-analytic tier), per the task
specification in `CLAUDE.md`.

**Status & caveats.**
- Framed deliberately as complementary/independent, not "first-ever": Bae et al. (2025,
  *Nucl. Fusion* 65 086051) independently published a polarized (a,b,c) source in OpenMC
  (tokamak, precomputed points). This work's differentiators are: native/differentiable,
  local-B̂ (not precomputed), and stellarator-capable via composition with `StellaratorSource`.
- An earlier conformal-QA smoke study (separate from the tiered verification plan) found an
  indicative, non-obvious directional sign: the **parallel** mode (B/C shape) raises coil fast
  flux (+13.5% free-streaming, +8.5% with scattering vs unpolarized) and **perpendicular**
  (A) lowers it (−15.3%/−11.3%) — the opposite of a documentation typo that had circulated.
  Its scatter/free-streaming estimate η = 0.63 ± 0.22 has large MC uncertainty. Treat the
  **direction of change** as the smoke-study result; the numerical percentages and η are not
  headline-ready without the fuller tiered verification.
- Scope explicitly excludes (this prototype): alpha-channel emission, spatially varying
  polarization, depolarization, non-convex geometry, energy-angle coupling beyond a
  decoupled monoenergetic/Ballabio-broadened option.
- A planned second, purely analytic/numerical-methods paper (perturbative field-period
  series + closed-form singularity-subtraction quadrature for off-axisymmetric NWL) is
  considered the stronger-novelty half of the SPF work, since Schwartz/anarrima is
  axisymmetric-only; open framing decisions are tracked separately.

---

## 3. Automated, response-aware shielding placement

**Goal metric:** peak magnet/coil neutron loading (not first-wall NWL) — the practical
engineering constraint. **Central question:** can non-uniform shield placement beat uniform
bulk shielding at fixed mass?

**What was built / method.** A pipeline on real QUASR-derived QH ParaStell/DAGMC geometry:
a direct-forward transport optimizer (periodic-BC one-field-period model, validated
essential — an open-boundary wedge washes out the redistribution lever entirely, 1.03x vs
the real 4.2x), adjoint-importance attribution (contributon = source × adjoint flux on the
coil guide curve), a 12-patch (4×4 sector) finite-difference "ground truth" response matrix,
an independent 40M-sample `VolumeCalculation`-based exact-mass pipeline, and a clearance
(coil-exclusion) field.

**Key results, with honest status:**
- **Kill-shot (done):** attribution-guided non-uniform shielding at fixed envelope cuts peak
  coil fast flux **4.2×** vs uniform — but the peak **migrates** from coil 20 to coil 25.
  Multi-coil placement (spread and greedy variants) both **lost** to single-coil placement;
  coil load is collective/inboard-dominated, not per-coil-local.
- **Patch-worth n=12 verdict (done, negative result):** the clean n=4 signed-worth
  correlation (Spearman +0.80) **did not generalize** to n=12 (all correlations n.s.,
  p≈0.3–0.6). Scalar contributon worth is in fact **anti-correlated** with true
  finite-difference worth on some patches (one patch that attribution ranks highest actually
  *backfires*, +41% dose). Verdict: **signed worth is a backfire guardrail, not a ranker.**
- **Structured 12-patch response matrix (done, scoped):** 98% of response entries exceed 2σ
  (real, not noise); inboard sectors robustly reduce peak coil load, wrong sectors robustly
  backfire (the *sign pattern* is robust) — but the absolute ~5× magnitude is
  **baseline-confounded** (patch builds used a slightly different shield-outside offset than
  the reference baseline) and the hotspot still migrates. A rung-2 attempt to correlate
  adjoint importance with finite-difference worth via a circular-torus proxy for the true
  VMEC surface gave an implausible anti-correlated result and was correctly flagged and
  **blocked** (not reported) per the project's standing integrity rule, pending an exact
  surface projection.
- **Mass method (done, authoritative):** a 40M-sample direct `VolumeCalculation` on the real
  DAGMC shows the *same* +20 cm shield patch costs **163–400 tonnes depending on sector**
  (2.45× spread) — "equal thickness ≠ equal mass." Every individually-built patch is
  *lighter* than the uniform +2.55 cm control (129–366 t lighter, up to 8σ), which
  quantitatively **refutes** an earlier "patch ≈ uniform at equal mass" claim. A cheaper
  ParaStell surface-area×thickness proxy was found to be shape-accurate but biased
  +6–19% high (a proposed correction made it worse, +20%, and was reverted) — the lesson
  generalized into the mass-method reference: use MC for absolute mass, surface method only
  for shape/field structure.
- **Clearance field (done):** minimum available radial growth toward the coils is only
  **4.7 cm**, located inboard (poloidal ≈173°) — exactly where coil-20 loading peaks. This is
  the project's concrete, quantified statement of "attribution ≠ intervention worth ≠
  feasible controllability": the place shielding is most wanted is also the place with least
  room to add it.
- **FW-CADIS / variance-reduction side investigation (done, pivotal negative-then-positive
  finding):** an earlier "coils are unreachable even at 24M histories" result that motivated
  a multi-week FW-CADIS/random-ray effort turned out to be an **artifact of the
  surface-source-replay pipeline** (black-inside cells killing coil-reaching neutrons), not
  real deep-penetration physics. In the nominal high-flux benchmark, a plain analog
  2M-history forward run resolves all 8 coils/period at 2–4% relative error in ~6 minutes.
  Heavily shielded designs can leave low-flux coils noisier; the matched-mass comparison
  therefore gates the **peak objective** separately (the 4M `famU` peak is 3.5%, while only
  5/8 low-to-high coils are below 10%). **Consequence: the shield optimizer uses direct-forward
  transport, no VR/SSR needed for the peak coil-load metric.** (The native
  FW-CADIS random-ray path was still successfully debugged and works for CSG and DAGMC
  smoke tests — miss rates 19–43% on the one-period model — but none of the three candidate
  weight-window files (native FW-CADIS, reciprocity-CADIS, MAGIC) were judged production-ready
  for the SSR-replay coil calculation; the direct-forward finding made this moot for the
  main optimizer.)
- **DONE, four-family direct-transport comparison:** matched-mass design *families* (uniform,
  clearance-driven, coil-20-guided/response-aware, and a linear-minimax stress test over the full 12×20
  sensitivity operator) built as design inputs at a common target mass M\*=2.16×10⁶ kg
  (0.00% mass-matching error, recalibrated to the MC-authoritative scale). At 4M histories,
  clearance-driven is +0.2% ± 4.8% and response-aware is +0.6% ± 4.9% relative to uniform;
  both are statistically unresolved and retain the c23 hotspot. This refutes the expected
  benefit of these two heuristics at equal mass. Linear-surrogate screening had predicted the
  response-aware family would backfire by +112% on a different coil, so transport has already
  demonstrated that this extrapolated magnitude and migration do not generalize. A follow-up
  LP audit makes the minimax failure sharper: at the common mass, adding `R_pred≥0` constraints
  makes the minimax surrogate infeasible even with a 20 cm per-sector cap, although the mass
  budget itself is feasible. Direct transport then falsified the extrapolation decisively: the
  stress-test design worsened the peak by **+891% ± 36% (9.91×; z=+88)**, again at c23. No
  nonuniform candidate improves uniform at this mass. The result supports bounded, sequential
  direct-forward relinearization or nonlinear design data—not one-shot linear extrapolation.

**Status & caveats (project-wide, from the standing integrity rule).** This arc has
repeatedly found that a "surprising" result (anti-correlation, huge multiplier,
byte-identical results across supposedly different configs) traced back to a pipeline bug
(cwd races serving the wrong DAGMC geometry to concurrent SLURM tasks — recurred twice;
a circular-torus proxy standing in for the real VMEC surface; a substring material-name
match inflating a volume). The team's discipline is to treat such results as bugs first and
either fix-and-reproduce or mark the finding BLOCKED rather than report it — several results
above (rung-2 attribution projection, the earlier "5×" absolute magnitude) are explicitly
held at that status rather than being overclaimed.

---

## What's ready to claim vs. in progress

**Ready to claim (transport-verified or independently cross-checked):**
- StellaratorSource runs real QUASR QH/QA fluxmap-driven transport in OpenMC (Arc 1).
- SPF polarized birth-direction physics and the B/C-share-shape correction are derived. The
  corner identities (2/3, 1, 2/3, 1/3 in σ₀) provide the scalar oracle, while the independent
  analytic NWL comparison remains part of the fuller tiered verification (Arc 2). An early
  conformal smoke run indicates the device-level sign (parallel raises coil flux,
  perpendicular lowers it), but not a headline-ready magnitude.
- Periodic-BC one-period modeling is *mandatory* for coil-load studies (quantified 1.03x vs
  4.2x with/without) (Arc 3).
- Attribution-guided single-coil placement beats uniform by 4.2× but migrates the peak; both
  multi-coil variants lose to single-coil (Arc 3).
- Patch-worth does not generalize past n=4; scalar/attribution worth anti-correlates with
  true finite-difference worth on at least one patch (Arc 3).
- Authoritative per-sector added-shield mass (40M MC), the 2.45× equal-thickness-mass-spread
  result, and the "every patch is lighter than uniform" refutation (Arc 3).
- The clearance-vs-response-need collision at 4.7 cm inboard (Arc 3).
- Coils are not VR-starved; direct-forward transport suffices for the coil-load metric,
  eliminating the need for VR/SSR in the main optimizer loop (Arc 3).

**In progress / explicitly not yet claimed:**
- Full tiered SPF sampler verification (χ²/KS tables, full-stack analytic NWL comparison,
  scalar-oracle ±43%/±22% reproduction) — physics spec and plan exist (`CLAUDE.md`); not yet
  executed/reported in this checkout.
- SPF composition with the real stellarator fluxmap source (as opposed to the analytic
  tokamak test problem) — a stated integration target, not yet built.
- Rung-2 exact adjoint-importance-vs-worth correlation — BLOCKED pending an exact VMEC/ParaStell
  surface projection (no circular-torus proxy allowed).
- Native FW-CADIS on the one-period DAGMC model is debugged and runs, but no candidate
  weight-window file has cleared the production-readiness bar for the SSR coil calculation;
  this is now secondary given the direct-forward finding, but the infrastructure remains for
  any future genuinely-deep quantity (e.g., total nuclear heating).
- Paper framing decisions (lead paper, Paper-1 headline choice, whether to include a
  ripple/tokamak QI arm) are open per the SPF paper-direction notes.
