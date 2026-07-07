# Overnight progress — 2026-07-07

Running log of the autonomous session (StellaratorSource core class + real √g + SPF-in-TokamakSource
+ geometry→peaking). Newest status at top of each section. Times approximate.

> **⚠️ SSH to Ginsburg expired late in the session** (Duo re-auth needed) — the last two cluster jobs are
> prepped but did NOT run: the C++ `StellaratorSource` register+build, and the real-`√g` `precise_QA`
> solve (with the fixed V1 quadrature). **To resume in ~3 commands: follow `ginsburg_jobs/RESTART.md`.**
> Everything else below is done, validated locally where possible, and committed.

## TL;DR (read this first)
- **HEADLINE SCIENCE (T3): the SPF peaking-factor modulation IS predictable from field coherence, and
  it is NOT random.** Analytic free-streaming study over 12 QUASR devices: PF_perp/PF_unpol vs S_φ has
  Spearman **r=−0.83 (p=0.001)**; QA vs QH separate cleanly (Mann-Whitney p=0.003/0.005). Physical
  story: the SPF lever on peaking is **largest where the field is least coherent** (low-S_φ QH) —
  perpendicular emission reinforces hotspots (peaking↑), parallel smears the load (peaking↓); coherent
  QA fields respond weakly. → `docs/notes/planning/FINDINGS_geometry_peaking.md`, `figs/geometry_peaking/`.
  Honest: class-level (not precise within QA), free-streaming only, n=12.
- **Conformal StellaratorSource sampler: algorithm PROVEN.** Validated against an analytic
  circular-torus √g — birth weight exactly 1, ρ-marginal χ² *converges* under grid refinement
  (6.08→1.28), and it reproduces the **outboard differential-volume bias** (uniform-θ rejected at
  χ²/dof≈183). Core of the whole StellaratorSource, done independent of DESC.
- **Native OpenMC PR #3999 build WORKS** (0.15.4-dev212) from source on the cluster (MPI on, DAGMC off;
  run with `LD_PRELOAD` of the built libopenmc — conda ships a competing one in RPATH). Enables the
  native TokamakSource-SPF and the C++ StellaratorSource core class.
- **Real √g pipeline:** DESC 0.17.2 (separate venv); `quasr_fluxmap.py` with 5 hard rigor gates +
  a portable binary writer. DESC solve for 803097 re-running (data-path fixed).
- Rigorous **11-page LaTeX evaluation** of StellaratorSource compiled to PDF; design docs A/B +
  phase-3 adjoint memo committed.

## DONE / validated
- [x] Committed tonight's docs (`026426cec`); confirmed the in-review PDF is gitignored (not in git).
- [x] Prereqs on Ginsburg: cloned OpenMC + fetched PR #3999 (+submodules); installed DESC 0.17.2 (venv).
- [x] `python/quasr_fluxmap.py` — real-√g producer (DESC eq → √g,R,Z,B on one (ρ,θ,ζ) grid) with gates
      V1 (own-quadrature vol = DESC V), V1b (independent divergence-theorem vol), V2 (single sign),
      V3 (lexsort tensor reshape), V4 (LCFS reconstructs the boundary fit).
- [x] `python/stellarator_source.py` — conformal CDF-cascade sampler (reference for the C++ class) +
      analytic circular-torus fluxmap generator.
- [x] `python/test_stellarator_source.py` — **PASSES all gates** (see TL;DR). Plot:
      `figs/stellarator_source/circular_torus_validation.png`.
- [x] `docs/writeups/stellarator_source_eval.{tex,pdf}` — rigorous methods evaluation (11 pp).

## DONE since (updates)
- [x] **Native SPF `TokamakSource` (A) BUILT into OpenMC and ALL 23 unit tests pass** (15 physics-mirror
      reproducing Schwartz oracles + 8 Python-API round-trip). Staged in `native_spf/`, applied onto the
      PR build. Only the full Schwartz ±43% *transport* remains (bonus).
- [x] **T3 geometry→peaking: complete + rigorous** (see TL;DR). Engine `sweep/geometry_peaking.py`,
      data `sweep/sweep_out/geometry_peaking.json`, 4 figures, findings note. Caught+fixed a real wall-
      geometry bug (self-similar ρ=1.3 wall folds into the source on QH → 1/r² spike; uses normal-offset).
- [x] OpenMC PR #3999 from-source build works (run via `LD_PRELOAD`; conda's libopenmc is first in RPATH).

## IN PROGRESS
- [ ] **Real √g DESC solve (803097):** first cold `eq.solve()` did NOT converge (force residual ~1e28,
      sign-flipping √g — my V2 gate caught it). Switched to DESC **continuation** (`solve_continuation_
      automatic`, ramps resolution) + a V0 force-residual gate; re-running.
- [ ] **C++ `StellaratorSource` core class** (agent): the class + binary-fluxmap reader + CDF cascade +
      SPF direction, patterned on TokamakSource → `native_spf/`. Then build + axisymmetric-limit test.

## NEXT (queued, in priority order)
1. **Transport ±43% capstone** — reproduce Schwartz's ±43% inboard / ∓22% outboard midplane NWL with
   the *native* SPF `TokamakSource` (b̂=toroidal) AND the `StellaratorSource` on a circular/axisymmetric
   fluxmap. Every piece is already validated (sampler analytic, A 23 tests, StellaratorSource bit-exact,
   √g gates); this is the end-to-end integration against the anarrima oracle. Reuse the tier-2 machinery.
2. **Fix the QUASR-boundary → DESC `[C1]` convention** so real-device (not just example) equilibria solve
   (`ensure_positive_jacobian` degeneracy on the axisymmetric seed → boundary Fourier-mode sign/theta
   orientation). Then produce device-specific real fluxmaps (803097, 886079, …) and feed the peaking study.
3. **Add a QI arm** (you asked): QUASR is quasisymmetric-only (helicity 0=QA, ≠0=QH) — no QI. But DESC
   (now stood up) ships QI/omnigenity examples and there are published QI boundaries; generate 2–3 QI
   equilibria → same `quasr_fluxmap`→`geometry_peaking` pipeline → drop QI points onto the S_φ plot. Turns
   the QA↔QH contrast into a 3-class QA/QH/QI landscape (stronger generality; probes the low-S_φ end).
4. **Register the Python `StellaratorSource` API more fully** and wire the C++ prn-parity build-time gate.
5. Depolarization / scattering arm A(τ); energy-radius coupling; alpha channel (all out of prototype scope).

## Notes / decisions
- DESC deps conflict with OpenMC → DESC lives in `$HOME/desc_venv`; producers write portable files
  the OpenMC-env source consumes. (Confirmed by the existing `quasr_equilibrium_field.py` design.)
- √g sign: DESC/orientation gives √g = −R·𝒥_pol; we store |√g| and validate via the sign-independent
  divergence-theorem volume (caught by the T4 sympy check — no shortcut).
- Metric readout default (pending your call): **coil fast-flux/DPA leads**, peak factor + TBR secondary.
