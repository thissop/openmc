# SPF documentation — master index (`key.md`)

The single map of every SPF doc/paper/deck. Paths are relative to `spf_prototype/`
(this file lives at `spf_prototype/docs/key.md`); a few coupled files live one level
up at the repo root and are shown with `../`.

## Conventions for future agents (read before adding a doc)

1. **Every new doc must be filed under the right subfolder AND get a bullet here** in
   `key.md` with a one-line summary. A doc that is not in this index does not exist.
2. **Historical notes are condensed into the `RETROSPECTIVE_*.md` files** in
   `docs/notes/`. If you write a new results/planning/handoff note that later becomes
   history, fold its substance (every number, verbatim key passages as attributed
   blockquotes) into the matching retrospective and `git mv` the original into the
   right `docs/notes/<sub>/` folder — never delete it.
3. **Active design docs and operational runbooks stay in place** (e.g. the `sweep/`
   study docs incl. its self-contained Ginsburg runbook `sweep/docs/README.md`). Only
   *index* them here; do not move them.
4. **Taxonomy:** `docs/notes/{orientation,results,handoffs,planning}/` = condensed
   historical notes; `docs/reference/` = standing reference (parameters, scope,
   prior work); `docs/papers/` = external reference PDFs; `docs/decks/` = figure/slide
   decks; `docs/teaching/` = teaching deck; `docs/archive/` = superseded (kept, never
   deleted; incl. `legacy_precise_qa/`). The historical smoke-run deck is coupled to
   `../figures.py` and stays at the repo root — see below.
5. **Never move** any `.py`, `src/`, `CLAUDE.md`, `AGENTS.md`, or a file coupled to a
   `.py` build (the repo-root smoke-run deck). Non-destructive only.

---

## docs/ — papers (ACTIVE; left in place)

- `analytic upgraded.tex` — flagship analytic manuscript (Kiker & Peterson): a
  field-period-perturbation analytic neutron wall load for spin-polarized plasmas in
  quasisymmetric stellarators, reduced to Legendre F/E elliptic integrals, JAX-
  differentiable, validated to ~1% vs quadrature.
- `Analytic_Neutron_Wall_Loading_for_Spin_Polarized_Plasmas_in_Quasisymmetric_Stellarators_by_Field_Period_Perturbation.pdf`
  — compiled 21-page PDF of the analytic paper (readable without pdflatex).
- `mc_spf_nwl.tex` — companion MC manuscript: the OpenMC spin-polarized fusion neutron
  *source* (exact stratified inverse-CDF sampler + Gram-Schmidt rotation), cross-
  validated against the analytic NWL in the free-streaming limit (≤1.8%, 4–6% inboard
  residual = analytic point-patch limit).
- `mc_spf_nwl.pdf` — compiled 10-page PDF of the MC paper.

## docs/notes/ — retrospective syntheses (condensed history)

- `RETROSPECTIVE_tiers.md` — the Tier 1→8 arc: sampler verification → analytic NWL →
  OpenMC recovery → scattering/TBR reactor tiers → pitched/3-D field source; combines
  all `RESULTS_*.md` + `JOINT_VALIDATION.md`, numbers preserved.
- `RETROSPECTIVE_stellarator_and_smoke_run.md` — the stellarator extension (toy QA-low
  and real QUASR 59509 free-streaming MC↔analytic checks) and the precise_QA conformal
  smoke-run deck (η≈0.63 survival, the coil-flux sign convention + typo warning).
- `RETROSPECTIVE_orientation_planning_handoffs.md` — the project scaffolding: §0
  API/env ground-truth, the verification-first roadmaps, the compaction-surviving
  state dumps, the Bae-2025 positioning, and the still-live data/scope docs.

### docs/notes/orientation/

- `NOTES_orientation.md` — §0 pre-implementation ground truth for this checkout:
  verified `Source`/`SourceSite`/`ParticleType`/`prn` signatures, the custom-source
  CMake build recipe, the bare-environment finding, and the corner-cross-section
  re-derivation.

### docs/notes/results/

- `RESULTS_tier1.md` — Tier-1 direction sampler (no transport): sympy identities,
  C++↔Python parity 6.1e-16, per-mode cosθ χ²/dof≈1, rotation isotropy/steering, KS.
- `RESULTS_tier2a.md` — Tier-2a independent analytic NWL: three engines agree to
  1.8e-15; reproduces Schwartz Fig.2 oracles; corrects spec's "iso≡C" to **B≡C**.
- `RESULTS_tier2b.md` — Tier-2b full-stack OpenMC vs analytic: free-streaming
  directionality residual mean~0 (σ~0.68), oracles reproduced, B≡C / linearity.
- `RESULTS_tier3.md` — Tier-3 presentation artifacts (postprocessing): peaking factor
  RISES with polarization (iso 1.33 → A 1.65, B/C 1.61); B/C center-stack relief.
- `RESULTS_tier5.md` — Tier-5 scattering + tungsten first wall: thin-wall recovers
  free-streaming, real scattering dilutes steering ~3×; heating split (FLiBe 91%).
- `RESULTS_tier6.md` — Tier-6 FLiBe TBR: 1.149 at 30% Li-6 (ARC/LIBRA band),
  modestly polarization-dependent (~0.6%); enrichment peaks ~50%.
- `RESULTS_tier7.md` — Tier-7 absolute rate (postprocessing): η restored gives A
  +50% power / C −50%; B = free steering; TBR-per-neutron ~mode-independent.
- `RESULTS_tier8_angled.md` — Tier-8/Part-A pitched-field rigor anchor: quad matches
  anarrima angled kernels to 3.6e-14; β=0 recovers toroidal bit-identically.
- `RESULTS_tier8_3d.md` — Tier-8/Part-B end-to-end 3-D-field stellarator run
  (placeholder Helios rig): per-mode TBR ~1.37, first-wall steering, η=0.89.
- `RESULTS_quasr.md` — free-streaming MC↔analytic on real QUASR QA 59509 (ε_eff~0.43):
  ≤1.8% on 3 walls, bounded 4–6% inboard = point-patch limit (shaping² scaling).
- `RESULTS_toy_qalow.md` — free-streaming MC↔analytic on a toy low-ε_eff (~0.08) QA
  source: worst-wall discrepancy 1.8%; pure-numpy analytic engine (==anarrima to 3e-14).
- `JOINT_VALIDATION.md` — the analytic↔MC cross-check protocol: shared precise_QA
  equilibrium, audited φ/polarization conventions, and the "not yet numbers-level" caveat.

### docs/notes/handoffs/

- `CONTEXT_REPORT.md` — read-only `file:line`-cited snapshot for a collaborator scoping
  the 3-D phase: repo map, sampler internals, validation status, gaps, reusable assets.
- `HANDOFF.md` — compaction-surviving authoritative state dump (env, commits, physics/
  API facts, per-tier numbers, Bae positioning) as of the Bae-comparison session.

### docs/notes/planning/

- `DATA_NEEDED.md` — **still live**: ranked real-Helios inputs to replace public
  placeholders, each tagged to its `INJECT(helios)` code seam.
- `NEXT_STEPS.md` — three ranked directions beyond free-streaming (angled field →
  scattering → blanket/TBR), cheap-extension vs new-tier.
- `PLAN_next_phase.md` — reviewed proposal for the scattering/blanket phase (Tiers 4–6),
  organized around "every new physics step ships a replacement verification".

## docs/reference/ — standing reference

- `BLANKET_PARAMS.md` — fully-cited ARC/LIBRA-class parameter table (FLiBe breeder,
  Li-6 enrichment, W/steel/Be stack) grounding the reactor model.
- `LIMITATIONS.md` — the core prototype's scope boundaries (scattering off = geometric
  NWL, toroidal/constant field, single energy, convex square cross-section).
- `DEFERRED.md` — **still live**: what the 3-D-field tier does NOT support (axisymmetric
  transport geometry, synthetic field map, circular plasma), each with its code seam.
- `PRIOR_WORK_BAE2025.md` — full-text comparison to Bae et al. 2025: airtight mode
  mapping, the (a,b,c) symbol-collision warning, and the honest contribution statement.

## docs/archive/ — superseded (kept, never deleted)

- `README.md` — the old front-door README (physics summary, layout, reproduce commands)
  covering only Tiers 1–2b; stale as a front door. **A fresh top-level README should be
  written** (see report); this copy is preserved for its reproduce recipes.

## docs/teaching/ — teaching deck (ACTIVE; authored by another agent)

- `spf_order_parameter_teaching_deck.tex` / `.pdf` — from-scratch teaching Beamer deck
  on the SPF field-coherence order parameter (the `sweep/` study's C). Indexed only.

## docs/papers/ — external reference literature (gitignored: large/copyright, kept local)

- `Analytic neutron wall loading from spin-polarized fusion in axisymmetric geometries.pdf`
  — Schwartz 2025 (arXiv:2507.11758): the axisymmetric NWL analytic our field-period
  work expands around; C=1 is the tokamak-limit validation anchor.
- `Neutronics analysis of spin-polarized fuel in spherical tokamaks.pdf` — Bae et al. 2025
  (Nucl. Fusion 65 086051): the closest prior MC SPF neutronics (spherical tokamak).

## docs/decks/ — figure / slide decks

- `peterson_figures.tex` / `.pdf` — short figure summary (7 figs + captions) of the SPF
  source verification + blanket neutronics, prepared for E. Peterson (.pdf gitignored).

## docs/archive/ — superseded / bulky artifacts (kept, never deleted)

- `legacy_precise_qa/` — the superseded single-equilibrium precise_QA Ginsburg workflow
  (`SETUP_GINSBURG.md`, `RUN_ON_GINSBURG.md`, `GINSBURG_CLAUDE_HANDOFF.md`,
  `ginsburg_job.sh`), replaced by the QUASR sweep (`sweep/docs/README.md`,
  `sweep/sweep.sbatch`); see its README. The driving CODE (`run_ginsburg.py`,
  `run_conformal.py`) stays in `python/` — the sweep reuses it.
- `tier2b_statepoints.tgz` (gitignored) — saved OpenMC statepoints from the Tier-2b
  free-streaming-vs-analytic cross-check; kept locally for re-analysis.
- `README.md` — the prototype's original root README, superseded by this index.

---

## Left in place elsewhere (indexed, NOT moved)

### spf_prototype/ root — active + operational

- `slides/spf_talk.tex` (ACTIVE) — from-scratch teaching Beamer deck (SPF physics,
  (a,b,c) bookkeeping, Schwartz analytic, OpenMC source, Tiers 1–7, honest findings).
- `ginsburg_preflight.sh` (OPERATIONAL) — read-only Ginsburg env checks (account,
  storage, cross sections, SLURM partitions); workflow-agnostic. The Ginsburg runbook is
  now the self-contained `sweep/docs/README.md`; the legacy precise_QA Ginsburg docs are
  archived under `docs/archive/legacy_precise_qa/`.

### spf_prototype/sweep/ — the QS-stellarator coherence sweep (ACTIVE/OPERATIONAL)

- `sweep/NOTES_orientation.md` (ACTIVE) — ground-truth inventory of reusable modules in
  this checkout + the `spf_fieldmap_v1` interop contract + local-cylindrical-frame C.
- `sweep/docs/THEORY.md` (ACTIVE) — the physical rationale: field-direction coherence C
  as the order parameter, the direction tensor T, the universality/factorization hypotheses.
- `sweep/docs/FACTORIZATION.md` (ACTIVE) — the η(C,blanket) ≈ η_source(C)·A(τ_scatter)
  separability hypothesis and its three break modes.
- `sweep/docs/EXPERIMENTAL_DESIGN.md` (ACTIVE) — the fixed-blanket control, C-first
  config selection, the two observables, and the two-paper plan.
- `sweep/docs/SCOPE_AND_CAVEATS.md` (ACTIVE) — honest calibration: pymoab-free DAGMC as
  throughput contribution, uniform-offset blanket, coil-fit selection effect.
- `sweep/docs/README.md` (OPERATIONAL) — the exact Ginsburg runbook for the sweep
  (rsync, env preflight, single-config bring-up, SLURM array, pull-back, analysis).

### repo root (../) — the precise_QA smoke-run deck (HISTORICAL; left in place)

Coupled to `../figures.py` (a `.py` — must not move) + `../figs/` + `../spf_export/` via
repo-root-relative paths; moving them would break the documented build. Narrative folded
into `docs/notes/RETROSPECTIVE_stellarator_and_smoke_run.md`.

- `../README_spf_deck.md` — build/reproduce hub for the smoke-run deck package; lists
  headline numbers, the honest-framing caveats, and the "parallel reduces coil flux" typo flag.
- `../talk_script.md` — first-person speaker notes (one block per slide) for the deck.
- `../qa.md` — 11 anticipated expert Q&A (free-streaming TBR~0, sampler verification,
  DAGMC watertightness, VR plan, the coil-flux sign convention).
- `../slides/talk.tex` — the 12-frame Beamer deck ("Does SPF neutron steering survive a
  real stellarator wall?"), headline η = 0.63 ± 0.22.
