# CONTEXT RECOVERY — adjoint magnet-shielding program (tonight, 2026-07-25)

Read this FIRST after compaction to continue seamlessly. It captures the live state, the plan, the
strategic direction, the just-arrived Kappel paper, and the concrete next actions. Cross-references the
other docs (all under `spf_prototype/`): `docs/EXECUTION_PLAN_magnet_adjoint.md` (the 6-phase plan +
adjustments), `docs/SCIENCE_STRATEGY.md`, `docs/DATABASES_AND_RELATIONAL_DESCRIPTORS.md`,
`docs/COIL_DESIGN_AND_FILTERING.md`, `paper/PAPER_DIRECTION.md`, `paper/METHODS_DRAFT.md`,
`paper/RESULTS_adjoint.md`. Branch `spf-prototype`, NEVER pushed. Ginsburg via `ssh -o BatchMode=yes ginsburg`.

## 0. What this program is
Overarching: geometry-first stellarator neutronics (Ethan Peterson steer, 2026-07-23). The methodological
centerpiece = an **adjoint through-shield coil-importance map**: a localized coil-response (kerma/flat)
adjoint source solved with OpenMC random-ray ADJOINT; the adjoint flux over the plasma IS the importance
map psi_dagger(r) = how much a neutron born at plasma point r drives that coil's load THROUGH ~1 m of
shield. **psi_dagger is REACTIVITY-INDEPENDENT** (geometry + coil response only) -> re-weight S(r) for free.
Contributon C(r)=S(r)*psi_dagger(r) = attribution. Feeds a closed-loop breeder<->shield placement optimizer.

## 1. VALIDATED + COMMITTED (with numbers)
- **Adjoint pipeline works on real DAGMC** (QH Landreman-Paul reactor, dagmc_corr_w20/w35f). Files:
  `shield_opt/fwcadis/adjoint_importance.py` (the run: filamentary coil-guide-curve adjoint source,
  --response flat/kerma, --scatter-order 0/1, --no-vacuum for QA), `shield_opt/adjoint_placement.py`
  (contributon + (theta,phi) placement priority; emissivity uniform/bosch_hale, bosch_hale DEFAULT),
  `shield_opt/adjoint_closed_loop.py` (adjoint-informed shield optimizer),
  `shield_opt/fwcadis/reciprocity_check.py`, `shield_opt/fwcadis/compare_p0_p1.py`.
- **RECIPROCITY validated on BOTH devices** (adjoint INT(S*psi) vs forward per-coil flux, matched geometry,
  n=coils, uniform emissivity to match the uniform-source forward StellaratorSource):
  - QH: Spearman +0.76 (raw) / +0.70 (length-norm), n=20, w35f.
  - **QA: Spearman +0.83 / +0.85, n=16** (STRONGER; QA=nfp2 precise-QA reactor, dagmc_qa_base). Cross-device
    method validation across QS type. QA maps in `~/adjoint_test/qa_coil{11..26}_flat/`.
- **Closed-loop benefit (QH coil-15, flat, SURROGATE objective):** 19% Bosch-Hale / 55% uniform /
  -1% no-attribution baseline. CRITICAL: this is a SURROGATE number (exp(-k*delta), UNCALIBRATED MFPs) and
  the coil dose SATURATES with shield thickness (leak-around, not exponential) -> the surrogate is the WRONG
  FORM and CANNOT test the thesis. Must be replaced by a DIRECT transport eval (= Step 1 below).
- **M7 angular (P1) built into the Ginsburg DAGMC build** (`~/src/GitHub/openmc-m7`, angular C++ patched +
  compiled). P1 vs P0 on the real coil: P1/P0 median ~1.006, placement phi0 UNCHANGED -> SCALAR ADJOINT
  SUFFICES for coil placement. (M6 slab: P1/P0 ->1.7x at 5 mfp deep in scattering bulk, ~3% in clean gaps.)
- **Reactivity robustness:** peak coil + hot-5 set UNCHANGED uniform->bosch_hale (free re-weight).
- Key committed hashes (branch spf-prototype): reciprocity+free-ranking `4b80fd662`; Phase A docs+plan
  `ccf88ce58`; QA bring-up `6f5db25e1`; science strategy+reprioritization `e2cd69b22`; rigor pass
  (realistic reactivity, near-void fix, P1 wiring) `851c762d7`; self-reporting emissivity `d5cf20002`.

## 2. LIVE STATE / IMMEDIATE NEXT (what was mid-flight at compaction)
- **STEP 1 KILL-SHOT (the falsification experiment) -- NEEDS A FIX then re-run.** Goal: build adjoint-PLACED
  vs equal-material-budget UNIFORM shield DAGMC on QH, run FULL-TRANSPORT coil dose, compare (placed beats
  uniform? if placed~uniform even at high concentration -> thesis FALSIFIED, still publishable). delta fields
  computed (placed concentrates 67% of shield material in top-25% cells): `~/qh_step1_deltas.npz`. Two DAGMC
  built (`dagmc_qh_placed.h5m`, `dagmc_qh_uniform.h5m`) via `build_qh_placed.py`/`build_qh_uniform.py`.
  **BUG:** those used build_qh.py = a SIMPLIFIED 5-layer radial build; `coil_run_v2.py` expects the FULL
  9-layer corrected build (multiplier/gap/thermal_shield) -> coil-dose job 9191387 aborted (material
  mismatch). FIX: rebuild placed/uniform from `pstl_test/corrected/build_corrected.py` (full radial build),
  feed the delta field to its shield/breeder thickness_matrix, then coil dose via
  `pstl_test/corrected/coil_run_v2.py MODE BATCHES PARTICLES DAGMC OUT_TAG LI6 [WW_FILE]` (StellaratorSource,
  tallies coil_fast_flux + coil_heating; deep penetration -> may need weight windows). Build env: SLURM
  `conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl` (parastell), account=astro partition=short.
- **QA cross-device metrics (Phase E/Step 3):** QA reciprocity done (+0.83). Next: compute the 4 comparable
  metrics for QH+QA (attribution concentration Gini/PR via `concentration.py`; benefit from Step 1 transport;
  angular P1/P0; reactivity robustness) -> `shield_opt/magnet_decomposition.py` (TO WRITE) -> honest n=2
  framework + comparison, standoff-labeled, HYPOTHESIS not claimed law.
- Two research agents completed tonight: DATABASES_AND_RELATIONAL_DESCRIPTORS.md, COIL_DESIGN_AND_FILTERING.md.
  A THIRD (extended science vision, based on the Kappel paper) is being launched now.

## 3. THE KAPPEL PAPER (just read -- TRANSFORMATIVE) 
Kappel, Landreman & Malhotra, PPCF 66 (2024) 025018, "The magnetic gradient scale length explains why
certain plasmas require close external magnetic coils." DOI 10.1088/1361-6587/ad1a3e. Data: Zenodo
**10.5281/zenodo.8349408**.
- **L_grad_B** (magnetic gradient scale length) = sqrt(2)*B_coils / ||grad B_coils||_F (Frobenius norm of the
  3x3 grad-B tensor), computed from B + grad B at a point on the LCFS (B_coils via VIRTUAL_CASING in
  simsopt). L*_gradB = min over the LCFS. CHEAP, LOCAL, reactivity-independent, no coil design needed.
- **L*_gradB predicts the minimum achievable plasma-coil separation (L_REGCOIL) with R^2 = 0.941**, slope
  1.585, across a database of **>40 configs** (QA, QH, QI, tokamaks, misc), ALL scaled to ARIES-CS
  (a=1.704 m, B=5.865 T). Insensitive to the exact K/Brms constraints.
- KEY TREND: **larger plasma-coil separation is possible with SMALL number of field periods** (nfp 1-2 ->
  large separation -> reactor-attractive). L*_gradB varies by ~10x across the database.
- L_gradB is smallest on the INSIDE of the bean cross-section (where coil current density K is largest).
- Appendix D = the full table (~40 configs: NCSX, W7-X, LHD, ARIES-CS, HSX, precise QA/QH, Wechsung QA,
  Goodman QI nfp1/2/3, Jorge QI nfp1, tokamaks TJ-II/ITER-like, CTH, ATF, ESTELL, CFQS...) with nfp, beta,
  L_REGCOIL, L*_gradB. THIS IS the ~40-device engineering-relevant curated dataset -- on Zenodo, reactor-
  scaled, spanning QA/QH/QI+tokamak. Coils via REGCOIL (Landreman 2017; simsopt has the pieces -- simsopt
  1.10.6 is ALREADY installed on the base env per COIL_DESIGN doc).
- Kappel proposes L_gradB as a DIFFERENTIABLE stage-1 objective to MAXIMIZE plasma-coil separation.

### Why it's transformative for us (the mechanistic chain is now complete + grounded)
- **L_gradB (plasma field geometry, cheap) -> min plasma-coil separation [Kappel R^2=0.94] -> available
  shield thickness -> magnet neutron load / shieldability [OUR adjoint].** The whole geometry->magnet-
  protectability law now has a VALIDATED mechanistic upstream variable instead of a confounded 2-point corr.
- **Dissolves the standoff confound properly:** standoff is not a lurking variable, it is the MEDIATOR; L_gradB
  is the root geometric CAUSE. The chain is causal and grounded, not spurious.
- **Kappel stops at GEOMETRIC clearance (can the blanket fit?). We add the NEUTRONICS (does the blanket
  actually protect the magnets?).** L_gradB is a proxy for "room for shielding"; our adjoint gives the actual
  magnet load through the achievable shield. We are the neutronics extension of an adopted optimization proxy.
- **Unified thesis:** geometry sets the STRUCTURE of loading -- plasma SHAPE (elongation) -> first-wall
  concentration (xi=0.96, committed); field-gradient scale length (L_gradB) -> coil standoff -> magnet load.
  Different descriptors, one meta-claim. L_gradB is the magnet-side analog of elongation.
- **The Zenodo-8349408 database is the engineering-relevant curated set** (n~40, QA/QH/QI+tokamak, reactor-
  scaled, L_gradB + L_REGCOIL pre-computed). Pull it -> real n for the predictive law, coils via REGCOIL.

## 4. STRATEGIC DIRECTION (the framing to build toward)
- The field has NO engineering-relevant QUASR-across-all-symmetries because (a) two-stage optimization death
  loop, (b) multi-physics software silos, (c) filament vs finite-build errors. It is shifting to SINGLE-STAGE
  differentiable optimization (simsopt/DESC) + engineering proxies + AI surrogates (StellFoundry/STELLAR-AI).
- **CRITICAL GAP we fill: every proposed engineering proxy is GEOMETRIC/STRUCTURAL (clearance, curvature,
  forces) -- NONE is NEUTRONIC.** The 1.3 m clearance is a stand-in for "magnets protected" that is never
  verified. Our adjoint magnet-shieldability IS the missing neutronics engineering proxy.
- **User-endorsed reframe:** the ParaStell->DAGMC->adjoint AUTOMATION goes from "glue" to **"the neutronics
  evaluator for the single-stage / generative-design future."** L_gradB is already used/adopted; we provide
  the neutronics it is a proxy for. Honest limit: making the adjoint a DIFFERENTIABLE objective term needs a
  shape gradient (differentiate the transport) -- we are a fast EVALUATOR/filter, not (yet) an autodiff term.
- **Layered data strategy:** (1) engineering-relevant backbone READY NOW = Kappel/Zenodo ~40 + simsopt-bundled
  (W7-X/NCSX/HSX) + our QH/QA; (2) QUASR filtered (aspect + full-radial-build ~1 m standoff cut -> ~15-30%
  survive) for the observational free-streaming coil law; (3) simsopt/FOCUS single-stage at CONTROLLED
  standoffs = the confound-breaker (same plasma, N standoffs); (4) QI from ConStellaration boundaries + coils.
- Compute: DAGMC build (~30-60 min, license-free cad_to_dagmc) is the bottleneck, not coil design (minutes).
  n~20-50 reachable overnight via SLURM array (~25 concurrent). MVP = method + Step-1 transport benefit +
  zoo/relational law + honest n=2..few comparison. Ambitious (gated) = the L_gradB->standoff->magnet law
  across the Zenodo-40 + QI.
- Prior-art to pre-empt: ParaStell WISTELL-D 3-D shield/TBR (differentiate: adjoint says WHERE); coil-non-
  planarity<->plasma-surface arXiv:2604.26763 (stay on NEUTRONICS axis); HELIAS ADVANTG adjoint VR
  (differentiate on ATTRIBUTION). Kappel itself: cite as the standoff-predictor; we add the neutronic outcome.

## 4b. EXTENDED-VISION AGENT VERDICT (docs/SCIENCE_VISION_extended.md) -- read it, key points:
- **The chain is TWO links, asymmetric.** Link 1 (L_gradB -> plasma-coil separation) is ALREADY VALIDATED
  by Kappel (n>40, R^2=0.94) AND **our QH+QA are named points in his Appendix D**: precise QH L_REGCOIL~1.52 m
  / L*_gradB~3.19 m; precise QA ~2.87 m / ~5.30 m; measured as-built standoffs 1.63/3.11 m (ratio 1.91) sit on
  his predicted separations. So "geometry->standoff" is NO LONGER an n=2 line we must defend -- it's a
  published 40-device law with QH/QA as named points (table lookup). Link 2 (separation->shield room->magnet
  load) is OURS to establish, still n-limited to the DAGMC+adjoint subset.
- **BIGGEST RISK (test it first):** Link 2 is trivial -- magnet load tracks standoff by plain 1/r^2 geometric
  dilution, so the adjoint earns nothing over Kappel + inverse-square; OR leak-around saturation makes shield
  room past ~1 m irrelevant (available-shield-thickness is the wrong operative variable). Either collapses the
  chain to geometry alone (= the Step-1 null in a new guise).
- **SHARPEST NEW EXPERIMENT: the controlled-standoff test at fixed plasma / fixed L_gradB.** Hold one boundary,
  simsopt `CurveSurfaceDistance` to build coils at ~3 standoffs, measure whether magnet load AND attribution
  concentration track standoff; decompose the response into geometric dilution (1/r^2) vs shield attenuation vs
  leak-around. THE test of whether the neutronics adds signal beyond 1/r^2. Cheap variant: spatial-coincidence
  -- does the adjoint magnet-load hotspot land at Kappel's L_gradB minimum (bean-inside)?
- **PAPER-1 scope (reframe of the existing MVP):** label the QA-vs-QH comparison BY standoff/L_gradB, cite
  Kappel, note both devices are in his DB -> the "confound" dissolves into a named, root-caused axis. Deliver:
  adjoint attribution + reciprocity (+0.76/+0.83); the Step-1 placed-vs-uniform transport kill-shot (confirm OR
  null, both publishable); the free-streaming zoo coil-concentration law (real-n companion); the
  L_gradB->standoff->magnet chain stated as the organizing HYPOTHESIS with P2 as its test; SPF supplementary.
  NO through-shield predictive law claimed in P1 (that's P2 = the Zenodo-40 study).
- HONESTY FLAG: QH/QA preserve Kappel's ordering + the ~1.9x gap but sit somewhat BELOW his regression (L_gradB
  ratio 1.66 under-predicts the measured 1.91) -> quote the ordering + "large gap", not a spuriously exact ratio.

## 5. KEY PATHS / ENVS / JOBS
- Ginsburg openmc builds: adjoint runs use `~/src/GitHub/openmc-spf` (PYTHONPATH=$SRC + PATH=$SRC/build/bin +
  LD_PRELOAD=$SRC/build/lib/libopenmc.so, env spf-stellarator for deps, OPENMC_CROSS_SECTIONS=
  ~/openmc_data/endfb-viii.0-hdf5/cross_sections.xml). Angular P1 build: `~/src/GitHub/openmc-m7`.
  find_cell/discovery: plain conda openmc (no LD_PRELOAD). ParaStell builds: SLURM, env
  `/ginsburg/astro/users/tjk2147/spf_work/envs/pstl`, account=astro partition=short.
- QA: dagmc `~/pstl_test/qa/dagmc_qa_base.h5m`, coils `~/pstl_test/qa/Wechsung_QA24.coils` (x100 m->cm),
  centroids `~/qa_coil_centroids.npz` (cells 11-26, filament find_cell-verified mapping), fluxmap
  `~/qa_reactor_fluxmap.npz` (nfp2), forward `~/pstl_test/corrected/percoil_unpol_qa_base.npz`. QA peak=cell23.
- QH: dagmc `~/pstl_test/corrected/dagmc_corr_w20.h5m` / `_w35f.h5m`, coils `~/pstl_test/coils_qh` (x100),
  centroids `~/pstl_test/corrected/coil_centroids_corr.npz` (cells 11-30), fluxmap `~/qh_freestream_fluxmap.npz`,
  forward `~/pstl_test/corrected/percoil_unpol_corr_w35f.npz`. QH adjoint-peak=cell20 (forward-hottest=15,
  an outlier). Committed maps: `spf_prototype/data/adjoint/coil15_fil_flat.npz`, `m7/coil15_P0.npz`,`P1.npz`.
- Multi-thickness (saturation evidence): percoil_unpol_corr_{w20,w35f,w60,w90} (shield ~flat 7e-3 w35->w90);
  breed_{d5,d15} (breeder scan). NO clean exponential -> Phase C must be direct transport, not lambda fit.
- Reciprocity length-norm + engineering filters computable from filament polylines (`_parse_coils`).

## 6. IMMEDIATE NEXT ACTIONS (resume order)
1. **FIX + RE-RUN STEP 1 (kill-shot)** on the FULL radial build: rebuild placed/uniform from
   build_corrected.py, coil dose via coil_run_v2.py, compare (weight windows if stats thin). THE result that
   gates the thesis.
2. **Phase E metrics:** write `shield_opt/magnet_decomposition.py`, compute 4 metrics for QH+QA, honest n=2.
3. **Pull the Kappel/Zenodo-8349408 dataset** (~40 reactor-scaled configs w/ L_gradB + L_REGCOIL); compute
   L_gradB ourselves (from VMEC B+gradB, or use their table) -> the observational geometry(L_gradB)->standoff
   spine; then the extension standoff->magnet-load via our free-streaming coil concentration (cheap) and
   adjoint (DAGMC subset).
4. **Step 2 zoo law:** free-streaming per-coil concentration vs relational descriptors (standoff, non-
   planarity, gap-streaming) + L_gradB, over the filtered set. Real n, dissolves confound.
5. Stand up simsopt stage-2 controlled-standoff coil design (~1 day) for the confound-breaker + QI.
6. The extended-science-vision agent's output -> `docs/SCIENCE_VISION_extended.md` (folds Kappel + neutronics-
   evaluator framing into the roadmap). Read it and reconcile with EXECUTION_PLAN.

## 7. STYLE / CONSTRAINTS (persistent)
Hard-gate on correctness; strict phase ordering; honest non-cherry-picked diagnostics; no tuning-to-target;
n=2 is n=2 (framework+comparison+hypothesis, not a law); label surrogate-vs-transport; every result stamps
device/geometry/emissivity/response/objective/relerr/n. Branch spf-prototype, never push. Commit messages end
`Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. Collaboration circle: Ethan Peterson, E. Paul,
Sophia (MIT), Jack Fletcher. Confidential PDF plasma_source_sampling_paper.pdf must not be distributed.
