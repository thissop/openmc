# analysis_status — 2026-08-24

## Action: reconstructed full finite-difference response matrix from existing patch study
- Inputs: cluster_pull/coil_step1b_uniform.npz (baseline) + 12 coil_patch_t*p*.npz + delta_qh_patch_*.npz
- Script: reconstruct_b_fd.py (rerunnable, repo-relative). Outputs: patch_manifest.json, *.csv, results.md, 3 figs.

## Numerical result (with uncertainty caveats)
- 98% of dR entries > 2sigma (independent-error bound) -> effects are RESOLVED, not noise.
- SVD (whitened E, 12x20): effective rank 4.66; 2 modes = 91% energy; top mode common-sign 0.75.
- 5 patches REDUCE the peak (inboard sectors: t3p3,t0p2,t0p3,t3p2,t1p3), 6 TRADEOFF, 1 BACKFIRE (t2p3, +41% on coil 20).
- Best single patch t3p3: peak 5.1x lower (BASELINE-CONFOUNDED number).

## OUTCOME = #1 HIDDEN STRUCTURED CONTROL (placement-dependent), scoped
- ROBUST: placement structure — inboard sectors reduce the peak & the whole inboard cluster; wrong sectors backfire.
- NOT ROBUST: the absolute 5x. Confounds: (a) patch shield-outside 40 vs baseline 42.55 (~2.55cm); (b) patches used
  FW-CADIS WW, baseline did NOT -> correlated-seed cancellation likely broken; (c) cell-integrated flux != density.

## Interpretation
Reproduces/generalizes the kill-shot (adjoint-placed = 4.2x; best inboard sector patch = 5x). The prior "flat worth"
was the INTENSIVE per-added-cell worth (~uniform); the EXTENSIVE peak reduction is large and placement-structured.
Superposition of patches NOT assumed (disjoint single sectors, built one-at-a-time).

## Next gate: A (structure is resolved) -> clean-baseline nonlinear validation before any absolute claim.

## 2026-08-24 (cont) — scope + Gate A launched
- SCOPE: this 12-patch B_FD is the free launch point, NOT the paper. Full program = source-attribution A vs
  intervention B vs feasible controllability -> response-aware 3D blanket/shield co-design -> Pareto fronts
  (TBR/heating/mass/build) vs uniform/clearance/wall-loading/attribution/worth baselines -> SPF coupling -> 2nd config.
- SPATIAL-RESOLUTION CAVEAT (critical): the 12-patch matrix is only the projection of the continuous design->response
  operator onto a COARSE 4x4 basis. A flat/coarse B must NOT be generalized to the continuous field. Next: build a
  hierarchical geometry-aware basis (uniform + broad inboard/outboard + low-order stellarator-symmetric Fourier in
  theta,phi + overlapping periodic B-splines; ~30-60 modes), test multiple spatial widths, set resolution from the
  sensitivity correlation length (not patch count). Twisted field-period-compatible modes may reveal control the
  coarse patches average away.
- GATE A launched: job 9471397 reruns baseline(step1_uniform)/t3p3(reducer)/t2p3(backfire) NO-WW, same seed,
  full-torus -> clean correlated dR removing the WW-vs-noWW + baseline confound. Confirms whether structured
  placement (inboard reduces, wrong sector backfires) survives a consistent-settings comparison.

## 2026-08-24 (correction) — the B_FD is already clean; my WW-confound caveat was WRONG
- slurm_patch_worth.sh: WW="" (NO weight windows), 12x4M, seed=1 -- IDENTICAL settings to coil_step1b_uniform.
  So patches & baseline are same-settings, correlated. The "patches used WW / baseline no-WW" confound DOES NOT EXIST.
- => existing dR is a clean apples-to-apples finite difference. Only residual: ~2.55cm uniform shield offset
  (patch-outside 40 vs baseline 42.55), which cancels in patch-to-patch structure.
- => OUTCOME #1 (structured placement-dependent control, ~5x on inboard sectors, backfire on wrong sector) STANDS,
  clean, from existing data. No clean-baseline rerun needed. Cleanval job 9471397 cancelled as redundant.
- Practical note: no-WW full-torus coil run ~16h (not 3h). One node hit a transient XS-mount error (not systematic).

## RUNG 2 (attribution/importance vs worth) — INCONCLUSIVE (projection unreliable)
- adjoint_coil20_flat_P0.npz = coil-20 ADJOINT IMPORTANCE psi_dagger on 30x30x20 cartesian mesh (flat response),
  NOT contributon or exact worth.
- Projected onto 12 patch sectors via CIRCULAR-TORUS approx -> ANTI-correlates with worth (Spearman dR20 -0.62,
  dpeak -0.78), 0/3 top-3 overlap. This is physically implausible (best reducers t0/t3 = phi=0/90 seam = coil-20's
  toroidal position, should be HIGH importance) -> the projection is misaligned, NOT a real A!=B result.
- OUTCOME: mapping too ambiguous. Do NOT claim A!=B or A=B. A defensible test needs projection onto the true VMEC
  shield surface + verified coil-20 toroidal angle. (Consistent with prior scalar W_FD-vs-A n.s.: UNRESOLVED.)
- Artifacts: rung2_attribution_vs_worth.py/.json, fig_rung2_attribution_vs_worth.png.

## RUNG 3 (validate the large coarse effect) — design
Minimal cases (independent seed; identical no-WW settings; full 20-coil vectors; full-torus):
  (a) baseline (step1_uniform), (b) t3p3 reducer, (c) t2p3 backfire  -> reproducibility of the 5x/+41% with NEW seed;
  (d) half-amplitude t3p3 (+10cm) -> linearity;  (e) EQUAL-MASS uniform shield addition (+~2cm everywhere, same
  added mass as t3p3's 500 cell-cm) -> tests placement-vs-bulk (the decisive engineering control).
Cost: no-WW 48M ~16h; use 12M (~4h, resolves a 5x effect) on short. (d),(e) need new ParaStell builds.
Prereq: confirm XS path stable (one node hit a transient mount error).

## RUNG 3 — minimal validation SUBMITTED (job 9472669)
- coil_run_v3_seed.py (copy of coil_run_v3, adds s.seed=RUNSEED env; original preserved).
- Array: baseline(step1_uniform)/t3p3/t2p3, NEW seed=7, 12M no-WW full-torus (~4h) -> reproducibility of 5x/+41%.
- XS path confirmed present (earlier failure was transient node-mount).
- PENDING BUILDS (more decisive, queue next): (d) half-amplitude t3p3 (+10cm) linearity; (e) clean zero-delta
  (shield 40) baseline for the equal-mass placement-vs-bulk comparison.
- PRELIMINARY placement-vs-bulk (needs zero-baseline to finalize): coil_step1b_uniform is ~an equal-mass uniform
  addition (+2.55cm, 630 cell-cm) and barely helps (peak 0.036->0.0345, ~4%), while t3p3 CONCENTRATED (+20cm/500
  cell-cm) gives 5x. Suggests placement >> bulk at equal mass -- the key engineering result, to be confirmed clean.
- NEXT: analyze 9472669 when done; build (d),(e); then Rung 4 (metrics/design-variable definitions).

## RUNG 3 attempt 1 INVALID (cwd race) -> relaunched
- All 3 seed=7 tasks shared cwd rung3/ -> raced on model.xml -> all loaded t3p3 geometry (identical coil20=6.40e-3). INVALID.
- Relaunched (rung3b) with ISOLATED per-task cwd. seed=1 originals remain valid (separate jobs).

## STANDING INTEGRITY RULE adopted 2026-08-24 (see memory spf-integrity-rule) — retroactive reclassification
- RUNG 2 is now BLOCKED, not "done". The importance->patch projection used a CIRCULAR-TORUS proxy for the
  VMEC/ParaStell shield surface = an UNCONTROLLED, UNVALIDATED approximation (the standing rule forbids this).
  Its anti-correlated/0-overlap result is treated as a mapping error, not physics. To UNBLOCK: project the coil-20
  adjoint importance onto the EXACT ParaStell shield reference surface using the actual (theta,phi)->(x,y,z) mapping
  from wout_qh.nc + the ParaStell build, with a geometry-overlay + reconstruction check. No A-vs-worth claim until then.
- Preflight is now MANDATORY before: (i) the exact Rung-2 re-projection; (ii) any "equal-mass/equal-cost" design
  comparison (compute mass from actual t_shield/t_breeder volumes x densities, not cell counts); (iii) volume-
  normalized heating density (Rung 4) using real cell/mesh volumes; (iv) any A vs B correlation.

## RUNG 2 (exact) UNBLOCK — in progress (integrity-rule compliant)
- Preflight written (rung2_exact_PREFLIGHT.md). No torus proxy.
- INTEGRITY FINDING: existing importance (adjoint_coil20_flat_P0) is on the w35f build; patches are on step1.
  Importance depends on intervening shield -> mismatch -> NOT used. Per-coil importance on disk only for coils 11-20.
- FIX launched (job 9489463): adjoint_importance.py on dagmc_qh_step1_uniform, coil 20, response flat, SAME
  (30,30,20) full-device mesh (-1800..1800, -550..550) -> step1-consistent importance covering the shield.
- REMAINING for the exact projection (after importance lands): (2) extract ParaStell step1 shield surface
  (theta,phi)->(x,y,z) via the SAME construct_invessel_build as build_step1 (pstl env); (3) sample importance over
  the shield shell per patch; VALIDATE with openmc.lib.find_cell geometry-overlay (points must land in shield cells)
  before any A-vs-worth correlation. If the ParaStell surface can't be extracted exactly, STOP and report (no proxy).

## LOCAL COMPUTE UNBLOCKED via lima (2026-08-24) — cluster contention no longer blocks GEOMETRY work
- lima VM "spf" HAS a full stack: OpenMC 0.15.3 (DAGMC=True) in conda env `spf` (bin ~/opt/openmc/bin,
  lib ~/opt/openmc/lib), ENDF/B-8.0 at ~/xs/endfb80/cross_sections.xml. Mac home mounts at /Users/tkiker (rw).
  I was wrong that it couldn't run OpenMC; it can (4 cores -> heavy transport impractical, but GEOMETRY is fast).
- RUNG 4 (critical-path) DONE LOCALLY: VolumeCalculation on dagmc_qh_step1_uniform (31s, no cluster) ->
  cell_volumes_step1.npz. Coils 11-30 volumes 6.58e6-7.63e6 cm3 (~16% spread -> "hottest coil" IS volume-confounded
  in principle). Dummy H1 materials (composition irrelevant for volume; names match DAGMC tags).
- VOLUME-NORMALIZATION RESULT (integrity check passed): baseline coil density <phi>=F/V still peaks at COIL 20
  (both integrated and density). So the peak-coil finding is ROBUST to volume normalization here; the ~13x flux
  spread dominates the ~16% volume spread. Reported as density (1/cm^2), not integrated -- metric now correct.
- IMPLICATION: geometry-only critical-path (volumes, exact added-mass, clearance field, ParaStell surfaces) can
  run LOCALLY on lima while the cluster queue is contended. Transport (coil runs, importance) still needs cluster.
- Cluster jobs still queued (volcalc 9489710 now redundant -> can cancel; rung3b 9489690; adj20 9489689 secondary).

## RUNG 5 step (d) hierarchical basis DONE (local): A=20 DOF (coarse sectors), B=12 (m<=3,n<=2, lam_pol 454-1363cm),
C=35 (m<=6,n<=4, finest lam_pol ~227cm ~= 1 coil). Stellarator-symmetric cosines. Files: hierarchical_basis_*.
## CLUSTER-DEPENDENCY (per integrity rule, NOT faked locally): ParaStell/cadquery/pymoab/netCDF4 NOT in lima ->
exact surfaces + clearance (steps b,c) and the zero-delta reference build for exact masses (step a) require a LIGHT
critical-path cluster job (pstl env). Local lima CAN do: VolumeCalculation on any transferred DAGMC (geometry-only),
the basis (done), and design-family delta construction (step e math).
## Cluster still queued: rung3b 9489690 (transport), adj20 9489689 (attribution, secondary).

## 2026-08-25 RESUME (claude-control)
- rung3b transport (9489690) RUNNING (indep-seed reproducibility); adj20 not in queue.
- masses cluster job 9490687 RAN+FAILED (pymoab not in spf-stellarator env). BUT zero-delta ref t0p0 built OK.
- FIX: pymoab-free local mass calc (openmc.lib cell->material), lima. Transferring 12 patches+t0p0 to lima_share.
- BASIS CORRECTED (rung5_nested_basis.md): invalid 20/12/35 mixed-family -> nested Fourier F0..F4 (5/13/32/59/124),
  nesting PROVEN (subset + lstsq resid 6.3e-15). Feasibility cap: grid Nyquist m<=9 pol, n<=6 tor -> can't resolve
  8-coil toroidal (needs n=8); needs finer ParaStell mesh. basis!=convergence (needs opt+transport stabilization).
- ParaStell NOT on Mac host either -> exact surfaces/clearance remain a genuine cluster dependency.

## 2026-08-25 EXACT MASSES + CLEARANCE DONE (RESULTS_masses_clearance.md)
- Surface extraction (cluster pstl, light) SUCCEEDED: get_loci() per-layer surfaces, 16-sector EXACT added-shield
  mass (222-457 t, 2.06x spread for same +20cm -> "equal thickness != equal mass" proven). Cross-check exact vs MC
  V_shield 7.19e8 vs 7.04e8 (2%). MC method bug fixed (shield vs thermal_shield substring).
- Clearance FIELD: min 29.7cm INBOARD (tor3,pol173), median 314; avail-to-grow min 4.7cm inboard. Design domain:
  I1 swap delta in[-30,+35]; I2 growth g in[0,avail]. KEY: tightest clearance is inboard = where coil-20 peaks ->
  response-need vs feasibility collide (the thesis, quantified).
- Refining: 40M MC rerun (lima, 2/14); indep-seed rung3b still transporting (2:31, no output). find_cell overlay TODO.

## 2026-08-25 DESIGN FAMILIES (matched-mass, design INPUTS not yet transport-evaluated)
- Masses CORRECTED to true normal-offset hex-vol (area*t was +9% biased vs true-hex/MC): per-sector 225-491t,
  spread 2.18x (robust). RESULTS_masses_clearance.md updated.
- design_families.py: uniform / clearance-driven / response-aware(B_FD c20) built on 4x4 sectors, ALL matched to
  M*=2.63e6 kg (uniform d=10cm) at 0.00% err -> equal-mass comparison VALID. Peaks: clearance t3p0, response t3p3.
- PENDING families: attribution-driven (Rung-2 exact importance->shield projection, NOW unblocked by surface_extract);
  wall-loading-driven (needs first-wall neutron-flux map). Next: build family DAGMCs -> transport -> equal-mass Pareto.
- Background still running: 40M MC mass rerun (lima, ground-truth cross-check ~6/14); rung3b indep-seed (cluster, ~2.5h).

## 2026-08-25 MASS TABLE FINALIZED + surface method CORRECTED (integrity)
- 40M MC finished -> AUTHORITATIVE per-patch masses (RESULTS_masses_clearance.md sec1): added mass 163t(t2p1)-
  400t(t3p3), spread 2.45x; uniform +2.55cm=529t; ALL patches lighter than uniform (t3p3 -129t, 8sigma) ->
  "t3p3=uniform equal-mass" REFUTED. swapResid<=6% (fixed-envelope swap confirmed).
- OVER-CLAIM CORRECTED: surface area*t OVERESTIMATES MC dV_shield by +6-19% (hex-vol "fix" was WORSE +20%, reverted);
  boundary patches (I_TOR 0,3) have a build_patch SEAM-extension artifact (t3p2 surface -32% vs MC). Surface loci
  geometry is correct (total-shell true-hex -2.6% vs MC) but area*t != swept vol for thin delta. USE MC for abs mass.
- design_families recalibrated to MC: M*=2.16e6 kg (was 2.63e6, r~0.85); uniform/clearance/response matched at 0.00%.

## 2026-08-31 EQUAL-MASS FAMILY TRANSPORT — COMPLETE + ANALYSIS GATED
- The minimax family was re-solved at the common 16-sector budget after an audit found its original 12-sector
  candidate was only 1.528e6 kg. All four build inputs now target **M*=2.163095e6 kg**.
- Cluster job **9569111** builds and evaluates `famU/famC/famR/famM` sequentially. Serialization is intentional:
  the ParaStell builders share a STEP export directory, so concurrent builds would recreate the known cwd race.
- `family_transport_analysis.py` + standard-library tests added. It verifies each JSON's declared peak against the
  coil vector, requires equal histories/cell ordering, propagates independent Monte Carlo errors, applies a 2sigma
  improvement/degradation gate, and reports hotspot migration only when the peak beats the runner-up by >=2sigma.
- First result `famU` completed at 4M histories: peak `2.5449e-2` on c23 with 3.5% relative error. It reproduces the
  prior 2M uniform result (`2.5372e-2` on c23) within statistics. Remaining family rankings are deliberately withheld
  until all JSONs exist and pass the same gate.
- `famC` completed next: peak `2.5509e-2` on c23, **+0.2% ± 4.8% vs `famU` (z=+0.05)**. The clearance-driven
  family is statistically indistinguishable from uniform.
- `famR` completed next: peak `2.5597e-2` on c23, **+0.6% ± 4.9% vs `famU` (z=+0.12)**. The response-aware
  family is also statistically indistinguishable from uniform, and the resolved hotspot remains c23. This is a
  negative result for the single-patch response heuristic, not evidence that placement never matters. `famM`
  remains in flight, so no four-family ranking is claimed yet.
- The original serialized job ended after the `famM` CAD layer split into 31 triangle sets while ParaStell supplied
  25 material tags. The wrapper printed `BUILD_FAIL_famM` but returned zero, so scheduler `COMPLETED` was not a
  valid family-success signal. Recovery job **9574368** reuses the repo's guarded split-solid tag expansion,
  preserves the eight magnet tags, and runs fail-closed (`set -euo pipefail`) before the same 4M forward solve.
- The guarded exporter reconciled the split exactly: the variable breeder became seven disconnected solids,
  expanding 25 tags to 31 volumes while retaining 8 magnet volumes. That topology shifts actual magnet cell IDs
  from 19–26 to 25–32, so the first recovery forward was cancelled before completion when its hard-coded IDs were
  detected. `forward_eval.py` now accepts an H5M-derived cell manifest and writes both actual tally IDs and stable
  physical coil labels; `family_transport_analysis.py` compares the stable labels and has a regression test for
  shifted geometry IDs. Fail-closed forward-only job **9574427** completed the corrected 4M evaluation.
- **FINAL EQUAL-MASS VERDICT:** `famM` peaks at `2.5222e-1` on physical c23 (actual DAGMC cell 29),
  **+891% ± 36% vs `famU` (z=+88)**, with 1.0% peak relative error and all 8 coils below 10%. All four
  peak-value and hotspot-label gates pass. Thus clearance and response are unresolved versus uniform, while the
  extrapolating linear-minimax stress test is decisively worse. No nonuniform candidate improves uniform at this
  mass. Production optimization needs bounded/sequential direct-forward relinearization or nonlinear design data.
- Minimax surrogate audit added (`minimax_surrogate_audit.py`, `MINIMAX_SURROGATE_AUDIT.md`). The 2.163e6 kg
  budget is reachable under a 20 cm sector cap, but the unconstrained superposition still predicts 13/20 negative
  coil responses; adding `R_pred>=0` makes the LP infeasible at both 20 and 35 cm caps. `famM` is therefore an
  extrapolation stress test. A production minimax needs sequential direct-forward re-linearization/nonlinear data.
- Meeting artifacts: `fig_peterson_summary.png`, `PETERSON_TALKING_POINTS.md`, and
  `../../../talk/peterson_checkin_2026-08-31.pptx`. The deck's family slide is regenerated from
  `family_transport_results/family_transport_summary.json` and switches to a quantitative chart only when complete.
