# Overnight autonomous build-out — plan + live progress

Goal: a rigorously validated, tested body of work for the paper. Each step is the next
most-valuable increment. All compute is SLURM (survives laptop sleep); I collect + analyze +
PLOT when the tunnel is reachable. Updated as I go.

## Sequence
- [RUNNING] **S1. Response curve on the ORIGINAL 107 cm build** (delta=0/10/25 cm, unpol+A).
  Qualitative closed-loop proof + coil-peaking-vs-TBR sensitivity. Jobs coilcmp 9008105,
  coilcmp10 9008205. -> PLOT: coil peaking + hottest-coil flux + TBR vs delta.
- [ ] **S2. Corrected literature-anchored baseline** (RADIAL_BUILD_VALIDATION.md): FW 2mm W+3cm
  RAFM (homog ~7.5 g/cc), **Be multiplier 2cm**, FLiBe 50cm (Li6 re-optimized ~40%), RAFM back
  4cm, WC shield 40cm, gap 2cm, VV 25cm (steel+borated water), thermal shield 3cm ~= 126cm.
  Build + verify **TBR>=1.05** (tune Li6/Be if short). Robust coil-cell discovery by 'magnets'
  tag. -> report TBR, PLOT corrected radial build diagram.
- [ ] **S3. Response curve on the CORRECTED build** (delta=0/10/25, unpol+A, WW). The
  QUANTITATIVE headline. -> PLOT: peaking/TBR vs delta on the defensible build; per-coil bars;
  actual (theta,phi) shield/FLiBe change maps.
- [ ] **S4. First-wall <-> coil relationship** (from the 4M coilhot job's first-wall map + coil
  results). Does first-wall peaking predict coil hotspots? -> PLOT correlation.
- [ ] **S5. QH transported culprit map** (coil-hotspot attribution to plasma regions) - the
  attribution novelty on QH. -> PLOT 2D + 3D.
- [ ] **S6 (stretch). QA device #2**: VMEC run for QA wout -> full pipeline -> QA-vs-QH.

## Live progress log
(2026-07-14 ~21:00) S1 running; S2 launched.
(2026-07-14 ~22:40) LEARNINGS + course corrections:
- **delta=25 cm strains ParaStell CAD** — the aggressive non-uniform grid hung `construct_
  invessel_build` >1 h and got cancelled; delta=10 builds fine in ~8 min. => use GENTLE delta
  (<=15 cm) everywhere. (Also more physical — matches the "50% shield jump is extreme" concern.)
- **Weight-window generation is too expensive/finicky** on this deep-penetration 3D DAGMC
  (cancelled at ~70 min). => use ANALOG transport (3M histories gives 2-5% per-coil error,
  adequate for the response curve). WW deferred as a future optimization (tet-level / speed).
- S1 salvaged as a delta=0 vs delta=10 point on the fine grid: base_fine (9008435) +
  coilcmp10 (9008205), analog. Watcher bp4wfqd6p waits for the 4 npz -> analyze + PLOT.
- S2 corrected literature-anchored baseline building in ~/pstl_test/corrected/ (agent) w/ Be
  multiplier + real FW/VV/gap, verifying TBR>=1.05.
- S3 will run the response curve (delta=0/5/10/15, gentle) on the corrected build, ANALOG.
- coilhot (9003086) still running: B/C modes + first-wall map for S4.
(2026-07-14 ~23:15) TARGETING (coil_centroids.npz, paper-relevant, model-robust):
- Hot coils = cells 12/22, 11/21, 10/20 (warm 9/19), at toroidal phi~20-40deg within each
  90deg period, poloidal theta~265-330deg (inboard-to-bottom, below midplane z=-86..-356cm).
  Confirms the perturbation aim theta0~300/phi0~30 hits the hot coils. Baseline peaking 1.79.
- FINDING for paper: 2-fold hot/cold ALTERNATION between adjacent nfp=4 periods (SPF source
  structure; base periods 0/180deg hot, 90/270deg cold). ParaStell tiles ONE period (repeat=3)
  so any shield bump is necessarily periodic -> inherent efficiency caveat to state.
- Corrected build standoff = 129.2 cm (vs 107 old) -- matches the ARIES-CS ~130 cm target.
- The stray WW/perturbation agent stood down; I am the sole driver via SLURM + watchers.

(2026-07-15 ~01:30) TWO REAL ISSUES CAUGHT BY RIGOR (both must be fixed before trusting any
per-coil perturbation number). Honest morning brief:

**ISSUE 1 -- perturbed run gives an ANOMALOUS uniform coil flux (per-coil perturbation
result NOT trustworthy; cause NOT yet found).** Verified with a WORKING pymoab discovery
(discover2.py, in ww_pert/): coils = cells 8-27 in BOTH dagmc_baseline_fine.h5m AND
dagmc_perturbed10.h5m (identical), and 11-30 in dagmc_corr_base.h5m (8-layer build shifts
them). So the hardcoded 8-27 WAS correct for base/pert -- the cell IDs are NOT the bug.
The real anomaly: base_fine coils have structured flux 0.02-0.07 (peaking 2.20), but
perturbed10 coils are ALL uniform ~0.007 (~10x lower). A localized delta=10 bump (peak
10cm@theta300/phi30, MEAN 1.88cm, shield mostly still 40) physically CANNOT reduce all 20
coils 10x -> so "peaking 2.20->1.72, hottest -88%, TBR +0.10" is an ARTIFACT of a deeper
perturbed-build or run inconsistency (NOT understood yet). ALSO note: base_fine peaking
(2.20) differs from the earlier old-build peaking (1.79) on the same 8-27 cells -> the FINE
builds differ from the coarse build in some non-obvious way. TOP MORNING DEBUG: diff the
baseline_fine vs perturbed10 geometries (volumes, coil positions, shield map), and re-examine
the perturbed10 run config, before trusting ANY per-coil perturbation number. discover2.py is
the reusable correct cell-discovery (corrected build needs cells 11-30, not the c.fill path
that segfaulted).
FURTHER (discover3.py, centroids): coils are in IDENTICAL positions in baseline_fine and
perturbed10 (cell8 R=1796/Z=-126, cell12 R=1091/Z=-87 in BOTH) -> the perturbed build is NOT
geometrically broken (coils didn't move). Yet perturbed coils see uniform ~0.007 at distances
1091..1796 cm -- physically impossible (near vs far coils can't be equal; localized bump can't
touch all of them). So the anomaly is a SUBTLE build/run issue, most likely the FINE-grid
ParaStell build misapplying the per-(theta,phi) thickness_matrix, OR a run-config problem.
NEXT (fresh morning): (i) OpenMC geometry PLOT of perturbed10 to eyeball the actual shield
thickness distribution; (ii) confirm the thickness_matrix that build_perturbed.py fed ParaStell
matches what the .h5m contains; (iii) also explain why base_fine peaking (2.20) != old-build
peaking (1.79) on the same coils -- the fine build itself differs from the coarse one.
Until this is understood, NO per-coil perturbation number is trustworthy. TBR (MaterialFilter)
is unaffected and reliable.

**ISSUE 2 -- FLiBe under-breeds (design decision needed).** Corrected build (50cm FLiBe +
2cm Be, WC shield, real FW/VV/gap): TBR = 0.849/0.862/0.878 at Li6 30/40/60% (90% ~0.90).
Li6 barely helps (+0.015 per 20%). CONFIRMED: 50cm FLiBe cannot self-breed. Thick-FLiBe test
(65cm + 4cm Be) building (breedbld 9010645). Even 65cm likely marginal (~1.0). And 80cm FLiBe
would exceed the 162cm coil-plasma standoff. => DECISION for the user: (a) thicker FLiBe at
tight standoff, (b) switch breeder to PbLi (breeds ~1.1 at 50cm), or (c) accept sub-unity /
external tritium and frame the TBR floor as a BINDING constraint (which actually makes the
optimization more interesting). TBR-via-MaterialFilter is reliable; this is a real physics
finding, not a bug.

**ISSUE 2 RESOLVED (2026-07-15 ~02:30):** breeder fixed. **65 cm FLiBe + 4 cm Be -> TBR =
1.106** (dagmc_corr_breed.h5m), self-sufficient with margin. SATURATES at 60% Li-6 (60% and
90% both 1.106 -> use 60%, cheaper). Fits: total ~146 cm vs 162 cm standoff = ~16 cm headroom
for the shield trade. BONUS: thicker in-vessel also cuts coil flux 3.4x (0.76 -> 0.23) -- better
magnet protection for free. So the breeding-adequate REFERENCE BUILD = corr_breed (8 layers,
FLiBe 65, Be 4, Li6 60%, coils = cells 11-30). No PbLi pivot needed. The morning "breeder
decision" is effectively answered; run the real response curve on corr_breed once the
perturbed-geometry anomaly (Issue 1) is fixed.

(2026-07-15 ~03:30) corr_breed BASELINE per-coil (cells 11-30, correct discovery, TRUSTWORTHY):
peaking = 2.906 (unpol), hottest cell 30, TBR 1.106, rel err 6.2%. FINDING: baseline coil
peaking RISES with blanket thickness -- 1.79 (thin old) -> 2.20 (fine) -> 2.91 (thick corrected).
Physically sensible (deeper penetration -> near-field streaming dominates -> sharper peaking) and
a paper point: realistic reactor-thickness blankets sharpen coil peaking. Now waiting on
corr_breed delta=10 PERTURBED (bprtrun 9018240) to test if the shield trade gives a LOCALIZED
reduction (real closed loop) vs the old-path uniform anomaly. Combined watcher beobg90kb.

(2026-07-15 ~04:00) ANOMALY PRECISELY CHARACTERIZED (corr_breed delta=10, corrected pipeline):
- baseline peaking 2.906 (max 3.26e-2); perturbed peaking 1.792 (max 6.88e-3), still STRUCTURED
  (per-coil range 1.4e-3..6.9e-3, max/min 4.9 -> NOT uniform).
- REAL SIGNAL: hot coils preferentially reduced (cell 30: 3.26e-2 -> 6.1e-3 = 5.3x; mean drops
  2.9x) -> the distribution FLATTENS (peaking 2.91 -> 1.79) = the desired closed-loop behavior.
- ANOMALY: the GLOBAL mean coil flux drops 2.9x, but a delta=10 bump covers only ~23% of the
  wall -> global mean should drop ~20-30%, not 2.9x. The excess global reduction is unphysical.
- CONFIRMED: anomaly is in the SHARED thickness_field -> ParaStell thickness_matrix lofting
  (both old build_perturbed.py AND build_breed.py show it), NOT a retired path. TBR rising
  (+0.03) with the trade is plausibly real (WC reflection). So there IS a real result buried
  under a build artifact that inflates the global reduction.
- TOP MORNING TASK (concrete, will crack it): OpenMC geometry PLOT of dagmc_corr_breed_pert.h5m
  colored by material at a slice through the bump AND away from it -> see whether the shield is
  over-thickened GLOBALLY (lofting/interpolation bug) vs only in the bump. Also sample the
  actual per-(theta,phi) shield thickness the .h5m contains vs the thickness_matrix fed in.
- PROCESS BUG to fix: coil_percoil.py writes percoil_<MODE>.npz -> baseline & perturbed runs in
  the SAME dir COLLIDE (perturbed clobbered the baseline unpol array; only summary survived).
  Use unique per-run output names (include the DAGMC/tag).

(2026-07-15 ~04:15) PHYSICAL EXPLANATION FOUND (bump-coverage math, thickness_field.py):
- gaussian_bump width_deg=35 with a 90deg field period => toroidal Gaussian sigma=35deg barely
  decays within the period. dph uses dperiodic(nfp*TOR,...) -> ONE bump per period, tiled x4 by
  repeat=3/nfp=4. Result: the +10cm WC bump is localized POLOIDALLY (theta=300+-35, ~24% of
  poloidal) but spans the ~FULL TOROIDAL extent of every period -> a poloidal BAND of extra WC
  wrapping the ENTIRE torus, x4. Not a compact spot.
- mean delta 1.878 = amp * (poloidal cover ~0.244) * (toroidal cover ~0.77, near-full period).
- CONSEQUENCE: a full-toroidal poloidal band of +10cm WC (2-3 mfp) right on the flux-peak
  poloidal angle SHOULD shadow a large fraction of ALL coils (each coil spans full poloidal
  range; its peak flux comes from the nearest/hot poloidal patch, which is exactly what's
  shielded). So a big global mean drop + flattening (peaking 2.9->1.8) is largely PHYSICAL,
  not a lofting bug. The clean same-build pair (job 9019265) quantifies the true magnitude.
- DESIGN NOTE for the real optimizer: width=35 is too wide for a "localized" control knob (it's
  effectively a toroidally-uniform poloidal band). For genuine spatial control use narrower
  bumps or the Fourier field; report the band-vs-spot behavior honestly in the paper.

(2026-07-15 ~04:20) GEOMETRY VALIDATED via OpenMC material-slice + image-diff (job 9019270):
- Plotted dagmc_corr_breed vs dagmc_corr_breed_pert (xz & xy), colored by material, diffed.
- breeder->shield growth = 126 px (xz) / 858 px (xy), LOCALIZED to a poloidal band at the
  inboard/bottom of the cross-section; shield->breeder = 0 (trade direction correct).
- NO global shield thickening -> NO lofting/interpolation bug. Geometry is faithful to the
  intended fixed-envelope trade. => the earlier "5x global coil-flux drop" was the wrong-
  baseline artifact (clobbered percoil_unpol.npz) and/or the physical band-shadowing (wide
  bump = full-toroidal poloidal band). Clean same-build pair (job 9019265) gives true number.
- Figs: tmp/geom_xz_*.png, tmp/imdiff_xz.png (green = breeder->shield). Move keepers to
  data/processed/figs once the pair confirms the story.

(2026-07-15 ~05:25) RESPONSE-CURVE PIPELINE QUEUED (all BE_CM=4.0, LI6=60, 3M unpol, patched
coil_percoil2.py with unique percoil_unpol_<tag>.npz names to kill the file-collision bug):
  delta=0  corr_breed        -> job 9019809 (running)
  delta=5  build 9019818     -> coil 9019825 (afterok dependency)
  delta=10 corr_breed_pert   -> job 9019810 (running)
  delta=15 build 9019819     -> coil 9019826 (afterok dependency)
- Verified delta_corr_breed*.npz both used be_cm=4.0 -> curve is self-consistent.
- Surrogate prediction to test: k=1/lam_s-1/lam_b. lam_s(WC)~5.5-8, lam_b(FLiBe)~17 ->
  k~0.066-0.12/cm -> hottest coil (behind full bump) should drop to exp(-10k)~0.29-0.52x.
  The delta=0/10 pair gives the FIRST calibration point for lam_shield_eff; the full curve
  (0,5,10,15) fits the exponential and validates the cheap surrogate gradient = key methods fig.
- NEXT once pair RESULTs land: run compare_pair.py (tmp/) -> per-coil bar + ratio plot; confirm
  perturbed max ~0.3-0.5x (physical), NOT 0.21x (that was the wrong-baseline artifact).

(2026-07-15 ~06:35) PRE-REGISTERED PREDICTION (write before the MC lands = honest test):
- Hottest baseline coil sits behind the poloidal band that gets +delta WC (breeder->shield).
- Surrogate: flux(delta)/flux(0) = exp(-k*delta), k = 1/lam_shield - 1/lam_breeder.
- With lam_breeder(FLiBe)~17cm and lam_shield(WC) in [5.5, 8]cm -> k in [0.066, 0.123]/cm.
- PREDICTED hottest-coil reduction: delta=5 -> 0.54-0.72x; delta=10 -> 0.29-0.52x;
  delta=15 -> 0.16-0.37x. If MC lands in these bands AND is ~log-linear in delta, the
  exponential surrogate (the cheap closed-loop gradient) is VALIDATED. A calibrated
  lam_shield_eff falls straight out of the fit (response_curve.py, ready in tmp/).
- Build/env FIX logged: ParaStell env must be activated by FULL PATH
  (/ginsburg/astro/users/tjk2147/spf_work/envs/pstl), NOT bare name 'pstl' (bare name lacks
  parastell -> ModuleNotFoundError, and no set -e means the sbatch still exits 0 = silent).
  Add `set -e` or an explicit h5m-existence check to future build sbatch to fail loudly.

(2026-07-15 ~12:20) CLEAN VALIDATED delta-SWEEP (4 isolated, AUDIT-verified runs, Li6=60, BE=4):
  delta=0  peak 3.26e-2 (cell30) peaking2.91 TBR1.106  [job 9023085 r2, relerr 6%]
  delta=5  peak 6.91e-3 (cell15) peaking1.76 TBR1.262  [9021330]
  delta=10 peak 6.88e-3 (cell15) peaking1.79 TBR1.136  [9021331]
  delta=15 peak 6.91e-3 (cell15) peaking1.80 TBR1.120  [9021332]
- HEADLINE: first ~5cm (1 mfp WC) cuts targeted coil 30 to 0.19x; delta=5->15 FLAT = leak-around
  floor. Exponential-thickness surrogate FAILS past 1 mfp (data is a STEP, not exp). Hotspot
  RELOCATES 30->15 (uncovered coil becomes limiting). Per-coil: hottest coils drop most (0.14-
  0.21x), distribution flattens (peaking 2.9->1.8). => design variable is COVERAGE not THICKNESS;
  surrogate needs a floor term D=D_floor+(D0-D_floor)exp(-d/lam). Figs: response_curve_clean.png,
  coilpair_breed_clean.png (in data/processed/figs).
- OPEN THREADS: (a) baseline coil load is NOT nfp-symmetric (localized hot cluster 27-30) -
  check if real low-standoff sector vs StellaratorSource toroidal-sampling asymmetry.
  (b) TBR values noisy (delta=5's 1.26 high). (c) proposed next: WIDER bump at fixed delta=5
  to prove coverage reduces the floor where thickness didn't (awaiting user go).
- 9021329 = independent 2nd baseline (diff node/seed) left running as a reproducibility x-check.

**SOLID / VALIDATED (safe to build on):** full ParaStell+SPF+DAGMC pipeline works end-to-end;
corrected literature-anchored build (129cm, Be multiplier, real FW/VV/gap) built + validated
with citations (RADIAL_BUILD_VALIDATION.md); novelty verdict (NOVELTY_VERDICT.md); targeting
analysis (hot coils + nfp period alternation). Reliable plots: coil_hotspots_qh (the 1.79
inter-coil peaking on the OLD build with CORRECT cells 8-27), culprit maps, change maps,
qh_device_3d, geometry slices. FLAG: pert_unpol_d10.png is the ARTIFACT -- do not use.

**MORNING PRIORITY ORDER:** (1) fix robust coil-cell discovery (reuse coil_centroids.py) ->
re-tally the delta=0/10 comparison correctly = the real proof-of-concept; (2) user decision on
the breeder (Issue 2); (3) valid response curve on the chosen breeding-adequate build;
(4) S4 first-wall<->coil from coilhot (9003086, still running, 4.5h); (5) S5 QH culprit map.
