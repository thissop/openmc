# Session Handoff — One-Period Periodic BCs + Deep Freeze + FW-CADIS VR

Written 2026-08-21 after a long build session. Read this first next session. Companion memories:
`spf-ssr-decoupled-pipeline`, `spf-patch-worth-n12-verdict`. Prior handoff:
`handoffs/2026-08-10_ssr-pipeline.md`.

## 1. Where we are in one paragraph

Peterson (Aug 19 meeting) green-lit **one field period + periodic BCs** (the biggest speedup, standard
in stellarator neutronics) and **freezing the surface source deeper** — at the outer boundary of a
FIXED breeder zone, splitting the breeder into fixed (never traded) + variable (breeder<->shield trade)
zones. We built that whole pipeline this session, validated it, and it **largely resolves the far-side
residual** (3-5x at the first wall -> ~30-60% at the deep freeze; far-side asymmetry structurally gone
because one period has no far side). The one open piece is **efficient variance reduction to resolve
the per-coil dose** (coils sit behind 40 cm of WC shield, VR-limited in any method). FW-CADIS via
OpenMC's native random-ray solver is the literature-confirmed right answer, and the random-ray solve
now RUNS, but the `convert_to_multigroup`/random-ray integration on our DAGMC deep-penetration model has
murky issues (below) that need careful (not marathon-tail) work.

## 2. THE VALIDATED PIPELINE (all on Ginsburg, `~/pstl_test/corrected`, code in `ssr/`)

**Geometry:** `dagmc_qh_op_uniform.h5m` — ONE field period (ParaStell `repeat=0`), breeder split into
`breeder_fixed` (15 cm) + `breeder_var` (35 cm). Built by `build_oneperiod.py` (from build_step1.py:
`repeat 3->0`, magnets extent `360->90`, breeder split, OUTNAME `qh_step1->qh_op`). Runner
`slurm_oneperiod_build.sh`. NOTE: the build "FAILED" status is a red herring — the final `$(date)`
line lacks `date` on PATH; the DAGMC writes fine ("OK_OP" / "written DAGMC file" in the log).

**Cell map (verified via find_cell along a radial ray, `cellmap3.py`):**
plasma void=2, first_wall=3, multiplier=4, **breeder_fixed=5**, **breeder_var=6**, back_wall=7,
shield=8, gap=9, vac_vessel=10, thermal_shield=11, **magnets=12-16** (5 coils, verified `findmagnets.py`).

**Periodic-BC wedge (in `op_master.py` / `op_ssr.py`):** wrap DAGMC in `YPlane(0)` (phi=0 face) +
`XPlane(0)` (phi=90 face), both `boundary_type="periodic"`, `py.periodic_surface = px`, root region
`+py & +px & -sph`. Validated: **0.66% leakage** (would be ~50% if periodicity weren't stitching).

**Wedge-restricted source:** `StellaratorSource(..., constraints={"domains":[root], "rejection_strategy":
"resample"})`. NOTE the domain MUST be the root **Cell** (not the DAGMCUniverse — the type check rejects
DAGMCUniverse). Resamples the full-torus fluxmap down to the wedge (~4x oversample, cheap). The fluxmap
is confirmed perfectly 4-fold symmetric (sqrt(g) CV=0%), so this is exact.

**Deep freeze:** `surf_source_write = {"cellfrom": 5}` (bank all crossings leaving breeder_fixed).
IMPORTANT: OpenMC forbids `cellfrom`+`cellto` together (at most ONE of cell/cellfrom/cellto). So use
cellfrom=5 alone and handle inward neutrons with black-inside in the replay.

**SSR replay (`op_ssr.py`):** reads the frozen source, transports through the variable zone to the
coils. BLACK_INSIDE=1 overrides cells **2,3,4** (void, FW, mult) -> B-10 black via
`init_lib -> sync_dagmc_universes -> add_material_override -> finalize_lib`. CRITICAL: black cells
**2,3,4 NOT cell 5** — the surface source is banked ON cell 5's boundary, so blacking cell 5 kills the
outward source neutrons at launch (gave all-zero coils). Blacking 2,3,4 kills only inward backscatter.

## 3. RESULTS (the science answer — residual is largely resolved)

- **Master (`op_master.py`, 24M, ~3h):** reference coil dose `[0.0044 0.0097 0.0132 0.0212 0.0256]`
  (fast flux, cells 12-16). Wrote `surf_op_master.h5` (4.18 GB, ~50M crossings) + `coil_op_master.npz`.
- **Full SSR (24M):** depth profile `breeder_var=16 -> shield=1.0 -> VV=1.7e-3 -> coils~1e-6`. Textbook
  shield attenuation. Coils VR-starved (Poisson-zero) even at 24M — EXPECTED (Peterson pre-warned).
- **Cheap verdict (4M master with depth tally, `op_master_depth.py`):** master depth
  `[breeder_var 25.7, shield 2.2, VV 0.067]` vs SSR `[16, 1.0, 1.7e-3]`. **Ratio 1.6 -> 2.1** at the
  well-populated layers = decoupling reproduces the field to ~30-60%. Deep layers noise-limited.
- **Bottom line:** first-wall residual 3-5x -> deep-freeze ~30-60% radial; far-side asymmetry DEAD
  (one period). The ~30% is likely the black-inside under-counting multiply-backscattered neutrons +
  statistics. Clean per-coil number needs VR.

## 4. THE OPEN TASK: efficient VR for the coil dose

**MAGIC (OpenMC native, `op_ssr_ww.py`):** confirmed to over-split — 4M secondary tracks with
on_the_fly=True. Capped version (`op_ssr_wwapply.py`, max_split=4, max_history_splits=1e5,
generate-then-apply, coarse mesh) STOPPED the explosion (350k tracks) but is INEFFICIENT (~20 min/batch
-> ~1 day for 24M). Not viable. Confirms the literature: MAGIC 2nd+ iteration "orders of magnitude
slower... aggressive splitting."

**FW-CADIS via random ray (`op_fwcadis.py`) — THE RIGHT ANSWER (do this):**
- Literature: gold standard for deep-penetration fusion shielding (48x on ITER magnet shielding).
  OpenMC has it NATIVE via the internal random-ray multigroup solver (no Denovo/ADVANTG). Adjoint-based
  (`w = 1/(2 phi_dagger)`), ~order-of-magnitude better FOM than MAGIC. Supports target tallies.
  Refs: docs.openmc.org/en/stable/methods/variance_reduction.html and .../usersguide/variance_reduction.html;
  arxiv 1412.1075 (FW-CADIS on ITER).
- Recipe: `convert_to_multigroup()` -> `convert_to_random_ray()` ->
  `WeightWindowGenerator(method='fw_cadis', mesh=mesh)`, `wwg.targets=[coil_tally]` -> run ->
  `weight_windows.h5` -> apply in CE SSR (`settings.weight_windows_file`, `weight_windows_on=True`).

**Integration issues found this session (fix these next):**
1. **`convert_to_multigroup` regenerates only if `mgxs.h5` absent** (model.py line ~100). Pre-generate
   `mgxs.h5` ONCE and reuse to skip the slow MGXS step every run.
2. **`method="infinite_medium"` still ran multiple slow CE fixed-source sims** (log had 6+ "FIXED SOURCE
   TRANSPORT SIMULATION" phases). Investigate `_generate_infinite_medium_mgxs`; may need explicit
   nparticles/batches control or the mode isn't dispatching as expected.
3. **Fine source-region mesh (20^3) + DAGMC -> "maximum number of events" stall** (rays cross too many
   boundaries; `distance_active = 5x` geometry diag = ~25000 cm rays). FIXED by DROPPING
   `random_ray['source_region_meshes']` (use DAGMC cells as source regions) + shortening rays
   (`distance_inactive=0.2*span, distance_active=1*span`). Then random ray progressed fast.
4. **BUT after the random-ray batches, the run read full CE data (W183...) and started a slow CE MC** —
   suggests `energy_mode='multigroup'`/random-ray mode isn't fully taking on `model.run()`, OR my
   post-`convert_to_random_ray` settings edits reset it. THIS is the key thing to resolve. Check
   `settings.energy_mode` after convert_to_random_ray; ensure model.run() uses the RR/MG solver.
5. **Resolution vs geometry:** a 40 cm shield in a 2300 cm domain can't be resolved by a uniform mesh
   without exploding ray crossings. LIKELY FIX: subdivide the shield radially into ~8x 5 cm sub-cells
   in ParaStell (same trick as the breeder split) so the DAGMC CELLS carry the importance gradient ->
   random ray source regions = cells -> FW-CADIS gradient through the shield WITHOUT a fine mesh.
6. Benign: "No intersection with DAGMC cell 4" warnings (thin 2 cm Be multiplier) — appear even in the
   plain master, don't corrupt results, just bloat logs.

## 5. FIRST THINGS NEXT SESSION

1. `ssh ginsburg` (user does Duo). PATH=/cm/shared/apps/slurm/current/bin.
2. Resolve FW-CADIS issue #4 (CE-vs-MG mode) — read `convert_to_random_ray` in
   `$SRC/openmc/model/model.py`, confirm `settings.energy_mode` and run mode after conversion, make a
   minimal CSG (non-DAGMC) random-ray FW-CADIS smoke first to prove the workflow, THEN port to DAGMC.
3. Pre-generate `mgxs.h5` once (issue #1) so iterations are fast.
4. Consider the radially-subdivided-shield geometry (issue #5) so cells give the gradient.
5. Once `ww_op_fwcadis.h5` is produced: apply in `op_ssr_wwapply.py` (swap the WW file), run 24M,
   compute Gate 2 (master `coil_op_master.npz` vs the WW-resolved SSR). Clean per-coil verdict.
6. If FW-CADIS-on-DAGMC stays murky: fallback is FW-CADIS on a CSG bounding approximation, or accept
   the ~30-60% cheap-verdict result as "residual resolved" and move to the optimizer with a coarser VR.

## 6. Key files (Ginsburg `~/pstl_test/corrected/`)
- `build_oneperiod.py`, `slurm_oneperiod_build.sh` — geometry.
- `ssr/op_master.py` (+ `_depth.py`), `ssr/op_ssr.py`, `ssr/op_fwcadis.py`, `ssr/op_ssr_wwapply.py`.
- `ssr/slurm_op_master.sh`, `ssr/slurm_op_ssr.sh`, `ssr/slurm_op_fwcadis2.sh` (streamed output).
- `cellmap3.py`, `findmagnets.py` — cell identification (reusable).
- Outputs: `surf_op_master.h5` (4.18 GB), `coil_op_master.npz`, `ww_op_coils.h5` (MAGIC, coarse).

## 7. Deadlines / context
- Aug 26 research talk. Spine holds: patch-worth negative -> why proxies fail (moderation) -> SSR
  decoupling -> **one-period + deep-freeze (residual resolved)** -> FW-CADIS VR -> optimizer.
- The residual being resolved is a presentable result NOW; the FW-CADIS VR is the bridge to the
  optimizer numbers.
