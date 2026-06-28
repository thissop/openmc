# GINSBURG_CLAUDE_HANDOFF.md — mission brief for Claude Code running on Ginsburg

You are Claude Code on **Columbia Ginsburg** (the dev work was done on a separate
aarch64 sandbox that **cannot run DAGMC/OpenMC-with-DAGMC**, so the x86 transport
path here has **never been executed** — treat this as bring-up). Your mission:
**verify this spin-polarized-fusion (SPF) stellarator OpenMC package is set up
correctly, then run smoke SLURM jobs to confirm it works end-to-end.** Work on
branch `spf-prototype`. You MAY commit to `spf-prototype` (never to `develop`/main).

## Read these first (authoritative; don't duplicate, follow them)
- `spf_prototype/SETUP_GINSBURG.md` — the setup sequence + storage guidance.
- `spf_prototype/RUN_ON_GINSBURG.md` — the conceptual pipeline + smoke ladder.
- `spf_prototype/CONTEXT_REPORT.md` — what the whole project is.
- `spf_prototype/DEFERRED.md` / `JOINT_VALIDATION.md` — scope limits + the analytic↔MC plan.

## Preflight findings already gathered (don't re-discover; verify if unsure)
- Account `astro`; `module load anaconda/3-2023.09` (default) →
  `source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh` → `conda activate`.
- `HOME=/burg-archive/home/tjk2147` is the **cold archive tier, ~37 GB free**; **no
  scratch mount**; **no `/burg/astro/users/$USER`** yet. 37 GB IS enough for the
  ~13 GB footprint (env ~6 + cross sections ~5 + outputs ~2) — not blocked, but
  prefer warm `/burg` if you find roomy space (`df -h`).
- **No shared cross-section library exists** (the one real setup blocker — see below).
- SLURM: default partition `short` (12 h max) covers the 8 h baseline; you also have
  `burst` (14-day) for long/scan jobs. git/cmake/g++ present.
- `spf_prototype/ginsburg_preflight.sh` re-runs all these checks read-only.

## Setup sequence (login node; the ONLY place with internet)
1. `conda env create -f spf_prototype/environment.yml` (pins `nodefaults` →
   conda-forge only; Ginsburg's base conda has only `defaults`, whose ToS can break
   env creation). `conda activate spf-stellarator`. Confirm:
   `python -c "import openmc; assert hasattr(openmc,'DAGMCUniverse')"`.
   - **First, dry-run the solve** (catches any DAGMC-variant issue before the big install):
     `conda create -n _t --dry-run -c conda-forge "openmc=0.15.*=dagmc*" | tail`.
2. **Cross sections (REQUIRED — none shared):** download ENDF/B-VIII.0 HDF5 once
   (~2-3 GB) from https://openmc.org/official-data-libraries/ to roomy storage,
   `export OPENMC_CROSS_SECTIONS=.../cross_sections.xml`, and put the SAME path in
   `ginsburg_job.sh`.
3. Build the compiled source against the conda OpenMC (once):
   `cd spf_prototype/src && cmake -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" . && cmake --build build`
   (`run_ginsburg.py` will reuse this `.so`; it skips rebuilding when present.)

## Verification checklist (do these in order; STOP at the first failure and fix)
1. `python -m pytest spf_prototype/tests -q` → **expect 113 passed** (sampler
   C++/Python parity, anarrima angled-kernel <1e-6, field-map parity, geometry
   watertight/simple). This is architecture-independent and must pass cleanly.
2. **Tiny end-to-end on the LOGIN node** (the real bring-up gate; ~minutes):
   `OMP_NUM_THREADS=4 python spf_prototype/python/run_ginsburg.py --stem equil_precise_qa --scale 10 --particles 20000 --batches 5 --results-dir /tmp/spf_smoke`
   It bootstraps the geometry + DAGMC `.h5m` + scaled field map, runs 6 configs, and
   must write `/tmp/spf_smoke/RESULTS_tier8_conformal.md`. **Open that file and
   confirm the numbers are NOT `nan`** (a metrics blocker was fixed but is untested).
3. Then a **small SLURM job** to confirm batch mechanics + persistence: copy
   `ginsburg_job.sh`, set its paths (REPO, XS), drop `--particles` to ~100000, and
   `sbatch` it. Confirm results land in `$SLURM_SUBMIT_DIR/results_*` (NOT `/tmp`).
4. Only then the **baseline** (`ginsburg_job.sh` as-is, 2M×20).

## Watch-items / known risks (THIS is the high-value part — my checks can't cover these)
1. **Conformal-wall toroidal self-intersection (most likely to bite).**
   precise_QA is geometrically strongly shaped — its cross-section CENTER swings
   from R≈120→81 cm (native; ~39 cm, *larger than the ~17 cm minor radius*),
   rotating a tall bean (φ=0) into a wide ellipse (φ=π/2). My geometry guard checks
   only **poloidal** self-intersection per φ-slice (passes at `--scale 10`); it does
   **NOT** check **toroidal** self-intersection (adjacent φ-slices of the offset wall
   crossing on the inboard side), and `watertight` (edge-manifold) won't catch that
   either. **In the tiny run (step 2), watch for `lost particles`, negative-volume,
   or overlap warnings from OpenMC/DAGMC.** If they appear, the conformal offset
   folded. Fixes, in order: raise `--scale` (thinner build relative to minor
   radius), thin `stellarator_geometry.DEFAULT_LAYERS`, or use ParaStell (the
   production conformal mesher; `DEFERRED.md`). Verify the `.h5m` with
   `python -c "import openmc; u=openmc.DAGMCUniverse('stellarator.h5m'); print(u.material_names)"`.
2. **`stl_to_h5m` API drift** — `build_dagmc.py` has the single call site flagged;
   if it errors, fix against the installed version (the dry-run/tiny run surfaces it).
3. **Coil fast-flux is deep behind blanket+shield** → large MC error without
   variance reduction. Lead with the trend / report the ±; a precise coil η needs
   weight windows (not yet wired).
4. **DESC** (only for NEW equilibria / the QA→QH scan) must be a SEPARATE venv
   (`pip install desc-opt`) — its jax/numpy pins conflict with OpenMC. The baseline
   precise_QA run needs no DESC (field map + surface are committed).

## Honest scope — what these runs DO and DON'T validate
- **DO:** that the pipeline runs end-to-end and produces sane per-mode TBR / coil
  flux / heating + a valid conformal DAGMC geometry on the real precise_QA field —
  i.e. the **machinery**.
- **DON'T (yet):** the free-streaming **numbers** cross-check vs the analytic
  stellarator NWL — that is blocked on the analytic side (it currently runs a toy
  config, not precise_QA) and on the conformal wall, per `JOINT_VALIDATION.md`. So
  do not claim "validated against analytic" from these runs; claim "machinery
  validated, sane physics, geometry valid."

## Report back / hand back
Commit your findings to `spf-prototype` (a short `RESULTS_ginsburg_bringup.md`:
what passed, the tiny-run RESULTS numbers, any geometry/DAGMC issues + fixes, the
env/XS paths you used). The user can sync the branch back. If you hit the conformal
self-intersection, document the `--scale` (or layer thicknesses) that produced a
clean DAGMC model so the baseline uses it.
