# SETUP_GINSBURG.md — turnkey offline SLURM runs on Columbia Ginsburg

This is the **operational** runbook (do this on Ginsburg; you cannot use Claude
Code there). `RUN_ON_GINSBURG.md` has the conceptual pipeline; this file is the
copy-paste sequence for a no-internet SLURM batch run.

## The offline model (what runs where)

SLURM **compute nodes have no internet**. So every network step happens **once on
the login node**; the **batch job is 100% offline**.

| Phase | Where | Internet? | What |
|---|---|---|---|
| Setup | login node | yes | conda env, cross sections, build the `.so`, smoke test |
| (scan only) pre-gen | login node | yes* | field maps for NEW equilibria via DESC |
| Run | compute node (sbatch) | **no** | geometry → DAGMC → transport → metrics → RESULTS |

\* DESC's example equilibria (`precise_QA`, `precise_QH`) ship **inside** the
package, so even DESC needs no live internet once installed — but do it on the
login node anyway. **The committed `precise_QA` field map + geometry mean the
baseline run needs no DESC at all.**

---

## 0.5 Storage: put everything in GROUP space, not home

The two heavy items are the **conda env (~5–8 GB)** and the **cross sections
(~1–10 GB)** — these, not the outputs, are what blow a small home quota. Per-run
output is ~0.1 GB; a 10-equilibrium scan ~1 GB. So **work in your `astro` group
allocation** (`/burg/astro/...` — confirm your writable path), not `/burg/home`.

```bash
checkquota                      # confirm home vs group quota (or gpfsquota/mmlsquota)
WORK=/burg/astro/users/$USER    # <-- your writable group dir (adjust to your allocation)
```

- **conda env in group space** (so it doesn't eat home quota): create/activate by
  PATH — `conda env create -p $WORK/envs/spf-stellarator -f spf_prototype/environment.yml`
  then `conda activate $WORK/envs/spf-stellarator`.
- **cross sections**: reuse a SHARED cluster ENDF/B library if one exists (ask RCS /
  your group) and just point `OPENMC_CROSS_SECTIONS` at it — don't download your own
  multi-GB copy. Otherwise put it under `$WORK`.
- **repo + job output under `$WORK`**, and submit from there so `$SLURM_SUBMIT_DIR`
  (where results land) is group space. Set `REPO=$WORK/.../openmc` in
  `ginsburg_job.sh`.
- statepoints already avoid node-local `/tmp` (purged at job end) — they go to the
  persistent `--results-dir`.

## 1. Login-node setup (once; needs internet)

> **Library note (vs your GPU jobs):** OpenMC+DAGMC needs **no CUDA modules and no
> `LD_LIBRARY_PATH`/`XLA` hacks** — the entire stack (openmc, dagmc, moab, libs)
> comes self-contained from one conda-forge env. You only use the system
> `anaconda` module to get `conda`, then activate the env.

```bash
# clone (or copy) the repo into your /burg working space
cd /burg/home/$USER/src/GitHub && git clone <your_repo_url> openmc && cd openmc
git checkout spf-prototype

# get conda from the system anaconda module (your Ginsburg pattern)
module load anaconda/3-2023.09                                  # adjust to available version
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh

# (a) transport env (OpenMC+DAGMC + analytic stack). ~10-20 min, login node (has net).
conda env create -f spf_prototype/environment.yml
conda activate spf-stellarator
python -c "import openmc; assert hasattr(openmc,'DAGMCUniverse'); print('DAGMC OK')"

# (b) cross sections: point at an ENDF/B HDF5 library already on Ginsburg, or
#     download one once here. Then EXPORT THE SAME PATH in ginsburg_job.sh.
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
python -c "import openmc,os; print('XS:', os.path.exists(os.environ['OPENMC_CROSS_SECTIONS']))"

# (c) build the compiled spin-polarized source against the CONDA OpenMC.
#     build_so() auto-uses $CONDA_PREFIX; this just pre-builds it so the job
#     doesn't spend time (or hit a surprise) compiling.
cd spf_prototype/src && cmake -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" . && cmake --build build && cd ../..
ls spf_prototype/src/build/libpolarized_fusion_source.so   # must exist

# (d) smoke test the architecture-independent gates (sampler, analytic, geometry):
python -m pytest spf_prototype/tests -q          # expect: all pass

# (e) ONE tiny end-to-end DAGMC+transport check (catches stl_to_h5m / DAGMC / XS
#     issues in minutes, before submitting an 8-hour job). Runs on the login node:
OMP_NUM_THREADS=4 python spf_prototype/python/run_ginsburg.py \
    --stem equil_precise_qa --scale 10 --particles 20000 --batches 5 \
    --results-dir /tmp/spf_smoke
#   EXPECT: it builds <stem>.h5m, runs 6 configs, writes
#   /tmp/spf_smoke/RESULTS_tier8_conformal.md . If this fails, fix it HERE
#   (login node), not in a queued job. Common fixes: stl_to_h5m API drift
#   (build_dagmc.py single call site), XS path, conda env.
```

If step (e) writes a RESULTS file, the pipeline is good end-to-end and a full
batch run is just the same command with more particles.

---

## 2. Submit the baseline run (offline batch job)

Edit the env block at the top of `ginsburg_job.sh` (the `anaconda` module version,
its `profile.d` path, `REPO`, `OPENMC_CROSS_SECTIONS`; `-A astro` is already set),
then:

```bash
cd /burg/home/$USER/src/GitHub/openmc/spf_prototype
sbatch ginsburg_job.sh
```

Output lands in `$SLURM_SUBMIT_DIR/results_qa_<jobid>/`:
`RESULTS_tier8_conformal.md` (TBR per mode, coil fast-flux per mode free vs
scatter, the steering-survival η, per-material heating) + `statepoints/` (the raw
OpenMC `.h5` for any further postprocessing, incl. the φ-resolved wall mesh).
Statepoints are written to this **persistent** dir, never `/tmp`.

Bump `--particles` for production. The coil tally is deep behind the
blanket+shield, so its η is indicative without variance reduction — for a precise
coil η, add weight windows (a documented next step).

---

## 3. The QA→QH quasisymmetry scan (the science sweep)

The headline result correlates the steering benefit η with quasisymmetry. Each
equilibrium needs a field map, which uses DESC — so **pre-generate all field maps
on the login node**, then submit an array job that is fully offline.

Stems must be **lowercase** (`equil_precise_qa`) to match the committed files on
case-sensitive `/burg`; the DESC example NAME is mixed-case (`precise_QA`), so map
it down explicitly:

```bash
# (login node, desc env — DESC pins conflict with OpenMC, keep it separate)
python -m venv $HOME/desc_venv && $HOME/desc_venv/bin/pip install desc-opt
for EQ in precise_QA precise_QH ; do          # add QUASR devices for a denser scan
  STEM="equil_$(echo "$EQ" | tr 'A-Z' 'a-z')"   # -> equil_precise_qa / equil_precise_qh
  $HOME/desc_venv/bin/python spf_prototype/python/desc_to_fieldmap.py "$EQ" "$STEM"
  conda run -n spf-stellarator python spf_prototype/python/stellarator_geometry.py "$STEM" 10.0
done
# now data/equil_precise_q{a,h}.{meta,bin}, *_surface.npz, *_geom/ exist (lowercase)
```

Then an array job (one task per equilibrium): copy `ginsburg_job.sh` to
`ginsburg_scan.sh`, add `#SBATCH --array=0-1`, select the lowercase stem by index,
and **build the `.so` once on the login node first** (array tasks share one `REPO`,
so a concurrent `cmake --build` would race — `run_ginsburg` skips the build when the
`.so` already exists):

```bash
EQS=(equil_precise_qa equil_precise_qh)        # lowercase
STEM=${EQS[$SLURM_ARRAY_TASK_ID]}
python spf_prototype/python/run_ginsburg.py --stem "$STEM" --scale 10 \
    --particles 2000000 --batches 20 \
    --results-dir "$SLURM_SUBMIT_DIR/results_${STEM}_${SLURM_JOB_ID}"
```

> **Deploy note:** if you `rsync`/copy the working tree (instead of `git clone`),
> first `rm -rf spf_prototype/src/build` — a stale `CMakeCache.txt` from another
> machine breaks the `cmake` build. A fresh `git clone` is clean (`build/` + `*.so`
> are gitignored).

Note: QH has larger field-period excursion → larger plasma-shaping; if the
conformal build self-intersects, `stellarator_geometry.py` refuses it — raise
`--scale` (thinner build relative to minor radius) until it validates on the login
node before submitting.

---

## 4. Honest caveats (read before trusting numbers)

- **The x86 transport path was never executed on the aarch64 dev box** (no DAGMC
  there). Step 1(e) is the real bring-up gate — run it before any long job. The
  version-sensitive spots are isolated and flagged (`build_dagmc.py` `stl_to_h5m`
  call; the OpenMC build prefix).
- **Coil tallies** are statistics-starved without variance reduction; lead with
  the trend and the ± errors.
- **Geometry** is a conformal radial-build rig on `precise_QA` (a public QA
  equilibrium); **machinery/validation, not a Helios physics claim.** Device-
  specific values are `INJECT(helios)` (see `DATA_NEEDED.md`).
- **Joint analytic↔MC cross-check** (free-streaming OpenMC vs the analytic NWL) is
  a postprocessing step that needs your companion analytic code; the run saves the
  free-streaming statepoints for it. See `JOINT_VALIDATION.md`.
