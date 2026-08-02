#!/bin/bash
# Patch-worth experiment driver (Ginsburg). REWRITTEN to the REAL cluster recipe
# (coil_run_v3.py positional interface + StellaratorSource + two conda envs), not the
# CompiledSource/--flag assumptions of the first draft.
#
#   MODE=fwd   sbatch slurm_patch_worth.sh                 -> Tier-0 forward mesh flux (once)
#   MODE=patch sbatch --array=0-N slurm_patch_worth.sh     -> Tier-2 per-patch build + coil-20 dose
#
# Lessons baked in: burst partition (14-day walltime); per-array-task cwd so parallel builds
# don't race; FW-CADIS weight windows so the deep coil dose is affordable; correlated seed
# (coil_run_v3 uses a fixed internal seed given fixed args) so Delta R has low variance.
#SBATCH --job-name=patchworth
#SBATCH --partition=burst
#SBATCH --account=astro
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/burg-archive/home/tjk2147/pstl_test/corrected/patch_worth/slurm_pw_%A_%a.log

set -uo pipefail
MODE="${MODE:-patch}"

WORK=/burg-archive/home/tjk2147/pstl_test/corrected
PW=$WORK/patch_worth
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src   # transport tree (StellaratorSource)
ADJ=/burg-archive/home/tjk2147/adjoint_test/recip_w35f_coil20_flat/adjoint_importance_flat_P0.npz
# NO weight windows: the fwcadis WW over-splits on this geom (2M primaries -> 68M+ secondary
# tracks, batch 1 never finishes). The kill-shot proved brute-force no-WW resolves the deep
# coil dose to 0.5-2.4% relerr, so we match it exactly -> the patch doses become CORRELATED
# with the existing kill-shot uniform baseline (coil_step1b_uniform.npz), same seed=1.
WW=""
BASE_DAGMC=$WORK/dagmc_qh_step1_uniform.h5m            # baseline (delta=0 reference)
LI6=60
N_TOR=4; N_POL=4; DELTA_CM=20.0
BATCHES=12; PARTICLES=4000000                          # match kill-shot (48M no-WW, seed=1)
FWD_BATCHES=12                                         # 48M no-WW (~16h, safe under 24h wall)

act_transport() {
  source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
  conda activate spf-stellarator
  export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC
  export LD_PRELOAD=$SRC/build/lib/libopenmc.so
  export PATH=$SRC/build/bin:$PATH
  export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" PYTHONUNBUFFERED=1
}
act_build() {
  source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
  conda deactivate 2>/dev/null || true
  conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl
  export PYTHONNOUSERSITE=1
  unset LD_PRELOAD PYTHONPATH
}

if [ "$MODE" = "fwd" ]; then
  cd "$WORK"
  act_transport
  echo "=== Tier-0: forward mesh flux (co-located with adjoint, NO WW) $(date) ==="
  python -u "$PW/fwd_meshflux.py" "$FWD_BATCHES" "$PARTICLES" "$BASE_DAGMC" "$ADJ" \
      "$WORK/qh_fwd_meshflux.npz" "$LI6"
  echo "=== DONE fwd meshflux rc=$? -> $WORK/qh_fwd_meshflux.npz $(date) ==="
  exit 0
fi

if [ "$MODE" = "base" ]; then
  # Tier-2 finite-difference BASELINE: uniform (delta=0) coil-20 dose with the SAME
  # settings (no-WW, 12x4M, seed) as the patches -> apples-to-apples Delta R.
  cd "$WORK"
  act_transport
  # OPTIONAL: usually skip -- coil_step1b_uniform.npz (kill-shot, 12x4M no-WW seed=1) is already
  # a correlated baseline at identical settings. Run this only to refresh it.
  echo "=== Tier-2 baseline: uniform coil-20 dose (NO WW) $(date) ==="
  python -u "$WORK/coil_run_v3.py" unpol "$BATCHES" "$PARTICLES" \
      "$BASE_DAGMC" "patchbase" "$LI6" "$WORK/cells_step1_uniform.npz"
  echo "=== DONE baseline rc=$? -> coil_patchbase.npz $(date) ==="
  exit 0
fi

if [ "$MODE" = "doseonly" ]; then
  # Re-dose an ALREADY-BUILT patch DAGMC in an ISOLATED cwd (coil_run_iso.py chdir's to
  # $DOSE_CWD) so concurrent doses do NOT collide on OpenMC's default statepoint.12.h5 in
  # the shared WORK -- the race that made two different patches return identical flux.
  IDX="${SLURM_ARRAY_TASK_ID:-0}"
  read -r I_TOR J_POL < <(sed -n "$((IDX+1))p" "$PW/patches_dose.txt")
  TAG="t${I_TOR}p${J_POL}"
  RUNDIR="$WORK/doserun_${TAG}"; mkdir -p "$RUNDIR"
  export DOSE_CWD="$RUNDIR"
  act_transport
  echo "=== doseonly ($I_TOR,$J_POL): isolated coil-20 dose (NO WW) $(date) ==="
  # coil_run_iso can exit non-zero on the benign teardown double-free -> key on the npz, not rc.
  python -u "$WORK/coil_run_iso.py" unpol "$BATCHES" "$PARTICLES" \
      "$WORK/dagmc_qh_patch_${TAG}.h5m" "patch_${TAG}" "$LI6" \
      "$WORK/cells_patch_${TAG}.npz" || echo "coil_run rc=$? (may be benign exit double-free)"
  if [ -f "$RUNDIR/coil_patch_${TAG}.npz" ]; then
    mv -f "$RUNDIR/coil_patch_${TAG}.npz" "$WORK/"
    echo "=== DONE doseonly ($I_TOR,$J_POL) -> coil_patch_${TAG}.npz $(date) ==="
  else
    echo "DOSE_TRULY_FAILED ($I_TOR,$J_POL)"; exit 1
  fi
  exit 0
fi

# ---- Tier-2: one patch per array task. SLURM_ARRAY_TASK_ID -> line of patches.txt ("I J") ----
IDX="${SLURM_ARRAY_TASK_ID:-0}"
read -r I_TOR J_POL < <(sed -n "$((IDX+1))p" "$PW/patches.txt")
TAG="t${I_TOR}p${J_POL}"
RUNDIR="$WORK/patchrun_${TAG}"
mkdir -p "$RUNDIR"

echo "=== Tier-2 patch ($I_TOR,$J_POL): build DAGMC (pstl env, isolated to RUNDIR) $(date) ==="
cd "$RUNDIR"
act_build
# export to RUNDIR (6th arg) -> parallel builds don't clobber shared STEP/gmsh intermediates
python -u "$PW/build_patch.py" "$I_TOR" "$J_POL" "$N_TOR" "$N_POL" "$DELTA_CM" "$RUNDIR" \
    || { echo "BUILD_FAIL"; exit 1; }
# move the two final products to WORK (uniquely named per patch)
mv -f "$RUNDIR/dagmc_qh_patch_${TAG}.h5m" "$RUNDIR/delta_qh_patch_${TAG}.npz" "$WORK/" \
    || { echo "MOVE_FAIL"; exit 1; }

echo "=== Tier-2 patch ($I_TOR,$J_POL): magnet cells (pymoab) $(date) ==="
python -u "$WORK/step1_cells.py" "$WORK/dagmc_qh_patch_${TAG}.h5m" "$WORK/cells_patch_${TAG}.npz" \
    || { echo "CELLS_FAIL"; exit 1; }

echo "=== Tier-2 patch ($I_TOR,$J_POL): coil-20 dose (NO WW, transport env) $(date) ==="
cd "$RUNDIR"                                            # per-task cwd -> no HDF5 race
act_transport
python -u "$WORK/coil_run_v3.py" unpol "$BATCHES" "$PARTICLES" \
    "$WORK/dagmc_qh_patch_${TAG}.h5m" "patch_${TAG}" "$LI6" \
    "$WORK/cells_patch_${TAG}.npz" || { echo "DOSE_FAIL"; exit 1; }
# coil_run_v3 writes coil_patch_${TAG}.npz in cwd -> move to WORK for analysis
mv -f "coil_patch_${TAG}.npz" "$WORK/" 2>/dev/null || true
echo "=== DONE patch ($I_TOR,$J_POL) -> coil_patch_${TAG}.npz $(date) ==="
