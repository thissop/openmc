#!/bin/bash
# Patch-worth experiment driver (Ginsburg). Two modes:
#   MODE=fwd   sbatch slurm_patch_worth.sh            -> Tier-0 forward mesh flux (once)
#   MODE=patch sbatch --array=0-N slurm_patch_worth.sh -> Tier-2 per-patch build + coil dose
# Lessons baked in: burst partition (14-day walltime, kill-shot timed out on short=12h);
# per-array-task cwd so parallel builds/HDF5 don't race; FW-CADIS weight windows so the
# deep coil dose is affordable; correlated seed across patches (SEED fixed) so Delta R has
# low variance vs the baseline.
#SBATCH --job-name=patchworth
#SBATCH --partition=burst
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=slurm_patchworth_%A_%a.log

set -euo pipefail
MODE="${MODE:-patch}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"

WORK=/burg-archive/home/tjk2147/pstl_test/corrected
PW=~/src/GitHub/openmc-spf/spf_prototype/shield_opt/patch_worth   # or wherever synced
OMC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src            # transport tree (StellaratorSource)
SRCLIB=$OMC/build/libstellarator_source.so                      # adjust if named differently
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/cross_sections.xml

ADJ=~/adjoint_test/recip_w35f_coil20_flat/adjoint_importance_flat_P0.npz
WW=$WORK/weight_windows.h5           # FW-CADIS ww (baseline); regenerate if over-split
BASE_DAGMC=$WORK/dagmc_qh_step1_uniform.h5m   # baseline geometry (delta=0 ~ uniform ref)
MATS=$WORK/materials.xml
SEED=1                                # SAME seed across all patches (correlated sampling)

N_TOR=4; N_POL=4; DELTA_CM=20.0

if [ "$MODE" = "fwd" ]; then
  cd "$WORK"
  echo "=== Tier-0: forward mesh flux (co-located with adjoint) ==="
  python "$PW/fwd_meshflux.py" \
    --adjoint "$ADJ" --dagmc "$BASE_DAGMC" --materials "$MATS" \
    --source-lib "$SRCLIB" --source-params "emissivity=uniform" \
    --weight-windows "$WW" --particles 4000000 --batches 20 --seed "$SEED" \
    --out "$WORK/qh_fwd_meshflux.npz"
  echo "=== DONE fwd meshflux -> $WORK/qh_fwd_meshflux.npz ==="
  exit 0
fi

# ---- Tier-2: one patch per array task. Map SLURM_ARRAY_TASK_ID -> (I_TOR,J_POL) from
#      a patch list file (patches.txt: one "I J" per line, chosen from tier1 corners). ----
IDX="${SLURM_ARRAY_TASK_ID:-0}"
read -r I_TOR J_POL < <(sed -n "$((IDX+1))p" "$PW/patches.txt")
TAG="t${I_TOR}p${J_POL}"
RUNDIR="$WORK/patchrun_${TAG}"
mkdir -p "$RUNDIR"; cd "$RUNDIR"                 # per-task cwd -> no HDF5 race

echo "=== Tier-2 patch ($I_TOR,$J_POL): build DAGMC ==="
python "$PW/build_patch.py" "$I_TOR" "$J_POL" "$N_TOR" "$N_POL" "$DELTA_CM"

echo "=== Tier-2 patch ($I_TOR,$J_POL): coil-20 dose (WW, seed=$SEED) ==="
# reuse the debugged coil-dose runner on the patch DAGMC; tally cell-20 fast flux.
python "$WORK/coil_run_v3.py" \
    --dagmc "$WORK/dagmc_qh_patch_${TAG}.h5m" \
    --materials "$MATS" --weight-windows "$WW" \
    --source-lib "$SRCLIB" --source-params "emissivity=uniform" \
    --particles 4000000 --batches 20 --seed "$SEED" \
    --out "$WORK/coil_patch_${TAG}.npz"
echo "=== DONE patch ($I_TOR,$J_POL) -> coil_patch_${TAG}.npz ==="
