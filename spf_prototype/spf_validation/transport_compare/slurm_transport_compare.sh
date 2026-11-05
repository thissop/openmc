#!/bin/bash
#SBATCH -A astro
#SBATCH -J src-rep
#SBATCH -p short,burst,rent
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-7
#SBATCH -o /ginsburg/astro/users/tjk2147/advancedagents/spf_transport_compare/logs/%x_%A_%a.out
#SBATCH -e /ginsburg/astro/users/tjk2147/advancedagents/spf_transport_compare/logs/%x_%A_%a.err

set -euo pipefail

CASES=(
  qa_unpol_grid qa_unpol_spectral qa_A_grid qa_A_spectral
  qh_unpol_grid qh_unpol_spectral qh_A_grid qh_A_spectral
)
CASE=${CASES[$SLURM_ARRAY_TASK_ID]}
IFS=_ read -r DEVICE MODE REPRESENTATION <<< "$CASE"

REMOTE_ROOT=${GINSBURG_REMOTE_ROOT:-/ginsburg/astro/users/tjk2147/advancedagents}
WORK=${COMPARE_WORK:-$REMOTE_ROOT/spf_transport_compare}
GRID_SRC=${GRID_SRC:-/ginsburg/astro/users/tjk2147/spf_work/openmc_src}
SPECTRAL_SRC=${SPECTRAL_SRC:-$WORK/openmc-spf-pol}
DRIVER=${COMPARE_DRIVER:-$WORK/run_transport_compare.py}
PSTL_ENV=${PSTL_ENV:-/ginsburg/astro/users/tjk2147/spf_work/envs/pstl}
XS=${OPENMC_CROSS_SECTIONS:-/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml}
BATCHES=${COMPARE_BATCHES:-2}
PARTICLES=${COMPARE_PARTICLES:-10000}
SEED=${COMPARE_SEED:-20260901}

if [[ "$REPRESENTATION" == grid ]]; then
  SRC=$GRID_SRC
else
  SRC=$SPECTRAL_SRC
fi

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate "$PSTL_ENV"
export PYTHONNOUSERSITE=1
export PYTHONPATH=$SRC
export PATH=$SRC/build/bin:$PATH
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export OPENMC_CROSS_SECTIONS=$XS
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export PYTHONUNBUFFERED=1

if [[ "$DEVICE" == qa ]]; then
  DAGMC=${QA_DAGMC:-/burg-archive/home/tjk2147/pstl_test/qa/dagmc_qa_base.h5m}
  FLUXMAP=${QA_FLUXMAP:-/burg-archive/home/tjk2147/qa_reactor_fluxmap}
  WOUT=${QA_WOUT:-$WORK/wouts/wout_LandremanPaul2021_QA.nc}
  CELL_ARGS=(--coil-cells 11-26)
  # The staged canonical QA wout is already reactor-scaled (R spans roughly
  # 6.2--13.2 m), so from_vmec's m-to-cm conversion needs no extra scale.
  LENGTH_SCALE=${QA_LENGTH_SCALE:-1.0}
else
  DAGMC=${QH_DAGMC:-/burg-archive/home/tjk2147/pstl_test/corrected/dagmc_qh_step1_uniform.h5m}
  FLUXMAP=${QH_FLUXMAP:-/burg-archive/home/tjk2147/pstl_test/qh_fluxmap}
  WOUT=${QH_WOUT:-/burg-archive/home/tjk2147/pstl_test/wout_qh.nc}
  CELL_ARGS=(--cells-npz /burg-archive/home/tjk2147/pstl_test/corrected/cells_step1_uniform.npz)
  LENGTH_SCALE=${QH_LENGTH_SCALE:-1.0}
fi

OUT=$WORK/results/$CASE
mkdir -p "$OUT" "$WORK/logs"
python -u "$DRIVER" \
  --device "$DEVICE" \
  --representation "$REPRESENTATION" \
  --mode "$MODE" \
  --dagmc "$DAGMC" \
  "${CELL_ARGS[@]}" \
  --fluxmap "$FLUXMAP" \
  --wout "$WOUT" \
  --length-scale "$LENGTH_SCALE" \
  --batches "$BATCHES" \
  --particles "$PARTICLES" \
  --seed "$SEED" \
  --threads "$OMP_NUM_THREADS" \
  --output-dir "$OUT"
