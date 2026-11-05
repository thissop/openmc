#!/usr/bin/env bash
#SBATCH --job-name=recip-cadis
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=reciprocity_cadis_%j.out
#SBATCH --error=reciprocity_cadis_%j.err

set -euo pipefail

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator

openmc_source=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
if [[ ! -x "$openmc_source/build/bin/openmc" ]]; then
    printf 'ERROR: source-built OpenMC is missing at %s\n' "$openmc_source" >&2
    exit 2
fi

export PYTHONNOUSERSITE=1
export PYTHONPATH="$openmc_source"
export LD_PRELOAD="$openmc_source/build/lib/libopenmc.so"
export PATH="$openmc_source/build/bin:$PATH"
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

job_dir="$SLURM_SUBMIT_DIR/openmc_fwcadis/reciprocity_cadis"
run_dir="$job_dir/run_${SLURM_JOB_ID}"
mkdir -p "$run_dir"

printf 'HOST=%s\n' "$(hostname)"
printf 'OPENMC_SOURCE=%s\n' "$openmc_source"
openmc --version
python -u "$job_dir/run_reciprocity_cadis.py" \
    --workdir "$run_dir" \
    --dagmc /burg-archive/home/tjk2147/pstl_test/corrected/dagmc_qh_op_uniform.h5m \
    --coils /burg-archive/home/tjk2147/pstl_test/coils_qh \
    --batches "${RECIP_BATCHES:-20}" \
    --particles "${RECIP_PARTICLES:-20000}" \
    --mesh "${RECIP_MESH:-16,16,8}" \
    --threads "$SLURM_CPUS_PER_TASK"
