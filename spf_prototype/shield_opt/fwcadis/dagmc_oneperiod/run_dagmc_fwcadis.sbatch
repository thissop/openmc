#!/usr/bin/env bash
#SBATCH --job-name=dag-fwcadis
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=dagmc_fwcadis_%j.out
#SBATCH --error=dagmc_fwcadis_%j.err

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

job_dir="$SLURM_SUBMIT_DIR/openmc_fwcadis/dagmc_oneperiod"
run_dir="$job_dir/run_${SLURM_JOB_ID}"
mgxs_path="$job_dir/mgxs_casmo2.h5"
mkdir -p "$run_dir"

printf 'HOST=%s\n' "$(hostname)"
printf 'OPENMC_SOURCE=%s\n' "$openmc_source"
openmc --version
python -u "$job_dir/dagmc_fwcadis.py" \
    --workdir "$run_dir" \
    --dagmc /burg-archive/home/tjk2147/pstl_test/corrected/dagmc_qh_op_uniform.h5m \
    --fluxmap /burg-archive/home/tjk2147/pstl_test/qh_fluxmap \
    --mgxs "$mgxs_path" \
    --mgxs-particles "${FWC_MGXS_PARTICLES:-200}" \
    --source-points "${FWC_SOURCE_POINTS:-2000}" \
    --rays "${FWC_RAYS:-2000}" \
    --batches "${FWC_BATCHES:-12}" \
    --inactive "${FWC_INACTIVE:-3}" \
    --mesh "${FWC_MESH:-10,10,10}" \
    --threads "$SLURM_CPUS_PER_TASK"
