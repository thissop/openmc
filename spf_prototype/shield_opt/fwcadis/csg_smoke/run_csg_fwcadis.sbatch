#!/usr/bin/env bash
#SBATCH --job-name=fwc-csg
#SBATCH --partition=short
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=csg_fwcadis_%j.out
#SBATCH --error=csg_fwcadis_%j.err

set -euo pipefail

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"

openmc_source=""
for candidate in \
    /ginsburg/astro/users/tjk2147/spf_work/openmc_src \
    /burg-archive/home/tjk2147/src/GitHub/openmc-spf \
    /ginsburg/astro/users/tjk2147/spf_work/openmc; do
    if [[ -x "$candidate/build/bin/openmc" && -f "$candidate/openmc/__init__.py" ]]; then
        openmc_source="$candidate"
        break
    fi
done

if [[ -z "$openmc_source" ]]; then
    printf 'ERROR: no existing source-built OpenMC installation was found\n' >&2
    exit 2
fi

export PYTHONNOUSERSITE=1
export PYTHONPATH="$openmc_source"
export LD_PRELOAD="$openmc_source/build/lib/libopenmc.so"
export PATH="$openmc_source/build/bin:$PATH"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# The protocol submits this script from its configured scratch root, so Slurm's
# submit directory is the root rather than this file's directory.
job_dir="$SLURM_SUBMIT_DIR/openmc_fwcadis/csg_smoke"
run_dir="$job_dir/run_${SLURM_JOB_ID}"
mkdir -p "$run_dir"

printf 'HOST=%s\n' "$(hostname)"
printf 'OPENMC_SOURCE=%s\n' "$openmc_source"
printf 'OPENMC_EXEC=%s\n' "$(command -v openmc)"
openmc --version
python "$job_dir/csg_fwcadis_smoke.py" \
    --workdir "$run_dir" \
    --threads "$SLURM_CPUS_PER_TASK"
