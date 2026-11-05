#!/usr/bin/env bash
#SBATCH --job-name=ssr-vr-analysis
#SBATCH --partition=short
#SBATCH --account=astro
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=ssr_analysis_%j.out
#SBATCH --error=ssr_analysis_%j.err

set -euo pipefail

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PYTHONNOUSERSITE=1
export HDF5_USE_FILE_LOCKING=FALSE

job_root="${BENCH_ROOT:-$SLURM_SUBMIT_DIR/openmc_fwcadis/ssr_benchmark}"
label="${BENCH_LABEL:?BENCH_LABEL must name a completed result set}"
result_root="$job_root/results/$label"

python -u "$job_root/analyze_benchmark.py" \
    "$result_root" \
    --master /burg-archive/home/tjk2147/pstl_test/corrected/coil_op_master.npz \
    --surface-source /burg-archive/home/tjk2147/pstl_test/corrected/surf_op_master.h5 \
    --output-dir "$result_root/analysis"
