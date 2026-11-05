#!/usr/bin/env bash
#SBATCH --job-name=ssr-vr-pilot
#SBATCH --partition=short
#SBATCH --account=astro
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --array=0-3
#SBATCH --output=ssr_benchmark_%A_%a.out
#SBATCH --error=ssr_benchmark_%A_%a.err

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
export OPENMC_SOURCE_LIB="$openmc_source/build/lib/libopenmc.so"
unset LD_PRELOAD
export PATH="$openmc_source/build/bin:$PATH"
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

job_root="${BENCH_ROOT:-$SLURM_SUBMIT_DIR/openmc_fwcadis/ssr_benchmark}"
input_root="$job_root/inputs"
label="${BENCH_LABEL:-pilot}"

case "$SLURM_ARRAY_TASK_ID" in
    0) mode=analog; weight_windows= ;;
    1) mode=magic; weight_windows="$input_root/magic.h5" ;;
    2) mode=native; weight_windows="$input_root/native_fwcadis.h5" ;;
    3) mode=reciprocity; weight_windows="$input_root/reciprocity_cadis.h5" ;;
    *) printf 'ERROR: unsupported array index %s\n' "$SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

run_dir="$job_root/results/$label/${mode}_${SLURM_JOB_ID}"
mkdir -p "$run_dir"

weight_args=()
if [[ -n "$weight_windows" ]]; then
    weight_args=(--weight-windows "$weight_windows")
fi

printf 'HOST=%s\n' "$(hostname)"
printf 'MODE=%s RUN_DIR=%s\n' "$mode" "$run_dir"
LD_PRELOAD="$OPENMC_SOURCE_LIB" openmc --version
set +e
python -u "$job_root/benchmark_ssr.py" \
    --mode "$mode" \
    --workdir "$run_dir" \
    --dagmc /burg-archive/home/tjk2147/pstl_test/corrected/dagmc_qh_op_uniform.h5m \
    --surface-source /burg-archive/home/tjk2147/pstl_test/corrected/surf_op_master.h5 \
    --batches "${BENCH_BATCHES:-4}" \
    --particles "${BENCH_PARTICLES:-50000}" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --seed "${BENCH_SEED:-104729}" \
    --max-split "${BENCH_MAX_SPLIT:-4}" \
    --max-history-splits "${BENCH_MAX_HISTORY_SPLITS:-100}" \
    --max-lower-bound-ratio "${BENCH_MAX_LOWER_BOUND_RATIO:-5}" \
    "${weight_args[@]}"
runner_rc=$?
set -e

# The source-built Python extension can abort during interpreter teardown after
# a fully successful run. Treat that known teardown as benign only when the
# atomic result envelope was already written with PASS status.
result_json="$run_dir/result.json"
if [[ -f "$result_json" ]] && python -c \
    'import json,sys; raise SystemExit(json.load(open(sys.argv[1]))["status"] != "PASS")' \
    "$result_json"; then
    if [[ "$runner_rc" -ne 0 ]]; then
        printf 'WARNING: runner exited %s after writing a PASS result (known OpenMC teardown)\n' "$runner_rc" >&2
    fi
    exit 0
fi

if [[ "$runner_rc" -eq 0 ]]; then
    runner_rc=1
fi
exit "$runner_rc"
