#!/bin/bash
#SBATCH -A astro                            # <-- your Ginsburg account (from your example)
#SBATCH --job-name=spf_conformal
#SBATCH --output=spf_conformal_%j.out       # written to $SLURM_SUBMIT_DIR (persistent)
#SBATCH --error=spf_conformal_%j.err
#SBATCH -t 0-08:00                          # OpenMC CPU job (no GPU); raise for production
#SBATCH -N 1
#SBATCH -c 16                               # OpenMC OpenMP threads
#SBATCH --mem-per-cpu=4G
#
# Turnkey OFFLINE conformal stellarator SPF run. Submit with:
#     sbatch ginsburg_job.sh
# Everything network-dependent (conda env, installs, cross sections, field maps for
# NEW equilibria) is done ONCE on the login node first -- see SETUP_GINSBURG.md.
# This job touches the network for NOTHING. CPU only -- no CUDA / LD_LIBRARY_PATH /
# XLA hacks: the whole OpenMC+DAGMC stack is self-contained in the conda-forge env.
set -euo pipefail
echo "==== SLURM job on $(hostname); $(date) ===="

# --- 1. environment (system anaconda module -> conda env; your Ginsburg pattern) ---
module load anaconda/3-2023.09                          # <-- adjust to the available version
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator                          # created on the login node (env yml)
REPO="/burg/home/$USER/src/GitHub/openmc"               # <-- path to the cloned repo
export OPENMC_CROSS_SECTIONS="/burg/home/$USER/endf_b_viii/cross_sections.xml"  # <-- your XS

# --- 2. threads + persistent output (NOT /tmp: it is node-local on SLURM) ---
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
RESULTS="${SLURM_SUBMIT_DIR:-$PWD}/results_qa_${SLURM_JOB_ID:-local}"

# --- 3. sanity: fail fast (in the job, before the long run) ---
python - <<'PY'
import openmc, os
assert hasattr(openmc, "DAGMCUniverse"), "OpenMC has no DAGMC support (wrong build/env)"
xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
assert xs and os.path.exists(xs), f"OPENMC_CROSS_SECTIONS not found: {xs!r}"
print("preflight OK: DAGMC-enabled OpenMC + cross sections present")
PY

# --- 4. the run (single offline entrypoint; bootstraps .so/geometry/.h5m offline) ---
cd "$REPO"
python spf_prototype/python/run_ginsburg.py \
    --stem equil_precise_qa \
    --scale 10 \
    --particles 2000000 \
    --batches 20 \
    --results-dir "$RESULTS"

echo "DONE -> $RESULTS/RESULTS_tier8_conformal.md"
