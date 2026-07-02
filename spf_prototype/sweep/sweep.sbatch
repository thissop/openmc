#!/bin/bash
#SBATCH -A astro                              # Ginsburg account
#SBATCH --job-name=spf_sweep
#SBATCH --output=spf_sweep_%A_%a.out          # -> $SLURM_SUBMIT_DIR (persistent)
#SBATCH --error=spf_sweep_%A_%a.err
#SBATCH -p burst                              # 14-day; a C-scan can be long
#SBATCH -t 1-00:00
#SBATCH -N 1
#SBATCH -c 16                                 # OpenMP threads per config
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-0%8                          # <-- set 0-(Nconfigs-1); %8 = 8 concurrent
#
# SLURM ARRAY sweep: one array task per config. Submit with the array range set to
# the number of configs in the manifest, e.g. for 24 configs, 8 at a time:
#     sbatch --array=0-23%8 sweep.sbatch
# All network steps (conda env, cross sections, the .so PREBUILD) are done ONCE
# beforehand -- see docs/README.md runbook. This job touches the network for nothing.
set -euo pipefail
echo "==== array task ${SLURM_ARRAY_TASK_ID:-?} on $(hostname); $(date) ===="

# --- environment (exact activation; no silent base-python fallback) ----------
unset PYTHONPATH                                        # avoid a stray PYTHONPATH shadow
module load anaconda/3-2023.09                          # adjust to the available version
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"                   # env python wins, not base
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

REPO="${SPF_REPO:-$HOME/src/GitHub/openmc}"             # <-- edit to your clone path
export OPENMC_CROSS_SECTIONS="${OPENMC_CROSS_SECTIONS:-$HOME/endfb-viii.0-hdf5/cross_sections.xml}"
MANIFEST="${MANIFEST:-$REPO/spf_prototype/sweep/configs/production.json}"
OUT="${SLURM_SUBMIT_DIR:-$PWD}/sweep_out"

# neutral cwd: the repo root contains an openmc/ SOURCE tree; running from there
# would let it shadow the conda openmc. Run from the submit dir instead.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# --- preflight (fail fast, per task): conda openmc + DAGMC + XS + the .so ------
python - <<'PY'
import openmc, os, sys
assert hasattr(openmc, "DAGMCUniverse"), "OpenMC has no DAGMC support (wrong env)"
assert "site-packages" in openmc.__file__, f"imported the repo fork, not conda: {openmc.__file__}"
xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
assert xs and os.path.exists(xs), f"OPENMC_CROSS_SECTIONS not found: {xs!r}"
so = os.path.join(os.environ["SPF_REPO"] if "SPF_REPO" in os.environ else
                  os.path.expanduser("~/src/GitHub/openmc"),
                  "spf_prototype/src/build/libpolarized_fusion_source.so")
assert os.path.exists(so), f"prebuild the source .so first (missing {so})"
print("preflight OK:", openmc.__file__)
PY

# --- run this task's config ---------------------------------------------------
python "$REPO/spf_prototype/sweep/sweep.py" \
    --manifest "$MANIFEST" \
    --out "$OUT" \
    --index "${SLURM_ARRAY_TASK_ID:-0}"

echo "==== task ${SLURM_ARRAY_TASK_ID:-?} done; record in $OUT ===="
