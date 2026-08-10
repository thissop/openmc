#!/bin/bash
# SSR shielding run (reads a frozen surface source). Env: TAG DAGMC BATCHES PARTICLES SURF_SRC CELLS
#SBATCH -J ssrdose
#SBATCH -p burst
#SBATCH --account=astro
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=96G
#SBATCH -o /burg-archive/home/tjk2147/pstl_test/corrected/ssr/slurm_%x_%j.log
set -uo pipefail
WORK=/burg-archive/home/tjk2147/pstl_test/corrected
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
SSR=$WORK/ssr
export RUNDIR=$WORK/ssrrun_${TAG}; mkdir -p "$RUNDIR"
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh; conda activate spf-stellarator
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC LD_PRELOAD=$SRC/build/lib/libopenmc.so PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=32 PYTHONUNBUFFERED=1
echo "=== ssr_dose $TAG $(date) ==="
python -u "$SSR/ssr_dose.py" "$BATCHES" "$PARTICLES" "$WORK/$DAGMC" "$TAG" 60 "$WORK/$CELLS" "$SURF_SRC" \
  || echo "rc=$? (benign teardown?)"
[ -f "$WORK/coil_${TAG}.npz" ] && echo "OK_${TAG} $(date)" || { echo "FAILED_${TAG}"; exit 1; }
