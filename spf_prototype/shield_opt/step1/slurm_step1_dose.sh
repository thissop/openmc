#!/bin/bash
#SBATCH -J s1dose
#SBATCH -t 04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --account=astro
#SBATCH --partition=short,burst,rent
#SBATCH -o /burg-archive/home/tjk2147/pstl_test/corrected/slurm_s1dose_%j.log
set -x
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
# CORRECT SPF openmc tree (has StellaratorSource) -- matches slurm_cwff35.sh, NOT ~/src/GitHub/openmc-spf
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=32 PYTHONUNBUFFERED=1
cd /burg-archive/home/tjk2147/pstl_test/corrected
# SAME source + settings for both; ONLY the DAGMC (shield placement) differs. unpol, Li6=60%.
for VAR in placed uniform; do
  DAG=/burg-archive/home/tjk2147/pstl_test/corrected/dagmc_qh_step1_${VAR}.h5m
  CELLS=/burg-archive/home/tjk2147/pstl_test/corrected/cells_step1_${VAR}.npz
  echo "######### STEP1 coil dose: $VAR  $(date) #########"
  [ -f "$DAG" ] || { echo "MISSING $DAG"; continue; }
  python -u coil_run_v3.py unpol 50 10000000 "$DAG" step1_${VAR} 60 "$CELLS"
  echo "DONE_${VAR} rc=$?"
done
echo "=== STEP1_DOSE_DONE $(date) ==="
