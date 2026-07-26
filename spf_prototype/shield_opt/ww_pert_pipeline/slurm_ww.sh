#!/bin/bash
#SBATCH --job-name=wwgen
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:59:00
#SBATCH --output=/burg-archive/home/tjk2147/pstl_test/ww_pert/wwgen.log
set -x
SPF=/ginsburg/astro/users/tjk2147/spf_work
SRC=$SPF/openmc_src
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PYTHONPATH=$SRC:$PYTHONPATH
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
cd /burg-archive/home/tjk2147/pstl_test/ww_pert
echo "START WWGEN $(date)"
python gen_ww.py 12 100000 12
echo "END WWGEN $(date)"
