#!/bin/bash
#SBATCH --job-name=pstlbld
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=03:59:00
set -x
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl
export PYTHONNOUSERSITE=1
cd /burg-archive/home/tjk2147/pstl_test/ww_pert
echo "START BUILD AMP=$AMP OUT=$OUTNAME $(date)"
python build_perturbed.py $AMP $OUTNAME
echo "END BUILD $(date)"
