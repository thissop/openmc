#!/bin/bash
#SBATCH -J s1build
#SBATCH -t 03:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH -o /burg-archive/home/tjk2147/pstl_test/corrected/slurm_s1build_%j.log
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl
export PYTHONNOUSERSITE=1
cd /burg-archive/home/tjk2147/pstl_test/corrected
echo "=== BUILD placed  $(date) ==="; python -u build_step1.py placed  && echo "OK_PLACED"  || echo "FAIL_PLACED"
echo "=== BUILD uniform $(date) ==="; python -u build_step1.py uniform && echo "OK_UNIFORM" || echo "FAIL_UNIFORM"
echo "=== STEP1_BUILDS_DONE $(date) ==="
ls -la dagmc_qh_step1_placed.h5m dagmc_qh_step1_uniform.h5m
