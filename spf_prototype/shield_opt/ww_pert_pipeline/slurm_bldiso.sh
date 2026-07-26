#!/bin/bash
#SBATCH --job-name=bldiso
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:59:00
#SBATCH --output=/burg-archive/home/tjk2147/pstl_test/kiso/bldiso.log
set -x
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl
export PYTHONNOUSERSITE=1
export WORKDIR=/burg-archive/home/tjk2147/pstl_test/kiso
cd $WORKDIR
echo "START BLDISO $(date)"
python build_pert_iso.py 25 perturbed
echo "END BLDISO $(date)"
