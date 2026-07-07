#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_vmecflux
#SBATCH -p short
#SBATCH -t 0-00:30
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
module load anaconda/3-2023.09 2>/dev/null
export OMP_NUM_THREADS=8
cd /ginsburg/astro/users/tjk2147/spf_work/spf_pp/python
"$HOME/vmec_venv/bin/python" vmec_fluxmap.py 803097 2>&1 | tail -20
echo DONE_VMECFLUX
