#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_vmecbatch
#SBATCH -p short
#SBATCH -t 0-00:30
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=0-9
set -uo pipefail
module load anaconda/3-2023.09 2>/dev/null
export OMP_NUM_THREADS=8
IDS=(803097 923602 932746 59509 886079 1960314 2588494 2630694 2632243 1190023)
ID=${IDS[$SLURM_ARRAY_TASK_ID]}
cd /ginsburg/astro/users/tjk2147/spf_work/spf_pp/python
echo "=== VMEC fluxmap device $ID ==="
"$HOME/vmec_venv/bin/python" vmec_fluxmap.py "$ID" 2>&1 | tail -8
echo "DONE_$ID"
