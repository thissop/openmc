#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_wallmap
#SBATCH -p short
#SBATCH -t 0-01:30
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=0-4
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src; PP=$SPF/spf_pp
export PYTHONPATH="$SRC:${PYTHONPATH:-}"
export LD_PRELOAD="$SRC/build/lib/libopenmc.so"
export PATH="$SRC/build/bin:$PATH"
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=8
IDS=(803097 886079 932746 59509 1960314)
ID=${IDS[$SLURM_ARRAY_TASK_ID]}
cd "$PP/python"
echo "=== (theta,phi) wall-map transport: device $ID ==="
python realdev_wallmap.py "$ID" 2>&1 | tail -14
echo "DONE_WM_$ID"
