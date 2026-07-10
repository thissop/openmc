#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_cs59509
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src; PP=$SPF/spf_pp
export PYTHONPATH="$SRC:${PYTHONPATH:-}"; export LD_PRELOAD="$SRC/build/lib/libopenmc.so"
export PATH="$SRC/build/bin:$CONDA_PREFIX/bin:$PATH"
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=8
cd "$PP/python"
echo "=== conformal + scattering, device 59509, StellaratorSource on multi-layer blanket ==="
python conformal_scatter_59509.py
echo "EXIT_CODE=$? (benign double-free at exit is expected; check for SAVED above)"
ls -la "$PP/data/quasr59509_conformal_scatter.npz"
