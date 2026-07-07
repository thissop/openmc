#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_capstone
#SBATCH -p short
#SBATCH -t 0-01:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src; PP=$SPF/spf_pp
export PYTHONPATH="$SRC:${PYTHONPATH:-}"                # PR openmc (native SPF TokamakSource)
export LD_PRELOAD="$SRC/build/lib/libopenmc.so"        # built lib (conda's is first in RPATH)
export PATH="$SRC/build/bin:$PATH"
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=8
python -c "import openmc; print('openmc', openmc.__version__, '->', openmc.__file__)" 2>&1 | head -2
cd "$PP/python"
echo "=== Schwartz capstone: native SPF TokamakSource, square torus, midplane NWL ==="
python schwartz_capstone.py 2>&1 | tail -30
echo DONE_CAPSTONE
