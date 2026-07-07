#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_rebuild
#SBATCH -p short
#SBATCH -t 0-01:00
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=4G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src
echo "=== DAGMC/MOAB in conda env? ==="
ls $CONDA_PREFIX/lib/libdagmc* $CONDA_PREFIX/lib/libMOAB* 2>&1 | head -6
find $CONDA_PREFIX -name "DAGMCConfig.cmake" -o -name "MOABConfig.cmake" 2>/dev/null | head
echo "=== current build DAGMC flag ==="
grep -i "DAGMC\|dagmc" $SRC/build/CMakeCache.txt 2>/dev/null | head
echo "=== reconfigure build with DAGMC=ON ==="
cd $SRC/build
cmake -DOPENMC_USE_DAGMC=ON -DOPENMC_USE_MPI=ON -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" .. 2>&1 | tail -18
echo "=== rebuild ==="
make -j16 2>&1 | tail -12
echo "=== verify DAGMC now compiled in ==="
grep -i "OPENMC_USE_DAGMC" CMakeCache.txt
echo DONE_REBUILD
