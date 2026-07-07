#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_ompc_build
#SBATCH -p short
#SBATCH -t 0-01:30
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
cd "$SRC"
echo "=== configuring OpenMC (PR #3999) from source, DAGMC OFF (CSG suffices for source validation) ==="
rm -rf build && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPENMC_USE_DAGMC=OFF \
      -DOPENMC_USE_MPI=ON \
      -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
      -DHDF5_ROOT="$CONDA_PREFIX" \
      .. 2>&1 | tail -20
echo "=== building (make -j8) ==="
make -j8 2>&1 | tail -30
echo "=== result ==="
ls -la ./bin/openmc 2>/dev/null && ./bin/openmc --version 2>&1 | head -3
echo DONE_BUILD
