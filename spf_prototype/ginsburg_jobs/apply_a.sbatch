#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_apply_a
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
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src; NS=$SPF/spf_pp/native_spf

echo "=== apply native SPF edits onto the PR source ==="
cp "$NS/source.h"   "$SRC/include/openmc/source.h"
cp "$NS/source.cpp" "$SRC/src/source.cpp"
cp "$NS/source.py"  "$SRC/openmc/source.py"
grep -q spin_fractions_to_abc "$SRC/openmc/__init__.py" || \
  echo "from openmc.source import spin_fractions_to_abc" >> "$SRC/openmc/__init__.py"

echo "=== incremental rebuild (only source.cpp + relink) ==="
cd "$SRC/build" && make -j8 2>&1 | tail -14

echo "=== run env: PR python + built lib ==="
export PYTHONPATH="$SRC:${PYTHONPATH:-}"
export LD_PRELOAD="$SRC/build/lib/libopenmc.so"
export PATH="$SRC/build/bin:$PATH"
python -c "import openmc; print('openmc', openmc.__version__, '->', openmc.__file__)" 2>&1 | head -2

echo "=== SPF TokamakSource unit tests (physics-mirror + Python API) ==="
cd "$NS" && python -m pytest test_spf_tokamak.py -q 2>&1 | tail -22
echo DONE_APPLY_A
