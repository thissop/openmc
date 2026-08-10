#!/bin/bash
# Step (1) fix: build + dose the CONCENTRATE multi-coil placement (delta_multi_greedy) at the SAME
# 630.7 cm budget as the spread variant that lost. Greedy bang-bang dumps the budget on the
# collective inboard band (poloidal 120-200 deg), stellarator-symmetric across periods -> protects
# every coil's inboard sightline, testing whether concentrate beats both the spread multicoil
# (peak 2.33e-2) and the single-coil placed (peak 8.28e-3, which let the peak migrate coil20->25).
#
#   sbatch slurm_multicoil_greedy.sh
#SBATCH -J mcgreedy
#SBATCH -p burst
#SBATCH --account=astro
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G
#SBATCH --exclude=g002,g074,g080,g197,g088,g179,g232,g267
#SBATCH -o /burg-archive/home/tjk2147/pstl_test/corrected/slurm_mcgreedy_%j.log
set -uo pipefail
WORK=/burg-archive/home/tjk2147/pstl_test/corrected
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
RUNDIR=$WORK/mcgreedyrun; mkdir -p $RUNDIR
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh

echo "=== build multicoil_greedy DAGMC (pstl env) $(date) ==="
conda activate /ginsburg/astro/users/tjk2147/spf_work/envs/pstl
export PYTHONNOUSERSITE=1; unset LD_PRELOAD PYTHONPATH
cd $WORK && python -u build_step1.py multicoil_greedy || { echo BUILD_FAIL; exit 1; }

echo "=== magnet cells $(date) ==="
python -u step1_cells.py $WORK/dagmc_qh_step1_multicoil_greedy.h5m $WORK/cells_step1_multicoil_greedy.npz \
    || { echo CELLS_FAIL; exit 1; }
if [ "${SKIP_DOSE:-0}" = "1" ]; then echo BUILD_ONLY_DONE; exit 0; fi

echo "=== dose ALL coils (isolated, no-WW 48M) $(date) ==="
conda deactivate; conda activate spf-stellarator
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC LD_PRELOAD=$SRC/build/lib/libopenmc.so PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=32 PYTHONUNBUFFERED=1 DOSE_CWD=$RUNDIR
cd $RUNDIR
python -u $WORK/coil_run_iso.py unpol 12 4000000 $WORK/dagmc_qh_step1_multicoil_greedy.h5m \
    step1b_multicoil_greedy 60 $WORK/cells_step1_multicoil_greedy.npz || echo "rc=$? (maybe benign)"
[ -f $RUNDIR/coil_step1b_multicoil_greedy.npz ] \
    && mv -f $RUNDIR/coil_step1b_multicoil_greedy.npz $WORK/ && echo "DONE_multicoil_greedy $(date)" \
    || { echo DOSE_TRULY_FAILED; exit 1; }
