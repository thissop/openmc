#!/bin/bash
#SBATCH --job-name=cmpiso
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=07:59:00
#SBATCH --output=/burg-archive/home/tjk2147/pstl_test/kiso/compare.log
set -x
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PYTHONPATH=$SRC:$PYTHONPATH
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export WORKDIR=/burg-archive/home/tjk2147/pstl_test/kiso
cd $WORKDIR
BASE=$WORKDIR/dagmc_baseline_fine.h5m
PERT=$WORKDIR/dagmc_perturbed.h5m
WW=weight_windows.h5
B=${B:-20}; P=${P:-100000}
for f in $BASE $PERT; do [ -f "$f" ] || { echo "MISSING $f"; exit 2; }; done
echo "===== COMPARE base vs pert (WW), B=$B P=$P ====="
python coil_run_iso.py unpol $B $P $BASE base_unpol $WW
python coil_run_iso.py unpol $B $P $PERT pert_unpol $WW
python coil_run_iso.py A     $B $P $BASE base_A     $WW
python coil_run_iso.py A     $B $P $PERT pert_A     $WW
python analyze_iso.py
echo "COMPARE_ISO_DONE"
