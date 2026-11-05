#!/bin/bash
#SBATCH -J opfwcadis
#SBATCH -p short
#SBATCH --account=astro
#SBATCH -t 4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=96G
#SBATCH -o /burg-archive/home/tjk2147/pstl_test/corrected/ssr/slurm_%x_%j.log
set -uo pipefail
WORK=/burg-archive/home/tjk2147/pstl_test/corrected
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
export RUNDIR=$WORK/fwcadis_${TAG:-smoke}; mkdir -p "$RUNDIR"
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh; conda activate spf-stellarator
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC LD_PRELOAD=$SRC/build/lib/libopenmc.so PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=32 PYTHONUNBUFFERED=1 HDF5_USE_FILE_LOCKING=FALSE
echo "=== op_fwcadis rays=$PARTICLES ==="
python -u "$WORK/ssr/op_fwcadis.py" "$PARTICLES" 60 || echo "rc=$?"
echo "END"
