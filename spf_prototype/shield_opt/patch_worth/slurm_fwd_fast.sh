#!/bin/bash
# Multigroup step (2): FAST-group forward mesh flux, co-located with the coil-20 adjoint.
# The coil dose response is FAST flux; the single-group signed worth used the TOTAL forward flux,
# which over-weights the (thermalized) breeder band. Recomputing the contributon C = phi_fast*psi
# with the FAST forward flux is the first refinement toward a physical multigroup signed worth
# (see SIGNED_WORTH_ANALYSIS.md). fwd_meshflux.py argv[8]="1" tallies only 0.1-25 MeV.
#
#   sbatch slurm_fwd_fast.sh
#SBATCH --job-name=fwdfast
#SBATCH --partition=burst
#SBATCH --account=astro
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/burg-archive/home/tjk2147/pstl_test/corrected/patch_worth/slurm_fwdfast_%j.log

set -uo pipefail
WORK=/burg-archive/home/tjk2147/pstl_test/corrected
PW=$WORK/patch_worth
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src
ADJ=/burg-archive/home/tjk2147/adjoint_test/recip_w35f_coil20_flat/adjoint_importance_flat_P0.npz
BASE_DAGMC=$WORK/dagmc_qh_step1_uniform.h5m
LI6=60; FWD_BATCHES=12; PARTICLES=4000000

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" PYTHONUNBUFFERED=1

cd "$WORK"
echo "=== FAST-group forward mesh flux (NO WW, 48M) $(date) ==="
# args: BATCHES PARTICLES DAGMC ADJOINT OUT_NPZ LI6 WW FAST
python -u "$PW/fwd_meshflux.py" "$FWD_BATCHES" "$PARTICLES" "$BASE_DAGMC" "$ADJ" \
    "$WORK/qh_fwd_meshflux_fast.npz" "$LI6" "" "1"
echo "=== DONE fwd fast rc=$? -> $WORK/qh_fwd_meshflux_fast.npz $(date) ==="
