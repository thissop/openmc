#!/bin/bash
# Queue the SSR pipeline: a fast SMOKE chain (tiny stats) gates the REAL overnight chain via afterok,
# so the expensive runs only fire if the smoke pipeline works end to end.
#   smoke_master -> smoke_gate2 -> real_master -> real_gate2
# Run on a Ginsburg login node from $WORK/ssr.
set -euo pipefail
export PATH=/cm/shared/apps/slurm/current/bin:$PATH
WORK=/burg-archive/home/tjk2147/pstl_test/corrected
UNI=dagmc_qh_step1_uniform.h5m
CELLS=cells_step1_uniform.npz
CF=2                                   # cellfrom = plasma void (freezes the inward first-wall source)
cd "$WORK/ssr"

sub() { sbatch --parsable "$@"; }

# --- smoke chain (minutes) ---
S1=$(sub --export=ALL,TAG=smoke_master,DAGMC=$UNI,BATCHES=2,PARTICLES=20000,MAXPART=500000,CELLFROM=$CF,CELLS=$CELLS slurm_ssw_master.sh)
S2=$(sub --dependency=afterok:$S1 --export=ALL,TAG=smoke_gate2,DAGMC=$UNI,BATCHES=2,PARTICLES=20000,SURF_SRC=$WORK/surf_smoke_master.h5,CELLS=$CELLS slurm_ssr.sh)

# --- real chain (overnight), gated on smoke success ---
R1=$(sub --dependency=afterok:$S2 --export=ALL,TAG=real_master,DAGMC=$UNI,BATCHES=12,PARTICLES=4000000,MAXPART=80000000,CELLFROM=$CF,CELLS=$CELLS slurm_ssw_master.sh)
R2=$(sub --dependency=afterok:$R1 --export=ALL,TAG=real_gate2,DAGMC=$UNI,BATCHES=12,PARTICLES=4000000,SURF_SRC=$WORK/surf_real_master.h5,CELLS=$CELLS slurm_ssr.sh)

echo "smoke_master=$S1  smoke_gate2=$S2  real_master=$R1  real_gate2=$R2"
squeue -u tjk2147 -o "%.12i %.12j %.8T %.10M %R" | grep -E "sswmaster|ssrdose|JOBID" || true
