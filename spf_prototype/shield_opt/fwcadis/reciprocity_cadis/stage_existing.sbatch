#!/usr/bin/env bash
#SBATCH --job-name=fwc-stage
#SBATCH --partition=short
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=stage_existing_%j.out
#SBATCH --error=stage_existing_%j.err

set -euo pipefail

source_root=/burg-archive/home/tjk2147/pstl_test/corrected
stage_root="$SLURM_SUBMIT_DIR/openmc_fwcadis/reciprocity_cadis/source_snapshot"
mkdir -p "$stage_root/ssr"

for name in \
    op_fwcadis.py \
    op_ssr_wwapply.py \
    op_ssr.py \
    op_master.py \
    slurm_op_fwcadis2.sh \
    slurm_op_ssr.sh; do
    if [[ -f "$source_root/ssr/$name" ]]; then
        cp "$source_root/ssr/$name" "$stage_root/ssr/$name"
    fi
done

for name in cellmap3.py findmagnets.py build_oneperiod.py; do
    if [[ -f "$source_root/$name" ]]; then
        cp "$source_root/$name" "$stage_root/$name"
    fi
done

{
    printf 'SOURCE_ROOT=%s\n' "$source_root"
    for path in \
        "$source_root/dagmc_qh_op_uniform.h5m" \
        "$source_root/coil_centroids_corr.npz" \
        "$source_root/ssr/surf_op_master.h5" \
        "$source_root/ssr/coil_op_master.npz" \
        "$source_root/ssr/mgxs.h5" \
        "$source_root/ssr/ww_op_fwcadis.h5"; do
        if [[ -f "$path" ]]; then
            stat -c '%n|%s|%y' "$path"
        else
            printf '%s|MISSING\n' "$path"
        fi
    done
} > "$stage_root/manifest.txt"

printf 'STAGED_FILES\n'
find "$stage_root" -maxdepth 2 -type f -printf '%P|%s\n' | sort
printf 'DONE_STAGE_EXISTING\n'
