#!/bin/bash
# SLURM ARRAY scaffold for the magnet-side neutronics suite (200-300 devices).
#
# Each array task launches the FULL dependency-chained per-device pipeline for ONE device
# via run_device.py --submit (fetch inline -> equil -> dagmc -> cells -> adjoint -> finalize).
# The array task itself is lightweight (it submits the chain and exits); the heavy work runs
# as the chained sub-jobs, right-sized per stage. Resumable: run_device skips completed
# stages, so re-running the array picks up where each device left off. One JSON per device.
#
# GATE (do NOT launch the full 200-300 until BOTH hold):
#   (a) the automation passed end-to-end on >=1 FRESH device (adjoint gate PASS, JSON written);
#   (b) the adaptive-shielding kill-shot (job 9193026 / step1) verdict validates the shield step.
#
# Usage:
#   1. build the device list (subset of engineering_relevant_devices.csv):
#        awk -F, 'NR>1 && $7==1 {print $1}' ../data/engineering_relevant_devices.csv > devices.txt
#      (or a small validation batch: head -5 devices.txt > devices_batch.txt)
#   2. sbatch --array=1-$(wc -l < devices_batch.txt)%5 slurm_array.sh devices_batch.txt
#      (%5 = at most 5 device-chains launching concurrently)
#
#SBATCH -A astro
#SBATCH -J magnet_suite
#SBATCH -p short
#SBATCH -t 0-00:20
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=4G
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
set -uo pipefail

DEVLIST="${1:?usage: sbatch --array=1-N%C slurm_array.sh devices.txt}"
SUITE="${SUITE:-$(cd "$(dirname "$0")" && pwd)}"
WORK="${MAGNET_SUITE_WORK:-$HOME/pstl_test/magnet_suite}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DEVLIST" | tr -d '[:space:]')
if [[ -z "$ID" ]]; then echo "no device on line $SLURM_ARRAY_TASK_ID"; exit 0; fi
echo "=== array task $SLURM_ARRAY_TASK_ID -> device $ID  (suite=$SUITE work=$WORK) ==="

source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator          # numpy + internet for the inline fetch stage

python "$SUITE/run_device.py" "$ID" --work "$WORK" --suite "$SUITE" --submit
echo "=== launched device $ID chain ==="
