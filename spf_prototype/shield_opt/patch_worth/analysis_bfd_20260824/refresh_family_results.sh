#!/usr/bin/env bash
# Pull any completed equal-mass family results from Ginsburg and rerun the
# local, standard-library uncertainty gates. This script does not submit,
# cancel, or modify cluster work.
set -euo pipefail

analysis_dir=$(cd "$(dirname "$0")" && pwd)
results_dir="$analysis_dir/family_transport_results"
remote_host=ginsburg
remote_root=/burg-archive/home/tjk2147/pstl_test/corrected
common_mass_kg=2163095.2194506074

mkdir -p "$results_dir"

for family_tag in famU famC famR famM; do
  remote_file="$remote_root/opt/result_${family_tag}.json"
  if ssh -o ConnectTimeout=20 "$remote_host" "test -f '$remote_file'"; then
    scp -q -o ConnectTimeout=20 \
      "$remote_host:$remote_file" \
      "$results_dir/result_${family_tag}.json"
    printf 'pulled result_%s.json\n' "$family_tag"
  fi
done

python3 "$analysis_dir/family_transport_analysis.py" \
  --results-dir "$results_dir" \
  --output-json "$results_dir/family_transport_summary.json" \
  --output-markdown "$analysis_dir/FAMILY_TRANSPORT_RESULTS.md" \
  --common-mass-kg "$common_mass_kg"

ssh -o ConnectTimeout=20 "$remote_host" \
  'export PATH=/cm/shared/apps/slurm/current/bin:$PATH; squeue -j 9569111 -o "%.10i %.9T %.10M %.10L %R" 2>/dev/null || true'
