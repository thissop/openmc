#!/usr/bin/env bash
# Local batch runner: iterate candidate device IDs, run each at production fidelity,
# log per-device status + wall seconds. Resumable (sweep.py skips existing records).
set -u
source /Users/tkiker/Documents/GitHub/openmc/spf_prototype/local_worker/activate_worker.sh >/dev/null 2>&1
cd /Users/tkiker/Documents/GitHub/openmc/spf_prototype/sweep

# Leaky-recovery: raise OpenMC's default max_lost_particles=10 to a MODEST absolute 100
# (=2.5e-5 of 4M histories, negligible vs the ~0.2% coil-flux stat error), with the
# relative cap disabled so the absolute count governs, and per-particle dump files
# capped at 5 so a grossly-leaky build can't flood the disk. Devices with a trivial
# leak recover to 'done'; devices with a real (>100) leak still abort -> stay 'error'.
export SPF_MAX_LOST=100
export SPF_REL_MAX_LOST=1.0
export SPF_MAX_WRITE_LOST=5

MANIFEST=configs/production_scaleup.json
OUT=sweep_out_local
LIST=sweep_out_local_candidates.txt
LOG=sweep_out_local/batch_progress.log
mkdir -p "$OUT"
echo "=== batch start $(date -u +%FT%TZ) ===" >> "$LOG"

while read -r ID; do
  [ -z "$ID" ] && continue
  REC="$OUT/config_${ID}_baseline.json"
  if [ -f "$REC" ]; then
    st=$(python -c "import json;print(json.load(open('$REC')).get('status'))" 2>/dev/null)
    echo "$(date -u +%FT%TZ) $ID SKIP_EXISTS $st" >> "$LOG"
    continue
  fi
  t0=$(date +%s)
  python sweep.py --manifest "$MANIFEST" --only "$ID" --out "$OUT" >/dev/null 2>&1
  t1=$(date +%s)
  st=$(python -c "import json;print(json.load(open('$REC')).get('status'))" 2>/dev/null || echo NO_RECORD)
  echo "$(date -u +%FT%TZ) $ID $st $((t1-t0))s" >> "$LOG"
done < "$LIST"
echo "=== batch done $(date -u +%FT%TZ) ===" >> "$LOG"
