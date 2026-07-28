#!/bin/bash
# One-shot orchestrator: wait for host load to settle, relaunch pre-audit ONCE at 3
# procs (leave it alone), wait for all 220 pre-audit records, then build the passer
# manifest. The separate worker-launcher (waits for configs/grow_passers.json) then
# runs the sweep. Emits milestone lines for monitoring.
set -u
cd /Users/tkiker/Documents/GitHub/openmc/spf_prototype/sweep
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python3

# 1) wait for load to settle
while [ "$($PY -c 'import os;print(1 if os.getloadavg()[0]<9 else 0)')" != "1" ]; do sleep 10; done
echo "[orch] load settled -> $(uptime | awk -F'load averages:' '{print $2}')"

# 2) relaunch pre-audit ONCE at 3 procs if not already running and not complete
have=$(ls sweep_out_grow/preaudit/ 2>/dev/null | wc -l | tr -d ' ')
running=$(pgrep -f preaudit_grow.py | wc -l | tr -d ' ')
if [ "$have" -lt 220 ] && [ "$running" = "0" ]; then
    nohup $PY preaudit_grow.py 3 >> sweep_out_grow/preaudit.log 2>&1 &
    echo "[orch] pre-audit relaunched at 3 procs (resume from $have)"
fi

# 3) wait for all 220
while [ "$(ls sweep_out_grow/preaudit/ 2>/dev/null | wc -l | tr -d ' ')" -lt 220 ]; do sleep 20; done
echo "[orch] PRE-AUDIT COMPLETE 220/220"

# 4) build manifest
$PY build_grow_manifest.py
echo "[orch] MANIFEST BUILT"
