#!/bin/bash
# Run the pre-audit PASSER manifest through sweep.py on the Lima 'spf' worker,
# resumable, into sweep_out_grow/. One device at a time (sweep.py loops the manifest;
# each config JSON is written atomically so a crash loses at most one). Emits a
# per-device wall-clock line so we can measure real done/hr.
#
# Launch from the Mac:  nohup bash run_grow_worker.sh > sweep_out_grow/worker.log 2>&1 &
set -u
SWEEP=/Users/tkiker/Documents/GitHub/openmc/spf_prototype/sweep
MAN=configs/grow_passers.json
OUT=sweep_out_grow

limactl shell spf bash -lc "
source /Users/tkiker/Documents/GitHub/openmc/spf_prototype/local_worker/activate_worker.sh >/dev/null 2>&1
cd $SWEEP
export OMP_NUM_THREADS=4
python - <<'PY'
import json, subprocess, time
from pathlib import Path
man = json.load(open('$MAN'))
ids = [c['id'] for c in man['configs']]
outdir = Path('$OUT')
print(f'[worker] {len(ids)} passers to run (resumable)', flush=True)
t_start = time.time(); ndone_new = 0
for i, cid in enumerate(ids):
    rec = outdir / f'config_{cid}_baseline.json'
    if rec.exists():
        continue
    t0 = time.time()
    subprocess.run(['python','sweep.py','--manifest','$MAN','--only',str(cid),
                    '--out','$OUT'], check=False)
    dt = time.time() - t0
    st = 'unknown'
    if rec.exists():
        try: st = json.load(open(rec)).get('status','?')
        except Exception: st = 'parse_err'
    if st == 'done': ndone_new += 1
    el = time.time() - t_start
    rate = ndone_new/(el/3600.0) if el>0 else 0
    print(f'[worker] {i+1}/{len(ids)} id={cid} status={st} wall={dt:.0f}s '
          f'| new_done={ndone_new} net_done/hr={rate:.1f}', flush=True)
print(f'[worker] FINISHED: {ndone_new} new done in {(time.time()-t_start)/3600:.2f} h', flush=True)
PY
"
