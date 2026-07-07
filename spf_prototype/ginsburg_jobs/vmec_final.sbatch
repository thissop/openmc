#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_vmecfin
#SBATCH -p short
#SBATCH -t 0-00:30
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
module load anaconda/3-2023.09 2>/dev/null
export OMP_NUM_THREADS=8
cd /ginsburg/astro/users/tjk2147/spf_work/spf_pp/python
"$HOME/vmec_venv/bin/python" - <<'PY'
import sys, os, glob, numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "../sweep")
import vmecpp, simsopt
sbase = os.path.dirname(simsopt.__file__)
inputs = sorted(glob.glob(sbase + "/**/input.*", recursive=True))
print("simsopt inputs:", [os.path.basename(i) for i in inputs][:15], flush=True)
tmpl = None
for c in inputs:
    try:
        t = vmecpp.VmecInput.from_file(c)
        if getattr(t, "rbc", None) is not None:
            tmpl = (c, t); break
    except Exception as e:
        print("  skip", os.path.basename(c), type(e).__name__, flush=True)
assert tmpl, "no simsopt VMEC template loaded"
path, inp = tmpl
rbc0 = np.asarray(inp.rbc)
print(f"TEMPLATE {os.path.basename(path)}: nfp={inp.nfp} mpol={inp.mpol} ntor={inp.ntor} rbc.shape={rbc0.shape}", flush=True)

import quasr_equilibrium_field as qef
nfp, Rm, Zm, meta = qef._quasr_boundary_modes(803097, flip_theta=True)
# match the template's array shape convention. rbc0 shape is (mpol_t+1, 2*ntor_t+1).
mpol_t = rbc0.shape[0] - 1; ntor_t = (rbc0.shape[1] - 1) // 2
mpol = min(6, mpol_t); ntor = min(6, ntor_t)
rbc = np.zeros_like(rbc0); zbs = np.zeros_like(rbc0)
for (m, n), (rc, rs) in Rm.items():
    if m <= mpol and abs(n) <= ntor: rbc[m, n + ntor_t] = rc
for (m, n), (zc, zs) in Zm.items():
    if m <= mpol and abs(n) <= ntor: zbs[m, n + ntor_t] = zs
R00 = Rm[(0, 0)][0]; a = abs(Rm.get((1, 0), (0.1, 0))[0]) or 0.1
rax = np.zeros_like(np.asarray(inp.raxis_c)); rax[0] = R00
upd = dict(nfp=nfp, lasym=False, rbc=rbc, zbs=zbs, raxis_c=rax,
           zaxis_s=np.zeros_like(np.asarray(inp.zaxis_s)),
           ns_array=np.array([16, 51, 99]), ftol_array=np.array([1e-8, 1e-10, 1e-12]),
           niter_array=np.array([2000, 3000, 6000]), phiedge=float(np.pi * a * a),
           ncurr=1, curtor=0.0, pres_scale=0.0)
try:
    inp2 = inp.model_copy(update=upd)
except Exception:
    inp2 = inp
    for k, v in upd.items():
        try: setattr(inp2, k, v)
        except Exception as e: print("  set fail", k, type(e).__name__, e)
try:
    out = vmecpp.run(inp2, verbose=False)
    g = np.asarray(out.wout.gmnc)
    print(f"  *** VMECPP CONVERGED (803097) *** volume={out.wout.volume_p:.4f} "
          f"aspect={out.wout.aspect:.2f} gmnc{g.shape} sign={np.sign(np.median(g)):+.0f}", flush=True)
except Exception as e:
    import traceback; traceback.print_exc()
    print(f"  vmecpp FAILED: {type(e).__name__}: {e}", flush=True)
print("DONE_VMEC_FINAL")
PY
