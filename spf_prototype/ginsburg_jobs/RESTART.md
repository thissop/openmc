# Ginsburg restart — resume the overnight run

The SSH ControlMaster socket expired mid-session (Duo re-auth needed), so the **last two cluster
jobs did not run**: (1) register+build the C++ `StellaratorSource`, (2) the real-`√g` DESC solve on
`precise_QA` (with the fixed V1 quadrature). Both scripts are here; the code up to the V1 fix and
the register step did **not** get synced (that sync failed on the dead socket), so re-sync first.

Everything else is done and committed (`git log`, `MORNING_SUMMARY.md`).

## Step 0 — re-auth (you do the Duo)
```
ssh ginsburg          # complete the Duo push, then you can leave it / close it
```

## Step 1 — re-sync the latest code + batch scripts (from the Mac)
```
LOCAL=/Users/tkiker/Documents/GitHub/openmc/spf_prototype
SPF=/ginsburg/astro/users/tjk2147/spf_work ; PP=$SPF/spf_pp
( cd "$LOCAL" && COPYFILE_DISABLE=1 tar czf - python native_spf ) | ssh ginsburg "cd $PP && tar xzf -"
( cd "$LOCAL/ginsburg_jobs" && tar czf - *.sbatch )                | ssh ginsburg "cd $SPF && tar xzf -"
```

## Step 2 — submit the two jobs
```
ssh ginsburg "cd $SPF && sbatch register_stellarator.sbatch && sbatch fluxmap.sbatch"
```
Expected:
- **register_stellarator**: prints `CMake/source.cpp/source.py patched`, rebuilds, `nm ... | grep -c
  StellaratorSource` > 0, and the numpy-mirror `test_stellarator_source_cpp.py` passes on the cluster.
  ⇒ the StellaratorSource core class is compiled into OpenMC.
- **fluxmap** (`precise_QA`): prints `V0 … V1 … V1b … V2 … V3 … V4 skipped … ALL RIGOR GATES PASSED`
  and writes `spf_pp/data/precise_QA_fluxmap.{npz,meta,bin}`. ⇒ real `√g` validated.

## Step 3 — validate the sampler on the real equilibrium (Mac)
```
scp ginsburg:$PP/data/precise_QA_fluxmap.npz "$LOCAL/data/"
python "$LOCAL/python/validate_real_fluxmap.py" "$LOCAL/data/precise_QA_fluxmap.npz"
```
Expect: unit weight, inside-plasma, |b̂|=1, rho-marginal χ²/dof < 3 → sampler OK on real `√g`.

## Step 4 — the capstone: transport ±43% for the native sources (needs a model)
Reproduce Schwartz's +43% inboard / ∓22% outboard midplane NWL with the native SPF `TokamakSource`
(b̂=toroidal) and `StellaratorSource` (axisymmetric/circular fluxmap). Reuse the tier-2 machinery +
anarrima (`/Users/tkiker/Documents/GitHub/anarrima`). This is the last integration step; write it
carefully (the inboard/outboard reduction is the fiddly part that bit the coil-flux η earlier).

## Reference — the run environment on Ginsburg
```
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh ; conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SRC=$SPF/openmc_src
export PYTHONPATH="$SRC:$PYTHONPATH"                 # PR-#3999 python openmc (has native TokamakSource+SPF)
export LD_PRELOAD="$SRC/build/lib/libopenmc.so"      # conda ships a competing libopenmc first in RPATH
export PATH="$SRC/build/bin:$PATH"                    # the from-source openmc exe
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml
# DESC (separate venv, deps conflict with openmc):  $HOME/desc_venv/bin/python
```

## Known issue to fix (not blocking the above)
QUASR-boundary → DESC (`quasr_equilibrium_field._desc_surface`, `[C1]`): device-specific solves die in
`ensure_positive_jacobian` (degenerate axisymmetric seed) — a boundary Fourier mode-sign / theta-
orientation convention. `precise_QA` (a DESC example) is used meanwhile to validate the `√g` machinery.
