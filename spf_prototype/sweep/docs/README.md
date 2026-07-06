# SPF quasisymmetric-stellarator sweep: overview + runbook

Test whether the cheap, field-only coherence C predicts the expensive neutronics
directional efficiency eta of spin-polarized fusion (SPF) neutron steering, across a
family of QUASR quasisymmetric stellarators, and whether eta factorizes into a field
part and a blanket part.

## Hypotheses (falsifiable outputs, not assumptions)
1. eta_source collapses onto a single curve eta_source(C) across symmetry classes
   (QA/QH/QI). See THEORY.md.
2. eta factorizes: eta(C, blanket) ~ eta_source(C) * A(tau). See FACTORIZATION.md.

## eta definition
Single source of truth in EXPERIMENTAL_DESIGN.md and analyze_sweep.py's header.
eta_source = delta_free/delta_free(anchor), from the NEAR-SOURCE first-wall directional
contrast (free stream; low variance, G4); eta_coil = delta_scatter/delta_scatter(anchor),
from the deep coil fast flux (scatter). A = eta_coil/eta_source. Anchor = highest-C config
(eta_*(anchor)=1). The PARITY-CORRECT predictor is the nematic order S_phi=(3*lambda_phi-1)/2,
NOT the first-moment C (C is a proxy, valid on cap-shaped fields where lambda_phi~C^2);
analyze_sweep fits eta against BOTH and reports the winner. See THEORY.md (parity argument).

## Modules (all built; unit-tested locally in seconds)
- coherence_metrics.py : C, direction tensor, angular std.  DONE.
- biotsavart_field.py  : exact b_hat from coil filaments + spf_fieldmap_v1 writer.  DONE.
- quasr_loader.py      : QUASR coils (simsopt-serial JSON, no simsopt) + boundary + meta.  DONE.
- field_audit.py       : HARD-GATE unit / div-free / B.n-on-LCFS. DONE.
- sweep.py             : per-config driver (field->audit->C->geometry->DAGMC->transport->record). DONE.
- select_configs.py    : compute C + S_phi for candidates, tile the axis, class x C
                         coverage + aliasing/reversal report, --field-source. DONE.
- analyze_sweep.py     : aggregate, fit eta_source(C) AND eta_source(S_phi), C-vs-S_phi
                         comparison, universality, A(tau), figures. DONE.
- sweep.sbatch         : SLURM array on burst. DONE.
Reused (debugged, from the Ginsburg copy): dagmc_writer.build_from_stls (pymoab-free
DAGMC), run_conformal.build_model, run_ginsburg._mat_score (array tally reader),
export_for_plots.py.

## Local smoke test (no cluster; runs everything except transport)
```
cd spf_prototype/sweep
python3 -m pytest tests -q                                  # 36 tests, seconds
python3 field_audit.py 59509                                # coil-field hard-gate audit
python3 sweep.py --manifest configs/smoke.json --out /tmp/smoke --force   # field->...->DAGMC->record
python3 analyze_sweep.py --records <records-with-transport> --out /tmp/an  # (after Ginsburg)
```
Locally, transport auto-skips (no openmc), so records land as `staged_no_transport`
with C, the audit, tensor metrics, tau, scale, and a valid DAGMC .h5m.

---

# GINSBURG RUNBOOK (exact commands)

Ginsburg facts (from ../../GINSBURG_CLAUDE_HANDOFF.md, ../../SETUP_GINSBURG.md): account
`astro`; env `spf-stellarator`; no shared cross sections (download once); LOGIN NODE
= network only, ALL compute on a compute node / sbatch; the repo root has an
`openmc/` SOURCE tree that must NOT shadow the conda openmc.

Set these once in your shell:
```
MAC_REPO=/Users/tkiker/Documents/GitHub/openmc          # this clone (has the sweep + quasr_*)
GINS=tjk2147@ginsburg.rcs.columbia.edu                  # <-- your login
GREPO='$HOME/src/GitHub/openmc'                          # <-- clone path on Ginsburg (single-quoted)
```

## (a) PUSH to Ginsburg
rsync the sweep engine + the QUASR data + the coil serials it already downloaded.
```
# code: the new sweep engine and the incorporated debugged pipeline files
rsync -avz --exclude '__pycache__' \
  "$MAC_REPO/spf_prototype/sweep/" "$GINS:$GREPO/spf_prototype/sweep/"
rsync -avz \
  "$MAC_REPO/spf_prototype/python/dagmc_writer.py" \
  "$MAC_REPO/spf_prototype/python/export_for_plots.py" \
  "$MAC_REPO/spf_prototype/python/run_ginsburg.py" \
  "$GINS:$GREPO/spf_prototype/python/"
rsync -avz "$MAC_REPO/spf_prototype/src/polarized_fusion_source.cpp" \
  "$GINS:$GREPO/spf_prototype/src/"
# QUASR inputs + the coil serials select_configs already fetched (so the job is offline)
rsync -avz "$MAC_REPO/spf_prototype/data/quasr/" "$GINS:$GREPO/spf_prototype/data/quasr/"
```

## (b) SET UP + CLEAN the environment  (LOGIN node for the env; COMPUTE node for builds/tests)
Login node (network only; the env should already exist per SETUP_GINSBURG.md):
```
ssh $GINS
unset PYTHONPATH
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator || conda env create -f "$GREPO/spf_prototype/environment.yml"
export OPENMC_CROSS_SECTIONS=$HOME/endfb-viii.0-hdf5/cross_sections.xml   # download once if absent
```
Grab a COMPUTE node for everything that computes (never the login node):
```
salloc -A astro -N 1 -c 8 -t 2:00:00
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh && conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
cd "$GREPO"
# CLEAN + VERIFY (each is a hard gate; stop at the first failure):
python -c "import openmc,sys; assert 'site-packages' in openmc.__file__, openmc.__file__; \
           assert hasattr(openmc,'DAGMCUniverse'); print('conda openmc OK:', openmc.__file__)"
python -c "import os; p=os.environ['OPENMC_CROSS_SECTIONS']; assert os.path.exists(p); print('XS OK:', p)"
python spf_prototype/python/dagmc_writer.py --selftest        # pymoab-free writer: structure OK
python spf_prototype/python/dagmc_writer.py --selftest-stl    # end-to-end STL->h5m OK
python -m pytest spf_prototype/sweep/tests -q                 # 24 sweep tests
# PREBUILD the compiled source .so ONCE (array tasks must not race-build it):
cd spf_prototype/src && cmake -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" . && cmake --build build && cd "$GREPO"
# optional: rebuild the config family on the cluster (else use the pushed production.json)
# python spf_prototype/sweep/select_configs.py --auto 30 --tiles 12 \
#     --out spf_prototype/sweep/configs/production.json
```

## (b.5) BRING-UP: run ONE config end-to-end FIRST (do not skip)
Transport for a QUASR config has never actually run (only precise_QA ran before;
the sweep's local smoke validated everything EXCEPT OpenMC transport). Shake out the
transport integration on ONE cheap config on the compute node before the array:
```
cd "$GREPO"
OMP_NUM_THREADS=8 python spf_prototype/sweep/sweep.py \
    --manifest spf_prototype/sweep/configs/smoke.json --out /tmp/bringup --only 59509
python -c "import json; r=json.load(open('/tmp/bringup/config_59509_baseline.json')); \
           print(r['status']); print('delta_free_parallel', r.get('delta_free_parallel'))"
```
Expect status `done` and a FINITE delta_free_parallel. If not, check (in order): the
DAGMC .h5m loads (`python spf_prototype/python/dagmc_writer.py --load-test \
spf_prototype/data/quasr59509_baseline_geom/quasr59509.h5m` -> material_names
W/steel/Be/FLiBe/shield/coil); the .so reads the biotsavart field map (no NaN
births); the source R0/aminor sit inside the shaped wall. Then run a config with
`"streams": "both"` to exercise the scatter path before trusting eta_coil.

## (c) LAUNCH the sweep array (burst partition)
The manifest `configs/production.json` has N configs (the C-scan + the blanket
subset). Set the array range to 0..N-1. Current production.json has 12 configs.
```
cd "$SUBMIT_DIR"                 # a scratch/results dir you own (NOT the repo root)
export SPF_REPO="$GREPO"
export OPENMC_CROSS_SECTIONS=$HOME/endfb-viii.0-hdf5/cross_sections.xml
export MANIFEST="$GREPO/spf_prototype/sweep/configs/production.json"
N=$(python -c "import json;print(len(json.load(open('$MANIFEST'))['configs']))")
sbatch --array=0-$((N-1))%8 "$GREPO/spf_prototype/sweep/sweep.sbatch"
```
Per-config statistics live in the manifest defaults (particles=200000, batches=20 for
the C-scan; streams='free' so only eta_source is computed cheaply). The blanket
subset entries carry streams='both' for eta_coil + A. eta_coil is VR-pending: set
`"vr": true` on those entries (and expect longer runs) for a converged coil eta.

## (d) MONITOR
```
squeue -u $USER                                  # array task states
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS   # per-task exit/status
ls $SUBMIT_DIR/sweep_out/config_*.json | wc -l   # per-config records appearing (resumable)
# a record's status field: done | staged_no_transport | audit_failed | geometry_folded | error
grep -l '"status": "error"' $SUBMIT_DIR/sweep_out/config_*.json   # any failures
grep -H lost_particles $SUBMIT_DIR/sweep_out/config_*.json        # lost-particle guard (auto-retries scale)
```
Each array task writes ONE `config_<id>_<blanket>.json`; a crash loses at most that
config, and re-submitting the array skips existing records (idempotent; add --force
to recompute).

## (e) PULL results back to the Mac
Only the small JSON/CSV are needed for analysis (statepoints stay on Ginsburg).
```
rsync -avz "$GINS:$SUBMIT_DIR/sweep_out/" "$MAC_REPO/spf_prototype/sweep/sweep_out/"
```

## (f) ANALYZE locally -> scaling-law figures + verdicts
```
cd "$MAC_REPO/spf_prototype/sweep"
python3 analyze_sweep.py --records sweep_out --out analysis
# -> analysis/sweep.csv, analysis/report.json, analysis/figs/*.png
```
report.json states: the eta_source(C) fit (slope/R^2, does cheap C predict expensive
eta), the universality verdict (do QA/QH/QI collapse, or use the tensor), the A(tau)
factorization, and which separability failure mode (if any) is active. Figures:
fig_C_predicts_eta.png (the headline), fig_eta_vs_nfp.png, fig_A_of_tau.png,
fig_direction_spread.png.

## Honest status labels
- eta_source: PRIMARY, fast-converging (free stream), reported as the scaling law.
- eta_coil / A: VR-PENDING. Without weight windows the deep coil flux is noisy; lead
  with eta_source and the trend, quote coil eta with its error. Enable `--vr` per
  config for a converged coil eta.
- The whole run is machinery + method validation on public equilibria, not a reactor
  claim. See SCOPE_AND_CAVEATS.md.
