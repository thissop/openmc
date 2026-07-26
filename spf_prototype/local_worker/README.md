# Local OpenMC + DAGMC worker (arm64 Lima VM `spf`)

A persistent Linux (Ubuntu 24.04, aarch64) OpenMC **0.15.3 with DAGMC** build,
living in a Lima VM on the M2 Mac, for running the **cheap SPF transport tier**
(conformal free-streaming sweep + pipeline dev/validation) with **zero cluster
queue latency**. The heavy magnet DAGMC+adjoint job stays on Ginsburg.

Built 2026-07-26. Path built once; re-usable across sessions.

## TL;DR — run a device
```bash
limactl start spf                      # if not already Running
limactl shell spf
source /Users/tkiker/Documents/GitHub/openmc/spf_prototype/local_worker/activate_worker.sh
cd /Users/tkiker/Documents/GitHub/openmc/spf_prototype/sweep
python sweep.py --manifest configs/smoke.json --only 59509 --out /tmp/out
# -> /tmp/out/config_59509_baseline.json  (C, field_audit, delta_free_*, dagmc, ...)
```
The Mac repo `/Users/tkiker/Documents/GitHub/openmc` is mounted **writable** at the
same path inside the VM, so you edit code on the Mac and run it in the VM.

## What is installed (all inside the VM, under `~`)
| Component | Version | Location |
|-----------|---------|----------|
| Miniforge / conda env `spf` | py3.11 | `~/miniforge3/envs/spf` |
| MOAB (conda-forge) | 5.6.0 | `$CONDA_PREFIX/lib` |
| DAGMC (source, no double_down/embree) | 3.2.4 | `~/opt/dagmc` |
| OpenMC C++ (source, `-DOPENMC_USE_DAGMC=ON`, vendored xtensor) | 0.15.3 | `~/opt/openmc` |
| OpenMC Python (`pip install .`) | 0.15.3 | env `spf` |
| SPF compiled source `.so` | rebuilt vs `~/opt/openmc` | `spf_prototype/src/build/libpolarized_fusion_source.so` |
| Cross sections (ENDF/B-VIII.0 subset, neutron) | 26 nuclides, 386 MB | `~/xs/endfb80/cross_sections.xml` |

`openmc --version` reports **`DAGMC support: yes`** — this is the make-or-break
that conda-forge could not provide on aarch64.

## Environment (what activate_worker.sh sets)
- `conda activate spf`
- `PATH`   += `~/opt/openmc/bin`
- `LD_LIBRARY_PATH` = `~/opt/openmc/lib:~/opt/dagmc/lib:$CONDA_PREFIX/lib`  (required — RPATHs are not self-contained)
- `OPENMC_CROSS_SECTIONS=~/xs/endfb80/cross_sections.xml`
- `OMP_NUM_THREADS=4`

## Cross-section subset
Nuclides the conformal blanket uses (W, steel = Fe/Cr/Ni, Be, FLiBe = Li/Be/F,
WC shield = W/C, vacuum = H1), fetched per-nuclide with `openmc_data_downloader`:
```bash
openmc_data_downloader -l ENDFB-8.0-NNDC -p neutron -e H Li Be C F Cr Fe Ni W -d ~/xs/endfb80
```
Neutron only; no S(alpha,beta) (source is monoenergetic 14.1 MeV, materials near-void
for free-streaming), no photon.

## Rebuilding the SPF `.so` (if OpenMC is rebuilt)
The pipeline's `build_and_run.build_so()` cmake will FAIL against a conda HDF5 because
OpenMC's exported cmake target lists the HDF5 imported targets (`hdf5-shared`,
`hdf5_hl-shared`) but `OpenMCConfig.cmake` does not re-`find_package(HDF5)`. Inject it:
```bash
echo "find_package(HDF5 CONFIG REQUIRED)" > /tmp/hdf5inject.cmake
cmake $REPO/spf_prototype/src -G Ninja \
  -DCMAKE_PREFIX_PATH="$HOME/opt/openmc;$HOME/opt/dagmc;$CONDA_PREFIX;$CONDA_PREFIX/cmake" \
  -DCMAKE_PROJECT_spf_source_INCLUDE=/tmp/hdf5inject.cmake -DCMAKE_BUILD_TYPE=Release -B /tmp/sob
ninja -C /tmp/sob
cp /tmp/sob/libpolarized_fusion_source.so $REPO/spf_prototype/src/build/
```

## VM lifecycle
- `limactl list` / `limactl start spf` / `limactl stop spf`
- VM specs: 4 CPU, 8 GB RAM, 28 GB disk (leaves native cores/RAM for simsopt).
- Persistent: the VM disk and all builds survive stop/start.

## Scope
This worker runs the **free-streaming / near-void conformal tier only**. It does NOT
hold the full-reactor magnet adjoint (16 GB is too small) — that stays on Ginsburg.
