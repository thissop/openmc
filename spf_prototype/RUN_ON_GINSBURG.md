# RUN_ON_GINSBURG.md — conformal stellarator SPF neutronics on x86

This package was developed on a throwaway aarch64 sandbox, which **cannot** run
the compiled DAGMC/OpenMC-with-DAGMC stack (conda-forge has no `openmc`/`dagmc`
for linux-aarch64). Everything architecture-independent (the spin-polarized
source, the field map, the conformal geometry, the analytic cross-check) is built
and validated; the **conformal transport runs on Ginsburg (x86)** via the conda
env below. Nothing in the code needs editing — only the env and cross-section path
resolve at run time.

> **Pipeline at a glance** (each step is one script):
> `desc_to_fieldmap.py` (field map + LCFS surface, *desc env*) → `stellarator_geometry.py` (conformal STLs, *any env*) → `build_dagmc.py` (`.h5m`, *conda*) → `run_conformal.py` (transport + metrics, *conda*).

## 0. One-time setup

```bash
# (a) transport env — OpenMC+DAGMC + analytic cross-check stack
conda env create -f spf_prototype/environment.yml
conda activate spf-stellarator

# (b) DESC env — SEPARATE (its jax/numpy pins conflict with OpenMC); producer only
python -m venv ~/desc_venv
~/desc_venv/bin/pip install desc-opt

# (c) cross sections (ENDF/B HDF5; any reachable library)
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml   # the code honors this

# (d) build the spin-polarized compiled source against the CONDA OpenMC
cd spf_prototype/src
cmake -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" .
cmake --build build          # -> build/libpolarized_fusion_source.so
cd ../..
```

## 1. Smoke-test the env BEFORE any long job

```bash
conda activate spf-stellarator

# (i) architecture-independent gates: the SPF sampler + angled-kernel analytic
#     validation (Part A rigor anchor) + the conformal geometry math.
python -m pytest spf_prototype/tests -q
#     EXPECT: all pass (sampler C++/Py parity, anarrima angled-kernel <1e-6,
#     field-map parity, geometry watertight/simple/nested).

# (ii) DAGMC round-trip on a tiny build (confirms moab/dagmc/openmc-with-dagmc):
python spf_prototype/python/desc_to_fieldmap.py precise_QA equil_precise_qa   # in desc env if needed
python spf_prototype/python/stellarator_geometry.py equil_precise_qa 10.0
python spf_prototype/python/build_dagmc.py spf_prototype/data/equil_precise_qa_geom stellarator.h5m
#     EXPECT: "wrote stellarator.h5m with material volumes: [W, steel, Be, FLiBe, shield, coil]"
#     If add_stl_file errors, the cad_to_dagmc API moved — fix the single call site
#     flagged in build_dagmc.py against the installed version.

# (iii) low-history FREE-STREAMING conformal run -> must reproduce the analytic
#       free-streaming NWL of the companion paper (the JOINT analytic<->MC check).
OMP_NUM_THREADS=4 python spf_prototype/python/run_conformal.py stellarator.h5m equil_precise_qa 10.0
#     (lower `particles` in run_conformal for the smoke pass.)
```

The free-streaming smoke run is the **joint-validation anchor**: in the near-void
limit the OpenMC per-patch, per-mode wall load must match the analytic
free-streaming NWL (anarrima + the field-period-perturbation series of
`docs/Analytic_Neutron_Wall_Loading...`). Agreement there validates the source
sampling, the polarized angular distribution, and the line-of-sight geometry
*before* trusting the collided (scattering) result.

## 2. The science runs

`run_conformal.py` runs **unpolarized / perpendicular(A) / parallel(B/C)**, each in
**free-streaming** and **scattering**, at fixed neutron rate, per source neutron
(Bae 2025 convention), and reports: φ-resolved first-wall load, per-material
heating/damage/TBR, coil fast flux (>0.1 MeV), and the directional efficiency η
(steering retained). The headline question — **does SPF steering survive
non-axisymmetric conformal smearing, and how much does scattering dilute it** — is
answered by comparing free-streaming vs scattering steering on the conformal wall.

### The QA→QH quasisymmetry scan (the main result)
Loop the whole pipeline over a controlled family of equilibria and correlate η
against a quasisymmetry / non-axisymmetry metric:

```bash
for EQ in precise_QA precise_QH ; do          # add QUASR devices for a denser scan
  ~/desc_venv/bin/python spf_prototype/python/desc_to_fieldmap.py $EQ equil_$EQ
  python spf_prototype/python/stellarator_geometry.py equil_$EQ 10.0
  python spf_prototype/python/build_dagmc.py spf_prototype/data/equil_${EQ}_geom ${EQ}.h5m
  OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK python spf_prototype/python/run_conformal.py ${EQ}.h5m equil_$EQ 10.0
done
```

Note QH (precise_QH) has larger field-period excursion → larger ε_eff and a
non-negligible field-direction ripple; the conformal offset may need a larger
`scale` (thinner build relative to minor radius) to stay self-intersection-free
(`stellarator_geometry.py` flags this — it refuses invalid builds). Variance
reduction (weight windows / MAGIC) is recommended for the deep-coil tally; without
it the coil fast-flux carries ~25–30% MC error (lead with the resolved first-wall η).

## 3. Device scale / Helios injection

`scale` maps the normalized precise_QA (R0≈1 m, a≈0.17 m) to a reactor (×10 ≈
Helios a≈1.7 m). The field map's B̂ direction is scale-invariant (`run_conformal`
rescales the grid bounds automatically). To inject real Helios numbers, replace
the equilibrium (`desc_to_fieldmap.py <name>` or a real `wout` via `desc.VMECIO`)
and the radial build (`stellarator_geometry.DEFAULT_LAYERS`) — every device-specific
slot is tagged `INJECT(helios)`. See `DATA_NEEDED.md`.

## Gotchas
- OpenMC+DAGMC needs a **graveyard / vacuum boundary**; `run_conformal.make_geometry`
  uses `DAGMCUniverse(...).bounded_universe()` to add one automatically.
- The `.so` must be built against the **same** OpenMC it runs with (the conda one
  here), not the sandbox build — step 0(d).
- Run the DESC producer in the **desc env**, everything else in **spf-stellarator**.
