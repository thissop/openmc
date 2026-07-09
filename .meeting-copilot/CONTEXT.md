# OpenMC — Quick-Look Context for Meeting Copilot

## What It Does

**OpenMC** is a fully-featured Monte Carlo particle transport code for nuclear reactor analysis and fusion neutron transport. It performs continuous-energy simulations using HDF5-based cross-section libraries (ENDF/B), models arbitrary 3D geometries via Constructive Solid Geometry (CSG), tracks neutrons/photons/electrons through materials via ray-tracing, scores results using flexible tallying, and supports both fixed-source and eigenvalue (k-eff) modes with MPI and OpenMP parallelization.

---

## Key Modules & Entry Points

| Component | Path | Role |
|-----------|------|------|
| **Executable entry** | `src/main.cpp` | Reads args, calls `openmc_init()`, dispatches to `openmc_run()` (MC) or `openmc_run_random_ray()` (RR) |
| **Python API** | `openmc/__init__.py` | Exports ~50 modules for geometry, materials, sources, tallies, settings, execution |
| **Core particle transport** | `src/particle.cpp`, `include/openmc/particle.h` | `Particle` class; holds position, direction, energy, weight; tracks state through geometry |
| **Geometry engine** | `src/geometry.cpp`, `include/openmc/geometry.h` | CSG ray-tracing: cell location, distance-to-boundary, lattice/universe nesting |
| **Cross-section lookup** | `src/cross_sections.cpp`, `src/nuclide.cpp` | Load HDF5 ENDF libraries; sample reactions (elastic, inelastic, fission, absorption) |
| **Collision physics** | `src/physics.cpp` | Sample reaction types, handle neutron/photon/electron secondaries, scattering kernels |
| **Source sampling** | `src/source.cpp`, `include/openmc/source.h` | Abstract `Source` class; built-in distributions (spatial, angular, energy); user-extensible via `CompiledSource` |
| **Tallying** | `openmc/tallies.py`, `openmc/filter.py` | Python API to define what to score; C++ backend bins by cell/surface/energy/angle/time |
| **Simulation control** | `include/openmc/simulation.h` | Batch/generation loops, keff estimation, tally accumulation, eigenvalue convergence |
| **Settings & I/O** | `openmc/settings.py`, `src/initialize.cpp` | Parse XML input, set run mode, particle count, geometry/tallies write-out |

---

## Main Techniques & Algorithms

### Geometry & Tracking
- **CSG ray-tracing**: cells defined by Boolean combinations of surfaces (planes, spheres, cylinders, cones, tori, etc.); find particle location via coordinate hierarchy; compute distance to next boundary.
- **Lattice nesting**: rectangular/hexagonal repeated structures; cross-lattice bookkeeping for periodic tiling.
- **DAGMC (optional)**: faceted CAD geometry via Moab; alternative to CSG.

### Cross-Section Physics
- **Continuous-energy HDF5 libraries** (ENDF/B-VII.1, JEFF, JENDL): temperature-dependent interpolation; resonance unresolved (URR) probability-table sampling for energy-dependent widths.
- **Scattering kernels**: free-gas elastic, thermal (S(α,β)), Kalbach-Mann uncorrelated emission, correlated angle-energy products.
- **Fission**: prompt/delayed neutrons; Watt spectrum and energy dependence.

### Particle Transport (Monte Carlo Loop)
- **LCG RNG** (`prn(seed)`) for reproducibility & thread safety; unique seed per particle via stride.
- **Collision loop**: transport to next collision → sample reaction (type, angle, energy) → handle secondary particles → accrue weight-dependent scores.
- **Eigenvalue mode** (k-eff): cycle neutrons between source & fission bank across generations; estimate k-eff as (fission sites born) / (neutrons tracked).
- **Fixed-source mode**: single batch of independent histories from external source.

### Tallying
- **Flexible binning**: cell/surface tally; filters by energy, angle, time, material, universe level.
- **Reaction scoring**: reaction rates, heating, flux, etc.; per-nuclide & per-delayed-group resolution.
- **Statistics**: accumulate mean & variance per batch; combine across batches for final error estimate.

### Parallelization
- **MPI**: distribute particle histories across ranks; synchronize tally results & keff at batch end.
- **OpenMP**: per-rank thread parallelism within a batch; particle-loop level (each thread gets independent RNG seed).
- **Event-based**: (optional) all particles share a global collision queue; better load-balance on sparse geometries.

---

## Topic → File Index

| Question / Topic | Where to Look |
|------------------|----------------|
| **Particle state & data layout** | `include/openmc/particle_data.h` (struct `SourceSite`, `ParticleData`), `include/openmc/particle.h` (class `Particle`) |
| **Geometry definition (Python)** | `openmc/geometry.py` (class `Geometry`), `openmc/cell.py`, `openmc/surface.py`, `openmc/universe.py` |
| **Compiled custom source** | `include/openmc/source.h` (class `Source`), `openmc/source.py` (factory & `CompiledSource`); examples: `examples/custom_source/`, `examples/parameterized_custom_source/` |
| **Material & cross-section setup** | `openmc/material.py`, `openmc/nuclide.py`, `src/cross_sections.cpp` |
| **Running a simulation** | `openmc/executor.py` (`_process_CLI_arguments`, `_run`), `openmc/model.py` (class `Model.run()`) |
| **Collision & reaction sampling** | `src/physics.cpp` (functions `collision()`, `sample_neutron_reaction()`, `scatter()`, etc.) |
| **RNG & reproducibility** | `include/openmc/random_lcg.h` (`prn()`, `init_seed()`), `src/random_lcg.cpp` |
| **Source bank initialization** | `src/simulation.cpp` (functions `allocate_banks()`, `initialize_batch()`, `sample_source()`) |
| **Tally accumulation** | `openmc/tallies.py` (class `Tally`, `Filter`), `src/tallies/` (C++ backend scoring) |
| **Settings & XML I/O** | `openmc/settings.py`, `src/settings.cpp`, `src/initialize.cpp`, `src/xml_interface.cpp` |
| **Eigenvalue & keff estimation** | `src/eigenvalue.cpp`, `src/simulation.cpp` (finalize-batch logic) |
| **Parallel execution (MPI/OpenMP)** | `include/openmc/message_passing.h`, `src/message_passing.cpp`, `src/openmp_interface.cpp` |
| **Version & build info** | `include/openmc/version.h`, `CMakeLists.txt` |
| **Test suite & examples** | `tests/`, `examples/` (e.g., `examples/pincell/`, `examples/c5g7/`) |

---

## Common Call Chains

- **User runs Python model** → `openmc.Model.run()` → subprocess calls C++ executable → `main()` → `openmc_init()` → reads XML → `openmc_run()` → batch loop → transport particles → score tallies → write statepoint HDF5
- **Custom source sampling** → user defines `CompiledSource` C++ class → Python loads `.so` and calls `openmc_create_source(param_string)` → factory returns heap object → each particle birth calls `sample(seed)` method
- **Collision physics** → `transport_to_collision()` → `collision()` → `sample_nuclide()` → `sample_reaction()` → call appropriate physics routine (e.g., `scatter()`, `absorption()`, `fission()`) → create secondaries or terminate

---

**Key design principle**: Particle immutability in the hot loop (no mutable shared state except atomic counters); deterministic RNG seeding for reproducibility; tallying orthogonal to transport (flexible filter composition).
```

---

This document is now ready for the meeting copilot to ingest. It covers what OpenMC does, the major components and where to find them, the algorithms in use, and a quick lookup table for common questions about specific features. Use it with `/deep-research` or just point the copilot to this file when you need fast answers about the codebase.