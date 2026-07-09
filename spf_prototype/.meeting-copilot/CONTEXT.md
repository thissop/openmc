# Quick-Look Context: SPF Neutron Source Prototype

## What This Is

A research prototype implementing **spin-polarized DT neutron birth-direction sampling** in OpenMC, validated against closed-form analytic results. The core is a stratified inverse-CDF sampler (C++17 header + Python mirror) parameterized by Schwartz collision-mode fractions `(a,b,c)`, which controls anisotropic angular emission via two basis functions (`sin²θ` and `1/4+3/4cos²θ`). Verification tiers proceed from direction-sampler statistics (Tier 1) → free-streaming NWL geometry (Tier 2, validated against Schwartz/anarrima) → realistic source physics (Tier 4, planned) → scattering & blanket (Tier 5–6, planned). Current scope: **axisymmetric only, free-streaming, monoenergetic 14.1 MeV, toroidal or constant-direction B-field.**

---

## Entry Points & Key Modules

| Purpose | Files |
|---------|-------|
| **Pure-math sampler core** | `src/spf_sampler.hpp` (264 lines) — header-only, zero OpenMC deps, defines `sample_local_direction()`, inverse-CDF cubics, frame rotation |
| **OpenMC CompiledSource plugin** | `src/polarized_fusion_source.cpp` — wraps header, parses param string, implements `sample()` override, dlopen factory `openmc_create_source()` |
| **Python mirror + RNG** | `python/spf_mirror.py` — 1:1 port of sampler + PCG-RXS-M-XS RNG (bit-exact to OpenMC) |
| **Mode algebra & identities** | `python/abc_modes.py` — Schwartz mode decomposition, total-rate factor η, sympy verification |
| **Analytic NWL ground truth** | `python/analytic_nwl.py` (260 lines) — three independent engines: scipy quadrature, anarrima package (Schwartz's own code), direct Eq.5 evaluation |
| **Tier 1 direction tests** | `tests/test_sampler_stats.py` — χ² histogram, KS tests, rotation isotropy, C++/Python parity |
| **Tier 1 mode identities** | `tests/test_abc_modes.py` — sympy integral verification, corner cross sections |
| **Tier 2b NWL vs OpenMC** | `python/verify_nwl.py` — builds square-torus geometry, runs free-streaming, residual normality |
| **Tier 3–7 postprocessing** | `python/verify_reactor.py`, `python/tier3_artifacts.py`, `python/tier7_absolute.py` — wall-load maps, peaking factors, absolute-rate scaling |
| **Spatial source models** | `python/reactor_model.py` — box-torus, W/steel blanket, tallies (current, heating, TBR) |
| **Build & run harness** | `python/build_and_run.py` — CMake invocation, square-torus runner, surface-current extraction |
| **Unit test harness** | `tests/conftest.py` — standalone C++ driver build (g++, no libopenmc), fixtures |

---

## Core Algorithms & Design

### Angular sampling (stratified inverse-CDF)

- **Two basis shapes:** `sin²θ` (A mode, perpendicular to B̂) and `1/4 + 3/4cos²θ` (B/C modes, parallel to B̂).
- **Stratified:** pick shape with probability `P_perp = a/η`, then inverse-CDF within the chosen shape.
- **Cubic solvers:** depressed-cubic for both; roots derived in code comments (not magic constants).
  - sin²θ: solves `x³ − 3x + 2(1−2u) = 0` (trigonometric method).
  - par: solves `x³ + x + (2−4u) = 0` (Cardano).
- **Cross-check:** rejection-sampling variants exist (KS-tested against CDF, used for unit tests only).
- **Solid-angle integrals:** `∫sin²θ dΩ = 8π/3`, `∫(1/4+3/4cos²θ) dΩ = 2π` — sympy-verified, hardcoded with derivation comments.

### Frame rotation (Gram-Schmidt, no gimbal lock)

- Pick reference axis (x, y, or z) with smallest `|B̂·ê|`.
- Gram-Schmidt: `x̂_local = normalize(ref − (ref·B̂)B̂)`, `ŷ_local = B̂ × x̂_local`, `ẑ_local = B̂`.
- Tested on degenerate cases (B̂ = ±ẑ).

### Spatial sampling

- **Ring (filamentary):** uniform φ, fixed (r₀, z₀).
- **Plasma (parabolic):** `(1−ρ²)` rejection-sampled over poloidal disk; volume Jacobian `p` accounted (see `CONTEXT_REPORT.md:169-172`).

### Rate handling (crucial design point)

- **Sampler sets `wgt = 1.0` always** — all modes emit at constant rate.
- Total-rate factor `η = a + 2b/3 + c/3` computed but **never applied in C++**. Applied only in post-processing (Tier 7).
- Allows mode-wise rate manipulation at the OpenMC config level, not in source.

---

## Verification Tiers (What's Done vs. Planned)

| Tier | What | Status | Entry point |
|------|------|--------|-------------|
| **1** | Direction sampler alone (stats, C++/Py parity, frame rotation) | ✅ Green | `tests/test_sampler_stats.py`, `python/verify_sampler.py` |
| **2a** | Analytic NWL (3 independent engines: quad, anarrima, Eq.5) | ✅ Green | `tests/test_analytic_nwl.py`, `python/verify_nwl_analytic.py` |
| **2b** | Free-streaming NWL: OpenMC vs analytic (residual normality, isotropic≡C) | ✅ Green | `python/verify_nwl.py`, `RESULTS_tier2b.md` |
| **3** | Wall-load peaking factors, inboard/outboard fractions | ✅ Green | `python/tier3_artifacts.py`, `RESULTS_tier3.md` |
| **4a** | Ballabio energy broadening (planned) | — | — |
| **4b** | Angled/real field (Peterson bridge, planned) | — | — |
| **5** | Scattering + wall heating/dpa (planned) | — | — |
| **6** | Blanket + TBR (planned) | — | — |

---

## Quick Reference: Topic → File

- **"What's the angular pdf for mode X?"** → `python/abc_modes.py:24–50` (integral forms), `src/spf_sampler.hpp:15–31` (physics doc)
- **"How are cosθ values sampled?"** → `src/spf_sampler.hpp:141–168` (inverse-CDF), `python/spf_mirror.py:84–114` (Python equivalent)
- **"Does the C++ match Python on fixed seeds?"** → `tests/test_sampler_stats.py::test_cpp_vs_python_fixed_seed`
- **"How does OpenMC config pass the mode fractions to C++?"** → `src/polarized_fusion_source.cpp:83–97` (parse), `python/build_and_run.py:58–75` (invocation)
- **"What's the validation against Schwartz?"** → `python/analytic_nwl.py:227–233` (test configs), `python/verify_nwl.py:__main__` (full stack), `RESULTS_tier2b.md` (numbers)
- **"Why is Tier 7 needed?"** → `docs/notes/results/RESULTS_tier7.md:1–30` (absolute-rate effects), `python/tier7_absolute.py` (code)
- **"What are the known limitations?"** → `docs/reference/LIMITATIONS.md` (6 major categories, each justified)
- **"How do I run everything?"** → `RESULTS_tier{1,2a,2b,3,7}.md` each have a "Reproduce" section; see also `tests/conftest.py:__main__`, `python/verify_sampler.py:__main__`
- **"What's next (Tier 4–6)?"** → `docs/notes/planning/PLAN_next_phase.md:7–70` (detailed roadmap with verification strategy)
- **"Parallelism in the source?"** → `src/polarized_fusion_source.cpp:1–20` (thread-safety note), `CONTEXT_REPORT.md:Orientation` (concurrency model)

---

## Key Code Patterns & Gotchas

| Pattern | Location | Notes |
|---------|----------|-------|
| **Mode-to-weights mapping** | `spf_sampler.hpp:98–126` | `make_mode_weights(a,b,c)` returns struct with `w_perp`, `w_par`, `P_perp`, and **asserts internal invariants** (sum checks, probability bounds) |
| **RNG interface** | `spf_sampler.hpp:196–210` | Sampler is templated on `Rng` callable; `openmc::prn(seed)` at tier 2, PCG port in Python. **Reproducibility depends on fixed seed sequence.** |
| **Analytic dispatcher** | `python/analytic_nwl.py:109–114` | Tries anarrima first; falls back to scipy quad; both tested to 1e-6 agreement |
| **Tallied NWL extraction** | `python/build_and_run.py:88–115` | Extracts surface-current mesh tally; normalizes to unity for shape comparison |
| **Absolute vs. relative rates** | `python/tier7_absolute.py:__main__` | Multiplies each mode's tally by η to show true rate; shows A mode is +50% vs. non-polarized (i.e., η_A=1, η_iso=2/3, ratio=3/2) |

---

## Dependencies & Environment

- **C++:** C++17, only `<cmath>/<cassert>/<stdexcept>` in the header.
- **OpenMC API:** `CompiledSource`, `SourceSite`, `openmc::prn(seed)`, `Position`, `Direction` types.
- **Python:** scipy (quad), anarrima (GitHub, not PyPI), numpy, matplotlib.
- **Build:** CMake + `find_package(OpenMC)`, or standalone `g++` for tier-1 driver.
- **Data:** ENDF/B HDF5 library (assumed pre-loaded; tests assume `OPENMC_CROSS_SECTIONS` set).