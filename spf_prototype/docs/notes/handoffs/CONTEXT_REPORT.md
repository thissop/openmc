# CONTEXT REPORT — Spin-Polarized Fusion (SPF) neutron source prototype

Snapshot for a collaborator scoping the next phase (target: 3D stellarator,
Helios-class, quasi-axisymmetric). Read-only investigation; nothing in the code
was modified to produce this. All nontrivial claims cite `file:line`. Every
claim below was independently re-read/fact-checked against the source; the few
caveats and corrections found are called out explicitly.

- Repo: `/Users/tkiker/Documents/GitHub/openmc` (a full OpenMC checkout); the
  prototype lives entirely in `spf_prototype/`.
- Branch: `spf-prototype`; latest prototype commit `bbf80ada1`.
- This is a **research prototype**, honest about scope: it builds and *verifies*
  an anisotropic birth-direction source and recovers the Schwartz analytic NWL;
  it is **axisymmetric (tokamak-limit) only**. There is **no stellarator / 3D /
  general-field code yet** (see §5, §6).

---

## 1. Repo map

Tree (depth ~3, build/cache/data blobs omitted):

```
spf_prototype/
  src/
    spf_sampler.hpp              # header-only pure-math angular sampler (the core)
    polarized_fusion_source.cpp  # OpenMC CompiledSource wrapper (position + B̂ + glue)
    standalone_driver.cpp        # tier-1 CLI driver (sampler math only, no libopenmc)
    CMakeLists.txt               # builds the .so via find_package(OpenMC)
  python/
    spf_mirror.py                # 1:1 Python mirror of the header + bit-exact PCG prn port
    abc_modes.py                 # (a,b,c) mode algebra + sympy identities
    analytic_nwl.py              # independent analytic NWL (quad + anarrima + Eq.5)
    spf_stats.py                 # chi-square / KS helpers
    build_and_run.py             # tier-2b: build .so, square-torus model, wall extraction
    reactor_model.py             # tier-5/6: layered box-torus (W/steel/Be/FLiBe), tallies
    verify_sampler.py            # tier-1 diagnostics -> RESULTS_tier1.md
    verify_nwl_analytic.py       # tier-2a diagnostics -> RESULTS_tier2a.md
    verify_nwl.py                # tier-2b OpenMC-vs-analytic -> RESULTS_tier2b.md
    verify_reactor.py            # tier-5/6 runs + figures -> RESULTS_tier5/6.md
    tier3_artifacts.py           # tier-3 postprocessing (no new runs) -> RESULTS_tier3.md
    tier7_absolute.py            # tier-7 absolute-rate postprocessing -> RESULTS_tier7.md
  tests/
    conftest.py                  # builds the standalone C++ driver (g++), fixtures
    test_abc_modes.py            # sympy identities + corner cross sections (7 tests)
    test_analytic_nwl.py         # quad vs anarrima vs Eq.5, oracles (33 tests)
    test_sampler_stats.py        # C++/Py parity, chi2, KS, steering (16 tests)
  figs/  *.png                   # 8 result figures (committed)
  slides/spf_talk.tex            # 52-slide teaching deck
  peterson_figures.tex/.pdf      # 4-page figure summary for the PI
  RESULTS_tier{1,2a,2b,3,5,6,7}.md
  README.md NOTES_orientation.md LIMITATIONS.md NEXT_STEPS.md PLAN_next_phase.md
  HANDOFF.md PRIOR_WORK_BAE2025.md BLANKET_PARAMS.md
  (2 reference PDFs, gitignored)
```

**LOC / languages.** ~3,170 lines of code total: Python ~2,668 (`python/` +
`tests/`), C++17 ~506 (`src/*.cpp` + `*.hpp`); plus ~837 lines LaTeX. Primary
languages: **Python** (analysis, verification, transport driving) and **C++17**
(the sampler core + the OpenMC source plugin). Largest modules:
`tier3_artifacts.py` (361), `verify_reactor.py` (269), `spf_sampler.hpp` (264),
`analytic_nwl.py` (260), `reactor_model.py` (235).

**Build / run / test entry points (what actually exists):**
- **No** `pyproject.toml`, `setup.py`, `Makefile`, or `requirements.txt` anywhere
  in `spf_prototype/` (verified). Dependencies are assumed pre-installed in a venv.
- **C++ `.so` build**: CMake (`src/CMakeLists.txt`), invoked programmatically by
  `build_and_run.build_so()` (`build_and_run.py:36-41`: `cmake -B build
  -DCMAKE_PREFIX_PATH=$HOME/spf_venv` then `cmake --build`).
- **C++ tier-1 driver build**: a *separate* path — `g++ -O2 -std=c++17` in
  `tests/conftest.py:19-29`, linking `src/random_lcg.cpp` (not libopenmc).
- **Tests**: `pytest spf_prototype/tests/`.
- **Run entry points**: the `verify_*.py` and `tier{3,7}*.py` scripts each have a
  `__main__` and write a `RESULTS_tier*.md` + figures; `build_and_run.py` has a
  `_selftest` `__main__` (`build_and_run.py:116-133`). No CLI/notebooks.
- Built artifacts currently present on disk: `src/build/libpolarized_fusion_source.so`
  and `build/spf_driver`.

---

## 2. The sampler (the core asset)

### Where it is
- **Pure-math core**: `src/spf_sampler.hpp` — header-only, `namespace spf`,
  depends only on `<cmath>/<cassert>/<stdexcept>`, **no OpenMC types** (it defines
  its own `spf::Vec3` at `spf_sampler.hpp:52`). This deliberate isolation lets it
  be (a) compiled into the no-libopenmc tier-1 driver and (b) mirrored in Python.
- **OpenMC plugin**: `src/polarized_fusion_source.cpp` — `class
  PolarizedFusionSource : public openmc::Source` (`:74`), `sample()` override
  (`:118`), dlopen factory `extern "C" openmc_create_source(std::string)` (`:185`).
- **Python mirror**: `python/spf_mirror.py` — 1:1 port of the header **plus a
  bit-exact port of OpenMC's PCG-RXS-M-XS `prn`** (`Rng` class, `spf_mirror.py:27-41`;
  constants verified against `src/random_lcg.cpp`).

### Polarization schemes / parameterization
Parameterized by the three **Schwartz collision-mode fractions** `(a,b,c)` (Eq. 1).
The birth intensity on the sphere (`spf_sampler.hpp:15-20`):

```
I(θ) = w_perp·sin²θ + w_par·(¼ + ¾cos²θ),   w_perp = ¾a,   w_par = ⅔b + ⅓c
```

θ is the angle to the **local** field B̂. Modes used throughout:
`iso=(⅓,⅓,⅓)`, `A=(1,0,0)`, `B=(0,1,0)`, `C=(0,0,1)`, `mixed=(0.5,0.3,0.2)`
(`analytic_nwl.py:227-233`). Note there are only **two distinct angular shapes**:
`sin²θ` (A) and `¼+¾cos²θ` (the B/C "par" shape); B and C share the shape and
differ only in rate. This is **not** the Bae/Hupin–Navrátil vector/tensor `(a,b,c)`
parameterization — it is the Schwartz collision-mode form. The total-rate factor
`η = a + ⅔b + ⅓c` (`spf_sampler.hpp:119`).

### How θ is sampled — stratified inverse-CDF (closed forms derived in-code)
With probability `P_perp = a/η` pick the `sin²θ` shape, else the `¼+¾cos²θ` shape
(`sample_local_direction`, `spf_sampler.hpp:196-210`; `P_perp` set in
`make_mode_weights:117-126`). Each shape uses an **inverse-CDF whose cubic and root
selection are derived in the comments**, not just pasted:

```cpp
// spf_sampler.hpp:141-152  (sin²θ shape: p(x)∝1-x²; solves x³-3x+(4u-2)=0)
inline double sample_costheta_perp_invcdf(double u) {
  double x = 2.0 * std::cos(std::acos(1.0 - 2.0 * u) / 3.0 - 2.0 * PI / 3.0);
  ... clamp to [-1,1] ...
}
// spf_sampler.hpp:154-168  (¼+¾cos²θ shape: solves x³+x+(2-4u)=0, Cardano)
inline double sample_costheta_par_invcdf(double u) {
  const double q = 2.0 - 4.0 * u;
  const double disc = q*q/4.0 + 1.0/27.0;          // > 0 always
  double x = std::cbrt(-q/2.0 + std::sqrt(disc)) + std::cbrt(-q/2.0 - std::sqrt(disc));
  ...
}
```

Rejection-sampling variants of both shapes exist (`spf_sampler.hpp:172-189`) but
are used **only as a statistical cross-check** (KS test); the production path is
inverse-CDF (`use_rejection=false` default, `spf_sampler.hpp:197,244`). The
solid-angle integrals behind `P_perp` (`∫sin²θ dΩ=8π/3`, `∫(¼+¾cos²θ)dΩ=2π`) are
verified symbolically in the tests, not hardcoded blind.

### Local B̂ and lab-frame construction
B̂ is obtained **in the plugin**, not the header (`polarized_fusion_source.cpp:146-153`):
- `bmode="toroidal"` (default): `B̂ = φ̂ = (-py, px, 0)/rxy` — a real
  *position-dependent* field, but only in the trivial axisymmetric-toroidal sense.
- `bmode="constant"`: a fixed user vector `(bx,by,bz)`.
- **Only these two** are accepted; anything else throws (`:96-97`). **There is no
  general 3D field / field-map / equilibrium path.**

The lab-frame triad is built by **Gram–Schmidt off the least-aligned world axis**
(`build_frame`, `spf_sampler.hpp:220-234`), which avoids the `B̂×ẑ` degeneracy when
`B̂∥ẑ`; `sample_global_direction` asserts the result is unit (`:248`). Important for
the next phase: **the header math already accepts an arbitrary B̂** (tested with
off-axis and degenerate ±ẑ fields, §4), so only the *plugin's field model* is the
tokamak-limit bottleneck, not the rotation math.

### Spatial sampling (the parabolic ring/plasma)
Two shapes (`polarized_fusion_source.cpp:101-144`):
- `shape="ring"`: filamentary — fixed `r0,z0`, φ uniform (`:124-128`).
- `shape="plasma"`: a **parabolic `(1-ρ²)` profile, volume-weighted**, rejection
  sampled over the poloidal disk with φ uniform (`:131-144`):

```cpp
// polarized_fusion_source.cpp:133-139
for (;;) {
  p   = (R0_ - aminor_) + 2.0*aminor_*rng();
  z0w = -aminor_ + 2.0*aminor_*rng();
  rho2 = ((p-R0_)*(p-R0_) + z0w*z0w)/(aminor_*aminor_);
  if (rho2 <= 1.0 && rng() <= (1.0 - rho2) * p / pmax_) break;  // (1-ρ²)·p weight
}
```

The extra `·p` is the volume Jacobian (`dV = p dp dz dφ`). The analytic side
mirrors this with a grid of rings weighted by `(1-ρ²)` (`analytic_nwl._ring_grid`,
`analytic_nwl.py:197-206`); the docstring notes the `p` of `dV` cancels the `1/p`
in the kernel so no extra `p` weight appears there (`analytic_nwl.py:19-20`).

### Energy
**Monoenergetic, default 14.1 MeV** (`polarized_fusion_source.cpp:115,167`). There
is **no Ballabio broadening and no D-D channel** anywhere in the C++ (grep-verified).

### Rate handling (important design point)
`sample()` **always sets `wgt = 1.0`** (`polarized_fusion_source.cpp:168`). The
rate factor `η` is computed in `ModeWeights` but is **never applied inside the
source** — every mode emits at constant total rate. Absolute-rate effects (the
+50%/−50% A/C story) are applied **only in Tier-7 postprocessing**
(`tier7_absolute.py`). A collaborator wanting absolute neutron yields must add η
to `strength`/weighting; the source does not do it.

### Transport connection
In-process **OpenMC `CompiledSource`** via the dlopen factory; used at
`build_and_run.py:58` and `reactor_model.py:95`. **No `FileSource`/MCPL** path
exists (grep-verified). Samples are generated on the fly by OpenMC.

---

## 3. Schwartz recovery / analytic validation (what is actually verified)

### Tier 2a — analytic, three independent engines
`python/analytic_nwl.py` implements three independent computations of the
single-ring reduced intensity `g`:
1. `g_ring_quad` — our vectorized 256-node Gauss–Legendre quadrature (`:61,:64`).
2. `g_ring_anarrima` — **Schwartz's published closed forms via the `anarrima`
   package** (the reference; `:117`, dispatch `_an_fn:109-114`).
3. `g_RiI_closed_eq5` — our own elliptic-integral transcription of Eq. 5
   (inboard-isotropic case only; `:147`).

`anarrima` **is installed (0.1.0) and actually used** (not just hand-coded).

### Geometry / source / cases checked
- **Square cross-section torus**, `R0=1`, `a=0.5`, walls `u=0.4`, `w=1.6`,
  `z=±0.6` (`analytic_nwl.py:169-176`). Aspect ratio is **2.0 to match the
  `anarrima` example**; the Schwartz paper *text* says 2.5 — this discrepancy is
  documented in both `analytic_nwl.py:164-171` and `RESULTS_tier2a.md:3`.
- Sources: **single central ring** (tightest check) and the **parabolic
  multi-ring plasma**.
- Walls: **inboard, outboard, floor** (+ ceiling by up-down symmetry). Modes:
  `iso`, `A`, `cos2`, `BC` and combined `iso/A/B/C/mixed`.

### Quantitative agreement (actual numbers, from `RESULTS_tier2a.md`)
- Single ring, three engines: **max |Δ| ≈ 1.8e-15** across all walls/factors.
- Parabolic plasma, quad vs anarrima: **max relative diff ≈ 5.2e-08**.
- Schwartz §2 scalar oracles reproduced (our quad == anarrima; paper text rounds):
  A inboard **+40.6%** (paper "+43%"), A outboard **−21.2%** ("−22%"), iso
  outboard/inboard **+14.6%** ("+12%").
- Tests enforce: quad-vs-anarrima `<1e-9`, Eq.5 `<1e-9`, plasma `<1e-6`
  (`test_analytic_nwl.py:24,31,49`).

### Tier 2b — full OpenMC transport recovery
`python/verify_nwl.py`: square torus, **interior void (free-streaming)**, cylindrical
mesh surface-current tally, 4e6 histories/config, 50×50 bins. It compares
**directionality (mode/iso ratios)**, which is normalization-free, against the
exact constant-rate analytic `D = bracket/η`.

Result (`RESULTS_tier2b.md`): the standardized residual `(MC−analytic)/σ` has
**mean ≈ 0** (−0.022…+0.019) but **stdev ≈ 0.68, not 1.0**; the Schwartz scalar
oracles are reproduced from OpenMC (A inboard +39.2% vs analytic +40.6% vs paper
+43%); B≡C and linearity hold.

> **⚠ Honesty flags on Tier 2b (verified):**
> 1. The σ≈0.68 means the residual is **under-dispersed** (the MC agrees with the
>    analytic *better* than the tally error bars imply). The "mildly conservative
>    error bars" explanation exists, but in `LIMITATIONS.md:54-57` and
>    `RESULTS_tier3.md` — **not** in `RESULTS_tier2b.md`/`verify_nwl.py`, which
>    instead state the pass criterion as "stdev≈1" (`verify_nwl.py:123`,
>    `RESULTS_tier2b.md:18`) without reconciling the observed 0.68. The "N(0,1)"
>    framing is therefore loose; the real statement is "mean≈0, σ≈0.68".
> 2. The **B≡C check reports mean/stdev exactly 0.000** (`RESULTS_tier2b.md:32`)
>    because it is **trivially exact, not an independent statistical result**: pure
>    B and pure C both have `a=0 ⇒ P_perp=0`, so both *always* take the par branch
>    and consume the RNG in the identical order; with **no explicit seed set**
>    (OpenMC default master seed = 1; `build_and_run.make_model`/`verify_nwl` set
>    none) the two runs produce **bit-identical particle streams**. It confirms the
>    shape identity but not statistical convergence.

### What is NOT validated against analytics
- **Scattering NWL (Tier 5)**: no analytic ground truth. It is anchored only by a
  **thin-wall limit** (all densities ×1e-4) that must recover the Tier-2b
  free-streaming result (`RESULTS_tier5.md:8-17`); the scattering effect itself is
  explicitly "the kind of correction the analytic theory cannot provide" (`:34`).
- **3D / non-axisymmetric**: nothing.
- **Angled / general-B fields**: not validated. `anarrima` *ships* the angled-field
  kernels (`g_HAa/g_Hca/g_HBa`, `g_VAa/...`) — confirmed present in the installed
  package — but they are **referenced only in planning docs and never called** in
  any `.py` (`_an_fn` wires only the toroidal kernels, `analytic_nwl.py:109-114`).
- A vs B/C: both covered, but only in the free-streaming regime.

---

## 4. Tests & reproducibility

**Suite: 56 tests, all pass, ~27 s** (re-run during this review).
Command: `PATH=$HOME/spf_venv/bin:$PATH python -m pytest spf_prototype/tests/ -q`.
Breakdown: `test_abc_modes.py` **7**, `test_analytic_nwl.py` **33**,
`test_sampler_stats.py` **16**. **No skips, no smoke tests** — every test asserts a
physical identity, a numerical agreement, or a bit-parity:

- `test_abc_modes.py`: sympy-verified integrated cross section `= σ0(a+⅔b+⅓c)`;
  corner totals (2/3,1,2/3,1/3); B/C share shape; weight-map closed forms; Eq.1
  collision-mode fractions.
- `test_analytic_nwl.py`: quad-vs-anarrima single ring (`<1e-9`, 6 walls × 4
  factors); Eq.5 (`<1e-9`); free identities (`g_iso=g_A+g_cos2`,
  `g_BC=¼g_iso+¾g_cos2`); plasma quad-vs-anarrima (`<1e-6`); Schwartz scalar
  oracles; B≡C and iso≢C; linearity.
- `test_sampler_stats.py`: **C++/Python bit-parity ≤1e-12** (`:45`); per-mode cosθ
  chi-square vs analytic pdf; φ uniformity; rotation isotropy; off-axis steering
  `⟨(u·B̂)²⟩`; degenerate `B̂=±ẑ`; inverse-CDF-vs-rejection KS.

**Coverage gap**: the OpenMC **transport tiers (2b/5/6/7) have no pytest** — they
are driven by the `verify_*.py` / `tier*.py` scripts that write `RESULTS_tier*.md`
(confirmed by reading those scripts: e.g. `verify_nwl.py` writes
`RESULTS_tier2b.md`, `verify_reactor.py` writes tier 5/6, `tier7_absolute.py`
writes tier 7). So the headline transport results are **not under CI-style
regression protection**.

**RNG / seeding.** OpenMC's PCG `prn`, ported bit-for-bit in `spf_mirror.Rng`. Unit
tests use **fixed seeds** and assert C++/Python parity ≤1e-12 → reproducible. The
**transport scripts set no explicit seed** and rely on OpenMC's default master seed
(=1); they are reproducible run-to-run but this is implicit, and it is the reason
the B≡C check is bit-exact (§3).

**Dependencies installed in this environment** (venv `~/spf_venv`):
numpy 2.5.0, scipy 1.18.0, sympy 1.14.0, matplotlib 3.11.0, pandas 3.0.3,
h5py 3.16.0, **openmc 0.15.4.dev182+g608a1c338**, **anarrima 0.1.0**, jax 0.10.2,
jaxlib 0.10.2, smplotlib 1.0.0, pytest 9.1.1. Required for the core: numpy/scipy/
sympy; for the analytic reference: anarrima+jax; for transport: openmc+pandas+h5py;
matplotlib/smplotlib for figures.
**Cross sections**: `OPENMC_CROSS_SECTIONS` is **unset**, but the scripts set
`openmc.config["cross_sections"] = ~/nndc_hdf5/cross_sections.xml`
(`build_and_run.py:25`, `reactor_model.py:26`) and that file **exists**. So XS are
configured *in-script*, hardcoded to one user path.

---

## 5. Geometry & transport status

- **Geometry: OpenMC CSG only.** Tier-2b is a single **void cell** square torus
  from `ZCylinder`/`ZPlane` (`build_and_run.py:46-50`). Tier-5/6 is a layered
  box-torus: void cavity + 4 material layers (W → Fe/Cr steel → Be → FLiBe),
  `reactor_model.py:62-87`. **No DAGMC, no CAD.**
- **No stellarator / DESC / VMEC / equilibrium code anywhere** (grep-verified).
  The only `DAGMC` hit is a caption string that *disclaims* CAD geometry
  (`tier3_artifacts.py:350`); the only `desc` hits are CMake build-artifact
  substrings. Everything is **axisymmetric**; the field is toroidal or constant.
- **Transport runs do exist** (this is not source-construction-only): Tier-2b
  free-streaming surface currents, and **Tier-5/6 full scattering** with materials
  and tallies for `heating`, `damage-energy`, `H3-production`, and a FLiBe flux
  spectrum (`reactor_model.py:110-125`). Tier-7 is pure postprocessing of the
  Tier-5/6 statepoints. Reactor scale: `R0=200, a=100 cm`, cavity `U=80, W=320,
  ZC=120 cm`, layers `0.1/1/1/40 cm`, FLiBe `Li2BeF4` @ 30% ⁶Li, vacuum outer
  boundary (`reactor_model.py:29-58`).

---

## 6. State, gaps, and risks (candid)

**What works end-to-end today.** Given `(a,b,c)`, the prototype builds a compiled
OpenMC source that emits Schwartz-anisotropic 14.1 MeV neutrons about a toroidal (or
constant) B̂ from a parabolic ring/plasma, runs it through an axisymmetric CSG torus,
and tallies wall load / heating / damage / TBR. The **angular sampler is verified to
machine precision** (C++/Python parity ≤1e-12; symbolic identities) and the
**free-streaming NWL is recovered against an independent analytic** (anarrima to
~1e-6 analytically, and to statistics in transport). A FLiBe-blanket TBR and a
W-first-wall heating study exist with a thin-wall free-streaming anchor.

**Top gaps to a 3D stellarator SPF Monte Carlo study:**
1. **Real B̂(x).** Only `toroidal`/`constant` exist (`polarized_fusion_source.cpp:92-153`).
   A stellarator needs a general 3D field from an equilibrium (DESC/VMEC/Boozer or a
   field map). The **header rotation already handles any B̂** (§2), so the work is a
   field-model interface in the plugin + the spatial source region — not the math.
2. **3D geometry.** Everything is axisymmetric CSG; stellarators are inherently 3D
   and non-axisymmetric. The whole CSG geometry **and** the axisymmetric analytic
   NWL framework (`analytic_nwl.py`) **do not transfer** to 3D. A DAGMC/CAD route is
   needed; there is none today.
3. **No analytic ground truth off-axis/3D.** The recovery anchor is axisymmetric and
   free-streaming. The angled-field anarrima kernels are a *ready but unused* hook
   for pitched-field axisymmetric checks; there is **no** 3D analytic check at all.
4. **Tallies vs. real device.** Heating/damage/TBR tallies exist but on a toy square
   blanket; there is **no coil model** for coil-damage, and the blanket is single-zone.

**Bugs / hacks / hardcoded values / stability (all verified):**
- **No `TODO`/`FIXME`/`HACK`/`XXX` markers** in `python/`, `src/`, `tests/`.
- **Tier-2b residual σ≈0.68 framed as "stdev≈1"** in the tier-2b code/results; the
  honest explanation lives in `LIMITATIONS.md`/`RESULTS_tier3.md` (§3). Not a bug,
  but a loose pass-criterion a collaborator should not take literally.
- **B≡C "exact" agreement is a bit-identical-RNG artifact**, not statistics (§3).
- **Hardcoded environment paths**: `VENV=~/spf_venv`, `XS=~/nndc_hdf5/cross_sections.xml`
  (`build_and_run.py:24-25`). Will break on any other machine.
- **Statepoints written to ephemeral `/tmp/spf_*`** (`verify_nwl.py:54`,
  `verify_reactor.py:50`, etc.) and the postprocessing/figure scripts **read** those
  paths (`tier3_artifacts.py:44`, `tier7_absolute.py:43`, `verify_reactor.py
  --replot`). If `/tmp` is cleared, **figures cannot be regenerated without re-running
  the transport** — there is no cached intermediate.
- **Aspect-ratio 2.0 vs paper's 2.5** — documented, deliberate (matches anarrima),
  but a surprise if compared to the paper text.
- **Simplified "steel"**: Fe(0.91 wt)/Cr(0.09 wt) only, labeled "Eurofer-like,
  data-safe" (`reactor_model.py:44-46`) — omits W/Mn/V/Ta/C.
- **Temperature inconsistency**: 294 K cross sections used with a 900 K FLiBe density
  (stated caveat, `RESULTS_tier5.md:66`, `RESULTS_tier6.md:42`).
- **Fragile build**: the `.so` needs OpenMC **`make install`ed into the venv prefix**
  (`find_package(OpenMC)` + `OpenMC::libopenmc`, `CMakeLists.txt:10-11`); a build-tree
  config is not consumable. Two separate C++ build paths exist (CMake `.so` and g++
  driver) for the same header.

**Surprising / nonstandard, worth knowing:**
- "isotropic" in the figures/results = **unpolarized fuel baseline**, not a modeling
  error (the polarized cases are A/B/C).
- The sampler is intentionally OpenMC-free + Python-mirrored with a bit-exact PCG
  port — unusually rigorous for a prototype; lean on it.
- No packaging/Makefile; the "build system" is a Python function calling CMake.

---

## 7. Reusable assets

**Build on directly (verified, well-tested):**
- `src/spf_sampler.hpp` — the angular sampler. **Already handles an arbitrary B̂**
  (Gram–Schmidt, tested off-axis and degenerate). The single most reusable piece for
  the 3D phase: feed it a real B̂(x) and it works unchanged.
- `python/spf_mirror.py` — Python sampler + bit-exact PCG port; the **oracle** for
  validating any future sampler change.
- `python/abc_modes.py` — `(a,b,c)` algebra + sympy identities (`total_xs`,
  `int_*_solid_angle`, `angular_shape_b_minus_c`).
- `python/spf_stats.py` — generic chi-square / KS helpers.
- `python/reactor_model.py` extraction helpers `tbr()` (`:199`), `cell_table()`
  (`:176`), `midplane_intensities()` (`:160`) — tally-reading utilities, geometry-
  agnostic-ish.

**Refactor / replace before the next phase:**
- `polarized_fusion_source.cpp` **field model** — only toroidal/constant; needs a
  general B̂(x) / field-map / equilibrium interface. The **spatial sampler** is also
  tied to a circular ring/plasma; a stellarator needs a 3D source region.
- **Geometry** (`build_and_run.py`, `reactor_model.py`) — axisymmetric CSG; replace
  with DAGMC/3D for a stellarator.
- `python/analytic_nwl.py` — reusable **only** for axisymmetric free-streaming; the
  anarrima angled kernels are a ready hook for pitched-field (still axisymmetric)
  checks, but there is no 3D analytic and this module won't validate 3D.
- **Hardcoded paths + `/tmp` statepoints + no packaging** — make configurable and
  add `pyproject`/`requirements` if this is to be shared or run in CI.
- The `.so` build's `make install` dependency — document or vendor for portability.

---

## Open questions for the collaborator (could not be determined from the code)

1. **Field representation**: what will the stellarator equilibrium hand us for
   B̂(x) — a DESC output object, a Boozer-coordinate field, a sampled field map?
   This sets the plugin's new field interface (the header math is ready either way).
2. **Source spatial model in 3D**: the current parabolic ring/plasma is axisymmetric;
   what defines the 3D emissive region (flux-surface profile on a stellarator
   equilibrium)? Out of scope of anything present.
3. **Does any analytic check survive into 3D**, or is the Schwartz/anarrima recovery
   purely a tokamak-limit unit test to be retained for regression? (No 3D analytic
   ground truth exists.)
4. **Energy spectrum scope**: is mono-14.1 MeV acceptable, or are Ballabio broadening
   / D-D needed? (Mentioned in planning docs; not implemented.)
5. **Absolute rate**: should the source itself carry η (so absolute neutron yields /
   power come straight out), or keep the current "directionality-only, η in
   postprocessing" split?
6. **Transport geometry route**: DAGMC-from-CAD, or hand-built 3D CSG approximations?
7. **Coil/structure models**: needed for coil-damage tallies — none exist yet.
```
