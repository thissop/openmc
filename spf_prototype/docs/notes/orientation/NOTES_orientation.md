# NOTES — §0 Orientation (ground truth for THIS checkout)

Written before any implementation. Every signature below was read directly from
this checkout **and** independently re-verified by a 6-agent search workflow;
the two agree. File:line references are to this tree.

## 0.1 Version / repo state

| Item | Value |
|---|---|
| `git rev-parse HEAD` | `608a1c3386db6af88d70d35cc2202d3955aea30b` |
| `git describe --tags` | `v0.15.3-182-g608a1c338` (OpenMC **0.15.3** + 182 commits) |
| Branch at start | `develop` (working now on **`spf-prototype`**) |
| `git status` | only `CLAUDE.md` modified (the task spec) |
| version header | `include/openmc/version.h.in` (generated → `version.h` at build) |

## 0.2 The C++ source API (verified)

### `openmc::Source` — abstract base (`include/openmc/source.h:65-120`)
```cpp
class Source {
public:
  virtual ~Source() = default;
  virtual double strength() const { return strength_; }      // optional override
  virtual SourceSite sample(uint64_t* seed) const = 0;       // <-- THE method
  ...
protected:
  virtual bool constraints_applied() const { return false; } // optional override
  double strength_ {1.0};
};
```
- The **one** method we must override: `SourceSite sample(uint64_t* seed) const` —
  returns `SourceSite` **by value**, takes a `uint64_t*` seed, is `const`.
  (`source.h:88`)
- We *may* override `strength()` (used to weight this source's spatial rate vs
  other sources) and `constraints_applied()`. We will override `strength()` to
  carry the polarization total-rate factor η when comparing absolute rates.

### Custom-source plug-in path (dlopen) — verified end-to-end
- `CompiledSourceWrapper` (`source.h:193-213`) `dlopen`s the `.so`, `dlsym`s a
  symbol named **exactly `openmc_create_source`**, calls it, and delegates
  `sample()` to the returned object. (`src/source.cpp:543-587`, `RTLD_LAZY`.)
- Factory typedef (`source.h:215`):
  ```cpp
  typedef unique_ptr<Source> create_compiled_source_t(std::string parameters);
  ```
- So our `.so` **must** export:
  ```cpp
  extern "C" std::unique_ptr<openmc::Source>
  openmc_create_source(std::string parameters);
  ```
  (`extern "C"` is required or `dlsym` won't find it; returning
  `unique_ptr<Derived>` also works via implicit upcast, but we'll return
  `unique_ptr<openmc::Source>` to match the typedef exactly.)
- `parameters` is the **verbatim string** from the Python `parameters=...` /
  XML attribute. In-tree example parses it as comma-separated `key=value`
  (`examples/parameterized_custom_source/parameterized_source_ring.cpp:15-32`) —
  we will reuse that exact parsing idiom.

### `SourceSite` — fields we must populate (`include/openmc/particle_data.h:41-58`)
```cpp
struct SourceSite {
  Position r;                       // birth position [cm]
  Direction u;                      // birth direction (unit)
  double E;                         // energy [eV]      -- NO default, MUST set
  double time {0.0};
  double wgt {1.0};
  int delayed_group {0};
  int surf_id {SURFACE_NONE};       // SURFACE_NONE == 0  (constants.h:379)
  ParticleType particle;            // defaults to neutron, but we set explicitly
  // ... parent_id/progeny_id/wgt_born/... : leave defaults
};
```
- **Must set**: `r`, `u`, `E`, `particle`. `wgt` defaults to 1.0 (the spec's
  requirement). `E` and `particle` have *no* useful default for us — set both.

### `Position` / `Direction` (`include/openmc/position.h`)
- `using Direction = Position;` (`position.h:231`) — same struct; members
  `double x,y,z` default 0.
- Ctors: `Position()`, `Position(x,y,z)`. Brace-init works: `u = {1.,0.,0.};`.
- Has `.dot(o)`, `.cross(o)`, `.norm()`, unary `-`, and free `+ - * /` with
  scalar or `Position`. **Everything the Gram-Schmidt rotation needs is built in.**

### `ParticleType` (`include/openmc/particle_type.h`) — **API differs from spec**
- It is a **class**, not a plain enum. Set neutron via the static factory
  **`openmc::ParticleType::neutron()`** (note the parentheses) — `particle_type.h:103`.
  The bare `ParticleType::neutron` the spec assumed does **not** exist.
- In-tree examples do exactly `particle.particle = openmc::ParticleType::neutron();`.

### RNG (`include/openmc/random_lcg.h`, `src/random_lcg.cpp`)
- `double prn(uint64_t* seed);` (`random_lcg.h:26`), namespace `openmc`, returns
  **[0, 1)** via `ldexp(result, -64)` (`random_lcg.cpp:43`). PCG-RXS-M-XS engine.
- Source sampling uses a dedicated stream **`STREAM_SOURCE = 1`**
  (`random_lcg.h:14`; `src/source.cpp:683` uses `init_seed(id, STREAM_SOURCE)`).
- **Thread model**: each particle carries its own `seeds_[N_STREAMS]`; `sample()`
  is called concurrently with a per-particle `seed`. ⇒ **no mutable members** in
  our source; everything set in ctor, all methods `const`, only `prn(seed)` for
  randomness. (This is the data-race hard rule from the spec — confirmed real.)

## 0.3 Python side (`openmc/source.py`)
```python
openmc.CompiledSource(library, parameters=None, strength=1.0, constraints=None)
```
- `library` (PathLike → `.so`), `parameters` (str, passed verbatim), `strength`
  (float). Serializes `library`/`parameters` as XML attributes
  (`source.py:715-762`). Attach via `model.settings.source = src` (single or list;
  `settings.py:662-667`).
- Example: `examples/parameterized_custom_source/build_xml.py:21-23`
  ```python
  source = openmc.CompiledSource()
  source.library = 'build/libparameterized_source.so'
  source.parameters = 'radius=3.0, energy=14.08e6'
  ```

## 0.4 The custom-source BUILD RECIPE (the fragile step) — verified in-tree

There are **working in-tree examples** that match our geometry almost exactly:
- `examples/parameterized_custom_source/` — a *parameterized ring source*
  (`parameterized_source_ring.cpp` + `CMakeLists.txt` + `build_xml.py`).
- `examples/custom_source/` — simpler ring.
- Regression tests `tests/regression_tests/source_dlopen/` and
  `…/source_parameterized_dlopen/` build the `.so` **at test time** and run it.

**Canonical CMake recipe** (from the regression test, `source_parameterized_dlopen/test.py:19-26`):
```cmake
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(openmc_sources CXX)
add_library(source SHARED <our_source>.cpp)
find_package(OpenMC REQUIRED HINTS <OPENMC_REPO>/build)
target_link_libraries(source OpenMC::libopenmc)
```
then `mkdir build && cd build && cmake .. && make`.

Key implications (these define the tier-2 environment gate):
- `find_package(OpenMC … HINTS <repo>/build)` requires OpenMC to have been
  **built** so that `<repo>/build/OpenMCConfig.cmake` + `OpenMCTargets.cmake`
  exist. The config also `find_dependency(fmt CONFIG)` and
  `find_dependency(pugixml CONFIG)` (`cmake/OpenMCConfig.cmake.in`) — these are
  satisfied by OpenMC's own build tree, so **a full `make install` is NOT
  required**; pointing at `build/` is what the test harness does.
- `libopenmc` is `OUTPUT_NAME openmc` → file is `libopenmc.so`; it links HDF5 +
  `fmt::fmt` + `dl` (`CMakeLists.txt:492-526`). Headers in `include/openmc`
  (138 of them) are self-sufficient for the out-of-tree compile.
- Required includes in our `.cpp`: `openmc/source.h`, `openmc/particle.h`,
  `openmc/random_lcg.h` (matches all in-tree examples).
- Docs: `docs/source/usersguide/settings.rst:378-507` ("Compiled Sources").

## 0.5 Cross-section data, near-void, tallies (for tier-2)

- **`OPENMC_CROSS_SECTIONS` is UNSET.** No `cross_sections.xml` is reachable.
  Download recipe (`tools/ci/download-xs.sh`): NNDC HDF5 (~800 MB) from
  `https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz` →
  `$HOME/nndc_hdf5/cross_sections.xml`. **Only needed for tier-2.** For the
  analytic comparison the wall is near-void, so the actual nuclide data barely
  matters — but OpenMC still needs *a* valid `cross_sections.xml` to start.
- **Near-void / free-streaming** idioms in-tree: void cell `fill=None`
  (`test_cylindrical_mesh.py:124`) — neutrons stream with zero collisions, ideal
  for isolating source-direction physics; or ultra-low density
  `mat.set_density('g/cm3', 1e-5)` (`test_lost_particles.py:13`).
- **Wall-load tally**: `MeshSurfaceFilter` + `score='current'`
  (`test_mesh.py:863-866`); poloidal binning via `CylindricalMesh(phi_grid=…,
  z_grid=…, r_grid=…)` (`test_cylindrical_mesh.py:36-40`) or `SurfaceFilter` +
  `current` (`test_surface_flux.py:66-69`).
- **Fixed source run**: `settings.run_mode='fixed source'`, `batches`,
  `particles`, `inactive=0`. Read with `openmc.StatePoint(...).tallies[id].mean`.

## 0.6 ENVIRONMENT STATE — the big finding (the spec assumed this was done)

This checkout is **not built and has no Python environment**:
- ✗ No `build/` dir, no `libopenmc.so` anywhere, no `openmc` executable on PATH.
- ✗ `import openmc` fails — **numpy/scipy/h5py/… not installed**; the system
  `python3` (3.13.7) has **no `pip`** (`No module named pip`).
- ✗ No cross-section data; `OPENMC_CROSS_SECTIONS` unset.
- ✓ Toolchain present: `cmake 3.31.6`, `g++ 15.2` (C++17 OK).
- ✓ Network reachable (PyPI OK) — pip installs / data download are feasible.
- ✓ 76 GB free on `/`.
- ⚠ The only HDF5 on the box lives inside a sibling project:
  `…/OpenFUSIONToolkit/builds/hdf5-1_14_6/` (1.14.6), with its own
  `oft_venv` (py3.13, has h5py). No system HDF5.
- ⚠ The `.mcp.json` RAG server (`openmc_rag_search`/`openmc_rag_rebuild`) is
  **registered but not connected** this session — the tools don't surface, so I
  can't call them (and the "first-call rebuild widget" can't trigger). I used
  direct `Read`/`grep` instead (higher fidelity for known-file API work).

**Consequence (drives the plan):** Tier 1 (sampler-only) needs only a small
Python venv + a standalone C++ compile (no `libopenmc`). Tier 2 needs the full
OpenMC build (HDF5 wiring) + the `.so` linked against it + cross-section data.
Strict tier ordering ⇒ stand up the minimal env, get tier-1 green, *then* do the
heavy build for tier-2.

## 0.7 Physics re-derivation (confirming the spec's oracle, before any code)

Integrating Schwartz Eq. 2 over 4π (to be redone rigorously in sympy in tier-1):
`∫sin²θ dΩ = 8π/3`, `∫(¼+¾cos²θ) dΩ = 2π` ⇒
`σ_tot/σ₀ = a + (2/3)b + (1/3)c` ≡ **η** (total-rate factor). Corners confirmed:
non-pol(⅓,⅓,⅓)→**2/3**, A(1,0,0)→**1**, B(0,1,0)→**2/3**, C(0,0,1)→**1/3**. ✓
A ∝ sin²θ (perp to B̂); **B and C share** the ¼+¾cos²θ shape (parallel to B̂) and
differ only in rate — i.e. **C is NOT isotropic** (the spec's correction holds).

Angular sampling weights: `w_perp=(3/4)a`, `w_par=(2/3)b+(1/3)c`; stratified pick
`P_perp = a/η`, `P_par = ((2/3)b+(1/3)c)/η`.

**Re-derivation catch (the spec invited this):** for the sin²θ shape, with
`x=cosθ`, `p(x)∝(1−x²)`, inverse-CDF gives the depressed cubic
`x³ − 3x + (4u − 2) = 0` (verified: u=0→x=−1, u=1→x=+1). The spec wrote
`x³ − 3x + 2(1−2u) = 0`; that is the *same* sampler with `u→1−u` (also uniform),
so it is **equivalent, not wrong** — but I'll use the endpoint-verified form and
cross-check both against the rejection sampler (KS) in tier-1. For the ¼+¾cos²θ
shape: `p(x)=¼+¾x²` (already normalized on [−1,1]); inverse-CDF cubic
`x³ + x + (2 − 4u) = 0` (monotone, single real root; u=0→−1, u=1→+1). ✓

## 0.9 Tier-0 environment — VERIFIED working recipe (all steps passed)

Platform note: this box is **aarch64 (ARM64) Linux**; `$HOME = /home/tkiker.guest`
(working tree under `/Users/tkiker/...`). All of the below ran green.

```bash
# 1. System HDF5 (serial, has HL libs); sudo is passwordless here
sudo apt-get install -y libhdf5-dev            # -> HDF5 1.14.5, /usr/bin/h5cc

# 2. Dedicated Python venv (system python3.13 had no pip) + OpenMC (editable) + deps
python3 -m venv "$HOME/spf_venv"
"$HOME/spf_venv/bin/pip" install -e /Users/tkiker/Documents/GitHub/openmc sympy pytest
#   -> import openmc == 0.15.4.dev182+g608a1c338 ; numpy/scipy/sympy/h5py/pytest/mpl OK

# 3. Build OpenMC C++ and INSTALL into the venv prefix (build-tree config is NOT
#    consumable; install also lays down fmt/pugixml cmake configs find_package needs)
cd /Users/tkiker/Documents/GitHub/openmc && mkdir -p build && cd build
cmake .. -DOPENMC_USE_MPI=OFF -DOPENMC_USE_OPENMP=ON \
         -DCMAKE_BUILD_TYPE=RelWithDebInfo -DOPENMC_BUILD_TESTS=OFF \
         -DCMAKE_INSTALL_PREFIX="$HOME/spf_venv"
make -j"$(nproc)" && make install
#   -> $HOME/spf_venv/bin/openmc (0.15.4-dev182), lib/libopenmc.so,
#      lib/cmake/{OpenMC,fmt,pugixml}/*.cmake, include/openmc/*

# 4. Cross-section data (NNDC HDF5, 966 nuclides incl. H1)
#    already at $HOME/nndc_hdf5/cross_sections.xml  (set per-run:)
export OPENMC_CROSS_SECTIONS="$HOME/nndc_hdf5/cross_sections.xml"
#    (or in Python: openmc.config['cross_sections'] = '.../cross_sections.xml')

# 5. Build a custom source .so against the installed prefix  (VERIFIED)
#    CMakeLists.txt:  find_package(OpenMC REQUIRED)
#                     target_link_libraries(<tgt> OpenMC::libopenmc)
cmake .. -DCMAKE_PREFIX_PATH="$HOME/spf_venv" && make
#   -> .so exports `openmc_create_source` (nm -D shows 'T'); ldd resolves libopenmc.so

# 6. Run (venv bin on PATH so openmc.run() finds the exe; cap threads per AGENTS.md)
PATH="$HOME/spf_venv/bin:$PATH" OMP_NUM_THREADS=2 "$HOME/spf_venv/bin/python" model.py
```

**End-to-end smoke test (in-tree `parameterized_custom_source` ring + void sphere):**
Leakage = **1.0** (pure free-stream, zero collisions); surface current =
**1.0 ± 0.0 per source particle** (each neutron crosses once). The compiled-source
→ transport → tally pipeline is proven before writing any of our own code.

- **anarrima**: NOT on PyPI (`pip install anarrima` → "No matching distribution").
  Per the agreed fallback, Tier-2 uses an independent `scipy.quad` reimplementation
  of the Schwartz integral as the oracle; will attempt a git-repo install in Tier 2.

## 0.8 Deviations from the spec's assumptions (flag list)
1. **`ParticleType::neutron()` is a function call**, not an enum value `::neutron`.
2. **Environment is bare** — OpenMC unbuilt, no Python deps, no `pip`, no XS data.
   The spec assumed a working OpenMC; a real env bring-up is the first gate.
3. **RAG MCP tools not connected** this session (used direct reads instead).
4. **HDF5 only inside OpenFUSIONToolkit** — must point OpenMC's CMake at it
   (`HDF5_ROOT=…/builds/hdf5-1_14_6`) or install system HDF5 for the tier-2 build.
5. The sin²θ inverse-CDF cubic in the spec is the `u→1−u` mirror of the
   endpoint-verified form (equivalent). Will use the verified form + KS cross-check.
</content>
</invoke>
