# Native `StellaratorSource` — conformal spin-polarized stellarator source

A native OpenMC-core `Source` subclass (patterned on Ethan Peterson's
`TokamakSource`, PR #3999) that samples birth **positions** from a real
DESC/VMEC-equilibrium flux map via a rejection-free CDF cascade, and birth
**directions** from the Schwartz P2 spin-polarized-fusion (SPF) mixture about the
local field b̂ read from the *same* grid. It is a line-for-line C++ port of the
validated Python reference `ConformalStellaratorSampler`
(`spf_prototype/python/stellarator_source.py`).

Files in this directory:

| File | Role |
|------|------|
| `stellarator_source.h`            | class decl → goes to `include/openmc/stellarator_source.h` |
| `stellarator_source.cpp`          | impl → goes to `src/stellarator_source.cpp` |
| `stellarator_source.py`           | Python API (`StellaratorSource(SourceBase)`) → merge into `openmc/source.py` |
| `test_stellarator_source_cpp.py`  | NumPy-mirror unit test (cascade math vs the reference) |
| `README_stellarator.md`           | this file |

The SPF direction sampler (`sample_polarized_direction`, the `spf_costheta_*`
inverse-CDF helpers, the `(a,b,c) → w_perp/w_par/P_perp` derivation) is **reused
verbatim** from the native SPF `TokamakSource` work (`native_spf/source.cpp`); it
is not re-derived here.

---

## 1. The portable `spf_fluxmap_v1` binary format

Mirrors the style of `src/spf_fieldmap.hpp`'s `spf_fieldmap_v1`: a small text
`.meta` header plus a raw little-endian `float64` `.bin` payload (HDF5-free, so
the C++ reader is a plain `std::ifstream`; producer-agnostic).

**Writer:** `python/quasr_fluxmap.py:write_fluxmap_bin(...)`, now called inside
`build()` alongside the existing `np.savez` — the DESC solve emits **both** the
`.npz` and `quasr<ID>_fluxmap.{meta,bin}`.

`<stem>.meta` (text, `key value` per line):

```
magic spf_fluxmap_v1
nR <nR>
ntheta <nθ>
nzeta <nζ>
nfp <nfp>
sign_sqrtg <±1.0>
```

`<stem>.bin` (little-endian `f8`, C-order, concatenated in **this exact order**):

```
rho[nR], theta[nθ], zeta[nζ],
sqrtg, R, Z, phi, BR, Bphi, BZ            (each nR*nθ*nζ)
```

All 3-D arrays are indexed `[i_rho, j_theta, k_zeta]` (C-order) — the **same**
layout as `ConformalStellaratorSampler`. `theta`/`zeta` span the full torus
`[0, 2π)` (endpoint excluded); `rho` runs `(drho .. 1]` (the axis `rho=0` is
excluded because `sqrt(g) → 0` there). `sqrtg` is stored as the non-negative
magnitude `|sqrt(g)|`; `sign_sqrtg` records the sign of the raw value for
provenance. `nfp` is informational — the sampler treats the grid as the full
torus with `2π` periods (matching the reference).

---

## 2. Sampling algorithm & conventions (must match the reference line-for-line)

`sample(uint64_t* seed)` draws, **in this RNG order** (`openmc::prn` only):

1. **Cascade cell selection** — `searchsorted(side="right") − 1`, clamped
   (`std::upper_bound`), on:
   - `cdf_rho` (marginal `p(rho)`), then
   - `cdf_zeta[i]` (conditional `p(zeta|rho)`), then
   - `cdf_theta[i,:,k]` (conditional `p(theta|rho,zeta)`).
   Draw order: **rho → zeta → theta**. Returned cell is `(i_rho, j_theta,
   k_zeta)`.
2. **In-cell jitter**, uniform, draw order **rho → theta → zeta**
   (`x = grid[cell] + (prn − 0.5)·d`).
3. **Trilinear interpolation** of `R, Z, phi` → cylindrical → Cartesian. Term
   order is **theta (inner) → zeta → rho (outer)**. `rho` bracket uses
   `lower_bound` (`side="left"`) and is **clamped** (non-periodic; jitter may
   extrapolate a half-cell past the axis/edge, exactly as the reference does —
   the bias → 0 under grid refinement). `theta`/`zeta` brackets use the periodic
   `_wrap_index` with period `2π` and a non-negative modulo (correct for jitter
   just below `grid[0]`).
4. **b̂**: `field_model='fluxmap'` (default) interpolates `(BR,Bphi,BZ)`, rotates
   `rpz→Cartesian` (`bx = br·cosφ − bp·sinφ`, `by = br·sinφ + bp·cosφ`) using the
   **interpolated** `phi`, then normalizes. `field_model='toroidal'` overrides
   with pure `φ̂ = (−sinφ, cosφ, 0)`.
5. **Direction**: isotropic if unpolarized; else the SPF P2 mixture about b̂
   (RNG order selector → shape-u → φ, matching `spf_sampler.hpp`).
6. **Energy**, then **time**.

> Period-wrap note: `phi` is trilinearly interpolated with the same periodic
> `zeta` wrap as every other grid field (inherited from the reference). Across
> the last `zeta` cell (`k=nζ−1 → 0`) `phi` seams from `~2π` back to `0`; this is
> the reference's convention and is reproduced faithfully rather than "fixed".

**Thread safety / invariants.** All state is set in the constructor and `const`
at sample time (no mutable members). Hard-`assert`ed invariants: `sqrt(g) ≥ 0`,
`cdf_rho` monotone and ends at 1, `P_perp ∈ [0,1]`. User-input errors
(`polarization < 0`, missing `fluxmap`, bad `field_model`) are `fatal_error`,
and a non-unit `a+b+c` is renormalized (with a Python-side warning).

---

## 3. Register the class in an OpenMC build

Copy the two C++ files into the tree and add three small hooks:

```bash
cp stellarator_source.h   <openmc>/include/openmc/stellarator_source.h
cp stellarator_source.cpp <openmc>/src/stellarator_source.cpp
```

**`CMakeLists.txt`** — add the new translation unit next to `src/source.cpp`
(around line 426):

```cmake
  src/source.cpp
  src/stellarator_source.cpp      # <-- add
```

**`src/source.cpp`** — include the header and add a dispatch branch in
`Source::create` (right after the `"tokamak"` branch):

```cpp
#include "openmc/stellarator_source.h"
// ...
    } else if (source_type == "stellarator") {
      return make_unique<StellaratorSource>(node);
    }
```

**`openmc/source.py`** — merge the `StellaratorSource` class and the
`spin_fractions_to_abc` helper from `stellarator_source.py` into `openmc/source.py`
(or keep it a submodule and import it there), and add a dispatch branch in
`SourceBase.from_xml_element`:

```python
            elif source_type == 'stellarator':
                return StellaratorSource.from_xml_element(elem)
```

Also export the helper into the package namespace (next to the other `source`
imports in `openmc/__init__.py`):

```python
from openmc.source import spin_fractions_to_abc
```

Then build:

```bash
cd <openmc>/build && cmake --build . -j
pip install -e ..
```

The class adds **no** new dependencies: only `<fstream>`, `<cmath>`,
`<algorithm>`, `<cassert>`, `openmc::prn`, `openmc::Direction`, and the existing
XML/distribution helpers.

## 3b. Apply onto an `openmc_src/pr3999` checkout

`StellaratorSource` reuses the SPF direction sampler that the PLAN-A patch adds to
`TokamakSource`, but is otherwise a **self-contained new translation unit** — it
does not depend on that patch being applied first (it carries its own copy of the
`spf_costheta_*` helpers). Order of operations on the cluster:

```bash
cd openmc_src/pr3999
git checkout -b spf-stellarator        # never commit onto main/develop
cp <this>/stellarator_source.h   include/openmc/
cp <this>/stellarator_source.cpp src/
# apply the three hooks from §3 (CMake line, Source::create branch, Python dispatch)
cd build && cmake --build . -j && pip install -e ..
```

---

## 4. Validation gates

**Producer (already gated, `python/quasr_fluxmap.py`).** The `√g` volume
identity and equilibrium quality are asserted at write time: `V1` own-quadrature
volume vs DESC `V` (`<2e-3`), `V1b` independent LCFS divergence-theorem volume
(`<1e-2`), `V2` single-signed `√g` (nested surfaces), `V3` clean tensor-product
node set, `V4` LCFS reconstructs the QUASR boundary, plus `V0` force-balance
convergence. A fluxmap that fails these never reaches the sampler.

**Sampler (this directory).** `test_stellarator_source_cpp.py` is a NumPy mirror
of the C++ cascade math checked, on a small synthetic circular-torus fluxmap
(known analytic `√g`), against the validated reference:

- `spf_fluxmap_v1` `.meta/.bin` round-trips bit-exactly;
- the C++ (sequential-loop) CDF construction reproduces the reference's CDFs to
  `~1e-12`;
- **same seeds → same cells**: identical uniform draws give identical `(i,j,k)`;
- interpolated positions and grid b̂ match the reference to `~1e-12` (validates
  the interp term order and the `theta/zeta` period wrap);
- birth weight is exactly 1 (rejection-free) and every birth is inside the
  plasma (`rho ≤ 1 + jitter`);
- on the circular torus (`BR=BZ=0`), the grid b̂ equals `φ̂`, so
  `field_model='fluxmap'` and `'toroidal'` agree.

Run: `python native_spf/test_stellarator_source_cpp.py` (pure NumPy; no build) or
`pytest -q native_spf/test_stellarator_source_cpp.py`.

**Build-time gates (after building on the cluster):**

1. **bit-for-bit prn parity** — dump sampled positions from a tiny standalone
   driver and confirm the C++ `openmc::prn` stream reproduces the NumPy mirror on
   an identical seed sequence (same draw order §2). The math is transcribed
   line-for-line; this proves the transcription (analogous to the SPF
   `spf_sampler.hpp` parity gate in `README_apply.md`).
2. **axisymmetric limit vs `TokamakSource`** — feed an axisymmetric (`nfp=1`,
   circular/Miller) fluxmap with `field_model='toroidal'` and confirm the
   free-stream neutron-wall-loading pattern matches the native SPF
   `TokamakSource` (`field_model='toroidal'`, the `B̂=φ̂` case), reproducing
   Schwartz's ±43% inboard / ∓22% outboard midplane numbers. This ties the
   conformal sampler back to the analytic/anarrima oracle.

---

## 5. Scope / limitations

Birth distribution only (no depolarization in transport); energy decoupled
(single distribution for all radii in this prototype); volumetric source drawn
from one solved equilibrium fluxmap; b̂ from the gridded field or the pure
toroidal override; no alpha channel; convex free-stream NWL validation. The
`(a,b,c)` → direction/rate split and the SPF physics match the native SPF
`TokamakSource` exactly.
