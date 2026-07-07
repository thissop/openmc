# Native SPF for `openmc::TokamakSource` — staging + apply guide

This directory is a **staging area** that adds spin-polarized-fusion (SPF)
anisotropic neutron emission to Ethan Peterson's native `TokamakSource`
(PR #3999). It is implemented as edits on top of the three extracted PR files so
a clean patch can be generated and `git apply`-ed onto an `openmc_src/pr3999`
checkout on the cluster.

Design source of truth:
`spf_prototype/docs/notes/planning/PLAN_A_native_spf_tokamak.md`.
Validated math ported from `spf_prototype/src/spf_sampler.hpp`,
`spf_field.hpp`, `polarized_fusion_source.cpp`.

## What is here

| File | Role |
|------|------|
| `source.h`        | edited — `include/openmc/source.h` with SPF members + helper decls |
| `source.cpp`      | edited — `src/source.cpp` with the SPF branch, samplers, ctor wiring |
| `source.py`       | edited — `openmc/source.py` with the polarization API + XML round-trip |
| `*.orig`          | pristine PR baselines (byte-for-byte the extracted PR files) |
| `test_spf_tokamak.py` | unit tests (sampler-math oracle + Python-API layer) |
| `README_apply.md` | this file |

The `.orig` files are the exact PR sources; the non-`.orig` files are those same
files **plus** the SPF changes. That pairing is what lets you produce a minimal
patch (below).

## Backward compatibility (the whole point)

`polarization=None` (the default) writes **no** new XML, so C++ sees
`polarized_ = false` and takes the **unchanged** isotropic draw
(`site.u = angle_->sample(seed).first;`). Existing `test_source_tokamak.py`
passes unmodified; existing input files are byte-identical.

## Generate the patch

```bash
cd spf_prototype/native_spf

# Per-file unified diffs (labels rewritten to the real repo paths).
git diff --no-index source.h.orig   source.h   \
  | sed -e 's#source.h.orig#a/include/openmc/source.h#' \
        -e 's#b/source.h#b/include/openmc/source.h#'   > /tmp/spf_source_h.patch
git diff --no-index source.cpp.orig source.cpp \
  | sed -e 's#source.cpp.orig#a/src/source.cpp#' \
        -e 's#b/source.cpp#b/src/source.cpp#'         > /tmp/spf_source_cpp.patch
git diff --no-index source.py.orig  source.py  \
  | sed -e 's#source.py.orig#a/openmc/source.py#' \
        -e 's#b/source.py#b/openmc/source.py#'        > /tmp/spf_source_py.patch
```

(`git diff --no-index` works even outside a repo and ignores `.gitignore`.)

## Apply onto the PR checkout

```bash
cd openmc_src/pr3999            # a checkout that already contains PR #3999
git checkout -b spf-tokamak     # never commit onto main/develop

git apply --check /tmp/spf_source_h.patch /tmp/spf_source_cpp.patch /tmp/spf_source_py.patch
git apply         /tmp/spf_source_h.patch /tmp/spf_source_cpp.patch /tmp/spf_source_py.patch
```

If `git apply` complains about context drift (the PR moved), fall back to
`git apply --3way` or copy the three edited files in wholesale — they are the
complete PR files, not fragments:

```bash
cp source.h   ../pr3999/include/openmc/source.h
cp source.cpp ../pr3999/src/source.cpp
cp source.py  ../pr3999/openmc/source.py
```

Also export the helper into the package namespace so `openmc.spin_fractions_to_abc`
resolves (one line in `openmc/__init__.py`, next to the other `source` imports):

```python
from openmc.source import spin_fractions_to_abc
```

## Build

```bash
cd openmc_src/pr3999/build      # or wherever you configure
cmake --build . -j
pip install -e ..               # rebuild the Python extension bindings
```

The SPF code adds **no** new dependencies and **no** new translation units — it
is confined to `source.{h,cpp,py}`. It uses only `<cmath>`, `<algorithm>`,
`<cassert>` (all already included by `src/source.cpp`), `openmc::prn`, and
`openmc::Direction`'s existing operators.

## Test

```bash
# Sampler-math oracle (no build needed; validates the ported physics):
pytest -q spf_prototype/native_spf/test_spf_tokamak.py

# After building + installing, the Python-API layer un-skips automatically.
# Also run the PR's own regression to confirm zero behavior change:
pytest -q openmc_src/pr3999/tests/unit_tests/test_source_tokamak.py
```

The unit test has two layers:

1. A self-contained NumPy mirror of the exact `source.cpp` SPF math
   (`field_direction` + `sample_polarized_direction`). It checks the physics
   oracles: `(1/3,1/3,1/3)` → flat cosθ (isotropic), pure-A → concentrates ⊥ b̂
   (`<cos²θ>=1/5`), pure-B/C → concentrates ∥ b̂ (`<cos²θ>=7/15`), B/C share
   shape, Gram-Schmidt rotation preserves `u·b̂`, off-axis isotropy, degenerate
   `b̂=±ẑ`, toroidal `b̂=φ̂`, and pitched → toroidal as `q→∞`. Runs today.
2. Python-API tests (`polarization` setter, spin-fraction Eq. 1 map, XML
   round-trip, `pitched` requires `safety_factor`). Skipped until `openmc` with
   the patched `TokamakSource` is importable.

## Remaining build-time gate (per PLAN A §4.1)

Before trusting the native port in production, run the **bit-for-bit parity
test**: the native `sample_polarized_direction` must reproduce
`spf::sample_global_direction` (`spf_prototype/src/spf_sampler.hpp`) on an
identical seed sequence (same RNG draw order: selector → shape-u → φ). The math
is transcribed line-for-line, but the parity test is what proves the
transcription. Then the headline axisymmetric-NWL comparison vs anarrima
(`field_model='toroidal'`, the `B̂=φ̂` case) reproduces Schwartz's
±43% inboard / ∓22% outboard midplane numbers.

## What the patch changes (summary)

**`include/openmc/source.h`** — in `class TokamakSource`:
- private helper decls: `field_direction`, `sample_polarized_direction`,
  `interp_q`;
- private const data: `polarized_`, `w_perp_`, `w_par_`, `p_perp_`,
  `field_model_`, `field_sign_`, `q_r_over_a_`, `q_values_` (all set in ctor).

**`src/source.cpp`**:
- `#include <cassert>`;
- anon-namespace inverse-CDF helpers `spf_costheta_perp` / `spf_costheta_par`
  (ports of `spf_sampler.hpp`);
- constructor block reading `polarization` / `field_model` / `q_*` /
  `field_sign` and deriving `w_perp/w_par/η/P_perp` (no pasted magic constants);
- `interp_q`, `field_direction` (toroidal + pitched), and
  `sample_polarized_direction` (stratified P2 sampler + Gram-Schmidt rotation);
- the one-line hook: the isotropic draw is now branched on `polarized_`.

**`openmc/source.py`**:
- module helper `spin_fractions_to_abc` (Schwartz Eq. 1);
- `TokamakSource.__init__` kwargs `polarization` / `field_model` /
  `safety_factor` / `field_sign` + a cross-field check (pitched ⇒ safety_factor);
- setters for each (polarization accepts `(a,b,c)` or a spin-fraction dict,
  warns + renormalizes if `Σ≠1`);
- `populate_xml_element` / `from_xml_element` round-trip (only emitted when
  polarized).

## Scope / limitations (see PLAN A §"out of scope")

Birth distribution only (no depolarization in transport), energy decoupled,
volumetric/weighted-ring plasma source, toroidal or `q_cyl`-pitched `b̂` (not a
solved equilibrium), no alpha channel, convex free-stream NWL validation.
