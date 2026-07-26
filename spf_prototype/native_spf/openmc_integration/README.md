# Native SPF source — integrated into OpenMC core (canonical, version-controlled copy)

This directory is the **authoritative, git-tracked snapshot** of the spin-polarized-fusion (SPF)
neutron source *as integrated into the OpenMC core* — the version the Ginsburg transport pipeline
actually runs. It was previously living **only** as uncommitted edits in a working source tree on
the cluster (`/ginsburg/astro/users/tjk2147/spf_work/openmc_src`), which is exactly why it once
looked "lost." It is now preserved here so it can never disappear again.

## What this is
A parametric fusion source with polarized emission (Schwartz 2025, Eq. 1–2), added to OpenMC as:

- **New files** (mirrored under `src/`, `include/openmc/`, `openmc/` here at their real repo paths):
  - `src/stellarator_source.cpp` — the compiled `StellaratorSource` (fluxmap birth sampling ∝ √g
    over the true VMEC plasma; B̂ from the map; polarized angular sampler ported line-for-line from
    `spf_prototype/src/spf_sampler.hpp`).
  - `include/openmc/stellarator_source.h`
  - `openmc/stellarator_source.py` — the Python `openmc.StellaratorSource` wrapper (280 lines).
- **Registration hooks** in 5 existing core files, captured in `spf_core_integration.patch`
  (601-line diff): `CMakeLists.txt` (+1, builds the new .cpp), `include/openmc/source.h` (+43,
  SPF helpers + const members on the parametric source), `src/source.cpp` (+208, the Eq.2
  coefficient derivation + `TokamakSource` polarized path), `openmc/__init__.py` (+2, exports
  `StellaratorSource` and `spin_fractions_to_abc`), `openmc/source.py` (+181, `spin_fractions_to_abc`
  + `polarization`/`field_model`/`safety_factor` kwargs on `TokamakSource`).

## Base
Built on **PR #3999** (`parametric_tokamak_source`), base commit `23572ccc9` (see `BASE.txt`).
The hooks patch the parametric `TokamakSource` from that PR, so it applies on a pr3999-based tree,
**not** cleanly on `develop`. To fully re-home into this fork's main `src/`·`include/`·`openmc/`,
branch from PR #3999 (or cherry-pick it in) and then copy the 3 new files + `git apply` the patch.

## Reproduce / rebuild on Ginsburg (the working transport env)
```bash
SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src   # the tree with these files already applied
conda activate spf-stellarator
export PYTHONNOUSERSITE=1 PYTHONPATH=$SRC
export LD_PRELOAD=$SRC/build/lib/libopenmc.so PATH=$SRC/build/bin:$PATH
export OPENMC_CROSS_SECTIONS=/ginsburg/astro/users/tjk2147/spf_work/xs/endfb-viii.0-hdf5/cross_sections.xml
python -c "import openmc; print(hasattr(openmc,'StellaratorSource'))"   # -> True
```
Python usage (as every coil-dose / percoil run does it):
```python
src = openmc.StellaratorSource(fluxmap="qh_fluxmap", polarization=(a,b,c),
                               field_model="fluxmap",
                               energy=openmc.stats.Discrete([14.06e6],[1.0]), strength=1.0)
```
`fluxmap` is an `spf_fluxmap_v1` stem (`.meta`+`.bin` from `wout_to_fluxmap.py`). Note the two
openmc trees on the cluster — transport uses `spf_work/openmc_src`; adjoint uses
`~/src/GitHub/openmc-spf` (a bundle clone that does NOT have this source). Do not confuse them.

## Provenance
Snapshotted 2026-07-26 from `spf_work/openmc_src` (working-tree state: 3 untracked new files + 5
modified core files, all uncommitted). The prototype/standalone `CompiledSource` version remains at
`spf_prototype/src/polarized_fusion_source.cpp` (box-torus analytic validation); this integrated
version is what supersedes it for real 3-D DAGMC transport.
