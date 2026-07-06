# DATA_NEEDED — real Helios inputs to replace the placeholders

Tier 8 / Part B runs end-to-end on **public, literature-grounded placeholders**
(Helios overview arXiv:2512.08027; ARIES-CS radial build, El-Guebaly UWFDM-1336).
Every device-specific value is tagged `# INJECT(helios):` in the code at the seam
where a real Thea-supplied number drops in. This file lists those inputs, ranked.
Nothing here is a Helios physics claim until these are filled.

> **If any input is proprietary** (not from the public papers), keep it in a
> gitignored `spf_prototype/data/` path and out of commits.

## 1. Equilibrium / 3D field — the single biggest item
The source's spin-steering axis is the **local field direction B̂(x)**, and the
birth region is the **plasma boundary / flux surfaces**. Both come from the
equilibrium. In order of usefulness:
- **Best: a full equilibrium** — VMEC `wout_*.nc` or a DESC `.h5`. Gives the 3D
  B-field, the LCFS, nested flux surfaces, and profiles in one file. Drops into
  `python/make_standin_field.py`'s producer seam (writes the same `.meta`+`.bin`
  field map that `src/spf_fieldmap.hpp` / `python/fieldmap.py` already read), and
  into a real flux-surface source. **Code seam:** `make_standin_field.write_map`
  (replace the analytic field with equilibrium evaluation).
- **Good: boundary Fourier fits at φ=0, π/4, π/2, 3π/4** + plasma params. Fixes the
  3D source REGION shape and the wall shape (offset). Needed convention: VMEC
  `RBC(m,n)/ZBS(m,n)` + `nfp` + stellarator-symmetry flag, OR per-cut
  `R(θ)=Σ aₘcosmθ+bₘsinmθ`, `Z(θ)=...` with the m-range and units. **Note:** the
  boundary alone does NOT give B̂(x); the field map stays the stand-in (or is
  approximated from the cuts) until a field is supplied. **Code seam:** the
  `shape=plasma` spatial sampler in `src/polarized_fusion_source.cpp` (replace the
  circular `(1-ρ²)` disk with the real boundary) + `make_standin_field`.

## 2. Exact radial build at the ≥1.2 m standoff
Per-layer cm thicknesses + materials within the standoff: FW armor/coolant
fractions, blanket coolant channels, the full shield stack (WC / B₄C / 316L /
borated water / HDPE), vacuum vessel, coil case, winding-pack location.
**Placeholder now:** W 0.1 / steel 1.9 / Be 5 / FLiBe 50 / WC-shield 40 / coil 15 cm
(`python/stellarator_model.py:LAYERS`). **Code seam:** `stellarator_model.LAYERS`
+ `make_materials`.

## 3. Breeder spec
PbLi vs FLiBe vs ceramic; exact ⁶Li enrichment (Helios paper states ~65% — confirm);
coolant, structure fractions, operating temperature (sets the data temperature).
**Placeholder:** FLiBe (Li₂BeF₄), 65% ⁶Li, 294 K data (`reactor_model.make_materials`,
`stellarator_model.make_materials`). **Code seam:** `make_materials`.

## 4. Coil / conductor spec + radiation limits
HTS (ReBCO) tape-stack composition, winding-pack geometry, and the ReBCO
**radiation limits** Thea uses (fast-fluence threshold, dose, dpa proxy). The
ARIES-CS LTS limits (≤1e19 n/cm² >0.1 MeV, ≤0.006 dpa Cu, ≤1e11 rad insulator) are
placeholders. **Placeholder:** a steel coil block (`stellarator_model.make_materials`
`coil`). **Code seam:** `coil` material + the `coil_fast`/`coil_damage` tallies.

## 5. Source rate / fusion power
Confirmed fusion power + neutron source rate (Helios public ~958 MW) to convert
per-source-neutron tallies → absolute MW/m², W/cm³, dpa/FPY, fluence/plant-life.
Currently all results are **per source neutron**. **Code seam:** postprocessing
(multiply per-neutron tallies by the source rate).

## 6. Plasma profiles (reaction-rate weighting)
Density/temperature profiles n(ρ), T(ρ) for a physical reaction-rate weight on the
source (currently a `(1-ρ²)` stand-in). With Tᵢ(ρ) the Ballabio energy-broadening
hook (`polarized_fusion_source.cpp`, `INJECT(helios)`) can also be turned on.
**Code seam:** the `shape=plasma` weight in `polarized_fusion_source.cpp`.

## 7. Configuration scalars (cite to replace literature defaults)
field periods (placeholder **2**; ARIES-CS uses 3), aspect (**4.5**), R0 (**8 m**),
a (**1.8 m**), B0 (**6 T**), NWL avg/peak (ARIES-CS proxy 2.6/5.3 MW/m²), plant
lifetime/duty (40 FPY / 88%). **Code seam:** `stellarator_model` module constants
+ `make_standin_field.DEFAULTS`.

---
**Status:** the collaborator has offered (1) boundary Fourier fits at φ=0,π/4,π/2,3π/4
and (6/7) plasma params from the Helios papers — these fill items 1 (boundary),
6, and 7. A full VMEC/DESC equilibrium would additionally fill item 1's field map.
