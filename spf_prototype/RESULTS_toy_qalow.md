# RESULTS — toy QA-low free-streaming: OpenMC ↔ analytic cross-check

The clean low-ε_eff QA validation for the **Monte-Carlo (OpenMC) SPF paper**: in the
free-streaming limit, OpenMC reproduces the analytic (anarrima `ripple`) neutron
wall load on a genuinely 3-D, field-period-rippled QA stellarator source — the
regime where the analytic series/exact/hybrid all agree, so the comparison is
meaningful. Extends the axisymmetric Tier-2b (toroidal) and Tier-8a (angled) MC↔
analytic validations to a **stellarator source**. Runs on the dev box (square-torus
CSG wall, no DAGMC).

## Setup (identical problem on both sides — `python/toy_qalow_config.py`)
- **Plasma/source:** anarrima tier3 toy QA-low spectrum (Nfp=2, `qa_spectrum(0.6)`),
  flux surfaces ρ=0.6,1.0 weighted ρ², 14 poloidal loops/surface; toy model field
  B̂(α(φ),β(φ)) with field-period ripple. **ε_eff ≈ 0.08** (series regime).
- **Wall:** square-cross-section torus (R_in 0.74, R_out 1.30, z=±0.34; R₀=1
  normalized) — axisymmetric, shared, CSG. All non-axisymmetry is in the SOURCE.
- **Comparison metric:** per-wall **directionality** A/iso and B/iso (the
  polarization steering), which cancels source-position/normalization and isolates
  the emission about the local field.

## Result — agreement to MC statistics
OpenMC free-streaming (1×10⁶ pre-sampled births/mode via the bit-parity mirror
sampler, vectorized presampler) vs the analytic quadrature (384 nodes, 40 wall
patches/wall). Per-wall directionality is the **ratio of total loads** (ΣA/Σiso) on
both sides (apples-to-apples), and `corrA` is the per-bin Pearson correlation of the
poloidal A/iso profile:

| wall | A/iso MC | A/iso ana | rel | B/iso MC | B/iso ana | rel | corrA |
|---|---|---|---|---|---|---|---|
| inboard | 1.207 | 1.229 | 1.8% | 0.782 | 0.771 | 1.4% | 0.85 |
| outboard | 0.823 | 0.829 | 0.7% | 1.182 | 1.171 | 1.0% | 0.46 |
| floor | 1.065 | 1.069 | 0.4% | 0.937 | 0.931 | 0.7% | 0.96 |
| ceiling | 1.062 | 1.066 | 0.4% | 0.937 | 0.934 | 0.3% | 0.98 |

**Worst-wall discrepancy 1.8% (all four walls)**, poloidal shapes tracked. (Outboard's
lower corr just reflects that its analytic profile is nearly flat — little shape to
track.) Consistent with MC statistics.

## Analytic engine (rewritten; provably ≡ anarrima)
The analytic is now a fast pure-numpy free-streaming quadrature evaluating the loop
geometry and field DIRECTLY at the Gauss-Legendre nodes (identical to what the OpenMC
source samples), with correct visibility: the TRUE shaped-loop N>0 sub-arcs
intersected with the inner-cylinder occlusion (the central hole blocks far-side
source). It is gated at startup against `anarrima.free_streaming_quadrature` and
matches to **3e-14** here (1e-12 on QUASR). This replaced an earlier crude `front_arc`
(mean-radius, N0>0, cap=2.5) that biased the small inboard arc and double-counted the
outboard far side — errors that cancel in anarrima's own series-vs-quadrature tests
(same arc both sides) but not against ray-tracing MC. The OpenMC direction sampler was
separately verified to reproduce the normalized kernels exactly (A/iso=(3/2)sin²θ to
0.3%/bin, mean 1.0000).
The physics is the expected SPF steering: A (∝sin²θ) enhances the inboard / floor and
suppresses the outboard load; B/C (∝¼+¾cos²θ) is the mirror image.

Figures (`plot_toy_qalow.py` -> `figs/`): `toy_qalow_directionality` (A/iso, B/iso
poloidal profiles, analytic line vs OpenMC markers, all four walls),
`toy_qalow_agreement` (per-bin MC-vs-analytic scatter on y=x), `toy_qalow_geometry`
(rotating plasma cross-sections in the square wall).

## Convention note (state this when the two methods are shown together)
anarrima's quadrature uses **unnormalized** kernels (∫sin²θ dΩ = 8π/3,
∫(¼+¾cos²θ)dΩ = 2π); the OpenMC source samples emission PDFs **normalized to unit
total**. The exact reconciling factors are A ×3/2, B ×2 (= 4π/∫K dΩ). This is now a
first-class option in anarrima: `free_streaming_quadrature(..., normalize=True)`
divides each kernel by its solid-angle mean (`_KERNEL_MEAN`), so the analytic side is
directly comparable to a unit-emission MC source (verified A→×3/2, B→×2 exactly). The
toy analytic script uses `normalize=True`. Same rate(η)-vs-shape separation used
throughout the prototype. (This kwarg is a small, default-off addition to the anarrima
working tree — include it in the anarrima PR/JOSS so others can reproduce the MC
cross-check; it does not change any existing anarrima result.)

## Reproduce
```
# analytic (anarrima venv) -> /tmp/toy_qalow_analytic.npz
$HOME/ana_venv/bin/python spf_prototype/python/toy_qalow_analytic.py
# OpenMC (spf venv) -> the agreement table
OMP_NUM_THREADS=4 $HOME/spf_venv/bin/python spf_prototype/python/toy_qalow_openmc.py
```

## Scope / two-paper split (see also ana_precise_qa.py)
This is the **low-ε_eff QA** overlap where analytic and MC cross-validate (the MC
paper's credibility figure; also the analytic paper's home turf). The **high-ε_eff**
real equilibrium **precise_QA (ε_eff ≈ 2.3)** is MC-only: there the analytic series
diverges AND even the fixed-node exact quadrature breaks down (≈34% exact-vs-hybrid
disagreement and unphysical negatives, from the near-singular integrand / visibility
breakdown when the strongly-shaped loop nearly grazes the wall — see
`python/ana_precise_qa.py`). That is the MC-only regime the MC paper owns.
