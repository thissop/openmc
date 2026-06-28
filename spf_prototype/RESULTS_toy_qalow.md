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
OpenMC free-streaming (3×10⁵ pre-sampled births/mode via the bit-parity mirror
sampler) vs anarrima exact quadrature (384 nodes):

| wall | A/iso MC | A/iso ana | rel | B/iso MC | B/iso ana | rel |
|---|---|---|---|---|---|---|
| inboard | 1.221 | 1.236 | 1.2% | 0.793 | 0.764 | 3.8% |
| outboard | 0.818 | 0.839 | 2.5% | 1.174 | 1.161 | 1.1% |
| floor | 1.067 | 1.060 | 0.7% | 0.944 | 0.940 | 0.5% |
| ceiling | 1.059 | 1.057 | 0.2% | 0.931 | 0.943 | 1.3% |

**Worst-wall discrepancy 3.8%** (inboard B/iso — the lowest-current wall, so the
noisiest); all others 0.2–2.5%. Consistent with MC statistics. The physics is the
expected SPF steering: A (∝sin²θ) enhances the inboard / floor and suppresses the
outboard load; B/C (∝¼+¾cos²θ) is the mirror image.

## Convention note (state this when the two methods are shown together)
anarrima's quadrature uses **unnormalized** kernels (∫sin²θ dΩ = 8π/3,
∫(¼+¾cos²θ)dΩ = 2π); the OpenMC source samples emission PDFs **normalized to unit
total**. The exact reconciling factors are A ×3/2, B ×2 (= 4π/∫K dΩ), applied to the
analytic side in `toy_qalow_analytic.py`. This is the same rate(η)-vs-shape
separation used throughout the prototype.

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
