# Paper direction (reconstructed 2026-07-07 after a lost chat transcript)

> **Provenance / honesty note.** The live pre-paper direction discussion was lost when the
> chat transcript was cleared (hit agents menu + esc). This note is a *reconstruction* from the
> surviving artifacts on disk — `LIT_REVIEW_NOVELTY.md`, `LIT_REVIEW_SINGULARITY_SUBTRACTION.md`,
> `PLAN_A_native_spf_tokamak.md`, `PLAN_B_stellarator_conformal_source.md`,
> `PHASE2_realdevice_and_scaleup.md`, `FINDINGS_geometry_peaking.md`, `OPTIMAL_POLARIZATION.md`.
> **tkiker: correct anything below that misremembers what we actually decided.**

## The two-paper split

### Paper 1 — Source & transport (native SPF NWL in OpenMC on real stellarators)
- **Content:** native C++ `TokamakSource`/`StellaratorSource` with local-B̂ polarized (a,b,c)
  emission → free-streaming transport on real QUASR/VMEC **conformal DAGMC** first walls →
  per-device polarization-optimal `a₂` (min first-wall peaking).
- **Novelty verdict (from lit review):** PARTIALLY NOVEL. **Bae et al. 2025 (Nucl. Fusion 65 086051)
  already put a polarized (a,b,c) DT source in OpenMC** (spherical tokamak, Python-precomputed
  point-source cloud). This is contribution-1's Achilles heel.
- **Defensible delta — lead with these, NOT "first in OpenMC":**
  native compiled on-the-fly sampling / differentiable / **local position-dependent B̂** /
  **stellarator-capable** (3D VMEC flux surfaces) / a₂-optimization as a new observable.
- **Open question — what is Paper 1's headline?**
  - (A) the **per-device optimal a₂** result (`OPTIMAL_POLARIZATION.md`: e.g. 886079 → +16.1%
    peaking reduction) — novel-in-combination, concrete engineering lever; OR
  - (B) the **geometry→peaking coherence law** (`FINDINGS_geometry_peaking.md`: PF ratio vs field
    coherence S_φ, Spearman r=−0.83 p=0.001 across 12 devices; QA vs QH separate p≈0.003–0.005) —
    arguably the more striking *physics* story, but free-streaming-only + class-level + n=12.
  - Leaning: (B) as the intellectual hook, (A) as the actionable result. **CONFIRM WITH tkiker.**

### Paper 2 — Analytic / numerical methods (off-axisymmetry NWL + singularity subtraction)
- **Content:** perturbative field-period series for polarized NWL on weakly-shaped (non-axisymmetric)
  stellarators + the closed-form singularity-subtraction quadrature for near-singular trig-polynomial
  line integrals (12 nodes beats 6144 Gauss–Legendre; JAX-differentiable, grads match FD to 1e-8).
- **Novelty verdict:** STRONGEST of the four claims. Schwartz 2025 (anarrima) is **axisymmetric only**;
  no field-period perturbation series for NWL found in the literature.
- **Closest prior art to differentiate against:** af Klinteberg–Barnett **singularity-swap** quadrature
  (BIT 2021/2024) — same subtraction/factoring family, but they Newton-iterate for the complex
  preimage; our φ*, D_min, D''(φ*) are **closed-form** because D(φ) is a trig polynomial. Transform
  methods (sinh/Telles/Duffy) are a different family and do not subsume us.

## Cross-cutting decisions to make (were likely live in the lost chat)
1. **Lead paper / submission order** — Paper 1 (results-forward, timely w/ Bae et al.) vs Paper 2
   (methods-forward, stronger novelty). ?
2. **Paper 1 headline = coherence law (B) or a₂-optimization (A)?** (see above) ?
3. **Scope of device set** — 5 VMEC devices (transport, `OPTIMAL_POLARIZATION`) vs 83-device
   free-streaming scaleup (`geometry_peaking_scaleup.json`, now enriched w/ io-shift metrics). Which
   carries which paper's claims. ?
4. **QI arm** — add 2–3 QI equilibria (DESC) to turn QA/QH into a 3-class landscape → strengthens
   Paper 1(B) generality, probes the low-S_φ end. In/out of first submission? ?
5. **Bae et al. framing paragraph** — needs to be written carefully and early; it is the reviewer's
   first objection. Draft it before methods.

## Status of the underlying results (all on disk, mostly committed)
- Sampler tiers 1–2 green; native TokamakSource 23 tests; StellaratorSource bit-exact + √g gates.
- Transport capstone (±40%/−21% vs anarrima) closed for both native sources.
- Headline transport: `OPTIMAL_POLARIZATION.md` (committed).
- Coherence law: `FINDINGS_geometry_peaking.md` (committed); 83-device scaleup enriched w/
  io_shift/redist metrics — **uncommitted** in `geometry_peaking_scaleup.json`.
- Lit reviews: both docs written, **uncommitted**.
