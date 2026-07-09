# PAPER_DIRECTION — the single capstone paper

*Authoritative plan for the one paper that closes out this summer's SPF-neutronics work.
Written 2026-07-08. Supersedes the exploratory `docs/notes/planning/PAPER_DIRECTION.md`.*

## Thesis (one sentence)

> A native, equilibrium-rigorous, **spin-polarized** fusion neutron source for OpenMC — and a
> demonstration of **where and why Monte Carlo transport is necessary** for polarized stellarator
> wall loading (scattering dilutes the geometric steering ~3×; the filamentary/analytic method
> breaks in the strongly-shaped ε_eff ≳ 2 regime) — with an analytic companion (anarrima) that maps
> the boundary of its own validity.

This is a **methods + honest-findings** paper (same genre as Peterson's TokamakSource paper), not a
discovery paper. It does not need a headline physics surprise; it delivers a validated capability and
characterizes it honestly.

## Novelty framing (what to claim, what NOT to claim)

| Claim | Status | Prior art to cite & differentiate |
|---|---|---|
| First **native, on-the-fly, general** polarized fusion source in a production MC code | defensible | Bae 2025 = precomputed static point cloud, one fixed tokamak |
| First **polarized stellarator** source (b̂ from a real 3-D equilibrium) | first-of-kind | nobody |
| Equilibrium-rigorous **StellaratorSource** (√g volumetric, gated) in OpenMC | first *in OpenMC* | Lyytinen 2024 (Serpent2, ported-MCNP isotropic source, parametric *geometry*); ParaStell; Wu/W7-X |
| Free-streaming **overestimates steering ~3×** (scattering correction) | our result | Bae reports the MC benefit but did not isolate the free-streaming vs transport gap |
| **Where MC is necessary**: filamentary-analytic error ∝ shaping², quadrature breakdown at ε_eff≳2 | our result | — |

**Do NOT write:** "first polarized source in OpenMC" (Bae), "first MC stellarator source" (Lyytinen/ParaStell).
**DO write:** "first *native/general/on-the-fly* polarized source in a production MC code; first *polarized
stellarator* source; first equilibrium-rigorous stellarator source *in OpenMC*."

Peterson precedent: a native/validated *increment* is publishable methods-paper content. We clear that
bar more easily than Peterson (no predecessor stellarator source in OpenMC; we add polarization).

## The fidelity ladder (the coherence spine — state this explicitly)

Each model is used at exactly the fidelity needed to isolate one effect:

1. **Analytic square-torus (Schwartz)** — verification of the sampler & NWL (Tier 1–2).
2. **Free-streaming, real QUASR equilibrium (59509)** — geometric steering + the analytic-breakdown
   boundary (analytic↔MC agreement; bare-cylinder ∝s²; precise_QA quadrature breakdown).
3. **Layered ARC-class blanket (box-torus)** — the scattering correction (~3× dilution) and the
   rate-vs-steering trade-off (B is the sweet spot).
4. **3-D-field-driven Helios-class radial build** — end-to-end demonstration; directional efficiency
   η = 0.89 under a real pitched field.
5. **anarrima analytic companion** — the differentiable polarized analytic method whose validity
   boundary (2) rigorously characterizes.

## Structure (sections)

1. **Introduction** — SPF motivation (reactivity + component protection, Lyytinen's aim (2)); prior art
   up front (Schwartz, Bae, Lyytinen, Peterson); the gap.
2. **Polarized emission physics & the native sampler** — (a,b,c) modes, w(θ_B)∝1+a₂P₂, angular/rate
   separation, stratified exact sampler, thread-safe design. *(co-headline)*
3. **StellaratorSource: polarized emission from a verified 3-D equilibrium** — VMEC √g/geometry/b̂,
   rigor gates, TokamakSource-pattern extension. *(the vehicle)*
4. **Verification** — Tier 1 sampler stats; Tier 2 NWL-vs-Schwartz; the exact a₂-linearity identity.
5. **Where Monte Carlo is necessary** *(the thesis section)* — free-streaming QUASR 59509 analytic↔MC;
   bare-cylinder ∝s² (3× worse inboard); precise_QA quadrature breakdown at ε_eff≈2.3; **Tier 5 ~3×
   scattering dilution**; **Tier 7 rate-vs-steering (B sweet spot)**; **Tier 8b η=0.89**.
6. **Analytic companion (anarrima)** — singularity subtraction; differentiable polarized quadrature;
   ε_eff two-threshold; gradient closure (with XLA caveat); the gc₂ correction; validation triad;
   cross-validation of the weak-global-lever finding; brief field-period-perturbation aside (→ ripple).
7. **Discussion / Conclusion** — SPF as local actuator, not global flattener, agreed by two methods;
   native SPF now in OpenMC for both device classes.

**Appendices:** coherence-C predictor + honest n=83 collapse; nematic/C-term + exp(−N a/R₀) reason;
Carlson/Legendre engine; QUASR mining tooling; **Bae (a,b,c) symbol-collision footnote**.

## IN / OUT (locked)

- **IN:** §2–§6 above; the five ladder models; the a₂-linearity identity; the three "so-what" results
  (3× dilution, B sweet spot, η=0.89); the MC-necessity proof (bare-cylinder + precise_QA).
- **OUT (agreed):** Tier 6 TBR (one-line cross-validation of Bae only); the η=η_source·A(τ)
  factorization (A(τ) arm never run); the angular-adjoint memo (Discussion outlook only, unmerged
  dependency); the unrealized "precise_QA field-map" third-party cross-check (do not claim).

## Honest gaps (what bounds the claims — state in Limitations)

- Tiers 5–7 use **box-torus/layered-blanket** geometry, not the conformal QUASR wall — fidelity ladder
  must be stated or it reads as a grab-bag.
- **Polarized-at-grazing** (anarrima subtraction) is by-construction, **untested**.
- **XLA jit(grad)** through the anarrima engine hangs; eager grad works (~9s).
- anarrima **OpenMC uncollided cross-check** (HANDOFF §5) never run.
- First-flight / convex-patch for the analytic; single energy (14.1 MeV) for the steering results.
- The n=25 QUASR predictor batch (Ginsburg jobs 8885457→8885458) may upgrade the coherence-C /
  headroom predictor from "suggestive (n=5)" to a real correlation — hold that framing until it lands.

## Data/figure discipline

`paper/data/` holds **frozen** result snapshots (CSV/JSON extracted from the RESULTS_tierX.md and the
conformalmap npz) with `PROVENANCE.md` naming the source of every number. `paper/src/` figure scripts
read only from `paper/data/` and write to `paper/figures/`. smplotlib, minimalist. One script per figure,
fixed inputs, regenerable.
