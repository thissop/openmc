# LIMITATIONS — what this prototype does and does NOT model

This is a self-contained research prototype to get the spin-polarized DT
birth-direction physics right and validated before any integration. Read this
alongside any figure or number: several plainly-stated boundaries matter for how
far a result can be pushed.

## The model is a *geometric, free-streaming* NWL — not a dose, not a TBR
- **Scattering is OFF.** The vessel interior is near-void (Leakage = 1); neutrons
  travel in straight sightlines from the source to the first wall. The "neutron
  wall load" (NWL) here is therefore a **geometric** quantity (uncollided current
  per unit wall area). It is **not** a shielded dose, **not** a damage/dpa rate,
  and **not** heating in structure — all of which require transport through real
  materials.
- **No breeder blanket.** The inboard/outboard "current fractions" are a
  **geometric precursor** to a tritium-breeding argument (where do the source
  neutrons land), **not** a TBR. A real TBR needs a Li-bearing breeder blanket and
  full transport (scattering, multiplication, ⁶Li(n,α)t); none of that is here.

## Source physics
- **Purely toroidal field only** (plus a constant-direction option). The angled /
  pitched field (Schwartz §3.2, the α/β kernels) is **not exercised** in the
  source or the analytic comparison — even though the C++ sampler's Gram-Schmidt
  rotation already handles an arbitrary B̂ direction. So no magnetic pitch, no
  spatially varying field.
- **Polarization `(a,b,c)` is spatially uniform user input.** No depolarization
  physics during confinement; the fuel polarization is taken as a fixed global.
- **Neutron channel only.** No alpha-emission channel, no alpha directing, no
  knock-on/beam-target fusion, no center-of-mass motion of the reactants.
- **Single birth energy** (14.1 MeV, monoenergetic). The Ballabio thermal-
  broadened spectrum is **not** in the verification path (energy is decoupled from
  the polarized angular distribution to leading order and irrelevant to the NWL
  geometry).
- **Filamentary / parabolic-ring source.** The plasma is a `(1−ρ²)`-weighted sum
  of toroidal rings over a circular cross section — no realistic equilibrium, no
  flux-surface shaping, no profile beyond `(1−ρ²)`.

## Geometry
- **Rectangular, convex cross section only.** A square box torus (Schwartz §3).
  No D-shape, no divertor, no non-convex / self-occluding first wall. The
  visibility integrals assume the convex-cross-section sightline geometry.
- **Aspect ratio 2.0**, matching anarrima's `examples/square_torus.py` (which
  generates Schwartz Fig. 2), not the 2.5 stated in the paper text (see
  `RESULTS_tier2a.md`).

## Interpretation caveats surfaced by the results
- **Polarization RAISES the overall peaking factor in this FIXED geometry**
  (iso 1.33 → A 1.65, B/C 1.61): directional emission concentrates the load, so
  isotropic is the most uniform. The B/C benefit here is the **center-stack
  (inboard) load reduction** (inboard peak −41%, inboard current 8.7%→5.1%), NOT a
  global peaking reduction. Schwartz's peaking-factor *reduction* comes from
  re-shaping the first wall to exploit the steering — not modeled here.
  (`RESULTS_tier3.md`.)
- The Tier-2b per-bin directionality residuals have stdev ≈ 0.65–0.7 (just under
  1): the 10-batch tally σ is mildly conservative for the mode/iso ratio, so
  OpenMC agrees with the analytic even slightly better than the error bars imply.
  No coherent bias (means ≈ 0).
- Absolute (non-ratio) per-bin residual normality at very high statistics (≳1e8)
  is dominated by bin-center-vs-bin-integral discretization, not physics; we
  therefore compare directionality ratios at ~2–3%/bin statistics.

## Dependencies
- `anarrima` is **not on PyPI**; installed from GitHub
  (PrincetonUniversity/anarrima). Our independent `scipy.quad` engine is the
  fallback and agrees with anarrima to ~1e-6 (to 1.8e-15 on a single ring).
