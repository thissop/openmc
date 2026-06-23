# LIMITATIONS — what this prototype does NOT do

This is a self-contained research prototype to get the spin-polarized DT
birth-direction physics right before any integration. Deliberately out of scope
(matching Schwartz 2025's own caveats plus our verification choices):

- **No scattering in the verification.** Tier-2b runs a near-void interior
  (vacuum cells, Leakage = 1) so neutrons free-stream from source to wall,
  isolating source-direction physics. A real first-wall load needs transport
  through blanket/structure (multiple scattering, (n,2n), etc.).
- **Single birth energy.** Monoenergetic 14.1 MeV. The Ballabio thermal-broadened
  Gaussian is not implemented (energy is decoupled from the polarized angular
  distribution to leading order and irrelevant to the NWL geometry here).
- **Filamentary / weighted-ring source only.** The plasma is a parabolic
  (1−ρ²) sum of toroidal rings over a circular cross section; no realistic
  equilibrium, no flux-surface shaping, no profile beyond (1−ρ²).
- **Constant or purely toroidal B̂ only.** The angled-field (α,β) forms exist in
  anarrima and in the C++ header's reach, but the prototype's source + analytic
  comparison use only toroidal (and a constant-field option). No spatially
  varying / realistic equilibrium field.
- **Neutron channel only.** No alpha emission, no alpha directing, no knock-on or
  beam-target fusion, no center-of-mass motion of the reactants.
- **No depolarization.** Fuel polarization is treated as fixed input (a,b,c);
  no depolarization physics during confinement.
- **Convex, rectangular cross section only.** The square-cross-section torus
  (Schwartz §3). No divertor, no non-convex first wall, no stellarator geometry.
- **Geometry aspect ratio 2.0**, matching anarrima's `examples/square_torus.py`
  (which generates Fig. 2), not the 2.5 stated in the paper text. See
  `RESULTS_tier2a.md`.

## Known numerical caveats
- The Tier-2b per-bin directionality residuals have stdev ≈ 0.65–0.7 (slightly
  below 1): the 10-batch tally σ is mildly conservative for the mode/iso ratio,
  so OpenMC agrees with the analytic even a bit better than the error bars imply.
  No coherent bias (means ≈ 0).
- Absolute (non-ratio) per-bin residual normality at very high statistics (≳1e8)
  is dominated by bin-center-vs-bin-integral discretization, not physics; we
  therefore compare directionality ratios at the spec's ~2–3%/bin statistics.
- `anarrima` is **not on PyPI**; it is installed from its GitHub repo
  (PrincetonUniversity/anarrima). Our independent `scipy.quad` engine is the
  fallback and agrees with anarrima to ~1e-6.
