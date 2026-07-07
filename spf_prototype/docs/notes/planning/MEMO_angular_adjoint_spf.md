# MEMO — Angular-dependent adjoint / random-ray importance for SPF neutron steering

**Status:** Phase-3 (end-of-summer) idea, captured now, not to be implemented yet.
**Author context:** SPF prototype (this repo). Depends on an external feature.
**One-liner:** Pair a *direction-resolved* adjoint importance function `ψ†(r,Ω)` (from
OpenMC's random-ray solver, angular extension in progress by E. Peterson's group) with
our controllable-direction SPF source to turn "does polarization help?" into "**aim the
polarization optimally to protect the magnets.**"

---

## 1. Concept in plain terms

Spin-polarized DT fusion gives a source whose **emission direction is a control knob**.
About the local field `b̂(r)` the birth-angle kernel is anisotropic and tunable via the
fuel polarization:

- **perpendicular mode** (our "A") ∝ `sin²θ_B` — emission has a **true null along ±b̂**;
- **parallel mode** (our "B/C") ∝ `1/4 + 3/4 cos²θ_B` — peaks along `±b̂`, minimum (¼) ⟂ `b̂`;
- general mix: `w(θ_B) ∝ 1 + a₂ P₂(cosθ_B)`, with the polarization mapping to a single
  scalar `a₂ ∈ [−1, +1]` per point (`a₂=−1` pure perpendicular, `a₂=+1` pure parallel).

An **adjoint importance** `ψ†(r,Ω)` for a chosen damage response (fast-neutron fluence /
DPA / kerma at the toroidal-field magnets) answers the dual question: *a neutron born at
`r` heading in direction `Ω` — how much does it ultimately contribute to that damage?*
Ordinary (scalar) adjoints give only `ψ†(r)`; the **angular** adjoint resolves the
`Ω`-dependence, which is exactly the axis SPF controls.

**Synergy:** where `ψ†(r,Ω)` is large in some direction band, we want the SPF source to be
*small* there. Because `w` has an angular null/minimum we can orient (via `a₂`), we can
suppress birth flux precisely into the high-importance directions. The adjoint tells us
*which* directions matter at *which* locations; SPF lets us *avoid* them.

## 2. Math sketch

Response (damage) as the source–importance overlap:

```
R = ∫_V ∫_4π  S(r,Ω) ψ†(r,Ω)  dΩ dr .
```

SPF factorizes source strength from direction:

```
S(r,Ω) = S₀(r) · w(Ω·b̂(r); a₂(r)),     w = 1 + a₂ P₂(Ω·b̂),   (1/4π)∫ w dΩ = 1
```

(the P₂ term integrates to zero over the sphere, so the *shape* knob `a₂` leaves the
angle-normalization fixed; the overall rate factor `η(a,b,c)` scales `S₀` separately, as
in the sampler). Insert and split the importance into its angle-average plus its P₂
moment **about the local field**:

```
ψ̄†(r) = (1/4π) ∫ ψ†(r,Ω) dΩ                       (angle-averaged importance)
m₂(r) = (1/4π) ∫ P₂(Ω·b̂(r)) ψ†(r,Ω) dΩ            (P₂ moment of ψ† along b̂)

⇒  R = ∫_V S₀(r) [ ψ̄†(r) + a₂(r)·m₂(r) ] dr .
```

**Optimality is then explicit.** Minimizing `R` at fixed yield, pointwise over
`a₂(r) ∈ [−1,+1]`:

```
a₂*(r) = − sign( m₂(r) ) .
```

- `m₂ > 0` — importance is concentrated **along ±b̂** (damage flows along field lines):
  pick `a₂ = −1`, the **perpendicular mode**, whose null sits exactly on `±b̂`.
- `m₂ < 0` — importance is concentrated **perpendicular to b̂**: pick `a₂ = +1`, the
  **parallel mode**, minimum in the equatorial band.

The key physical constraint: `b̂(r)` is fixed by the equilibrium — we do **not** get to
aim the null anywhere, only to choose *which side* of the P₂ node to suppress. Note the
node structure: the perpendicular mode vanishes along `b̂` (a real zero), while the
parallel mode only dips to ¼ transverse to `b̂` (a minimum, not a zero) — so perpendicular
polarization is the sharper protective tool when high-`ψ†` directions align with the field.

## 3. The OpenMC path

1. **Adjoint solve (once):** run the random-ray solver in adjoint mode with an adjoint
   source = the magnet damage response (DPA/kerma weighting on the TF-coil region),
   scattering **on**. Output: `ψ†(r,Ω)` resolved in angle (spherical-harmonic / discrete-
   direction moments), on the random-ray spatial mesh. This is the external dependency —
   the angular extension Peterson's student is adding.
2. **Field map:** evaluate `b̂(r)` on the same mesh from the equilibrium we already carry
   (the SPF source consumes `b̂(r)` today), to form `Ω·b̂` and hence `m₂(r)`.
3. **Couple:** our `PolarizedFusionSource` is the forward `S(r,Ω)`. Compute `ψ̄†` and
   `m₂` moments per cell from the adjoint output, combine with `S₀(r)` per §2.
4. **Optimize:** pick `a₂*(r)` (or a single global `a₂*` under the realistic constraint
   that the fuel carries one global polarization: `a₂* = −sign(∫ S₀ m₂ dr)`), then confirm
   with one forward SPF transport run.

Reciprocity does the heavy lifting: with scattering on, `ψ†` already folds in *multi-
scatter* damage paths, so this captures the physics our free-streaming geometric NWL
picture (Tiers 1–2b) cannot, and directly addresses the ~3× steering dilution seen once
real scattering was turned on (Tier 5).

## 4. Payoff

- **Magnet-protection optimization.** A principled, per-location prescription for the
  polarization that minimizes coil damage — not a guess-and-check parameter scan.
- **Cheap screening quantity.** The overlap integral `∫ S·ψ† dΩ dr` and its sensitivity
  `dR/da₂ = ∫ S₀ m₂ dr` are **quadratures over precomputed fields**, not transport runs.
  One adjoint solve replaces a full forward SPF polarization sweep: evaluate any candidate
  `(a,b,c)` by integration, and read the *achievable* damage reduction straight off `m₂`
  and `ψ̄†` before running anything. This is the screening lever — rank polarizations, and
  even bound the maximum benefit, essentially for free.
- **Falsifiable "does SPF help here?"** If `m₂(r)` is tiny wherever `S₀` lives, SPF cannot
  help that response, and we learn it from one adjoint solve.

## 5. Dependencies / unknowns — why phase-3

- **External feature not yet in hand.** The angular-adjoint random-ray capability is being
  developed by Peterson's student; not merged/validated. Timeline and API are out of our
  control — this memo commits nothing until it lands.
- **Angular resolution.** The method only needs the P₂ moment `m₂`, so the adjoint must
  resolve at least quadratic angular structure (≥ P₂ moments / enough discrete directions);
  random-ray ray-effects and angular truncation must be checked to not alias the P₂ signal.
- **Response definition.** Need an agreed magnet damage metric (DPA vs fast fluence vs
  kerma) and its region as the adjoint source — a modeling choice, not yet fixed.
- **Realizable polarization.** `a₂ ∈ [−1,+1]` assumes ideal polarization and a single
  global fuel state; real depolarization shrinks the range and pointwise `a₂(r)` is likely
  *not* independently controllable — the practical optimum is the constrained global scalar.
- **Geometry / frame consistency.** `ψ†` lives in global `Ω`; `w` in the `b̂`-local frame.
  Both must share cross sections, geometry, and the same `b̂(r)` map — plumbing that only
  pays off once both halves independently exist and are trusted.
- **Sequencing.** This is a *synthesis* step: it presupposes (a) the merged angular adjoint
  solver, (b) our SPF source validated scattering-on in realistic geometry, and (c) a
  defined response. Hence end-of-summer, after the forward SPF neutronics is solid — not now.

---

*Cross-refs: SPF source physics and (a,b,c)↔(w_perp,w_par) mapping — `python/abc_modes.py`,
`docs/reference/PRIOR_WORK_BAE2025.md`; scattering-on steering dilution — Tier 5
(`docs/notes/results/RESULTS_tier5.md`); free-streaming geometric baseline — Tiers 1–2b.*
