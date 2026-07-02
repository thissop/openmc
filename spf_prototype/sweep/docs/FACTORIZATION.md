# FACTORIZATION: eta ~ eta_source(C) * A(tau)  (skeleton, completed in stage 7)

## Hypothesis
The directional efficiency factorizes into a field-geometry part and a blanket part:

    eta(C, blanket) ~ eta_source(C) * A(tau_scatter)

- eta_source(C): the field-geometry cancellation, measured from a near-source
  (pre-blanket) observable. PRIMARY, fast-converging.
- A(tau): a single scalar set by the scattering optical depth / scattering power
  between source and detector. Measured as A = eta_coil / eta_source per config, for
  free from data already generated.

## Why the blanket acts mainly through scattering optical depth
Prior precise_QA run: the directional signal fell from 13.5% (free-streaming) to
8.5% (with scattering), a ~0.63 survival. Scattering randomizes direction, so it
washes out the SPF directional bias roughly in proportion to how many scatters a
neutron takes reaching the detector. That is an optical-depth statement, so A should
depend on tau and little else if the field geometry and blanket separate.

## How separability is TESTED
1. Two eta observables per config: eta_source (near-source, pre-blanket) and
   eta_coil (deep, full blanket). A = eta_coil / eta_source.
2. A small blanket subset (2 to 3 recipes: baseline, lighter/thinner,
   heavier/more-hydrogenous). If eta_source(C) is INVARIANT across recipes while
   only A moves along an A(tau) curve, separability holds.

## The three ways separability can break (report which, if any, is active)
1. ENERGY COUPLING: source energy-angle correlation differs by config and the
   blanket attenuates energy-dependently. Signature: eta_source(C) invariant but
   A(tau) scatters by config in a way correlated with the source spectrum.
2. GEOMETRIC NON-SEPARABILITY: the residual steered-flux DIRECTION depends on C, and
   the nonuniform-offset blanket attenuates anisotropically. Signature: A depends on
   C at fixed tau.
3. NEAR-SOURCE OBSERVABLE CONTAMINATION: eta_source is not truly pre-blanket.
   Signature: eta_source(C) itself slopes with blanket recipe.

## Outcome
- If eta_source(C) is invariant across recipes: report eta ~ eta_source(C) * A(tau).
- If it slopes: report the C-blanket COUPLING as the finding (field geometry and
  blanket cannot be designed independently for SPF). Either way is publishable; the
  design DETECTS which world we are in and never assumes.
