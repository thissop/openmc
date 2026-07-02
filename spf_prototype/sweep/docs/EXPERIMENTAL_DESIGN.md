# EXPERIMENTAL DESIGN (skeleton, completed in stage 7)

## The control: one blanket, applied identically to every config
The baseline blanket recipe (thicknesses, materials, Li-6 enrichment) is fixed and
applied identically to every config in the C sweep. This uniformity is the
experimental control on eta: it isolates the field-geometry effect from any
per-config blanket variation. Note that C is blanket-INDEPENDENT by construction
(field plus source only), so the uniform blanket protects the clean MEASUREMENT of
eta, not C.

## Config family (stage 6 selects it)
Compute C for candidate QUASR devices FIRST, then tile the C axis roughly uniformly
(do not cluster). Span nfp and iota bands. Carry symmetry class QA/QH/QI as metadata.
QA is the core; include a QH arm and a QI arm to test whether eta collapses onto ONE
eta_source(C) curve regardless of class (universality), or whether the classes
separate (then the tensor structure, not scalar C, is the predictor). Include an
axisymmetric / tokamak-limit anchor with C ~ 1 for the eta normalization.

Note on symmetry class: QUASR's catalogue does not tag QA/QH/QI in its columns, so
the class is carried explicitly per device in the config manifest (the family is
chosen by a human who knows each device's class), not inferred by the loader.

## The eta definition (single source of truth)
eta is the fractional change of a directional observable from unpolarized to
polarized, normalized to the same fractional change at the high-coherence (C ~ 1)
anchor. Two observables:
- eta_source: near-source, pre-blanket (free-streaming first-wall directional
  anisotropy). PRIMARY for the scaling law.
- eta_coil: coil fast flux, full blanket, deep-penetration, VR-pending.
A = eta_coil / eta_source.

## Two-paper plan
- Paper 1 (validation and methods): verify OpenMC SPF in the FREE-STREAMING neutron
  wall loading (NWL) limit against the independent analytic/numeric NWL series on
  matched precise_QA patches (same boundary, source angular convention, patch
  definition, normalization). The pymoab-free DAGMC pipeline is the enabling method;
  the automated conformal blanket builder is a supporting capability.
- Paper 2 (physics): the C to eta scaling law, universality across QA/QH/QI, and the
  factorization completeness arm, with coils-as-detectors dose in an appendix.
