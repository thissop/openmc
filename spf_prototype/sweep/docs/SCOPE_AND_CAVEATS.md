# SCOPE AND CAVEATS (skeleton, completed in stage 7)

Honest calibration of what this study claims and does not claim.

## The conformal blanket builder is not novel geometry
Normal-offset surfaces and conformal stellarator blankets already exist (ARIES-CS,
HELIAS). What is defensible here is an automated, toolchain-independent
(pymoab-free) pipeline that emits DAGMC-ready conformal blankets uniformly across a
large device sweep, with the hand-written DAGMC writer as the enabling piece. The
claim is throughput and portability, not a new blanket concept.

## Blankets are representative, not design-optimal
Because the blanket is a uniform normal-offset shell applied identically to every
config and is NOT individually optimized per config, absolute TBR and dose are
REPRESENTATIVE, not design-optimal. The study isolates the field-geometry effect on
SPF steering; it does not optimize the reactor.

## Future work (scoped)
Per-configuration blanket optimization (breeding adequacy, shield sizing, coil
standoff) is scoped future work and a candidate for a fusion DESIGN COURSE (PI
intends to use this in fall 2027).

## QUASR coil caveat
QUASR coils are optimized for the magnetic objective, not for a blanket or a build,
so coil-to-plasma standoff may be tighter than a reactor. Coil dose is therefore a
RELATIVE cross-config comparison, not an absolute reactor number.

## The uniform shell quietly protects the factorization
A uniform offset shell attenuates more isotropically than a lumpy, ported blanket.
This makes the A(tau) factorization cleaner than it would be for a real engineered
blanket, so separability holding here is necessary but not sufficient evidence that
it holds for an engineered device. Stated so the completeness result is not
over-read.

## Field fidelity
b_hat comes from exact Biot-Savart on the QUASR coil filaments, audited by
div B ~ 0 and B.n ~ 0 on the LCFS. This removes the fitted-field risk. The remaining
field caveat is that we use the vacuum coil field (no plasma response); for the
low-beta QUASR equilibria this is a small correction, noted but not modeled here.

### Coil-fit quality varies by device, and the audit filters the family (methods note)
The LCFS residual B.n/|B| is not a code error but a real property of each QUASR coil
set: it is grid-converged and genuinely a few percent for aggressive (high-nfp,
high-iota) QH configs, and discretization-limited to ~1e-3 or below for gentle QA
configs (device 59509: rms 3e-4 at a fine grid, still falling). We HARD-GATE exclude
configs whose coils leave more than ~8% normal field on the nominal boundary, because
b_hat is then too unreliable to compute a trustworthy C. Two honest consequences to
state in the paper:
- This is a data-quality FILTER on the config family, so the exact exclusion count
  and threshold belong in the methods.
- The exclusion CORRELATES with the independent variable: the most aggressive QH
  devices have both the lowest C and the largest B.n, so they are the most likely to
  be cut. The low-C end of the family is therefore somewhat filtered toward the
  better-fit QH devices, a mild selection effect to acknowledge when reading
  eta_source(C) at low C. (C itself is an integral over the source and is robust to a
  few-percent local field error, so passing configs are trustworthy; the concern is
  only which devices survive to populate the low-C end.)
A fix that removes the caveat entirely is to take b_hat from the QUASR VMEC/Boozer
equilibrium field (which is a flux function by construction) instead of the vacuum
coil field, or to re-solve the coils to a tighter normal-field tolerance; either is
scoped future work, not required for the coherence study.
