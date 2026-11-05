# RUNG 5 — corrected nested design basis (replaces the invalid 20/12/35 sequence)

## Audit finding (what was wrong)
The earlier A/B/C "hierarchy" (20/12/35 DOF) was **not a nested resolution hierarchy**. It mixed basis
*families*: Level A was **indicator sectors** (discontinuous top-hats), Levels B/C were **Fourier modes**.
Sector indicators are not a subspace of the Fourier span, so `span(A) ⊄ span(B)` — the sequence had no nesting
relationship and could not support a spatial-resolution convergence argument. Calling basis *generation*
"spatial-resolution convergence" was also wrong: convergence is established only by **constrained optimization +
nonlinear transport evaluation** of successive spaces, showing both the objective value and the recovered
thickness field stabilize within Monte-Carlo uncertainty. That test has not been run yet.

## Corrected: ONE genuinely nested stellarator-symmetric 2-D Fourier ladder
Single family, `t(θ_pol, φ_tor) = Σ c_{m,n} cos(m·θ − n·nfp·φ)` on the ParaStell 13(tor)×19(pol) grid, nfp=4.
Stellarator symmetry `t(θ,φ)=t(−θ,−φ)` keeps cosines (VMEC convention; n±  for m≥1 carries both helicities,
m=0 keeps n≥0). Truncation index set `S(M,N)={(0,n):0≤n≤N} ∪ {(m,n):1≤m≤M,−N≤n≤N}`.

| level | (M,N) | DOF | λ_pol,min | (vs 227 cm coil) | λ_tor,min | sub-coil poloidal? |
|-------|-------|-----|-----------|------------------|-----------|--------------------|
| F0 | (1,1) | 5   | 1363 cm | 6.0× | 2147 cm | no |
| F1 | (2,2) | 13  | 682 cm  | 3.0× | 1074 cm | no |
| F2 | (4,3) | 32  | 341 cm  | 1.5× | 716 cm  | no |
| F3 | (6,4) | 59  | 227 cm  | 1.0× (≈1 coil) | 537 cm | no (marginal) |
| F4 | (9,6) | 124 | 151 cm  | 0.67× | 358 cm  | **yes** |

**Proof of nesting (both):** (1) literal — each level's mode-name set is a *strict subset* of the next
(F0⊂F1⊂F2⊂F3⊂F4; |·|=5<13<32<59<124). (2) numeric — least-squares fit of an F1 mode in the F2 basis has
max residual 6.3e-15 (machine zero). ⇒ `span(F_K) ⊂ span(F_{K+1})` exactly.

## Sub-coil level and the exact GEOMETRIC feasibility cap
F4 is the first clearly sub-coil-poloidal level (λ_pol 151 cm < 227 cm coil footprint). Going finer is **capped
by the grid**, not by choice: the 13×19 surface has poloidal Nyquist m=9 and toroidal Nyquist n=6. Consequences:
- **Poloidal:** sub-coil reachable but only to m=9 (λ_pol 151 cm) — "modestly" sub-coil, as flagged.
- **Toroidal:** the 8-coils-per-period structure needs n≈8 (λ_tor 268 cm); the grid maxes at n=6 (358 cm), so
  the current surface **cannot resolve individual coil toroidal spacing**. A genuinely sub-coil-in-BOTH design
  study requires **re-meshing ParaStell finer** (≥17 toroidal, ≥21 poloidal). This is an exact geometric limit.

## Separate basis FAMILIES (robustness checks, NOT members of F0..F4)
- **Sector indicators** (16 coarse patches): interpretable, discontinuous — a basis-family robustness check.
- **B-splines** (smooth, compact support): a *later, separate* family for the same convergence test.
These are run to confirm results are basis-family-independent; they are not a continuation of the Fourier ladder.

## Files
`nested_basis.py`, `nested_basis_modes.npz` (all mode arrays), `nested_basis_manifest.json` (levels, proofs,
feasibility). The old `hierarchical_basis*` (20/12/35) is superseded — kept only for provenance.
