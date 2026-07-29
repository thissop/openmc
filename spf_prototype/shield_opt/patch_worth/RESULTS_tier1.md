# Patch-worth Tier 1 — Attribution A vs Shielding-Worth W (QH coil-20)

Forward mesh flux: `qh_fwd_meshflux.npz` (48M no-WW, mesh [30,30,20], 18000 voxels,
median relerr 3.2%, p90 6.4%). Adjoint: `recip_w35f_coil20_flat` (coil_cells=[20]).
Plasma source S: uniform emissivity (matches the uniform forward source; edge-weighted).
Shield band: rho in [225,265] cm (a_minor~166, R0~1299 cm), 696 voxels.

## A-vs-W rank correlation is LOW and robust across binnings

| grid  | Spearman rho | p     | Chatterjee xi | top-k overlap |
|-------|--------------|-------|---------------|---------------|
| 4x4   | +0.263       | 0.33  | 0.153         | 0.25          |
| 5x5   | +0.107       | 0.61  | 0.014         | 0.33          |
| 6x6   | +0.127       | 0.46  | 0.004         | 0.44          |

**Reading:** attribution (where responsible neutrons are born, S*psi_dagger) and scalar
worth (where blocking helps, phi*psi_dagger in the shield band) are essentially
UNCORRELATED (rho ~ 0.1-0.26, all p > 0.3; Chatterjee xi ~ 0 => near-independence).
This is the "attribution != actionability" direction. Honest caveats: (1) underpowered
per-binning (16-36 patches, none individually significant) -- the strength is the
CONSISTENCY across binnings + xi~0, not any single p; (2) scalar (P0) leading-order worth
-- Tier 2 finite-difference is the confirmation that W (not A) predicts true dose drop;
(3) uniform emissivity for A (edge-weighted); shield band auto-estimated.

Corroborating physics: the kill-shot peak migrated coil20->coil25 (single-coil attribution
over-protects its target), consistent with A not tracking the true lever.

## Tier-2 corners selected (patch_worth_analyze --pick-corners, 4x4)
The decisive disagreement patches: (2,1) near-zero A but HIGHEST W; (2,3) highest A but
modest W. Finite-difference on these tests whether W_FD tracks W (novel) or A (proxy).
Tier-2 array 9227246 (8 patches) building + dosing now.
