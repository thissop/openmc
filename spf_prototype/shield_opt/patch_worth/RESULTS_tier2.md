# Patch-worth Tier 2 — finite-difference verdict (QH coil-20, 4 corner patches)

Baseline: coil_step1b_uniform.npz, coil-20 fast flux R0 = 3.4521e-2 +/- 0.50%.
Patches: +20 cm shield / -20 cm breeder in one (theta,phi) patch (fixed envelope), no-WW 48M,
correlated seed. Tier-1 A,W from the geometry-grounded band (R0=1367 cm, a_minor=217 cm,
band rho in [276,316] cm). W_FD = -(R_patch - R0) / (added_cells * delta_cm).

| patch | corner   | R(coil20) | dR       | dR/sigma | W_FD     | A (attr)  | W (scalar worth) |
|-------|----------|-----------|----------|----------|----------|-----------|------------------|
| (2,1) | loA-hiW  | 1.218e-2  | -2.23e-2 | -101.7   | +9.31e-5 | ~0        | low              |
| (1,0) | hiA-loW  | 9.837e-3  | -2.47e-2 | -124.8   | +8.23e-5 | 2.74e11   | low              |
| (1,2) | loA-loW  | 2.954e-2  | -4.98e-3 |  -17.9   | +1.66e-5 | 0         | mid              |
| (2,3) | hiA-hiW  | 4.868e-2  | +1.42e-2 |  +38.7   | -4.72e-5 | 2.35e12   | 1.81e6 (top)     |

**Spearman W_FD vs W = -1.000 ;  W_FD vs A = -0.400 ;  sign-match(W) = 0.75.**

## The result (surprising, and stronger than "attribution != worth")
Every dR is 18-125 sigma above noise (relerr 0.5-1.1%): the effects are REAL, not sampling scatter.
The patch with the HIGHEST attribution AND the HIGHEST scalar worth -- (2,3) -- is the ONE place
where adding shield made coil-20 WORSE (+41% dose, +38.7 sigma). So on these four corners the
scalar leading-order worth W = phi*psi_dagger is ANTI-correlated with the true finite-difference
worth, and attribution A is uncorrelated/negative too.

**Physics:** the fixed-envelope trade is +20 cm shield AND -20 cm breeder. The scalar worth only
counts the shield's blocking benefit (phi*psi_dagger over the shield band); it omits the breeder
that is removed. FLiBe breeder moderates/absorbs fast neutrons, so at (2,3) thinning it lets more
fast flux reach coil-20 than the added WC shield blocks -> net dose INCREASE. The breeder-for-
shield trade has a SIGN that flips across the device, which neither A (plasma contributon) nor the
leading-order scalar W captures. Only finite difference (or a SIGNED worth that debits the breeder
removal) ranks correctly.

## Caveats (honest)
- Only 4 patches. rho=-1 is a perfect rank reversal but n=4 is underpowered (the printed p=0 is a
  scipy/numpy-2.4.6 artifact; the true two-sided p for rho=+-1 at n=4 is ~0.08). 8-point run
  (patches (1,1),(1,3),(2,0),(2,2)) queued as job 9330876 to firm it up.
- P0 scalar worth. The sign-flip empirically supports the "need a signed/angular worth" branch.
- Only 4 of 8 corners built; (0,*)/(3,*) hit a parastell 29-vs-28 material-tag error (period-
  boundary patches) -- separate fix pending.

## Implication for the method
The corrected worth is not phi*psi_dagger but ~ (Sigma_shield - Sigma_breeder) * phi * psi_dagger
-- a SIGNED trade worth. Building/validating that signed worth against these finite differences is
the next methodological step and the likely paper headline: naive contributon placement (attribution
OR scalar worth) can be actively counterproductive under a fixed-envelope breeder-for-shield trade.
