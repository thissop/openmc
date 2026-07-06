# RESULTS — Tier 8 / Part A (angled/pitched field: the rigor anchor)

Square-torus free-streaming (geometry == analytic: R0=1.0, a=0.5, u=0.4, w=1.6, z=∓0.6); 4,000,000 histories/config, 50×50 poloidal bins/face. Pitched field B̂=cosβ·φ̂+sinβ(cosα·R̂−sinα·ẑ); test pitch (α,β)=(0.6,0.5) rad. Walls: inboard/outboard/floor (α≠0 breaks up-down symmetry → ceiling deferred).

## Rung 1 (analytic): generalized quad vs anarrima angled kernels g_*a

| check | max rel error | gate |
|---|---|---|
| single-ring quad vs anarrima (inb/out/floor × A/cos2/BC × 3 pitches) | 3.6e-14 | <1e-6 |
| parabolic-plasma quad vs anarrima_vec (α,β)=(0.6,0.5) | 2.4e-08 | <1e-6 |
| β=0 reduction: angled quad vs toroidal quad | 0.0e+00 | ≈0 |

The (α,β) convention (B̂ above) was pinned to anarrima's kernels to ~2e-14; β=0 is α-independent and machine-exact, so the regression below is unambiguous.

## Rung 0 (transport regression): bmode=angled,β=0 reproduces toroidal

- **max rel |I(angled,β=0) − I(toroidal)|** over iso/A/B × inb/out/floor = **0.0e+00** (bit-identical B̂ ⇒ identical RNG stream ⇒ exact recovery).
- toroidal A midplane oracle: inboard **+39.2%**, outboard **-21.0%** (Tier-2b: +39%/−21%; Schwartz text +43%/−22%).

## Rung 1 (transport): free-streaming OpenMC vs angled analytic

Per-bin standardized residual of the mode/iso directionality vs the angled analytic D_mode/D_iso (mean≈0, stdev≈O(1) ⇒ reproduced within statistics).

| config | residual mean | residual stdev | N bins |
|---|---|---|---|
| A/iso (angled) | -0.031 | 0.646 | 150 |
| B/iso (angled) | -0.002 | 0.818 | 150 |

| angled A directionality vs iso | analytic | OpenMC |
|---|---|---|
| inboard midplane | +21.9% | +21.4% |
| outboard midplane | -16.3% | -16.6% |

**Tier-8/Part-A gate: angled quad matches anarrima g_*a (<1e-6); β=0 recovers toroidal exactly; free-streaming OpenMC reproduces the angled analytic directionality within statistics. The verified pitched-field source is the bridge to a real equilibrium B̂(x).**
