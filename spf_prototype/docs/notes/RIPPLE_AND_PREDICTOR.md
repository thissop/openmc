# TF ripple NWL, the line-of-sight low-pass filter, and the right predictive parameter

*Recorded 2026-07-08. Companion to `LIT_REVIEW_NOVELTY.md`, `LIT_REVIEW_SINGULARITY_SUBTRACTION.md`,
and the anarrima `METHODS_RETROSPECTIVE.md`. This note captures two linked results so they are not
lost again: (1) the toroidal ripple NWL back-of-envelope, and (2) the reframe of "why a cheap
shape→neutronics predictor kept failing."*

## 1. Toroidal ripple NWL is negligible on a normal first wall

Back-of-envelope (`python/ripple_estimate.py`), ITER-like filamentary ring: R0=620 cm, wall minor
radius a=200 cm, N=18 TF coils, field ripple δ=0.5%, source-strength ripple ε_s = δ (any O(δ)
N-periodic source perturbation — strength, radial displacement, or SPF B̂-tilt — enters the same way).

| Quantity | Value |
|---|---|
| Poloidal in/out asymmetry | q_inboard/q_outboard = 0.83; **~19% swing** |
| Toroidal ripple NWL (outboard) | **0.007%** (from a 0.5% source ripple) |
| Toroidal ripple NWL (inboard) | 0.003% |
| source→wall suppression | **~100×** |
| poloidal / toroidal ratio | **~2800×** |

The toroidal NWL ripple is ~**50 ppm** — dead for any realistic first wall.

## 2. The mechanism is a closed-form low-pass filter — this is the keeper

The ~100× suppression is not incidental. A wall point integrates the source over a toroidal window of
half-width ~a/R0 with a Lorentzian 1/d² weight. The Fourier transmission of that kernel at toroidal
harmonic N is

    T(N) ~ exp(−N · a/R0).

For ITER, N·a/R0 ≈ 18·0.323 ≈ 5.8 → exp(−5.8) ≈ 0.003, matching the measured suppression.
**Free-streaming line-of-sight integration is a spatial low-pass filter with cutoff harmonic ~R0/a.**
The ripple harmonic N sits far above cutoff, so the wall never sees it. Channel-independent: strength
ripple, source displacement, and the SPF angular (B̂-tilt) channel are all O(δ) and N-periodic, so all
get the same exp(−N a/R0) haircut.

Consequences:
- Ripple-NWL as a headline application is **dead for a normal wall**; it only revives at tight standoff
  a/R0 ≲ 1/N (wall within ~R0/N of the plasma) — i.e. the **grazing regime the singularity subtraction
  was built for**. The two analytic pieces connect: ripple NWL matters only where you'd need the
  subtraction machinery anyway.
- exp(−N a/R0) is a **closed-form scaling the perturbation series hands you for free** — a citable
  "why the wall can't see the ripple" result even though the number is small.

## 3. Why the cheap shape→neutronics predictor kept failing (C-term, nematic Q)

Prior attempts to predict neutronics from a **pure-plasma global shape descriptor** — the C-term (wrong
order) and a nematic order parameter (tensor alignment of boundary vectors) — did not correlate. The
low-pass result explains *why*, and turns a dead end into a positive result:

- The nematic Q and C-term measure **high-harmonic content of the plasma boundary** (how strongly/
  coherently shaped the plasma is).
- But §2 proves the wall-load harmonic content is **exponentially damped above k ~ R0/a**. The wall
  *structurally cannot resolve* high-harmonic shaping. A predictor built from high-harmonic shape
  content is therefore decoupled from wall-load structure **by construction** — it was never going to
  work.

**ε_eff works because it is the right *class* of parameter.** ε_eff = (real-space excursion)/(source-to-
wall **standoff**) is standoff-based and local, not shape-harmonic and global. Standoff is exactly what
survives the low-pass filter and drives peaking (near-field grazing spike). The failed metrics measured
the *plasma*; ε_eff measures the *plasma-to-wall relationship*.

## 4. What actually predicts global SPF steering benefit

Chaining the project's findings:
- SPF is a **weak global flattener** — peaking is set by a narrow near-field grazing spike
  (0.3–2.8% of wall area) that the P₂ angular lever cannot move.
- The **broad, low-harmonic** part of the load (the ~19% poloidal in/out asymmetry) *is* steerable —
  it is the ±22–43% in/out lever SPF commands, and being low-harmonic it survives the filter.
- Therefore **SPF global benefit is predicted by standoff *uniformity*, not plasma shaping.**
  Generous, uniform standoff (broad in/out-dominated load) → SPF flattens well. A tight grazing spot →
  spike-dominated peaking → SPF can't touch it → weak benefit. The predictor lives in the
  standoff / ε_eff family, computable from geometry with no transport.

**One mechanism (the exp(−k a/R0) low-pass filter) explains three separate negative results:**
the nematic-Q failure, the C-term failure, and the weak global SPF lever. That unification is the
publishable insight; the small ripple number is the evidence for the filter.

## 5. Quantitative test — and what it refuted (`python/predictor_test.py`)

The §4 claim ("standoff predicts SPF benefit; plasma shape is decoupled") was tested against the 5
QUASR devices with local conformal-wall transport. **It did not cleanly validate — reported honestly.**

Y = SPF benefit = (PF_unpol − min_a2 PF)/PF_unpol, vs geometry-only predictors:

| Device | PF_unpol | a2_opt | Y | PF_geo | standoff_peak | spike_area | shaping_amp | nematic_Q |
|---|---|---|---|---|---|---|---|---|
| 803097 | 1.40 | −0.83 | 0.014 | 1.19 | 1.00 | 0.177 | 0.246 | 0.092 |
| 886079 | 1.78 | +0.83 | 0.160 | 1.60 | 1.03 | 0.043 | 0.376 | 0.037 |
| 932746 | 1.38 | +0.22 | 0.043 | 1.34 | 1.00 | 0.095 | 0.159 | 0.005 |
| 59509  | 1.37 | +0.21 | 0.044 | 1.26 | 1.00 | 0.120 | 0.246 | 0.025 |
| 1960314| 1.60 | +0.60 | 0.060 | 1.38 | 0.84 | 0.112 | 0.906 | 0.230 |

Spearman ρ vs Y: **PF_geo +0.90** · standoff_peak +0.30 · spike_area −0.70 · **shaping_amp +0.60** ·
nematic_Q +0.20.

**Findings (honest):**
1. **Best predictor is PF_geo (+0.90):** a transport-free isotropic 1/d² peaking proxy. SPF benefit ≈
   available peaking **headroom** — near-uniform devices gain little; the peaked one (886079) gains most.
   A legitimate cheap predictor, but a *headroom* story, not the standoff-*steerability* story of §4.
2. **The shape metric did NOT decouple (shaping_amp +0.60).** So "the low-pass filter proves plasma
   shape can't predict steering" is **unsupported**. The low-pass filter is validated for the *ripple*
   (§2, directly computed); extending it to "shape is a useless steering predictor" is an untested
   conjecture this data does not back. nematic_Q is weak (+0.20), consistent with the original failure
   but not a clean zero.
3. **The interesting physics is the outlier 803097:** peaked (PF=1.40) but **unsteerable** (Y=0.014,
   10× below the rest) and the **only perp-optimal device** (a2=−0.83). PF_geo flags it as *low* peaking
   (1.19) — so its transport peaking comes from an anisotropic/near-field structure the isotropic proxy
   can't see, and that structure is both perp-favorable and steering-resistant. This is where "when can
   SPF NOT help" lives.

**Dominating caveat: n=5.** Spearman ±0.6 at n=5 is within noise; only +0.90 stands out, and that is 5
points. This test **cannot** discriminate standoff-vs-shape predictors — no statistical power. It does:
(a) kill the premature "shape decoupled" claim, (b) show benefit tracks headroom, (c) flag 803097.
Deciding the predictor question needs a larger device set (batch more QUASR VMEC equilibria on Ginsburg).
