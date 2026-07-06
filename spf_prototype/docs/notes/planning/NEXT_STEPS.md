# NEXT_STEPS — ranked, honest about cost

Three real directions beyond this prototype. Ranked by recommended order
(readiness × value), with an honest note on whether each is a **cheap extension
of existing machinery** or a **new physics tier**.

## 1. Angled / spatially-varying B̂  (cheapest; the integration bridge)
**What:** replace the purely-toroidal field with a pitched/spatially-varying B̂
(toroidal + poloidal component, i.e. Schwartz §3.2 α/β). This is the natural
bridge to integrating with Peterson's physics-informed tokamak source, which
supplies a real equilibrium field.

**Cost: cheap extension of current machinery.**
- The C++ sampler already emits correctly for **any** B̂ — the Gram-Schmidt
  rotation in `spf_sampler.hpp` is tested with off-axis fields (Tier-1 §3.3). Only
  a field model needs adding to `polarized_fusion_source.cpp` (read α,β or a
  field map; the `bmode=constant` path already proves arbitrary directions work).
- The analytic cross-check is already available: anarrima ships the angled
  kernels `g_HAa/g_Hca/g_HBa` and `g_VAa/...` (the `a` = angled-field forms,
  Schwartz Eqs. 15–17, 21–23). `analytic_nwl.py` would call those instead of the
  toroidal ones.
- Stays in the validated free-streaming regime, so the analytic ground truth
  still applies — verification is straightforward.

## 2. Scattering ON with a real wall material  (moderate; first step toward dose)
**What:** give the vessel/first wall a real material and let neutrons scatter.
This is the **first result OpenMC produces that anarrima structurally cannot**
(anarrima is an uncollided-sightline integral), and the first honest step from a
geometric NWL toward an actual wall **dose / heating / dpa**.

**Cost: cheap machinery, but a new verification problem.**
- Mechanically small: swap the void cell for a material (e.g. steel/W first wall),
  add the wall as a real cell, switch the tally from surface current to wall
  heating / flux / reaction rate.
- The hard part is **losing the analytic ground truth**: once scattering is on
  there is no closed-form reference, so verification shifts to convergence
  studies, energy-spectrum sanity, and limiting cases (thin-wall → free-streaming
  recovers Tier-2b). Plan a verification strategy before trusting numbers.
- This is where the directionality result becomes physically actionable (does the
  steering survive scattering and reduce real center-stack heating?).

## 3. Breeder blanket + transport for an actual TBR  (new physics tier; biggest)
**What:** add a Li-bearing breeder blanket and compute a real tritium breeding
ratio, then ask how polarization steering changes it.

**Cost: a full new physics tier.**
- Needs blanket geometry + materials (Li enrichment, multiplier e.g. Be/Pb,
  structure, coolant), full neutron transport with the right cross sections, and
  ⁶Li(n,α)t / ⁷Li(n,n'α)t **tritium-production tallies** plus leakage accounting.
- The current inboard/outboard **current fractions are only a geometric precursor**
  to this — turning "where neutrons land" into "tritium bred per source neutron"
  requires all of the above. Highest effort, highest "real-reactor" payoff.

---
**Summary:** #1 is a near-term extension (sampler is ready, analytic check exists)
and the cleanest path to integration; #2 is the first physically-meaningful step
beyond free-streaming and the first OpenMC-only result; #3 is the headline reactor
metric and a genuine new tier of work.
