# Prior work: Bae et al. 2025 vs. this prototype (honest positioning)

**Bae, Borowiec, Badalassi, Parisi, Diallo, Menard, Khodak, Brown,
"Neutronics analysis of spin-polarized fuel in spherical tokamaks,"
*Nuclear Fusion* 65, 086051 (2025)** (ORNL + PPPL; DOI 10.1088/1741-4326/adf3c6).
This is the closest prior work and it predates and overlaps our Tiers 5–6. **Our
work is independent validation + tooling, not new physics.** Cite Bae up front.

> **Verified against the full PDF** (`Neutronics analysis of spin-polarized fuel
> in spherical tokamaks.pdf`, *Nucl. Fusion* 65 086051, 9 pp, pub. 1 Aug 2025).
> Everything below is confirmed from the text/tables/figures, not the abstract.

## What Bae et al. did
- **Code:** OpenMC (ENDF/B-7.1), 3D Monte Carlo, MAGIC variance reduction; CAD →
  Cubit → DAGMC geometry via the FERMI workflow.
- **Source:** 10 M pre-generated samples (position/angle/energy) per scheme from
  physics-based emission functions, exported to a **point source / fixed sample
  bank** via the OpenMC **Python API**; **40 M total particles** transported.
- **Machine:** STAR spherical tokamak — R=4 m, a=2 m (aspect=2), 800 MW,
  κ=2.3, δ=0.65, B_T,0=5.2 T; 90° reflective sector.
- **Blanket:** **Pb-Li**, 92% ⁶Li-enriched; homogenized single-shell layers.
- **Modes (Bae Table 1) and the mapping to OUR modes — make this airtight:**
  | Bae name | Bae W(θ) | Bae rel. σ | = our mode |
  |---|---|---|---|
  | unpolarized | 1 | σ₀ | **iso** |
  | parallel (anti-aligned) | (1+3cos²θ)/2 | σ₀ | **B/C** (∝1+3cos²θ) |
  | perpendicular (spin-aligned) | (9/4)sin²θ | 3σ₀/2 | **A** (∝sin²θ, +50% rate) |
  ⚠ **Symbol-collision warning:** Bae's (a,b,c) are the **Hupin–Navrátil vector/
  tensor polarizations** (a=d₁−d₋₁, b=t₁/₂−t₋₁/₂, c=1−3d₀), **NOT** Schwartz's
  collision-mode fractions (a,b,c) that we use. Never equate the symbols — compare
  by emission shape W(θ)/our-shape and by name only. Bae has no "mode C."
- **TBR:** unpol **1.082**, perpendicular **1.062** (−1.8%), parallel **1.111** (+2.7%).
  Inboard/outboard split (Table 3): parallel steers outboard (inboard 0.1777 vs
  0.1862 unpol; outboard 0.9338 vs 0.8954). Their stated **mechanism** for parallel
  giving the *highest total* TBR is **geometric and ST-specific**: the parallel
  scheme pushes flux outboard, and the ST outboard breeder is **~78.5% of total
  breeder volume** while the inboard is space-starved → more neutrons land where
  there is more blanket. (Our square torus is far more symmetric → same sign but
  much smaller magnitude; see below.)
- **Magnets:** parallel cuts inboard TF-coil flux **40%** → **+68% magnet lifetime**
  (linear-fluence assumption); perpendicular **+30%** flux (worse). CCFE-709 group
  spectra similar in shape across schemes, parallel lower in magnitude.
- **Held total neutron rate/power constant**; explicitly notes the perpendicular
  case "would have an approximately 50% increase in fusion reactivity," *not modeled*.
- **Did show:** a **pre-scatter "first-interaction location" poloidal distribution**
  (Fig. 4, first 200k neutrons) — i.e. a free-streaming-like view — *separately*
  from the full-transport TBR/flux; and a **sampled-pdf-vs-W(θ)sinθ overlay**
  (Fig. 3) validating their source bank (the analog of our Tier-1 χ² check).
- **Did NOT:** report wall heating/dpa/down-scattered spectrum; **tabulate the
  free-streaming→transport dilution *ratio***; optimize design; couple to
  plasma/burn-efficiency; do tritium transport; build an independent analytic NWL.

## Point-by-point vs. our prototype
| aspect | Bae et al. 2025 | ours | verdict |
|---|---|---|---|
| MC neutronics of SPF | OpenMC, sample bank, point source | OpenMC, **compiled C++ `Source` plugin** (on-the-fly, any N, any B̂) | **they're first**; ours = independent reimplementation (different route) |
| TBR per polarization | parallel **+2.7%**, perp −1.8% | B/C **+0.6%**, A −0.4% | **same sign/trend**; our smaller magnitude is *explained*: Bae's effect is amplified by the ST's inboard-starved blanket (outboard ~78.5% of volume), our square torus is ~symmetric → **independent cross-validation, not novel** |
| center-stack benefit | **+68% magnet life**, −40% inboard flux | B/C −15% inboard incident load (w/ scattering) | same physics; **they quantified magnet life, we the incident dilution** (different observable + location: magnet behind a *thin* ST inboard blanket vs our first wall behind 40 cm) |
| breeder | Pb-Li | **FLiBe** (Li₂BeF₄) | different material → cross-check |
| geometry | spherical tokamak (D-ish) | square cross-section torus | different → cross-check |
| analytic anchor | verified OpenMC vs other codes generally | **built independent analytic (Schwartz/anarrima) + showed MC recovers it** | our methodological addition (Bae has no independent analytic NWL) |
| sampler validation | Fig. 3: sampled pdf vs W(θ)sinθ overlay | Tier-1 χ²/dof≈1, C++↔Py parity 6e-16 | same idea, both validate the sampler against the analytic pdf |
| free-stream vs transport | shows pre-scatter Fig. 4 + transport **separately**; no dilution ratio | **explicit free-stream→transport ratio: ~3× dilution at the first wall** | **complementary, not a flaw in Bae**: we tabulate the ratio they don't; magnitude is geometry-dependent (their thin ST inboard preserves −40%, our thick symmetric blanket dilutes more) |
| wall heating/spectrum | not reported | Tier 5 heating split + down-scattered tail | minor delta |
| rate (un-normalized) | flagged as future | Tier 7 (η scaling, +50% A / −50% C) | we did the simple version they deferred |

## Honest contribution statement (use this with Ethan)
1. **No new physics.** Rate/anisotropy/steering = Kulsrud 1982, Ciullo 2016,
   Schwartz 2025; MC neutronics incl. TBR-per-polarization + magnet benefit = Bae 2025.
2. **Independent validation / cross-check.** We reproduce Schwartz analytically
   *and* independently reproduce Bae's TBR-polarization **trend** in a different
   code path (compiled source), different breeder (FLiBe), different geometry
   (square torus). Agreement across these is real corroboration.
3. **A verified, reusable, open tool** (compiled `openmc::Source` with an
   analytic↔MC recovery anchor; handles any B̂ → ready for a real equilibrium /
   Peterson's source).
4. **Thin complementary delta (verified against the full text, now softened):**
   we tabulate the explicit **free-streaming→transport dilution ratio** — the
   geometric steering contrast (±40%) is ~3× the realized transport contrast
   (±13–15%) *at the first wall behind a thick symmetric blanket*. Bae **does**
   show a pre-scatter first-interaction view (Fig. 4) and full-transport results
   *separately*, and even reports a large preserved −40% inboard *magnet* flux —
   but does **not** tabulate the dilution ratio, and their magnet sits behind a
   *thin* ST inboard blanket (less scattering → steering better preserved). So
   this is a **complementary quantification, not a correction to Bae**, and its
   magnitude is geometry-dependent. Do not oversell it.

## Doc/deck fixes (cite Bae) — ALL DONE
- RESULTS_tier6.md: softened "TBR ~independent" → "modestly polarization-dependent
  (±0.6% here; same sign as Bae's ±2.7%)". DONE.
- RESULTS_tier5.md: notes Bae's +68% magnet result; positions the ~3× dilution as
  a geometry-dependent complementary result. DONE.
- slides/spf_talk.tex: literature slide rewritten (removed the "open ground"
  overclaim, cites Bae); Tier-5 dilution, Tier-6 TBR, and Summary slides cite Bae
  and reposition as validation/tooling. DONE (recompiled).
- README.md: added a "Positioning" section + Bae/Schwartz/Kulsrud references. DONE.
