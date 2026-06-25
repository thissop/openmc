# Prior work: Bae et al. 2025 vs. this prototype (honest positioning)

**Bae, Borowiec, Badalassi, Parisi, Diallo, Menard, Khodak, Brown,
"Neutronics analysis of spin-polarized fuel in spherical tokamaks,"
*Nuclear Fusion* 65, 086051 (2025)** (ORNL + PPPL; DOI 10.1088/1741-4326/adf3c6).
This is the closest prior work and it predates and overlaps our Tiers 5–6. **Our
work is independent validation + tooling, not new physics.** Cite Bae up front.

## What Bae et al. did
- **Code:** OpenMC (ENDF/B-7.1), 3D Monte Carlo, MAGIC variance reduction.
- **Source:** 10 M pre-generated samples (position/angle/energy) from physics-based
  emission functions, loaded as an OpenMC source via the **Python API as point
  sources** (a fixed sample bank).
- **Machine:** STAR spherical tokamak — R=4 m, a=2 m, aspect ~2, 800 MW, 5.2 T.
- **Blanket:** **Pb-Li**, 92% ⁶Li-enriched; homogenized layers, 90° reflective sector.
- **Modes:** unpolarized / perpendicular / parallel (anti-aligned).
- **TBR:** unpol **1.082**, perpendicular **1.062** (−1.8%), parallel **1.111** (+2.7%).
  Inboard/outboard split: parallel steers outboard (inboard 0.178 vs 0.186 unpol).
- **Magnets:** parallel cuts inboard magnet flux 40% → **+68% magnet lifetime**
  (linear-fluence assumption); perpendicular +30% flux (worse).
- **Held total neutron rate constant** (explicitly notes the perpendicular/A case
  would give ~+50% reactivity, *not modeled*).
- **Did NOT:** report wall heating/spectrum; isolate free-streaming vs scattering;
  optimize design; couple to plasma/burn-efficiency; do tritium transport.

## Point-by-point vs. our prototype
| aspect | Bae et al. 2025 | ours | verdict |
|---|---|---|---|
| MC neutronics of SPF | OpenMC, sample bank, point source | OpenMC, **compiled C++ `Source` plugin** (on-the-fly, any N, any B̂) | **they're first**; ours = independent reimplementation (different route) |
| TBR per polarization | parallel **+2.7%**, perp −1.8% | B/C **+0.6%**, A −0.4% | **same sign/trend**, smaller magnitude (FLiBe sq-torus vs PbLi ST) → **independent cross-validation, not novel** |
| center-stack benefit | **+68% magnet life**, −40% inboard flux | B/C −15% inboard incident load (w/ scattering) | same physics; **they quantified magnet life, we the incident dilution** |
| breeder | Pb-Li | **FLiBe** (Li₂BeF₄) | different material → cross-check |
| geometry | spherical tokamak (D-ish) | square cross-section torus | different → cross-check |
| analytic anchor | verified OpenMC vs other codes generally | **built independent analytic (Schwartz/anarrima) + showed MC recovers it** | our methodological addition |
| free-stream vs transport | not isolated | **explicit: scattering dilutes steering ~3×** | **our clearest delta** (Bae didn't report it) |
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
4. **Possible thin delta:** the explicit **free-streaming-vs-transport** result —
   the analytic/geometric steering contrast (±40%) is ~3× the realized transport
   contrast (±13–15%); i.e. free-streaming over-predicts the steering benefit.
   Appears absent from Bae; verify against their full text before claiming.

## Doc/deck fixes required (cite Bae)
- RESULTS_tier6.md: soften "TBR ~independent / steering doesn't cost breeding" →
  "modestly polarization-dependent (±0.6% here; same sign as Bae's ±2.7%)". DONE.
- RESULTS_tier5.md: note Bae's +68% magnet result; position the ~3× dilution as
  our delta. DONE.
- slides/spf_talk.tex: literature + Tier-6 + summary slides must cite Bae and
  reposition as validation/tooling. **TODO (deck not yet updated).**
- README.md: add Bae to references / reposition. **TODO.**
