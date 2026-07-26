# Note: free-streaming first-wall attribution — QA/QH and the QUASR zoo

*2026-07-22. Local analytic (free-streaming point-kernel) study; no Monte-Carlo transport, no
magnets. All numbers reproducible from `shield_opt/` (culprit_59509 machinery + concentration.py +
emissivity.py + descriptor_correlation.py). Figures in `shield_opt/data/processed/figs/`.*

## What this is

First-wall neutron-load **attribution** (which plasma regions drive the peak load) and **load
concentration** for two canonical anchors (precise-QA, precise-QH) and a **24-device QUASR zoo**
(nfp 1–7, aspect 1.5–14.6). Because first-wall loading is dominated by uncollided source neutrons
streaming to a nearly-transparent-plasma-facing wall, the analytic point-kernel is valid here and
runs in seconds/device — so this regime is appropriate for the large QUASR sample (unlike the
magnet-side work, which needs transport + reactor-relevant devices).

Emissivity is a **pluggable weight** `r(rho)`: `uniform` (|√g| only — isolates geometry+polarization)
vs `bosch_hale` (core-peaked DT reactivity, profiles per Miralles-Dolz 2026: n∝(1−s⁵), T∝(1−s),
Bosch–Hale ⟨σv⟩). Concentration measured with Gini, participation ratio, and the culprit-map
participation (fraction of plasma volume supplying 50/90% of the peak load).

## Finding 1 — realistic emissivity nearly erases the QA–QH peaking gap

| first-wall peaking (max/mean) | uniform | Bosch–Hale |
|---|---|---|
| precise-QA | 1.38 | 1.56 |
| precise-QH | 1.91 | 1.58 |
| **QA–QH gap** | **0.54** | **≈0.02** |

Under uniform emissivity the strongly-shaped QH concentrates first-wall load far more than the
gentle QA. Weighting by a realistic **core-peaked** reactivity collapses the difference — because
uniform emissivity over-weights the shaped *edge*, while real fusion happens in the rounder *core*.

**Implication (honest, citable):** device-to-device first-wall peaking differences reported under a
uniform-emissivity assumption may be **substantially overstated**. The emissivity profile is a
first-class modelling knob, not a detail — which is the whole point of the pluggable-emissivity
decomposition (`geometry × polarization × emissivity`). This is exactly the kind of result the
analytic tool is *for*: it cleanly separates the three factors, where full transport would mix them.

## Finding 2 — iota drives the *peak*; elongation drives the *distribution*; nfp is weak

Full descriptor set (24 devices, Bosch–Hale, iota + QA/QH type + qs_error pulled from the live
QUASR records). Chatterjee ξ, **top predictor per concentration metric**:

| response | top predictor | ξ |
|---|---|---|
| **peaking (max/mean), unpol** | **iota** | **0.34** |
| NWL Gini (unpol) | **elongation** | 0.28 |
| culprit Gini | qs_error (QS quality) | 0.26 |
| NWL Gini, A-mode | **iota** | 0.35 |
| field periods (nfp) — any | (weak) | ≤0.12 |

The two convergent predictors flagged independently by **Miralles-Dolz §IV** (neutronics
sensitivity) and **QUASR-PCA** (geometry) both show up — and they split cleanly:

- **iota controls how hot the single peak gets** (ξ=0.34 for peaking; 0.35 for A-mode concentration)
  — matching Miralles-Dolz's finding that iota drives the deeper/peaked response.
- **elongation controls the overall distributional concentration** (Gini, ξ=0.28) — matching
  QUASR-PCA's finding that elongation is the dominant geometric axis.
- **QS quality (qs_error) also predicts concentration** (ξ≈0.26) — better-quasisymmetric configs
  concentrate differently; worth a second look.
- **nfp is weak throughout** (ξ≤0.12).

Honest caveat: ξ≈0.3 is a **real but moderate** signal — leading predictors, not determinants;
report with the scatter shown. Heatmap: `figs/zoo_descriptor_xi.png`.

## The polarization × geometry effect

SPF A-mode's effect on the culprit is **device- and emissivity-dependent**: under uniform emissivity
A *diffuses* the QA culprit (Gini 0.50→0.34) but *concentrates* the QH culprit (0.60→0.66); under
Bosch–Hale it slightly *raises* peaking on both (opposite to uniform). Some zoo devices show strong
A-mode peaking amplification (e.g. the nfp=7 device: 1.75→3.02). The polarized emission redistributes
source contributions differently depending on the local field geometry — a genuine SPF×geometry
coupling not previously shown.

## Caveats

- **Free-streaming, first-wall only** — valid here; the magnet-side attribution needs transport +
  the adjoint/reciprocity map (running separately on the cluster).
- **Geometry descriptors only so far** (nfp, aspect, elongation) — iota + QS-type incoming.
- **Moderate correlations** (ξ≤0.28) — a trend, not a law; report with the scatter shown.
- QUASR coils/scale are not reactor-engineering-relevant — irrelevant here (plasma-shape → first
  wall only), but this is why the magnet story stays on precise-QA/QH.

## Where it fits

This is the **first-wall / attribution spine of the science paper**, with real data: culprit maps +
concentration index + a 24-device trend + the geometry/polarization/emissivity decomposition showing
the emissivity assumption *changes the conclusion*. The magnet/shield/adjoint spine is the companion
(transport-dependent, cluster).

Figures: `bh_culprit_voxels_3d.png` (3D volumetric culprit, Bosch–Hale), `qa_qh_culprit_compare.png`
(uniform QA vs QH). Viewable artifact: the 2026-07-22 results page.
