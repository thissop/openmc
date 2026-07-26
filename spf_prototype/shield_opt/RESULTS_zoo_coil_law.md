# Zoo-wide free-streaming coil-load concentration → geometry law

**SCIENCE_STRATEGY.md §5 Step 2.** Cheap (no DAGMC, no OpenMC transport), real *n* from
the QUASR zoo, built to attack the standoff confound that muddies the n=2 through-shield
comparison. Publishable even if the through-shield law never closes.

- Repo commit: `8f0b7e002` (branch `spf-prototype`)
- Compute: `shield_opt/zoo_coil_law.py` → `shield_opt/data/zoo_coil_law.csv` (366 devices, 0 failures, 395 s)
- Analysis: `shield_opt/zoo_coil_law_analyze.py`
- Figures: `figs/zoo_coil_xi_heatmap.png`, `figs/zoo_coil_concentration_scatter.png`

---

## 1. Method (all geometric — no transport)

For every cached QUASR device we load the real coil filaments (`quasr_loader`) and the
(1−ρ²)-weighted plasma fusion source (`QuasrDevice.source_sample`, the same source the
compiled neutron source uses), then compute a **free-streaming, occlusion-free** neutron
load on each coil-filament point:

```
φ(p) = Σ_x  s(x) / |p − x|²        (isotropic point emission, 4π dropped; convex free-streaming limit)
```

This is the solid-angle / line-of-sight neutron load on the conductor — the same
free-streaming approximation the committed first-wall culprit map (`sweep/geometry_peaking.py`)
uses, extended from the conformal wall to the coil filaments as the strategy directed.
Resolution (source 6×16×48, coil 160 samples/curve) is convergence-checked: `point_gini`
and `coil_peak` change <1 % versus a 10×28×96 / grid.

**Two concentration views** (they answer different questions; both reported):

| metric | definition | note |
|---|---|---|
| `point_gini` **(primary)** | Gini over φ at **all** filament points, pooled across coils | non-degenerate; measures spatial hot-spot structure of the coil-side load — the reactor-relevant "are the leak paths few and targetable" question |
| `coil_peak` | max/mean of φ over all filament points | peak conductor load; rank-based ξ is robust to its near-field 1/r² tail |
| `coil_gini` | Gini over **per-coil mean** load (strategy-literal "across coils") | **degenerate**: stellarator-symmetry copies see identical load, so this collapses to the `nc_per_hp` distinct base coils and is exactly 0 for the 61 `nc_per_hp=1` devices. Reported with caveat; not trustworthy. |

---

## 2. Descriptor availability (Task 1)

All 366 cached devices have **every** descriptor (from `catalogue.csv.gz`, `coil_standoff.csv`,
or computed from filaments). Nothing is missing.

| descriptor | source | notes |
|---|---|---|
| coil–plasma standoff | `coil_standoff.csv` (`d_min_over_a`, `d_min_over_R0`) | precomputed for all 366 |
| coil non-planarity | **computed here** from filaments (SVD out-of-plane fraction `s₂/√(s₀²+s₁²)`) | not in catalogue; `nonplanarity` (mean) + `nonplanarity_max` |
| winding curvature | catalogue `max_kappa` (+`max_msc`) | max filament curvature |
| plasma elongation | catalogue `max_elongation`, `mean_elongation` | |
| iota | catalogue `mean_iota` | discrete-ish (many at 0.1) |
| nfp | catalogue `nfp` | integer 1–8 |
| QS type (QA/QH) | catalogue `helicity` (0=QA, ≠0=QH) | 157 QA / 209 QH |

---

## 3. Engineering-relevant subset (Task 2)

Reactor-relevance cuts (`coil_standoff.py` convention: scale QUASR R₀~1 to a 10 m reactor):

- **standoff**: `d_min/R0 × 1000 cm ≥ 101 cm` (the sol5+W0.2+steel3.8+Be2+FLiBe50+shield40
  blanket+shield stack to the coil) → **348 / 366** pass.
- **aspect ratio** in [4, 12] (reactor-plausible stellarator band) → **147 / 366** pass.
- **both** → **140 devices** (61 QA / 79 QH). This is the engineering subset used for the
  headline table.

The standoff cut is weak (almost everything clears 1 m at 10 m reactor scale); the aspect
cut does most of the filtering.

---

## 4. Descriptor → response ξ table (Task 4)

Chatterjee ξ (`descriptor_correlation.chatterjee_xi`). Uncertainty is **jackknife SE** +
**permutation p-value** — *not* a with-replacement bootstrap, which is corrupted here:
resampling rows makes duplicate X-values, and ξ with random X-tie-breaking is spuriously
inflated by ties. Null: ξ ~ N(0, 2/5n), so the 95 % one-sided threshold is
1.645·√(2/5n) ≈ **0.05** (full zoo) / **0.09** (engineering subset).

**Reliability caveat:** the *discrete* descriptors (`iota`, `nfp`, `qs_type`) have heavy
X-ties, which makes ξ itself seed-dependent (jackknife SE 0.2–0.7). Their ξ values are
**not trustworthy** — for QS type we use Mann-Whitney instead (§5). The *continuous*
descriptors (standoff, non-planarity, curvature, elongation; SE ≈ 0.06–0.09) are reliable.

### Engineering subset (n = 140), * = perm p < 0.05

| descriptor | `point_gini` (concentration) | `coil_peak` | `coil_gini` (degenerate) |
|---|---|---|---|
| standoff `d_min/a` | +0.06 ± 0.09 | +0.05 ± 0.09 | +0.02 ± 0.10 |
| non-planarity (mean) | −0.01 ± 0.10 | −0.07 ± 0.10 | +0.10 ± 0.09* |
| non-planarity (max) | +0.09 ± 0.08* | +0.16 ± 0.08* | +0.14 ± 0.08* |
| curvature | +0.02 ± 0.09 | **+0.17 ± 0.08*** | +0.04 ± 0.09 |
| **elongation** | **+0.15 ± 0.09*** | **+0.17 ± 0.08*** | +0.03 ± 0.10 |
| iota *(unreliable)* | +0.01 ± 0.23 | +0.17 ± 0.15* | +0.05 ± 0.34 |
| nfp *(unreliable)* | +0.07 ± 0.28 | +0.01 ± 0.40 | +0.00 ± 0.22 |
| qs_type *(→ §5)* | +0.02 ± 0.27 | +0.13 ± 0.24* | +0.00 ± 0.29 |

Full-zoo (n=366) values are consistent and are printed by the analysis script.

### Reading it

- **Everything is weak.** No descriptor reaches ξ > 0.2 for concentration. This directly
  **confirms the strategy doc's skeptical prediction (§1.2)**: the coil/magnet side is *not*
  a clean geometric projection the way the first wall is (first-wall elongation→concentration
  ξ = 0.96). The coil-side load is smeared over the whole conductor by 3-D geometry.
- **Elongation is the one robust geometry driver of concentration** (ξ ≈ 0.15, p < 0.05) —
  a faint echo of the strong first-wall elongation law. See §5 for its standoff-robustness.
- **`coil_peak` is more predictable than concentration**: curvature, elongation (and iota)
  all ξ ≈ 0.17. Physically sensible — the peak conductor load is set by the closest/curviest
  approach of the winding to the plasma, whereas the *shape* of the load spread (concentration)
  is nearly geometry-independent.

---

## 5. Standoff-confound decomposition (Task 4, the whole point)

Does concentration merely proxy standoff? **No — for concentration, standoff is not even a
driver.**

1. **ξ(standoff → point_gini) = +0.06 ± 0.09, perm p = 0.135 → not significant.** The
   free-streaming coil-load *concentration* is essentially **independent of coil standoff**.
   (Contrast: standoff *does* drive `coil_peak`, ξ = 0.15* in the full zoo — standoff sets
   peak *magnitude* but not concentration *shape*.) The scatter (`figs/…_scatter.png`, left)
   is a structureless cloud, confirming this.

2. **Standoff-stratified ξ(elongation → point_gini):** full 0.15, and within low/mid/high
   standoff tertiles **[0.10, 0.15, 0.16]** — stable across all three. Elongation's effect is
   *not* a standoff artifact. (Other continuous descriptors are unstable across tertiles:
   curvature [−0.13, 0.20, −0.11], non-planarity-max [0.00, 0.19, 0.00] — consistent with noise.)

3. **Standoff-residualized ξ (concentration − smooth standoff trend):** elongation survives at
   ξ = 0.11, perm p = 0.015; non-planarity-max marginal (0.09, p = 0.042); everything else null.

4. **QS type (Mann-Whitney, proper for a binary descriptor):** QH devices are slightly *more*
   concentrated than QA — median point_gini 0.338 (QH) vs 0.313 (QA), **p = 0.020**. A weak but
   significant effect, orthogonal to standoff.

**Conclusion of the decomposition.** The standoff confound that dominates the n=2 through-shield
story **does not dominate the free-streaming coil concentration** — concentration is nearly
standoff-independent. The only geometry descriptor with a genuine, standoff-robust link to
concentration is **plasma elongation (ξ ≈ 0.15)**, a weak coil-side echo of the ξ = 0.96
first-wall law. QH is marginally more concentrated than QA. This is an honest, real-*n*,
transport-free result and it lands as a **weak / partly-null law**, exactly the outcome the
strategy doc flagged as likely for the magnet side.

---

## 6. Limitations (honest)

- **Free-streaming proxy, not transport.** φ ∝ Σ s/r² with no scattering, no occlusion, no
  blanket/shield attenuation. This is the coil-side analog of the first-wall culprit map and is
  explicitly a **geometric proxy** for where neutrons would deposit *if* they streamed freely.
  Through-shield magnet concentration (Step 5) is smeared further and is expected to be *even
  less* geometry-predictable.
- **`coil_gini` (across-coil) is degenerate** under stellarator symmetry (0 for `nc_per_hp=1`);
  do not use it. `point_gini` (pooled filament points) is the trustworthy concentration metric.
- **`coil_peak` has a near-field 1/r² tail** (max 54.7 for the tightest coils) — magnitudes are
  distorted by the point-source singularity, but the rank-based ξ is unaffected.
- **Discrete-descriptor ξ (iota, nfp, qs_type) is unreliable** (heavy X-ties → seed-dependent);
  use the Mann-Whitney / stratified views for those.
- **QUASR is QS-only (QA/QH), R₀ ≈ 1 native.** No QI; reactor scaling is a uniform ×1000.
  Standoff cut is therefore aspect-independent by construction (`d_min/R0`).
- **Effect sizes are small.** The headline elongation ξ ≈ 0.15 is real (p < 0.05, standoff-robust)
  but weak; this is a trend, not a tight law. Do not overclaim.
