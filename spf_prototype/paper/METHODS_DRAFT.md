# Methods (draft)

*Draft for the geometry-first stellarator-neutronics paper. Covers the adjoint
magnet-shielding pipeline (§4 of `PAPER_DIRECTION.md`). Written from committed
results; every transport number cited here was produced on Ginsburg (x86) with the
DAGMC-enabled OpenMC 0.15.x build and frozen. Citation keys match
`paper/manuscript/refs.bib`.*

---

## 2.x Overview

The methodological center of this work is an **adjoint coil-kerma importance map**:
an adjoint transport solve whose source is the magnet response function and whose
solution, read over the plasma, is a field that attributes coil neutron loading to
its phase-space origins *through* the ~1 m of intervening blanket and shield. Unlike
the FW-CADIS / CADIS use of the adjoint flux as a variance-reduction weight
[differentiate from Miralles-Dolz, Thea random-ray FW-CADIS], we use the adjoint flux
directly as a coil→plasma *attribution* map and feed it to a closed-loop
breeder-for-shield optimizer. The pipeline has six stages, described in turn:
(i) geometry and plasma-source model; (ii) the adjoint importance map; (iii) the
contributon and its reactivity weighting; (iv) projection to a shield control grid and
the closed-loop optimizer; (v) forward↔adjoint reciprocity validation; and
(vi) an angular-sensitivity (P0 vs P1) check. All results below are for the
Landreman–Paul quasi-helically symmetric (QH) reactor with the Wiedman QH coil set;
the quasi-axisymmetric (QA) device is in progress.

## 2.x.i Geometry and plasma-source model

The device geometry is a watertight DAGMC model [`dagmc`] built from the real QUASR
equilibrium [`quasr`]: a conformal radial build offset outward from the last closed
flux surface (LCFS) along the poloidal normal into nested shell solids
(plasma-facing wall / structural steel / Be multiplier / FLiBe breeder / WC shield /
coil), plus the QH coil filaments. The coil cell is presently a lumped homogenized
region; per-coil spatially-resolved damage tallies are future work (see Limitations).

The plasma neutron source density is

```
S(r) = |sqrt(g)(r)| * r(rho(r)),      r(rho) = n(rho)^2 <sigma_v>_DT(T(rho))
```

where `|sqrt(g)|` is the VMEC flux-surface volume element carried on the equilibrium
grid, and `r(rho)` is a normalized DT reaction-rate profile. The reactivity uses the
Bosch–Hale (1992) DT `<sigma_v>` [`schwartz2025` cross-checks the same reaction] with
the Miralles-Dolz core-peaked profiles `n(s) = 4.8e20 (1 - s^5) m^-3`,
`T(s) = 11.5 (1 - s) keV`, `s = rho^2`. **Core-peaked Bosch–Hale is the default and
physical choice.** A uniform-emissivity option (`r = 1`, weighting only `|sqrt(g)|`)
is retained for geometry isolation and for matching uniform-source forward runs, but
it over-weights the strongly-shaped plasma edge and ~3x over-states the localizable
optimizer benefit on QH; the code emits a one-time warning whenever uniform emissivity
is used so it is never applied silently.

## 2.x.ii Adjoint through-shield coil-importance map

The importance map is the adjoint scalar flux `psi_dagger(r)` obtained by solving the
adjoint transport equation with a source equal to the **magnet response localized in
the coil**. Physically, `psi_dagger(r)` is the expected contribution to the chosen
coil response of one source neutron born at plasma point `r`, integrated over all paths
through the blanket and shield. We solve it with OpenMC's random-ray adjoint solver
(the Transient/Random Ray Method, TRRM) [same solver family as `peterson_tokamaksource`
uses forward], using the mainline localized-adjoint-source capability
(`FlatSourceDomain::set_local_adjoint_sources`); no feature-branch rebuild is required.

The workflow is: build the continuous-energy model with a plasma-mesh flux tally;
`convert_to_multigroup` on a 13-group structure with fast detail above 0.1 MeV;
`convert_to_random_ray`; then set `random_ray['adjoint'] = True` and
`random_ray['adjoint_source']` to the localized coil response. In adjoint mode the
plasma-mesh "flux" tally scores the adjoint flux, i.e. the importance map itself.

**Adjoint source geometry.** A QH coil is a non-planar curved winding; a point or box
at the coil-cell centroid is unreliable (the volume centroid of a curved coil can fall
off the conductor — verified for cell 15, whose centroid pointed away from the coil),
and DAGMC cells are not enumerable before the run so a per-cell domain constraint is
not available. We therefore localize the response *on the coil guide curve*: parse the
MAKEGRID coil filament, align it to the DAGMC frame (m→cm; the scaled filaments match
the DAGMC coil centroids bijectively to ~35–50 cm, one filament per coil cell),
subsample ~40 points along the winding, and emit a list of isotropic point sources that
random ray sums into one filamentary adjoint source lying on the conductor path.

**Response energy dependence.** Two options are supported: `flat` (equal weight to
groups in a fast band, default 0.1–15 MeV, capturing the fast neutrons that dominate Cu
heating and damage) and `kerma` (weighted by the Cu KERMA coefficient per group, the
physically correct magnet nuclear-heating response). Adjoint sources in random-ray mode
require a discrete multigroup energy distribution (a continuous distribution aborts the
solve), so both responses are evaluated at group midpoints. The committed QH coil-15
closed-loop result uses the `flat` fast-band response.

**Key property — reactivity independence.** `psi_dagger(r)` depends only on geometry and
the coil response; it contains no plasma-source information. Consequently, switching the
plasma reactivity profile (uniform ↔ Bosch–Hale, or any user scenario) is a *free
re-weight* of the attribution and requires no new transport solve. This is the property
that makes the emissivity study and the physical per-coil ranking cheap.

A structure gate hard-asserts that the map is non-empty (`total > 0`), localized
(nonzero-voxel fraction `< 0.98`), structured (`peak/mean > 1.5`), and converged
(median relative error `< 0.5`) before it is saved for downstream use.

## 2.x.iii Contributon and reactivity weighting

The bare adjoint flux peaks in the optically thin central vacuum, which is physically
irrelevant because no source neutrons are born there. The attribution map is the
**contributon**, the product of the plasma source and the importance:

```
C(r) = S(r) * psi_dagger(r).
```

`C(r)` answers "which plasma region drives this coil's load," restricted to where
coil-bound neutrons are actually born — the through-shield analog of the free-streaming
first-wall culprit map. Because `psi_dagger` is reactivity-independent, `C(r)` for any
emissivity model is obtained by re-weighting the single committed adjoint solve. On QH,
switching from uniform to Bosch–Hale emissivity moves the contributon toward the core
and compresses the localizable benefit (see §2.x.iv).

## 2.x.iv Placement and closed-loop breeder-for-shield optimizer

**Projection to the control grid.** The optimizer's degree of freedom is a smooth added
shield-thickness field `delta(theta, phi)` on the (poloidal, toroidal) shield surface.
Each plasma voxel is assigned a toroidal angle `phi = atan2(y, x)` and a poloidal angle
`theta = atan2(z, R - R0)` about the magnetic-axis major radius `R0` (auto-estimated as
the contributon-weighted mean cylindrical radius), and its contributon weight is
Gaussian-splatted (periodic in both angles, default width 15°) onto the control grid to
form a placement-priority field `P(theta, phi)`. The peak of `P` is the data-driven
target `(theta0, phi0)`. As a sanity check, `phi0` for the QH coil-15 contributon
coincides with the coil's own toroidal angle (`|phi0 - phi_coil| < 45°`), i.e. the map
points shield at the coil it was built for.

**Fixed-envelope breeder↔shield trade.** At each `(theta, phi)` the shield is thickened
by `delta` and the breeder thinned by the same `delta`, holding the outer envelope
fixed. The coil dose along a sightline follows a near-exponential attenuation surrogate,

```
D(delta) = D0 * exp(-k * delta),     k = 1/lambda_shield - 1/lambda_breeder,
```

with `k > 0` because the dense WC/steel shield removes fast neutrons in fewer cm than
FLiBe, so trading breeder for shield over a hot coil always reduces its dose until the
tritium-breeding ratio (TBR) floor bites. `D0` is the adjoint-informed baseline
coil-dose field (the placement priority `P`); `lambda_shield`, `lambda_breeder` are
effective fast-neutron removal mean free paths (literature-scale defaults 8 cm / 17 cm,
**to be calibrated to OpenMC uniform-thickness scans** — see Limitations).

**Objective and optimizer.** The design vector is the coil standoff scalar(s) plus the
Fourier coefficients of `delta(theta, phi)`. The objective to minimize is a
Kreisselmeier–Steinhauser soft-max of the coil-dose field (the hot tail, not a single
noisy cell) plus a one-sided TBR-floor penalty and light smoothness/material/standoff
regularizers:

```
J = KS_softmax( D(delta) )
    + penalty * max(0, TBR_floor - TBR(delta))
    + lam_smooth * smoothness(delta) + lam_material * volume(delta) + lam_standoff * standoff,
```

with the TBR linearized about the baseline,
`Delta_TBR = -sum_cells (dTBR/d(breeder_cm)) * delta * area_frac`. Because MC coil-dose
gradients are too noisy to descend cell-by-cell, we optimize with Bayesian optimization
(Gaussian process + Expected Improvement; scikit-optimize when available, otherwise a
self-contained GP/EI fallback): most evaluations use the cheap noise-free surrogate,
and an MC-verify hook occasionally corrects a promising candidate. (In the committed
offline run the MC-verify hook defaults to the surrogate; a cluster hook is future work.)

**Result.** With the adjoint-informed baseline field, the optimizer localizes shield on
the coil-facing sightlines and cuts the peak coil dose by **~19%** on QH coil-15 under
realistic Bosch–Hale emissivity, at a small breeder sacrifice and with TBR held above
floor. A no-attribution uniform baseline (flat `D0`, no spatial signal) has nothing to
target, spreads shield thin, and achieves only **~-1%**. We report **19% as the honest
number.** Uniform emissivity would inflate the adjoint-informed reduction to ~55%
because it over-weights the shaped plasma edge where the contributon is most localizable;
that figure is non-physical and is reported only to bound the emissivity sensitivity.

## 2.x.v Forward↔adjoint reciprocity validation

Reciprocity predicts that the adjoint contributon integral for a coil equals the coil's
neutron load from the full plasma source computed by an independent forward transport
run:

```
INT S(r) psi_dagger(r) dr   ~   forward per-coil coil flux.
```

We evaluated the left side from the committed per-coil adjoint maps (kerma sweep) with
the uniform-emissivity `S(r)` that matches the uniform-source forward run, and correlated
it against the committed forward per-coil flux (matched geometry `w35f`, matched response
band) across n = 20 coils. Because a longer winding spreads the fixed 40-point source
over more plasma, the raw integral over-weights long coils relative to the forward flux
density; dividing by coil arc length removes that geometric bias. The rank correlation is
**Spearman +0.76**, confirming that the adjoint attribution reproduces the forward
per-coil loading order. **Honest caveat:** coils 15 and 28 are per-coil outliers, most
likely from filament source-matching on highly non-planar windings; the overall ranking
holds but individual coils can be misranked. The *direction* of attribution (that
`phi0` points at the coil) was validated separately in §2.x.iv; reciprocity here
validates the *magnitude/ranking*.

## 2.x.vi Angular sensitivity (P0 vs P1)

The random-ray adjoint can run scalar (P0, transport-corrected) or with explicit
first-order anisotropic scatter (P1), the latter requiring the l=1 scatter moment
(`NU_SCATTER_FMU` → `Sigma_s1`) in the exported multigroup library; we force
`legendre_order = 1` via a targeted patch to the MGXS library for the P1 path. Both
orders use an untransport-corrected total cross section so the comparison isolates the
explicit angular effect (and because the P0 transport correction drives near-void
materials to negative total XS, which random ray rejects).

On a controlled deep-penetration slab (test M6), the P1/P0 importance ratio grows with
depth into scattering shield — up to ~1.7x at 5 mfp in the scattering bulk — but stays
within ~3% across clean gaps. On the **real QH coil** (test M7), comparing the P0 and P1
adjoint importance maps for the same coil, the median P1/P0 ratio is **~1** and the
shield-placement target is **unchanged** (`phi0` identical for P0 and P1). We conclude
that **a scalar adjoint suffices for coil placement**: the plasma-side attribution the
optimizer actually uses lives in the optically thin region where angular effects are
small, even though angular treatment matters deep inside a scattering slab. We report
this honestly and do not over-claim angular insensitivity in general.

---

## 2.x.z Assumptions and limitations

- **Emissivity model.** Bosch–Hale core-peaked reactivity is the default and physical
  choice; uniform emissivity is retained only for geometry isolation and for matching
  uniform-source forward runs, and it ~3x over-states the localizable benefit (55% vs
  the honest 19% on QH coil-15). We report the realistic number as headline.
- **Single QS device demonstrated.** The full adjoint→optimizer→validation chain is
  demonstrated on the QH reactor only; QA is in progress. The raw QA-vs-QH magnet
  comparison is additionally confounded by coil standoff (the QA coil set sits ~1.9x
  further from the plasma, ~148 cm extra), not by quasisymmetry; a fair comparison
  requires standoff-normalized or same-standoff coil sets.
- **Per-coil outliers.** Reciprocity is Spearman +0.76 over 20 coils, but coils 15 and
  28 are outliers, likely from filament source-matching on strongly non-planar windings.
  Trends are solid; individual coil magnitudes carry that caveat. QA statistics are
  additionally loose (relerr ~23–25%; exact factors ±15%) while QH is tight (~4–5%).
- **Angular treatment is P1-only.** The angular check covers P0 vs P1; higher moments are
  not tested. The conclusion "scalar suffices for coil placement" is specific to the
  plasma-side attribution and does not extend to deep in-shield fields, where M6 shows
  angular growth to ~1.7x at 5 mfp.
- **Filament source approximation.** The adjoint response is a ~40-point filamentary
  source on the coil guide curve, not a volumetric coil-cell response; this is the source
  of the non-planar-winding matching error and the per-coil outliers.
- **Surrogate attenuation objective is uncalibrated.** The closed-loop objective uses an
  exponential-attenuation surrogate with literature-scale removal MFPs (8 cm shield,
  17 cm breeder) and a linearized TBR; these must be calibrated to OpenMC
  uniform-thickness scans before the absolute dose-reduction and TBR numbers are
  quantitatively trusted. The adjoint importance map itself is a full transport solve;
  only the optimizer's fast inner objective is surrogate, with an MC-verify hook that is
  presently a stub. The 19% figure is the surrogate-model peak reduction given the
  adjoint-informed vs uniform baseline *shape*.
- **Lumped coil, placeholder materials.** The coil is a single homogenized cell with a
  generic Fe-Cr-Ni block and bare-WC shield; there is no ReBCO/Cu/steel winding-pack or
  HTS-vs-casing split, so absolute magnet-damage/lifetime numbers are not yet physically
  meaningful. Per-coil unstructured-mesh damage tallies are future work.
- **Plasma importance mesh envelope** is presently hardcoded to the QH bore and flagged
  for per-device confirmation.

## Prior art and differentiation (for the Introduction; noted here for the reviewer)

Parametric stellarator radial-build sweeps (ParaStell, Davis et al.) are manual and
per-device with no coil→plasma attribution [contrast]. FW-CADIS/CADIS use the adjoint
flux as a variance-reduction weight, and random-ray FW-CADIS (Miralles-Dolz, Thea) for
weight-window generation, not as a closed-loop shield-design attribution map [contrast].
Bae et al. 2025 [`bae2025`] published a precomputed static polarized point-cloud source
for one fixed tokamak; the SPF lever here (supplementary, §6) is a native, general,
on-the-fly stellarator source and is explicitly *not* claimed as the novelty. The novel
contribution is the adjoint coil-kerma importance map plus the closed-loop
breeder-for-shield optimizer on a 3-D stellarator.

---

## Supplementary: the SPF lever (for §6 cross-reference)

Where the paper layers spin polarization on top of the geometry/adjoint spine, the native
source samples the Schwartz differential cross section [`schwartz2025`]

```
dsigma/dOmega = (sigma_0 / 2pi) [ (3/4) a sin^2(theta)
                                  + ((2/3) b + (1/3) c) (1/4 + (3/4) cos^2(theta)) ],
```

with `theta` measured from the local field direction, separating the angular sampling
(weights `w_perp = 3a/4`, `w_par = 2b/3 + c/3`) from the total-rate factor
`eta = a + 2b/3 + c/3`. On the magnet side (identical radial build, each device compared
to its own unpolarized case), parallel-emission (B/C) polarization flattens QA coil
peaking ~17% and lowers peak flux ~15% while barely moving QH (~3%), i.e. a beneficial,
strongly device-dependent magnet-protection lever; A-mode (perpendicular) does the
opposite on QA (+10%). B≡C on both devices confirms the sampler (shared angular shape,
differing only in total rate, which normalizes out of a unit-strength direction-sampled
run). These are within-device results and are independent of the QA-vs-QH standoff
confound.
