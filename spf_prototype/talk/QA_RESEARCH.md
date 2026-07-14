# Ethan Peterson meeting — reference definitions + anticipated Q&A

*Prep for the SPF `StellaratorSource` review with Ethan Peterson (TokamakSource / PR #3999 author).
Grounded in the actual code under `spf_prototype/`. File references are absolute-relative to that dir.
Anything not verbatim in the repo is marked **[UNCERTAIN]**.*

---

# PART A — precise reference definitions (verbatim-accurate for slides)

## A1. `PF_geo` — geometry-only isotropic-flux peaking proxy
Source: `python/predictor_test.py:geometry()` (lines 69–80).

**Formula.** For wall points `w_k` (conformal wall = LCFS offset outward by `0.30·a` along the
poloidal normal) and source voxels `r_v` (the VMEC fluxmap nodes), with per-voxel weight
`w_v = |sqrt(g)|_v / Σ|sqrt(g)|` (normalized volumetric reactivity support):

```
G(w) = Σ_v  w_v / |w − r_v|²                 # transport-free isotropic 1/d² flux proxy
PF_geo = max_k G(w_k) / Σ_k G(w_k) dA_k       # dA_k = area-fraction weights, Σ dA_k = 1
```

i.e. `PF_geo = max(G) / area-weighted-mean(G)`.

**Plain English.** The neutron wall load a *point isotropic* source cloud would deposit, using only
inverse-square geometry — no scattering, no angular physics, no polarization. Its peaking factor
(hot-spot / mean) is the **peaking headroom the geometry alone bakes in**: how non-uniform the wall
load is before any transport or steering. It is "transport-free headroom" because SPF can only flatten
peaking that exists to be flattened; `PF_geo` measures that budget from geometry in milliseconds.
(Best single predictor of SPF benefit in the n=15 set, but only borderline: Spearman ρ=+0.525.)

## A2. `nematic_Q` — LCFS wall-normal director alignment
Source: `python/predictor_test.py:geometry()` (lines 84–86).

**Formula.** For the LCFS poloidal *outward* normals `(nR, nZ)` (from
`stellarator_geometry.poloidal_outward_normals`), let `ψ = atan2(nZ, nR)` be each normal's angle in
the (R,Z) plane. Then

```
Q = | mean( exp(2i·ψ) ) |          # 2i because a director is a headless axis (period π)
```

**Plain English.** The **nematic order parameter** of the wall-normal *director* field. The `2ψ`
(not `ψ`) makes it orientation-blind: `n` and `−n` count the same, as for liquid-crystal directors.
`Q=0` ⇒ the boundary normals point every which way (isotropic orientation, strongly/asymmetrically
shaped cross-section); `Q=1` ⇒ all normals aligned (a slab/ellipse-like boundary). It is a **pure
plasma-shape metric** — computed from the LCFS alone, no wall standoff, no transport. In the predictor
test it is the representative "does global plasma shape forecast steering?" descriptor (answer: weakly,
ρ=+0.396, n.s.).

## A3. `ε_eff` (epsilon_eff) — helical excursion / source-to-wall standoff
**This is a coined name for the anarrima-companion controlling parameter.** Exact definition:
`docs/analytic upgraded.tex`, Eq. `\label{eq:epseff}` (line 416); computed in
`python/quasr_geom.py:eps_eff()` (lines 116–133).

**Formula (paper form, Eq. eq:epseff):**
```
ε_eff ~ max_φ sqrt(dR² + dZ²) / min_φ |Δ₀|
```
= (non-axisymmetric boundary excursion) / (smallest source-to-wall standoff). `Δ₀` is the
source-to-wall separation; the closest approach controls the bound. **Truncation error of the analytic
series at order p scales as `ε_eff^(p+1)`.** For a wall conformal to the plasma, `min|Δ₀|` is the
plasma–wall gap and `ε_eff` = ratio of helical excursion to that gap. It is a property of the
**boundary geometry, not of |B| quasisymmetry** (a QA device need not be geometrically axisymmetric).

**Formula (code form, `quasr_geom.eps_eff`, per constant-θ filament loop):**
```
e_geom      = max_θ ( max_φ |loop − <loop>| ) / a         # intrinsic, aspect-invariant excursion
ε_eff(gap)  = max_θ ( excursion / standoff )              # standoff to a square wall gap·a outside bbox
   with standoff so = min(p−rin, rout−p, ztop−zc, ztop+zc),  gap = 0.25 default
```

**Regime thresholds (from `docs/analytic upgraded.tex` §Convergence, lines 411–427):**
- `ε_eff ≲ 0.15` — second-order series already ≤1.5% (Schwartz's uniform-NWL residual level); QA
  regime, `g⁰+g¹+g²` sufficient.
- `ε_eff < 1` — series still convergent (need higher orders `g³` or the hybrid quadrature as it grows).
- `ε_eff ≳ 2` — **series breaks down**; characteristic spurious double-peak in the inboard valley;
  must fall back to differentiable hybrid quadrature. (`precise_QA ≈ 2.3` is MC-only.)

**Per-device / per-config values found in the repo:**
- **59509** (the scatter/validation device, nfp=3, aspect 6.7, R0=1, a≈0.175): `e_geom ≈ 0.29`;
  at the square validation wall (gap 0.4·a) **ε_eff ≈ 0.43**. (`RESULTS_quasr.md` lines 18–20.)
- Toy-QA-low patch: `ε_eff ≈ 0.08`. Published QAs at a tight wall: precise_QA/reactor_QA **2.3**,
  ESTELL 1.0, NCSX 1.3, ARIES-CS 1.6, WISTELL-A 1.6 (`RESULTS_quasr.md` lines 8–10).

**[UNCERTAIN] on your "two thresholds ~1/2 and ~2":** the repo explicitly states the **~2** grazing/
breakdown threshold (precise_QA 2.3) and the `<1` series-convergence bound and the `≲0.15`
second-order-accuracy bound. I did **not** find a verbatim "~1/2 (series)" threshold; the closest
lower crossover stated in the files is `0.15`. Recommend citing `≲0.15` (accurate) and `≳2`
(breakdown), or verify the `1/2` number before it goes on a slide.

## A4. How the 15 QUASR devices were selected
The device set is layered; the *provenance* is a stratified mine of the full QUASR catalogue plus a
"gentle-QA" shortlist, then whatever survived the VMEC + conformal-wall pipeline.

**Pool construction (`sweep/select_candidate_pool.py`):**
- Start from the full QUASR catalogue: **371,701 devices** (`data/quasr/catalogue.csv.gz`).
- Can't compute the field-coherence metric `S_φ` for all 371k, so **stratify on the catalogue columns
  that drive it**: symmetry class (QA=helicity 0 / QH), `nfp`, quasisymmetry-error quintile, aspect
  ratio. Sample **6 per (class × nfp × qs_error-quintile) stratum** (seed 20260707) → 354-device
  `candidate_pool.csv`. Explicitly chosen to **span** the coherence range, not cluster.
- (`compute_sphi_pool.py` then computes `S_φ`; `stratified_select.py` down-selects ~120 to a **flat
  `S_φ` histogram** per class — the anti-bias step, since QA piles up near `S_φ≈1`.)

**The 15 transport devices (`predictor_test.py:DEV`, PROVENANCE.md, RESUME.md):**
Ginsburg VMEC batch (`vmec_batch30`, job 8885457): **30 QUASR devices → 22 converged** fixed-boundary
VMEC fluxmaps (QH devices self-intersect in fixed-boundary VMEC and fail — expected). Chained conformal
batch (`conformal_batch30`, job 8885458) built conformal walls + free-streaming NWL maps → **15
devices** with usable maps, **plus the original 5** from the earlier conformal transport run
(1960314, 59509, 803097, 886079, 932746).

**Span of the 15** (from `candidate_pool.csv`; 5 originals not in that file):
- **Class:** ~13 QA + 2 QH (1642553, 1854549). QH under-represented because fixed-boundary VMEC
  rejects strongly-shaped QH.
- **Aspect ratio:** **3.33 → 24** (883496=3.33, 11357=5, 24285/59509≈6.7, 1642553=8, 180790/66633/
  262171/1090019=10–20, 1854549=24). Wide.
- **nfp:** 2–5. **qs_error (log10):** ≈ −2.1 → −4.1 (i.e. spanning modest to good quasisymmetry).
- **Shaping (`shaping_amp` in predictor_n15.csv):** 0.16 → 1.00 — from near-round to strongly shaped.

Devices: `1052272 1090019 11357 1642553 180790 1854549 1960314 24285 262171 59509 66633 803097
883496 886079 932746`.

---

# PART B — anticipated Q&A

## B1. Lineage & design (vs TokamakSource)

**Q1. Why a new `StellaratorSource` class instead of extending `TokamakSource`?**
- TokamakSource generates positions from an **analytic parametric equilibrium** (Miller R/Z with
  elongation/triangularity/Shafranov, per-flux-surface energy spectra). A stellarator has **no such
  closed form** — the plasma boundary and √g are fully 3-D `(ρ,θ,ζ)`.
- `StellaratorSource` samples positions from a **real gridded equilibrium** (DESC/VMEC `spf_fluxmap_v1`
  file) via a rejection-free CDF cascade; there is no parametric surface to reuse.
- The field b̂ must come from the **same 3-D grid** (BR,Bφ,BZ) so position and polarization axis are
  self-consistent — not derivable from Miller parameters.
- They share the *polarized-direction* machinery (see Q3), not the *position* machinery. Two position
  models, one direction model.

**Q2. What did you reuse vs rewrite from PR #3999?**
- **Reused verbatim** (copied, not re-derived) into both `native_spf/source.cpp` (the SPF-patched
  Tokamak) and `native_spf/stellarator_source.cpp`: the SPF direction sampler
  `sample_polarized_direction`, the per-shape inverse-CDF roots `spf_costheta_perp/par`, and the
  `(a,b,c)→w_perp/w_par/P_perp` derivation block. Diff-confirmed identical code.
- **Reused pattern** from TokamakSource: `public openmc::Source` subclass; constructor takes a
  `pugi::xml_node`; `SourceSite sample(uint64_t*) const override`; Python `SourceBase` subclass with
  `populate_xml_element`/`from_xml_element`; `distribution_from_xml` for energy/time; registration via
  `Source::create` dispatch branch.
- **Rewrote** (stellarator-specific, new translation unit): the `spf_fluxmap_v1` reader, the
  marginal→conditional CDF cascade `build_cdf_cascade`, trilinear interp with θ/ζ period wrap, and the
  b̂-from-grid path.

**Q3. Could the polarized-direction hook be upstreamed straight into TokamakSource?**
- Yes — that is exactly what `native_spf/source.cpp` does: it adds an optional `<polarization>` node
  to TokamakSource (absent ⇒ unchanged isotropic behavior, backward compatible), plus a `field_model`
  (`toroidal` default, or `pitched` with a q-profile) to build b̂ = φ̂ + s·λ·p̂.
- The clean upstreaming story: a small self-contained `spf_sampler` header (the `w_perp/w_par/P_perp`
  algebra + two inverse-CDF roots + Gram-Schmidt rotation) that any `Source` can call with a local b̂.
  TokamakSource supplies b̂ from Miller; StellaratorSource supplies b̂ from the grid.

**Q4. Why port Python → C++ line-for-line rather than just call the Python?**
- The source `sample()` runs **inside the OpenMC transport loop, per history, multithreaded** — it must
  be native C++ (no Python GIL, no per-call overhead).
- The Python reference (`ConformalStellaratorSampler`) is the **validated oracle**; the line-for-line
  port (same cumsum axis order, same `searchsorted(side="right")`, same interp term order, same period
  wraps, same RNG draw order) lets a **bit-for-bit parity test** (C++ vs Python on an identical seed
  sequence) catch any transcription bug. Correctness is inherited, then *proven* preserved.

**Q5. Why is the fluxmap a bespoke `spf_fluxmap_v1` binary instead of HDF5?**
- Deliberately **HDF5-free** so the C++ reader is a plain `std::ifstream` — no new build dependency
  (mirrors `spf_fieldmap_v1` style). `.meta` text header (`magic/nR/ntheta/nzeta/nfp/sign_sqrtg`) +
  raw little-endian f8 `.bin` payload in fixed order (`rho,theta,zeta, sqrtg,R,Z,phi,BR,Bphi,BZ`).
- Producer-agnostic: a real DESC/VMEC solve or a synthetic test map writes the identical format;
  round-trips bit-exactly (tested).

## B2. Correctness of the sampler

**Q6. How do you know the direction distribution is exact (not just approximately right)?**
- **Stratified exact sampling, not rejection-into-histogram:** with prob `P_perp=a/η` draw the `sin²θ`
  shape, else the `¼+¾cos²θ` shape, each via a **closed-form inverse-CDF** (exact cubic roots). No
  binning approximation enters the sampler.
- Tier-1 gates: χ² of sampled `cosθ` vs the analytic marginal `pdf_costheta` for {non-pol, A, B, C,
  mixed}; **inverse-CDF vs rejection two-sample KS** (independent method, same distribution) for each
  shape; φ-uniformity; rotation isotropy under a random b̂.

**Q7. What's the inverse-CDF actually solving?** (`spf_costheta_*`, spf_mirror.py 101–114)
- `sin²θ` shape, `p(x)∝1−x²`: CDF inversion solves `x³ − 3x + (4u−2)=0`; in-range trig root
  `x = 2cos(⅓·arccos(1−2u) − 2π/3)`.
- `¼+¾cos²θ` shape, `p(x)=¼+¾x²`: solves `x³ + x + (2−4u)=0` (discriminant >0 always ⇒ single real
  Cardano root). Endpoints verified `u=0→x=−1`, `u=1→x=+1`; both `std::clamp`ed against fp drift.

**Q8. The a₂-linearity identity — what is it and why does it prove the mode weights are right?**
- The emission intensity `I(θ) = w_perp·sin²θ + w_par·(¼+¾cos²θ)` is exactly a **P₂ multipole**:
  `I ∝ 1 + a₂·P₂(cosθ_B)`, with `a₂ = (−⅔w_perp + ½w_par)/(⅔w_perp + ½w_par)`.
  → pure-A `a₂ = −1` (oblate ⟂ B̂), pure-B/C `a₂ = +1` (prolate ∥ B̂), unpolarized `a₂ = 0`.
- Because the wall load is **linear in the emission multipole**, the NWL map of any mixed `(a,b,c)` must
  equal the same linear combination of the pure-mode maps. `predictor_test.py:benefit()` sweeps
  `a₂∈[−1,1]` as `unpol + t·(par−unpol)` and finds `a₂_opt`; the linearity of that interpolation is the
  identity. If mode-selection probabilities were wrong, linearity would fail.

**Q9. Thread safety?**
- **No mutable members.** Every field (grids, CDFs, `w_perp/w_par/p_perp/polarized_/field_model_`) is
  set in the constructor and `const` at `sample()` time. Each thread owns its `seed`. Verified by
  construction and asserted in the header comment.

**Q10. RNG discipline & reproducibility?**
- `openmc::prn(seed)` only — never `rand()`/`<random>`. Fixed draw **order** (cell: ρ→ζ→θ; jitter:
  ρ→θ→ζ; direction: selector→shape-u→φ; then energy, then time).
- The Python mirror (`spf_mirror.Rng`) is a **bit-exact port of OpenMC's PCG-RXS-M-XS** `prn`
  (`src/random_lcg.cpp`), including `future_seed`/`init_seed`, so C++ vs Python is bit-for-bit on an
  identical seed sequence — the transcription proof.

**Q11. Gram-Schmidt degeneracy — how do you avoid the b̂∥ẑ blow-up?**
- Do **not** build the frame from `b̂×ẑ`. Pick the world axis **least aligned** with b̂ as the
  reference (`ref = argmin |b̂·ê|` over x̂,ŷ,ẑ), then `ex = normalize(ref − (ref·b̂)b̂)`, `ey = b̂×ex`.
  The reference is never near-parallel to b̂, so `ex` never degenerates. Tier-1 explicitly tests
  b̂=±ẑ.

**Q12. Birth weight — always 1?**
- Yes, spatially: the cascade is **rejection-free** (marginal→conditional CDF), so every position is an
  exact draw with weight 1. `site.wgt = E_wgt·t_wgt`, which is 1 for the default monoenergetic/delta
  distributions and carries importance weights only for biased energy/time.

**Q13. Jitter can push ρ a half-cell past the axis/edge — is that a bug?**
- No: ρ is **clamped** (non-periodic) in the interp bracket; the half-cell extrapolation is a bounded
  discretization bias that **→0 under grid refinement** (documented, matches the reference exactly).
  θ/ζ wrap periodically with a non-negative modulo (correct even for jitter just below `grid[0]`).

## B3. Physics

**Q14. Derive w_perp, w_par, η from Schwartz Eq. 2.**
- Eq. 2: `dσ/dΩ = (σ₀/2π)[ (¾)a·sin²θ + (⅔b+⅓c)(¼+¾cos²θ) ]`, θ measured from **local B̂**.
- Angular weights: `w_perp = ¾a` (the sin²θ piece), `w_par = ⅔b + ⅓c` (the ¼+¾cos²θ piece).
- Total-rate factor: `η = ∫dσ/dΩ dΩ / σ₀ = a + ⅔b + ⅓c`. Verified symbolically in
  `python/abc_modes.py`/`tests/test_abc_modes.py`.

**Q15. Why is C-mode NOT isotropic (a common mistake)?**
- Set `a=0` in Eq. 2: both the b and c terms multiply the **same** angular function `¼+¾cos²θ`. So B
  and C are **identical in shape** (peak ∥ B̂) and differ **only in total rate**: pure-B `η=⅔`, pure-C
  `η=⅓` (C is half the rate of B, same directionality). C is emphatically not isotropic.
  `abc_modes.angular_shape_b_minus_c()` proves `bracket_B − 2·bracket_C ≡ 0` symbolically.

**Q16. Where does `P_perp = a/η` come from?**
- Solid-angle integrals: `∫sin²θ dΩ = 8π/3`, `∫(¼+¾cos²θ) dΩ = 2π`. So the integrated weights are
  `W_perp = w_perp·8π/3 = 2π·a` and `W_par = w_par·2π`. Then
  `P_perp = W_perp/(W_perp+W_par) = 2πa/(2πa + 2π·w_par) = a/η`. Derived in-code, no magic constants.

**Q17. Corner totals (your oracle)?**
- Non-pol `(⅓,⅓,⅓)`: `η=⅔`. Pure-A `(1,0,0)`: `η=1` (**+50%** vs 2/3). Pure-B `(0,1,0)`: `η=⅔`.
  Pure-C `(0,0,1)`: `η=⅓` (**half** non-pol). All verified as corner cases in `abc_modes.CORNERS`.

**Q18. Angular vs rate separation — why is it the key design decision?**
- The **direction pdf** depends only on `(w_perp, w_par)` → `P_perp`; the **overall rate** `η` is a
  separate scalar that multiplies the spatial source strength. This lets A and C emit different
  *numbers* of neutrons but each emit correctly *shaped* directions — mirroring how Schwartz separates
  angular factors from rate factors. In code: the sampler only ever uses `P_perp`; `η` is exposed for
  rate/TBR accounting (Tier 7).

**Q19. (a,b,c) from spin fractions?**
- Schwartz Eq. 1: `a = d₊t₊ + d₋t₋`, `b = d₀` (`= d₀t₊ + d₀t₋`), `c = d₊t₋ + d₋t₊`. Implemented as
  `spin_fractions_to_abc` (Python) / documented in the header. Non-polarized fuel ⇒ `a=b=c=⅓`.
- **Footnote for Ethan:** Bae et al. use `(a,b,c)` for *different* symbols — flag this collision
  (already an appendix footnote in the paper plan).

**Q20. Energy handling — monoenergetic 14.06 MeV?**
- Default `Discrete([14.06e6],[1.0])` — energy is decoupled from the polarized angular distribution to
  leading order, and the steering/NWL verification is energy-independent. TokamakSource carries
  per-flux-surface Ballabio spectra; the stellarator prototype uses **one** distribution for all ρ
  (documented limitation). A per-r spectrum is a drop-in extension (the XML already loops `<energy>`).

## B4. Geometry & equilibrium

**Q21. Where does the flux map come from?**
- `python/quasr_fluxmap.py` in a DESC venv: solves a fixed-boundary DESC equilibrium from a QUASR
  boundary (or loads a converged DESC example), evaluates `sqrt(g),R,Z,phi,B,V` on a full-torus
  `(ρ,θ,ζ)` tensor grid (default 16×64×192), and writes both `.npz` and the portable `spf_fluxmap_v1`
  `.meta/.bin`. The Ginsburg batch used VMEC equilibria (job 8885457).

**Q22. What rigor gates guard the fluxmap? (Ethan will care about √g correctness.)**
- Asserted at write time (`quasr_fluxmap.py`, not warnings):
  **V1** own-quadrature volume `Σ|√g| dρdθdζ` vs DESC `V` (<3e-2, trapezoid w/ axis point);
  **V1b** independent LCFS volume via the **divergence theorem** `V=⅓∮X·(X_θ×X_ζ)` (<1e-2);
  **V2** √g single-signed on ρ>0 (nested surfaces);
  **V3** DESC node set is a clean `(ρ,θ,ζ)` tensor product after explicit lexsort;
  **V4** ρ=1 surface reconstructs the fitted QUASR boundary (<5e-2 of extent);
  **V0** force-balance convergence. A map that fails never reaches the sampler.

**Q23. Why √g volumetric weighting, and why does it matter?**
- Birth-position density is `p(ρ,θ,ζ) ∝ S(ρ)·|√g|` — `√g` is the metric Jacobian (real volume element
  of the flux coordinates). Without it you'd sample uniformly in coordinate space and **over-weight the
  compressed inboard flux surfaces**, biasing the source (and hence the wall load) toward the wrong
  poloidal/toroidal distribution. `S(ρ)` is the emissivity profile (default parabolic `1−ρ²`). Uniform
  `(ρ,θ,ζ)` spacing means the constant cell volume cancels in the normalized CDFs, so only `S·√g`
  enters.

**Q24. The field model — fluxmap b̂ vs toroidal override?**
- `field_model='fluxmap'` (default): trilinear-interp `(BR,Bφ,BZ)` from the **same grid** as position,
  rotate rpz→Cartesian using the interpolated φ, normalize → the real 3-D field direction.
- `field_model='toroidal'`: override with pure φ̂ = (−sinφ, cosφ, 0) — the **axisymmetric limit** used
  to cross-check against TokamakSource's B̂=φ̂ case and the Schwartz/anarrima analytic oracle.
- (The SPF-patched TokamakSource instead offers `toroidal` vs `pitched` = φ̂ + s·λ·p̂ from a q-profile.)

**Q25. How is the conformal wall / blanket built?**
- `python/stellarator_geometry.py`: take the LCFS grid `R(θ,φ),Z(θ,φ)`, offset **outward along the
  poloidal cross-section normal** by cumulative layer thickness, emit each layer as a **watertight
  closed-shell STL**; `build_dagmc.py` (Ginsburg/conda) → DAGMC `.h5m`; `run_conformal.py` transports.
- Guards: edge-manifold watertightness + a `_poly_is_simple` **self-intersection check** per φ
  cross-section (naive normal-offset folds on the concave inboard side beyond the curvature radius),
  and a signed-volume orientation guard (G6) so outer normals point outward for either θ-handedness.

**Q26. The DEFAULT_LAYERS thicknesses — where from?**
- `DEFAULT_LAYERS` (`stellarator_geometry.py`, cm, plasma→out): **sol 5.0** (scrape-off vacuum gap),
  **W 0.2** (tungsten first-wall armor), **steel 3.8** (RAFM structural FW), **Be 2.0** (multiplier),
  **FLiBe 50.0** (breeder), **shield 40.0** (WC/steel shield), **coil 15.0** (coil block, damage
  tally). ~111 cm total build.
- Provenance: **ARC-class / ARIES-CS reactor radial-build proportions** (reused material identities
  from `reactor_model`); flagged `INJECT(helios)` — placeholder public values, swap for real Helios
  numbers. The scatter run scales the device ×10 so the ~110 cm blanket fits a reactor-size plasma
  (b̂ is scale-invariant).

**Q27. Where is the conformal-wall offset gap set?**
- `0.30·a` outward from the LCFS in `predictor_test.geometry()` / `conformal_wallmap` (`a` = half the
  R-extent). The 59509 free-stream validation used `0.40·a`.

## B5. Results & interpretation

**Q28. Schwartz recovery is +39/−21% (MC) vs his +43/−22% analytic — why not exact?**
- The MC numbers are the free-streaming StellaratorSource on the square torus; the ~4/1-point gap is
  the **ε_eff≈0.43 shaping residual**, not a sampler error. On the toy-QA-low patch (ε_eff≈0.08) MC↔
  analytic agree ≤1.8% on **all** walls; at ε_eff≈0.43 the gap appears **only on the inboard wall**
  (closest, most curved), 4.1%/5.7% (`RESULTS_quasr.md`). Both methods were eliminated as the cause:
  sampler exact (A/iso=(3/2)sin²θ to 0.3%/bin), analytic ≡ anarrima to 1.2e-12, insensitive to
  visibility, field, and near-field gap.

**Q29. What exactly causes the inboard residual?**
- The **point-patch approximation** of the filamentary free-streaming *analytic* on a strongly-curved
  wall — a limitation of the analytic, with the sampler-exact/geometry-exact MC as the reference.
  `bare_cylinder_demo.py` isolates it: single shaped filament, analytic point-patch vs exact
  ray/surface intersection; error grows **≈ s²** with shaping amplitude and is **~3× larger on the
  convex inboard cylinder** (|Δ| 0.002→0.098 as s:0→1). This is the paper's "where MC is necessary."

**Q30. Scattering dilutes steering ~3× — physical reason?**
- Free-streaming, A-mode steers inboard **+40% / −20%** outboard; with the full W/steel/Be/FLiBe/shield
  blanket ON it drops to **+13.2% / −9.2%** (`steering_dilution.csv`, Tier 5). Physically: the
  polarized anisotropy is a **birth-direction** effect; every scatter **randomizes direction**, so the
  wall load is a mix of uncollided (steered) + scattered (isotropized) neutrons. Roughly a factor ~3
  of the first-flight lever survives to the wall. **This is the headline MC-necessity result** — a
  free-streaming/analytic method overstates the benefit ~3×.

**Q31. FLiBe = 84% of heating — is that expected?**
- Yes: the 50 cm FLiBe breeder is the thickest layer and the primary neutron-energy sink. The
  volumetric-flux run (device 59509 ×10, scattering ON, mesh (nr,nφ,nz)=(39,32,39)) gives FLiBe = 84%
  of the 1.1e7 eV/src heating (`PROVENANCE.md`). Near-wall fast flux modulates **perp −10% / parallel
  +10%** vs unpolarized — the volumetric echo of the surface steering.
- Caveat to volunteer: the volumetric mesh shows heating, **not** inboard/outboard steering — the
  near-field steering asymmetry is a **surface** effect a volumetric mesh can't resolve (would need a
  first-wall surface-current tally).

**Q32. Free-streaming cross-check ON the conformal wall (not just the box torus)?**
- `conformal_freestream_xcheck.py`, device 59509: analytic point-patch vs ray-traced OpenMC
  StellaratorSource on the **real 3-D conformal wall**. Agree **0.5–2.3% by region** (inboard A/iso
  1.6%, B/iso 2.3%; outboard 0.5%/0.6%; global 0.9%/1.2%), per-bin r≈0.91 (`conformal_xcheck.csv`).
  Validates the source free-streaming on real geometry, not only the square torus.

**Q33. Directional efficiency η = 0.89 — definition?**
- Fraction of the axisymmetric steering **retained under the real pitched 3-D field**. Tier 8b: parallel
  mode reduces the inboard load **2.2%** with the 3-D QA field vs **2.4%** with pure toroidal ⇒
  η = 2.2/2.4 ≈ **0.89** (`RESULTS_tier8_3d.md`). The pitched field tilts B̂ off φ̂, costing ~11% of the
  ideal-toroidal lever. Emission moments confirm the sampler: `<(u·B̂)²>` = 0.201/0.333/0.467 (perp/
  unpol/par) vs analytic 0.200/0.333/0.467 (`emission_moments.csv`).

**Q34. The n=15 predictor NULL — what did it show?**
- SPF global benefit `Y=(PF_unpol−min_a2 PF)/PF_unpol` is **weak: 0–16%, median ~4%** across the 15
  devices (`predictor_n15.csv`).
- Best predictor is the **transport-free headroom proxy PF_geo (ρ=+0.525, borderline)**; **plasma-shape
  metrics do not forecast it** (nematic_Q ρ=+0.396 n.s., shaping_amp ρ=−0.186). The strong n=5
  correlation (PF_geo ρ=+0.90) **regressed to noise** at n=15 — a properly-powered null.

**Q35. Why is the SPF benefit weak for stellarators / why doesn't plasma shape forecast it?**
- Peaking is dominated by a **narrow near-field grazing spike** (0.3–2.8% of wall area) that the P₂
  angular lever **cannot move**; SPF only steers the **broad, low-harmonic** in/out asymmetry (the
  ±22–43% lever). So benefit tracks **available headroom**, not shaping.
- Mechanistically, free-streaming line-of-sight integration is a **spatial low-pass filter**
  `T(N)~exp(−N·a/R₀)` (`RIPPLE_AND_PREDICTOR.md`): the wall **structurally cannot resolve** high-
  harmonic plasma shaping, so a shape-harmonic predictor is decoupled by construction. (Honest caveat:
  the low-pass is directly validated for the *ripple*; extending it to "shape can't predict steering"
  is supported by the n=15 null but not independently proven — `shaping_amp` did not cleanly decouple.)

**Q36. Device 803097 is the interesting outlier — why?**
- Peaked (PF_unpol=1.40) but **unsteerable** (Y=0.014, ~10× below the rest) and the **only perp-optimal
  device** (a₂_opt=−0.83). PF_geo flags it as *low* peaking (1.19), so its transport peaking comes from
  an anisotropic near-field structure the isotropic proxy can't see — both perp-favorable and
  steering-resistant. This is "when can SPF **not** help."

**Q37. Rate-vs-steering trade — why is B the "sweet spot"?**
- Tier 7 (`rate_steering.csv`): A gives the biggest steering (+70% center-stack) but **η=1.0** requires
  full alignment; C halves the rate (η=⅓). **B keeps the non-polarized rate (η=⅔, fusion_rate=1.00)
  while still steering** (center-stack 0.85, outboard 1.09) — best rate-preserving directional lever.
  TBR/neutron is nearly mode-flat (~1.15).

## B6. Comparison to prior art

**Q38. How does this differ from Bae et al. 2025?**
- Bae = **precomputed static point-cloud** polarized source, one fixed tokamak, tabulated and replayed.
  Ours is **native / on-the-fly / general**: samples position+polarized direction inside the transport
  loop from a real equilibrium, so it carries **local-field polarized emission** (b̂ varies in space) —
  a tabulated cloud cannot. One-liner: *"Bae replays a frozen tokamak source; we generate a
  local-field-aware polarized source on the fly, for any equilibrium."*
- Do **not** claim "first polarized source in OpenMC" (Bae has that); **do** claim "first *native/
  general/on-the-fly* polarized source in a production MC code."

**Q39. vs Lyytinen 2024?**
- Lyytinen = **Serpent2**, HELIAS, a **ported MCNP isotropic** source with parametric *geometry* — no
  polarization, different code. One-liner: *"Lyytinen ports an isotropic stellarator source into
  Serpent2; we add spin-polarized, local-field emission natively in OpenMC."* Don't claim "first MC
  stellarator source."

**Q40. vs Peterson's own TokamakSource?**
- TokamakSource = analytic axisymmetric tokamak geometry, isotropic (unpolarized) emission, per-surface
  spectra. We (a) add the **SPF polarized direction** hook (upstreamable into TokamakSource itself) and
  (b) add a sibling **StellaratorSource** for real 3-D equilibria with grid-field b̂. One-liner: *"same
  Source pattern; we add polarization + a 3-D-equilibrium position/field model."*

**Q41. What's genuinely first-of-kind here?**
- First **polarized stellarator** source (b̂ from a real 3-D equilibrium) — nobody has this. First
  **equilibrium-rigorous StellaratorSource in OpenMC** (gated √g volumetric). Plus the honest-findings
  results (3× scattering dilution, MC-necessity boundary, B sweet spot).

## B7. Limitations & likely pushback

**Q42. Monoenergetic — isn't that unphysical for heating?**
- Yes for absolute heating; but the **steering/NWL physics is energy-independent** to leading order,
  and the verification isolates source direction. Per-r Ballabio spectra are a drop-in (the XML already
  loops `<energy>`). Stated as a limitation.

**Q43. No depolarization physics in transport?**
- Correct — this is a **birth-distribution** model; once born, neutrons transport unpolarized (standard
  cross sections). Depolarization/spin-dependent scattering is out of scope (as in Bae). The 3× scatter
  dilution already captures the dominant "steering washes out" effect via ordinary isotropizing
  scatter.

**Q44. Convex-ish wall assumption?**
- The analytic companion assumes a **convex cross-section / single visible arc** (inherited from
  Schwartz); non-convex features (divertor legs) can split the visible region — flagged as future work.
  The **MC** side has no such restriction (it ray-traces the true DAGMC wall); the assumption bounds
  only the analytic oracle. The conformal builder also refuses self-intersecting (over-offset concave)
  walls.

**Q45. Steering is "just" a source/free-stream effect that scattering dilutes — so what?**
- That **is** the paper's point, stated honestly: SPF is a **local actuator, not a global flattener**,
  and quantifying the ~3× dilution is *why you need transport* rather than the analytic. Two independent
  methods (MC + anarrima) agree on the weak-global-lever finding.

**Q46. Single device for the scatter run (59509) — is that enough?**
- The **3× dilution** headline currently rests on the box-torus Tier-5 blanket + the 59509 conformal
  volumetric run. Broadening to more conformal *scattering* devices is the acknowledged to-do (the
  script `run_conformal.py` exists, wired for precise_QA, not yet looped over QUASR devices). The
  15-device set is free-streaming NWL maps, not scattering. State the fidelity ladder explicitly.

**Q47. The box/square torus — is that your physics geometry?**
- **No — validation only.** Per Ethan's own earlier guidance: the square/box torus is for sampler +
  NWL verification and Schwartz axisymmetric recovery **only**; the physics belongs on the conformal
  QUASR wall. The fidelity ladder is explicit in `PAPER_DIRECTION.md` (5 rungs, each isolating one
  effect).

**Q48. Are the blanket materials/thicknesses real Helios values?**
- No — **public ARC/ARIES-CS placeholders**, every device-specific value marked `INJECT(helios)`. The
  machinery is validated; the specific radial build is a stand-in pending real Helios numbers. Not a
  Helios physics claim.

**Q49. anarrima cross-checks you did NOT do (be upfront)?**
- **[from PAPER_DIRECTION honest-gaps]** polarized-at-grazing (singularity subtraction) is
  by-construction, **untested**; `jit(grad)` through the anarrima engine hangs (eager grad works ~9s);
  the anarrima↔OpenMC uncollided cross-check was never run. These are Limitations, not claims.

**Q50. What would you want from Ethan specifically?**
- Co-author sign-off on §7–§10; a sanity check on the TokamakSource upstreaming path for the
  `<polarization>` hook; and whether the `spf_fluxmap_v1`-vs-HDF5 choice is acceptable for an eventual
  upstream (or whether it should read the existing OpenMC HDF5 mesh/field infrastructure).

---

### Quick numbers cheat-sheet
- Corner η: non-pol ⅔, A 1, B ⅔, C ⅓. `P_perp=a/η`. `∫sin²=8π/3`, `∫(¼+¾cos²)=2π`.
- a₂: A=−1, unpol=0, B/C=+1.
- 59509: nfp3, aspect 6.7, e_geom 0.29, ε_eff 0.43. precise_QA ε_eff 2.3 (MC-only).
- Scatter dilution: A inboard +40%→+13.2%, outboard −20%→−9.2% (~3×).
- FLiBe 84% heating; near-wall flux perp −10%/par +10%. η_dir=0.89. TBR/n ~1.15.
- n=15 benefit: 0–16%, median ~4%; PF_geo ρ=+0.525 (best, borderline); nematic ρ=+0.396 (n.s.).
- Conformal free-stream xcheck: 0.5–2.3% by region, per-bin r≈0.91.
