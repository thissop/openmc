# RESULTS — exact masses + clearance field + design domain (2026-08-25)

Two independent methods, agreeing → the mass table is trustworthy. No "equal-mass" claim is made where the
computed costs disagree beyond tolerance.

## 1. AUTHORITATIVE per-patch added-shield mass — clean 40M-sample VolumeCalculation (direct DAGMC)
The definitive masses are the **direct `openmc.VolumeCalculation` on each real DAGMC** (40M samples, cells by exact
material match). This measures the actual transport geometry, so it is primary over any surface approximation.
Added mass = M(config) − M(zero-δ); the fixed-envelope swap is confirmed (dV_shield ≈ −dV_breeder, swapResid ≤ 6%).

| patch | ΔV_shield (cm³) | **ΔM_net vs zero-δ (kg)** | ΔM_net vs uniform (kg) |
|-------|-----------------|---------------------------|------------------------|
| t0p2 | 2.001e7 | 274,406 ± 15,842 | −255,027 |
| t0p3 | 2.800e7 | 384,029 ± 15,887 | −145,404 |
| t1p0 | 1.504e7 | 207,138 ± 15,813 | −322,295 |
| t1p1 | 1.256e7 | 172,216 ± 15,799 | −357,216 |
| t1p2 | 1.512e7 | 207,193 ± 15,814 | −322,240 |
| t1p3 | 2.069e7 | 284,595 ± 15,846 | −244,838 |
| t2p0 | 1.832e7 | 251,578 ± 15,832 | −277,855 |
| t2p1 | 1.184e7 | **163,086 ± 15,795** (min) | −366,347 |
| t2p2 | 1.418e7 | 193,956 ± 15,808 | −335,476 |
| t2p3 | 1.385e7 | 190,154 ± 15,806 | −339,279 |
| t3p2 | 2.503e7 | 343,095 ± 15,870 | −186,338 |
| t3p3 | 2.906e7 | **399,947 ± 15,893** (max) | −129,486 |
| **uniform (+2.55cm)** | 3.861e7 | **529,433 ± 15,948** | 0 (ref) |

**Per-patch added mass varies 2.45× (163 t → 400 t) for the same +20 cm δ** — "equal thickness ≠ equal mass",
and *every* patch is lighter than the uniform +2.55 cm control (−129 to −366 t vs uniform). The prior "t3p3 ≈
uniform at equal mass" claim is **quantitatively refuted** (t3p3 is 129 t lighter than uniform, 8σ).

## 1b. ParaStell-surface cross-check (shape-accurate; corrects an earlier over-claim)
The surface `get_loci()` area×thickness estimate (see §earlier drafts) reproduces the per-patch SHAPE but is
**biased +6–19% high in absolute volume** vs the MC ground truth (interior patches), and a normal-offset hex-volume
"correction" was WORSE (+20%) — so it was reverted. Boundary patches (I_TOR∈{0,3}) additionally carry a **seam
artifact**: `build_patch` extends the delta across the field-period seam, which a plain 4×4 sector re-binning of the
surface area misses (e.g. t3p2 surface 1.71e7 vs MC 2.50e7, −32%). Lesson (integrity rule): the surface area×t is a
geometry-consistent SHAPE field, not an exact mass; **use the MC for absolute mass**. The total-shell hex volume does
match MC (true-hex 6.45e8 = −2.6% vs MC 6.62e8), confirming the loci geometry itself is correct — the bias is purely
the area×thickness vs true-swept-volume approximation for the thin δ increment plus the seam binning.

## 2. Method provenance
The §1 table is the 40M-sample MC (errors ±15.8 t net, swapResid ≤ 6%). A 5M pilot first flagged the
`shield`↔`thermal_shield` substring bug (fixed via exact material match) and established the swap + not-equal-mass
findings at coarser precision; the 40M rerun tightened them 3×. The §1b surface method is retained only as a
shape-consistent field for the design families (§4), calibrated to the MC absolute scale.

## 3. Clearance field C(θ,φ) and design domain
C = min 3D distance from thermal-shield OUTER surface(θ,φ) to nearest coil filament (coils_qh, m→cm, nfp-replicated;
scale validated: coil R0≈1366 cm matches blanket). Field: **min 29.7 cm, median 314.5, max 783.8 cm; tightest at
tor=3°, pol=173° (INBOARD).** Available growth (C − 25 cm coil half-width): **min 4.7 cm inboard**, median 289.5.

**Key coupling:** the tightest clearance is inboard (pol≈180°) — exactly where coil-20 loading peaks. Response-need
and feasibility collide inboard: the growth family (I2) has almost no room where shielding is most wanted. This is
the concrete "attribution ≠ intervention worth ≠ feasible controllability" statement, now with numbers.

Design domain:
- **I1 fixed-envelope swap** (design var δ_shield): δ ∈ [−30, +35] cm (breeder floor 15, shield floor 10);
  clearance-independent (envelope fixed). Uniform bound across (θ,φ).
- **I2 growth into clearance** (extra build g toward coils): g(θ,φ) ∈ [0, avail(θ,φ)], avail 4.7–743 cm — the
  spatially-varying, coil-exclusion-limited bound. This is where the clearance field bites.
- Manufacturability: |∇t| ≤ 15 cm/step; stellarator symmetry: 1 period × nfp=4.

## Files
`surface_extract.npz` (loci, areas, exact masses), `clearance_domain.npz` (C field, avail, t_min/t_max),
`masses_local.json` (MC). Pending: full DAGMC `find_cell` overlay of shield-surface points (geometry gate); 40M MC
rerun cross-check; coil-placement validation vs DAGMC coil cells.

## 5. Corrections & claim status (2026-08-25, per review)
**Repeated-geometry bug ROOT-CAUSED (not relaunched blindly):** `coil_run_v3_seed.py:15` hardcoded
`os.chdir(~/pstl_test/corrected)` (shared). Concurrent array tasks exported model.xml/materials.xml into that one
dir; C++ init read whichever won the race -> ALL tasks transported one geometry (tags stayed correct as filename
labels). Fixed in `coil_run_v3_iso.py`: physically-isolated per-task RUNDIR, absolute immutable DAGMC paths,
pre-transport manifest (config, task id, seed, rundir, DAGMC+model.xml sha16, bbox). Local smoke gates PASS: 4
distinct DAGMC hashes, distinct V_shield, find_cell differs at the patched sector. Cluster low-history smoke
(job 9494080) verifies statepoint isolation + non-byte-identical response vectors before the real rerun.

**Filament-clearance (renamed):** the §3 field is min Euclidean distance thermal-shield-outer -> coil FILAMENT
centerline (min 29.7 cm, avail-to-grow 4.7 cm). This is NOT engineering clearance: real t_max needs the full
winding-pack swept/exclusion solid (coil cross-section + structural margin) and a collision/geometry check on final
candidates. The 29.7 cm / 4.7 cm are filament distances, not exact usable growth. [engineering t_max: TODO]

**Matched-mass caveat:** only the 12 built patches have DIRECT MC volumes; sectors t0p0,t0p1,t3p0,t3p1 are
surface-derived, globally calibrated (ratio r~0.85). "0.00% matched" = matched under the calibrated sector-cost
MODEL, not exact equality in final DAGMC geometries. Each generated candidate must get its own VolumeCalculation +
thickness trim to hit the mass target. The seam artifact (I_TOR 0,3 cross the field-period seam) must be represented
identically in the intervention basis, B_FD row mapping, and cost map -- do not compare a seam-crossing transport
patch with a non-seam analytic sector.

**coil-20-guided vs response-aware (renamed):** the B_FD field that follows the coil-20 column is now
`coil-20-guided`. The true response-aware candidate is the constrained MINIMAX over the full 12x20 operator
(`minimax_candidate.py`): min peak-coil s.t. matched mass + bounds. Result (SCREENING only -- linear surrogate,
known non-generalizing, baseline confounded): coil-20-guided BACKFIRES on the full vector (predicted new peak c30
+112% vs base) while minimax migrates the limiting coil c20->c30. Magnitudes are linear-extrapolation-unreliable;
the candidate MUST be transport-verified.
