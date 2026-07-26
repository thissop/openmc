# Auto coil design, QUASR engineering-relevance filtering, and the controlled-standoff experiment

*Feasibility + design memo for the magnet-adjoint program. Written 2026-07-25. Grounds:
`SCIENCE_STRATEGY.md` (Step-2 free-streaming coil law; Step-4 QI spike; standoff confound),
`DATABASES_AND_RELATIONAL_DESCRIPTORS.md` (QUASR = spine, standoff = dominant relational descriptor),
`EXECUTION_PLAN_magnet_adjoint.md` (ParaStell→DAGMC→adjoint per-device cost). Every tool-availability
claim below was **probed on this machine** (base conda, simsopt 1.10.6) and is flagged as verified or
unverified. Honesty rule: where a capability is paper-only, absent from the installed build, or
untested at scale, it is marked as such — a described capability is not an installed one.*

---

## 0. Report-back (read this first)

**Recommended coil-design tool + effort.** **simsopt stage-2 filamentary optimization** is the right
near-term engine, and it is **already installed and importable on the group's machine** (base conda,
`simsopt 1.10.6` — the memory note that the group parses QUASR "without simsopt" is about their *loader*,
not about availability). Every class the standoff-controlled recipe needs is present and imported clean:
`BiotSavart`, `Coil`, `Current`, `coils_via_symmetries`, `CurveXYZFourier`,
`create_equally_spaced_curves`, `SurfaceRZFourier`, `CurveLength`, `CurveCurveDistance`,
**`CurveSurfaceDistance`** (the standoff knob), `MeanSquaredCurvature`, `SquaredFlux`. Standoff is
controlled directly by the `CurveSurfaceDistance(target=d)` penalty. Effort to a working
boundary→coils→ParaStell-ready filament file: **~1 day** to stand up and validate on one QUASR boundary
(against its shipped coils as ground truth), **~2–4 hours per additional (plasma, standoff) point**
thereafter, mostly weight-tuning babysitting. **REGCOIL is NOT the near-term path** — the current-potential
classes (`CurrentPotential`, `CurrentPotentialSolve`) are **absent from the installed simsopt build**
(verified: `'CurrentPotential' in dir(simsopt.field)` → `False`; the build ships without the compiled
regcoil backend). REGCOIL lives cleanly in **DESC** (pure-Python `run_regcoil` + `CoilSet`) if you want
a winding-surface offset knob, but DESC is **not installed here** and adds a JAX/DESC dependency; treat it
as the fallback, not the default.

**Surviving fraction + reachable n (the two numbers asked for).** Under the sharp **radial-build fit**
cut (min coil-plasma standoff ≥ the ~1 m reactor radial build), **~95% of QUASR survives** — the cut
barely bites because QUASR's baked-in 0.1 m min-distance ≈ 1 m at reactor scale, so the devices are
*already* build-clearing by construction (computed on the group's real 366-device standoff cache;
requiring shield *headroom* beyond the minimum, gap ≥ 150 cm, cuts to ~70%). The **binding** engineering
cut is **aspect ratio**: radial-build-fit **AND** compact aspect (4 ≤ A ≤ 10) leaves **~15–30%** of
QUASR. On the compute side, the pipeline is embarrassingly parallel and the **DAGMC build (~30–60 min) is
the bottleneck, not coil design (minutes)** — so a modest SLURM array (~25 concurrent tasks, license-free
`cad_to_dagmc`) clears **n ≈ 20–50 overnight**; compute is *not* the limiter for the n a defensible law
needs (coil-design convergence and watertight-`.h5m` success rate are). Details in §2.3 and §4.

**QUASR filtering vs auto-design-at-controlled-standoff — which first?** **Both, sequenced, and they are
complementary, not competing.** QUASR bulk-filtering is the cheap, real-*n*, no-new-tooling path and should
go first (it is SCIENCE_STRATEGY Step 2, uses only code the group already has). But filtering **cannot break
the standoff confound** — QUASR ships one coil set per boundary, so you get *one* standoff per shape and
cannot separate "shape effect" from "standoff effect." Auto-design-at-controlled-standoff is the **only**
way to get multiple standoffs at *fixed* shape, and it is exactly what the confound-breaking experiment
needs. So: **filter QUASR now to get the real-*n* free-streaming coil law (Step 2); stand up simsopt stage-2
next to run the small controlled-standoff experiment (the confound-breaker); use QI auto-design as the
Step-4 spike.** Filtering feeds the observational law; auto-design feeds the controlled experiment.

**Single biggest feasibility risk.** **stage-2 non-convergence / unrealistic coils at tight standoff.**
Filamentary stage-2 is non-convex; at a small `CurveSurfaceDistance` target the optimizer trades Bnormal
against the distance penalty and can converge to wavy, high-curvature, or coil-coil-touching geometries —
or fail to hit the target standoff at all. That directly threatens the "tight" arm of the controlled
experiment (the most scientifically valuable one). Mitigation: warm-start each standoff from the previous
(looser) solution, cap curvature/length with penalties, and **hard-gate** each output on
achieved-standoff, max-curvature, and min-coil-coil before it enters ParaStell. A secondary risk is that
the standoff knob (`CurveSurfaceDistance`) controls the *minimum* distance, not a *uniform* offset, so a
"tight" coil may be tight only locally — which is physically fine (min standoff is the reactor-binding
quantity per `DATABASES_AND_RELATIONAL_DESCRIPTORS.md` §2) but must be **measured and reported**, not
assumed uniform.

---

## 1. Coil-design tool assessment (QUESTION 1)

### 1.1 Assessment table

| Tool | What it does | Automation | Cost / device | Standoff control | Output format | Availability **here** | Main risk |
|---|---|---|---|---|---|---|---|
| **simsopt stage-2** (filament) | Fix plasma; optimize N `CurveXYZFourier` coils to min `SquaredFlux` (Bnormal) + length + curvature + coil-coil + coil-plasma penalties | **High** — scriptable; QUASR ships the exact recipe (`stage_two_optimization*.py`) | Minutes (L-BFGS on ~10² DOFs); +human weight-tuning | **Direct**: `CurveSurfaceDistance(target=d)` sets the min coil-plasma distance = the standoff knob | Filament polylines → **MAKEGRID/`.coils` or point-locus** = exactly what ParaStell + the group's adjoint filament source consume | **INSTALLED & VERIFIED** (simsopt 1.10.6, all classes import) | Non-convex: non-convergence / wavy coils at tight standoff; weight-tuning babysitting |
| **REGCOIL** (current potential on winding surface) | Solve a *linear* regularized least-squares for a sheet current on a winding surface offset from the plasma; cut contours → filaments | **High & robust** (convex linear solve) | **Seconds** (linear solve) + contour-cutting | **Direct & clean**: standoff = **winding-surface offset distance** (the natural REGCOIL knob) | Current-potential contours → filament coils (extra cut step) | **simsopt: ABSENT** (`CurrentPotential*` not in build — verified). **DESC: present in code, NOT installed here** | Coils can be unrealistic at small offset (high current density / tight bends); needs a winding-surface constructor; contour-cutting adds a step |
| **DESC coil tools** (`run_regcoil` + stage-2 `CoilSet`) | Pure-Python REGCOIL **and** filamentary stage-2; reads DESC/VMEC equilibria | High | Seconds–minutes (JAX-accelerated) | REGCOIL offset **or** stage-2 min-distance penalty | DESC `CoilSet` → export to MAKEGRID `.coils` | **NOT installed here** (adds DESC+JAX dep) | New dependency stack; less battle-tested filament export into ParaStell than simsopt |
| **NESCOIL / BNORM** (STELLOPT) | Original (un-regularized) current-potential method | Medium | Seconds | Winding-surface offset | Fortran; `.coils` via cut | Not installed; heavyweight Fortran build | Un-regularized → noisy coils (REGCOIL exists precisely to fix this); legacy tooling |
| **pyoculus** | Field-line / island analysis | — | — | — | — | Not a coil-design tool | Wrong tool (analysis, not synthesis) — do not use for this |

### 1.2 REGCOIL, specifically (as asked)

- **How automatic/robust:** REGCOIL is the *most* robust of the family because the core solve is a
  **regularized linear least-squares** (current potential on a fixed winding surface), not a non-convex
  shape optimization. Given a boundary + a winding surface, one regularization sweep (`lambda`) is a
  handful of seconds and essentially always converges. That is its selling point over stage-2.
- **Standoff knob:** **Yes, cleanly** — the standoff *is* the winding-surface offset distance, and
  building the winding surface as a uniform normal-offset of the plasma boundary makes standoff a single
  explicit scalar. This is arguably a *better* standoff knob than stage-2's min-distance penalty because
  the offset is (approximately) uniform by construction, whereas `CurveSurfaceDistance` pins only the
  minimum.
- **Output:** current-potential contours, which must be **cut into discrete filament coils** (an extra,
  well-established step). simsopt/DESC both provide the cutter; the cut filaments then export to
  MAKEGRID `.coils` for ParaStell.
- **Compute cost:** seconds per device per `lambda`; cheaper than stage-2.
- **Failure modes:** (i) at **small offset** the sheet current develops sharp features → cut coils are
  high-current-density, tightly-bent, possibly non-buildable (this is the tight-standoff failure mode,
  same physics as stage-2's); (ii) you must *supply* a winding surface — a bad/self-intersecting offset
  surface at large offset or high curvature breaks it; (iii) contour cutting can produce awkward coil
  counts/topologies that need manual arbitration.
- **Availability verdict:** the **winding-surface offset knob is attractive for the controlled-standoff
  experiment**, but it is **not installed in this simsopt build** and would come via **DESC** (new dep).
  Recommendation: default to simsopt stage-2 (installed, direct ParaStell path); adopt DESC-REGCOIL only
  if stage-2's min-distance-only standoff proves too non-uniform for a clean 2-way separation.

### 1.3 simsopt stage-2, specifically (as asked)

- **Cost/robustness/automation:** minutes per device (L-BFGS on ~10²–10³ DOFs); fully scriptable. QUASR
  itself was **built** with this machinery, and the group already has a copy of the recipe on disk:
  `spf_prototype/shield_opt/data/raw/QH_unzipped/zenodo/optimization_scripts/stage_two_optimization_planar_coils_iterate.py`.
  So the recipe is not only available — it is the *provenance* of the coils the group already parses.
- **Standoff control:** the objective is
  `J = SquaredFlux + Σ wᵢ·penaltyᵢ`, with
  `CurveSurfaceDistance(curves, plasma_surface, target=d_standoff)` the coil-plasma min-distance penalty
  (the **standoff knob**), `CurveCurveDistance(target=d_cc)` the coil-coil gap, `CurveLength`/`LpCurve*`
  the length, and `MeanSquaredCurvature` the buildability cap. Raising `d_standoff` and re-optimizing =
  moving the coils out = the controlled-standoff sweep.
- **Output format:** filament curves → **MAKEGRID `.coils` / point-locus**, which is *exactly* what
  ParaStell ("VMEC + coil filament point-locus data") and the group's adjoint filament source
  (`reciprocity_check._parse_coils` / `coil_adapter.py`) consume. No format bridge needed.
- **Verified importable here** (probe): all of
  `BiotSavart, Coil, Current, coils_via_symmetries, CurveXYZFourier, create_equally_spaced_curves,
  SurfaceRZFourier, CurveLength, CurveCurveDistance, CurveSurfaceDistance, MeanSquaredCurvature,
  SquaredFlux` → OK.

### 1.4 QUASR download is also already wired in simsopt

`simsopt.configs` exposes **`download_ID_from_QUASR_database(ID, return_style=...)`** and
`get_giuliani_data` — verified present. It pulls
`https://quasr.flatironinstitute.org/simsopt_serials/{ID[:4]}/serial{ID:07d}.json` (the *same*
`simsopt_serials` JSON the group's `quasr_loader` already parses), returns base curves + currents + nfp +
full coil set, and caches locally. So a boundary+coils pull for any device is a one-liner — no scraping.

### 1.5 Verdict for Q1 (10–20 standoff-controlled coil sets from existing boundaries)

- **Best tool:** **simsopt stage-2** (installed, direct ParaStell output, standoff via
  `CurveSurfaceDistance`). Keep **DESC-REGCOIL** as the fallback if a *uniform* offset knob is needed.
- **Realistic effort:** **~1 day** to build+validate the driver on one QUASR boundary (validate by
  re-deriving coils close to the shipped QUASR coils — a free ground-truth check the group uniquely has);
  then **~2–4 h per (plasma, standoff) point**. 10–20 sets ≈ **3–6 focused days** including hard-gating
  and ParaStell hand-off, dominated by weight-tuning and the tight-standoff arm.
- **Main risk:** tight-standoff non-convergence / non-buildable coils (see §0). Warm-start from looser
  standoff; hard-gate on achieved-standoff + curvature + coil-coil before ParaStell.
- **Design constraint (ties to §2.2):** set the coil-design standoff target ≥ the **radial-build offset
  (~1 m at reactor scale)** so every generated device is fit-feasible by construction — with
  `CurveSurfaceDistance(target ≥ radial_build)` in stage-2, or a winding surface offset ≥ radial_build in
  REGCOIL. The controlled-standoff sweep (§3) then varies this target *above* the ~1 m floor, so every arm
  houses the blanket and the knob is pure shield-headroom.

---

## 2. QUASR bulk access + engineering-relevance filter (QUESTION 2)

### 2.1 Bulk access — four channels, ranked for this group

1. **Metadata dataframe (best for bulk *selection*):** `https://quasr.flatironinstitute.org/QUASR.pkl` —
   a pandas pickle with **all valid IDs + per-device metadata** (aspect, nfp, QS class, rotational
   transform, and coil-length *threshold* columns; note the Zenodo erratum that `total_coil_length` /
   `coil_length_per_hp` should read `..._threshold`). Pull this **once**, filter in-memory, then download
   only the surviving serials. This is the "query metadata for many configs without pulling all 370k"
   answer.
2. **Per-device serial JSON (what the group already parses):**
   `.../simsopt_serials/{ID[:4]}/serial{ID:07d}.json` — coils + boundary together, one HTTP GET each,
   via `simsopt.configs.download_ID_from_QUASR_database` or the group's `quasr_loader`. Chunk by
   downloading only filtered IDs.
3. **Zenodo bulk deposit:** `10.5281/zenodo.10050655` (concept) / `10.5281/zenodo.13717741`
   (~12.7 GB tar) — the whole zoo offline when you want everything.
4. **Web navigator** (`quasr.flatironinstitute.org`) — interactive, for eyeballing, not bulk.

**Recommended flow:** `QUASR.pkl` → filter on cheap metadata (aspect band, nfp, QS class) → download the
surviving `simsopt_serials` → compute the coil-referenced filters (§2.2) from the filaments the group
already parses. No REST/GraphQL is needed (none is documented); the pickle + serial-URL pattern *is* the
API.

### 2.2 Engineering-relevance filters + concrete thresholds

**Primary filter — "can the reactor radial build physically fit?"** The sharp, physically-grounded
definition of engineering-relevant is: **the minimum coil-plasma standoff must be ≥ the full reactor
radial-build thickness** that has to be stacked between plasma boundary and coils —
`first_wall + multiplier + breeder + back_wall + shield + gap + vacuum_vessel + thermal_shield`. The
group's actual layer sum (`coil_standoff.py` stack SOL+W+steel+Be+FLiBe+shield; matching the ParaStell
radial build first_wall/breeder/back_wall/shield/vacuum_vessel in `PARASTELL_SIMSOPT_CAPABILITY.md`) is
**≈ 100 cm**. So the **hard cut is: min coil-plasma standoff ≥ ~1 m at reactor scale** — i.e. at R₀=10 m,
**d_min/R₀ ≥ 0.10**. This is cleaner than any margin heuristic: it is literally "can you house the
blanket+shield at all." Below it, the device is disqualified regardless of quasisymmetry. **Above it,
the *excess* standoff (standoff − radial-build) is the neutronics-relevant knob** — it is the available
*extra* shield thickness / headroom that sets how much you can locally re-route shield to leak paths
(the placement lever in `EXECUTION_PLAN` ADJ-1). So standoff enters twice: as a **binary fit/no-fit
gate** (this filter) and, beyond the minimum, as the **continuous shield-headroom predictor** the law
correlates against.

Reactor-practice grounding (cited): ARIES-CS requires **1.5–2 m plasma-to-coil-winding-pack** midplane
spacing for first-wall+blanket+shield, with a **min standoff ~1.3 m** achieved by swapping in a WC-based
shield at the critical minimum (which *cut* the required standoff **20–30%**) — i.e. ARIES-CS is a
~1.3–2 m radial-build device, consistent with the ~1 m floor here being a *minimum*, not a comfortable
build. QUASR bakes in engineering constraints already: **min coil-surface & coil-coil distance 0.1 m,
max curvature 5 m⁻¹, max mean-squared curvature 35 m⁻²** at the dimensionless R₀=1 m scale — so the
0.1 m@R₀=1 m min-distance ≈ **1 m@R₀=10 m ≈ exactly the radial-build floor**, which is why (see §2.3) the
fit/no-fit gate barely bites on QUASR.

| Filter | Compute from | Threshold (suggested) | Grounding |
|---|---|---|---|
| **Min coil-plasma standoff ≥ radial build** (PRIMARY, fit/no-fit) | filaments + LCFS (group's `coil_standoff.py`, KD-tree) | at R₀=10 m, **gap ≥ 100 cm** = full radial build; i.e. **d_min/R₀ ≥ 0.10**; *excess* standoff beyond this = shield-headroom knob | group radial build ≈1 m; ARIES-CS 1.3–2 m |
| **Aspect ratio band** | metadata | **4 ≤ A ≤ 10** (compact-reactor band; ARIES-CS A ≤ 6) | ARIES-CS A ≤ 6; HELIAS ~10 |
| **Max coil curvature (buildability)** | filaments (Frenet κ) | κ·R₀ ≤ 5 m⁻¹ → at R₀=10 m, **κ ≤ 0.5 m⁻¹ (R_curv ≥ 2 m)** | QUASR built-in cap 5 m⁻¹ @R₀=1 m |
| **Coil-coil min distance / R₀** | filaments (KD-tree) | **≥ 0.10** (≥ 1 m @R₀=10 m) for winding-pack + structure | QUASR built-in 0.1 m; stage-2 standard penalty |
| **Coil length per period** | metadata (`coil_length_threshold_per_hp`) | prefer shorter for cost/complexity; use as a tiebreak, not a hard cut | Giuliani et al. complexity penalties |

Because QUASR **already enforces** curvature and 0.1 m distances, the curvature and coil-coil filters are
**nearly non-binding on QUASR** (they *will* bind on self-generated coils — apply them there). The
**binding** cuts for QUASR are **standoff** and **aspect ratio**.

### 2.3 Surviving-fraction estimate — grounded in the group's real data

Computed on the group's **cached 366-device standoff table**
(`spf_prototype/data/quasr/coil_standoff.csv`; 157 QA + 209 QH, nfp 1–8). Distribution:
`d_min/R₀` median **0.166** (min 0.007, max 0.418); aspect median **12** (min 2.9, max 24).

| Cut | Surviving | Fraction |
|---|---|---|
| **radial-build fit: gap ≥ 100 cm @R₀=10 m (PRIMARY, standoff alone)** | 348/366 | **95%** |
| gap ≥ 120 cm (build + small headroom) | 334/366 | 91% |
| gap ≥ 150 cm (build + generous shield headroom) | 256/366 | 70% |
| radial-build fit **AND** 4 ≤ A ≤ 10 | 102/366 | **28%** |
| radial-build fit **AND** 5 ≤ A ≤ 8 | 55/366 | **15%** |

**Key finding (honest):** the **radial-build fit/no-fit cut alone barely bites** (~**95%** pass) because
QUASR's baked-in 0.1 m min-distance constraint at R₀=1 m ≈ **1 m reactor gap** — i.e. QUASR devices are
*already* approximately radial-build-clearing by construction. Only ~5% fail to house the ~1 m stack.
Requiring **shield headroom beyond the minimum** (gap ≥ 150 cm) cuts harder (→70%). The **binding
constraint is aspect ratio** (this sample's median A=12 is far from compact-reactor relevance). So the
realistic surviving fraction under **radial-build-fit + compact-aspect** lands at **~15–30%** of QUASR.
**Caveat:** the 366-device
cache is a *sub-sample* and its aspect distribution (median 12) is almost certainly not representative of
the full 370k; treat 15–30% as an order-of-magnitude, and recompute on the full `QUASR.pkl` before
quoting a number in the paper. The **cut is cheap** — standoff from filaments the group already parses
(`coil_standoff.py`, O(N log N) KD-tree), aspect straight from metadata; the whole 370k is a metadata
filter + a KD-tree pass on the survivors.

---

## 3. Controlled-standoff experiment design (the confound-breaker)

**Goal.** Separate *plasma-shape* effect from *standoff* effect on the coil-load concentration/shieldability
metric — the confound flagged in `SCIENCE_STRATEGY.md` §1.1 (QA vs QH differ ~1.9× in standoff, so a 2-point
concentration→benefit trend is inseparable from standoff→everything).

**Design = a 2-way grid: {plasma shape} × {standoff}.**

- **Plasma shapes (fixed boundaries): 3** — one QA + one QH (from QUASR, coils re-derivable as
  ground-truth) + one QI (from ConStellaration, coils *must* be generated — the Step-4 spike). Three
  spans the QS-type axis and includes the reactor-relevant QI point.
- **Standoffs per shape: 3** — `standoff/a ∈ {tight, nominal, generous}`, e.g. **{0.7, 1.0, 1.4}** in
  minor-radii (tight ≈ ARIES-CS-aggressive with local WC shield; nominal ≈ the 1.01 m group stack;
  generous ≈ comfortable blanket). Generate with simsopt stage-2 at three `CurveSurfaceDistance` targets,
  warm-started tight←nominal←generous.
- **Total: 3 × 3 = 9 (plasma, standoff) coil sets.** This is the minimum for a **clean 2-way separation**:
  ≥3 standoffs per shape gives a within-shape standoff *slope* (not just a difference), and ≥3 shapes at
  each standoff gives a shape *contrast* at fixed standoff. 9 points let you fit
  `metric ~ shape + standoff` and read off each main effect. (A 2×2 would give differences but no slope and
  no test of monotonicity; 3×3 is the smallest grid that separates the two axes with a shape of curve.)
  If budget forces it, a **2 shapes × 3 standoffs = 6** still separates the axes but weakens the
  shape-generality claim.

**Cheap vs expensive split.**

- **Cheap (free-streaming proxy, NO DAGMC) — run on all 9 first:** per-coil free-streaming flux
  concentration (Gini/PR via committed `concentration.py`) and peak, with the group's existing
  ring/volumetric plasma source streaming to the parsed filaments. This is SCIENCE_STRATEGY Step 2's
  machinery. It **cheaply shows** whether concentration moves with standoff at fixed shape and with shape
  at fixed standoff — the whole 2-way question — for essentially free (no new DAGMC builds). What it
  **cannot** show: through-shield magnet dose, leak-around saturation, or the *shieldability* (placed-vs-
  uniform) benefit — those are not free-streaming quantities.
- **Expensive (DAGMC + adjoint) — run on a 2–3-point subset only:** the through-shield contributon
  concentration + the transport placed-vs-uniform benefit (Step-1 kill-or-confirm) needs a ParaStell→DAGMC
  build per (plasma, standoff) point (days each per `EXECUTION_PLAN` Phase B evidence). Do it on the
  **diagonal** (tight/QA, nominal/QH, generous/QI) or on one shape's full standoff column, gated by whether
  the cheap proxy already shows a standoff effect.

**Gate (pass/fail).** The experiment *succeeds as a design* (regardless of outcome sign) if it can answer
**both**:
1. **Does concentration vary with standoff at fixed shape?** — within-shape slope of concentration vs
   standoff/a, with CI, per shape. Non-zero, consistent-sign slope ⇒ standoff is a real driver.
2. **Does concentration vary with shape at fixed standoff?** — spread across the 3 shapes at each standoff,
   with CI. Non-zero spread ⇒ shape carries signal *beyond* standoff.

The scientifically decisive outcomes: if (1) is strong and (2) is null → the "geometry→concentration" law
is really "standoff→concentration" (confound confirmed, law demoted to a standoff statement — a clean,
publishable negative). If (2) survives at fixed standoff → shape carries independent signal and the
relational-descriptor thesis (`DATABASES_AND_RELATIONAL_DESCRIPTORS.md`) holds. Either way the confound is
*measured*, not lurking. **Honesty gate:** report both slopes with CIs even if one is null; do not tune the
standoff grid to manufacture a clean separation.

---

## 4. Parallelized compute budget on Ginsburg (reaching n ≈ 20–50)

The whole pipeline is **embarrassingly parallel across devices** — each (plasma, standoff) point is
independent — so the question is not per-device speed but **how many devices a SLURM job array clears per
wall-clock window**, and **which stage is the bottleneck**.

### 4.1 Per-device cost, broken down

| Stage | Tool | Wall time (single device) | Parallelism | Notes |
|---|---|---|---|---|
| (a) Coil design | simsopt stage-2 (L-BFGS on ~10² DOFs) | **~2–15 min** | multi-core Biot-Savart per task; also parallel *across* devices | REGCOIL alternative is **~seconds** (linear solve) — coil design is **never** the bottleneck |
| (b) ParaStell → DAGMC build | **`cad_to_dagmc`** (CadQuery/OCC + MOAB; **Cubit-free, license-free** — the group's verified path) | **~30–60 min** | **poorly threaded**: CAD lofting + imprint/merge/facet is largely serial (OCC single-thread); MOAB partly threads | **THE BOTTLENECK.** The group's `rebuild_dagmc.sbatch` budgets `-t 0-01:00` on `-N1 -c16 --mem-per-cpu=4G`, i.e. ~1 h/16-core envelope per build |
| (c) Adjoint solve | random-ray / MC adjoint + VR | **~10 min** on a node (reactivity-independent → solved once, re-weighted free) | threads well across cores | cheap; runs after the `.h5m` exists |
| **Total** | | **~45–90 min/device**, DAGMC-dominated | | |

**Bottleneck = the DAGMC build (b), not coil design.** Coil design (a) is minutes; the adjoint (c) is
minutes; the CAD→watertight-`.h5m` faceting is the tens-of-minutes, weakly-threaded step. Crucially the
group uses **`cad_to_dagmc` (license-free)**, so **there is no Coreform Cubit seat contention** — the array
can run as wide as the cluster allocation allows. (If they were on the Cubit path, concurrent tasks would
be capped by license seats — a real reason to stay Cubit-free for a bulk campaign.) The honest per-device
risk is not time but **`cad_to_dagmc` robustness on thin, highly-shaped shield layers** (imprint/merge of
near-coincident surfaces — flagged in `PARASTELL_SIMSOPT_CAPABILITY.md`): a failed facet is a *retry*, not
a slow success, and should be caught by a watertightness gate before the adjoint stage.

### 4.2 Throughput as a SLURM job array

Model: one device per array task; DAGMC build sets the ~1 h/task wall time; `T` = number of concurrent
tasks the allocation sustains. Devices cleared ≈ `T × (window / 1 h)`.

| Window | With T = 25 concurrent tasks | With T = 50 |
|---|---|---|
| **A few hours** (~3 h) | ~75 device-slots → **n ≈ 20–40 realistic** (after coil-design + watertight gating losses) | ~150 slots → n ≈ 40–60 |
| **Overnight** (~12 h) | ~300 slots → **n ≈ 20–50 trivially cleared**; bounded by coil/gating, not compute | ~600 slots |
| **A weekend** (~60 h) | ~1500 slots → **hundreds** | thousands |

**Verdict: compute is NOT the binding constraint for n ≈ 20–50.** Even a modest array (T≈25) clears the
n a defensible law needs **overnight**, and a few-hours window already reaches n≈20–40. What actually
limits n is (1) **coil-design convergence + the hard-gates** (tight-standoff failures, curvature/coil-coil
rejects) and (2) **`cad_to_dagmc` faceting robustness** (watertight-`.h5m` retries) — both are
success-rate, not throughput, problems.

### 4.3 Recommended array design

- **Pre-step (single small job):** generate all coil sets — a cheap serial/loop job over the filtered
  device list producing MAKEGRID `.coils` per (plasma, standoff), each hard-gated on achieved-standoff,
  max-curvature, min-coil-coil. (Minutes-scale; do *not* burn array slots on this.)
- **Array 1 — DAGMC builds (the bottleneck):**
  `sbatch --array=0-49%25 -p short -N1 -c8 --mem-per-cpu=4G -t 0-02:00 build_dagmc_device.sbatch`
  — 50 devices, **25 concurrent**, 2 h wall each (generous vs the 1 h envelope for retries), 8 cores/task
  (memory-bound faceting; more cores buy little given weak threading). Emit a **watertightness gate** per
  task; log failures for retry rather than blocking.
- **Array 2 — adjoint solves (dependent):** `--array=... --dependency=aftercorr:<array1>` so each adjoint
  starts when its `.h5m` lands; `-c 16 -t 0-00:30`. Reactivity-independence means one solve per device
  serves all emissivity re-weights for free.
- **Scale knob:** raise `%25 → %50` if the `short` partition allocation allows; throughput scales linearly
  in concurrency until you hit the queue/QOS cap. **n ≈ 20–50 is an overnight run**; n ≈ 100+ is a weekend.

---

## Sources

**Tools / availability (probed on-machine, base conda simsopt 1.10.6):** simsopt stage-2 classes all import
(`BiotSavart, Coil, Current, coils_via_symmetries, CurveXYZFourier, create_equally_spaced_curves,
SurfaceRZFourier, CurveLength, CurveCurveDistance, CurveSurfaceDistance, MeanSquaredCurvature, SquaredFlux`);
`simsopt.field.CurrentPotential*` **absent** (REGCOIL not in build);
`simsopt.configs.download_ID_from_QUASR_database` present. Group recipe on disk:
`shield_opt/data/raw/QH_unzipped/zenodo/optimization_scripts/stage_two_optimization_planar_coils_iterate.py`;
group standoff code `sweep/coil_standoff.py`; cache `data/quasr/coil_standoff.csv` (366 devices).
- simsopt — https://github.com/hiddenSymmetries/simsopt ; coil example
  https://simsopt.readthedocs.io/en/stable/example_coils.html (`CurveSurfaceDistance`, `CurveCurveDistance`,
  stage-2 objective).
- DESC coil tools (REGCOIL + stage-2, pure-Python) —
  https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/coil_optimization_REGCOIL.html ;
  .../coil_stage_two_optimization.html
- REGCOIL method — Landreman & Hanson, *Nucl. Fusion* 57, 046003 (2017), arXiv:1609.04378;
  surface-current / coil-cutting — arXiv:2508.09321.

**QUASR access + constraints:**
- Giuliani et al., *J. Plasma Phys.* 2025, arXiv:2409.04826 (engineering constraints: min coil-surface &
  coil-coil 0.1 m, max curvature 5 m⁻¹, max MSC 35 m⁻²); QA-only method arXiv:2310.19097.
- Navigator/metadata — https://quasr.flatironinstitute.org/ ; `QUASR.pkl` metadata dataframe;
  `simsopt_serials/{ID[:4]}/serial{ID}.json`; Zenodo `10.5281/zenodo.10050655` / `10.5281/zenodo.13717741`
  (~12.7 GB; coil-length-threshold erratum noted in the deposit).

**Reactor coil-plasma thresholds:**
- ARIES-CS — Najmabadi et al., IAEA FEC2006 FT/P5-26 (1.5–2 m plasma–coil; min ~1.3 m; WC shield at min
  standoff cut it 20–30%); https://www-pub.iaea.org/MTCD/Meetings/FEC2006/ft_p5-26.pdf
- ParaStell input (VMEC + coil filament point-locus / MAKEGRID `.coils`) — Frontiers Nucl. Eng. 2024,
  10.3389/fnuen.2024.1384788 ; https://github.com/svalinn/parastell ; MAKEGRID
  https://princetonuniversity.github.io/STELLOPT/MAKEGRID.html
- Compute budget grounded on-machine: group DAGMC path is **Cubit-free `cad_to_dagmc`** and radial build
  (first_wall/breeder/back_wall/shield/vacuum_vessel) — `spf_prototype/talk/PARASTELL_SIMSOPT_CAPABILITY.md`
  (VERIFIED: Cubit optional; layers); SLURM envelope from `spf_prototype/ginsburg_jobs/rebuild_dagmc.sbatch`
  (`-t 0-01:00 -N1 -c16 --mem-per-cpu=4G`, partition `short`); radial-build ~1 m layer sum from
  `sweep/coil_standoff.py`. Adjoint reactivity-independence — `EXECUTION_PLAN_magnet_adjoint.md`.
- ConStellaration (QI boundaries, coils self-generated) — arXiv:2506.19583.

*Uncertain / not verified here: full-370k surviving fraction (computed only on a 366-device sub-sample —
recompute on `QUASR.pkl`); DESC install cost on the group's stack (not installed, not tried); whether
`CurveSurfaceDistance`'s min-distance target yields a standoff uniform enough for a clean 2-way separation
(the §0 risk — measure achieved standoff per coil set).*
