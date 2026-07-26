# Extended science vision — the L∇B → standoff → magnet-load causal chain, and the neutronics evaluator for single-stage stellarator design

*Written 2026-07-25. A critical, grounded multi-paper vision that folds the Kappel–Landreman–Malhotra
magnetic-gradient-scale-length result (PPCF 66 (2024) 025018) into the group's built adjoint
magnet-shielding infrastructure. Grounds: `CONTEXT_RECOVERY_tonight.md`, `SCIENCE_STRATEGY.md`,
`DATABASES_AND_RELATIONAL_DESCRIPTORS.md`, `COIL_DESIGN_AND_FILTERING.md`,
`EXECUTION_PLAN_magnet_adjoint.md`, `paper/PAPER_DIRECTION.md`, `paper/RESULTS_adjoint.md`. This is a
strategy/vision document — WRITING/RESEARCH ONLY, no code, no transport. It does not supersede
`EXECUTION_PLAN_magnet_adjoint.md`; it sits above it as the "why" and the "where this goes" and should
be reconciled into the roadmap. Honesty rules from the group style apply throughout: n=2 is n=2;
surrogate ≠ transport until Step 1 closes; label every device by its standoff/L∇B, not by QS folklore;
where the Kappel connection is a hopeful analogy rather than a validated link, it is flagged as such.*

---

## 0. Report-back (read this first)

**Does the L∇B chain make the predictive law genuinely reachable, vs the earlier "over-reached at
n=2"?** *Partially, and asymmetrically — it is a real upgrade, not a rescue.* The thesis is a two-link
causal chain. **Link 1 (L∇B → minimum plasma–coil separation) is already validated by Kappel at n>40,
R²=0.94 — and, decisively, the group's own two devices are IN Kappel's database** (the Landreman–Paul
precise QH and precise QA, Appendix D). So the "plasma geometry → standoff" half is no longer a
two-point line the group must defend; it is a 40-device published law the group inherits, with QH/QA as
named points on it. **Link 2 (separation → available shield thickness → magnet neutron load /
shieldability) is the group's to establish, and it is still n-limited** to the handful of DAGMC+adjoint
devices the pipeline can build. So the chain converts a *fully confounded* n=2 correlation into *one
validated upstream link plus one link the group must earn at modest n*. That is genuinely reachable as
"an illustrative predictive law with a mechanistically grounded, independently validated upstream
cause" — **not** a tight universal law, and **still gated** on the Step-1 transport kill-shot and on
the ParaStell→DAGMC→adjoint automation.

**Single biggest risk.** Link 2 turns out to be *trivial or absent*: magnet load tracks standoff by
plain 1/r² geometric dilution and collimation, so the adjoint attribution earns nothing on top of
Kappel's geometry ("farther coils see less flux" needs no adjoint); OR the group's own leak-around
**saturation** finding means added shield room past the radial-build floor doesn't help, so "available
shield thickness" is the wrong operative variable. Either way the neutronics does not justify itself
over the geometric proxy, and the chain collapses to Kappel + inverse-square. This is the Step-1 null
in a new guise, and it must be probed before any n≈5 QI investment.

**Sharpest new experiment the chain enables.** The **controlled-standoff experiment at fixed plasma /
fixed L∇B**: hold one plasma boundary, use simsopt stage-2 (`CurveSurfaceDistance`) to generate coil
sets at 3 standoffs, and measure whether magnet load *and its attribution concentration* track standoff
as predicted — **decomposing the response into geometric dilution vs shield attenuation vs
leak-around.** This isolates the mediator (standoff) from its root cause (L∇B), and it is the one test
that shows whether the neutronic attribution adds signal beyond 1/r². A cheaper, equally sharp variant:
the **spatial-coincidence test** — Kappel shows L∇B is minimum on the *inside of the bean* cross
section; does the adjoint magnet-load attribution hotspot land there too? If yes, one figure ties root
cause to neutronic outcome; if no, the "same geometry drives both" unification weakens.

**Recommended Paper-1 scope.** The methods + first-comparison paper, essentially as
`PAPER_DIRECTION.md` already frames it, with **one reframe**: label the QA-vs-QH comparison by
**standoff / L∇B, citing Kappel, and note both devices sit in his 40-config database** — this dissolves
the standoff "confound" into a named, root-caused axis instead of an apology. Deliverables:
(i) adjoint coil-kerma attribution + reciprocity (QH +0.76, QA +0.83); (ii) the Step-1 transport
placed-vs-uniform kill-shot (confirm or null — both publishable); (iii) the zoo-wide free-streaming
coil-concentration law as the real-n companion; (iv) the L∇B→standoff→magnet-load chain stated as the
*organizing hypothesis*, with P2 as its test. Claim the method + the neutronics-evaluator framing +
the first honest comparison + the free-streaming law. **Do not claim the through-shield predictive
law** in P1.

---

## 1. The unified thesis, sharpened

### 1.1 The causal chain, stated precisely

```
   L∇B                    min plasma–coil          available shield        magnet neutron load /
(plasma field   ──────▶   separation        ─────▶ thickness        ─────▶ shieldability
 geometry;                [Kappel R²=0.94,          (= separation −         [group adjoint attribution;
 cheap, local,             slope 1.585, n>40]        fixed radial build]     reciprocity-validated
 reactivity-indep]         LINK 1 (validated)        LINK 2a                 Spearman +0.76 QH / +0.83 QA]
                                                                             LINK 2b (to establish)
```

- **L∇B** (Kappel eq. 12): `L∇B = √2·B_coils / ‖∇B_coils‖_F`, the Frobenius norm of the 3×3 gradient
  tensor of the *coil* field, evaluated on the LCFS; `L*∇B = min over the LCFS` (eq. 27). It is
  **cheap** (local `B` + `∇B` at boundary points via virtual-casing in simsopt), **reactivity- and
  coil-design-independent** (needs only the plasma field geometry, no coil set), and it has a clean
  physical meaning: for an infinite straight wire `L∇B = R`, the distance to the wire (eq. 17); the √2
  is chosen to make that identity exact. So L∇B ≈ "how far can a source of this field be" — a *local*
  read-out of the *global* difficulty of placing coils far from this plasma.

- **Link 1 — L∇B → separation — is Kappel's result and it is validated:** across >40 reactor-scaled
  configs (QA, QH, QI, tokamaks, misc; all scaled to ARIES-CS a=1.704 m, B=5.865 T), L*∇B predicts the
  REGCOIL minimum achievable plasma–coil separation with **R²=0.941, slope 1.585**, insensitive to the
  `‖K‖∞` / `B_RMS` constraints. L_REGCOIL varies by ~10× across the database. Kappel's own headline
  reactor-relevant trend: **large separation is available only at small field-period number** (the nine
  most-separated configs are all nfp 1–2).

- **Link 2a — separation → available shield thickness — is near-tautological but bounded.** The
  reactor radial build (first wall + multiplier + breeder + back wall + shield + gap + vessel) is
  ≈1 m; separation beyond that floor is the *shield-headroom* knob. The group's own `COIL_DESIGN`
  survey shows the fit/no-fit cut barely bites on QUASR (~95% clear the ~1 m stack, because QUASR bakes
  in a 0.1 m@R₀=1 m min-distance ≈ 1 m at reactor scale); the *excess* standoff is where the physics
  lives.

- **Link 2b — shield thickness / standoff → magnet load & shieldability — is the group's contribution,
  and it is the only link not yet established.** This is exactly where the adjoint coil-kerma
  attribution map, reciprocity (+0.76 QH / +0.83 QA), and the closed-loop placement lever live.

**The one-sentence meta-thesis:** *Kappel stops at geometric clearance — can the blanket physically
fit between plasma and coils. The group adds the neutronics — does the blanket, in the room that L∇B
permits, actually protect the magnets. L∇B sets the ceiling on protectability; the coil design and the
neutron transport determine how much of that ceiling is realized.*

### 1.2 How Kappel dissolves the standoff "confound"

`SCIENCE_STRATEGY.md §1.1` correctly flagged the QA-vs-QH comparison as fatally confounded: the two
devices differ mostly in coil standoff (1.9×), which plausibly drives *both* attribution concentration
*and* shieldability, so a 2-point concentration→benefit correlation cannot be separated from
standoff→everything. The critique treated standoff as a **lurking common cause** sitting outside the
model.

Kappel reframes standoff from a *lurking variable* to the **mediator on a causal path with a measured,
independent root cause**. In DAG terms: it was `standoff → {concentration, benefit}` with standoff
unexplained; it is now `L∇B → standoff → {load, shieldability}`, with L∇B computable from the plasma
alone and validated as the driver of standoff at n=40. Three consequences:

1. **Standoff is no longer a nuisance to apologize for; it is a named axis with a root cause.** Label
   every device by (L∇B, standoff), not by "QA vs QH." The QA-vs-QH gap *is* the L∇B axis, not a
   quasisymmetry effect — which is exactly the correction `PAPER_DIRECTION.md` already insists on ("the
   QA-vs-QH magnet gap is a coil-standoff confound, not a quasisymmetry advantage"). Kappel supplies
   the upstream variable that makes that statement mechanistic instead of merely cautionary.

2. **The group's two devices are already on Kappel's validated line.** From Appendix D (reactor-scaled):
   Landreman–Paul QH (nfp 4) has **L_REGCOIL ≈ 1.52 m, L*∇B ≈ 3.19 m**; precise QA (nfp 2) has
   **L_REGCOIL ≈ 2.87 m, L*∇B ≈ 5.30 m**. The group's *measured as-built* standoffs — QH 1.63 m, QA
   3.11 m (ratio 1.91) — sit right on Kappel's predicted separations (ratio 1.89), and L∇B correctly
   orders QA above QH by a large margin. **This means Link 1 for QH/QA requires zero new computation —
   it is a table lookup.** (Honest caveat: these two points preserve ordering and the ~1.9× ratio but
   are not perfectly on the R²=0.94 regression line — both sit somewhat below it; L∇B's *ratio* 1.66
   slightly under-predicts the measured standoff ratio 1.91. Cite the ordering and the "large gap,"
   not a spuriously exact number.)

3. **What the confound-dissolution does NOT do.** It does not establish Link 2b. Proving L∇B drives
   standoff says nothing about whether standoff drives magnet load in a way the adjoint captures beyond
   trivial 1/r². And at *fixed* L∇B, whether residual plasma shape or coil-design choices still move
   magnet load is an open question — that is precisely what the controlled-standoff experiment (§4.2)
   tests. Kappel removes the confound from the *upstream* half; the *downstream* half still needs the
   group's own controlled experiment.

### 1.3 What would falsify each link (be explicit)

| Link | Status | What falsifies it | Cost to probe |
|---|---|---|---|
| **1. L∇B → separation** | **Validated (Kappel, n>40, R²=0.94); QH/QA in the table** | L∇B computed across the group's target device set fails to correlate with the group's *achieved* standoffs; or QH/QA fall far off the line once re-derived in the group's frame | Cheap — read Kappel's table; recompute L∇B from VMEC `B+∇B` for QUASR devices |
| **2a. separation → shield room** | Near-tautological, bounded by radial build | Leak-around makes shield room irrelevant past the floor (group's own saturation finding) → "available thickness" is the wrong variable | Already half-shown (ADJ-1 saturation); Step-1 confirms |
| **2b-i. standoff → magnet load** | **To establish** | Magnet load tracks standoff by pure 1/r² dilution, so attribution adds nothing over Kappel + inverse-square (the neutronics doesn't earn its keep) | Controlled-standoff experiment (§4.2) |
| **2b-ii. concentration → shieldability** | **To establish** | Placed shield ≈ uniform shield even at high attribution concentration (leak-around irreducible) → "concentration ⇒ shieldable" is false | **Step-1 kill-shot** (already the top priority) |
| **Unification: same geometry drives both** | Hopeful | L∇B's bean-inside minimum does NOT coincide with the adjoint attribution hotspot; magnet-side concentration keys off coil descriptors independent of the plasma-shape descriptors that drive the first wall | Spatial-coincidence test (§4.4), cheap |

**Net:** the chain is now *grounded* (mechanism + one validated link) rather than *spurious* (a
confounded 2-point line). It is not yet *proven end-to-end*. Link 2b is the whole scientific bet, and
its two failure modes (trivial 1/r²; irreducible leak-around) are both cheap to probe and must be
probed first.

---

## 2. The neutronics-evaluator framing, honestly bounded

### 2.1 The gap the group fills

The stellarator-design field is shifting away from two-stage optimization toward **single-stage
differentiable optimization** (simsopt / DESC) and **AI-surrogate / generative design**
(StellFoundry-style, ML challenges on ConStellaration). To make that tractable, the field is adopting
cheap **engineering proxies** as objective terms or filters: plasma–coil clearance, coil curvature,
coil–coil distance, forces — and now Kappel's **L∇B** as a differentiable stand-in for "how much
separation is achievable." **Every one of these proxies is geometric/structural. None is neutronic.**
The 1.3–1.5 m clearance that a proxy enforces is a stand-in for "the magnets are protected" that is
*never actually verified against neutron transport*.

**The group's adjoint magnet-shieldability pipeline is the missing neutronic engineering proxy.** L∇B
tells you the blanket can *fit*; the adjoint tells you whether, in that space, the magnets are
*protected* — and *where* the load concentrates so you can re-route shield. This is the neutronic
ground truth that every geometric proxy is silently standing in for.

### 2.2 Evaluator/filter (reachable now) vs differentiable objective term (the named lift)

This distinction is load-bearing and must not be blurred:

- **EVALUATOR / filter — reachable now.** Given a device (boundary + coils), the
  ParaStell→DAGMC→adjoint pipeline returns: magnet neutron load, its phase-space attribution
  (contributon), and shieldability (placed-vs-uniform benefit). This is a *scoring function* that can
  be dropped into a generative or single-stage loop as a **ranker/rejecter**: score candidates by
  actual neutronic protectability instead of by a geometric clearance proxy, or use it to **calibrate
  and validate the cheap proxies** the loop optimizes against. The reactivity-independence of ψ†
  (solve once, re-weight the source for free) makes the evaluator cheap to sweep over scenarios. This
  is the honest, immediate contribution.

- **DIFFERENTIABLE objective term — the future lift, name it plainly.** To put magnet-shieldability
  *inside* a gradient-based single-stage optimizer as an objective term, you need
  `∂(magnet load)/∂(boundary Fourier coefficient or coil DOF)` — a **transport shape gradient**
  (differentiate the transport solve). The reactivity-independence trick does **not** provide this: it
  re-weights the *source*, not the *geometry*; it gives no geometric derivative (this is the same
  conclusion `SCIENCE_STRATEGY.md §4(b)` reached for adjoint-shape steering). Building it is a whole
  adjoint-shape-optimization program overlapping groups whose comparative advantage *is* shape
  optimization. **Do not promise it.** Frame it as: the contributon is the neutronics half of a coupled
  shape gradient; wiring it to a differentiable equilibrium/coil model is future work.

### 2.3 Positioning against the differentiable-L∇B proposal and single-stage optimization

- **vs Kappel's differentiable-L∇B (their §5 proposal).** Kappel explicitly proposes L∇B as a
  differentiable stage-1 objective to *maximize* plasma–coil separation, precisely because it is cheap,
  local, and easy to differentiate — a *geometric proxy for shield room*. The group's relationship to
  this is complementary and honest: **we are the neutronic ground truth that L∇B is a proxy for.** The
  strongest near-term contribution is not to compete with L∇B as a differentiable term but to
  **calibrate L∇B (and standoff, conformity, streaming solid angle) against actual magnet load via the
  adjoint evaluator** — so the single-stage community can keep optimizing the cheap differentiable
  proxy *knowing what neutronic quantity it does and does not predict*. That reframes "we can't
  differentiate transport yet" from a weakness into the correct division of labor.

- **vs single-stage (simsopt/DESC) and generative (StellFoundry / STELLAR-AI / ConStellaration ML).**
  These optimize plasma+coils jointly against geometric/confinement objectives. Their neutronics is,
  at best, a clearance constraint. The group provides the evaluator that turns "clearance satisfied"
  into "magnets verified protected, and here is where the load concentrates." Prior-art hygiene
  (already in the group docs): differentiate from ParaStell's brute-force 3-D shield/TBR sweeps
  (*adjoint says WHERE*), from HELIAS/ADVANTG adjoint variance-reduction (*attribution, not VR*), and
  from arXiv:2604.26763 coil-nonplanarity↔plasma-surface (*stay on the neutron-loading axis*).

---

## 3. The multi-paper arc

Three papers, in dependency order. Each states its claim, the data/experiments it needs, its gate, its
prior-art differentiation, and honest reachability.

### P1 — Methods + first cross-device comparison (reachable THIS cycle)

- **Claim.** (a) An adjoint coil-kerma flux, solved with OpenMC random-ray adjoint, is a coil→plasma
  *attribution* map that drives a closed-loop breeder-for-shield lever — validated by reciprocity on
  two devices of different QS type (QH +0.76, QA +0.83). (b) In transport, adjoint-placed shield does /
  does not beat equal-budget uniform shield (Step-1). (c) A zoo-wide free-streaming coil-concentration
  law relates load structure to coil descriptors at real n. (d) The QA-vs-QH magnet difference is a
  standoff/L∇B effect, not a quasisymmetry effect. **Framing, not a new claim:** the
  L∇B→standoff→magnet-load chain is the organizing hypothesis; cite Kappel; note both devices are in
  his database.
- **Data / experiments.** All in hand or in the execution plan: QH+QA DAGMC+adjoint (built, validated);
  the Step-1 placed-vs-uniform transport run (top priority, `EXECUTION_PLAN` Phase C/Step 1); the
  free-streaming QUASR zoo (Step 2, no new DAGMC); SPF as supplementary (built).
- **Gate.** Step-1 returns a transport result (confirm: placed beats uniform, sign tracks
  concentration; OR null: leak-around irreducible). Both are publishable; the null becomes "placement
  can't beat leak-around," still a real finding.
- **Prior-art differentiation.** ParaStell (WHERE, not brute force); ADVANTG (attribution, not VR);
  Bae 2025 (SPF is native/on-the-fly/supplementary, not the headline); 2604.26763 (neutron loading,
  not plasma-surface geometry).
- **Reachability.** High. This is the MVP of `SCIENCE_STRATEGY.md` and `PAPER_DIRECTION.md`, plus the
  L∇B *framing* which costs a table lookup and a paragraph. **No new physics required to write it.**

### P2 — The L∇B → standoff → magnet-protectability predictive law (reachable IF Step-1 confirms + pipeline automates)

- **Claim.** Across a reactor-scaled device set, plasma-field geometry (L∇B) predicts achievable
  plasma–coil separation (Kappel, inherited), and separation predicts magnet-load protectability
  (group), with L∇B as the spine and the adjoint attribution as the mechanism that explains *how much
  of the L∇B-permitted room actually protects the magnets*. Reactor-relevant corollary: **small-nfp
  configs, which L∇B/Kappel show permit large separation, are intrinsically more magnet-protectable.**
- **Data / experiments.** (i) L∇B for the Kappel/Zenodo-8349408 ~40 configs (read the table) + compute
  from VMEC `B+∇B` for a QUASR sub-sample and for QI (ConStellaration boundaries + self-generated
  coils). (ii) Free-streaming per-coil concentration across the whole set (cheap, real n). (iii)
  DAGMC+adjoint magnet load on a standoff-spanning subset (n≈5–7, automation-gated). (iv) The
  controlled-standoff experiment (§4.2) as the confound-breaker. (v) The nfp→protectability test
  (§4.3).
- **Gate.** Step-1 must have confirmed the mechanism (concentration ⇒ shieldable) in P1; the
  ParaStell→DAGMC→adjoint per-device cost must be affordable (QI feasibility spike, `SCIENCE_STRATEGY`
  Step 4); the standoff-controlled slope of magnet load must be non-trivial (not pure 1/r²).
- **Prior-art differentiation.** Kappel (we add the *neutronic outcome* to his *geometric clearance*);
  Miralles-Dolz symbolic-regression neutronics surrogates (ours is *mechanistic + adjoint-attribution*,
  not black-box); 2604.26763 (neutron loading axis). The ConStellaration ML challenge is filling the
  geometry→neutronics surrogate space fast — P2's defensible edge is the *mechanism* (L∇B as validated
  root cause + adjoint as the mediator's neutronic read-out), not a fitted surrogate.
- **Reachability.** Conditional. The L∇B→separation half is n=40 (borrowed, strong). The
  separation→magnet-load half is n≈5–7 DAGMC = **illustrative trend + hypothesis, not a tight law**
  (be explicit; below n≈20 no ξ is a law). Honest landing: "a mechanistically grounded predictive
  *framework* with one validated upstream link and a first neutronic test of the downstream link across
  a few standoff-spanning devices." Still a strong Nuclear Fusion paper if Step-1 is positive.

### P3 / vision — The neutronics evaluator for single-stage / generative stellarator design

- **Claim.** Automated ParaStell→DAGMC→adjoint magnet-shieldability is the missing *neutronic*
  engineering proxy for single-stage and generative stellarator optimization: it converts the field's
  geometric clearance proxies (L∇B, standoff, conformity) into verified magnet protectability + a
  targeting map, and it calibrates which cheap differentiable proxies are neutronically trustworthy.
- **Data / experiments.** The automation itself (a repeatable per-device job) is the infra
  contribution; a demonstration coupling the evaluator into a candidate-ranking loop over a generative
  boundary set; a calibration curve of L∇B / standoff vs actual adjoint magnet load across the P2 set.
- **Gate.** P2's pipeline automates to n≈20–50 overnight (the `COIL_DESIGN` compute analysis says
  compute is not the binding constraint — coil-design convergence and watertight-`.h5m` success rate
  are); the evaluator demonstrably changes a ranking that a geometric proxy alone would get wrong.
- **Prior-art differentiation.** The differentiable-objective community (Kappel-L∇B, simsopt/DESC
  single-stage) optimizes cheap geometric proxies; we are the neutronic evaluator/calibrator they lack.
  Explicitly disclaim the differentiable-transport-term lift as future work (§2.2).
- **Reachability.** Vision. The evaluator is reachable; framing it as *the* neutronics layer for
  generative design is a positioning bet that depends on the field continuing to move single-stage
  (it is) and on the automation maturing.

---

## 4. Concrete experiments the chain newly enables

These were not obvious before the Kappel connection because there was no cheap, validated upstream
variable to anchor the magnet-side story.

### 4.1 L∇B vs the group's magnet-load concentration across the Zenodo-40 (+ QUASR sub-sample)

Correlate **L∇B** (read from Kappel's Appendix D for the 40; compute from VMEC `B+∇B` via virtual
casing for the QUASR sub-sample and QI) against the group's **magnet-load concentration** (free-
streaming per-coil Gini/PR — cheap, no DAGMC — for all; adjoint through-shield contributon for the
DAGMC subset). This is the direct test of the full chain at the largest reachable n. **Cheap half:**
L∇B→separation is a table lookup; free-streaming concentration is committed machinery. **Expensive
half:** the through-shield magnet load needs the DAGMC+adjoint subset. Honest expectation: the
free-streaming concentration may key more on coil descriptors than on L∇B (magnet side is not a clean
projection); report that decomposition rather than forcing a single-variable law.

### 4.2 The controlled-standoff experiment (the confound-breaker AND the mechanism-decomposer)

Hold one plasma boundary (fixed shape, fixed L∇B). Use simsopt stage-2 `CurveSurfaceDistance` to
generate coil sets at 3 standoffs (`standoff/a ∈ {0.7, 1.0, 1.4}`, warm-started tight←nominal←generous;
design detail in `COIL_DESIGN §3`). Measure magnet load and its attribution concentration at each
standoff. **What makes this the sharpest experiment:** it separates the *mediator* (standoff) from its
*root cause* (L∇B, held fixed), and it decomposes the standoff response into:
- **geometric dilution** (1/r² + collimation — present even with no shield),
- **shield attenuation** (more standoff = more room = more shield, bounded by the leak-around
  saturation the group already found),
- **leak-around** (irreducible streaming that placement targets).

If magnet load drops with standoff purely via dilution, the adjoint adds nothing over Kappel +
inverse-square (**the biggest risk, §0, made testable**). If the *concentration* and the *placed-vs-
uniform benefit* move with standoff in a way dilution alone cannot explain, the neutronic attribution
earns its keep. Repeat across ≥2 shapes to test whether residual shape matters at fixed standoff (the
`COIL_DESIGN §3` 2-way grid). This is the experiment a reviewer will demand and the one that most
cleanly answers whether the whole program is more than geometry.

### 4.3 The nfp → magnet-protectability trend

Kappel's headline reactor-relevant trend is **small nfp → large achievable separation** (the nine
most-separated configs are all nfp 1–2). If Link 2b holds, this implies a **novel, reactor-design-
relevant corollary: small-nfp configurations are intrinsically more magnet-protectable** — a
neutronics reason (beyond confinement/ripple) to prefer low field-period counts. Test it directly:
does magnet load / shieldability improve monotonically as nfp decreases across the device set, and is
that trend *mediated* by separation (i.e. does it vanish when standoff is controlled, §4.2)? A clean
positive is a genuinely new, quotable engineering result; a null (nfp matters for separation but not
for realized magnet load) is itself publishable and sharpens the story.

### 4.4 Does L∇B's bean-inside minimum coincide with the adjoint attribution hotspot?

Kappel (Fig. 8) shows L∇B is smallest on the **inside of the bean** cross section, coinciding with the
largest coil current density K — "where coil engineering is hardest." The group's adjoint gives the
*neutronic* attribution hotspot (where on the plasma boundary the coil load originates). **Do they
coincide?** If the adjoint hotspot lands on the bean-inside where L∇B is minimum, one figure ties root
cause (tight field gradient → close coils → thin shield) to neutronic outcome (concentrated magnet
load) at a *spatial* level — the tightest possible statement of the unification. If they are spatially
decorrelated, the honest conclusion is that magnet-side concentration is governed by coil-set geometry
somewhat independently of the plasma-field minimum, and the unification is a meta-claim ("geometry sets
load *structure*") rather than a pointwise correspondence. Cheap: both fields already exist per device;
this is an overlay, not a new run.

---

## 5. Risks, what could kill it, and the MVP-vs-ambitious split

### 5.1 What could kill it

1. **Step-1 null (leak-around irreducible).** Placed shield ≈ uniform even at high concentration →
   "concentration ⇒ shieldable" (Link 2b-ii) is false. Kills the closed-loop payoff and P2's
   shieldability half. *Mitigation:* it is still publishable as a null, and it is the top-priority
   cheap probe; the free-streaming coil law (P1) survives regardless.
2. **Trivial-dilution collapse (Link 2b-i).** Magnet load tracks standoff by 1/r² alone → the adjoint
   adds nothing over Kappel + inverse-square. *Mitigation:* the controlled-standoff experiment (§4.2)
   is designed to detect exactly this; if concentration and placement benefit move beyond dilution,
   the neutronics is justified; if not, honestly demote to "L∇B + geometry predicts magnet load,
   attribution is diagnostic not predictive."
3. **Small-n stalls P2.** If ParaStell→DAGMC→adjoint automation doesn't mature, the downstream link
   stays n≈2–3 and P2 is an illustrative trend at best. *Mitigation:* `COIL_DESIGN §4` shows compute
   is not the binding constraint (n≈20–50 overnight); the binding constraints are coil-design
   convergence and watertight-`.h5m` success rate — engineering, not physics. Budget them honestly.
4. **Niche entry.** ParaStell/Proxima or the ConStellaration ML community add a neutronics proxy first.
   *Mitigation:* differentiate hard on *adjoint attribution + mechanism* (L∇B root cause, contributon
   mediator), not on "a geometry→neutronics surrogate," which is filling up.
5. **L∇B under-predicts magnet load because it is a plasma-only variable.** L∇B captures the
   separation the plasma *permits*; the achieved magnet load also depends on how the coil design
   *realizes* that room (conformity, streaming gaps). *Mitigation:* this is not fatal — it is the
   honest two-stage framing (L∇B → ceiling on protectability; coil design → realized fraction),
   mirroring the plasma/coil two-stage split. State it that way.

### 5.2 MVP vs ambitious

- **MVP (P1, this cycle).** Method + reciprocity (QH+QA) + Step-1 transport kill-shot + free-streaming
  zoo coil law + SPF supplementary, with the QA/QH comparison **re-labeled along the standoff/L∇B axis
  and Kappel cited**. The L∇B connection enters as *framing that dissolves the confound*, costing a
  table lookup and a paragraph. No new law claimed. **Reachable without new physics.**
- **Ambitious (P2, multi-quarter, gated).** The L∇B→standoff→magnet-protectability predictive law on
  the Zenodo-40 + QI, with the controlled-standoff experiment and the nfp→protectability test. Buy
  **only after** Step-1 confirms the mechanism and the QI feasibility spike shows an affordable
  per-device cost. Honest landing: framework + validated upstream link + first neutronic test of the
  downstream link across a few standoff-spanning devices — an "illustrative predictive law," not a
  tight universal one below n≈20.
- **Vision (P3).** The neutronics evaluator for single-stage/generative design — the automation as
  infra contribution, the evaluator as the field's missing neutronic proxy, the calibration of cheap
  geometric proxies against adjoint ground truth. Evaluator reachable; the positioning is a bet on the
  field's single-stage trajectory.

### 5.3 Honesty ledger (carry into every venue)

- **Link 1 is Kappel's, validated; Link 2b is the group's, not yet established.** Never present the
  chain as end-to-end proven.
- **The group's QH/QA sitting on Kappel's line preserves ordering and the ~1.9× gap, but is not
  perfectly on the R²=0.94 regression.** Quote the ordering and "large gap," not an exact ratio.
- **n=2 is n=2** for the through-shield comparison; the L∇B connection upgrades the *upstream* half to
  n=40, not the downstream half.
- **Surrogate ≠ transport** until Step-1 closes; the −19% closed-loop benefit is surrogate.
- **Label every device by (L∇B, standoff), not by QS type.** The QA-vs-QH magnet gap is standoff/L∇B,
  not quasisymmetry.
- **Evaluator now; differentiable transport term is a named future lift**, not a promise.

---

## 6. Reconciliation with the execution plan

This vision does not change the `EXECUTION_PLAN_magnet_adjoint.md` critical path — it re-motivates it:

- **Step 1 (kill-shot)** is unchanged and remains top priority; it now also tests Link 2b-ii of the
  chain (concentration ⇒ shieldable) and, run at varied standoff, feeds §4.2.
- **Step 2 (free-streaming zoo law)** is unchanged; the L∇B framing tells it to **include L∇B (from
  Kappel's table / VMEC) as a descriptor alongside standoff, non-planarity, curvature** — turning the
  standoff decomposition into an L∇B-rooted one.
- **Phase E (cross-device metrics)** gains L∇B as the organizing plasma-field variable and the
  §4.4 spatial-coincidence overlay as a near-free figure.
- **QI feasibility spike (Step 4)** is the gate for P2; the L∇B trend (small nfp → large separation)
  makes the QI points (nfp 1–2, high L∇B, large separation) the *most informative* end of the axis to
  reach — QI is not just "another QS type," it is the high-separation anchor the predictive law needs.
- **New, cheap, do-early:** read L∇B for QH/QA and the Zenodo-40 out of Kappel's Appendix D; compute
  L∇B for the QUASR free-streaming sub-sample from VMEC `B+∇B` (virtual casing, simsopt — installed per
  `COIL_DESIGN`). This anchors the entire magnet-side story to a validated upstream variable for
  essentially no cost and should precede the QI investment.

*End. Reconcile the P1 scope here with `paper/PAPER_DIRECTION.md` (they agree; this adds the L∇B
framing and the Kappel citation), and fold §4.1/§4.4 into `EXECUTION_PLAN` Phase E as cheap L∇B-anchored
figures.*
