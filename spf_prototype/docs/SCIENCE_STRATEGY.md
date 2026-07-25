# Science strategy — where the magnet-adjoint program should go next

*A critical, grounded assessment of the "geometry → magnet-protectability as a predictive law"
thesis, plus a prioritized science roadmap. Written 2026-07-25 against the committed infra
(`EXECUTION_PLAN_magnet_adjoint.md`, `RESULTS_adjoint.md`, `METHODS_DRAFT.md`,
`OVERNIGHT_2026-07-22.md`) and the `shield_opt/` code, verified fact-by-fact. Skeptic's stance:
the goal here is to name what could falsify the thesis and to separate what is publishable **now**
from what is a multi-quarter bet.*

---

## 0. Verdict, risk, next steps (read this first)

**Verdict (one paragraph).** The *method* — an adjoint coil-kerma flux used as a coil→plasma
**attribution** map on a 3-D stellarator, reactivity-re-weighted for free, feeding a placement-driven
breeder-for-shield trade, with a cheap free-streaming first-wall culprit map as a zoo-wide companion —
is the right bet and is publishable **now** as a methods + honest-first-comparison paper. The
**predictive law** ("attribution concentration predicts localized-shielding benefit, and concentration
is predictable from geometry") is *over-reached at n=2* and is confounded three ways: (i) the QA/QH pair
differs mostly in **coil standoff** (1.9×), which plausibly drives *both* concentration and
shieldability, so a 2-point concentration→benefit correlation cannot be separated from
standoff→everything; (ii) the benefit metric feeding the law is a **surrogate with the wrong functional
form** — the group's own finding is that coil dose **saturates** (leak-around-limited, not exponential),
which the `exp(−kδ)` surrogate cannot represent; (iii) the through-shield magnet side is **not** a clean
geometric projection the way the free-streaming first-wall ξ=0.96 is, so there is no strong prior that the
*same* plasma-shape descriptors predict magnet-side concentration. Treat the law as a **hypothesis to
test**, not a deliverable to claim. MVP = method + a zoo-wide *free-streaming* coil-concentration→geometry
relationship + an honest, standoff-controlled n=2 through-shield comparison. Ambitious = the through-shield
predictive law, which needs ParaStell→DAGMC→adjoint automation, n≈7+ standoff-spanning devices (QI from
ConStellaration/CIEMAT/Stellaris), and a **transport-grounded** benefit — and even then likely lands as an
"illustrative trend," not a tight law, below n≈20.

**Single biggest risk.** The transport-grounded benefit (once you actually build the placed-vs-uniform
shield and run coil dose) shows that leak-around is **irreducible** — placement buys little *even when the
contributon is concentrated* — so the core mechanism "concentration ⇒ shieldable" collapses; and in
parallel the standoff confound means the "geometry → concentration" half is really "standoff → everything."
Both failure modes are cheap to probe and must be probed **before** any QI-DAGMC investment.

**Top-3 concrete next steps** (each expanded in §5):
1. **Kill-or-confirm the mechanism (transport, QH+QA).** Build the adjoint-*placed* non-uniform-shield
   DAGMC and an equal-material-budget *uniform* shield DAGMC; run coil dose in full transport; report the
   real placement benefit and whether it is larger on the higher-concentration device. This is the one
   experiment that can falsify Claim 1 outright and it retires the surrogate.
2. **Zoo-wide free-streaming coil-concentration law (cheap, real n).** On the QUASR zoo you already have,
   correlate free-streaming per-coil-flux concentration (Gini/PR) against **coil** descriptors (standoff,
   non-planarity, winding curvature) *and* plasma descriptors (elongation, iota). This attacks the standoff
   confound with a real n, and is publishable even if the through-shield law never materializes.
3. **Freeze the 4-metric definitions, compute QH/QA honestly, and do a QI feasibility spike** — pull 1–2
   ConStellaration QI boundaries through coil-design → ParaStell → DAGMC → adjoint to *measure* the
   per-device cost before committing to n=5.

---

## 1. Critical assessment of the thesis

The thesis chains two claims. Keep them separate — they can fail independently.

> **Claim 1 (engineering).** Attribution concentration (Gini / participation ratio on the contributon
> map) predicts whether *localized* shielding works.
>
> **Claim 2 (geometry law).** That concentration is itself predictable from geometry (aspect ratio, QS
> type, coil non-planarity/standoff, |B| spectrum, iota) — the magnet analog of the first-wall
> ξ=0.96 elongation→concentration finding.

### 1.1 Is Claim 1 falsifiable and testable *with this infra*? Partly — but not as currently wired.

**The good news: the mechanism is coherent, and the group's own saturation finding is what makes it
interesting.** ADJ-1 in the execution plan established that the peak coil dose **saturates** with uniform
shield thickness (~7e-3 flat from w35 to w90) — it is *leak-around-limited*, not exponential. That kills
"add thickness" and elevates "re-route shield to the leak paths" as the only lever. A **concentrated**
contributon means the leak paths are few and targetable → placement should help; a **diffuse** contributon
means load leaks in from everywhere → placement can't help. So "concentration ⇒ shieldable-by-placement" is
a real, physical, falsifiable mechanism, and the saturation result is its motivation, not a nuisance.

**The fatal wiring problem: metric 2 (the benefit) cannot currently test the mechanism.** The closed-loop
benefit (`adjoint_closed_loop.py` → `attenuation_surrogate.py`) uses `D(δ)=D₀·exp(−kδ)` with uncalibrated
literature MFPs. That functional form is *exactly the regime the group proved wrong* (ADJ-1: saturation /
leak-around). An exponential surrogate is monotone and unbounded-in-benefit; it structurally *cannot*
represent "shield more and nothing happens because neutrons leak around." So the current metric-2 number
(−19% Bosch–Hale / −55% uniform / −1% control) is not a valid measurement of shieldability, and correlating
it against concentration would be correlating concentration against an artifact of the surrogate. **To test
Claim 1 you must replace metric 2 with a transport-grounded benefit**: build the adjoint-placed shield and a
uniform shield of equal material budget as two DAGMC models, run coil dose, and compare. That is per-device
and DAGMC-limited — feasible on the n=2 you have, expensive to scale. This is the central honest gap and it
is squarely in the execution plan (ADJ-1's revised Phase C: "a DIRECT OpenMC evaluation ... rather than
surrogate calibration"). Do it.

**Confounds on Claim 1:**

- **Standoff is a lurking variable for both concentration and benefit.** QA coils sit 1.9× further
  (3.11 m vs 1.63 m) → 25–40× lower flux (`OVERNIGHT_2026-07-22.md` follow-up). A farther coil sees a more
  collimated, more attenuated flux, which likely changes *both* the plasma-boundary concentration of the
  contributon *and* the leak-around geometry that sets shieldability. So with the only two devices differing
  primarily in standoff, any concentration→benefit trend is inseparable from standoff→(concentration,
  benefit). You cannot conclude concentration is the *predictor* rather than a *co-symptom*.
- **The optimizer's degrees of freedom bound the achievable benefit, independent of concentration.** The
  benefit also depends on how much breeder you may sacrifice before the TBR floor bites and on the shield
  Fourier basis (M=3, N=2). A concentrated contributon in a place where the TBR floor forbids thinning the
  breeder yields *low* benefit despite *high* concentration — a real decoupling of the two metrics that has
  nothing to do with geometry.
- **n=2 makes any "prediction" a two-point line.** Chatterjee ξ over 2 points is degenerate. The execution
  plan (E3) already concedes this and labels it "anecdotal / illustrative." Respect that in the paper.

### 1.2 Is Claim 2 ("concentration predictable from geometry") plausible? Weaker on the magnet side.

The first-wall ξ=0.96 is a **free-streaming** result: first-wall load is a near-clean line-of-sight
projection of the plasma shape, so elongation sets concentration *directly and geometrically*. The magnet
side is through ~1 m of blanket + shield — the contributon is **smeared by scattering and dominated by
leak paths and standoff**, i.e. by exactly the physics that *breaks* the clean geometric projection. So the
mechanism that makes elongation→concentration crisp on the first wall is the mechanism the shield destroys
on the magnet side. Concretely, magnet-side concentration is more likely a function of **coil-set**
descriptors (standoff, non-planarity, winding curvature, coil-plasma conformity) than of the plasma-shape
descriptors that dominate the first wall. That is not fatal to the program — "coil geometry sets coil-load
concentration" is a perfectly good law — but it **breaks the headline unification** ("the *same* geometry
sets *both* first-wall and coil loading"). The honest version is: *different* descriptors govern the two
surfaces; the unifying statement is only the meta-claim "geometry sets load *structure*," which is softer.

### 1.3 What single result would KILL the thesis?

- **Kills Claim 1 (cleanest).** One device where the contributon is highly concentrated (high Gini) but the
  **transport** placed-vs-uniform benefit is ≈ 0 within noise — i.e. leak-around is irreducible even for a
  concentrated load. That directly falsifies "concentration ⇒ shieldable." (This is why step 1 in §0 is the
  priority: it is the experiment most able to end the program early and cheaply.)
- **Guts Claim 1 by confounding.** Across the reachable device set, concentration turns out collinear with
  standoff (concentration is just a proxy for how far the coils sit). Then "concentration predicts benefit"
  adds nothing over "standoff predicts benefit," and the mechanistic story is vacuous.
- **Kills the unification (Claim 2 as stated).** Magnet-side concentration correlates with coil descriptors
  but is ~independent of the plasma-shape descriptors (elongation, iota) that drive the first wall — so
  "geometry sets *both*, via the same knobs" is false; you have two separate laws, not one.

**Net.** The mechanism is sound and testable; the *predictive law* is not reachable at n=2 and is
confounded by standoff and by an invalid benefit metric. Reframe the law as a hypothesis with a concrete
falsification test (step 1), and deliver the framework + first comparison + a *free-streaming* coil law
(step 2) as the paper.

---

## 2. What is genuinely novel vs incremental

| Contribution | Novelty | Prior art that threatens it | Honest framing |
|---|---|---|---|
| **Adjoint coil-kerma flux as a coil→plasma *attribution* map** on a 3-D stellarator, feeding a **closed-loop breeder-for-shield** trade | **Genuinely novel** — the *use* (attribution + design lever), not the tool | CADIS/FW-CADIS/ADVANTG use adjoint flux as a **variance-reduction weight**; ADVANTG was already applied to **HELIAS shielding**; Miralles-Dolz / Thea use random-ray FW-CADIS for **weight windows** | Claim the *attribution + closed-loop* framing, **not** "adjoint on a stellarator" (done) nor "adjoint importance map" (that *is* the adjoint flux by definition). The novelty is treating ψ† as a design-attribution field and closing a loop on it. |
| **Reactivity-independence trick** (solve ψ† once, re-weight S(r) for free) | **Real but not fundamental** | It is a textbook property of the adjoint (the adjoint solution is source-independent by construction) | Frame as a *practical* contribution — makes scenario/emissivity studies and per-coil ranking free — not a discovery. Verified in `adjoint_placement.contributon` (importance loaded once, S(r) re-weighted). Genuinely useful for the paper's cheap emissivity sweep. |
| **Two-tier design**: free-streaming culprit map zoo-wide (first wall) + adjoint DAGMC-limited (magnets) | **Defensible as an architecture** | — | Honest caveat: the two tiers measure **different physics** (free-streaming FW projection vs through-shield magnet attribution), so it is *complementary methods*, **not** one method at two fidelities. Sell it as "cheap screen + expensive confirm," and be explicit they are not the same quantity at two resolutions. |
| **Predictive law** concentration→shieldability, geometry→concentration | **Would be the headline — if reachable** | Miralles-Dolz symbolic-regression neutronics surrogates; the ConStellaration ML challenge; **arXiv 2604.26763** already links **coil non-planarity to surface geometry across a QI dataset** | The differentiator is *mechanistic + adjoint-attribution-based* (concentration as the mediating variable), not a black-box surrogate. But the "geometry → neutronics" surrogate space is filling up fast; at n=2 you cannot compete on the law, only on the mechanism. |
| **First-wall concentration law** (elongation→ξ=0.96, iota→peak) on a real n=83 zoo | **Solid, already committed** | ParaStell parametric sweeps (per-device, no decomposition); Lyytinen 2024 | Keep as spine part 1. This is the strongest empirical result in hand. |

**The two sharpest prior-art threats to flag in the manuscript:**

1. **ParaStell already demonstrated 3-D blanket/shield thickness variation for magnet shielding + TBR on
   WISTELL-D.** So "vary the shield in 3-D to protect stellarator magnets" is published. Your differentiator
   is that the adjoint tells you **where** (placement, data-driven target) versus their brute-force
   parametric sweep — and that you close a loop. Cite ParaStell (Frontiers Nucl. Eng. 2024) and draw that
   line explicitly, or a reviewer will.
2. **Coil-non-planarity ↔ geometry across QI is already being published** (arXiv 2604.26763). This directly
   overlaps adjacent direction (a). Their axis is coil-geometry↔*plasma-surface*-geometry; yours must be
   coil-geometry↔*neutron-loading*. Stay on the neutronics axis to keep clear water.

---

## 3. The gating constraint, honestly costed

**What a ParaStell→DAGMC→adjoint per device actually costs — evidence from the QA bring-up.** The plan's
Phase B and the ADJ-1/2/3 adjustments are the real cost data, and they are sobering:

- QA "standard setup segfaults" on a different material-tag set (9 tags, no `Vacuum`) → device-parameterized
  `build_materials` needed (B2).
- Coil-cell identification cannot use `Cell.fill` (NotImplementedError for DAGMC, ADJ-2) → centroid-matching
  workaround, and the **Wechsung×100 frame-match is poor** (221 cm offset vs QH's ~35 cm) — an unresolved
  scale/transform saga that must be fixed before the QA adjoint source can even be trusted (B3, ADJ-2).
- The demo coil itself was mis-identified (forward-hottest cell 15 is a reciprocity *outlier*; adjoint-peak
  is cell 20, ADJ-3).

So today, **a new DAGMC+adjoint device is a multi-day, human-in-the-loop bring-up**, not a button. The
one-button **ParaStell→DAGMC→adjoint automation is itself the publishable infra contribution** — precisely
*because* QA proved it is not free. Budget it as real engineering, not glue.

**Per-device cost breakdown (realistic):**

| Stage | QUASR QA/QH (coils shipped) | QI (ConStellaration etc.) | Notes |
|---|---|---|---|
| Equilibrium | in hand (QUASR) | free-boundary VMEC from boundary | ConStellaration ships boundaries + ideal-MHD equilibria + a forward model |
| **Coils** | shipped (simsopt) | **must design** (REGCOIL/simsopt or the ConStellaration coil benchmark) | QI coil design is a real extra step; "simple-to-build QI" is an *open benchmark*, not a solved supply |
| Radial build → watertight DAGMC | ParaStell (Cubit-licensed or CAD-to-DAGMC) | same | ParaStell automates it; QA showed material-tag + coil-cell identification still bite |
| Adjoint solve (once, reactivity-independent) | random-ray adjoint, ~1 cluster job | same | free re-weight across reactivity is the one genuine cost saving |
| Forward per-coil (reciprocity) | 1 cluster job | same | needed to validate each new device |
| **Human debug** | **days** (QA evidence) | **days+** (coil design on top) | the dominant cost until automation exists |

**Realistic n.** From n=2, with the automation built, n≈5 QI + resurveyed QA/QH → **n≈7** comparable
DAGMC+adjoint devices in a paper cycle is ambitious-but-plausible; without automation, n stays ≈2–3. A
*tight* ξ law wants n≥20–30 (Miralles-Dolz screened tens); that is a multi-quarter program, not this paper.
So: n≈7 buys an **illustrative trend + a hypothesis**, not a law. Plan the paper around that honestly.

**Where reactor-relevant QI equilibria actually live (concrete, public):**

- **ConStellaration** (Proxima Fusion + Hugging Face, arXiv:2506.19583; `huggingface.co/datasets/proxima-fusion/constellaration`;
  `github.com/proximafusion/constellaration`) — ~160k QI-like boundaries, ~7,500 QI equilibria with a
  stabilizing vacuum well, nfp 1–4, **with ideal-MHD equilibria, a forward model, and scoring code**. The
  single best public QI-boundary source. Caveat: **boundaries, not coil sets** — you supply the coil-design
  step. NeurIPS 2025 dataset track.
- **CIEMAT-QI4X** (Nuclear Fusion 2025, IOP `10.1088/1741-4326/ae54ad`) — a **reactor-relevant** QI config
  with fixed- and free-boundary VMEC equilibria *including coils*, island-divertor compatible. A concrete,
  named, reactor-scale QI point design.
- **Stellaris** (Proxima Fusion, *Fusion Eng. & Design* 2025, `S0920379625000705`) — high-field QI
  prototypical power plant, free-boundary VMEC. Reactor-relevant, high-field regime.
- **Thea Energy Eos / Helios** (4 *Nuclear Fusion* papers, Jan 2025) — **planar-coil-array** QI-like device.
  Very relevant to adjacent direction (a) because the coil topology is a dipole array, *not* filaments — but
  that also means your filament adjoint-source machinery does not transfer directly; treat as a stretch
  target / discussion, not a near-term DAGMC device.
- **Near-axis QI database** (JPP, Cambridge) and **HELIAS 5-B** (EUROfusion; DCLL blanket; DAGMC; note they
  used **ADVANTG** = adjoint/CADIS VR for shielding — prior art to differentiate on the *attribution* axis).

**Tooling reality:** **ParaStell** (svalinn/parastell, Frontiers Nucl. Eng. 2024) is the right automation
backbone — it already builds stellarator DAGMC (in-vessel + magnets) from equilibrium + parametric radial
build and *demonstrated* 3-D shield/TBR variation. Lean on it; do not rebuild it. The gap you fill is
wiring its output into the adjoint-attribution + reciprocity pipeline as a repeatable per-device job.

---

## 4. The two adjacent directions

### (a) Coil-complexity → neutronics penalty — **do a free-streaming version now; high value**

**Question.** Does coil non-planarity worsen magnet loading, or does it *improve* it via standoff (a more
non-planar coil can sit farther / conform to a bay)? This is genuinely open and it is the **question that
dissolves the standoff confound** currently muddying the whole magnet story — turning the confound into a
result.

**Feasibility: high, and mostly already computable.** You have per-coil forward flux (QUASR zoo,
free-streaming), the adjoint attribution + concentration modules, and standoff is already computed
(QA 3.11 m, QH 1.63 m). Non-planarity and winding curvature are cheap geometric descriptors from the
filament polylines you already parse (`reciprocity_check._parse_coils`). So a **free-streaming** version —
correlate free-streaming per-coil-flux concentration/peak against (standoff, non-planarity, curvature,
elongation, iota) across the existing zoo — is reachable **without new DAGMC builds**. That gives a real n
and a first, honest coil-side law on a defensible proxy.

**The interesting hypothesis** is non-monotonicity: non-planarity and standoff trade off, so the net
neutronic penalty is U-shaped — exactly the nonlinear structure Chatterjee ξ is built to catch (and exactly
why ξ, not Pearson, is in `descriptor_correlation.py`). Publishable finding either way.

**Threat.** arXiv:2604.26763 (coil non-planarity ↔ surface geometry across QI) is adjacent. Differentiate:
their response is *plasma-surface geometry*; yours is *neutron loading*. Stay on the neutronics axis.

### (b) Attribution steers design (contributon → shape/coil optimization; "Paul's adjoint-shape world") — **outlook only; low near-term value, high competition**

**Feasibility: low for this paper; large lift; most competitive space.** The contributon gives
∂(coil load)/∂(birth location). Steering *design* needs ∂(coil load)/∂(boundary Fourier coefficient) or
∂/∂(coil DOF) — a **shape gradient** requiring differentiation of the equilibrium *and* the transport,
i.e. the full adjoint-shape-optimization machinery (simsopt/Landreman/Paul). Crucially, the
reactivity-independence trick does **not** give you this — it re-weights the *source*, not the *geometry*;
it provides no geometric derivative. So (b) is a whole new adjoint-shape program that overlaps heavily with
groups whose comparative advantage *is* shape optimization and whose advantage in *neutronics attribution*
is exactly where you are strong. **Verdict: keep (b) as a one-paragraph Discussion outlook** ("the
contributon is the neutronics half of a coupled shape gradient; wiring it to a differentiable equilibrium is
future work"), not an in-scope arm.

---

## 5. Prioritized implementation roadmap (science, not infra)

Ordered so the **cheap kill-shots come first** and no expensive QI infra is bought until the mechanism
survives. Each step: what to run, cost, gate, what it buys.

### Step 1 — Transport-ground the placement benefit on QH + QA (the kill-or-confirm) — **DO FIRST**
- **Build/run.** For QH (tight stats) and QA: build two DAGMC variants at *equal material budget* — (i)
  the adjoint-*placed* non-uniform shield from `adjoint_placement.placement_priority`, (ii) a *uniform*
  shield. Run full-transport peak coil dose (existing VR / weight windows) on both. Report placed-vs-uniform
  benefit **in transport**, and whether it is larger on the higher-concentration device.
- **Cost.** 2 DAGMC rebuilds + a few cluster runs each. Reuses the coil_run/percoil machinery. Days.
- **Gate.** Placed shield beats uniform beyond per-bin noise on QH; benefit correlates (right sign) with
  concentration across the 2 devices. If placed ≈ uniform even at high concentration → **Claim 1 falsified;
  stop and rewrite the paper around the null** (still publishable: "placement can't beat leak-around").
- **Buys.** Retires the invalid `exp(−kδ)` surrogate; converts metric 2 to transport; the single most
  decisive experiment. Directly executes ADJ-1's revised Phase C.

### Step 2 — Zoo-wide free-streaming coil-concentration law — **cheap, real n, attacks the confound**
- **Build/run.** On the existing QUASR zoo, compute free-streaming per-coil-flux concentration (Gini, PR via
  `concentration.py`) and peak; parse coil filaments for standoff, non-planarity, winding curvature; run
  `descriptor_correlation.chatterjee_xi` of {standoff, non-planarity, curvature, elongation, iota, nfp,
  QS type} → {coil concentration, peak}.
- **Cost.** Low — free-streaming, no DAGMC. Reuses committed modules.
- **Gate.** ξ estimates with honest CIs; explicitly test whether concentration is collinear with standoff
  (partial-dependence or standoff-stratified ξ). Report even if the top driver is standoff.
- **Buys.** A real-n coil-side "geometry → concentration" relationship on a defensible proxy; the
  standoff-vs-shape decomposition that the n=2 through-shield comparison cannot give; delivers adjacent
  direction (a) as a byproduct. **Publishable even if the through-shield law never closes.**

### Step 3 — Freeze the 4 metrics + compute QH/QA (framework) — parallel with 1–2
- **Build/run.** `shield_opt/magnet_decomposition.py` per plan E1: (1) contributon concentration
  (Gini/PR + plasma-boundary projection), (2) **transport** benefit from Step 1 (not surrogate), (3) angular
  P1/P0 shift (scalar suffices — R3), (4) reactivity robustness (free re-weight — R4). Tabulate QH vs QA.
- **Cost.** Low once Steps 1–2 land; each metric must reproduce the known QH value (gate).
- **Buys.** The framework contribution + the honest n=2 comparison, standoff-labeled.

### Step 4 — QI feasibility spike (measure the per-device cost before committing) — gates the ambitious arm
- **Build/run.** Pull **1–2 ConStellaration QI boundaries** → coil design (simsopt/REGCOIL or the
  ConStellaration coil benchmark) → ParaStell radial build → watertight DAGMC → adjoint + reciprocity.
  Instrument the wall-clock and the manual-debug points.
- **Cost.** Unknown — that is the point; QA suggests days-to-week+ with the added coil-design step.
- **Gate.** One QI device reaches a reciprocity-validated adjoint map. Record the true cost.
- **Buys.** A grounded estimate of reachable n, and the first non-QUASR (non-QS) device — the standoff/QS
  spread the law needs. **Decision point:** if a device costs > ~1 week, the tight law is out of scope; ship
  MVP.

### Step 5 (ambitious) — the through-shield predictive law
- **Build/run.** Add n≈5 QI devices (ConStellaration + CIEMAT-QI4X + Stellaris), standoff-spanning; repeat
  Steps 1+3 per device; fit concentration→transport-benefit and geometry→concentration with **standoff
  controlled**; report ξ with CIs.
- **Cost.** High — the ParaStell→DAGMC→adjoint automation (itself an infra contribution) + n≈5 bring-ups +
  transport benefit each.
- **Gate.** Standoff-controlled ξ is stable and the concentration→benefit link survives out-of-sample.
  **Honesty rule:** at n≈7 this is an *illustrative trend*, not a law; say so.
- **Buys.** The headline predictive law *if* it survives — the ambitious paper. Realistically an "n≈7
  framework + trend + hypothesis for a zoo-scale study," which is still a strong Nuclear Fusion paper.

### MVP vs ambitious
- **Minimum viable science paper (reachable this cycle):** the adjoint attribution + closed-loop **method**;
  the **transport-grounded** placement benefit on QH+QA (Step 1); the **zoo-wide free-streaming
  coil-concentration→geometry** relationship with the standoff decomposition (Step 2); the 4-metric
  framework + honest, standoff-labeled n=2 through-shield comparison (Step 3). Claims: *method* + *first
  comparison* + *free-streaming coil law* + *through-shield law as a stated hypothesis*. No overclaim.
- **Ambitious version (multi-quarter, gated on Step 4's cost):** the through-shield predictive law across
  n≈7 standoff-spanning QI+QS devices, enabled by the ParaStell→DAGMC→adjoint automation. Buy this **only
  after** Step 1 confirms the mechanism and Step 4 shows the per-device cost is affordable.

---

## Appendix — infra facts verified against the code (for provenance)

- **ψ† reactivity-independent → free re-weight.** Confirmed: `adjoint_placement.contributon` loads the
  committed importance map once and multiplies by `plasma_source_on_mesh(..., emissivity=...)`; emissivity
  enters *only* in S(r) (`emissivity.py`), never in ψ†. Re-weighting is free. ✓
- **Through-shield magnet attribution needs a DAGMC build per device.** Confirmed: Phase B + ADJ-2 (QA
  material-tag mismatch segfault, coil-cell centroid-matching workaround, `Cell.fill` NotImplemented,
  Wechsung frame-match saga). ✓
- **Free-streaming first-wall culprit map is zoo-wide cheap; ξ=0.96 found.** Confirmed:
  `PAPER_DIRECTION.md` §2 (n=83, elongation→concentration ξ=0.96). ✓
- **Reciprocity +0.76** (raw) / **+0.70** (length-normalized), n=20, coils 15/28 outliers. Confirmed:
  `reciprocity_check.py`, `RESULTS_adjoint.md` R1. ✓
- **Angular P1 built; scalar suffices for QH placement.** Confirmed: `compare_p0_p1.py`, R3 (median
  P1/P0≈1.006, φ₀ unchanged); P1/P0→~1.7× only deep in a scattering slab (M6). ✓
- **Benefit is surrogate; coil dose saturates (leak-around, not exponential); benefit is placement- not
  thickness-driven.** Confirmed: `attenuation_surrogate` `exp(−kδ)` in `adjoint_closed_loop.py`; ADJ-1
  (dose ~7e-3 flat w35→w90). This is why metric 2 cannot currently test the thesis (§1.1). ✓

## Sources (QI equilibria, tooling, prior art)

- ConStellaration — arXiv:2506.19583; `huggingface.co/datasets/proxima-fusion/constellaration`;
  `github.com/proximafusion/constellaration`; NeurIPS 2025 dataset track.
- CIEMAT-QI4X — *Nuclear Fusion* 2025, IOP `10.1088/1741-4326/ae54ad`.
- Stellaris (Proxima) — *Fusion Eng. & Design* 2025, `S0920379625000705`.
- Thea Energy Eos/Helios — 4 *Nuclear Fusion* papers, Jan 2025 (planar-coil QI-like).
- QUASR — `quasr.flatironinstitute.org` (QA/QH, ~370k configs, simsopt coils); Landreman & Paul.
- Coil non-planarity ↔ QI surface geometry — arXiv:2604.26763 (threatens adjacent direction (a)).
- ParaStell — Frontiers Nucl. Eng. 2024; `github.com/svalinn/parastell` (3-D shield/TBR sweep on WISTELL-D).
- HELIAS 5-B neutronics — EUROfusion (DCLL blanket, DAGMC, ADVANTG adjoint VR — differentiate on attribution).
- Miralles-Dolz — symbolic-regression neutronics surrogates + Chatterjee ξ (the descriptor-screening method
  reproduced in `descriptor_correlation.py`).
