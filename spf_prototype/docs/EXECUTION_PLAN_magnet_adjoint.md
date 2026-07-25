# Execution plan — magnet-side adjoint program (through-shield coil attribution → cross-device ξ)

Owner: Claude Code, executing autonomously with hard correctness gates. Style (per CLAUDE.md):
hard-gate on correctness, strict phase ordering (don't start N+1 until N's gate passes), honest
non-cherry-picked diagnostics, no warnings-where-an-assert-belongs. Everything on branch
`spf-prototype`, never pushed. Ginsburg runs via `ssh -o BatchMode=yes ginsburg` (user does 2FA).

**Goal.** Turn the (now built + validated) single-device adjoint magnet-shielding machinery into a
**cross-device finding** — the magnet-side analog of the first-wall ξ decomposition — completing the
geometry-first thesis (geometry → *both* first-wall and coil loading). Plus close the three
honesty gaps the methods draft surfaced.

**Current state (done, committed).** QH adjoint through-shield importance (runs on real DAGMC),
filamentary coil-response adjoint source, contributon C=S·ψ†, placement adapter, closed-loop
optimizer (19% *surrogate* benefit), reciprocity validation (Spearman +0.76, n=20), M7 angular
(P1 built into the DAGMC build; scalar suffices for placement), realistic-reactivity default,
near-void fix, self-reporting assumptions. Methods draft written.

---

## Phase A — Provenance & un-staling (low-risk, do first; unblocks a consistent paper)

**A1. Un-stale `paper/PAPER_DIRECTION.md`.** It still calls the adjoint map "NOT YET BUILT /
aspirational." Update §4 + the Limitations/Open-questions to reflect: built, runs on real DAGMC,
reciprocity-validated (+0.76), M7 angular done (scalar suffices), realistic-reactivity default.
Keep the honest caveats that remain (surrogate benefit pending calibration, single device→QA in
progress, per-coil outliers 15/28, P1-only angular). *Deliverable:* updated doc. *Gate:* no claim
in the doc contradicts a committed result or the honest caveats.

**A2. Freeze `paper/RESULTS_adjoint.md`.** A citable provenance doc: each headline number with its
run parameters, commit hash, geometry (QH LP reactor, w20/w35f), emissivity, response, rel.err, and
the honest caveat. Numbers: reciprocity +0.76 (n=20, w35f, uniform-matched); closed-loop 19%
Bosch-Hale / 55% uniform / −1% baseline (flat response, surrogate objective); M7 P1/P0 median ~1,
placement φ0 unchanged; peak coil reactivity-robust (hot-5 unchanged). *Gate:* every number traces
to a committed artifact or a logged run; mark surrogate-vs-transport explicitly.

*Phase A is documentation — safe to run in parallel with B/C.*

---

## Phase B — QA new-device bring-up (fix + execute) — PREREQUISITE for the cross-device finding

The QA DAGMC (`~/pstl_test/qa/dagmc_qa_base.h5m`) has a **different material-tag set** than QH
(9 tags, no `Vacuum`) → the standard setup segfaults. QA = 16 coils, nfp=2, cells to be identified.

**B1. Reproduce the working QA material/geometry setup.** Base it on the *known-good*
`~/pstl_test/qa/qa_coil_discovery.py` (it produced `qa_cell_hits.npz`). Identify the exact
material list + DAGMC load that does NOT segfault, and the **16 coil cells** (the `magnets`-filled
cells — not the "last-16 id" heuristic, which included an R≈0 central cell). *Gate:* `find_cell`
binning yields exactly 16 magnets-cells with sensible centroids (R in the coil shell, not the axis).

**B2. QA-specific `build_materials`.** `adjoint_importance.py` imports the QH 10-material
`build_materials`; make it device-parameterized (QA = the 9-tag set). *Gate:* the QA DAGMC loads +
MGXS generates without the material NotImplementedError/segfault.

**B3. Frame-match the QA coils.** Parse `Wechsung_QA24.coils` (16 filaments), find the scale
(likely ×100 m→cm) that maps filament centroids bijectively onto the B1 coil-cell centroids
(~tens of cm, like QH). *Gate:* 16↔16 bijection, mean offset < ~1 coil cross-section.

**B4. QA adjoint importance sweep.** Run `adjoint_importance.py` for the QA hot coils (identify the
QA peak coil from `percoil_unpol_qa_base.npz` first), flat + kerma responses, on `dagmc_qa_base`,
with the QA fluxmap (`~/qa_reactor_fluxmap.npz`, nfp=2) + QA centroids + QA filaments. *Gate:*
structure gate passes (localized, rel.err ≲0.25), placement points at the coil.

**B5. QA reciprocity + placement + closed-loop.** Reciprocity vs `percoil_unpol_qa_base.npz`
(uniform-matched); placement priority; closed-loop benefit (Bosch-Hale). *Gate:* QA reciprocity
Spearman > 0.6 (method reproduces forward on a 2nd device); benefit beats the no-attribution
baseline. *Deliverable:* QA maps committed, QA numbers logged in RESULTS_adjoint.md.

*Depends on nothing but the QA data (present). The #4 filament fix (Phase D) improves B4/B5 accuracy.*

---

## Phase C — Surrogate MFP calibration (transport-ground the benefit)

Make the closed-loop benefit an OpenMC-grounded number, not surrogate-only.

**C1. OpenMC coil-dose vs uniform shield thickness.** Run full-transport coil fast-flux/dose at
~3–4 uniform breeder→shield trades δ (e.g. 0, 10, 20, 30 cm) on QH (and QA), for the peak coil,
with the existing VR (weight windows) so the coil dose converges. *Gate:* per-δ coil dose rel.err
small enough to fit (≲10%); monotone decreasing in δ.

**C2. Fit λ_shield, λ_breeder.** Fit `D(δ)=D0·exp(−k δ)`, `k=1/λ_s−1/λ_b`, to the C1 curve (and a
breeder-only scan if needed to separate the two λ). Replace the literature-scale defaults
(8/17 cm) in `attenuation_surrogate.py`. *Gate:* the calibrated surrogate reproduces the C1 OpenMC
dose-vs-δ within a stated tolerance (report the residuals).

**C3. Re-run the closed loop with calibrated MFPs.** The transport-grounded benefit (replaces the
surrogate-only 19%), reported with the calibration uncertainty. *Deliverable:* calibrated benefit
in RESULTS_adjoint.md, figure regenerated. *Gate:* honestly report how much the calibration moved
the number (could be up or down).

*C runs on Ginsburg in parallel with B. C1 reuses the coil_run/percoil machinery.*

---

## Phase D — Filament source-matching fix (the #4 per-coil outliers)

Coils 15 & 28 are per-coil outliers (adjoint vs forward), traced to the filamentary adjoint source
for highly non-planar windings (coil-15 sweeps R 470→1775; 40 cm centroid offset).

**D1. Diagnose.** For the outlier coils, check filament↔cell alignment (does the matched filament
actually thread the cell?), source-point density along arc length, and whether more/segment-length-
weighted points change ∫S·ψ†. **D2. Fix.** Segment-length-weighted sampling and/or a per-coil
alignment check; re-run the outlier coils. *Gate:* the reciprocity Spearman improves (target >0.8)
OR a clear, documented reason the residual is physical (e.g. a real streaming path the point-kernel
mis-weights). Honest either way — do not tune to a target.

*Do D before finalizing B5/Phase E rankings so QA uses the improved source. Cheap (analysis + a
few re-runs).*

---

## Phase E — Cross-device magnet-side ξ decomposition (the science)

For each device with a DAGMC build (QH, QA now; framework generalizes), compute **four comparable
metrics** and correlate with geometry — the magnet analog of the first-wall ξ story.

**E1. Metrics machinery (`shield_opt/magnet_decomposition.py`).** Per device, per (peak) coil:
1. **Attribution concentration** — Gini + participation ratio on the contributon C(r) (and its
   plasma-boundary projection). Higher = more localized = more shieldable.
2. **Adjoint-informed shield benefit** — the *calibrated* (Phase C) closed-loop peak-dose reduction
   (%), realistic reactivity.
3. **Angular sensitivity** — P1-vs-P0 placement shift (deg) + median P1/P0 (M7 machinery).
4. **Reactivity robustness** — stability of the per-coil load ranking uniform→Bosch-Hale
   (rank-corr; and does the peak coil change).
Reuse committed code (`concentration.py` Gini/PR, `adjoint_placement`, `compare_p0_p1`,
`reciprocity_check --emissivity`). *Gate:* each metric reproduces the already-known QH value.

**E2. Run the metrics for QH + QA.** Tabulate the 4 metrics × 2 devices. *Deliverable:* a per-device
metrics table.

**E3. Correlate with geometry.** Descriptors: QS type, aspect ratio, coil-plasma standoff, coil
non-planarity, plasma elongation/shaping. Chatterjee ξ(descriptor → metric) using the committed
`descriptor_correlation.py`. **HONEST FRAMING:** n=2 devices → correlation is *anecdotal*; frame
E3 as (a) establishing the framework + metric definitions, (b) the first clean QA-vs-QH magnet
comparison, and (c) a hypothesis for the zoo-scale study (which needs more DAGMC builds). Do NOT
over-claim a correlation from 2 points. *Gate:* the QA-vs-QH differences are stated with their
uncertainty; any ξ over 2 points is labeled illustrative.

**E4. Report.** `spf_prototype/RESULTS_magnet_decomposition.md` + figures (per-device metric bars,
QA-vs-QH attribution maps, the framework schematic). Honest, non-cherry-picked. *Deliverable:*
committed results + figures.

*Depends on B (QA), C (calibrated benefit), D (clean per-coil). This is the synthesis.*

---

## Phase F — Paper integration

Fold E into `METHODS_DRAFT.md` (methods) + a new RESULTS section; update `PAPER_DIRECTION.md` with
the completed cross-device story (first-wall ξ zoo-wide + magnet ξ QA-vs-QH deeply). *Gate:* the
paper's claims match RESULTS_adjoint.md + RESULTS_magnet_decomposition.md exactly.

---

## Dependency graph & parallelization

```
A (docs) ──────────────────────────────────┐ (parallel anytime)
B (QA bring-up) ──┐                          │
C (calibration) ──┤ (both Ginsburg, parallel)│
D (filament fix) ─┘ (before B5/E)            │
                    └──> E (cross-device ξ) ──┴──> F (paper)
```

**Immediate parallel launches:** A (docs, me) ‖ B1–B3 (QA bring-up, me+Ginsburg) ‖ C1 (OpenMC
dose-vs-δ, Ginsburg background). D interleaves with B. E is the gated synthesis.

## ADJUSTMENTS (from kickoff, 2026-07-25)

**ADJ-1 (Phase C is now the critical scientific question, not a formality).** The existing per-coil
data at shield thicknesses w20/w35/w60/w90 shows the peak coil dose SATURATES (~7e-3 flat from w35
to w90; w20 is a known CAD-defect artifact) — it is **leak-around-limited, not exponential**. So the
closed-loop surrogate `D=D0·exp(−kδ)` is the WRONG functional form in the operating regime, and a
naive λ fit is invalid. Revised C: (a) confirm saturation cleanly; (b) recognize the benefit is
PLACEMENT-driven (re-route shield to leak paths), not total-thickness-driven — consistent with the
prior "thickness saturates → placement is the lever" finding; (c) the honest transport-grounded
benefit likely requires a DIRECT OpenMC evaluation of the adjoint-placed vs uniform shield (build the
non-uniform-shield DAGMC, run coil dose) rather than surrogate calibration. This means the 19% may be
optimistic for a different reason than reactivity, and C is a real result, not a rubber-stamp.

**ADJ-2 (Phase B — DAGMC c.fill is NotImplemented).** Cannot filter coil cells by material via
openmc.lib (`Cell.fill` raises NotImplementedError for DAGMC). Workaround: match discovered solid
cells to the coil FILAMENTS by centroid. QA coil cells = 11–26 (bijective). BUT the Wechsung×100
frame-match is POOR (mean 221 cm offset vs QH ~35 cm) → the QA coil scale/transform to the DAGMC
frame is NOT ×100 clean (QA reactor-scaling saga). B3 must resolve the correct Wechsung→DAGMC
transform before the QA adjoint source can be trusted; else per-coil assignment may be wrong.

**ADJ-3 (reconcile the demo coil).** The committed closed-loop 19% is on coil-15 = the FORWARD-hottest
but a reciprocity OUTLIER; the ADJOINT-hottest (true peak) is cell 20. Re-run the closed-loop on the
adjoint-peak coil (20) for the headline, and keep coil-15 only as a demo. Folds into Phase D (the
15/28 outlier fix) + E.

## Global honesty gates (apply throughout)
- Every reported number stamps: device, geometry (shield build), emissivity, response, objective
  (surrogate vs transport-calibrated), rel.err, n. (The self-reporting convention already added.)
- No tuning-to-target: if a gate fails, report and investigate, don't massage.
- n=2 is n=2: the cross-device "finding" is a *framework + first comparison*, not a zoo correlation.
- Surrogate ≠ transport until Phase C closes; label the 19% accordingly until then.
