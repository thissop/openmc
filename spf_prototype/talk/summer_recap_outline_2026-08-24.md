# Summer Research Recap — Slide-by-Slide Draft (v1)

**Talk:** Thaddaeus Kiker · advisor Ethan Peterson (MIT PSFC) · Wednesday group recap
**Working thesis of the talk:** *"Where fusion neutrons originate, where shielding is effective, and which magnet-load reductions are physically achievable are three distinct problems — and separating them is the key to optimizing stellarator magnet shielding."*

---

## Draft assumptions (please correct — see questions at the bottom)
- **Audience:** Peterson's group + fusion-neutronics-literate peers. Technical, but not everyone is an OpenMC/neutronics specialist → keep transport jargon lightly explained.
- **Time:** ~15–20 min talk + questions → **~16 slides** (target ~1 min/slide, a couple of lingering result slides).
- **Format:** content-first here; build as Beamer (your dove/serif style, smplotlib figures, plain titles, minimal bold) after the story is locked.
- **SPF role:** the *origin* and a supplementary knob — 1–2 slides, honestly framed as "how I got in," not the headline.
- **Centerpiece:** the **attribution ≠ worth ≠ controllability** framing, anchored by the concrete magnet-shielding results.

---

## Slide 1 — Title
- **Title:** *Optimizing magnet shielding in stellarator fusion reactors: attribution, worth, and controllability*
- Name, advisor, affiliation, date. Subtitle: "A summer of stellarator neutronics."
- **Visual:** one clean 3D render of the QH device — plasma boundary + a few real coils (you have the QH wireframe figure).
- **Notes:** One sentence of what the audience will leave with: "how much can non-uniform 3D shield design actually protect the magnets — and what limits it."

## Slide 2 — Why this matters: magnets are the thing you can't replace
- Fusion power plants need superconducting magnets; the blanket/first wall can be maintained, the magnets largely cannot.
- 14 MeV DT neutrons drive coil damage: nuclear heating (cryogenic load / quench), fast-neutron fluence (HTS degradation), insulation dose, dpa.
- The magnet **peak** load — a hot spot, not the average — sets machine lifetime.
- Shielding steals space and breeder from the blanket → a real **TBR-vs-protection tradeoff**.
- **Visual:** simple radial-build cartoon: plasma → first wall → breeder → shield → gap → vessel → **coil**, with the neutron arrow and a "hot spot" star on the coil.
- **Notes:** Frame the whole talk as a design-economics question: protect the magnets at the least cost in breeding and machine size.

## Slide 3 — Why stellarators make this hard
- No axisymmetry: the load varies strongly in **both** poloidal and toroidal angle.
- Coils are **non-planar and twisted** — plasma-to-coil clearance varies dramatically around the device (our QH coil windings sweep R ≈ 470→1775 cm).
- A uniform shield is either too thin where clearance is tight or wasteful where it's roomy → non-uniform design is *natural*, but the design space is huge.
- Device: reactor-scale **QH (Landreman–Paul), nfp = 4**, R₀ ≈ 1367 cm, a ≈ 217 cm, **8 coils per field period**.
- **Visual:** the QH coil set colored by plasma-coil clearance (or the twisted coil-15 winding highlighted).
- **Notes:** This is the "why it's interesting and underexplored" slide — tokamaks are ~axisymmetric; stellarators force the 3D question.

## Slide 4 — The summer arc (roadmap slide)
- **Started:** spin-polarized fusion (SPF) — building a polarization-aware neutron source in OpenMC.
- **Pivoted (July, w/ Ethan):** to broad **geometry-first stellarator neutronics**; magnet shielding / peaking as the key metric. SPF becomes a knob.
- **Built:** an end-to-end source → geometry → transport → attribution → optimization pipeline.
- **Now:** a structured study of how much 3D shield design freedom is worth.
- **Visual:** a horizontal arrow timeline: SPF source → pivot → pipeline → shielding-optimization science.
- **Notes:** Set expectations honestly: "a maturing framework and a strong pipeline, with the first real results — not a finished optimization."

## Slide 5 — Origin: a spin-polarized DT source sampler (SPF)
- Built a self-contained, polarization-aware DT-neutron **source** in OpenMC (compiled `StellaratorSource`).
- Emission depends on polarization mode fractions (a, b, c): aligned emission ∝ sin²θ; others ∝ ¼ + ¾cos²θ — anisotropic birth directions relative to B̂.
- **Validated** against (i) the unpolarized isotropic limit and (ii) Schwartz (2025) analytic axisymmetric results — corner integrated cross sections (2/3, 1, 2/3, 1/3 σ₀) recovered; A-mode ±43%/∓22% inboard/outboard wall-load pattern reproduced.
- Hard-gated on correctness (χ²/KS on sampled directions, analytic NWL comparison).
- **Visual:** the cosθ histogram vs analytic overlay (per-mode), or the NWL inboard/outboard pattern.
- **Notes:** Keep this crisp — it shows rigor and where you started; don't dwell.

## Slide 6 — The pivot: geometry-first stellarator neutronics
- SPF is one lever; the richer, less-explored question is **how 3D geometry shapes the coil load** — and how to design against it.
- Key metric shifts from wall loading to **magnet peaking / coil protection**.
- This reframes the project around a design-and-methods question, with SPF as a knob layered on top.
- **Visual:** a simple "SPF knob on top of geometry" schematic, or just a text pivot slide.
- **Notes:** Credit Ethan's framing; this is the intellectual turn of the summer.

## Slide 7 — The pipeline I built (capability overview)
- **Source:** VMEC-fluxmap `StellaratorSource` (real plasma emissivity), optional polarization.
- **Geometry:** ParaStell build — plasma → FW → multiplier → breeder → shield → gap → vessel → thermal shield → coils, with a **toroidal × poloidal thickness field** per layer.
- **Transport:** one-field-period + **periodic BCs** (≈4× cheaper, exact for nfp symmetry); direct analog forward evaluation.
- **Analysis:** adjoint magnet-importance + contributon attribution; FW-CADIS variance reduction (built/debugged); SLURM automation + CMA-ES.
- **Visual:** labeled radial-build render + a small "toolbox" icon row (source / geometry / transport / adjoint / optimizer).
- **Notes:** This is the "I built real infrastructure" slide — the platform enables everything after.

## Slide 8 — A useful early result: the coils aren't as buried as feared
- Naively, coils sit behind ~40 cm of shield → assumed deep, variance-reduction-limited.
- In practice, a plain **analog** forward run resolves **all 8 coils to 2–4 % in ~6 min** — no weight windows needed for this metric.
- Lesson: don't reach for fancy variance reduction before checking the analog cost — it reshaped the whole workflow (simpler, faster iteration).
- (FW-CADIS is built and kept in reserve for genuinely deep quantities, e.g. peak heating inside the winding pack.)
- **Visual:** per-coil flux bar chart with 2–4% error bars; time-to-solution callout.
- **Notes:** A nice "honest engineering" beat — a methods detour that turned into a clarifying result.

## Slide 9 — The optimization question
- At a **fixed machine envelope**, can we redistribute breeder ↔ shield to cut the **peak** coil load — for free (equal material)?
- More generally: how much does **non-uniform 3D** shield design improve the engineering tradeoff (TBR vs peak coil heating vs mass/build) vs a uniform design?
- The honest complication we'll see: the peak is a **moving target**.
- **Visual:** a cartoon Pareto axes sketch (TBR vs peak coil load), with "uniform" and "?" points.
- **Notes:** Pose it as the central engineering question the rest of the talk chips at.

## Slide 10 — Adjoint-guided placement: a strong lever, with a catch
- Localize an adjoint source on a coil → adjoint importance map ψ†(r) → **contributon** C = S·ψ† (which plasma regions drive that coil).
- Adjoint-**placed** shield cut the *target* coil **4.2×**.
- But the global **peak migrated** to a different coil — protecting one coil exposed another.
- Multi-coil placement (spread and greedy) both underperformed → the load is **collective / inboard-dominated**, not one local plasma→coil pathway.
- **Visual:** before/after per-coil bar chart showing the target coil crushed and the peak jumping to a neighbor.
- **Notes:** This is the pivotal physics observation — set up the "attribution isn't the whole story" message.

## Slide 11 — The reframe: three different maps
- **Attribution** A: which plasma regions cause a coil's load. (`S·ψ†`)
- **Intervention worth** B: where a material/thickness change actually *reduces* the load. (`−δR/δu`)
- **Feasible controllability:** which load reductions are *achievable* under mass, clearance, TBR, smoothness, and non-negative-thickness constraints.
- **These are not the same.** High attribution ≠ high worth ≠ a feasible fix. Conflating them explains peak migration and "flat" local worth.
- **Visual:** three side-by-side stylized θ–φ maps (A, B, feasible region) that visibly *disagree*.
- **Notes:** This is the intellectual core / likely the paper's spine. Slow down here.

## Slide 12 — First evidence: a resolved, placement-dependent response operator
- Built the full **12 × 20** finite-difference response matrix (patch → all-coil vector) from existing runs — the prior "flat worth" only ever looked at *one* coil.
- Result: effects are **resolved** (not noise); control is strongly **placement-dependent**.
- Concentrated inboard +20 cm patches cut the peak substantially and drop the whole inboard cluster together; the *wrong* sector **backfires (+41%)**.
- **Preliminary** signal: *concentrated* placement beats *uniform* addition at ~equal added mass (to be confirmed against a clean baseline).
- **Visual:** the fractional-response heatmap (patches × coils, blue = reduced) + peak-migration bar chart.
- **Notes:** Frame carefully — "coarse, first rung; establishes structure exists, not optimality." Don't overclaim the 5×.

## Slide 13 — What this does *not* yet show (honesty slide)
- Not optimality; not an equal-cost engineering benefit yet; not the *continuous* (fine-resolution) shielding field.
- Attribution-vs-worth is still **unresolved** — needs the importance mapped onto the true device surface before any "A ≠ B" claim.
- The load being collective is *why* naive local tricks migrate the peak — a feature to design around, not a dead end.
- **Visual:** a small "coarse basis vs continuous field" schematic (4×4 patches vs a smooth twisted mode).
- **Notes:** This slide buys credibility and pre-empts the obvious critiques. Reviewers/Ethan will respect it.

## Slide 14 — Where it's headed: a validated co-design ladder
- A fixed rung ladder: coarse operator → validate the large effect → define real metrics (heating **density**, TBR, mass, clearance) → **hierarchical geometry-aware basis** (twisted, field-period-compatible modes) → clearance-constrained optimization → **Pareto fronts**.
- Compare response-aware design against uniform, clearance-only, wall-loading, and attribution baselines **at equal engineering cost**.
- Contrasts: **SPF** (does polarization change the optimum?) and a **second configuration** (generality).
- **Visual:** the target Pareto figure (TBR vs peak coil heating), with the baseline ladder as points/curves — clearly marked "projected."
- **Notes:** Show ambition + rigor; make clear the plumbing to do this is already built.

## Slide 15 — SPF as a knob (calibrated)
- Polarization reshapes the *source* anisotropy → shifts where neutrons go → potentially shifts attribution and the limiting coil.
- The interesting question isn't "SPF lowers flux" but "**does SPF change the engineering optimum**" (shield mass / build / TBR recovered at equal protection)?
- Honest status: a planned contrast, not yet a result; some prior art exists (polarized sources in OpenMC), so we'll frame the *coupling to design* as the novel angle.
- **Visual:** the SPF emission-anisotropy sketch (sin²θ vs ¼+¾cos²θ) tied back to the coil map.
- **Notes:** Keep SPF secondary and honest; it's the origin story and a future lever, not the headline.

## Slide 16 — Summary & takeaways
- Built an end-to-end stellarator neutronics + shield-optimization pipeline (source → geometry → transport → attribution → optimizer).
- Key idea: **attribution ≠ worth ≠ controllability** — the frame that explains why naive shielding tricks fail and what to optimize instead.
- First evidence of resolved, placement-dependent magnet-load control; concentrated inboard shielding looks promising.
- Next: turn the coarse operator into a validated, clearance-constrained 3D co-design with real Pareto fronts + SPF/second-config contrasts.
- **Acknowledgments:** Ethan Peterson; collaborators (E. Paul, Sophia, Jack Fletcher); tools (OpenMC, DAGMC, ParaStell).
- **Notes:** End on the thesis line; invite the room to poke at the attribution-vs-worth framing.

### Optional backup slides
- B1: validation details (SPF χ²/KS tables; NWL residual normality).
- B2: FW-CADIS / variance-reduction methods note (why analog sufficed).
- B3: the periodic-BC one-period model + validation caveat (vs full torus).
- B4: patch-worth n=12 details / signed-worth-as-backfire-guardrail.

---

## Clarifying questions for Thaddaeus
1. **Room & level:** just Ethan's group, or a broader PSFC/Columbia audience? How neutronics-literate — can I assume "adjoint / TBR / dpa" land without definition?
2. **Exact time budget:** 12, 15, or 20 minutes? (Sets whether we cut to ~10 slides or keep ~16 + backups.)
3. **Format:** Beamer in your dove/serif style, or a shareable HTML artifact? (I'll build whichever once the story's locked.)
4. **SPF prominence:** keep it as the 1–2 slide origin/knob (current draft), or elevate it (it was where you started and may be the "distinctive" angle)?
5. **Centerpiece:** which single thing do you want people to remember — (a) the *pipeline/capability*, (b) the *4.2× kill-shot + peak migration*, (c) the *placement-dependent control* result, or (d) the *attribution-vs-worth-vs-controllability* framing? (Draft currently leads with (d), anchored by (c).)
6. **Figures on hand:** which of these do you already have plottable — QH coil render, per-coil flux bars, the 12×20 heatmap, SPF cosθ histogram, radial-build cartoon — vs. which should I generate (smplotlib)?
