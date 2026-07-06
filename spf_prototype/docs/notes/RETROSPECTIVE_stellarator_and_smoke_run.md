# Retrospective — stellarator validation and the precise_QA smoke-run deck

*A condensed, in-retrospect narrative of two threads that extended the axisymmetric
prototype toward real stellarators: (1) the free-streaming OpenMC↔analytic cross-checks
on genuinely 3-D quasisymmetric sources (`RESULTS_toy_qalow.md`, `RESULTS_quasr.md`,
now in `notes/results/`), and (2) the precise_QA conformal-blanket smoke-run talk
package (`README_spf_deck.md`, `talk_script.md`, `qa.md`, `slides/talk.tex` — left in
place at the repo root because they are coupled to `figures.py`/`figs/`/`spf_export/`).
Numbers preserved; blockquotes verbatim and attributed.*

---

## The two-paper split (the framing that organises both threads)

**At the time** a DESC survey showed a hard fact: no *standard published* QA is
natively low-ε_eff — precise_QA/reactor_QA sit at ε_eff ≈ 2.3, ESTELL 1.0, NCSX 1.3,
ARIES-CS 1.6. Real QAs are strongly shaped, which is the **Monte-Carlo regime**. This
forced a clean split that governs both the analytic and MC papers:

- **Low-ε_eff overlap** — where the analytic series/quadrature *and* the MC agree, so
  they cross-validate each other. This is the MC paper's credibility figure and the
  analytic paper's home turf.
- **High-ε_eff, MC-only** — the real strongly-shaped equilibria (e.g. precise_QA,
  ε_eff ≈ 2.3), where the analytic series diverges and even fixed-node quadrature
  breaks down (≈34% exact-vs-hybrid disagreement, unphysical negatives from the
  near-singular integrand). This is the regime the MC method owns.

> "the high-ε_eff real equilibrium precise_QA (ε_eff ≈ 2.3) is MC-only" —
> RESULTS_toy_qalow.md

---

## Thread 1 — free-streaming OpenMC ↔ analytic on 3-D QA sources

### Toy QA-low (ε_eff ≈ 0.08) — the clean series-regime check

**At the time** we wanted the MC↔analytic validation on a genuinely field-period-
rippled (Nfp=2) stellarator *source*, still on the shared square-torus CSG wall (no
DAGMC), in the regime where series/exact/hybrid all agree. **We achieved** worst-wall
directionality agreement of **1.8%** across all four walls, with poloidal shapes
tracked (inboard A/iso 1.207 MC vs 1.229 analytic; outboard 0.823 vs 0.829; floor
1.065 vs 1.069; ceiling 1.062 vs 1.066).

This run also produced a reusable engineering artifact: the analytic engine was
**rewritten as a fast pure-numpy free-streaming quadrature** evaluating the loop
geometry and field directly at the Gauss-Legendre nodes, with correct visibility (true
shaped-loop sub-arcs ∩ inner-cylinder occlusion), gated against
`anarrima.free_streaming_quadrature` to **3e-14**. It replaced an earlier crude
`front_arc` whose errors cancelled in anarrima's own series-vs-quadrature tests (same
arc both sides) but not against ray-tracing MC.

> "the exact reconciling factors are A ×3/2, B ×2 (= 4π / ∫K dΩ)" — RESULTS_toy_qalow.md

The **normalization convention** was pinned: anarrima's kernels are unnormalized
(∫sin²θ dΩ = 8π/3, ∫(¼+¾cos²θ)dΩ = 2π) while the OpenMC source samples unit-total
emission PDFs; `free_streaming_quadrature(..., normalize=True)` was added (default-off)
so the two sides are directly comparable — a small addition proposed for the anarrima
PR/JOSS.

### Real published QUASR QA 59509 (ε_eff ≈ 0.43) — the literature config

**At the time** we wanted the same check on a *real, published* equilibrium, not a toy,
so the MC paper rests on a literature config. Finding a gentle published QA required
mining the QUASR database **without the ~13 GB bulk download**:

- `quasr_catalogue.py` pulls the single 37 MB `database.json.gz` (**371,701 devices**),
  filters QA (helicity==0, 71,794 in pool), shortlists by `max_elongation`;
- `quasr_eps_eff.py`/`quasr_geom.py` download only per-device VMEC namelists (~18 kB)
  and compute geometric ε_eff in pure numpy. **corr(max_elongation, e_geom) = 0.94.**
- Gentlest published QA: **device 59509** (nfp=3, aspect 6.7, e_geom ≈ 0.29, ε_eff ≈
  0.43 at the wall). (`simsopt`'s compiled core SIGILLs on this CPU, so geometry comes
  from the namelist and the field from a pure-numpy flux-surface-tangent model.)

**We achieved** ≤1.8% agreement on outboard, floor, and ceiling (poloidal shape
tracked), with a **bounded 4–6% inboard residual** that was investigated exhaustively
and traced — *not* to the sampler (verified exact), *not* to the analytic (≡ anarrima
to 1e-12), *not* to visibility, field, or near-field — but to the **point-patch
approximation** of the filament flux integral on a strongly curved wall:

> "The point-patch integral is exact for an axisymmetric source (|Δ|=0.002 = statistical
> floor), and its error grows ≈ s² with shaping, ~3× larger on the convex inner
> (inboard) cylinder than the concave outer (outboard)." — RESULTS_quasr.md

A targeted `bare_cylinder_demo.py` reproduced the inboard-specific, shaping²-scaling
residual **by construction** (s=0.0 → |Δ| 0.002; s=1.0 → inner |Δ| 0.098 vs outer
0.029), confirming it is a **limit of the filamentary analytic**, with the
geometry-exact, sampler-exact MC as the reference there — *"exactly the kind of regime
the Monte-Carlo method is meant to own."* The gentle toy QA-low (ε_eff 0.08) shows the
residual on **no** wall; only the strongly-shaped 59509 (ε_eff 0.43) shows it, and only
inboard — consistent with the shaping-scaling explanation.

---

## Thread 2 — the precise_QA conformal-blanket smoke-run deck

**At the time** (a separate, earlier session, packaged at the repo root) we ran a
10,000-history OpenMC transport of spin-polarized DT neutrons through a **conformal
breeding blanket** on the `precise_QA` stellarator (nfp=2) — a *smoke run*, not a
converged study — to ask whether the neutron-steering signal *survives a real
stellarator wall*. The deliverables (`README_spf_deck.md`, `talk_script.md`, `qa.md`,
`slides/talk.tex`, and `figures.py` → `figs/`) are built only from a plain-data bundle
(`spf_export/`), numpy/pandas/json only — no OpenMC/DAGMC/cluster needed to reproduce.

**We achieved** an *indicative* directional-efficiency headline the deck itself flags
as unconverged:

> "η ≡ Δ_scatter / Δ_free = 8.5% / 13.5% ≈ 0.63 ± 0.22" — slides/talk.tex

i.e. roughly two-thirds of the ideal free-streaming steering survives the conformal
wall + scattering. The **robust** result was the *trend*, not the number:

> "The robust result is the trend, not the precise number — a coherent directional
> ordering appears in both streams." — talk_script.md

### The sign convention (the recurring gotcha — get this right)

In this geometry, **parallel (B/C, ∝ 1+3cos²θ) RAISES the coil load; perpendicular
(A, ∝ sin²θ) LOWERS it** (coil fast flux: A −15.3% free → −11.3% scatter; B/C +13.5%
free → +8.5% scatter). A magnet-protection benefit would come from the *perpendicular*
mode.

> "parallel (B/C, 1 + 3cos²θ) raises the coil load; perpendicular (A, sin²θ) lowers it.
> A magnet-protection benefit in this geometry would come from the perpendicular mode.
> (The bundle's RESULTS_tier8_conformal.md has a 'parallel reduces coil flux' typo — it
> is wrong; we do not repeat it.)" — README_spf_deck.md

This is important enough to have its own memory note: **RESULTS.md / the conformal
bundle carries a "parallel reduces coil flux" typo; it is wrong.** The `qa.md` Q&A
also records that the free-streaming TBR ≈ 0 is a *construction artifact* isolating pure
geometric steering (all real breeding, TBR ≈ 0.68–0.69/src, appears only in the scatter
stream), that `summary.json` in the bundle is empty (`{}`), and that `tallies_long.csv`
is the source of truth (the 1.9M-row `wall_current_phi` mesh is noisy and unused).

**Honesty framing** was enforced throughout: 10k-history smoke run, large MC error deep
behind the blanket → *"Machinery validation, not a reactor claim."*

---

## How the two threads relate

Thread 1 is the **rigorous, converged, analytic-anchored** validation (the two papers'
credibility figures): free-streaming, machine-to-1.8% MC↔analytic agreement on real and
toy QA sources. Thread 2 is an **earlier, indicative smoke run** with scattering and a
conformal blanket but only 10k histories — less rigorous, distinct in purpose, and the
origin of the η ≈ 0.63 "steering survives the wall" story and the coil-flux sign
convention. Both agree on the physics direction; only Thread 1 supports a
"validated-against-analytic" claim.
