# Retrospective — the tier arc (Tier 1 → Tier 8)

*A condensed, in-retrospect narrative of the verification tiers, synthesised from
the per-tier `RESULTS_*.md` files (now archived alongside this note in
`notes/results/`). Every substantive number is preserved; the original files hold
the full tables. Blockquotes are verbatim, attributed to their source file.*

The prototype was built as a **strict tier ladder**: no tier started until the one
below it was green. The organising discipline — inherited from the analytic anchor —
was that *every step past a closed-form ground truth must ship a replacement
verification*. That single rule shaped the whole arc.

---

## Tier 1 — the direction sampler alone (no transport)

**At the time** we did not yet trust any geometry; the question was only whether the
angular sampler emitted the Schwartz distribution exactly. **We built** a C++
standalone driver (`spf_driver`, using OpenMC's real `prn`) plus a 1:1 Python mirror
(`spf_mirror.py`) as the oracle, and gated the two against each other bit-for-bit.

**We achieved** an all-green gate:

- **sympy mode-algebra**: `∫dσ/dΩ dΩ = σ₀(a + ⅔b + ⅓c)` confirmed; corner total cross
  sections **2/3, 1, 2/3, 1/3** (nonpol, A, B, C) in σ₀ units; `∫sin²θ dΩ = 8π/3`,
  `∫(¼+¾cos²θ)dΩ = 2π`; and **B and C share the same angular shape** (bracket_B −
  2·bracket_C = 0).
- **C++ ↔ Python parity**: mixed mode, off-axis B̂, N=5000 → **max|C++ − Python| =
  6.106e-16**.
- **per-mode cosθ** (N=10⁶, B̂=ẑ): χ²/dof ≈ 1 for every mode (nonpol 0.972, A 1.324,
  B 0.747, C 1.370, mixed 0.698); ⟨cos²θ⟩ matched analytic (A 0.200, B/C 0.467,
  iso 0.333). φ uniform (χ²/dof 1.034).
- **rotation**: an isotropic source rotated by an off-axis B̂ stays isotropic
  (χ²/dof 0.985); pure-A steers ⊥ B̂; degenerate B̂=±ẑ handled (`max‖u‖−1` ≈ 4e-16).
- **inverse-CDF vs rejection** (two-sample KS): sin²θ D=3.1e-3 (p 0.102),
  ¼+¾cos²θ D=1.8e-3 (p 0.687) — the closed-form cubic root selection is correct.

> "Tier-1 gate: all checks pass." — RESULTS_tier1.md

This tier established the two design decisions the whole project rests on: **only two
distinct angular shapes** (`sin²θ` and `¼+¾cos²θ`), and a clean **direction-vs-rate
split** (the sampler emits at constant rate; η is applied later).

---

## Tier 2a — the independent analytic NWL

**At the time** we needed a trustworthy free-streaming ground truth before running any
transport. **We built** three *independent* analytic engines for the square-cross-
section torus (geometry exactly `anarrima`'s `square_torus.py`: R0=1.0, a=0.5, inboard
u=0.4, outboard w=1.6, floor/ceiling z=∓0.6 — **aspect 2.0**, matching the reference
implementation, though the paper *text* says 2.5): our own Gauss-Legendre **quad**,
Schwartz's published closed forms via **anarrima**, and our own transcription of the
inboard-isotropic closed form **Eq.5**.

**We achieved** machine-precision agreement across all three:

- single central ring: **max |Δ| = 1.8e-15** across all walls × factors;
- parabolic plasma, quad vs anarrima: **max rel diff 5.2e-08**;
- Schwartz §2 scalar oracles reproduced (paper text rounds): A inboard **+40.6%**
  (paper +43%), A outboard **−21.2%** (−22%), B/C inboard **−40.6%** (−43%), B/C
  outboard **+21.2%** (+22%), iso outboard/inboard **+14.6%** (+12%).

Crucially, this tier **corrected the spec**: the task's §4.3 claimed "iso ≡ C", but
the data shows **C/iso at inboard = 0.594 ≠ 1** — C is −40.6% vs iso, *like B*. The
valid cheap equivalence is **B ≡ C**, not iso ≡ C.

> "iso ≢ C (spec §4.3 says 'iso≡C' — this is the §1.2 error): C/iso at inboard =
> 0.594 (≠1; C is -40.6% vs iso, like B). The valid cheap equivalence is B≡C, not
> iso≡C." — RESULTS_tier2a.md

---

## Tier 2b — full-stack OpenMC recovers the analytic

**At the time** this was the headline free-streaming result. **We ran** OpenMC on the
same square torus with a near-void interior (Leakage=1, free-streaming), 4,000,000
histories/config, 50×50 poloidal bins/face, comparing **directionality** (mode/iso,
normalisation-free) against the exact constant-rate analytic `D = bracket/η`.

**We achieved** statistical agreement — the standardised residual `(MC−analytic)/σ_MC`
has **mean ≈ 0** across every config (A/iso −0.022, B/iso +0.019, C/iso +0.019,
mixed −0.009). The Schwartz oracles came out of *OpenMC itself* (A inboard **+39.2%**,
A outboard −21.0%, B/C inboard −40.9%, B/C outboard +21.1%, iso out/in +13.5%). B ≡ C,
iso ≢ C, and linearity all held.

> "OpenMC reproduces the analytic directionality to statistics; a coherent nonzero
> mean would signal a bug." — RESULTS_tier2b.md

**Honesty flag (recorded later, in LIMITATIONS/CONTEXT_REPORT, not here):** the
residual **stdev is ≈0.68, not 1.0** — the MC agrees with the analytic *better* than
the 10-batch tally error bars imply (mildly conservative error bars; no bias). Tier-2b
itself loosely states the pass criterion as "stdev≈1"; the honest statement is
"mean≈0, σ≈0.68". Separately, the B≡C check reads **exactly 0.000** because pure B and
pure C both have `a=0 ⇒ P_perp=0`, take the identical RNG branch, and (with no explicit
seed) produce **bit-identical particle streams** — it confirms the shape identity, not
statistical convergence.

---

## Tier 3 — presentation artifacts (postprocessing only)

**At the time** we needed engineering-facing figures/tables; **no new transport ran** —
this is pure postprocessing of the Tier-2b wall-current data + the Tier-2a analytic.

**We achieved** the peaking-factor finding that contradicts a naive expectation:

> "in this FIXED square geometry, both polarized modes RAISE the overall peaking
> factor vs isotropic (iso 1.33 → A 1.65, B/C 1.61). Directional emission concentrates
> the wall load, so the isotropic case is the most uniform." — RESULTS_tier3.md

The real B/C benefit here is **center-stack load reduction**, not a lower global
peaking factor: inboard peak −41%, and the inboard neutron-current fraction drops from
**8.7% (iso) → 5.1% (B/C)** while outboard gains 8.6 pts (42.6% → 51.2%); A mode does
the reverse. This was framed carefully:

> "This is the geometric precursor to a TBR argument — it is NOT a TBR: a real TBR
> needs a breeder blanket + scattering, absent here." — RESULTS_tier3.md

---

## The analytic anchor fades: Tiers 5–7 (scattering, TBR, absolute rate)

**At the time** we recognised the pivot point: *"The moment we add scattering, that
closed-form ground truth disappears."* (PLAN_next_phase.md). The design rule became
explicit — every new physics step ships a replacement verification (thin-wall recovers
free-streaming; TBR anchored by benchmarks + trends). Recommended order was **4b → 4a →
5 → 6**; in practice the delivered tiers are 5, 6, 7.

### Tier 5 — scattering + tungsten first wall

**We built** a layered box-torus (cavity → W 0.1 cm → Fe-9Cr steel 1 cm → Be 1 cm →
FLiBe 40 cm), ARC-class (minor radius ~1 m), scattering ON, 400k histories/config.

**The verification anchor held**: with all densities ×10⁻⁴ (thin wall) the layered
model recovered the Tier-2b free-streaming steering (A inboard +40.0% / outboard
−19.5%; B/C −40.5% / +19.2%). Then the key new physics:

> "free-streaming overestimates the steering benefit by ~3×. B/C still relieves the
> inboard center stack, but by ~15%, not the ~41% the geometric model suggested. This
> is exactly the kind of correction the analytic theory cannot provide." — RESULTS_tier5.md

With the *real* blanket, backscatter dilutes the steering: A inboard +40% → **+13.2%**;
B/C inboard −40% → **−15.0%**. Energy deposition: **FLiBe 91.4%**, Be 6.28%, steel
2.26%, W 0.09%; total heating = 1.22e7 eV/src = **0.87× the 14.1 MeV** (rest leaks the
unreflected boundary); a large down-scattered tail appears in the FLiBe flux — the
signature that scattering is active. Positioned as complementary to **Bae et al. 2025**
(+68% magnet lifetime), whose free-streaming-vs-transport ratio was not tabulated.

### Tier 6 — FLiBe blanket TBR

**We achieved** a defensible TBR anchored by trends (no analytic TBR exists):

- **TBR = 1.149 ± 0.002** at 30% Li-6 (ARC/LIBRA band; ARC ≥1.1, DEMO goal 1.15).
- Per mode (constant rate): iso 1.149, A 1.144 (−0.44%), B/C 1.155 (+0.57%).
- Enrichment: rises then plateaus — 7.5% → 1.067, 30% → 1.149, 50% → 1.154, 90% →
  1.135 (reproducing LIBRA's ~20–50% optimum, Peterson 2022).
- Li-6/Li-7 split (iso, 30%): 1.040 / 0.099.

The earlier strict "TBR-independent" claim was **softened**:

> "TBR is modestly polarization-dependent — here ~0.6% (B/C parallel-emitting higher,
> A perpendicular lower)… the strict 'TBR-independent' claim is wrong; Bae's +2.7% is
> the careful number." — RESULTS_tier6.md

This **reproduces the sign/trend of Bae et al. 2025** (parallel +2.7%, perpendicular
−1.8% in a Pb-Li spherical tokamak) in an independent code path, breeder, and geometry
— i.e. independent cross-validation, not a new result.

### Tier 7 — the absolute rate effect (rate restored)

Tiers 2b/3/6 held the fusion rate constant to isolate *directionality*. Tier 7
(postprocessing of the Tier-5/6 statepoints) **restores** the rate factor `η = a + ⅔b +
⅓c`: A → 1.50× power (+50%), B → 1.00×, C → 0.50× (−50%) — the Kulsrud (1982)
cross-section enhancement the constant-rate comparisons had normalised away.

At fixed fuel density: A gives +49% absolute tritium but steers inboard, so the
center-stack load compounds to **1.70×**; B is the **sweet spot** (same power, same
breeding, center-stack load **0.85×** — free steering, no rate cost); C halves
everything (0.42× center-stack). **Self-sufficiency is preserved**: TBR per neutron is
~mode-independent (1.144–1.155).

> "B (free steering): same power, same breeding, center-stack load 0.85x — the sweet
> spot." — RESULTS_tier7.md

---

## Tier 8 — pitched and 3D fields (the bridge to a real equilibrium)

**At the time** the goal was to make the source's steering axis a real, spatially
varying B̂(x) while keeping a verification anchor. Two parts.

### Part A — angled/pitched field (the rigor anchor)

Pitched field `B̂ = cosβ·φ̂ + sinβ(cosα·R̂ − sinα·ẑ)`, test pitch (α,β)=(0.6,0.5).
**We achieved** an exact analytic anchor:

- single-ring generalized quad vs anarrima's angled kernels g_*a: **3.6e-14**;
- parabolic plasma vs anarrima_vec: 2.4e-08;
- **β=0 recovers toroidal exactly**: max rel |I(angled,β=0) − I(toroidal)| = **0.0e+00**
  (bit-identical B̂ ⇒ identical RNG stream ⇒ exact recovery);
- free-streaming OpenMC reproduces the angled analytic (A inboard analytic +21.9% vs
  OpenMC +21.4%; outboard −16.3% vs −16.6%).

> "The verified pitched-field source is the bridge to a real equilibrium B̂(x)." —
> RESULTS_tier8_angled.md

### Part B — 3D-field-driven stellarator (end-to-end)

**At the time** the 3D physics lived *only in the source* (a 2-field-period field map
B̂(x) + flux-surface births) on an **axisymmetric** Helios-class radial-build rig with
**public placeholder** geometry. This was framed unambiguously:

> "Machinery validation with PUBLIC placeholder geometry — NOT a Helios physics claim."
> — RESULTS_tier8_3d.md

**We achieved** a validated end-to-end pipeline, anchored by:
- **emission moments** about the local 3-D B̂ (⟨(u·B̂)²⟩ = 0.3328/0.2007/0.4668 for
  unpol/perp/parallel vs analytic 0.3333/0.2000/0.4667);
- **B̂ genuinely varies** over the source (std of B̂_z across births = 0.092 ≠ 0);
- **reduce-to-axisymmetric**: field-map toroidal vs trusted bmode=toroidal agree to
  max **0.8 pts**.

End-to-end (scattering, rate fixed): global TBR 1.368/1.366/1.373; first-wall steering
(inboard ±~2.2%, outboard ±~7%); and a **directional efficiency η = 0.89** — the
fraction of the axisymmetric steering retained under the pitched 3-D field. The deep
coil fast-flux tally was left at ~30% error (variance reduction deferred).

---

## The joint analytic ↔ MC validation ladder (framing)

`JOINT_VALIDATION.md` records the protocol that ties the analytic companion paper to
this MC pipeline: both consume the *same* equilibrium; conventions (φ, polarization
K(cosθ) with cosθ = B̂·Δ̂) were audited against anarrima and confirmed to match.

> "Because the analytic kernel is exact in the free-streaming limit, it is the quantity
> a Monte Carlo neutronics calculation with a polarized source must reproduce, per wall
> patch and per polarization mode." — JOINT_VALIDATION.md

The standing caveat: a **numbers-level precise_QA** comparison is not yet possible —
the analytic Tier-3 driver still runs a toy spectrum + model field. *"Until that swap
is done, do NOT treat the analytic Tier-3 outputs as precise_QA ground truth."* Key the
comparison on each patch's (R, Z, n̂), not on index.

---

## What the tier arc adds up to

Machine-precision sampler (C++/Python parity ≤1e-12) → machine-precision analytic NWL
(three engines to 1e-15) → OpenMC recovers it to statistics → scattering dilutes the
geometric steering ~3× → a defensible FLiBe TBR reproducing Bae's polarization sign →
the absolute +50%/−50% rate trade-off → a pitched/3-D field source anchored to
anarrima and reduce-to-axisymmetric. The honest overall verdict, carried in the
handoffs: **independent validation + a reusable verified tool, not new physics.**
