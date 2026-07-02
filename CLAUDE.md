## OpenMC Codebase Tools

Read the FULL `AGENTS.md` in this directory before starting work. It contains
project context, coding conventions, and documentation of the RAG search tools
registered in `.mcp.json`.

# Claude Code task: prototype a spin-polarized DT neutron source sampler in OpenMC

You are operating inside my cloned OpenMC repository. Your job is to **plan, then implement, then verify** a prototype polarization-aware fusion neutron *source* (the birth-direction sampler), validate it against the unpolarized isotropic limit, and validate it against the analytic axisymmetric results of Schwartz (2025). This is a research prototype standing on its own — it does **not** depend on any private MIT codebase. Where my collaborator's group has an internal physics-informed tokamak source, that is a later integration target; here we build a self-contained `CompiledSource` and a worked analytic test problem so the physics is correct before any integration.

Work in the style I prefer: **hard-gate on correctness, no warnings-where-an-assert-belongs, honest non-cherry-picked diagnostics, strict tier ordering** (do not start tier N+1 until tier N passes). Do not over-engineer. The MVP is: unpolarized verification → single-mode angular samplers → mixture sampler → frame rotation → analytic NWL comparison in a square-cross-section torus. Alpha emission, spatially varying polarization, realistic equilibria, and stellarator geometry are explicitly out of scope for this prototype.

---

## 0. Orientation (do this first, write findings to a scratch file, do NOT write code yet)

Before touching anything, establish ground truth about *this* checkout, because the source API changed across versions and I need you anchored to what's actually here, not to your training data.

1. `git rev-parse HEAD`, `git status`, and read `include/openmc/version.h` (or equivalent) to record the exact version/commit.
2. Read these files and report back what the real signatures are in this checkout (do not assume — the API was refactored around 0.14/0.15):
   - `include/openmc/source.h` — confirm the `openmc::Source` abstract base class, the exact signature of the virtual `sample(uint64_t* seed) const` method, and the return type.
   - `include/openmc/particle_data.h` (or wherever `SourceSite` lives) — confirm the exact field names and types of `SourceSite`. I expect roughly `r` (Position), `u` (Direction), `E` (double, eV), `wgt`, `delayed_group`, `particle` (ParticleType), but **verify** field names and the `Position`/`Direction` types and their member access (`.x/.y/.z`).
   - `include/openmc/random_lcg.h` — confirm the RNG call is `openmc::prn(seed)` returning a double in [0,1).
   - `openmc/source.py` — confirm `CompiledSource(library, parameters=..., strength=..., constraints=...)` and how `parameters` (a string) reaches the C++ `openmc_create_source(std::string)` factory.
   - Find at least one in-tree example of a compiled source if present (search `tests/`, `examples/`, `docs/` for `openmc_create_source`). If you find the build pattern (CMake target, or a standalone `gcc`/`g++` compile line linking against `libopenmc`), record it verbatim. The custom-source build is the single most fragile step; get the real recipe from this repo, not from memory.
3. Confirm cross-section data availability: is `OPENMC_CROSS_SECTIONS` set, and is there an ENDF/B HDF5 library reachable? For the analytic comparison we will run with scattering effectively disabled, so we need *some* valid data file to load but the wall material will be near-void.
4. Record the OpenMP threading model note: the `sample()` method is called concurrently by multiple threads, each with its own `seed`. **Any mutable member state in the source class is a data-race bug.** All sampler state must be set at construction and then const.

Write all of this to `spf_prototype/NOTES_orientation.md`. Stop and show me this file's key findings (especially the exact `SourceSite` fields and the custom-source build recipe) as part of your plan before implementing.

---

## 1. The physics you must implement — and a correction you must apply

I have a personal foundation document on this project. **Trust the Schwartz primary-source formulas below over any half-remembered version**, including over my own notes, which contain at least one transcription error. Here is the corrected, authoritative physics. Re-derive anything you're unsure of; do not just copy.

### 1.1 Polarization mode parameterization (Schwartz 2025, Eq. 1)

Deuteron spin-projection fractions `d₊, d₀, d₋`; triton fractions `t₊, t₋`. The three collision-mode fractions are:

```
a = d₊·t₊ + d₋·t₋
b = d₀                = d₀·t₊ + d₀·t₋
c = d₊·t₋ + d₋·t₊
```

with `a + b + c = 1`. Non-polarized fuel: `a = b = c = 1/3`.

Note: `b` is just the deuteron transverse (m=0) fraction. This `(a,b,c)` grouping is **not** the naive "S=3/2 fully-aligned / S=3/2 mixed / S=1/2" grouping — be careful, because a common secondary-source mistake (and one in my own notes) is to mis-assign which mode is which and to claim mode c is isotropic. It is not. See below.

### 1.2 Differential cross section (Schwartz 2025, Eq. 2) — THIS IS THE AUTHORITATIVE FORM

```
dσ/dΩ = (σ₀ / 2π) · [ (3/4)·a·sin²θ  +  ( (2/3)·b + (1/3)·c )·( 1/4 + (3/4)·cos²θ ) ]
```

where `θ` is the angle between the emitted neutron direction and the **local magnetic field direction B̂**, and `σ₀` is the cross section in the S=3/2 combined-spin state.

Key consequences (use these as sanity checks and as expected diagnostic outcomes — **memorize these, they are your oracle**):

- **Total (integrated) cross section** as a function of mode:
  ```
  σ_tot / σ₀  =  a  +  (2/3)·b·(... )  ...   →  integrate Eq. 2 over 4π:
  σ_tot = σ₀ · [ a + (2/3)b + (1/3)c ]
  ```
  Verify this integral yourself symbolically (sympy) as part of the build. Check the corners:
  - Non-polarized `(1/3,1/3,1/3)`: `σ_tot = σ₀·(1/3 + 2/9 + 1/9) = σ₀·(2/3)`. **Non-polarized integrated cross section is 2σ₀/3**, NOT σ₀.
  - Pure **A mode** `(1,0,0)`: `σ_tot = σ₀`. That is the **50% enhancement** relative to 2σ₀/3 (σ₀ / (2σ₀/3) = 3/2). Emission ∝ `sin²θ` (peaks perpendicular to B̂, zero along B̂).
  - Pure **B mode** `(0,1,0)`: `σ_tot = σ₀·(2/3)` = same as non-polarized. Emission ∝ `1/4 + (3/4)cos²θ` (peaks parallel to B̂).
  - Pure **C mode** `(0,0,1)`: `σ_tot = σ₀·(1/3)` = **half** the non-polarized rate. Emission has the **same angular shape as B mode**, `1/4 + (3/4)cos²θ` — NOT isotropic. (The B and C modes share directionality; they differ only in total rate. This is the correction to my notes — verify it directly from Eq. 2 by setting a=0 and seeing that the bracket's θ-dependence is identical for the b and c terms.)

- **Angular-shape decomposition you will actually sample.** Because the b and c terms share the same angular function, the normalized birth-direction pdf depends only on two effective weights. Define
  ```
  w_perp  = (3/4)·a                      (the sin²θ piece)
  w_par   = (2/3)·b + (1/3)·c            (the 1/4 + 3/4 cos²θ piece)
  ```
  Then the **unnormalized** angular intensity is `I(θ) = w_perp·sin²θ + w_par·(1/4 + 3/4 cos²θ)`, azimuthally symmetric about B̂. The birth-direction pdf on the sphere is `p(θ,φ) = I(θ) / ∫I dΩ`, with `φ ~ U[0,2π)`. **Sampling the birth direction only needs (w_perp, w_par); the overall total-rate normalization (which distinguishes A from non-pol from C) is a separate scalar that multiplies the spatial source strength, not the direction pdf.** Keep these two concerns separate in the code — this separation is the single most important design decision and it mirrors how Schwartz himself separates the angular factors from the rate factors.

### 1.3 Energy (decoupled, standard)

Birth energy is decoupled from the polarized angular distribution to leading order. For the **analytic Schwartz comparison**, emit monoenergetic 14.1 MeV (`14.06e6`–`14.1e6` eV; pick one and be consistent) — Schwartz drops the energy spread entirely. For a more physical mode (optional, behind a flag), sample a Ballabio-broadened Gaussian at a user-given ion temperature; mean ≈ 14.08 MeV, width ≈ `177·sqrt(T_i[keV])` keV. Energy is irrelevant to the angular verification, so do not let it block you.

---

## 2. What to build (proposed structure — adjust to match what you find in §0, and propose your final layout in the plan)

Create a self-contained subdirectory `spf_prototype/` so nothing pollutes the OpenMC tree. Proposed contents:

```
spf_prototype/
  NOTES_orientation.md          # from §0
  src/
    polarized_fusion_source.cpp # the CompiledSource
    CMakeLists.txt              # or a documented standalone compile command
  python/
    abc_modes.py                # (a,b,c) <-> (w_perp,w_par), total-rate factor, all with sympy-verified identities
    analytic_nwl.py             # independent reimplementation of the square-torus NWL integral (see §4)
    build_and_run.py            # drives compilation + OpenMC runs for each config
    verify_sampler.py           # §3 standalone direction-sampler tests (no transport)
    verify_nwl.py               # §4 full-stack NWL comparison vs analytic
  tests/
    test_abc_modes.py           # pytest: identities, corner cases
    test_sampler_stats.py       # pytest: chi-square / KS on sampled directions
  README.md                     # how to reproduce every number and figure
```

### 2.1 The C++ sampler — required design points

- Class `PolarizedFusionSource : public openmc::Source`, constructed from a parsed parameter string. Parameter string format (parse it robustly): `"a=0.5,b=0.3,c=0.2,bmode=toroidal,r0=600,z0=0"` with lengths in **cm** (OpenMC convention). Support at minimum: the three weights, a B̂-field specification (start with two options: a global constant direction, and a "purely toroidal" field `B̂ = φ̂(position)`), and the ring-source geometry `(r0, z0)` for the test problem.
- At construction: validate `a,b,c ≥ 0`, renormalize to sum 1 (warn to stderr if the input sum deviates by >1e-6 — this is a legitimate warn, it's user input, not an invariant). Precompute `w_perp`, `w_par`, and the mode-selection probabilities. **Assert** internal invariants (e.g. selection probabilities sum to 1 within 1e-12) — that is an invariant, so hard-assert, don't warn.
- `sample(uint64_t* seed) const override`:
  1. Sample birth position. For the test problem this is a filamentary ring: `φ = 2π·prn(seed)`, `r.x = r0·cos φ`, `r.y = r0·sin φ`, `r.z = z0`. (Keep position sampling pluggable so a volumetric source can replace it later.)
  2. Determine local `B̂` at that position (constant, or toroidal `φ̂ = (−sin φ, cos φ, 0)`).
  3. Sample a direction in the **B̂-aligned local frame** (`ẑ_local ∥ B̂`) from the mixture (see §2.2).
  4. Rotate the local direction to global Cartesian (see §2.3). This MUST be unit-normalized; assert `|u|−1 < 1e-12`.
  5. Set energy (monoenergetic for the verification path).
  6. Return the `SourceSite` with `wgt = 1.0`, `particle = neutron`, `delayed_group = 0`.
- **RNG**: only `openmc::prn(seed)`. Never `rand()`, never `<random>` engines. Reproducibility and thread-safety depend on it.
- **Thread safety**: no mutable members. Everything set in constructor, all sampling methods `const`.

### 2.2 Angular sampling within the local frame (stratified, exact)

Implement **stratified sampling** over the two angular shapes — it gives exact samples from the mixture and makes per-mode verification trivial:

- With probability `P_perp = W_perp / (W_perp + W_par)` choose the `sin²θ` shape; else the `1/4+3/4 cos²θ` shape. Here `W_perp`, `W_par` are the **solid-angle-integrated** weights of each shape (i.e. `w_perp` and `w_par` times the integral of their respective angular functions over the sphere). Compute these integrals symbolically and write the closed forms into a comment block; do not hand-wave the constants. (`∫ sin²θ dΩ = 8π/3`; `∫ (1/4 + 3/4 cos²θ) dΩ = 4π·(1/4) + ... ` — derive it and put the number in.)
- **sin²θ shape**: marginal in `cosθ` is `∝ sin³θ`. Provide BOTH:
  - the closed-form inverse-CDF via the depressed-cubic / trig method (solve `x³ − 3x + 2(1−2u) = 0`, `x = cosθ`, pick the root in [−1,1]); and
  - a rejection fallback (envelope `∝ sinθ`, acceptance 2/3) used to **cross-check** the inverse-CDF in a unit test.
- **`1/4 + 3/4 cos²θ` shape**: marginal in `cosθ` is `∝ (1/4 + 3/4 x²)`. Inverse-CDF is a cubic in `x` too; implement and cross-check with rejection.
- `φ ~ U[0,2π)` for both.
- **Decouple direction from rate.** The stratified split above samples the *direction* correctly for any `(a,b,c)`. Separately compute the scalar total-rate factor `η = a + (2/3)b + (1/3)c` (relative to σ₀) and expose it via the source `strength`/spatial weighting so that, e.g., A mode and C mode produce different *numbers* of neutrons but each produces *correctly shaped* directions. For the pure-shape angular verification in §3 you normalize it out; for the NWL comparison in §4 it matters only if you compare absolute (not relative) wall loads — Schwartz mostly shows relative patterns and even rescales the A mode by 2/3 in his Fig. 2 to isolate directionality, so plan to compare normalized patterns first.

### 2.3 Local-frame → global rotation (robust, no gimbal pitfall)

Do NOT build the basis as `B̂ × ẑ_world` (degenerate when B̂ ∥ ẑ). Use the "least-aligned reference axis + Gram-Schmidt" construction:

- pick `ref` = whichever of x̂,ŷ,ẑ has the smallest `|B̂·ê|`;
- `x̂_local = normalize(ref − (ref·B̂)B̂)`, `ŷ_local = B̂ × x̂_local`, `ẑ_local = B̂`;
- `u_global = local.x·x̂_local + local.y·ŷ_local + local.z·B̂`.

Unit test the rotation in isolation (see §3.3).

---

## 3. Verification tier 1 — the direction sampler ALONE (no OpenMC transport)

Gate: **this whole section must pass before you build any OpenMC geometry.** These tests call the sampling math directly (either compile a tiny standalone C++ driver that prints sampled directions, or — cleaner — also implement the identical sampler in `python/` and test that, then assert the C++ matches the Python on a fixed seed sequence). Prefer the latter: a Python mirror of the sampler is the cheapest way to get a trustworthy oracle, and a "C++ vs Python on identical seeds" test catches transcription bugs.

3.1 **Mode identities (`tests/test_abc_modes.py`)**: symbolically (sympy) verify `∫ dσ/dΩ dΩ = σ₀(a + 2b/3 + c/3)`; verify the four corner total cross sections (2/3, 1, 2/3, 1/3 in σ₀ units); verify B and C share angular shape; verify `w_perp/w_par` and the stratified selection probabilities are self-consistent.

3.2 **Sampled-direction distribution (`tests/test_sampler_stats.py`)**: for each of {non-pol, pure-A, pure-B, pure-C, and a mixed `(0.5,0.3,0.2)`}, draw ~1e6 directions with B̂ = ẑ, histogram in `cosθ` over equal-solid-angle bins, and chi-square test against the analytic `p(cosθ) ∝ I(θ)`. Report the χ²/dof and p-value. **Show me the actual numbers**, not "it passed." Non-pol and the isotropic check must be flat; pure-A must vanish at `cosθ=±1`; pure-B/C must peak at `cosθ=±1`. Also confirm `φ` is uniform.

3.3 **Rotation isotropy + steering (`tests/test_sampler_stats.py`)**: (i) sample isotropic in local frame, rotate with a random non-axis B̂, confirm global distribution is still isotropic (χ² on the sphere) — catches rotation bias. (ii) sample pure-A with a chosen B̂ and confirm the global emission lobe is perpendicular to that B̂ (max density in the plane ⊥ B̂). (iii) Explicitly test the degenerate case B̂ = ±ẑ to confirm the Gram-Schmidt branch doesn't blow up.

3.4 **Inverse-CDF vs rejection cross-check**: both methods must produce statistically identical `cosθ` distributions for each shape (two-sample KS). This is how you know the closed-form cubic root selection is correct.

Write a short `spf_prototype/RESULTS_tier1.md` with the tables of χ²/dof, p-values, and KS statistics. No cherry-picking: if a test is marginal, say so and investigate.

---

## 4. Verification tier 2 — full-stack NWL vs Schwartz analytic (square-cross-section torus)

Gate: tier 1 green. This is the headline result. We reproduce Schwartz's square-cross-section-torus example (his Section 2 / Figure 2) with OpenMC running the compiled source, scattering disabled, and compare the wall-load pattern to an **independent** analytic computation.

### 4.1 Independent analytic NWL (`python/analytic_nwl.py`)

Do **not** trust a single implementation. Reimplement the filamentary-ring NWL integral yourself from Schwartz Eq. 3 specialized to a rectangular cross section:

```
Q_NWL(r_w) ∝ ∮ dφ  [ angular_factor(θ_s(φ)) ] · ( n̂ · (r_s−r_w) / |r_s−r_w|³ )
```

where for a ring at `(p, z0)` and a wall point `r_w`, `r_s(φ) = (p cosφ, p sinφ, z0)`, `θ_s` is the angle between the emission direction `(r_w − r_s)/|...|` and the local B̂ at the source (toroidal: `φ̂(φ)`), and `angular_factor` is `sin²θ` (A), `1/4+3/4cos²θ` (B/C), or `1` (isotropic). Integrate over only the toroidally-visible arc. Do this two ways and require agreement to ~1e-6:
  - (a) direct numerical quadrature in φ (scipy `quad`), and
  - (b) — **strongly preferred as the real oracle** — install and call **anarrima**, Schwartz's own open-source package (`pip install anarrima` or from its repo; it's referenced as ref [44] in arXiv:2512.09242 and §5 of arXiv:2507.11758). If anarrima installs cleanly, use it as the primary analytic reference and keep your own quadrature as a secondary cross-check. If it does not install, proceed with your quadrature and note it.

Reproduce, as scalar checks against the paper's text: outboard wall ~12% hotter than inboard for isotropic (pure geometry); A mode **+43%** on the inboard midplane and **−22%** on the outboard midplane; B/C modes the mirror image (−43% inboard, +22% outboard). These specific numbers are stated in Schwartz §2 and are your pass/fail oracle.

### 4.2 OpenMC model (`python/build_and_run.py`, `verify_nwl.py`)

- Geometry: square-cross-section torus. Use Schwartz's example proportions — circular plasma of aspect ratio 2.5 in a square vessel; place the **filamentary ring** at the plasma center for the first comparison (a single ring reproduces Fig. 2's "circular plasma" only after you sum rings weighted by a parabolic `(1−ρ²)` profile, so do BOTH: first a single central ring vs your single-ring analytic, then the parabolic-weighted multi-ring sum vs anarrima's circular-plasma result). Lengths in cm.
- **Scattering off**: set the wall/interior material to an extremely low total cross section (near-void), OR run a sampling-only / current-tally configuration so neutrons free-stream from ring to wall. The point is to isolate source-direction physics. Confirm in the log that essentially no collisions occur (check the collision tally / mean free path is ≫ device size).
- Tally: surface current on each of the four walls, binned poloidally (a mesh or surface filter giving ≥50 bins per wall). Normalize each pattern to its own integral so you compare *shapes* first (matching Schwartz's normalized presentation and his 2/3 rescale of A).
- Configs to run: {isotropic non-pol, pure-A, pure-B, pure-C, mixed `(0.5,0.3,0.2)`}. Particle count high enough that per-bin statistical error ≪ the effect (target ≤2–3% per bin; ~1e8 histories if it's affordable, fewer if not — report what you used and the resulting per-bin error).

### 4.3 Pass/fail and the linearity test

- **Primary criterion**: for each config, `(Q_MC − Q_analytic)/σ_MC` across wall bins is a standard normal — histogram it, report mean (should be ~0) and stdev (should be ~1). A coherent nonzero mean signals a real bug (wrong constant, wrong rotation), not noise.
- **Isotropic ≡ pure-C shape**: their normalized wall patterns must be identical within statistics (they share directionality). This is a strong, cheap correctness check and a direct test of the §1.2 correction.
- **Linearity**: the `(0.5,0.3,0.2)` normalized pattern must equal `0.5·A + 0.3·B + 0.2·C` of the pure-mode patterns within statistics. If linearity fails, the mode-selection probabilities are wrong.
- Reproduce the scalar oracles from §4.1 (the ±43%/±22% inboard/outboard midplane numbers) from the **OpenMC** result, not just the analytic.

Write `spf_prototype/RESULTS_tier2.md` with: the residual histograms (mean/stdev per config), the isotropic-vs-C equivalence, the linearity check, the scalar-oracle comparison table (paper vs analytic vs OpenMC), and the 1D outboard-equator line plots overlaying all configs. Save figures to `spf_prototype/figs/`.

---

## 5. Diagnostics / figures to produce (honest, non-cherry-picked)

1. Per-mode `cosθ` histograms with analytic overlay (tier 1).
2. 2D wall-load maps per config on each wall face (tier 2).
3. 1D outboard-equator and inboard-equator line plots, all configs overlaid, symmetric in Z.
4. Residual histograms `(MC−analytic)/σ` per config.
5. A single summary table: for each config — total-rate factor η, peak/avg per wall, inboard & outboard midplane % change vs isotropic, χ²/dof of the residuals. Put paper-stated values in the same table where they exist.

Every figure must be regenerable by one script with a fixed seed. State the OpenMC version/commit and the exact run parameters in `README.md`.

---

## 6. Deliverables and how to report back to me

Proceed in this order and **pause for my sign-off after the PLAN and after each tier**:

1. **PLAN** (no code yet): the §0 orientation findings (especially real `SourceSite` fields + the custom-source build recipe from this checkout), your final directory layout, and any place where what you found in the repo differs from the assumptions in this prompt. Flag anything that would change the approach.
2. **Tier 0+1**: the C++ sampler + Python mirror + tier-1 tests green, with `RESULTS_tier1.md`. Show me the actual χ²/p-value/KS tables.
3. **Tier 2**: NWL comparison green, with `RESULTS_tier2.md`, the residual-normality results, the isotropic≡C and linearity checks, and the scalar-oracle table.
4. A top-level `spf_prototype/README.md` that lets me (or a collaborator) reproduce every number from scratch, plus a short `LIMITATItONS.md` listing what this prototype does NOT do (no scattering in the verification, single energy, filamentary/weighted-ring source only, constant or toroidal B̂ only, no alpha channel, no depolarization, convex geometry only — matching Schwartz's own caveats).

Constraints / style reminders:
- Hard-gate on invariants (assert), warn only on user input. Strict tier ordering.
- No mutable state in the source class; `openmc::prn` only.
- Verify the integrals symbolically; do not paste magic constants without a derivation comment.
- Treat the Schwartz Eq. 1–2 above as authoritative and **re-derive** the integrated cross sections yourself to confirm the corner values (2/3, 1, 2/3, 1/3) before trusting any sampler output.
- Honest diagnostics: report marginal/failed tests plainly and investigate; do not tune particle counts or binning to make a plot look clean.
- If the custom-source build against this checkout fights you, surface the exact error and the build command rather than silently switching approaches.

Begin with §0 and the PLAN. Do not write the sampler until I sign off on the plan.

### Claude Code-specific: first-call behavior

The first `openmc_rag_search` call of each session returns an index status
message instead of search results. When this happens, you MUST use the
`AskUserQuestion` tool to present the rebuild/use-existing choice to the user.
Do not ask conversationally — always use the widget. Do not skip this step even
if the index looks current — the user may have uncommitted changes that warrant
a rebuild.

FINALLY: DO NOT PUSH COMMITS TO MAIN MAKE MY OWN BRANCH AND SAVE COMMITS ON THAT! DO YOU COPY!??