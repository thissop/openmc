# HANDOFF / CONTEXT SUMMARY — spin-polarized fusion source prototype

Written to survive a context compaction. Authoritative state of the project as of
this session. Read this first to continue.

## 1. What this is
A self-contained prototype of a **spin-polarized DT neutron *source* for OpenMC**,
validated against the unpolarized limit and Schwartz (2025), then extended to
scattering + a FLiBe breeder blanket. Lives entirely in `spf_prototype/`. Research
prototype; does NOT depend on any private codebase (Peterson-source integration is
a future target).

PI: **Ethan Peterson** (MIT PSFC neutronics). User collaborates/joining; will email
him results. User style (HARD CONSTRAINTS): honest non-cherry-picked diagnostics;
hard-gate invariants (assert), warn only on user input; strict tier ordering;
**commit to branch `spf-prototype` only, NEVER push, NEVER to main/develop, do not
commit `CLAUDE.md`**; commit msgs end `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

## 2. Environment (all working; recipe also in NOTES_orientation.md §0.9)
- Platform aarch64 Linux. `$HOME=/home/tkiker.guest`; repo `/Users/tkiker/Documents/GitHub/openmc` (host mount); branch **spf-prototype**.
- **venv**: `$HOME/spf_venv` — openmc (editable), numpy/scipy/sympy/h5py/pandas/matplotlib/pytest, anarrima(+jax), smplotlib.
- **OpenMC C++** built + installed INTO the venv prefix: `$HOME/spf_venv/bin/openmc`, `lib/libopenmc.so`, `lib/cmake/{OpenMC,fmt,pugixml}`. (Build-tree config is NOT consumable; must `make install`.)
- **Cross sections**: `$HOME/nndc_hdf5/cross_sections.xml` (NNDC, 966 nuclides). Set via `openmc.config['cross_sections']` or `build_and_run.XS`.
- **Custom source .so**: `spf_prototype/src/build/libpolarized_fusion_source.so` — build via `build_and_run.build_so()` (cmake `-DCMAKE_PREFIX_PATH=$HOME/spf_venv` + make). CMakeLists: `find_package(OpenMC) + target_link_libraries(... OpenMC::libopenmc)`.
- **Run pattern**: `PATH=$HOME/spf_venv/bin:$PATH OMP_NUM_THREADS=2 $HOME/spf_venv/bin/python <script>`.
- git identity (local): user.email tjk2147@columbia.edu, user.name tkiker.
- `.gitignore` (spf_prototype): build/, *.pyc, *.so, *.pdf, *.tgz, latex aux; **figs/*.png are force-tracked** (override repo *.png ignore).
- **/tmp statepoints are EPHEMERAL** (Tier 5/6 runs: `/tmp/spf_t56_*`, Tier 1 driver, etc.). Regenerate via the verify_* scripts. A bundle exists: `spf_prototype/tier2b_statepoints.tgz` (gitignored).

## 3. Commit history (branch spf-prototype, nothing pushed)
b8db20c09 Tier0+1 · ea6482aff tier1 fig · 3cf0adfec Tier2a · a3ffc027e Tier2b+README/LIMITATIONS ·
ab0070663 Tier3 · 7248b21eb beamer deck · 802ae9c69 deck rewrite (non-expert) · e44cdb786 Tier5+6 ·
5e24e88d4 Tier7 + deck tiers5-7.
**Uncommitted (user's in-progress, leave unless asked):** `import smplotlib` in python/tier3_artifacts.py + restyled figs/tier3_*.png. Also CLAUDE.md (M, task spec — never commit).

## 4. File map (spf_prototype/)
- NOTES_orientation.md — §0 OpenMC API ground-truth + verified env/build recipe (§0.9).
- src/: `spf_sampler.hpp` (pure-math core), `polarized_fusion_source.cpp` (the .so), `standalone_driver.cpp` (tier-1, no libopenmc), CMakeLists.txt.
- python/: `spf_mirror.py` (PCG prn port + sampler mirror), `abc_modes.py` (sympy identities), `spf_stats.py`, `analytic_nwl.py` (quad + anarrima + Eq.5; parabolic plasma), `build_and_run.py` (free-streaming torus, .so build, extract_walls; consts U=0.4,W=1.6,ZW=0.6,R0=1,AM=0.5), `reactor_model.py` (Tier5/6 layered reactor), `verify_sampler.py`/`verify_nwl_analytic.py`/`verify_nwl.py`/`verify_reactor.py`/`tier3_artifacts.py`/`tier7_absolute.py` (diagnostics→RESULTS+figs).
- tests/: conftest.py, test_abc_modes.py, test_sampler_stats.py, test_analytic_nwl.py (56 pass tier1+2a).
- RESULTS_tier{1,2a,2b,3,5,6,7}.md · figs/*.png · slides/spf_talk.tex (52-slide deck, compiles w/ pdflatex).
- README.md, LIMITATIONS.md, NEXT_STEPS.md, PLAN_next_phase.md, BLANKET_PARAMS.md (cited), HANDOFF.md (this).

## 5. Physics/API facts (don't re-derive)
- Eq1 modes: a=d+t+ + d-t-, b=d0, c=d+t- + d-t+ (sum 1; nonpol 1/3 each).
- Eq2: dσ/dΩ=(σ0/2π)[¾a sin²θ + (⅔b+⅓c)(¼+¾cos²θ)], θ to local B̂.
- η=a+⅔b+⅓c (=σtot/σ0); corners nonpol 2/3, A 1, B 2/3, C 1/3. w_perp=¾a, w_par=⅔b+⅓c, P_perp=a/η.
- **B≡C** (share ¼+¾cos²θ); **C is NOT isotropic** (key correction; spec's "iso≡C" wrong).
- inverse-CDF cubics: perp p(x)∝(1−x²) → x³−3x+(4u−2)=0, root x=2cos(acos(1−2u)/3−2π/3); par p=¼+¾x² → x³+x+(2−4u)=0 (Cardano).
- Gram-Schmidt rotation (least-aligned axis) → works for any B̂.
- OpenMC: override `SourceSite sample(uint64_t* seed) const`; `.so` exports `extern "C" unique_ptr<openmc::Source> openmc_create_source(std::string)`; `ParticleType::neutron()` (a FUNCTION); `prn(uint64_t*)`∈[0,1); STREAM_SOURCE=1. Python `CompiledSource(library,parameters,strength)`.
- Geometry: Tier2/3 normalized R0=1,a=0.5,u=0.4,w=1.6,z=0.6 (aspect 2.0 = anarrima example; paper text says 2.5). Tier5/6 reactor (cm): R0=200,a=100,U=80,W=320,ZC=120; layers W0.1/Fe9Cr-steel1/Be1/FLiBe40 cm, vacuum outer.
- Tritium score: **'H3-production'** (appears as '(n,Xt)' in pandas df; cell tally cols are FLAT: cell,nuclide,score,mean,std. dev.).
- anarrima angled-field kernels exist: g_HAa/g_Hca/g_HBa, g_VAa/g_Vca/g_VBa (args p,z,r,α,β,φ).
- FLiBe = Li2BeF4 (atoms Li:Be:F=2:1:4), ρ=1.94 g/cm³, Li-6 enrich via add_element('Li',2,enrichment=E,enrichment_target='Li6').

## 6. Results (key numbers)
- T1: C++↔Py parity 6e-16; per-mode cosθ χ²/dof≈1; ⟨cos²θ⟩ A 0.2/BC 0.467/iso 0.333; corners 2/3,1,2/3,1/3.
- T2a: quad=anarrima=Eq.5 to 1.8e-15; plasma quad-vs-anarrima 5e-8; A inboard +40.6%/outboard −21.2%, iso out/in +14.6% (paper text rounds 43/22/12).
- T2b: directionality residual mean≈0, σ≈0.69 (MC recovers analytic); B≡C; linearity holds.
- T3: peaking iso 1.33 → A 1.65 / B-C 1.61 (RAISES; isotropic most uniform); inboard current fraction iso 8.7% → B/C 5.1% (geometric precursor, NOT TBR).
- T5: thin-wall recovers free-stream steering (A +40%/B −40%); **scattering DILUTES steering ~3×** (A inboard +40%→+13%, B/C −40%→−15%); heating FLiBe 91%/Be 6%/steel 2%/W 0.1%; total 0.87× neutron E; down-scattered spectrum tail.
- T6: TBR=1.149 (30% Li-6, ARC/LIBRA band); per mode A 1.144/B,C 1.155 (±0.6%); enrichment peaks ~50%, drops 90% (reproduces 20-50% optimum); Li6 1.04 + Li7 0.10.
- T7: at fixed density A 1.5× rate/1.49× tritium/**1.70× center-stack load**; B free steering (0.85× center-stack, 1× power); C 0.5× all; **TBR/neutron unchanged** (self-suff preserved). Rate (A) vs steering (B/C) are different spin states.

## 7. LITERATURE / NOVELTY (critical — recalibrated this session)
**Most of our work is NOT new physics.** Key prior work:
- **Bae et al., Nucl. Fusion 65, 086051 (2025)** "Neutronics analysis of spin-polarized fuel in spherical tokamaks" — **OpenMC (40M particles)**, TBR vs polarization (parallel/anti-aligned +2.7% = 1.111 vs 1.082; perpendicular lower), inboard/outboard TBR variation, **+68% magnet (center-stack) lifetime**. **This scoops our Tier 5/6.** (Also OSTI 2583820.)
- Schwartz arXiv:2507.11758 (2025) — analytic free-streaming NWL (we reproduce).
- Kulsrud PRL 49,1248 (1982) — +50% rate, perpendicular emission. Ciullo (ed.) Springer 2016 — rate-OR-directionality. Segantin FED 154,111531 (2020) — 20-50% Li-6 optimum. arXiv:2502.15941 (2025) — power enhancement. Peterson LIBRA (FST 2022) — FLiBe immersion blanket.
- Verdict: our work = **independent validation + reusable verified tool**, not discovery. We independently reproduce Schwartz (analytic+MC) AND Bae's TBR-polarization trend (our parallel-emitting B/C higher TBR, perp A lower — same sign). Possible thin novelty: explicit **free-streaming-vs-transport "analytic overestimates steering ~3×"** quantification (not found in lit; verify vs Bae).

## 8. CURRENT TASK (in progress when this was written)
1. **Fetch Bae et al. 2025** (Nucl Fusion 65 086051 / OSTI 2583820 / any arXiv) and do a point-by-point comparison (what they did vs us; agreements; what — if anything — we add, esp. the free-stream-vs-scattering dilution).
2. **Fix docs/deck to cite Bae honestly**: soften Tier-6 "steering doesn't cost breeding / TBR ~independent" → "modestly polarization-dependent (±0.6% here; consistent with Bae's ±2.7%)"; add Bae citation to RESULTS_tier6.md, the deck (Part 7 + literature slide), and reposition the whole project as validation/tooling (not novelty) in README/deck.
3. Then user emails Ethan; reposition results as independent validation + tool.

## 9. Next steps (post-Bae)
Angled/realistic field → couple Peterson's physics-informed source (sampler + anarrima g_*a kernels ready); depletion/multi-zone blanket; D-shape geometry; plasma/systems coupling for net-electricity & burn-efficiency (the upstream SPF benefits, out of neutron-transport scope).
