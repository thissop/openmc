# Spin-polarized DT neutron source prototype (OpenMC)

A self-contained prototype of a **polarization-aware fusion neutron source**
(the birth-direction sampler) for OpenMC, validated against the unpolarized
isotropic limit and against the analytic axisymmetric NWL results of
**Schwartz 2025** (arXiv:2507.11758). It stands alone — no private codebase.

Physics: emission ∝ `dσ/dΩ = (σ0/2π)[¾a sin²θ + (⅔b+⅓c)(¼+¾cos²θ)]` where θ is
measured from the local B̂. Mode fractions `(a,b,c)` (Schwartz Eq. 1); total-rate
factor `η = a+⅔b+⅓c`; birth-direction pdf depends only on
`w_perp=¾a`, `w_par=⅔b+⅓c`.

## Positioning (read first): this is validation + tooling, not new physics
Monte-Carlo neutronics of spin-polarized fuel was already done by **Bae et al.
2025** (*Nucl. Fusion* **65**, 086051; OpenMC, Pb–Li spherical tokamak: TBR
parallel +2.7%, +68% magnet lifetime). This prototype's contribution is
**independent cross-validation + a reusable, verified tool**: a compiled
`openmc::Source` anchored by an analytic↔MC recovery check against Schwartz, which
independently reproduces Bae's TBR-vs-polarization *trend* via a different code
path (compiled C++ source), breeder (FLiBe), and geometry (square torus). The
naming maps as: Bae **"perpendicular"** = our **A** (∝sin²θ, +50% rate); Bae
**"parallel"** = our **B/C** (∝1+3cos²θ); Bae **"unpolarized"** = iso. Full
point-by-point comparison and the honest contribution statement are in
**`PRIOR_WORK_BAE2025.md`**.

## Layout
```
spf_prototype/
  NOTES_orientation.md     §0 ground truth: OpenMC API + the VERIFIED env/build recipe (§0.9)
  src/
    spf_sampler.hpp        pure-math core: abc→weights, inverse-CDF angular samplers, Gram-Schmidt
    polarized_fusion_source.cpp  openmc::Source wrapping the header (the .so)
    standalone_driver.cpp  tier-1 driver: header + random_lcg.cpp, no libopenmc
    CMakeLists.txt         builds the .so (find_package(OpenMC))
  python/
    spf_mirror.py          Python mirror incl. bit-exact port of OpenMC's PCG prn
    abc_modes.py           (a,b,c) algebra + sympy identities
    spf_stats.py           chi2/KS helpers
    analytic_nwl.py        independent NWL: quad + anarrima + Eq.5; parabolic plasma
    build_and_run.py       build .so + OpenMC model + wall-current extraction
    verify_sampler.py      tier-1 diagnostics  -> RESULTS_tier1.md, figs/
    verify_nwl_analytic.py tier-2a diagnostics -> RESULTS_tier2a.md, figs/
    verify_nwl.py          tier-2b diagnostics -> RESULTS_tier2b.md, figs/
  tests/                   pytest gates (test_abc_modes, test_sampler_stats, test_analytic_nwl)
  RESULTS_tier{1,2a,2b}.md · figs/ · LIMITATIONS.md
```

## Environment (one-time; full recipe + rationale in NOTES_orientation.md §0.9)
```bash
sudo apt-get install -y libhdf5-dev
python3 -m venv "$HOME/spf_venv"
"$HOME/spf_venv/bin/pip" install -e <openmc-repo> sympy pytest
"$HOME/spf_venv/bin/pip" install "git+https://github.com/PrincetonUniversity/anarrima.git"
# build + install OpenMC into the venv prefix:
cd <openmc-repo> && mkdir -p build && cd build
cmake .. -DOPENMC_USE_MPI=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo \
         -DOPENMC_BUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX="$HOME/spf_venv"
make -j"$(nproc)" && make install
# cross sections (NNDC HDF5):
export OPENMC_CROSS_SECTIONS="$HOME/nndc_hdf5/cross_sections.xml"
```

## Reproduce every number
All scripts use fixed seeds. From the repo root:
```bash
PY="$HOME/spf_venv/bin/python"

# Tier 1 — direction sampler (no transport). 23 checks; builds the C++ driver.
$PY -m pytest spf_prototype/tests/test_abc_modes.py spf_prototype/tests/test_sampler_stats.py -q
$PY spf_prototype/python/verify_sampler.py          # -> RESULTS_tier1.md, figs/tier1_*

# Tier 2a — independent analytic NWL (quad vs anarrima vs Eq.5). 33 checks.
$PY -m pytest spf_prototype/tests/test_analytic_nwl.py -q
$PY spf_prototype/python/verify_nwl_analytic.py     # -> RESULTS_tier2a.md, figs/tier2a_*

# Tier 2b — OpenMC full-stack vs analytic (builds .so, runs 5 configs).
PATH="$HOME/spf_venv/bin:$PATH" OMP_NUM_THREADS=2 \
  $PY spf_prototype/python/verify_nwl.py            # -> RESULTS_tier2b.md, figs/tier2b_*
```

## Results in one line each
- **Tier 1**: C++↔Python sampler parity 6e-16; per-mode cosθ χ²/dof≈1; corners
  σ_tot/σ0 = 2/3,1,2/3,1/3; rotation isotropy + steering correct.
- **Tier 2a**: quad = anarrima = Eq.5 to 1.8e-15; parabolic-plasma quad-vs-anarrima
  5e-8; Fig.2 oracles reproduced (A inboard +40.6%, outboard −21.2%, iso out/in +14.6%).
- **Tier 2b**: OpenMC reproduces the analytic directionality within statistics on
  all walls/modes (means≈0); A inboard +39%, outboard −21%; **B≡C**, **iso≢C**,
  linearity all confirmed.

## Two physics notes carried through the work
1. **C is not isotropic.** B and C share the `¼+¾cos²θ` shape; only their total
   rate differs. The valid cheap equivalence is **B≡C** (the spec's "iso≡C" is the
   error §1.2 warned about). Verified analytically and in OpenMC.
2. The paper text says aspect ratio 2.5, but its figure-generating code uses 2.0;
   we match the reference implementation (anarrima). See `RESULTS_tier2a.md`.

## References
- **Schwartz**, "Analytic neutron wall loading for spin-polarized fusion,"
  arXiv:2507.11758 (2025) — the analytic free-streaming NWL we reproduce; his
  `anarríma` package is our analytic ground truth.
- **Bae, Borowiec, Badalassi, Parisi, Diallo, Menard, Khodak, Brown**,
  "Neutronics analysis of spin-polarized fuel in spherical tokamaks,"
  *Nuclear Fusion* **65**, 086051 (2025), DOI 10.1088/1741-4326/adf3c6 — the
  closest prior work (OpenMC MC neutronics of SPF); see `PRIOR_WORK_BAE2025.md`.
- **Kulsrud, Furth, Valeo, Goldhaber**, PRL **49**, 1248 (1982) — original
  spin-polarized DT rate/anisotropy.

OpenMC: `0.15.4-dev` @ commit `608a1c338`.
