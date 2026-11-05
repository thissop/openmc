# SPF source validation — final product, on Ginsburg (2026-11-05)

Validation of the spin-polarized fusion (SPF) neutron source as split for upstream
PRs, run on Columbia's Ginsburg cluster against compiled OpenMC. This supersedes the
prototype `CompiledSource` validation: it exercises the **merged** `TokamakSource` /
`StellaratorSource` `polarization=` API that the two stacked PRs introduce.

**PRs under test** (local branches in `../openmc-ethan-pr`, based on `openmc-dev/develop@e43e16583`):
- **PR A** `spf-tokamak` — shared SPF kernel + `TokamakSource` polarization.
- **PR B** `spf-stellarator` — `StellaratorSource` polarization (VMEC + DESC), stacked on A.

## Summary — all green

| Validation | Build | Result |
|---|---|---|
| Tokamak Schwartz / anarrima (full-stack MC) | `openmc_src` (pr3999) | **PASS** — perp(A) +40.0% / −21.2%, par(B) −40.6% / +21.0% vs oracle +40.6% / −21.2% |
| PR-B compiled samplers (unit, 73 tests) | `openmc_prb` (PR B) | **PASS** |
| PR-B DESC `from_desc` (unit, 10 tests) | `openmc_prb` (PR B) | **PASS** against real DESC equilibria |
| Bae STAR spherical tokamak (poloidal steering) | `openmc_prb` (PR B) | **PASS** — reproduces Bae 2025 Fig. 4 |

Figures: `figs/fig_schwartz_capstone_bw.png`, `figs/fig_bae_star_poloidal_bw.png`.

## 1. Tokamak Schwartz / anarrima capstone (full-stack Monte Carlo)

Square-cross-section torus, parabolic `(1−ρ²)` plasma, toroidal field — the anarrima
oracle geometry, scaled to cm, free-streaming (near-void). `TokamakSource(polarization=,
field_model="toroidal")`. Net current on the inboard (r-min) and outboard (r-max) walls,
poloidally resolved; midplane ratios cancel area, giving the direct % steering.

```
            inboard        outboard
unpol       0.011572       0.052855
perp (A)    0.016200       0.041656   ->  +40.0% / -21.2%
par  (B)    0.006880       0.063929   ->  -40.6% / +21.0%
```

Oracle (Schwartz 2025 §2, reproduced locally with anarrima 0.1.0): A **+40.6% inboard /
−21.2% outboard**; B/C the mirror. **Match < 1%.** The full-stack MC therefore confirms
the analytic oracle for the final merged source.

## 2. PR-B compiled samplers — unit tests

`pytest tests/unit_tests/test_source_stellarator.py test_source_tokamak.py` on the PR-B
build: **73 passed**. Includes the compiled, transport-exercised checks:
- `test_tokamak_source_polarized_direction_sampling` (all 4 polarizations) and
  `…_pitched_direction_sampling`.
- `test_stellarator_source_polarized_direction_sampling`: ⟨cos²θ⟩ = 1/5, 1/3, 7/15 for
  A / unpolarized / B about an independently rebuilt B̂ (from VMEC λ + ι).
- `test_stellarator_source_from_vmec_uses_native_half_mesh`, polarization round-trips,
  input validation.

## 3. PR-B DESC path — `from_desc` unit tests

`pytest … -k desc`: **10 passed** (`from_desc`, `desc_file`, `desc_invalid`, each for
axisymmetric + stellarator, both symmetries) against equilibria built with **desc 0.17.2**.
This validates PR-B's rebased DESC adapter — λ via `convert('L')` + mode-union, ι via
`eq.compute('iota')` — which replaced the removed manual `_desc_*` helpers.

## 4. Bae STAR — polarized steering in a spherical tokamak

Comparison to **Jin Whan Bae et al., *Nucl. Fusion* 65 (2025) 086051**, "Neutronics
analysis of spin-polarized fuel in spherical tokamaks." STAR device: R=4 m, a=2 m
(aspect ratio 2), κ=2.3, δ=0.65, dominant toroidal field. Bae's emission schemes map 1:1
onto ours: *perpendicular* (spin-aligned) `W∝sin²θ` = our A; *parallel* (anti-aligned)
`W∝¼+¾cos²θ` = our B. (Bae's (a,b,c) are deuteron/triton vector & tensor polarizations,
not Schwartz's collision-mode fractions — different parameterization, identical shapes.)

Method (mirrors Bae §2.1): sample the STAR `TokamakSource`, ray-trace each birth neutron
to a conformal D-shaped first wall, and histogram the first-hit poloidal angle (0° =
outboard midplane, 180° = inboard). Result:

```
          inboard (135-225 deg)   outboard (<45, >315 deg)
unpol     0.0728                  0.3132
perp      0.1023  (+40.6%)        0.2512  (-19.8%)
par       0.0425  (-41.6%)        0.3763  (+20.1%)
```

Perpendicular drives the neutron load **inboard**, parallel drives it **outboard** —
exactly Bae's Fig. 4 trend, and the magnitudes agree with the Schwartz oracle (same
sin²/cos² steering physics, independent geometry).

## Environment / reproduce (Ginsburg)

Connect `ssh -o BatchMode=yes ginsburg` (Duo on an interactive `ssh ginsburg` first).
Login node forbids >few-second compute → everything via `sbatch` (`account=astro`,
`partition=short`). Common env:

```
export PYTHONNOUSERSITE=1                       # ~/.local shadows pandas otherwise
export PYTHONPATH=$SRC
export LD_PRELOAD=$SRC/build/lib/libopenmc.so
export PATH=$SRC/build/bin:$ENV/bin:$PATH        # ENV=$SPF/envs/spf-stellarator
export OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml   # ENDF/B-VIII.0
```

- `SPF=/ginsburg/astro/users/tjk2147/spf_work`
- Tokamak capstone build `SRC=$SPF/openmc_src` (pr3999, has polarized TokamakSource).
- PR-B build `SRC=$SPF/openmc_prb` — clean clone of `develop@e43e16583` + `$SPF/prb_combined.patch`
  (`git diff develop..spf-stellarator`), configured DAGMC OFF / MPI ON.
- DESC tests: `desc` lives only in `~/desc_venv`; its HDF5 clashes with the LD_PRELOAD'd
  `libopenmc.so`, so `from_desc` (pure Python) runs **without** LD_PRELOAD.

**Caveats.** Free-streaming (no scattering) isolates source-direction physics, matching
Schwartz's and Bae's normalized presentations. The Bae comparison uses a representative
`(1−ρ²)²` emission profile and `field_model="toroidal"` (not STAR's full equilibrium
n,T profiles or pitched field); it reproduces the steering *trend* and magnitude, not an
absolute point-match of Bae's TBR/fluence. The exit-time `__cxa_finalize` double-free in
the logs is benign (fires after statepoints/JSON are written).
