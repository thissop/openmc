# Overnight SPF Integration Master Plan

Last updated: 2026-09-01 22:20 EDT

## Objective

Deliver the smallest reviewable spin-polarized-fusion contribution on Ethan
Peterson's `stellarator_source` branch first, then extend the same verified
angular kernel to the spectral `StellaratorSource`. Preserve the gridded source
as an independent representation/physics oracle and turn the comparison into a
paper-quality subsection only if it includes end-to-end source or transport
consequences rather than timing alone.

## Authoritative workspace

- Worktree: `/Users/tkiker/Documents/GitHub/openmc-spf-pol`
- Branch: `spf-polarization`
- Base: `eepeterson/stellarator_source`
- Base commit: `cf616f4322864455ba46e46e4da04a87bd8d94b0`
- Original verified SPF implementation:
  `/Users/tkiker/Documents/GitHub/openmc/spf_prototype/native_spf/`
- Representation harness: `spf_validation/`
- Claude session: `638db22c-edea-421d-88fe-e7f5d5ecfd10`

Do not edit the original `spf-prototype` worktree as part of this integration.
Do not commit until all changes have been reviewed as one coherent patch.

## Single-writer ownership

| Window (America/New_York) | Owner | Rule |
|---|---|---|
| Now to 2026-09-01 00:30 | Codex | Claude session must remain unwritten. |
| 00:30 handoff | Codex controller | Stop at a safe boundary, inventory diff/tests, verify Claude is idle, send one uniquely tagged handoff using the full session ID, verify ACK. |
| After verified ACK to 05:30 | Claude | Claude owns edits in `openmc-spf-pol`; Codex does not write. |
| 05:30 onward | Codex controller | Read Claude messages since 00:30, verify Claude is idle, reconcile git/test state, then resume edits. |

If ownership cannot be verified, do not create a competing writer. Record the
blocker in the handoff ledger.

## Non-negotiable technical invariants

1. With no `<polarization>` element, source direction sampling and RNG draw
   order remain identical to Ethan's isotropic implementation.
2. All sampler state is constructed once and read-only inside `sample()`.
3. Pure A samples the `sin^2(theta)` shape; pure B and pure C share the
   `1/4 + 3/4 cos^2(theta)` shape. B and C differ only in total rate.
4. Non-polarized `(1/3, 1/3, 1/3)` is isotropic and has `eta = 2/3`.
5. Direction vectors and local magnetic-field directions are normalized.
6. VMEC/DESC coordinate conventions must be checked against an independent
   field oracle; do not accept a visually plausible field direction.
7. Preserve all unrelated and pre-existing user changes.

## Track A: reviewable Tokamak polarization PR

### A1. Python API

- [x] Add `spin_fractions_to_abc`.
- [x] Add `polarization`, `field_model`, `safety_factor`, and `field_sign`.
- [x] Add XML serialization and parsing.
- [x] Validate spin-fraction dictionary keys, component bounds, and the two
      species' normalization explicitly.
- [x] Validate safety-factor grids: finite, strictly increasing `r/a`, finite
      nonzero/positive `q` according to the adopted convention, and suitable
      radial coverage.
- [x] Add focused API and round-trip tests to Ethan's existing
      `tests/unit_tests/test_source_tokamak.py`.

### A2. C++ implementation

- [x] Add shared inverse-CDF functions for the perpendicular and parallel P2
      shapes. Source: `native_spf/source.cpp`.
- [x] Add polarization parsing and derived `p_perp` state.
- [x] Add q-profile interpolation and Miller-geometry pitched `B-hat`.
- [x] Add local-to-global polarized-direction sampling.
- [x] Gate the new direction path so the unpolarized branch is unchanged.
- [x] Prefer one shared internal SPF helper for Tokamak and Stellarator rather
      than duplicating the kernel between the two classes.

### A3. Gates

- [x] Python syntax/import gate.
- [x] Existing Tokamak unit tests.
- [x] New XML/API negative and round-trip tests.
- [x] Compiled source-sampling moments: `<cos^2 theta> = 1/5, 1/3, 7/15`
      for A, unpolarized, and B/C within statistical tolerance.
- [x] Confirm unpolarized XML omits all SPF elements.
- [x] Format changed C++ with the repository `.clang-format`; final Python
      style check remains part of packaging.

Tokamak completion condition: a compact, backward-compatible patch with actual
compiled sampling evidence, ready for review independently of stellarator work.

## Track B: VMEC spectral Stellarator polarization

### B1. Production data model

- [x] Add the polarization API to Ethan's `StellaratorSource`, reusing Track A
      validation instead of maintaining a second interpretation.
- [x] Extend VMEC reading to carry rotational transform and lambda coefficients
      with explicit radial-grid and Fourier-mode conventions (half-mesh handled;
      `from_vmec` polarization path + `test_..._native_half_mesh`).
- [x] Serialize/parse the additional spectral field data.
- [x] Retain accepted `(rho, theta, zeta)` coordinates through the C++ rejection
      loop so the field is evaluated at the sampled birth point.

### B2. Analytic field direction

- [x] Derive the physical field-line tangent from geometry derivatives, iota,
      and lambda derivatives in the exact VMEC convention
      (`B^theta = iota - lambda_zeta`, `B^zeta = 1 + lambda_theta`, coordinate
      basis form to avoid dividing by a small `1 + lambda_theta`).
- [x] Implement the derivation in a small independently testable helper
      (`StellaratorSource::field_direction`).
- [x] Handle the magnetic axis and near-degenerate tangent norms explicitly.
- [x] Rotate the cylindrical tangent into Cartesian coordinates and normalize.
- [x] Reuse Track A's verified P2 sampler around this `B-hat`.

### B3. VMEC gates

- [x] Finite-iota circular-torus limit (stellarator polarized-direction test
      rebuilds B-hat independently in numpy and recovers `<cos^2>` = 1/5, 1/3,
      7/15).
- [x] Pointwise angular error against **independently exported** VMEC
      contravariant components `bsupumnc/bsupvmnc` (not the lambda formula) on
      QA and QH: p99 = 1.3e-4 deg (QA), 3.8e-4 deg (QH), 0.022 deg (lower-res
      QUASR). See `vmec_bhat_validation.json` / `validate_vmec_bhat.py`.
- [x] Compiled source-direction moment checks using reconstructed local B-hat
      (includes a genuinely 3-D lambda harmonic so a toroidal-only bug fails).
- [x] Confirm unpolarized path leaves birth distribution and RNG draw order
      unchanged (gated `polarized_`; `..._unpolarized_xml_unchanged`).

VMEC completion condition MET: spectral field reconstruction agrees with VMEC's
own exported field to <1e-3 deg at source-relevant points on QA and QH, and the
unpolarized branch is byte-identical to Ethan's isotropic path.

## Track C: DESC and publication-grade representation study

### C1. DESC integration -- DONE and VALIDATED

- [x] **DESC-native iota/lambda extraction implemented in `from_desc`.** When
      polarization is requested (and field data not supplied), the DESC stream
      function `lambda` (`L_lmn` on `L_basis`) is Zernike-collapsed onto the rho
      grid and Ptolemy-folded into the SAME combined VMEC-style mode table as the
      geometry (the mode union now includes L modes), and `iota` is evaluated on
      rho (live eq via DESC `compute`; HDF5 via the stored PowerSeriesProfile).
      Helpers: `_desc_lambda_data`, `_desc_iota`.
- [x] **Validated against DESC's own `compute('B')`** (the independent field
      oracle, analogous to VMEC's `bsupumnc/bsupvmnc`), via
      `spf_validation/validate_desc_bhat.py`. Field-direction error is sub-degree
      and CONVERGES with n_rho (precise_QA median 0.031/0.014/0.010 deg at
      n_rho=16/33/49; precise_QH 0.048/0.023 deg) -- a sign/convention bug would
      stay near 90 deg, so the convention is confirmed. Integrated end-to-end
      through `from_desc(..., polarization=...)`: median 0.020 deg, p99 0.21 deg.
      Record: `desc_bhat_validation.json`.
- [x] **End-to-end COMPILED-sampling proof** (`validate_desc_sampling.py`): a
      DESC-built polarized source sampled through the real OpenMC engine
      (`from_desc -> XML -> C++ constructor -> C++ field_direction -> C++
      sampler`) yields the correct angular moments about the DESC field. Births
      are inverted back to flux coordinates by a vectorized Newton on the
      geometry (residual ~1e-13 cm), then `<cos^2>` is measured about the
      reconstructed B-hat: precise_QA perp/par/unpol = 0.1998/0.4680/0.3330;
      precise_QH = 0.1994/0.4651/0.3327 (analytic 0.2000/0.4667/0.3333). This
      closes the full chain -- DESC is not just extracted correctly, it samples
      correctly in the compiled source.
- [x] Docs (`settings.rst`) and the unit suite updated; the old
      NotImplementedError guard/test replaced (DESC now works). The heavy DESC
      validation runs in a separate `descenv312` interpreter (desc 0.17.3 needs
      py<=3.12; py3.13 hits a jaxtyping bug) and is not part of the core unit
      suite, which stays desc-free (56 tests green).

### C2. Paper subsection: maximum-value version

Working title: **Source representation tradeoffs for Monte Carlo fusion
neutronics**.

The subsection earns space in the paper if it answers this engineering question:

> When does the faster gridded representation change physically meaningful
> source and transport results enough that the compact spectral representation
> is necessary?

Required evidence -- STATUS 2026-09-02:

1. [x] 20-device QA/QH sweep, median/IQR (`agg_representations_20dev.json`).
2. [x] Accuracy vs resolution for R, Z, sqrt(g), field-direction angle; plus a
      high-res extension to 48³/243 MiB showing sqrt(g) does NOT refine away
      (`agg_representations_highres.json`).
3. [x] Speed / memory / serialized size (grid 2.6-3.1x faster/eval but 12-325x
      more memory; spectral 0.09 MiB exact).
4. [x] Source-distribution impact: source-weighted B-hat angular error + birth
      density; polarized moments verified.
5. [x] **End-to-end tally impact (the headline).** Free-streaming NWL wall-current
      bias vs spectral truth, with an N-scaling noise-vs-bias discriminator:
      coarse 8³ grid carries a coherent under-steering bias growing ~sqrt(N)
      (5.3 sigma at 20M), physically ~1% and NOT refining away; and across the
      full 20-device ensemble the coarse-grid bias correlates with device
      steering strength (**Pearson r = +0.83, N=20**) -- the grid error bites in
      the high-steering regime that motivates SPF. See
      `spf_validation/end_to_end/RESULTS_end_to_end.md`, figures
      `e2e_bias_vs_resolution.png`, `e2e_ensemble_correlation.png`. OpenMC
      `FileSource` cross-check job 9607486 completed all three 3-million-history
      arms successfully. It reproduced the coarse grid's directional signature
      (integrated current grid - spectral: -0.82% inboard, +0.79% outboard), while
      the fine grid was -0.004%/+0.22%. This validates the near-void projection
      with OpenMC's engine, but is not a scattering-blanket result.
6. [x] Decision map: spectral by default; the grid's error floor (~1% wall bias,
      near-axis sqrt(g)) is irreducible by uniform refinement, so "grid above
      resolution Y" is the wrong frame -- see conclusion in RESULTS_end_to_end.md.

### C3. Hybrid method (the stretch, now PROTOTYPED)

- [x] **Spectral-truth-driven adaptive radial grid** (`adaptive_grid.py`): error
      equidistribution against the spectral sqrt(g) concentrates radial nodes near
      the axis. At equal memory it is 3.2x more accurate in sqrt(g) and **6x** in
      the source-weighted B-hat (0.92 -> 0.15 deg at 30 MiB), i.e. equal accuracy
      at half the memory (`adaptive_frontier.png`).
- [x] **End-to-end closure** (`adaptive_e2e.py`): at identical 30 MiB the adaptive
      grid cuts the coherent wall bias from uniform's |z_mean| 0.32/0.45 (inb/out)
      to **0.04/0.17** -- ~8x, essentially to the noise floor. The hybrid grid
      reproduces the spectral wall pattern at the same memory. The completed
      48x96x288 uniform ceiling used 60.8 MiB yet retained |z_mean| 0.24/0.43
      (max |z| 3.2), so doubling radial storage did not match the adaptive grid.
- [ ] 3-D monitor (theta,zeta adaptivity) and multi-device frontier: next step.

High-value stretch result: use the spectral representation as the source of
truth to construct an adaptive or compressed grid whose refinement is driven by
Jacobian and `B-hat` error. That would turn the comparison from a benchmark into
a new hybrid method with a measurable accuracy/speed frontier.

Avoid a subsection based only on lookup timing and pointwise geometry error.
That is supporting-method material, not a strong scientific result.

## Verification commands and environment notes

Record commands and exact outcomes here as they are run. Do not hide dependency
or build failures.

- A task-local Python dependency prefix under ignored `build/` supplies `endf`
  and editable package metadata without changing the user environment.
- Configure/build: succeeded with `OPENMC_BUILD_TESTS=ON`, Debug, OpenMP off.
- `tests/unit_tests/test_source_tokamak.py`: **45 passed**. A dummy CE library
  manifest under ignored `build/` was sufficient because these tests sample a
  void source and never load nuclides.
- New compiled SPF gates: **5 passed** for A, explicit nonpolarized, B, C, and
  q-pitched A-mode sampling.
- C++ tests excluding optional `test_mcpl_stat_sum`: **11/11 passed**. The full
  CTest invocation reports only `test_mcpl_stat_sum`, whose own output says it
  was skipped because the optional MCPL library is unavailable.
- Final focused Python gate after reconciliation and packaging: **62/62
  passed** (45 Tokamak + 11 Stellarator + 6 representation-harness tests).
- `clang-format 22.1.3 --dry-run --Werror` is clean on the changed C++ files;
  CI's required version 18 is not installed locally, so its exact formatting
  check remains a PR/CI gate. `git diff --check` is clean.
- User-facing RST now documents SPF inputs, normalization/strength semantics,
  tokamak field models, VMEC extraction, and the guarded DESC limitation.
  Sphinx was not available in the local runtime, so a rendered docs build was
  not run.
- The 22-device representation campaign is separate evidence and does not block
  the Tokamak PR.

### 20-device grid-vs-spectral sweep (2026-09-01, `agg_representations_20dev.json`)

Ginsburg job 9585942 produced **20/22** converged QUASR equilibria (nfp 1-7;
failed: 803097 [known unsteerable perp-optimal outlier] and 948374). All 20
were pulled and converted to ignored `wout_arrays_*.npz` bundles under
`spf_validation/wouts/quasr20/`. Sweep ran on all 20, **0 failures**:

| grid (rho,theta,zeta) | grid speed | grid mem | spectral mem | sqrt(g) RMS (med, IQR) | B-hat p99 err (med, IQR) |
|---|---|---|---|---|---|
| 8x32x96   | 3.1x | 1.1 MiB  | 0.093 MiB | 0.211 [0.196,0.218] | 5.36 deg [1.28,11.46] |
| 12x48x144 | 3.0x | 3.8 MiB  | 0.093 MiB | 0.168 [0.157,0.174] | 0.98 deg [0.50,1.37] |
| 16x64x192 | 2.8x | 9.0 MiB  | 0.093 MiB | 0.145 [0.135,0.149] | 0.62 deg [0.40,0.83] |
| 24x96x288 | 2.6x | 30.4 MiB | 0.093 MiB | 0.119 [0.110,0.122] | 0.44 deg [0.28,0.58] |

Honest headline: the gridded representation is ~2.6-3.1x faster **per evaluation**
but costs 12-325x more memory and never becomes exact. The Jacobian sqrt(g)
(which sets the birth-density) keeps a ~12-21% median RMS error even at the
finest grid, tight IQR => systematic and axis-dominated, not noise. The
polarization-relevant metric -- the source-weighted field-axis angular error --
is 5.4 deg (IQR up to 11 deg) on a coarse grid and only drops below ~0.5 deg at
the 30 MiB grid. So a coarse grid materially mis-orients the SPF steering axis.
This is a real tradeoff tied to the paper's own physics, not just lookup timing.
**Still missing for a top-tier subsection: end-to-end tally impact** (does the
12% sqrt(g) / few-deg B-hat error move a wall-load or coil-flux tally?) -- needs
a transport run and is the C2 stretch item.

## Handoff ledger

| Time | From -> To | Delivery/session | Git state | Verified results | Next bounded step |
|---|---|---|---|---|---|
| 2026-08-31 22:00 EDT | User/Claude -> Codex | Claude `638db22c-edea-421d-88fe-e7f5d5ecfd10`; no new delivery yet | `openmc/source.py` modified; `spf_validation/` untracked; HEAD `cf616f432` | Python file previously parsed; focused pytest now blocked during import by missing `endf` | Finish Track A C++ graft, integrate tests, establish a usable test environment |
| 2026-08-31 22:10 EDT | Codex | No Claude delivery; Codex owns worktree | Modified `include/openmc/source.h`, `openmc/source.py`, `src/source.cpp`, and Tokamak tests; plan untracked | Build succeeded; 45/45 Tokamak tests and 11/11 non-MCPL C++ tests passed; compiled SPF moments passed | Begin Track B VMEC iota/lambda production data model and independent field convention gate |
| 2026-09-01 ~02:00 EDT | Codex -> Claude (user handed off after Codex hit usage limit) | Claude resumes in same worktree; Codex idle (limit until 2026-09-07) | Same 4 tracked files + untracked `spf_validation/`, `test_source_stellarator.py`, plan | Claude re-verified from scratch: **55/55** Tokamak+Stellarator unit tests pass against the built dylib (incl. compiled SPF direction-moment gates `<cos^2>` = 1/5, 1/3, 7/15 and the pitched-field and 3-D-lambda end-to-end checks). Inverse-CDF math re-derived and confirmed. | Guard DESC polarization; run 20-device sweep; verify Ginsburg generation |
| 2026-09-01 (Claude session) | Claude | Owns worktree per plan (05:30 window not yet reached; Codex idle) | Added `from_desc` NotImplementedError guard + test (56/56 now); pulled and converted 20 QUASR equilibria to ignored array bundles under `spf_validation/wouts/quasr20/`; launched 20-device sweep | Ginsburg genwouts job 9585942 COMPLETED: **20/22 devices converged** (nfp 1-7; 2 failed: 803097 known unsteerable outlier + 948374). `from_desc` guard test green. | Land 20-device sweep result into paper evidence; final packaging review |
| 2026-09-01 20:55 EDT | Claude -> Codex | Exact Claude session `638db22c-edea-421d-88fe-e7f5d5ecfd10` verified idle; no competing writer; Codex reclaimed ownership | HEAD remains `cf616f432`; 4 tracked source/test files modified plus untracked plan, validation harness, and Stellarator tests; generated wouts, logs, and caches are ignored | Reconciled Claude's DESC guard and 20-device aggregate against files; 20/20 data bundles present; aggregate says 20 devices/0 failures; rebuild succeeded; **56/56** focused Python tests and **11/11** non-MCPL C++ tests pass; `git diff --check` clean | Finish reviewable packaging, paper-subsection artifact, and meeting-deck update; defer native DESC conversion and transport-tally stretch unless independently validated |
| 2026-09-01 22:20 EDT | Codex final checkpoint | Exact Claude session remains idle; Codex is sole writer | Added two RST documentation edits plus `spf_validation/PAPER_SUBSECTION_DRAFT.md` and refreshed harness README; no commit; source HEAD remains `cf616f432` | **62/62** focused Python tests pass; build succeeds; **11/11** non-MCPL C++ tests pass; clang-format dry-run and `git diff --check` clean. Seven-slide meeting deck passes overflow and template-fidelity checks. | Human review and commit split: source/API/docs/tests in the OpenMC PR; keep `spf_validation/` and the paper study separate. Highest-value new science remains identical-history QA+QH transport tallies. |
| 2026-09-01 ~23:40 EDT | Codex local staging | Ginsburg bounded helper still blocked: old transfer host `128.59.124.72` times out and no reusable `~/.ssh/cm-tjk2147@ginsburg.rcs.columbia.edu:22` socket exists; no remote mutations or jobs submitted | Added untracked `spf_validation/transport_compare/` driver, analyzer, eight-case Slurm array, README, and standard-library-facing tests; all earlier dirty work preserved; no commit | New transport tests **6/6** pass; combined representation/transport helpers **12/12** pass; `py_compile`, `bash -n`, and `git diff --check` pass. Local VMEC bounds show the canonical QA wout is already reactor-scaled, so QA spectral length scale is correctly 1.0 (not 10.10). | Once a human-authenticated `ssh ginsburg` master is left connected: inspect exact QA/QH files and builds, verify source/geometry bounds and material tags, stage package, compile the spectral tree separately, then submit a 20k-history smoke array before any production histories. |

## Final packaging

- [x] Review full diff for accidental generated files and `__pycache__`.
- [x] Keep `spf_validation/` in a separate research branch/package; generated
      wouts, logs, caches, NetCDF, and NPZ files are ignored and must not enter
      the OpenMC source PR.
- [x] Update user-facing docs with inputs, defaults, normalization semantics,
      field models, and VMEC/DESC limitations.
- [x] Draft the paper subsection at
      `spf_validation/PAPER_SUBSECTION_DRAFT.md` with the 20-device table,
      interpretation, and explicit transport-claim boundary.
- [x] Produce the meeting deck at
      `/Users/tkiker/Documents/GitHub/openmc/spf_prototype/talk/peterson_checkin_2026-09-01_updated.pptx`;
      original deck preserved, overflow/fidelity checks pass.
- [ ] Produce one clean commit only after review; do not push or open a PR
      without explicit user authorization.
