# SPF / OpenMC — direction from Ethan (2026-10-14 meeting) + plan

## What Ethan asked for

1. **Tokamak SPF PR first, into openmc proper (upstream `develop`).** Pull the
   spin-polarized tokamak source out of the combined PR and land it as its own
   standalone PR against openmc-dev, not against the stellarator branch.
2. **Rebase the stellarator SPF work onto Ethan's StellaratorSource PR** (his
   StellaratorSource is itself a PR into openmc). The stellarator polarization
   stacks on that, not on a one-off fork branch.
3. **Then the stellarator spin-polarized portion.** Must work for **both VMEC and
   DESC**, with the field direction built from the **lambda stream function**
   (λ) + iota (what we already do — confirm this is the "streaming function" he
   meant).
4. **Validation / comparisons:**
   - Tokamak **analytic**: compare against **Schwartz** via **anarrima**
     (his open-source analytic NWL package).
   - **Spherical tokamak**: compare against **Bae et al. 2025** (polarized source
     in OpenMC; they published this route first — frame our contribution
     carefully).
   - "some other things" — TBD, confirm the full list with Ethan.
5. **Paper draft** — Thaddaeus writes it (no AI); these are supporting notes only.

## PR restructure

Current state: PR eepeterson/openmc#9 = combined tokamak+stellarator on his
`stellarator_source` branch. This gets split:

- **PR A — Tokamak SPF -> openmc-dev `develop`** (land first)
  - shared SPF kernel: `spin_fractions_to_abc`, the two P2 inverse-CDF samplers,
    `sample_spf_direction` (the source-agnostic C++ helper)
  - `TokamakSource` polarization (toroidal + q-pitched field)
  - `tests/unit_tests/test_source_tokamak.py` additions
  - tokamak user-guide + API docs
  - This PR carries the kernel, since tokamak lands first.

- **PR B — Stellarator SPF -> rebased on Ethan's StellaratorSource PR** (after A)
  - `StellaratorSource` polarization: `field_direction` from iota + λ
  - `from_vmec` (half/full mesh) and `from_desc` (Fourier-Zernike λ + iota) field
    extraction
  - `test_source_stellarator.py`, stellarator docs
  - Depends on PR A (kernel) + StellaratorSource. Once both are on develop it
    rebases cleanly.

- **PR #9**: close or repurpose as PR B once the split is done. (Decide with Ethan
  whether to keep the thread or open fresh.)

## Validation plan (what proves each piece)

### Tokamak, analytic (Schwartz / anarrima)
- Reproduce the Schwartz square-cross-section-torus NWL example (his Fig. 2):
  filamentary ring, scattering off, wall-current tally binned poloidally.
- Modes: unpolarized, A, B, C (+ a mixture). Compare normalized wall patterns.
- Scalar oracle from Schwartz text: A mode ~ +43% inboard / -22% outboard
  midplane; B/C the mirror (-43% / +22%); isotropic ~ 12% outboard-hotter from
  geometry alone.
- Primary check: `(Q_MC - Q_analytic)/sigma` across wall bins ~ N(0,1), coherent
  mean ~ 0. anarrima is the reference implementation; our own quadrature is the
  secondary cross-check.
- Note the §1.2 correction: B and C share the `1/4 + 3/4 cos^2` angular shape
  (C is NOT isotropic); isotropic == pure-C normalized pattern is a cheap, strong
  check.

### Spherical tokamak (Bae et al. 2025)
- Reproduce Bae's polarized-source spherical-tokamak case and show agreement.
- Framing: Bae published polarized-source-in-OpenMC first; our contribution is the
  in-tree source classes + the stellarator (VMEC/DESC) generalization + the
  representation study. State this explicitly, don't overclaim novelty on the
  tokamak side.

### Stellarator (our own work, supporting)
- Field reconstruction vs each code's own field: VMEC `bsupumnc/bsupvmnc`
  (~1e-4 deg), DESC `compute("B")` (sub-deg, converges with n_rho).
- Compiled angular moments `<cos^2> = 1/5, 1/3, 7/15`; end-to-end DESC->sampler on
  QA + QH.
- Representation study: analytic vs grid vs adaptive (speed/memory/accuracy;
  grid's √g error doesn't refine away; coarse-grid wall bias correlates with
  steering, r=0.83; adaptive reaches spectral accuracy at ~6x less memory).

## Paper notes (skeleton — for Thaddaeus to write)

Working frame: methods paper on an in-tree spin-polarized fusion source for
OpenMC, tokamak + stellarator, with analytic verification.

- **Intro / prior art:** Schwartz (analytic SPF NWL, anarrima); Bae (first
  polarized source in OpenMC, spherical tokamak); Kulsrud (depolarization
  context). Position: first in-tree, source-agnostic sampler spanning tokamak +
  stellarator with VMEC/DESC field reconstruction.
- **Physics:** Schwartz Eq. 1-2; (a,b,c) modes; `w_perp = 3/4 a`,
  `w_par = 2/3 b + 1/3 c`; rate factor `eta = a + 2/3 b + 1/3 c`; B/C share shape.
- **Method:** stratified P2 inverse-CDF sampling; local->global frame rotation;
  field models — tokamak toroidal/q-pitched; stellarator B̂ from iota + λ stream
  function (VMEC half/full mesh; DESC Fourier-Zernike + Ptolemy).
- **Verification:** unpolarized == isotropic (bit-identical); compiled moments;
  tokamak NWL vs Schwartz/anarrima; spherical-tokamak vs Bae; stellarator field
  vs VMEC/DESC native field.
- **Representation study (stellarator):** analytic vs interpolated vs adaptive;
  the speed/memory/accuracy frontier; end-to-end wall-load bias and its scaling
  with steering strength.
- **Numbers to cite:** in `/Users/tkiker/Documents/GitHub/SPF_POLARIZATION_REPORT.md`
  and `openmc-spf-pol/spf_validation/`.

## Open questions for Ethan
- Exact base for the tokamak PR: openmc-dev `develop`?
- The "some other things" — full validation list?
- Keep PR #9 open as the stellarator thread, or close and reopen after the split?
- Does the tokamak PR need the Schwartz NWL comparison in-tree (regression test),
  or is that paper-only?
