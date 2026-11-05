# SPF source validation: Thaddaeus' `StellaratorSource` vs. Ethan's

Origination/validation case for the spin-polarization work. Ethan's `stellarator_source`
branch already ships a full `StellaratorSource` (VMEC/DESC **Fourier flux-surface**
representation, `from_vmec`/`from_desc`) and `TokamakSource`, with **no polarization**. The
plan is to (1) validate the two stellarator sources against each other on the QUASR QA/QH set,
adopt his as the base, then (2) add the spin-polarization layer + a B-hat backend onto his
`TokamakSource`/`StellaratorSource`. This directory is step 1.

Two representations answer the same question — given `(s, theta, zeta)`, where is the point and
what is the birth weight `sqrt(g)`:

| | representation | eval cost | memory | error |
|---|---|---|---|---|
| **Thaddaeus (grid)** | precompute `R,Z,sqrt(g)` on a `(nrho,ntheta,nzeta)` grid, trilinear interp | O(1) per draw | the grid | interpolation |
| **Ethan (spectral)** | keep VMEC `rmnc,zmns,...` and Fourier-sum on the fly | O(modes) per draw | the coeffs | exact (to VMEC truncation) |

## Phase 1 — representation head-to-head (runs now, no compiled OpenMC)
`fluxrep.py` implements both as pure NumPy; `metrics.py` scores accuracy/speed/memory;
`compare_representations.py` sweeps grid resolution on a real `wout`.

```bash
python compare_representations.py /path/to/wout_qh.nc --out results_qh.json
pytest -q test_fluxrep.py     # correctness gates (sqrt(g)->volume, node-exactness, convergence)
```

The current aggregate covers **20 converged QUASR devices**. Across the grid
resolution sweep, the grid is 2.6--3.1x faster per batched evaluation but uses
12--327x more memory than the spectral representation. Median normalized-source
total variation falls from 6.16% to 1.94%, while the median source-weighted
`B-hat` p99 error falls from 5.36 to 0.44 degrees. See
`agg_representations_20dev.json` for the complete per-device record and
`PAPER_SUBSECTION_DRAFT.md` for the paper-ready interpretation and caveats.

## Phase 2 — end-to-end transport comparison (next)
Run identical-history grid and spectral sources on representative QA and QH
models and compare wall loading plus coil flux or heating with Monte Carlo
uncertainties. This is the evidence needed to decide whether the measured
representation errors are engineering-significant; the current results should
not be described as a transport-tally effect.

## Data
- Smoke: `spf_prototype/talk/geom/wout_qh.nc` (QH) + `spf_prototype/data/*_fluxmap.*` (QA).
- Canonical boundaries: `spf_prototype/shield_opt/data/raw/equilibria/input.LandremanPaul2021_{QA,QH}`.
- Full set: QUASR devices — local `spf_prototype/data/quasr/` (catalogue + serials) and Ginsburg
  `spf_work/spf_pp/data/quasr<ID>_wall_geom/`. Each device needs a `wout` (for `from_vmec`) and a
  fluxmap; generate `wout` from the QUASR boundary via the existing `quasr_fluxmap` pipeline.

## Note: iota + lambda are in the wout
`read_wout` also returns `iotaf` (rotational transform) and `lmns` (lambda stream function) — the
two free functions Ethan wants for computing B-hat on the fly. They are carried here so the B-hat
backend for the polarized stellarator source can use them directly (no separate file needed).
