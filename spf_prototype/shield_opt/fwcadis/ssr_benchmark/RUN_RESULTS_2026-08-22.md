# One-period SSR variance-reduction benchmark — 2026-08-22

## Verdict

None of the three current weight-window files is ready for a production
per-coil calculation.

- The 10 x 10 x 10 native FW-CADIS integration map is computationally
  nonviable when applied to the surface-source replay. It creates a deep
  secondary cascade even after source-weight rescaling and a conservative
  split ceiling.
- The 16 x 16 x 8 reciprocity-CADIS map is safe but effectively neutral. At
  2,000,000 source histories it reproduced the same single coil contribution
  as analog transport with nearly the same runtime.
- The 12 x 12 x 12 MAGIC map is safe under the common guard and produced three
  nonzero coil bins instead of analog's one, but every nonzero bin still had
  91.5--100% relative error. Its median nonzero FOM was not better than the
  analog result at this resolution.

Do not spend a 24-million-history production allocation on these files. The
next native FW-CADIS map should be generated only after radial shield
subdivision and a higher-quality random-ray solve reduce the 19.38% forward
and 43.39% adjoint source-region miss rates reported for job `9467364`.

## Controlled inputs

All completed comparisons used source-built OpenMC `0.15.4-dev212`, commit
`23572ccc928a9ab7be19e6d7d2868852f12c0f0e`, with the same:

- one-field-period DAGMC geometry and periodic planes;
- deep-freeze bank `surf_op_master.h5`;
- B-10 black overrides on cells 2, 3, and 4;
- coil cells 12--16 and 0.1--25 MeV fast-flux response;
- seed `104729`, 32 OpenMP threads, and 60 atom-percent Li-6;
- application limits `max_split=4`, `max_lower_bound_ratio=5`, and
  `max_history_splits=100`.

The surface bank contains 8,291,969 sites with total starting weight
8,226,805.9125. Dividing by the 4,000,000 master histories gives the replay
normalization 2.0567014781 source weight per master birth. Absolute-bias
metrics use this factor; pattern metrics normalize the five-coil vector to
unit sum.

The candidate files retain their generated energy ranges:

| method | mesh | active bins | energy range | SHA-256 |
|---|---:|---:|---:|---|
| MAGIC | 12 x 12 x 12 | 251 / 1,728 | 0--25 MeV | `f772b4967da17812fb0b45e0570ef71f5b5f8a450122e74c4b109a732688397c` |
| native FW-CADIS | 10 x 10 x 10 | 595 / 1,000 | all energies | `61b61b88b7c57814eb0cb2dc5ad5e374996f4b63630ae965e48d03c802878021` |
| reciprocity-CADIS | 16 x 16 x 8 | 1,992 / 2,048 | 0.1--25 MeV | `168be831a381bfbe625d4327b7f74008d76508b9816cb06aa949ce370c0f18ab` |

## Guarded 200,000-history pilot

Jobs `9467410` and `9467414` established that the common dynamic-rescaling
guard prevents the immediate runaway seen in the unguarded array `9467405`.
The pilot was too small for a dose verdict:

| method | positive coil bins | OpenMC simulation time | tracks/source |
|---|---:|---:|---:|
| analog | 1 / 5 | 26.1 s | 1.00 |
| MAGIC | 1 / 5 | 30.3 s | 1.15 |
| native FW-CADIS | 0 / 5 | 105.3 s | 4.16 |
| reciprocity-CADIS | 1 / 5 | 28.5 s | 1.07 |

## Two-million-history screen

The completed cases used 20 batches of 100,000 histories. The first analog
attempt ran on a slow node and was replaced by job `9467424` on `g231` before
analysis. MAGIC and reciprocity were array tasks from job `9467416`.

| method | positive coil bins | simulation time | tracks/source | coil fast flux, per bank-source weight | relative error of positive bins |
|---|---:|---:|---:|---|---|
| analog | 1 / 5 | 292.8 s | 1.00 | `[0, 0, 0, 5.6101e-6, 0]` | `[1.000]` |
| MAGIC | 3 / 5 | 319.3 s | 1.16 | `[0, 0, 6.5786e-6, 6.1244e-6, 2.0778e-5]` | `[1.000, 0.915, 1.000]` |
| reciprocity-CADIS | 1 / 5 | 287.2 s | 1.06 | `[0, 0, 0, 5.6101e-6, 0]` | `[1.000]` |

The native task was stopped after nine minutes rather than allowed to consume
the one-hour allocation. Four 100,000-primary batches had already reported at
least 3,677,449 secondary tracks; individual secondary generations peaked at
183,119 tracks. That is a lower-bound amplification of 10.2 total tracks per
primary over the reported portion, with no completed result envelope.

The sparse scores cannot support a physics-bias conclusion: all methods are
far below the master coil response and zero of five bins agree with the master
within combined 2-sigma uncertainty. They do support the variance-reduction
verdict above. In particular, the extra MAGIC hits are not yet a convergence
result, and the identical analog/reciprocity score shows that the current
reciprocity map does not materially alter rare-event sampling.

## Artifacts

- Guarded pilot envelopes and analysis:
  `cluster_pull/pilot_guarded_20260822/`
- Two-million envelopes, analysis, and the cancelled native log:
  `cluster_pull/screen_2m_20260822/`
- Reusable runners: `benchmark_ssr.py`, `run_ssr_benchmark.sbatch`,
  `analyze_benchmark.py`, and `analyze_ssr_benchmark.sbatch`.
- The launcher now preloads the source-built shared library only into the
  external OpenMC process. Job `9467423` verified that this removes the former
  duplicate-library teardown abort.

## Next calculation

Subdivide the 40 cm shield into approximately eight 5 cm radial ParaStell
layers. Regenerate native FW-CADIS windows with enough forward and adjoint rays
to make the source-region miss rate acceptably small, then repeat this exact
200,000-history split-safety pilot before scaling. A regenerated map must show
bounded track amplification and nonzero responses in more coil bins than
analog before receiving a multi-million-history allocation.
