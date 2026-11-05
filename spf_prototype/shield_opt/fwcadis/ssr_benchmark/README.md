# Surface-source variance-reduction benchmark

This directory compares four ways to score the five one-period coil responses
from the validated deep-freeze surface source:

1. analog transport;
2. the existing coarse MAGIC map;
3. the native random-ray FW-CADIS map from job `9467364`; and
4. the continuous-energy reciprocity-CADIS map from job `9467350`.

The array runner holds the geometry, source bank, black-inside cells, tally
definitions, source histories, random seed, OpenMC build, and split-safety caps
fixed. Each case receives a unique run directory and never overwrites the
canonical files under `pstl_test/corrected`.

At application time all three maps use `max_split=4`,
`max_lower_bound_ratio=5`, and a pilot `max_history_splits=100`. The common
dynamic-rescaling guard prevents a source-bank particle whose starting weight
sits far above a generated lower bound from immediately spawning an
uncontrolled secondary cascade.

The three candidate maps deliberately retain their generated energy ranges:
MAGIC covers 0--25 MeV, native FW-CADIS covers all neutron energies, and the
reciprocity map covers 0.1--25 MeV, matching the fast-flux response it targets.

## Interpretation

The default 4 x 50,000 array is only a wiring and split-amplification pilot.
`estimated_tracks_per_source_history` is the OpenMC track rate divided by its
source-particle rate. The conventional per-coil figure of merit is
`1 / (R^2 T)`, where `R` is relative standard error and `T` is OpenMC's
simulation time.

The analyzer multiplies each replay response by the surface bank's total
starting weight divided by the 4,000,000 master histories. This restores the
per-master-birth normalization. Absolute bias then includes the intended
black-inside approximation; the separately reported pattern error removes the
overall response magnitude.

A production conclusion requires all five coil bins to be positive, stable
results across independent seeds, no unacceptable track amplification, and
agreement with the master within combined uncertainty. The 10 x 10 x 10 native
map remains an integration map because its random-ray source-region miss rates
were 19.38% forward and 43.39% adjoint.

Run `analyze_ssr_benchmark.sbatch` with `BENCH_LABEL` set to a completed result
set. The source-bank scan is submitted through Slurm rather than performed on
the login node.

The source-built shared library is preloaded only into the external `openmc`
process. Python loads the matching library shipped beside the source-tree
bindings, avoiding two independent library finalizers in one interpreter.

See `RUN_RESULTS_2026-08-22.md` for the guarded pilot and two-million-history
screen. The current native, reciprocity, and MAGIC files all failed the stated
production-readiness test for different reasons.
