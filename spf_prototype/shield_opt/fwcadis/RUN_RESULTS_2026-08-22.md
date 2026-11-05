# FW-CADIS path results — 2026-08-22

All cluster work used the source-built OpenMC `0.15.4-dev212` executable at
commit `23572ccc928a9ab7be19e6d7d2868852f12c0f0e` on Ginsburg. Jobs ran through
Slurm in the AdvancedAgents scratch root; no work ran on the login node.

## Core isolation

Job `9467343` completed the canonical three-region CSG test. Both the forward
and adjoint random-ray solves ran and exported a valid 343-bin
`weight_windows.h5` with every upper bound positive. This rules out the OpenMC
build and core FW-CADIS implementation as the source of the archived failure.
The apparent continuous-energy fallback in earlier logs was the expected
per-material infinite-medium MGXS generation.

## Path A — native one-period DAGMC FW-CADIS

Job `9467364` completed the real one-period DAGMC forward and adjoint
random-ray solves and exported one 10 x 10 x 10 weight-window object:

- Active bins: 595 / 1000.
- Weight-window file: 28,184 bytes.
- SHA-256: `61b61b88b7c57814eb0cb2dc5ad5e374996f4b63630ae965e48d03c802878021`.
- Cached CASMO-2 MGXS SHA-256:
  `728e29b623e2674a3bef2bd4d17f1c07045dc7c73e999d156802b3e67348924e`.
- Local artifacts: `dagmc_oneperiod/cluster_pull/run_9467364/`.

The controlled runner fixes the archived prototype's integration defects:

1. The local target is passed as an `openmc.Tallies` collection instead of an
   invalid list of `Tally` objects.
2. The requested ray count is restored after `convert_to_random_ray()` rather
   than silently accepting OpenMC's geometric default.
3. The MGXS library is generated once and reused. The negligible-density
   `gap` and `Vacuum` placeholders receive explicit near-zero two-group data,
   avoiding pathological maximum-event histories in the 1000-m reflective
   infinite-medium generator.
4. Runtime macro names retain the ID suffixes introduced by DAGMC
   synchronization.
5. The native `StellaratorSource` fluxmap is represented by 2,000 weighted,
   first-period point sources because random ray accepts only isotropic
   `IndependentSource` point or domain sources.

This is an integration pass, not a production window map. The forward average
source-region miss rate was 19.38%; the adjoint solve reported 43.39%. Increase
ray density and active distance before using these windows for final scoring.

## Path B — continuous-energy reciprocity-CADIS

Job `9467350` completed 400,000 continuous-energy histories from 320
equal-arc-length samples on the first five QH coil filaments. It constructed a
16 x 16 x 8 reciprocal importance map and manual CADIS windows using
volume-normalized inverse importance:

- Positive sampled bins: 2000 / 2048.
- Reliable bins at relative error <= 0.5: 1992 / 2048 (97.27%).
- Median reliable relative error: 0.0350.
- Maximum accepted relative error: 0.4686.
- Weight-window file: 44,952 bytes.
- SHA-256: `168be831a381bfbe625d4327b7f74008d76508b9816cb06aa949ce370c0f18ab`.
- Local artifacts: `reciprocity_cadis/cluster_pull/run_9467350/`.

This is a spatial reciprocity approximation rather than a fully
energy-transposed adjoint solve. Its advantage is that the importance map uses
the real continuous-energy DAGMC transport and avoids MGXS/random-ray
dependencies. Unsampled and high-error bins are disabled using OpenMC's
negative-bound convention.

## Production next step

Subdivide the 40-cm shield into roughly eight 5-cm ParaStell radial layers,
then rerun Path A with enough rays to eliminate the high miss rate. In
parallel, apply both candidate window files to identical bounded SSR replays
and compare coil-response figure of merit, bias, split behavior, and runtime
against analog and existing MAGIC windows.

## SSR application follow-up

The bounded comparison is complete. The present native map over-splits, the
reciprocity map is effectively neutral, and MAGIC remains statistically
unresolved. Exact inputs, jobs, metrics, and pulled artifacts are recorded in
`ssr_benchmark/RUN_RESULTS_2026-08-22.md`. The shield-subdivision and improved
native random-ray solve above remain the next calculation; none of the current
maps warrants a 24-million-history production replay.
