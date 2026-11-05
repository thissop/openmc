# QA/QH grid-versus-spectral transport gate

This directory turns the representation-level comparison into a controlled
transport experiment.  It runs QA and QH with the gridded and VMEC-spectral
`StellaratorSource` implementations while holding geometry, materials, source
emission, source energy, tallies, history count, and RNG seed fixed.

The default Slurm array includes unpolarized emission (isolates spatial-source
differences) and mode A polarized emission (also exercises the field-axis
representation), for eight cases total.  Smoke defaults are 20,000 histories
per case.  Increase `COMPARE_BATCHES` and `COMPARE_PARTICLES` only after every
smoke case writes `result.npz` without lost-particle or material-tag failures.

The shared seed makes each run reproducible, but the two samplers consume random
numbers differently.  It therefore does not create a particle-by-particle
common-random-number coupling.  `analyze_transport_compare.py` conservatively
combines the two tally standard deviations in quadrature.

The canonical QA VMEC file staged with this package is already reactor-scaled:
its sampled plasma spans roughly R=6.2--13.2 m before OpenMC's m-to-cm
conversion.  `QA_LENGTH_SCALE` therefore defaults to 1.0.  If a native-scale
Landreman--Paul wout is substituted, set the verified reactor scale explicitly;
only spectral R/Z coefficients are rescaled, while iota and lambda remain
dimensionless.  Always compare the staged wout bounds with the grid metadata
and QA DAGMC bounding box before production.  QH also defaults to unity.

Cluster launch prerequisites:

- an existing gridded-source build at `GRID_SRC`;
- a separate compiled spectral-source tree at `SPECTRAL_SRC`;
- the QA/QH DAGMC, cell-discovery, fluxmap, and VMEC files named in the Slurm
  script (override any path through its environment variable);
- a successful smoke array before production histories are submitted.

After results land, run:

```bash
python analyze_transport_compare.py \
  results/qa_unpol_grid/result.npz,results/qa_unpol_spectral/result.npz \
  results/qa_A_grid/result.npz,results/qa_A_spectral/result.npz \
  results/qh_unpol_grid/result.npz,results/qh_unpol_spectral/result.npz \
  results/qh_A_grid/result.npz,results/qh_A_spectral/result.npz \
  --json comparison.json --markdown comparison.md
```

The smoke comparison is a pipeline gate, not a paper result.  A publication
claim also needs adequate tally precision and stability under an independent
seed.
