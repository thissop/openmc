# Source-representation tradeoffs for polarized stellarator neutronics

## Draft subsection

The stellarator source admits two useful numerical representations. A spectral
source retains the VMEC Fourier coefficients and evaluates the geometry,
Jacobian, and magnetic-field direction analytically at every sampled point. It
is exact up to VMEC's own radial interpolation and spectral truncation. A
gridded source tabulates the same quantities and uses trilinear interpolation,
exchanging interpolation error and storage for constant-cost evaluation. This
tradeoff matters more for spin-polarized fuel than it does for an isotropic
source: the Jacobian controls where particles are born, while the interpolated
field direction controls the axis of the anisotropic birth-direction law.

We compared the representations on 20 converged QUASR equilibria spanning one
to seven field periods. Two attempted equilibria did not converge and were
excluded before the representation comparison (QUASR 803097, a known
unsteerable outlier, and 948374). For each retained device and resolution, the
grid was constructed directly from the spectral representation. Accuracy was
then evaluated at 20,000 common off-grid points. We report medians and
interquartile ranges across devices because the near-axis Jacobian error is
skewed. The source-distribution metric is total-variation distance after
normalizing the two Jacobian-weighted discrete distributions. The field metric
is the 99th percentile of the unoriented-axis angular error, weighted by the
spectral source density; identifying `B-hat` and `-B-hat` is appropriate because
the implemented P2 angular laws are invariant under field reversal. Throughput
is the best of three vectorized NumPy trials on the local development machine,
so the absolute rates are implementation- and hardware-specific.

| Grid `(rho, theta, zeta)` | Evaluation speed vs spectral | Stored grid | Source TV, median [IQR] | Source-weighted field-axis p99, median [IQR] |
|---|---:|---:|---:|---:|
| 8 x 32 x 96 | 3.14x | 1.125 MiB | 6.16% [5.85%, 6.29%] | 5.36 deg [1.28, 11.46] |
| 12 x 48 x 144 | 2.96x | 3.797 MiB | 4.01% [3.81%, 4.10%] | 0.98 deg [0.50, 1.37] |
| 16 x 64 x 192 | 2.83x | 9.000 MiB | 2.95% [2.81%, 3.02%] | 0.62 deg [0.40, 0.83] |
| 24 x 96 x 288 | 2.61x | 30.375 MiB | 1.94% [1.84%, 1.98%] | 0.44 deg [0.28, 0.58] |

The spectral representation stored a median 0.093 MiB. Thus, within this
prototype, the grid was 2.6--3.1 times faster per batched evaluation but used
12--327 times more memory. The accuracy frontier was not captured by geometry
error alone. At the coarsest resolution the normalized birth distributions
differed by 6.2% in total variation, and one percent of the source weight saw a
field-axis error exceeding a median 5.4 degrees across devices. Increasing the
grid to 24 x 96 x 288 reduced those errors to 1.9% and 0.44 degrees,
respectively, but required 30.4 MiB per source. The raw Jacobian RMS error also
fell only from 21.1% [19.6%, 21.8%] to 11.9% [11.0%, 12.2%]; this metric is
deliberately secondary because it is dominated by the difficult near-axis
region, whereas total variation measures the normalized birth distribution
more directly.

These results support a practical division of labor. The spectral form is the
fidelity and portability default: it is compact, preserves the equilibrium
truncation, and supplies an analytic field direction for polarized emission.
A grid is attractive when repeated evaluation is demonstrably dominant, but
its resolution should be selected against both a source-distribution tolerance
and a field-axis tolerance rather than against position error alone. An
adaptive grid refined by Jacobian and field-axis error could recover most of
the speed advantage without paying the full uniform-grid memory cost.

This comparison is representation-level evidence, not yet an end-to-end
transport claim. The next publication-grade test is an identical-history
spectral-versus-grid run on at least one QA and one QH configuration, reporting
wall loading and coil flux or heating with Monte Carlo uncertainty. That test
will determine whether the measured 1.9--6.2% birth-distribution difference
and 0.44--5.36 degree field-axis error produce a resolvable engineering-tally
bias. Until then, this subsection should be presented as a quantified method
tradeoff and a motivation for transport validation, not as proof that either
representation changes a device conclusion.

## Evidence and reproduction

- Aggregate: `agg_representations_20dev.json`
- Per-device driver: `compare_multidevice.py`
- Representation implementation: `fluxrep.py`
- Metrics: `metrics.py`
- Unit gates: `test_fluxrep.py`
- Ginsburg equilibrium-generation job: `9585942` (20/22 converged)

The generated equilibrium arrays, caches, and run logs are intentionally
ignored by `spf_validation/.gitignore`; the aggregate JSON retains the complete
per-device numerical record needed to audit the table.
