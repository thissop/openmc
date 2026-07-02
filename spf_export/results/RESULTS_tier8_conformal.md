# RESULTS — Tier 8 conformal stellarator SPF (Ginsburg)

Equilibrium `equil_precise_qa`, device scale 15.0 (R0~1545cm, a~258cm), conformal DAGMC wall, 10,000 histories/config. Modes: unpolarized=iso, perpendicular=A(sin^2), parallel=B/C(1+3cos^2). Per source neutron, rate held fixed.

> Machinery validation on a public equilibrium (precise_QA), NOT a Helios physics claim. Coil tallies are deep behind the blanket+shield; without variance reduction expect large MC error -- read the +/- and prefer the trend (see DEFERRED.md).

## Global TBR per mode (FLiBe H3-production)

| mode | free-streaming | scattering |
|---|---|---|
| unpolarized | 0.000 ± 0.000 | 0.680 ± 0.008 |
| perpendicular | 0.000 ± 0.000 | 0.687 ± 0.007 |
| parallel | 0.000 ± 0.000 | 0.693 ± 0.015 |

## Coil fast flux (>0.1 MeV) per mode -- magnet-lifetime proxy

| mode | free-streaming [/cm²/src] | scattering [/cm²/src] |
|---|---|---|
| unpolarized | 3.168e+01 (±2%) | 1.145e+01 (±2%) |
| perpendicular | 2.684e+01 (±1%) | 1.016e+01 (±3%) |
| parallel | 3.596e+01 (±1%) | 1.243e+01 (±2%) |

**Directional efficiency η (coil fast flux):** parallel reduces coil flux by -13.5% free-streaming → -8.5% with scattering ⇒ **η = 0.63** (fraction of the steering benefit that survives non-axisymmetric conformal walls + scattering). The headline answer to 'does SPF steering survive'. (Large coil MC error without variance reduction -- treat as indicative; rerun with weight windows for a precise η.)

## Per-material heating (scattering, per source neutron) [eV]

| layer | unpolarized | perpendicular | parallel |
|---|---|---|---|
| W | 1.007e+04 | 9.000e+03 | 1.135e+04 |
| steel | 3.653e+05 | 3.516e+05 | 3.839e+05 |
| Be | 4.351e+05 | 4.387e+05 | 4.351e+05 |
| FLiBe | 5.633e+06 | 5.763e+06 | 5.587e+06 |
| shield | 3.012e+05 | 3.046e+05 | 2.970e+05 |
| coil | 4.145e+05 | 3.720e+05 | 4.579e+05 |

## Not in this table (postprocessing hooks)

- **φ-resolved first-wall load** (toroidal/poloidal map): the `wall_current_phi` mesh tally is saved in each statepoint; extract the poloidal map there.

- **Joint analytic↔MC free-streaming cross-check**: compare the free-streaming per-patch wall load to the analytic stellarator NWL (companion field-period-perturbation code / anarrima angled kernels) on the same precise_QA patches. See JOINT_VALIDATION.md.

