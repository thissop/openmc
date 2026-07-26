# magnet_suite — per-device magnet-side neutronics automation

Parameterizes the hand-run QH/QA adjoint-shielding pipeline into a single per-device driver so
it can scale to the 200–300 engineering-relevant QI/QA/QH devices. Nothing here is rewritten
from scratch — the stages call the existing, debugged scripts (`equil` reuses
`python/quasr_fluxmap.py`; `dagmc` mirrors `pstl_test/corrected/build_corrected.py`; `cells`
reuses `step1/step1_cells.py`; `adjoint` reuses `fwcadis/adjoint_importance.py`; `finalize`
reuses `adjoint_placement.py` + `concentration.py` + `emissivity.py`).

## Pipeline (one QUASR device ID -> one JSON record)

| stage | script | env | output artifact |
|-------|--------|-----|-----------------|
| fetch | `make_coils.py` | numpy + internet | `<label>.coils` (reactor-scaled MAKEGRID, cm), `<label>_coil_centroids.npz`, `<label>_reactor_scale.txt` |
| equil | `equil_device.py` | `~/desc_venv` (DESC 0.17.2) | `<label>_fluxmap.{npz,meta,bin}`, `wout_<label>.nc` (VMEC, for ParaStell) |
| dagmc | `build_device_dagmc.py` | conda `…/envs/pstl` (ParaStell) | `dagmc_<label>.h5m` (8-layer, ~1.29 m) |
| cells | `step1/step1_cells.py` | conda `spf-stellarator` (pymoab) | `<label>_cells.npz` (magnet cell IDs = GLOBAL_ID+1) |
| adjoint | `fwcadis/adjoint_importance.py` | `openmc-spf` build | `adjoint/adjoint_importance_{flat,kerma}_P0.npz` |
| finalize | inline (`run_device.finalize`) | numpy | `record.json` |

`run_device.py` orchestrates: **fetch** runs inline (login node — needs QUASR internet and
writes the shared reactor-scale factor), then **equil→dagmc→cells→adjoint** are SLURM jobs
chained with `--dependency=afterok`, then **finalize** runs inline once the adjoint maps exist.

### The wout link (why fresh devices now reach DAGMC)
QH/QA had hand-obtained VMEC wouts; a zoo device has only a QUASR boundary+coils. `equil_device.py`
solves a DESC fixed-boundary equilibrium from the boundary and `desc.vmec.VMECIO.save` writes the
VMEC wout ParaStell consumes — the same equilibrium also produces the transport/adjoint fluxmap, so
DAGMC geometry and the source geometry are consistent. **Reactor scaling:** the boundary Fourier
modes and the coils are both scaled by `f = 1.704 / a_native` (catalogue `minor_radius`), the single
authoritative factor in `<label>_reactor_scale.txt`, so the ~1.29 m absolute radial build is
proportioned to a reactor plasma and coils share the wout frame.

## Bosch-Hale emissivity (real DT reactivity, consistent across attribution AND transport)
`bosch_hale_emission.py`. The native `StellaratorSource` births neutrons with per-cell weight
`w = S(rho)·|sqrt(g)|` (C++ `build_cdf_cascade`); `S(rho)` defaults to the crude `1-rho^2` and is
overridden by `emission=(rho, density)`. Setting `S(rho)=r_BH(rho)=n(rho)^2<σv>_DT(T(rho))`
(from `emissivity.py`) makes the transport birth density the **real DT reaction-rate density** —
the *same* `r(rho)` the adjoint contributon already folds into `S(r)` via
`adjoint_placement.plasma_source_on_mesh(..., emissivity='bosch_hale')`. So one emissivity drives
both attribution and transport.

    from bosch_hale_emission import bosch_hale_emission
    rho, dens = bosch_hale_emission("<label>_fluxmap.npz", "bosch_hale")
    src = openmc.StellaratorSource(fluxmap=stem, emission=(rho, dens), ...)

Verified core-peaking (`python bosch_hale_emission.py`, device 112734 fluxmap):
`<rho>` births move 0.558 (uniform) → 0.405 (1-rho²) → **0.307 (Bosch-Hale)**; core fraction
(rho<0.5) 48% → 70% → **88%**. BH is strongly core-peaked, as physically expected.

## Reproduce

    # one fresh device (dry plan; inline fetch + finalize only):
    python run_device.py 1328722 --work /scratch/msuite --suite $PWD --dry-run
    # real run on Ginsburg (submits the SLURM chain):
    python run_device.py 1328722 --work ~/pstl_test/magnet_suite --suite $PWD --submit
    # small validation batch (measure per-device wall-clock before the full array):
    awk -F, 'NR>1 && $7==1 {print $1}' ../data/engineering_relevant_devices.csv > devices.txt
    head -5 devices.txt > devices_batch.txt
    sbatch --array=1-5%5 slurm_array.sh devices_batch.txt

## Per-device JSON (`record.json`) fields
Device descriptors (nfp, aspect, qs_class, d_min_over_a, gap_m, blanket_fit) + per-response
(`flat`, `kerma`) attribution: `gini`, `participation_ratio`, `peaking`, `nnz_frac`,
`median_relerr`, placement `(phi0,theta0,R0)`, `peak_coil_index` + `peak_coil_dphi_deg` +
`placement_points_at_coil`, and the structure gate (`gate_localized`, `gate_structured`,
`gate_relerr_ok`, `gate_pass`).

## HARD GATES (honest, no hand-tuning)
- Adjoint structure gate: localized (`nnz_frac<0.98`), structured (`peaking>1.5`),
  converged (`median_relerr<0.25`), placement points at the coil (`dphi<45°`).
- DAGMC watertight (ParaStell export); coil-cell bijection (pymoab GLOBAL_ID+1).
- DESC volume gate `relV<3e-2` (equil).
- Do **not** launch the full 200–300 array until (a) the automation passes end-to-end on a fresh
  device and (b) the adaptive-shielding kill-shot verdict validates the shield step.

## Selection-effect caveat
The device list is the blanket-fit-filtered engineering-relevant set. That filter is NOT
population-neutral (as-shipped → high aspect; min-achievable → low nfp/QI — see
`../RESULTS_engineering_relevance.md`). **Stratify/control every result by aspect and nfp.**
