# Gil 2026 coil sets → ParaStell/OpenMC neutronics pipeline

Turn the [Gil et al. 2026](https://doi.org/10.1103/PhysRevLett.137.065101) augmented-Lagrangian
stellarator coil archive (Zenodo [18497939](https://doi.org/10.5281/zenodo.18497939)) into
DAGMC-ready coil geometry for magnet-shielding neutronics, and pick a controlled coil-surface-distance
(standoff) sweep from it. Everything here is tested **locally** against synthetic simsopt data — no
481 MB download and no cluster needed to run the test suite.

## Why these coils
86 optimized coil sets across QA (ncoils 3/5) and QH (ncoils 4/5), each a simsopt `BiotSavart` JSON.
The archive varies penalty knobs including **`cs` = coil-surface distance = our fixed blanket
envelope**, so it is an ideal *controlled* input for a standoff sweep, complementary to QUASR's broad
statistical variety. See memory `spf-gil2026-coils`.

## Environment
- **simsopt required** to load the JSONs. Locally: system `python3` has simsopt 1.10.6. On the
  cluster: `~/vmec_venv` (NOT the `spf-stellarator` / `pstl` envs).
- **ParaStell required** only for the final DAGMC build — cluster-only (`pstl` env). Everything up to
  the MAKEGRID file runs anywhere simsopt is present.

## Pipeline
```
Zenodo archive (simsopt BiotSavart JSONs)
   │  batch_gil_to_makegrid.py            (simsopt env)
   ▼
MAKEGRID coils files + manifest.csv       (ParaStell format; per-set R0, cs/cc/…, nfp)
   │  standoff_sweep_planner.py           (pure python)
   ▼
sweep plan (which sets span the cs knob)
   │  ParaStell construct_magnets_from_filaments  (pstl env, CLUSTER)
   ▼
DAGMC magnet geometry → OpenMC transport (coil_run_v3 pipeline)
```

## Steps

**1. Convert one set** (`gil_to_makegrid.py`)
```bash
python gil_to_makegrid.py IN.json coils_out --nfp 4 --target-major-radius 7.75
```
`--target-major-radius R0[m]` measures the set's own R0 (mean cylindrical R of coil centroids) and
auto-computes the rescale — "put this set at ARIES-CS scale" in one flag. Or `--scale S` for a raw
ratio. ParaStell needs only geometry, so the current is irrelevant to the DAGMC build.

**2. Convert the whole archive + manifest** (`batch_gil_to_makegrid.py`)
```bash
python batch_gil_to_makegrid.py ~/gil_coils OUT_DIR --nfp-map QA=2,QH=4 [--target-major-radius 7.75]
```
Writes one `coils_*` per set and `manifest.csv` (symmetry, ncoils, measured R0, applied scale, cs,
cc, length, curvature, msc, force, status). A bad set becomes `status=ERROR` and is skipped, never
fatal. `--nfp-map` sets field periods per symmetry (verify QA/QH nfp against the archive metadata).

**3. Plan a controlled standoff sweep** (`standoff_sweep_planner.py`)
```bash
python standoff_sweep_planner.py OUT_DIR/manifest.csv --symmetry QH --n-points 5 --out plan.csv
```
Finds the (symmetry, ncoils) family with the widest usable `cs` span and picks `n_points` sets whose
cs values best cover the range (endpoints kept). The plan names exactly which `coils_*` files to feed
ParaStell for a sweep where the standoff is the only thing changing.

**4. Inspect / reuse a layout** (`read_makegrid.py`)
```bash
python read_makegrid.py coils_out          # prints per-coil R, phi, Z
```
Also drives the surrogate placement demo on a **real** coil arrangement:
```bash
python ../shield_opt/multicoil_placement_demo.py --gil-coils coils_out
```
(The dose model there stays synthetic; only the coil *positions* become real.)

**5. Build DAGMC** — cluster/`pstl` env, ParaStell `construct_magnets_from_filaments(coils_out, ...)`,
then the standard `coil_run_v3` transport. (Not covered by the local test suite.)

## Tests
```bash
python -m pytest .        # 19 tests, all offline (synthetic simsopt round-trips)
```
Covers converter round-trip + format, auto-scale, batch config-token parsing + bad-file isolation,
MAKEGRID reader (angle/Z/count recovery), and sweep-family selection + cs coverage.

## Caveats
- `--nfp-map` and the config-token parsing follow the archive's own directory scheme; **verify nfp
  and cs units against the archive metadata** before a production run (parsing never blocks a
  conversion — unrecognised tokens land in `manifest.raw_tokens`).
- Coils are ~1 m major radius as archived; rescale to reactor size before building.
- Local tests use synthetic BiotSavart sets; the real archive lives on the cluster at `~/gil_coils/`.
