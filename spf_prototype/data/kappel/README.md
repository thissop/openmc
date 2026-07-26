# Kappel / Landreman / Malhotra — Magnetic Gradient Scale Length dataset

Reactor-scaled stellarator database underpinning the L∇B → plasma–coil-separation
law used as the "geometry → standoff spine" (real n≈40) for SPF Paper-2.

- **Paper:** J. Kappel, M. Landreman, D. Malhotra, *The magnetic gradient scale
  length explains why certain plasmas require close external magnetic coils*,
  Plasma Phys. Control. Fusion **66** (2024) 025018.
- **Zenodo record:** 8349408 (concept) → resolves to version record **8349409**,
  DOI **10.5281/zenodo.8349409**, "Magnetic Gradient Scale Length Auxillary Dataset",
  published 2023-09-15.
- **Downloaded:** 2026-07-26, on the Ginsburg cluster to `~/kappel_zenodo/`.

## File manifest (Zenodo)

Single archive:

| file | size | md5 |
|------|------|-----|
| `2023-09-14.MagneticGradientScalength.tar.gz` | 392,326,242 B (392 MB) | `158030bc5147b722ac04dfe609962397` |

Archive contents (145 entries):
- `finalData_RMSB_0.01_targetMaxK_17.16e6_REGRES.csv` — **the paper dataset** (45 configs).
- `scaleLengths.py`, `MagneticGradient.py`, `infrastructure.py`, `regcoil_in.data`,
  `readme.txt`, `job.mutter` — the pipeline (REGCOIL commit 4813f74, simsopt dae548e).
- `configs/` — 45 VMEC `input.*` files.
- `configswout/` — 45 VMEC `wout_*.nc` equilibria (source of nfp, aspect, betatotal).
- `bnorm/`, `vcasing/` — normal-field / virtual-casing inputs for finite-β configs.

Raw archive is **not** committed (large). Only the parsed CSV + this README live here.
Ginsburg copy of the parsed CSV:
`~/kappel_zenodo/2023-09-14.MagneticGradientScalength/kappel_configs.csv`.

## Parsed table — `kappel_configs.csv`

45 devices, columns:

- `config` — device name (VMEC input basename).
- `nfp` — field periods (read from `wout_*.nc`).
- `qs_type` — heuristic QS class from the config name (QA/QH/QI/tokamak-axisym/
  stellarator-other). **Heuristic only** — a name-based guess, not the paper's own
  taxonomy; do not treat as authoritative. Coverage-check use only.
- `aspect` — VMEC aspect ratio (from `wout`).
- `beta_wout` — VMEC `betatotal` (actual equilibrium β).
- `beta_csv` — β field from the source dataset CSV (target/input β; 0 for vacuum runs).
- `L_REGCOIL_m` — min achievable plasma–coil separation (source column `Regcoil_sep`), metres.
- `L_gradB_m` — min magnetic gradient scale length L*_gradB (source column `L_gradB`), metres.

Mapping note: the source CSV also carries alternate scale lengths
(`L_gradBScalar`, `Curvature`, `L_gradSVEuc`, `L_gradSVMax`); the paper's L*_gradB is
the `L_gradB` column (the min-over-volume ∇B scale length), confirmed by the anchor
check below. All configs are reactor-scaled ("aScaling"), so lengths are directly in metres.

## Verification against paper anchors

| check | expected (paper) | measured | status |
|-------|------------------|----------|--------|
| regression slope, L*_gradB = m·L_REGCOIL + b | 1.585 | **1.5850** | MATCH |
| regression R² | 0.941 | **0.9413** | MATCH |
| precise QH — L_REGCOIL | 1.52 m | **1.5206 m** | MATCH |
| precise QH — L*_gradB | 3.19 m | **3.1882 m** | MATCH |
| precise QA — L_REGCOIL | 2.87 m | **2.8748 m** | MATCH |
| precise QA — L*_gradB | 5.30 m | **5.2977 m** | MATCH |

- Regression is `L_gradB ~ L_REGCOIL` (ordinary least squares over all 45 points):
  slope 1.5850, intercept 0.4857, R² 0.9413 — reproduces the paper's headline law.
- **precise QH** = `20210728-01-026_QH_nfp4_A8_magwell_aScaling` (nfp=4, aspect 8.0).
- **precise QA** = `new_QA_aScaling` (nfp=2, aspect 6.0).
- Note: `20220124_qfm_well0_length24_aScaling` (2.875 / 5.304) and
  `20220124_qfm_well1_length24_aScaling` mirror `new_QA` / `new_QA_magwell` almost
  exactly — QFM reconstructions of the same precise-QA family.
- **Wechsung QA:** no config is unambiguously named as the Wechsung 2022 QA; the
  precise-QA anchor here is `new_QA` (Landreman–Paul precise QA). Not forced to a match.

## Coverage

- n = 45 reactor-scaled configurations.
- nfp spread: 1 (×3), 2 (×13), 3 (×6), 4 (×11), 5 (×7), 6 (×1), 10 (×3), 12 (×1).
- QS/type coverage spans QA, QH, QI, and tokamak/axisymmetric (ITER-hybrid,
  reactor-scaled tokamak `r35x/r36x/r37x_g12538` cases) plus classic stellarators
  (W7-X, HSX, NCSX, TJ-II, CTH, ATF, QPS, CFQS, ESTELL, ARIES-CS, etc.).
- **L*_gradB spread: 0.969 m (`n4qh.b4.a79a`) → 7.536 m (`QI_NFP1_r1_test`) ≈ 7.8×.**
  (Paper describes an "~order of magnitude" / ~10× spread; the measured factor on the
  `L_gradB` column is 7.8×. Reported honestly, not forced.)
- L_REGCOIL spread: 0.406 m → 4.713 m.

## Reproduce

On Ginsburg (has internet, scipy; no netCDF4 — scipy.io.netcdf_file used instead):

```
mkdir -p ~/kappel_zenodo && cd ~/kappel_zenodo
curl -sL https://zenodo.org/api/records/8349409/files/2023-09-14.MagneticGradientScalength.tar.gz/content -o data.tar.gz
tar xzf data.tar.gz
cd 2023-09-14.MagneticGradientScalength
python3 build_kappel.py    # reads finalData CSV + configswout/*.nc -> kappel_configs.csv
```
