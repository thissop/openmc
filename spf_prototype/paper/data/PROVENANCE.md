# Data provenance

Every number, table, and figure in the manuscript traces to a frozen file here. Each was
extracted from a documented transport run recorded in `spf_prototype/docs/notes/results/`
(the RESULTS_tierX.md ledger) or produced by the conformal-wall pipeline. Figure scripts in
`../src/` read only these files; none recompute transport.

## Frozen tables (CSV)

| file | contents | upstream source | run config |
|---|---|---|---|
| `steering_dilution.csv` | free-streaming vs full-transport steering (mode/iso-1) | `RESULTS_tier5.md` | layered box-torus (W/Fe-9Cr/Be/FLiBe), ARC-class, scattering ON, 4e5 hist/config; thin-wall anchor 1e6 |
| `rate_steering.csv` | absolute rate/load/TBR per mode (rel. unpolarized=1) | `RESULTS_tier7.md` | postprocess of Tier 5/6 statepoints; eta = a+2b/3+c/3 restored |
| `bare_cylinder.csv` | point-patch analytic error vs shaping s | `RESULTS_quasr.md`, `python/bare_cylinder_demo.py` | single shaped filament, SPF kernel, analytic point-patch vs exact ray/surface |
| `quasr59509_validation.csv` | analytic<->OpenMC per-wall A/iso, B/iso | `RESULTS_quasr.md` | real QUASR QA 59509 (nfp=3, eps_eff~0.43), free-streaming, 1e6 births/mode; analytic gated to anarrima 1.2e-12 |
| `emission_moments.csv` | <(u.Bhat_local)^2> about the 3-D field | `RESULTS_tier8_3d.md` | 2e5 flux-surface births, QA field map |
| `predictor_n15.csv` | 15-device conformal-wall predictor table (Y, PF_geo, shape metrics) + Spearman rho | `python/predictor_test.py` | 15 QUASR devices; Ginsburg VMEC batch 8885457 -> conformal batch 8885458, plus original 5 |

## Frozen transport maps

`conformal_maps/quasr<ID>_conformalmap.npz` -- 15 devices. Each holds the (theta,phi) NWL
maps `unpol`, `perp`, `par` and the wall-cell areas `dA`, on the DAGMC conformal first wall
(VMEC LCFS + 0.30a gap). Produced by `python/conformal_wallmap.py` via the native
`StellaratorSource`, 3M wall crossings per mode, on Ginsburg (jobs 8885457 -> 8885458).
These regenerate `fig_conformal_map` and the peaking/benefit numbers (`predictor_n15.csv`,
`benefit()` in `predictor_test.py`).

Devices: 1052272 1090019 11357 1642553 180790 1854549 1960314 24285 262171 59509 66633
803097 883496 886079 932746.

## Not frozen here (large / derivable)

- The VMEC fluxmaps `quasr<ID>_vmec_fluxmap.npz` (~1.7 MB each, source geometry sqrt(g)/R/Z/
  Bhat) live in `spf_prototype/data/` and on Ginsburg at
  `/ginsburg/astro/users/tjk2147/spf_work/spf_pp/data/`. They are needed only to *recompute*
  the geometry predictors in `predictor_test.py`; the derived result is frozen in
  `predictor_n15.csv`.
- Raw OpenMC statepoints remain on Ginsburg; the RESULTS_tierX.md ledger records the numbers.

Regenerate figures: `cd ../src && python fig_*.py` (each writes to `../figures/`).

## Conformal + scattering flux tally

`conformal_maps/quasr59509_conformal_flux.npz` -- volumetric flux tally on the conformal
multi-material DAGMC blanket (device 59509, scaled x10 to reactor size, scattering ON), from
`python/conformal_scatter_59509.py` (Ginsburg job 8912742, native StellaratorSource, 3M/mode).
Holds the poloidal fast-flux mesh `flux_<mode>` (nr,nphi,nz)=(39,32,39) and per-layer heating
`heat_<mode>` for unpolarized/perpendicular/parallel, plus r/z/phi grids and R0/a. Figure
`fig_conformal_flux`. FLiBe = 84% of the 1.1e7 eV/src heating; near-wall fast flux modulates
perp -10% / parallel +10% vs unpolarized.
