# Data provenance

Every number in the manuscript traces to a frozen file here. Each was extracted from a documented
transport run recorded in `spf_prototype/docs/notes/results/` (the RESULTS_tierX.md ledger) or from a
saved `.npz` produced by the conformal-wall pipeline. No figure script recomputes transport; they read
only these files.

| file | contents | upstream source | run config |
|---|---|---|---|
| `steering_dilution.csv` | free-streaming vs full-transport steering (mode/iso−1) | `RESULTS_tier5.md` | layered box-torus (W/Fe-9Cr/Be/FLiBe), ARC-class, scattering ON, 4e5 hist/config; thin-wall anchor 1e6 |
| `rate_steering.csv` | absolute rate/load/TBR per mode (rel. unpolarized=1) | `RESULTS_tier7.md` | postprocess of Tier 5/6 statepoints; η = a+2b/3+c/3 restored |
| `bare_cylinder.csv` | point-patch analytic error vs shaping s | `RESULTS_quasr.md`, `python/bare_cylinder_demo.py` | single shaped filament, SPF kernel, analytic point-patch vs exact ray/surface |
| `quasr59509_validation.csv` | analytic↔OpenMC per-wall A/iso, B/iso | `RESULTS_quasr.md` | real QUASR QA 59509 (nfp=3, ε_eff≈0.43), free-streaming, 1e6 births/mode; analytic gated to anarrima 1.2e-12 |
| `emission_moments.csv` | ⟨(u·B̂_local)²⟩ about the 3-D field | `RESULTS_tier8_3d.md` | 2e5 flux-surface births, QA field map |
| `quasr886079_conformalmap.npz` | (θ,φ) NWL maps unpol/perp/par + dA on the DAGMC conformal wall | conformal-wall pipeline (`python/conformal_wallmap.py`) | VMEC LCFS + 0.30a gap, DAGMC first wall, StellaratorSource, 3M crossings/mode |

Regenerate figures: `python src/fig_*.py` (each writes to `figures/`).
