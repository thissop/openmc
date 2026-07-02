# NOTES_orientation.md — what actually exists in this checkout vs the brief

Ground truth for the SPF QS-stellarator sweep engine, established by reading the
committed code on branch `spf-prototype`. Read this before trusting the brief's
file names: the brief describes a **newer** pipeline (a later Ginsburg session)
that is **not fully synced into this local checkout**.

## Existing, trusted, REUSE (paths under spf_prototype/python/ unless noted)
- `stellarator_geometry.py::build_layers(stem, layers=DEFAULT_LAYERS, scale=1.0, outdir=None)`
  -> writes per-layer conformal shell STLs + `manifest.json`. Pure numpy. Manifest
  layers carry: name, inner_offset_cm, outer_offset_cm, stl, watertight,
  simple_cross_section, n_self_intersecting_phi, enclosed_volume_cm3. Reads the LCFS
  from `data/{stem}_surface.npz` (keys R(ntheta,nphi), Z, nfp). HARD self-intersection
  guard already built in (`cross_sections_simple`, `is_edge_manifold`).
- `run_conformal.py`
  - `MODES = {unpolarized:(1/3,1/3,1/3), perpendicular:(1,0,0), parallel:(0,1,0)}`
  - `build_model(abc, h5m_path, fieldmap_stem, R0_cm, a_cm, particles, batches, density_scale)`
    -> (openmc.Model, mats). Source is the compiled `.so` via
    `openmc.CompiledSource(SO, parameters="a=..,b=..,c=..,bmode=fieldmap,fieldmap=<path>,shape=plasma,R0=..,aminor=..")`.
    Tallies already built: `wall_current_phi` (CylindricalMesh surface current, 32
    toroidal bins — the phi-resolved FIRST-WALL signature), `cell_response`
    (MaterialFilter: heating/damage/H3-production), `coil_fast` (coil flux >0.1 MeV).
  - `scale_fieldmap(stem, scale, out_stem)` — rescales a field map's grid bounds.
  - `make_geometry(h5m)` — loads DAGMC, forces Material.name == layer tag.
  - Free stream = `density_scale=1e-4` (near-void, pre-blanket physics);
    scatter = `1.0`.
- `run_ginsburg.py` — the SLURM entrypoint. `bootstrap()` (build .so, geometry,
  .h5m, scaled fieldmap), `run_all()` (6 streams), `extract_metrics()`, writes
  `RESULTS_tier8_conformal.md`. Resumable-ish. **CAVEAT below.**
- `build_dagmc.py::build(geom_dir, out_h5m)` — STL -> DAGMC `.h5m` via the
  `stl_to_h5m` conda package. HARD gate: refuses non-watertight/self-intersecting
  layers.
- `fieldmap.py::FieldMapField(stem)` — Python reader for the `.so`'s field-map
  format `spf_fieldmap_v1`. **This is the interop contract** (see format below).
- `quasr_geom.py` — loads a QUASR device from its VMEC boundary namelist
  (`data/quasr/nml/input.<ID>`: RBC/ZBS harmonics + NFP). `Device.RZ(theta,phi,rho)`,
  `Device.bhat(theta,phi,rho)` (flux-surface-tangent field direction from boundary
  shape + iota), `Device.R0`. `eps_eff(dev,...)` geometric-gentleness metric.
- `quasr_config.py` — problem def for device 59509; `data/quasr/{catalogue.csv.gz,
  shortlist.csv,eps_eff_results.csv}` give per-device metadata.
- Compiled source: `spf_prototype/src/build/libpolarized_fusion_source.so` (present;
  source `.cpp` in `spf_prototype/src/`). `build_and_run.py` (`br.SO`, `br.build_so()`)
  builds/locates it.
- HPC: `environment.yml` (conda-forge `openmc=0.15.*=dagmc*`), `ginsburg_job.sh`
  (sbatch template), `ginsburg_preflight.sh`, `SETUP_GINSBURG.md`,
  `RUN_ON_GINSBURG.md`, `GINSBURG_CLAUDE_HANDOFF.md`.

## Field-map interop contract (spf_fieldmap_v1) — how a NEW field backend plugs in
`.meta` (text k v): magic `spf_fieldmap_v1`, ndim 3, nR, nphi, nZ, R_min, R_max,
phi_min, period, Z_min, Z_max, nfp, endian little, source <str>.
`.bin`: little-endian float64, length `3*nR*nphi*nZ`, three contiguous blocks
[bx][by][bz], each C-order reshape (nR, nphi, nZ). Read code normalizes to unit at
lookup (so storing raw B or b_hat both work). **Any field backend that writes this
pair drives the trusted `.so` unchanged.** `coherence_metrics` reads b_hat from the
SAME map via `fieldmap.FieldMapField.bhat`, so C and the neutronics source use an
identical field by construction.

## Ginsburg environment (for the runbook; from GINSBURG_CLAUDE_HANDOFF.md)
- Account `astro`; `source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh`;
  env `spf-stellarator`. HOME=/burg-archive/home/tjk2147 (~37 GB, cold); no scratch;
  prefer warm /burg if roomy.
- No shared cross-section lib — download ENDF/B-VIII.0 HDF5 once; set
  `OPENMC_CROSS_SECTIONS`. LOGIN NODE = network installs only; ALL compute on a
  compute node (salloc/srun) or sbatch. Partitions: `short` (12 h), `burst` (14 d).
- Repo-fork shadow: there is an `openmc/` dir at repo root (the C++ source tree).
  Must import the CONDA openmc, not the fork — run from a neutral cwd / set sys.path.

## DISCREPANCIES vs the brief (decisions needed — see below)
1. **`dagmc_writer.py` (pymoab-free `build_from_stls`) is NOT here.** This checkout
   has `build_dagmc.py` using the `stl_to_h5m` package. The brief's pymoab-free
   writer exists only on the Ginsburg copy. (Handoff even flags `stl_to_h5m` for API
   drift.) -> sweep will call a swappable `build_dagmc_h5m()` shim that prefers a
   `dagmc_writer.build_from_stls` if present, else `build_dagmc.build`.
2. **Array-based (no-pandas) `_mat_score` is NOT here.** This checkout's
   `run_ginsburg._mat_score` uses `get_pandas_dataframe()` — which the brief says is
   broken on Ginsburg's pandas. -> the sweep will NOT reuse that; it reads tallies
   via the openmc StatePoint array API directly (no pandas), independent of which
   `run_ginsburg` is on Ginsburg.
3. **No coil filaments, no Biot-Savart kernel present.** `quasr_geom` downloads only
   the VMEC BOUNDARY namelist and models b_hat as the flux-surface tangent (an
   O(shaping) approximation). QUASR devices DO have coils (catalogue has
   `nc_per_hp`, `Nfourier_coil`, `total_coil_length`), but the coil geometry (simsopt
   device JSON) is a SEPARATE download not yet fetched, and the brief's "HTS
   Biot-Savart kernel" is not in this checkout. -> exact-coil-BS (brief Path A) is
   buildable and locally testable, but needs the coil-JSON fetch. Decision below.
4. **`export_for_plots.py` is NOT here** (it produced the `spf_export/` bundle on
   Ginsburg). Not needed for the sweep; `analyze_sweep.py` reads per-config JSON.

## KEY INTERPRETATION the brief leaves implicit (frame of the order parameter C)
`C = |∫ s b_hat dV| / ∫ s dV` gives **C=1 for a tokamak** ONLY if b_hat is expressed
in the LOCAL CYLINDRICAL basis (e_R, e_phi, e_Z): a purely toroidal field is
b_hat_cyl=(0,1,0) everywhere -> resultant length 1. In the LAB (fixed xyz) frame a
toroidal field winds around the torus and integrates to ~0, so lab-frame C would be
~0 for tokamak AND stellarator and could not be the discriminator. **Decision: define
C, the direction tensor T, and the angular std on b_hat expressed in the local
cylindrical basis**, so the tokamak limit calibrates to C=1 and stellarator helical
/pitch variation drives C<1. This is documented in THEORY.md and is the single most
important definitional choice; flag if lab-frame or field-line-Frenet was intended.

## eta observables mapped onto EXISTING tallies (no new transport tallies needed)
- `eta_source` (PRIMARY, pre-blanket, fast): directional anisotropy of the
  FREE-stream first-wall current `wall_current_phi` (free stream is near-void, so it
  is genuinely pre-blanket). Fractional change unpol->polarized, normalized to the
  C≈1 anchor.
- `eta_coil` (deep, VR-pending): SCATTER-stream `coil_fast`. `A = eta_coil/eta_source`.
- Anchor for normalization: a synthetic purely-toroidal / axisymmetric-limit config
  with C=1 exactly, included in the family.

## Proposed layout (new work only; imports existing modules, does not fork them)
```
spf_prototype/sweep/
  NOTES_orientation.md   quasr_loader.py     biotsavart_field.py   field_audit.py
  coherence_metrics.py   sweep.py            analyze_sweep.py      sweep.sbatch
  configs/  tests/  docs/{THEORY,FACTORIZATION,EXPERIMENTAL_DESIGN,SCOPE_AND_CAVEATS,README}.md
```
