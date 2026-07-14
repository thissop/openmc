# CURRENT_PIPELINE.md — SPF Stellarator Neutronics Takeover Inventory

Brief section 15, task 1: a path-specific inventory that lets a new engineer **run the
existing pipeline** and identifies the **reuse targets** (do NOT reimplement) for the
hotspot-driven nonuniform-shield-optimizer project (S15.2).

All paths are relative to the repo root
`/Users/tkiker/Documents/GitHub/openmc/` unless noted. The project tree lives under
`spf_prototype/`. Repo commit at time of writing: `2b842f1c2` on branch `spf-prototype`.
Signatures below were read from the actual files, not guessed.

Cluster context (everything transport runs on **Ginsburg**, x86):
- conda env `spf-stellarator` (`spf_prototype/environment.yml`): `openmc=0.15.*=dagmc*`,
  `dagmc`, `python=3.11`. This is the PR-3999 (`TokamakSource`) OpenMC built **WITH**
  `-DOPENMC_USE_DAGMC=ON`, and patched with the native `StellaratorSource`.
- `OPENMC_CROSS_SECTIONS` -> ENDF/B-VIII.0 HDF5.
- Ginsburg working copy: `/ginsburg/astro/users/tjk2147/spf_work/spf_pp/` (mirror of
  `spf_prototype/`), with `data/` holding the VMEC fluxmaps.
- The native `StellaratorSource` is a **compiled OpenMC core class** — it must be
  registered into the OpenMC source tree and the stack rebuilt (see §1.4). It is NOT a
  runtime-loaded `CompiledSource.so`.

---

## Overview diagram (data + control flow)

```
QUASR portal (flatironinstitute.org)
   |  quasr_geom.fetch_nml (VMEC boundary RBC/ZBS)      -> data/quasr/nml/input.<ID>
   |  quasr_loader.fetch_serial (simsopt JSON coils)    -> data/quasr/serials/serial<ID>.json
   v
[EQUILIBRIUM]  quasr_fluxmap.build(ID)  (DESC venv, SEPARATE from OpenMC)
   |   solves fixed-boundary DESC equilibrium from the QUASR boundary
   |   -> data/quasr<ID>_vmec_fluxmap.npz   (rho,theta,zeta, sqrtg, R,Z,phi, BR,Bphi,BZ, nfp)
   |   -> <stem>.meta + <stem>.bin          (spf_fluxmap_v1 = the SOURCE input)
   v
[SOURCE]  native StellaratorSource(fluxmap=<stem>, polarization=(a,b,c), field_model=...)
   |   births: rejection-free CDF cascade over S(rho)*sqrt(g); weight == 1
   |   directions: Schwartz P2 spin-polarized mixture about local b-hat from same grid
   v                                          (matched pol/unpol PAIR = same fluxmap, abc swapped)
[GEOMETRY]  two independent builders:
   |   (A) conformal blanket : stellarator_geometry.build_layers -> per-layer STL shells
   |                           -> dagmc_writer.build_from_stls / build_dagmc.py -> <ID>.h5m
   |   (B) coil solids       : coil_geometry (RMF sweep) -> data/quasr/coil_geom/<ID>_*.stl
   |   (B') coil standoff    : coil_standoff.py -> data/quasr/coil_standoff.csv (filter)
   v
[TRANSPORT]  OpenMC + DAGMCUniverse(h5m), fixed source
   |   conformal_scatter_59509.py  (scattering ON, volumetric flux + per-layer heat)
   |   conformal_wallmap.py        (free-streaming NWL maps, surf_source_write)
   |   run_conformal.py            (box-torus radial-build rig, free vs scatter)
   |   geometry_peaking.py         (LOCAL analytic free-streaming peaking, no transport)
   v
[RESULTS]  frozen -> paper/data/*.csv, paper/data/conformal_maps/*.npz
           figures -> paper/src/fig_*.py -> paper/figures/*.pdf
```

Two distinct source implementations exist and must not be confused:
1. **native `StellaratorSource`** (`native_spf/stellarator_source.*`) — the current,
   real-equilibrium, conformal source. **This is the reuse target.** Consumes
   `spf_fluxmap_v1`.
2. **legacy box-torus `CompiledSource`** (`src/`, driven by `build_and_run.py`) — the
   Tier-2b analytic prototype; a runtime `.so` parsed from a parameter string. Used only
   by `reactor_model.py` / `stellarator_model.py` (box-torus rigs). Superseded for real
   devices, but still the driver of the axisymmetric validation rigs.

---

## 1. THE SPF SOURCE (reuse target — do NOT reimplement)

### 1.1 Native `StellaratorSource` (C++ core class)

| file | key symbols | what it does |
|---|---|---|
| `spf_prototype/native_spf/stellarator_source.h` | `class StellaratorSource : public openmc::Source`; ctor `explicit StellaratorSource(pugi::xml_node)`; `SourceSite sample(uint64_t* seed) const override` | Header. Declares grid members (`rho_,theta_,zeta_,sqrtg_,R_,Z_,phi_,BR_,Bphi_,BZ_,S_`), CDF cascade (`cdf_rho_,cdf_zeta_,cdf_theta_`), SPF weights (`w_perp_,w_par_,p_perp_,polarized_,field_model_`). All state set in ctor, const at sample time (thread-safe, `openmc::prn` only). |
| `spf_prototype/native_spf/stellarator_source.cpp` | `read_fluxmap`, `build_cdf_cascade`, `interp`, `upper_bound_strided`, `sample_polarized_direction`, `sample` | Impl. Line-for-line C++ port of the Python reference `ConformalStellaratorSampler`. |
| `spf_prototype/python/stellarator_source.py` | `class ConformalStellaratorSampler`; `circular_torus_fluxmap(...)` | **Python reference oracle** (pure numpy) that the C++ mirrors bit-for-bit. Use for validation and prototyping without a rebuild. |
| `spf_prototype/native_spf/stellarator_source.py` | `class StellaratorSource(openmc.SourceBase)`; `spin_fractions_to_abc(d_plus,d_zero,d_minus,t_plus,t_minus)` | **Python API** (the object you instantiate). Serializes to/from XML; the C++ core reads that XML + the fluxmap. |
| `spf_prototype/native_spf/test_stellarator_source_cpp.py` | pytest | NumPy-mirror unit test: cascade CDFs, same-seed->same-cell, interp, weight==1. Runs with no build. |
| `spf_prototype/native_spf/README_stellarator.md` | — | The authoritative build/registration + validation-gate doc. |

**Physics (Schwartz 2025, Eq. 2), constants derived in the ctor (not pasted):**
`dσ/dΩ = (σ₀/2π)[ (3/4)a·sin²θ + ((2/3)b+(1/3)c)(1/4+3/4cos²θ) ]`;
`w_perp = 0.75a`, `w_par = (2/3)b+(1/3)c`, `η = a+(2/3)b+(1/3)c`, `P_perp = a/η`.
θ measured from local **b-hat**. `∫sin²θ dΩ=8π/3`, `∫(1/4+3/4cos²θ)dΩ=2π`.

**Sampling algorithm** (`sample`, RNG draw order — must not be reordered):
1. cascade cell select rho→zeta→theta via `upper_bound_strided` = `searchsorted(side="right")−1`, clamped.
2. in-cell jitter rho→theta→zeta.
3. trilinear interp of R,Z,phi (theta inner→zeta→rho outer); rho non-periodic/clamped, theta/zeta `wrap_index` period 2π.
4. b-hat: `field_model_==0` (fluxmap) interpolates (BR,Bphi,BZ)→Cartesian→normalize; `==1` (toroidal) uses pure φ̂=(−sinφ,cosφ,0).
5. direction: isotropic if unpolarized else `sample_polarized_direction` (stratified P2: selector→shape-u→φ; Gram-Schmidt rotation via least-aligned reference axis in `spf_costheta_perp/par` inverse-CDF helpers).
6. energy (single dist for all radii) then time. `site.wgt = E_wgt*t_wgt` (==1 for defaults).

### 1.2 Instantiating from Python

```python
import openmc  # the StellaratorSource is registered into openmc.* on the patched build
src = openmc.StellaratorSource(
    fluxmap="data/quasr59509_cm_fluxmap",     # stem of spf_fluxmap_v1 .meta/.bin (.meta/.bin auto-stripped)
    polarization=(1/3, 1/3, 1/3),             # (a,b,c) collision-mode fractions; None => isotropic
    field_model="fluxmap",                    # or "toroidal" (axisymmetric-limit cross-check)
    energy=openmc.stats.Discrete([14.1e6], [1.0]),  # default 14.06e6 if omitted
    time=None, emission=None, strength=1.0, constraints=None)
```
Python-API signature (`native_spf/stellarator_source.py:94`):
`StellaratorSource(fluxmap, polarization=None, field_model='fluxmap', energy=None, time=None, emission=None, strength=1.0, constraints=None)`.
- `polarization`: `None` (isotropic, writes no XML — backward compatible), a 3-seq `(a,b,c)`,
  or a dict `{d_plus,d_zero,d_minus,t_plus,t_minus}` mapped via `spin_fractions_to_abc`.
  Fractions must be ≥0; renormalized to sum 1 with a warning if the input sum deviates >1e-6.
- `emission`: optional `(rho_grid, density)` overriding the default parabolic `S(rho)=1−rho²`.

### 1.3 The `spf_fluxmap_v1` format (writer + reader)

Writer: `spf_prototype/python/quasr_fluxmap.py:write_fluxmap_bin(stem, rho, theta, zeta,
sqrtg, R, Z, phi, BR, Bphi, BZ, nfp, sign_sqrtg=1.0)`. HDF5-free so the C++ reader is a
plain `std::ifstream`.
- `<stem>.meta` (text `key value`/line): `magic spf_fluxmap_v1`, `nR`, `ntheta`, `nzeta`, `nfp`, `sign_sqrtg`.
- `<stem>.bin` (little-endian f8, C-order, exact order): `rho[nR], theta[ntheta], zeta[nzeta],
  sqrtg, R, Z, phi, BR, Bphi, BZ` (each `nR*ntheta*nzeta`, indexed `[i_rho,j_theta,k_zeta]`).
- `theta,zeta` span `[0,2π)` (endpoint excluded); `rho` runs `(drho..1]`; `sqrtg` stored as
  `|sqrt(g)|`. The `write_fluxmap_bin` used by conformal_wallmap/conformal_scatter is the
  same function (also re-exported into `conformal_wallmap.py`).

### 1.4 Registering / building the source (from README_stellarator.md §3)

```bash
cp native_spf/stellarator_source.h   <openmc>/include/openmc/
cp native_spf/stellarator_source.cpp <openmc>/src/
# CMakeLists.txt: add "src/stellarator_source.cpp" next to src/source.cpp
# src/source.cpp: #include header + add Source::create branch:
#     } else if (source_type == "stellarator") { return make_unique<StellaratorSource>(node); }
# openmc/source.py: merge StellaratorSource + spin_fractions_to_abc; add from_xml_element branch
# openmc/__init__.py: from openmc.source import spin_fractions_to_abc
cd <openmc>/build && cmake --build . -j && pip install -e ..
```
`StellaratorSource` is self-contained (carries its own `spf_costheta_*` helpers); it does
NOT require the PLAN-A `TokamakSource` SPF patch to be applied first.

### 1.5 Matched polarized/unpolarized PAIR

A matched pair is the **same fluxmap + same emission profile + same energy + same
strength**, with only `polarization` changed. Every run script builds the three canonical
modes over one fluxmap:
```python
MODES = {"unpolarized": (1/3,1/3,1/3), "perpendicular"/"perp": (1.0,0.0,0.0),
         "parallel"/"par": (0.0,1.0,0.0)}   # Bae naming: unpol=iso, perp=A (sin²θ), par=B/C (1/4+3/4cos²θ)
```
Because the CDF cascade gives birth weight exactly 1 for every mode and the P2 kernel
integrates to the same neutron count, comparing modes isolates DIRECTION (rate factor η
restored separately). This is the normalization guarantee that makes pol vs unpol a fair
pair. Files: `conformal_scatter_59509.py:23`, `conformal_wallmap.py:32`, `run_conformal.py:41`.

### 1.6 SPF-patched `TokamakSource` (axisymmetric sibling, also a reuse target)

| file | key symbols | what it does |
|---|---|---|
| `spf_prototype/native_spf/source.cpp` | `TokamakSource::TokamakSource(pugi::xml_node)` (line 721); `sample_polarized_direction` (1304); `field_direction` (1259); `spf_costheta_perp/par` (702/710) | PR-3999 `TokamakSource` **patched** with the SPF direction sampler. `.orig` sibling is the unpatched upstream. Same `w_perp/w_par/p_perp` derivation as `StellaratorSource`. |
| `spf_prototype/native_spf/source.py` | `class TokamakSource(SourceBase)` (934); `__init__(..., polarization=None, field_model='toroidal', safety_factor=None, field_sign=1)` (1056); `spin_fractions_to_abc` (901) | Python API. `field_model ∈ {'toroidal','pitched'}`; `'pitched'` requires `safety_factor=(r_over_a,q)`. `polarization=None` => unpolarized (isotropic, no XML). |
| `spf_prototype/native_spf/README_apply.md` | — | How to apply the `TokamakSource` SPF patch onto a pr3999 checkout. |
| `spf_prototype/native_spf/test_spf_tokamak.py` | pytest | SPF `TokamakSource` tests. |

Use `TokamakSource(field_model='toroidal')` as the analytic axisymmetric oracle to
cross-check `StellaratorSource` (feed an axisymmetric fluxmap + `field_model='toroidal'`;
should reproduce Schwartz ±43% inboard / ∓22% outboard).

---

## 2. EQUILIBRIUM / VMEC / BOUNDARY LOADERS

| file | key symbols | what it does | run command |
|---|---|---|---|
| `spf_prototype/python/quasr_geom.py` | `fetch_nml(ID)`, `parse_nml(fn)`, `class Device(spec,iota)` with `RZ(theta,phi,rho)` & `bhat(theta,phi,rho)`, `load_device(ID,iota)`, `eps_eff(dev,...)` | Pure-numpy VMEC-boundary loader. Evaluates `R=ΣRBC cos(mθ−nNFPφ)`, `Z=ΣZBS sin(...)` from the QUASR namelist; `bhat` = flux-surface tangent winding with iota. No simsopt. | `python python/quasr_geom.py <ID>` |
| `spf_prototype/sweep/quasr_loader.py` | `fetch_serial(ID)`, `load_coils(ID,n_samples=256)` -> `[(poly(N,3),current)]`, `total_filament_length`, `class QuasrDevice` with `source_sample(...)`, `boundary_grid(ntheta,nphi,rho)`, `load_device(ID,n_samples=256,coils=True)` -> `QuasrDevice(coils, device, meta)` | Parses the simsopt-serialized JSON coil graph (CurveXYZFourier / RotatedCurve / ScaledCurrent) WITHOUT importing simsopt; wraps `quasr_geom` for the boundary. `meta` = nfp,iota,aspect,minor_radius,qs_error,nc_per_hp,total_coil_length,symmetry_class(QA/QH),R0. | `python sweep/quasr_loader.py <ID>` |
| `spf_prototype/python/quasr_fluxmap.py` | `build(ID, outdir=None, L=8, n_rho=16, n_theta=64, n_zeta=192, ftol=1e-2, maxiter=100, example=None)`; `write_fluxmap_bin(...)` | **DESC equilibrium producer.** Solves fixed-boundary DESC from the QUASR boundary (or uses a DESC example like `precise_QA`), computes REAL `sqrt(g),R,Z,phi,B` on one `(rho,theta,zeta)` grid, asserts rigor gates V1/V1b/V2/V3/V4, writes `.npz` + `spf_fluxmap_v1` `.meta/.bin`. | `$HOME/desc_venv/bin/python python/quasr_fluxmap.py <ID> [outdir]` (**DESC venv, NOT the OpenMC env**) |
| `spf_prototype/python/quasr_equilibrium_field.py` | `_quasr_boundary_modes(ID)`, `_desc_surface(nfp,Rmodes,Zmodes)` | Boundary->DESC-surface machinery reused by `quasr_fluxmap`. |
| `spf_prototype/python/vmec_fluxmap.py` | fluxmap builder variant | Produces the `quasr<ID>_vmec_fluxmap.npz` used by the wall-map scripts. |

**Cached data locations:**
- `spf_prototype/data/quasr/nml/input.<ID>` — VMEC boundary namelists.
- `spf_prototype/data/quasr/serials/serial<ID>.json` — simsopt coil serials (~366 cached).
- `spf_prototype/data/quasr/catalogue.csv.gz`, `database.json.gz` — device metadata.
- `spf_prototype/data/quasr<ID>_vmec_fluxmap.npz` — per-device equilibrium (~1.7 MB, ~25 devices frozen).
- `spf_prototype/data/quasr<ID>_fluxmap.{meta,bin}` / `..._cm_fluxmap.*` — the source-input form (cm-scaled variants written on-the-fly by the run scripts).
- `spf_prototype/data/quasr<ID>_surface.npz` — LCFS `R,Z,nfp` grid (input to geometry).
- `spf_prototype/data/equil_precise_qa_*` , `precise_QA_fluxmap.npz` — DESC-example anchors.

**How to obtain a boundary + fluxmap for a device:** `quasr_fluxmap.build(ID)` (DESC venv)
is the single entry point — it fetches the boundary via `quasr_equilibrium_field`, solves,
and writes both the `.npz` and the `spf_fluxmap_v1` files the source consumes.

---

## 3. GEOMETRY

| file | key symbols | what it produces | run command |
|---|---|---|---|
| `spf_prototype/python/stellarator_geometry.py` | `DEFAULT_LAYERS`, `load_surface(stem)`, `poloidal_outward_normals(R,Z)`, `offset_surface(R,Z,nR,nZ,delta)`, `build_layers(stem, layers=DEFAULT_LAYERS, scale=1.0, outdir=None)`, `write_stl_shell`, `_write_binary_stl`, `is_edge_manifold`, `signed_volume`, `cross_sections_simple` | Pure-numpy conformal radial build: offsets the LCFS outward along the poloidal normal into nested **watertight shell-solid STLs**, one per layer. Emits `manifest.json` with per-layer watertight/simple/volume checks. `DATADIR` is module-level (rebound by callers). | `python python/stellarator_geometry.py <stem> [scale]` -> `data/<stem>_geom/*.stl` |
| `spf_prototype/python/build_dagmc.py` | `build(geom_dir, out_h5m="stellarator.h5m")` | Ginsburg/conda: converts the per-layer STL shells to one watertight DAGMC `.h5m` via `stl_to_h5m`, tagging each volume with the LAYER name. Hard-gates on watertight+simple. SOL/plasma interior deliberately untagged -> DAGMC implicit complement -> vacuum. | `python python/build_dagmc.py <geom_dir> [out.h5m]` (conda) |
| `spf_prototype/python/dagmc_writer.py` | `write_h5m`, `assemble`, `build_from_stls(geom_dir, out)`, `build_from_geometry(stem, geom_dir, out)`, `verify_structure`, `_read_stl_tris`, `_dedup_surface`, `_orient_outward`, `selftest`/`selftest_stl` | **pymoab-free** STL->DAGMC `.h5m` writer (raw HDF5). The alternative to `build_dagmc.py`/`stl_to_h5m`; used by `conformal_wallmap.build_wall`. Welds coincident nested-shell surfaces. | `python python/dagmc_writer.py ...` |
| `spf_prototype/python/coil_geometry.py` | `rmf_frames(P)` (rotation-minimizing frame + holonomy), `cross_section(radius,section,m)`, `build_tube(P,r,s,radius,section,m)`, `sweep_device(ID, radius=None, radius_frac=0.06, scale=1.0, section='circle', m=16, n_samples=256)`, `tangents`,`curvature` | **NEW.** Sweeps a circular/square winding-pack cross-section along REAL QUASR coil filaments (RMF double-reflection frame, no Frenet twist) into a homogenized coil solid STL. Diagnostics: `holonomy_deg`, `max_pinch` (curvature*radius; >1 = self-intersection). | `python python/coil_geometry.py <ID> [--section circle|square] [--radius R] [--scale S] [--m 16] [--preview]` -> `data/quasr/coil_geom/<ID>_<section>.stl` |
| `spf_prototype/sweep/coil_standoff.py` | `coil_standoff(ID, n_samples=256, ntheta=64, nphi=192)`, `_min_point_to_set`, `main()` | Min coil-filament-to-LCFS distance per device (reactor-relevance filter): `d_min`, `d_min/a`, `d_min/R0`, per-coil closest approach, absolute `gap_cm` at a target reactor R0. | `python sweep/coil_standoff.py [ids...] [--reactor_R0_cm 1000] [--blanket_cm 101]` -> `data/quasr/coil_standoff.csv` |

**Material-tag convention:** each STL shell layer STL is named `<layer>.stl`; the DAGMC
volume group tag == the layer name (W/steel/Be/FLiBe/shield/coil). At transport time OpenMC
matches a volume's `mat:NAME` to the `openmc.Material` whose `.name == NAME` — so every run
script forces `Material.name = layer_key` (see `run_conformal.make_geometry`,
`conformal_scatter_59509.build`). Mismatch => geometry-load abort.
`DEFAULT_LAYERS` (cm, plasma->out): `sol 5.0(vac), W 0.2, steel 3.8, Be 2.0, FLiBe 50.0,
shield 40.0, coil 15.0`.

**Current coil representation caveat:** the conformal `coil` layer from `build_layers` is a
lumped uniform shell (one cell / one material). `coil_geometry.py` produces true swept coil
solids as STLs, but these are **not yet wired into a transport model** (no tally on them
beyond the lumped `coil` cell). See Gaps.

---

## 4. OPENMC MODELS + TALLIES

| file | model builder | tallies | run command / env |
|---|---|---|---|
| `spf_prototype/python/conformal_scatter_59509.py` | `build(abc, fluxstem, R0, a, particles=300000)`; `scaled_fluxmap()` (writes x10 cm fluxmap); `DAGMCUniverse(H5M).bounded_universe()` | `flux_mesh`: `CylindricalMesh(40 r × 32 φ × 40 z)` + `EnergyFilter([0.1e6,20e6])`, score `flux`. `heat`: `MaterialFilter([W,steel,Be,FLiBe,shield,coil])` score `heating`. | Ginsburg conda; native `StellaratorSource`, scattering ON. `python python/conformal_scatter_59509.py` -> `data/quasr59509_conformal_flux.npz` |
| `spf_prototype/python/conformal_wallmap.py` | `build_wall(ID)` (LCFS+0.30a shell -> `dagmc_writer.build_from_stls`); `run_mode(pol, fluxstem, h5m, cwd, sid)` | Free-streaming NWL: `surf_source_write` on the plasma-facing wall surface; crossings KD-tree'd to `(θ,φ)` -> `NWL(θ,φ)`. 3M particles/mode, near-void firstwall. | Ginsburg conda; `python python/conformal_wallmap.py <ID>` -> `data/quasr<ID>_conformalmap.npz` |
| `spf_prototype/python/run_conformal.py` | `make_geometry(h5m)`, `build_model(abc, h5m, fieldmap_stem, R0_cm, a_cm, particles=200000, batches=10, density_scale=1.0)`; runs free(dscale 1e-4) + scatter(1.0) × 3 modes | `wall_current_phi`: `MeshSurfaceFilter` (φ-resolved current). `cell_response`: `MaterialFilter` heating/damage-energy/H3-production. `coil_fast`: `MaterialFilter([coil])`+E>0.1MeV flux. `firstwall_dir`: `MeshFilter` fast flux (steering contrast). | Uses legacy box-torus `CompiledSource` (`br.SO`, `bmode=fieldmap`). `python python/run_conformal.py <h5m> <stem> <scale>` |
| `spf_prototype/sweep/geometry_peaking.py` | `analyze_device(ID, ...)`, `wall_load(...)`, `wall_mesh(...)`, `sample_source(...)`, `run_all(ids=ALL_IDS)` | **No OpenMC** — analytic free-streaming `Q(r_w)=Σ s(x)·incid·[1+a2·P2(n·b̂)]/r²`, peaking `max/area-mean`, inboard/outboard asymmetry, per-mode. b̂ from real coils (Biot-Savart). | Fully local. `python sweep/geometry_peaking.py` -> `sweep/sweep_out/geometry_peaking.json` |
| `spf_prototype/python/reactor_model.py` | `build_model(abc, li6_enrich=30, density_scale=1.0, particles=100000, ...)` — box-torus (ZCylinder+ZPlane), legacy `CompiledSource` `bmode=toroidal` | `wall_current` (MeshSurfaceFilter), `cell_response` (CellFilter heating/damage/H3), `flibe_spectrum`, `tbr_nuclide` (Li6/Li7). Helpers `tbr`, `midplane_intensities`, `wall_ratio_inboard_outboard`. | `python python/reactor_model.py` (self-test) |
| `spf_prototype/python/stellarator_model.py` | `build_model(abc, fieldmap_stem, li6_enrich=65, ...)` — Helios-class axisymmetric box-torus radial-build rig; `make_materials`; legacy `CompiledSource` `bmode=fieldmap` | `wall_current`, `cell_response` (CellFilter), `coil_fast` (CellFilter[coil]+E>0.1MeV), `tbr_poloidal`. | `python python/stellarator_model.py` (smoke) |

`DAGMCUniverse` usage pattern (all conformal models): `openmc.Geometry(openmc.DAGMCUniverse(h5m).bounded_universe())`;
`bounded_universe()` adds the vacuum bounding cell; the implicit complement (SOL + plasma +
exterior) is filled with a near-void `vacuum` material.

Legacy `CompiledSource` build helper: `spf_prototype/python/build_and_run.py` — `SO =
src/build/libpolarized_fusion_source.so`, `build_so(force=False)` (cmake against
`$OPENMC_PREFIX`/`$CONDA_PREFIX`), `XS = $OPENMC_CROSS_SECTIONS`. Parameter-string format:
`"a=..,b=..,c=..,bmode={toroidal|fieldmap},[fieldmap=STEM,]shape=plasma,R0=..,aminor=.."`.

---

## 5. MATERIALS

`make_materials` lives in two layers. Base stack — `spf_prototype/python/reactor_model.py:make_materials(li6_enrich=30.0, density_scale=1.0)`:

| key | `Material.name` | composition | density (g/cm³) |
|---|---|---|---|
| `W` | tungsten | `W` 1.0 | 19.3 |
| `steel` | FeCr_steel | Fe 0.91 / Cr 0.09 (wo) — Eurofer-like | 7.8 |
| `Be` | beryllium | `Be` 1.0 | 1.85 |
| `FLiBe` | FLiBe | Li:Be:F = 2:1:4 (ao), Li enriched to `li6_enrich`% Li-6 | 1.94 |

Extended stack — `spf_prototype/python/stellarator_model.py:make_materials(li6_enrich=65.0,
density_scale=1.0)` calls `rm.make_materials` then adds:

| key | `Material.name` | composition | density | note |
|---|---|---|---|---|
| `shield` | WC_shield | W 1.0 / C 1.0 (elemental) | 15.6 | INJECT(helios): real WC+B4C+steel stack TBD |
| `coil` | coil_block | Fe 0.70 / Cr 0.18 / Ni 0.12 (wo) | 8.0 | **GENERIC Fe-Cr-Ni PLACEHOLDER — no real HTS/ReBCO winding-pack or steel-casing split. Flagged INJECT(helios).** |

`conformal_scatter_59509.py` and `run_conformal.py` both consume `stellarator_model.make_materials`
(keys W/steel/Be/FLiBe/shield/coil). Density is scaled by `density_scale` (`<<1` gives the
free-streaming near-void anchor). No separate `reactor_model.py` material beyond the base 4.

---

## 6. RESULTS PIPELINE + DATA

**Frozen tables** (`spf_prototype/paper/data/`, read-only inputs to figures; provenance in
`paper/data/PROVENANCE.md`):

| file | contents |
|---|---|
| `steering_dilution.csv` | free-streaming vs full-transport steering per mode (Tier 5, box-torus, scattering ON) |
| `rate_steering.csv` | absolute rate/load/TBR per mode, η restored (Tier 7) |
| `bare_cylinder.csv` | analytic point-patch error vs shaping (single filament) |
| `quasr59509_validation.csv` | analytic↔OpenMC per-wall A/iso, B/iso (QUASR 59509, anarrima-gated) |
| `emission_moments.csv` | ⟨(u·b̂)²⟩ about the 3-D field |
| `predictor_n15.csv` | 15-device conformal-wall predictor table + Spearman ρ |
| `conformal_xcheck.csv`, `predictor_n15.csv` | conformal free-stream cross-check / predictor |

**Frozen transport maps** (`spf_prototype/paper/data/conformal_maps/` and `spf_prototype/data/`):
- `quasr<ID>_conformalmap.npz` (15 devices) — `(θ,φ)` NWL maps `unpol/perp/par` + wall-cell
  `dA` on the DAGMC conformal first wall (LCFS+0.30a). From `conformal_wallmap.py`, 3M
  crossings/mode, Ginsburg jobs 8885457→8885458.
- `quasr59509_conformal_flux.npz` — volumetric fast-flux mesh `flux_<mode>` (39,32,39) +
  per-layer `heat_<mode>`, scattering ON, device scaled ×10. From `conformal_scatter_59509.py`
  (job 8912742). FLiBe = 84% of ~1.1e7 eV/src heating; near-wall fast flux perp −10% / par +10% vs unpol.
- `quasr59509_conformal_xcheck.npz` — conformal free-streaming cross-check.
- `data/quasr/coil_standoff.csv` — 366-device standoff filter table.
- `data/quasr/coil_geom/*.stl` — swept coil solids (`1035368_circle`, `45395_square`, `59509_circle`).

**Figure scripts** (`spf_prototype/paper/src/fig_*.py`, each reads only frozen files, none
recompute transport; `python fig_*.py` -> `../figures/*.pdf`): `fig_conformal_flux.py`
(NWL/dilution flux map), `fig_conformal_map.py`, `fig_conformal_xcheck.py`, `fig_dilution.py`,
`fig_predictor.py`, `fig_rate_steering.py`, `fig_quasr_validation.py`, `fig_bare_cylinder.py`;
shared style `paper/src/_style.py`.

---

## Gaps for the shield-optimizer project (what does NOT exist yet)

1. **ParaStell integration — absent.** All geometry is home-grown: conformal shells
   (`stellarator_geometry.build_layers`) + a bespoke pymoab-free `dagmc_writer`, and a
   separate coil sweep (`coil_geometry.py`). There is no ParaStell `Stellarator` /
   `build_cubit_model` / `export_cad_to_dagmc` path, and the shield-opt
   `PUBLIC_DATA_MANIFEST.yaml` names ParaStell as the intended primary engine. This is
   greenfield.

2. **Variable-thickness radial build — absent.** `DEFAULT_LAYERS` / `LAYERS` are **uniform
   scalar thicknesses per layer** offset by a constant normal distance. There is no
   `thickness_matrix` `t_l(φ,θ)` (per-poloidal/toroidal-angle thickness field), which is the
   core degree of freedom the nonuniform-shield optimizer needs. `offset_surface` takes a
   single scalar `delta`.

3. **Unstructured/tet coil tallies — absent.** The coil is a **lumped cell/material number
   only** (uniform `coil` shell in `build_layers`; `coil_fast` = a single `MaterialFilter`/
   `CellFilter` flux over that one region). `coil_geometry.py` can produce real swept coil
   solids as STLs, but they are not meshed into an unstructured/tet mesh, not converted to
   DAGMC in a transport model, and carry no per-coil or spatially-resolved damage tally.

4. **Real HTS + casing coil materials — absent.** `coil_block` is a generic Fe-Cr-Ni
   placeholder (density 8.0, flagged INJECT(helios)); `shield` is bare WC. There is no ReBCO/
   Cu/steel winding-pack composition and no HTS-bulk vs steel-casing split, so magnet-damage
   / lifetime numbers are not yet physically meaningful.

5. **Closed-loop optimizer — absent.** Everything is a one-shot forward pipeline
   (equilibrium → source → geometry → transport → frozen maps). `geometry_peaking.py` is an
   analytic *predictor* study, not an optimizer; `optimize_polarization.py` tunes the a2/abc
   *emission*, not the *shield geometry*. There is no objective over hotspot NWL, no
   thickness-field parameterization, no finite-diff/surrogate loop over ParaStell+OpenMC.
   The manifest explicitly calls the closed loop "our novelty" (future work).
