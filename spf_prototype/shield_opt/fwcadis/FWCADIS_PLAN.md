# FW-CADIS Weight-Window Generation via Random Ray — Recipe & Job Plan

Target: deep-penetration **reciprocity** run. Emit neutrons FROM a coil (DAGMC
cell, box source constrained to that cell, isotropic 14.06 MeV), tally flux on a
`RegularMesh` over the plasma region to get the plasma-region importance map
(which plasma regions can reach the coil through ~1 m of shield+blanket). Flux
drops ~6 orders of magnitude, so it needs variance reduction. This replaces the
failed MAGIC attempt (particle-splitting explosion, ~1 h/batch).

Anchored to the OpenMC checkout at `/Users/tkiker/Documents/GitHub/openmc`:
**version 0.15.4.dev182+g608a1c338 (tag 0.15.3 + 182 commits, branch
spf-prototype)**. The Ginsburg fork is 0.15.1-dev; the FW-CADIS + random-ray API
below was already present at 0.15.1/0.15.2 (see release notes 0.15.1.rst /
0.15.2.rst), but the *convenience* methods `Model.convert_to_multigroup()` and
`Model.convert_to_random_ray()` may be newer — **confirm they exist on the
Ginsburg fork** (see "Version risk" at the end). Everything is verified against
the local package source, not from memory.

---

## 1. Does OpenMC expose FW-CADIS directly? — CONFIRMED YES, with a required random-ray solve

FW-CADIS is **not** a one-line switch on a CE Monte Carlo deck. It is exposed
through `openmc.WeightWindowGenerator(method='fw_cadis', ...)`, but that
generator **only functions when the model is run in random-ray mode**
(`settings.random_ray` populated, multigroup energy mode). Confirmed in the
source and docs:

- `openmc/weight_windows.py` L487–636: `WeightWindowGenerator` accepts
  `method` in `{'magic','fw_cadis'}` (checked at L613–614) and a `targets`
  attribute (`openmc.Tallies` or iterable of tally IDs) for **local** variance
  reduction (L622–636).
- `docs/source/usersguide/variance_reduction.rst` L77–207 ("Generating Global or
  Local Weight Windows with FW-CADIS and Random Ray"):
  > "The only procedural difference [vs MAGIC] is that the code must be run in
  > random ray mode." … "OpenMC will automatically run a **forward solve
  > followed by an adjoint solve**, with a `weight_windows.h5` file generated at
  > the end."
- So: you do **not** hand-run separate forward/adjoint jobs and build WW
  yourself. You set up ONE random-ray fixed-source model with a `fw_cadis`
  `WeightWindowGenerator`; the solver internally does forward + adjoint and emits
  `weight_windows.h5`. The adjoint source is derived automatically from the
  `targets` tallies (FW-CADIS local VR) — this is exactly our reciprocity case:
  target = the plasma-region mesh tally.

Regression tests that pin the API (all in `tests/regression_tests/`):
- `weightwindows_fw_cadis/test.py` — global FW-CADIS: `WeightWindowGenerator(
  method="fw_cadis", mesh=ww_mesh, max_realizations=batches)` +
  `settings.random_ray['volume_estimator']='naive'`.
- `weightwindows_fw_cadis_local/test.py` — **local** FW-CADIS with
  `targets=model.tallies` (this is the pattern we use).
- `random_ray_adjoint_fixed_source/`, `random_ray_fixed_source_domain/` — pin the
  fixed-source random-ray source spec (see §3).

**Verdict: FW-CADIS-via-random-ray is feasible and directly supported in this
OpenMC version.** It is the intended replacement for MAGIC on deep-penetration
problems: the adjoint (importance) solve is deterministic-flavored (random ray),
so there is no MAGIC-style iterative MC splitting explosion during generation.

---

## 2. Multigroup cross sections (random ray is MG-only)

Random ray cannot use CE data; every material needs MG macroscopic XS in an
`mgxs.h5` library. Two paths:

### 2a. Recommended: auto-generate with `Model.convert_to_multigroup()` — CONFIRMED in-repo
`openmc/model/model.py` L2687–2810. One call runs CE MC, tallies material-wise
MGXS, writes `mgxs.h5`, and rewrites every material to `add_macroscopic(name)` /
`set_density('macro',1.0)` / `energy_mode='multi-group'`:

```python
model.convert_to_multigroup(
    method="material_wise",   # default; runs the real geometry+source, best fidelity
    groups="CASMO-2",         # or a custom openmc.mgxs.EnergyGroups / edge list
    nparticles=... ,          # per batch of the CE MGXS-gen run; raise for fidelity
    mgxs_path="mgxs.h5",
    overwrite_mgxs_library=False,  # reuse if present
)
```

**DAGMC is handled**: L2758–2773 detect `DAGMCUniverse`, init the lib in
`volume` mode, call `sync_dagmc_universes()` to attach materials to DAGMC cells,
finalize, restore run mode — so `convert_to_multigroup` works on the DAGMC coil
model. This is the single most important confirmation for us.

Method trade-offs (`docs/source/usersguide/random_ray.rst`, "Generating MGXS"):
- `material_wise` (default, best): runs the **full geometry+source**. Caveat that
  matters for us — *"If a material is only present far from the source and
  doesn't get tallied to in the CE simulation, the MGXS will be zero for that
  material."* In the reciprocity model the source is IN the coil, so the coil,
  magnets, vac_vessel, shield all get tallied; but a material buried where almost
  no CE particle reaches could come out with zero/garbage XS. **Mitigation:**
  either use a lot of `nparticles`, or fall back to `stochastic_slab`.
- `stochastic_slab` (medium, robust): simplified 1D geometry, guarantees non-zero
  XS for every material regardless of distance from source. **This is the safer
  choice for a deep-penetration deck** where far materials won't be reached by an
  analog CE gen source. Recommended fallback / primary if material_wise gives any
  zero-XS material.
- `infinite_medium` (low): one infinite-medium run per material; may hang if a
  material has k_inf>1 (not a risk for these non-fissile shield materials, but
  lowest fidelity).

For fusion 14 MeV, `CASMO-2` (2-group) is coarse. Prefer a finer fast-resolved
structure — e.g. a CCFE/VITAMIN-J-style or at least a several-group custom
structure with fast bins above 0.1 MeV. `openmc.mgxs.GROUP_STRUCTURES` lists
built-ins; pass a name or an explicit eV edge list. Start with a modest custom
structure (e.g. ~7–11 groups with fast detail) and refine if the WW look wrong.
For MGXS energy sampling in stochastic_slab/infinite_medium, pass
`source_energy=openmc.stats.delta_function(14.06e6)` so the collapse spectrum is
the DT source (docs "Generating MGXS", note on `source_energy`).

### 2b. Existing public MG library instead of generating?
Possible in principle (`materials.cross_sections = <some mgxs.h5>` +
`add_macroscopic`), but public MPACT/CASMO libraries are **LWR/thermal-spectrum**
libraries keyed to specific reactor nuclide names and group structures; they do
**not** contain macroscopics for FLiBe/WC/Be/HTS-magnet compounds at a 14 MeV
fusion spectrum, and OpenMC MG macroscopics are looked up by **material name**.
**Do not rely on an off-the-shelf library** — generate the MGXS from our own CE
model with `convert_to_multigroup`. That is the honest, correct path here.

---

## 3. Random-ray fixed source setup for the reciprocity problem

Random-ray fixed source is a normal `openmc.IndependentSource` with a **domain
constraint** and MG energy (verified `random_ray_fixed_source_domain/test.py` and
`examples.py` L1271–1276):

```python
src = openmc.IndependentSource(
    energy=openmc.stats.Discrete([14.06e6], [1.0]),
    constraints={'domain_type': 'cell', 'domain_ids': [<coil_cell_id>]},
    strength=1.0)
settings.source = [src]
```

`convert_to_random_ray()` (`model.py` L2812–2870) auto-fills `ray_source`
(uniform Box over bbox), `distance_inactive` (=max(chord,30 cm)),
`distance_active` (=5×), and a starting `particles` guess. Then overlay a
mesh so the DAGMC domain is subdivided into finer flat source regions AND so the
WW mesh does not subdivide any source region (the doc warning at
variance_reduction.rst L150–155 — **use the same mesh for
`random_ray['source_region_meshes']` and for the WW generator**):

```python
model.settings.random_ray['source_region_meshes'] = [(ww_mesh, [root_universe])]
```

FW-CADIS generator + target (local VR = drive flux to the plasma mesh tally):

```python
wwg = openmc.WeightWindowGenerator(
    method='fw_cadis', mesh=ww_mesh,
    targets=[plasma_mesh_tally],          # the reciprocity deliverable tally
    max_realizations=settings.batches)
model.settings.weight_window_generators = wwg
model.settings.random_ray['volume_estimator'] = 'naive'   # per fw_cadis tests
```

**Gotcha (documented, enforced at export):** the `targets` tallies must also be
present in `model.tallies`, and you must export via `Model.export_to_model_xml`
/ `Model.run` (NOT standalone `Settings.export_to_xml`) or the target check is
skipped (variance_reduction.rst L191–201).

---

## 4. The full workflow (what `reciprocity_fwcadis.py` implements)

```
CE model (DAGMC + real materials + coil box source + plasma mesh tally)
   │
   ├── deepcopy ──► convert_to_multigroup(material_wise|stochastic_slab)  ► mgxs.h5
   │                convert_to_random_ray()
   │                random_ray['source_region_meshes'] = [(ww_mesh,[root])]
   │                WeightWindowGenerator(fw_cadis, ww_mesh, targets=[plasma tally])
   │                model.run()   ──►  forward RR solve + adjoint RR solve
   │                                    ──►  weight_windows.h5
   │
   └── CE model (unchanged) + weight_windows.h5 applied
                    settings.weight_windows = openmc.hdf5_to_wws('weight_windows.h5')
                    settings.weight_windows_on = True
                    model.run()  ──► reciprocity flux on plasma mesh (the deliverable)
                    ──► structure gate (non-empty, localized, structured)
```

Two distinct OpenMC executions: (1) MG random-ray WW generation, (2) CE
fixed-source reciprocity tally with WW applied. The MGXS gen is a sub-step of (1)
(an internal CE run driven by `convert_to_multigroup`).

---

## 5. Ginsburg job sequence, runtimes, failure points

| Step | Job | What runs | Rough wall time | Notes |
|---|---|---|---|---|
| 0 | (interactive) | import check on fork; confirm `convert_to_multigroup`/`convert_to_random_ray` exist | minutes | **hard gate** — see Version risk below |
| 1 | `mgxs_gen` | CE MC to collapse MGXS (inside `convert_to_multigroup`) | ~10–40 min, `nparticles`×batches; DAGMC init overhead | material_wise runs full geometry; far materials may get zero XS → use stochastic_slab |
| 2 | `rr_ww` | random-ray forward+adjoint FW-CADIS → `weight_windows.h5` | ~minutes–1 h; scales with rays × groups × FSRs, NOT with penetration depth | random ray is deterministic-flavored: **no MAGIC splitting blow-up** |
| 3 | `reciprocity` | CE fixed-source, WW applied, plasma-mesh flux tally | target FOM-dependent; deep penetration but WW-accelerated; expect far better than the ~1 h/batch MAGIC run | can be MG instead of CE for speed, but CE recommended for the physics deliverable |

Suggested: steps 1+2 in one Slurm job (they share the MG model in-process);
step 3 a separate job that only reads `weight_windows.h5`. Persist `mgxs.h5` and
`weight_windows.h5` between jobs. Use a per-job cwd (HDF5 race note from memory).

### Where this can fail (flag these honestly)
1. **Version risk (highest):** if the Ginsburg 0.15.1-dev fork predates
   `Model.convert_to_multigroup` / `convert_to_random_ray`, the "easy way" is
   unavailable and you must build the MGXS manually with `openmc.mgxs.Library`
   (the "Hard Way" in random_ray.rst) and hand-populate `settings.random_ray`
   (ray_source, distances, source_region_meshes). The `fw_cadis`
   `WeightWindowGenerator` itself is older and should be present. **Verify first.**
2. **MGXS for DAGMC:** `convert_to_multigroup` needs a DAGMC init/finalize in
   volume mode to sync materials (model.py L2758–2773). Requires the OpenMC build
   to have **DAGMC enabled** (`-DOPENMC_USE_DAGMC=ON`, MOAB available) — the
   memory note flags that a prior PR build had DAGMC OFF. If DAGMC is off, this
   step errors immediately.
3. **Zero-XS materials:** material_wise gives zero MGXS for any material the CE
   gen source never reaches. Deep-penetration geometry makes this likely for the
   outermost magnets. Detect (inspect `mgxs.h5` for all-zero total XS) and switch
   to `stochastic_slab`.
4. **Random-ray domain / FSR setup on DAGMC:** random ray needs finite source
   regions; DAGMC cells are the base FSRs, subdivided by
   `source_region_meshes`. If a mesh cell straddles vacuum + material or the
   outer sphere void, the flat-source approximation degrades. Keep the WW/SRM
   mesh inside the device envelope; give near-void materials a real (tiny) Sigma_t
   (the existing `Vacuum` mat uses 1e-10 — OK, but confirm RR tolerates it; RR has
   a `random_ray_void`/`random_ray_low_density` test, so voids are supported).
5. **WW mesh subdividing source regions:** must reuse the SAME mesh for
   `source_region_meshes` and the WW generator (doc warning). The script does.
6. **Target tally not registered:** FW-CADIS `targets` must be in `model.tallies`
   and exported via `Model.run`/`export_to_model_xml` (doc warning L191–201).
7. **MG source spatial extent:** the coil box source is constrained to a DAGMC
   cell; RR samples the physical source only where that cell overlaps a source
   region. Make sure the coil cell is covered by the SRM mesh.
8. **Group structure too coarse:** CASMO-2 is thermal-reactor-oriented; a 14 MeV
   shielding penetration wants fast-group detail. Wrong structure → WW that
   under/over-split at the shield. Use a fusion-appropriate multi-group structure.

---

## 6. Concrete API summary (verified against this checkout)

- `openmc.WeightWindowGenerator(mesh, energy_bounds=None, particle_type='neutron',
  method='magic', targets=None, max_realizations=1, update_interval=1,
  on_the_fly=True)` — `weight_windows.py` L539–561.
- `method='fw_cadis'` requires random-ray mode; produces `weight_windows.h5`.
- `settings.random_ray` keys: `distance_inactive`, `distance_active`,
  `ray_source`, `volume_estimator` ∈ {naive, simulation_averaged, hybrid},
  `source_shape` ∈ {flat, linear, linear_xy}, `volume_normalized_flux_tallies`,
  `adjoint`, `sample_method`, `source_region_meshes`,
  `diagonal_stabilization_rho`, `adjoint_source` — `settings.py` L193–246,
  L1406–1456. (FW-CADIS drives the adjoint solve automatically; you don't set
  `adjoint=True` yourself for WW gen.)
- `Model.convert_to_multigroup(method, groups, nparticles, overwrite_mgxs_library,
  mgxs_path, correction, source_energy, temperatures, temperature_settings)` —
  model.py L2687; DAGMC-aware.
- `Model.convert_to_random_ray()` — model.py L2812; auto-fills ray source +
  distances + particles.
- Apply WW later: `settings.weight_windows = openmc.hdf5_to_wws('weight_windows.h5')`;
  `settings.weight_windows_on = True` (pattern already used in
  `../ww_pert_pipeline/coil_run.py` L73–76).

## References
- OpenMC docs: `docs.openmc.org` → User's Guide → Variance Reduction (FW-CADIS +
  Random Ray) and Random Ray (MGXS generation, quick start). In-repo:
  `docs/source/usersguide/variance_reduction.rst`,
  `docs/source/usersguide/random_ray.rst`,
  `docs/source/methods/random_ray.rst`,
  `docs/source/methods/variance_reduction.rst`.
- Tramm et al., random ray method in OpenMC (Ann. Nucl. Energy); FW-CADIS:
  Wagner & Haghighat, "Automated variance reduction of Monte Carlo shielding
  calculations using the discrete ordinates adjoint function" (CADIS) and
  Wagner, Peplow, Mosher (FW-CADIS). Cite these for the forward-weighted adjoint
  importance methodology the random-ray solve reproduces.
- Reuse geometry/materials/source from `../ww_pert_pipeline/coil_run.py` and
  `gen_ww.py` (the MAGIC version).
