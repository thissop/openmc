# ParaStell + SIMSOPT Capability & Novelty Brief

Scope: feasibility of building a closed-loop coil-shield optimizer **on top of** ParaStell
(geometry) + SIMSOPT (smooth fields / optimization) so the novelty is the optimizer, not the CAD.

Legend: **[VERIFIED]** = quoted from repo/docs/paper. **[INFERENCE]** = my reading. **[UNVERIFIED]** =
could not confirm; flagged for you to check.

Bottom line up front:
- **ParaStell already supports a spatially-varying (poloidal + toroidal) per-layer thickness matrix** — exactly your `t_shield(theta,phi)` control field. This is the single most important find; you do **not** need to modify ParaStell's radial-build core.
- **ParaStell builds real magnet coil solids from filament data** (QUASR-compatible), and **has a Cubit-free DAGMC path** (`cad_to_dagmc` / PyDAGMC / MOAB). Cubit is optional. Adoption gate is clear.
- **SIMSOPT can drive an external black-box objective** (an OpenMC tally) via `make_optimizable`, and its `SurfaceRZFourier` Fourier(theta,phi) parameterization is directly reusable for a smooth `t_shield` field.
- **Novelty holds, but narrowly.** The closest group (Moreno/Bader/Wilson, UW-Madison — the ParaStell authors) is *actively* working on numerical optimization of stellarator blanket/shield systems (Moreno's PhD topic). Their **published** work is a manual 3-D parameter sweep + global-scalar surrogate formulas, **not** a closed MC-tally → spatial-thickness-field → rebuild loop. The specific closed-loop tool you describe appears unpublished — but frame carefully; this is their backyard.

---

## 1. PARASTELL (svalinn/parastell — Moreno, Bader, Wilson; UW-Madison CNERG + Type One Energy)

Repo: https://github.com/svalinn/parastell
Paper: Moreno, Bader, Wilson, "ParaStell: parametric modeling and neutronics support for stellarator
fusion power plants," *Front. Nucl. Eng.* 3 (2024). https://doi.org/10.3389/fnuen.2024.1384788

### 1a. What it builds — [VERIFIED]
Core capabilities (repo README):
- "Model **in-vessel components of uniform or non-uniform thickness** using plasma equilibrium VMEC
  data or custom first wall data, and a user-defined radial build."
- "Model **magnet coils** using coil filament point-locus data and a user-defined cross-section."
- "Generate tetrahedral meshes of in-vessel components and magnets."
- "Generate DAGMC geometries" + neutron wall-loading support.

Construction model (paper): in-vessel components are built by lofting outward from the plasma/first-wall
surface. First wall, breeder/blanket, back wall, shield, vacuum vessel are stacked layers of a
**radial build**.

### 1b. CRITICAL — spatially-varying radial build (poloidal + toroidal)? — [VERIFIED: YES]
This is exactly your `t_shield(theta,phi)`.

Paper (quote): *"In-vessel components are defined by providing a 3-D radial build that assembles a set of
1-D radial builds at each vertex in a user-defined grid of (phi, theta) locations."* Layer thickness is
`t_l(phi, theta)`; cumulative offset `o_l(phi, theta)` = running sum of layer thicknesses at each
(phi,theta); points are spline-lofted into solids. The paper's demonstration varies breeder+shield
thickness in 3-D from 30 cm up to 148.7 cm, constrained by plasma-to-coil space.

**Concrete API** (from `examples/parastell_cad_to_dagmc_example.py`) — [VERIFIED]:
```python
# angular grid on which the radial build is defined (degrees)
toroidal_angles = [0.0, 11.25, 22.5, 33.75, 45.0, 56.25, 67.5, 78.75, 90.0]   # phi
poloidal_angles = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0] # theta

# each layer gets a 'thickness_matrix' of shape (len(toroidal_angles), len(poloidal_angles)) = 9x9
radial_build = {
    "first_wall": {"thickness_matrix": uniform_unit_thickness * 5},
    "breeder":    {"thickness_matrix": <9x9 array, values ~25.0 .. 75.0>},  # spatially varying
    "back_wall":  {"thickness_matrix": uniform_unit_thickness * 5},
    "shield":     {"thickness_matrix": uniform_unit_thickness * 50},
    "vacuum_vessel": {"thickness_matrix": uniform_unit_thickness * 10},
}
```
Each layer's `thickness_matrix` is a 2-D array over the (phi, theta) grid, interpolated between vertices.
So **"locally thicken shield / remove FLiBe at fixed envelope" is a native operation**: you decrement the
`breeder` matrix and increment the `shield` matrix at the same (phi,theta) cells. Lengths in cm.

Note the paper's phrasing "1-D radial builds at each (phi,theta) vertex" = a per-layer thickness field
sampled on a coarse angular grid and lofted. Your control field just becomes: given SIMSOPT Fourier dofs,
evaluate `t_shield(phi,theta)` and `t_breeder(phi,theta)` on this grid, write the matrices, rebuild.
- [INFERENCE] The grid resolution (9x9 here) is user-chosen; finer grids = smoother field, more CAD cost.
- [UNVERIFIED] Exact dict key spelling/nesting may differ slightly by release (I read one example on `main`);
  confirm against your pinned version. Behavior (per-(phi,theta) matrix) is solid.

### 1c. Magnet coil geometry from filaments — [VERIFIED: YES]
Paper: *"For each coil filament, the point-locus data is connected via spline interpolation to create a
closed loop. A cross-section, defined by the user as either a parametric rectangle or a parametric circle,
is then swept along the filament spline to generate the magnet solid."*

API (example):
```python
stellarator.construct_magnets_from_filaments(
    coils_file, width=40.0, thickness=50.0, toroidal_extent=90.0, sample_mod=6
)
```
Input is a coil filament point-locus file (`coils.example`) + rectangular/circular cross-section.
- [INFERENCE] This is directly compatible with **your QUASR simsopt-serials coils**: QUASR gives filament
  point loci; you feed those as the coils file. This lets you replace the "uniform conformal magnet shell"
  with realistic QA/QH coils as planned. (You'll need to export QUASR/SIMSOPT curves to ParaStell's coil
  file format — a small adapter, not a physics problem.)

### 1d. DEPENDENCY GATE — Cubit required? — [VERIFIED: NO, Cubit is optional]
README dependency list:
- **Required:** CadQuery, PyDAGMC, MOAB, CAD-to-DAGMC (`cad_to_dagmc`), OpenMC, NumPy, SciPy, PyYAML.
- **Optional:** Coreform Cubit.

Two DAGMC paths exist (paper: *"ParaStell is able to automate the generation of DAGMC neutronics
geometries via Coreform Cubit or CAD-to-DAGMC"*):
- **Cubit-free (your path):** `stellarator.export_cad_to_dagmc(filename="dagmc", export_dir=...)`
  → produces a DAGMC `.h5m` via the open-source `cad_to_dagmc` + MOAB. Examples:
  `parastell_cad_to_dagmc_example.py`, `parastell_pydagmc_example.py`, `nwl_pydagmc_example.py`.
- Cubit path: `parastell_cubit_example.py` (needs a Coreform license on PYTHONPATH; the authors note Cubit
  seat count is the scaling bottleneck and they are moving *away* from it, incl. magnet CAD via CadQuery).

**Adoption verdict:** a fully open-source, Cubit-free CAD→`.h5m` pipeline exists and is a supported example.
This is the key gate and it passes. [INFERENCE] `cad_to_dagmc` faceting quality/robustness on thin,
highly-shaped shield layers is the practical risk to test early (imprinting/merging of coincident surfaces
between adjacent lofted layers). Budget time to validate watertightness of the `.h5m` before trusting tallies.

### 1e. OpenMC hand-off — [VERIFIED: YES]
ParaStell emits a DAGMC `.h5m` that loads directly as an OpenMC `DAGMCUniverse` for transport, and ships
NWL/tally examples (`nwl_pydagmc_example.py`, `custom_source_example.py`). So: ParaStell builds geometry →
DAGMC `.h5m` → OpenMC transport + coil-heating tally. That is your forward model, ready to wrap.

---

## 2. SIMSOPT (hiddenSymmetries/simsopt — Landreman et al.)

Docs: https://simsopt.readthedocs.io

### 2a. SurfaceRZFourier for a smooth angular control field — [VERIFIED representation; INFERENCE on reuse]
`SurfaceRZFourier` represents a toroidal surface as truncated Fourier series in (theta, phi) with optional
stellarator symmetry (the standard R(theta,phi), Z(theta,phi) VMEC convention). The dofs are the Fourier
amplitudes; SIMSOPT exposes them as named, individually fix/free-able Optimizable dofs with box constraints.

- [INFERENCE] Your scalar `t_shield(theta,phi)` can reuse the **same** double-Fourier basis (a scalar field
  is simpler than a surface's R,Z pair). Cleanest approach: define your own small `Optimizable` whose dofs
  are cos/sin Fourier amplitudes `t_mn`, enforce stellarator symmetry by keeping only the symmetric terms,
  evaluate `t_shield` on ParaStell's (phi,theta) grid. You get smoothness (bandlimited), a small dof count,
  and symmetry for free. You don't strictly need `SurfaceRZFourier` itself — you need its parameterization
  pattern, which is easy to mirror. Reusing the class directly is possible but it's built for R,Z surfaces,
  not arbitrary scalars.

### 2b. Optimization framework + external black-box objective — [VERIFIED: YES]
- `Optimizable` base class: manages dofs (fix/free, names, box bounds) and a dependency graph (re-evaluates
  only what changed). Global state vector auto-assembled for the solver.
- **Wrapping a non-simsopt forward model** — `make_optimizable`:
  ```python
  from simsopt import make_optimizable
  def myfunc(v):
      v.run()
      return v.wout.DMerc[-2]
  myopt = make_optimizable(myfunc, vmec)   # myopt.J() returns the value
  ```
  [INFERENCE] Swap `myfunc` for "write thickness matrices → ParaStell rebuild → OpenMC run → return coil
  hotspot tally." `make_optimizable(objective, tshield_field)` yields an Optimizable whose `.J()` is your
  OpenMC-derived scalar. This is the documented, intended mechanism for external/black-box objectives.
- **Solvers:** `LeastSquaresProblem.from_tuples([(fn, goal, weight), ...])` +
  `least_squares_serial_solve()` / `least_squares_mpi_solve()` (MPI-parallel **finite-difference** gradients;
  cost linear in #dofs). SIMSOPT explicitly supports both derivative-free and derivative-based optimization
  and integrates scipy minimizers.
- [VERIFIED caveat] Analytic/adjoint derivatives (`.dJ()`) exist **only** for native simsopt.geo/simsopt.field
  objects (coils, fields, geometric penalties). **An OpenMC tally is a black box → no analytic gradient.**
  You are limited to finite-difference or derivative-free (Nelder-Mead, etc.). Each FD gradient step costs
  (#dofs + 1) full ParaStell-rebuild + OpenMC transport runs — so **keep the Fourier dof count low** (a big
  reason to use a bandlimited field, not per-cell thicknesses). MC tally noise also corrupts FD gradients;
  plan for variance reduction / high history counts or a noise-tolerant / surrogate-assisted optimizer.

### 2c. Coil representation (relevant to QUASR serials) — [VERIFIED]
SIMSOPT models coils as `CurveXYZFourier` (Fourier curve in x,y,z) wrapped in `Coil` (curve + current),
with `BiotSavart` for the field. QUASR serials are simsopt objects, so you can load them, evaluate/plot,
and **export filament point loci to feed ParaStell's `construct_magnets_from_filaments`**. [INFERENCE] The
adapter is: sample each `CurveXYZFourier` at N points → write ParaStell coil file. Straightforward.

---

## 3. NOVELTY CHECK — closed-loop MC-tally → spatial shield field → rebuild, for a stellarator

Question: has anyone published a closed loop that maps Monte-Carlo (OpenMC/MCNP/Serpent) coil/component
dose or hotspot tallies back onto a **spatially-varying (poloidal + toroidal) shield-thickness field** for a
**stellarator**, trading breeder for shield at fixed outer envelope?

Findings:
- **ParaStell paper (Moreno et al. 2024) — closest, but NOT closed-loop.** [VERIFIED] It varies breeder+shield
  thickness in 3-D and tallies TBR + magnet heating in OpenMC, but the study **manually sampled 99
  configurations**; results were **not** fed back to iteratively update geometry. The paper explicitly states
  the loop is future work: *"coupling ParaStell to machine-driven optimization is planned in order to
  establish a FWBS neutronics optimization framework for stellarators."* → As of that publication the
  closed loop was unbuilt.
- **Connor Moreno's PhD dissertation topic** — [VERIFIED, and the main risk to your novelty claim]:
  *"investigating the application of numerical optimization methods to stellarator fusion power plant design,
  focusing on the modification of plant breeding blanket and shielding systems and the corresponding effects
  on neutronics performance."* This is the **same problem space**. Watch for his thesis / follow-on papers.
- **Published UW follow-on = surrogate formulas, not a spatial closed loop.** [VERIFIED] Their newer work
  derives *global-scalar power-law surrogates* (via sensitivity analysis + symbolic regression) relating TBR,
  nuclear heating, and dpa to **global** parameters (fusion power, plasma surface area, elongation) for
  systems-code use. That is the opposite abstraction level from a `t_shield(theta,phi)` field — global
  scalars, not a spatial control field, and not an iterative geometry-rebuild loop.
- **European stellarator (HELIAS/W7-X) power-plant studies** — [VERIFIED intent only]: reviews describe
  parametric meshed-CAD blanket models converted to MCNP/OpenMC and an *aim* to "close the optimization loop
  by providing feedback to the physics-based stellarator optimization framework." This is a stated goal /
  coupling to the *equilibrium* optimizer, not a published spatial breeder-for-shield trade at fixed envelope.
- **Tokamak variable-shield / adjoint dose work** exists in the broader literature [INFERENCE, not deeply
  scanned here], but it is axisymmetric (poloidal-only) and not the stellarator 3-D (poloidal+toroidal) field
  you're targeting.
- **AI-driven stellarator coil optimization with in-the-loop FE** (arXiv 2603.15240) — [VERIFIED] exists but
  is **structural/coil** optimization, not neutronics shield.

**Verdict:** The specific tool — closed MC-tally → smooth `t_shield(theta,phi)` field → ParaStell rebuild →
re-tally, trading FLiBe for shield at fixed envelope on a stellarator — **appears unpublished; the novelty
claim holds.** But it is *narrow* novelty: (i) the constituent geometry ability (3-D shield field, coils,
DAGMC) is exactly what ParaStell already provides, so novelty must live in the **optimizer + the
breeder-for-shield-at-fixed-envelope objective**, precisely as you intend; and (ii) the ParaStell authors have
publicly flagged this closed loop as their planned next step and it is a PhD student's dissertation topic —
so **move deliberately, cite Moreno et al. explicitly, and frame your contribution as the smooth-field
SIMSOPT-parameterized optimizer with the fixed-envelope FLiBe/shield trade**, not as "first neutronics shield
optimization for stellarators." Check for a 2025-2026 Moreno thesis / TOFE/ANS proceedings before claiming
priority.

---

## 4. Decision summary

| Requirement | Status | Evidence |
|---|---|---|
| Realistic QA/QH coils from QUASR | **Supported** | `construct_magnets_from_filaments` (filament loci + cross-section); QUASR serials are simsopt curves → export loci [VERIFIED + INFERENCE] |
| Spatial `t_shield(theta,phi)` field | **Native** | per-layer `thickness_matrix` over (phi,theta) grid, lofted [VERIFIED] |
| Trade breeder↔shield at fixed envelope | **Native** | decrement `breeder` matrix, increment `shield` matrix same cells [VERIFIED mechanism / INFERENCE on your use] |
| Cubit-free DAGMC `.h5m` | **Yes** | `export_cad_to_dagmc` via `cad_to_dagmc`+MOAB+PyDAGMC; Cubit optional [VERIFIED] |
| OpenMC transport hand-off | **Yes** | DAGMC `.h5m` → `DAGMCUniverse`; NWL/tally examples shipped [VERIFIED] |
| Smooth Fourier(theta,phi) field param | **Reuse pattern** | mirror `SurfaceRZFourier` double-Fourier + stellarator symmetry [INFERENCE] |
| SIMSOPT drives external OpenMC objective | **Yes** | `make_optimizable(myfunc, ...)` black-box wrap; scipy/least-squares solvers [VERIFIED] |
| Gradients for MC objective | **FD / derivative-free only** | no adjoint for black-box tally; FD cost ~(#dofs+1) runs/step, MC noise sensitive → keep dofs low [VERIFIED caveat] |
| Closed-loop stellarator spatial-shield optimizer published? | **No (novelty holds, narrowly)** | ParaStell paper: loop is "planned"; Moreno PhD = same space; published follow-on = global scalar surrogates [VERIFIED] |

### Risks / things to verify early
1. **`cad_to_dagmc` robustness** on thin, curved, adjacent shield/breeder layers — watertight `.h5m`,
   correct surface imprint/merge. Test before trusting tallies. [INFERENCE]
2. **FD-gradient cost × MC noise.** Each optimizer step = many full rebuild+transport runs. Bandlimit the
   Fourier field (few dofs), and/or use a noise-tolerant / surrogate-assisted / derivative-free optimizer.
   [VERIFIED caveat + INFERENCE]
3. **Exact ParaStell radial-build dict schema** on your pinned release (key names/nesting may drift). [UNVERIFIED]
4. **QUASR→ParaStell coil-file adapter** format. Small but real. [INFERENCE]
5. **Prior-art watch:** Moreno 2025-2026 thesis / ANS-TOFE proceedings could pre-empt the novelty. [VERIFIED risk]

## Sources
- ParaStell repo: https://github.com/svalinn/parastell (README, `examples/parastell_cad_to_dagmc_example.py`)
- Moreno, Bader, Wilson, *Front. Nucl. Eng.* 3 (2024): https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2024.1384788/full
- OSTI record: https://www.osti.gov/pages/biblio/2333714
- Connor Moreno (CNERG): https://cnerg.github.io/community/people/cam/ ; CNERG pubs: https://cnerg.github.io/pubs/journals.html
- UW-Madison feature (fusion startups use ParaStell): https://engineering.wisc.edu/blog/grad-student-develops-software-used-by-fusion-startups-to-automate-stellarator-design/
- SIMSOPT docs — Optimizable / defining problems: https://simsopt.readthedocs.io/en/latest/optimizable.html ; solve pkg: https://simsopt.readthedocs.io/en/latest/simsopt.solve.html
- European stellarator power-plant review: https://www.sciencedirect.com/science/article/pii/S0920379624002394
- AI-driven coil optimization (structural, not shield): arXiv:2603.15240
