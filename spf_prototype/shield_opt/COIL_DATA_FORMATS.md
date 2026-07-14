# Public stellarator coil datasets — format characterization

Purpose: characterize two public coil datasets so a **simsopt-free loader adapter**
(matching the pure-JSON style of `sweep/quasr_loader.py`) can be written. This doc
does NOT implement the adapter — it gives file paths, exact reconstruction recipes,
symmetry-expansion rules, currents, licenses, and sha256.

Environment note: `simsopt` **1.10.6 is importable** in this environment. Both
recipes below were cross-checked against simsopt, but each is documented so a
pure-numpy loader is possible.

All cached raw files live under:
`/Users/tkiker/Documents/GitHub/openmc/spf_prototype/shield_opt/data/raw/`

---

## DATASET 1 — Precise QA (Wechsung, Landreman, Giuliani, Cerfon, Stadler 2022)

- Repo: https://github.com/florianwechsung/CoilsForPreciseQS
- Paper: *Precise stellarator quasi-symmetry can be achieved with electromagnetic coils*
  (PNAS 119(13) 2022; arXiv:2108.03711). Targets the Landreman–Paul QA and QA+Well
  precise-QS configurations.
- **4 unique modular coils, order-16 Fourier, 16 total after symmetry (nfp=2, stellsym).**

### 1.1 Where the coil geometry actually lives

The repo itself contains **no coil files** — they are inside `archive.zip`, a
**git-LFS** object (the working-tree `archive.zip` in a normal clone is a 134-byte
LFS pointer). git-lfs is not installed here; the real object was fetched via the
GitHub media endpoint:

```
curl -L https://media.githubusercontent.com/media/florianwechsung/CoilsForPreciseQS/main/archive.zip
```

Cached at: `data/raw/CoilsForPreciseQS/archive.zip`
(317,990,178 bytes uncompressed store; sha256
`ce2a95762b930950e348b6dfb66b85c6ad88a02b588dc12dd92f4656fc6f8971`)

Inside `archive.zip`, coil geometry is under `output/<config>/xmin.txt`. The config
directory naming encodes the paper configuration:

- `well_False_...` = **QA**, `well_True_...` = **QA+Well**
- `lengthbound_{18,20,22,24}` = the `[18]/[20]/[22]/[24]` configs in the paper
- `ig_N` = which random initial guess; the paper's chosen (best) minimizer per config,
  from `eval_geo.py` / `eval_find_best.py`:
  - QA:      lb18→ig7, lb20→ig6, lb22→ig4, **lb24→ig2**
  - QA+Well: lb18→ig6, lb20→ig4, lb22→ig7, **lb24→ig5**

**Canonical files** (the four most useful, extracted to
`data/raw/CoilsForPreciseQS/qs_xmin/`):

| Config      | source dir (in archive.zip)                                              | file            | sha256 |
|-------------|--------------------------------------------------------------------------|-----------------|--------|
| QA[24]      | `output/well_False_lengthbound_24.0_..._ig_2_order_16_expquad/xmin.txt`   | `xmin_QA24.txt` | `29c7187109802b79bb8feadc05cc6c3b16b861bc82e5cc56778c1c0269fd358a` |
| QA+Well[24] | `output/well_True_lengthbound_24.0_..._ig_5_order_16_expquad/xmin.txt`    | `xmin_QAWell24.txt` | `5fc0090fc457b1e60190b93ac894fe377dad39f8ce5828577bf4e7559ad93760` |
| QA[18]      | `output/well_False_lengthbound_18.0_..._ig_7_order_16_expquad/xmin.txt`   | `xmin_QA18.txt` | `fa6e311041b19a667e32e55f8afb4033eb068a5186fe8da14ff1eea08b086f38` |
| QA+Well[18] | `output/well_True_lengthbound_18.0_..._ig_6_order_16_expquad/xmin.txt`    | `xmin_QAWell18.txt` | `41f071f886082afb481a806850eb0b928a18e7271d276fc1409feadaab7e9977` |

(Any of the 32 `output/*/xmin.txt` dirs are loadable the same way.)

### 1.2 File format: `xmin.txt`

Plain ASCII, **399 floats, one per line** (`np.loadtxt`). This is the simsopt
`BiotSavart(coils).x` free-DOF vector for an **order-16** configuration.
It is NOT a JSON/BiotSavart serialization and NOT a MAKEGRID file — it is the raw
optimizer state vector, meaningful only against the fixed coil-construction recipe below.

**DOF layout (verified empirically against simsopt 1.10.6):**

```
index   0..2    : 3 free current DOFs   -> Current2:x0, Current3:x0, Current4:x0
index   3..101  : CurveXYZFourier1  (99 dofs)
index 102..200  : CurveXYZFourier2  (99 dofs)
index 201..299  : CurveXYZFourier3  (99 dofs)
index 300..398  : CurveXYZFourier4  (99 dofs)
```
399 = 3 + 4×99. The **first coil's current is FIXED** (not in the vector).

Per-curve 99 dofs = three coordinate blocks of 33, in order **x-block, y-block, z-block**.
Each 33-dof block for order 16 is laid out as:
```
[ c0, s1, c1, s2, c2, ..., s16, c16 ]    (1 + 2*16 = 33)
```

**Pure-numpy centerline reconstruction (validated to 2.2e-15 vs simsopt):**
```python
order = 16 ; N = 160
t  = np.linspace(0, 1, N, endpoint=False)
th = 2*np.pi*t
def coord(block33, th):           # block33 = [c0, s1,c1, s2,c2, ...]
    v = block33[0]*np.ones_like(th)
    for m in range(1, order+1):
        s = block33[1 + 2*(m-1)]; c = block33[2 + 2*(m-1)]
        v += c*np.cos(m*th) + s*np.sin(m*th)
    return v
# base curve k (k=0..3): dofs = x[3 + k*99 : 3 + (k+1)*99]
X = coord(dofs[0:33], th); Y = coord(dofs[33:66], th); Z = coord(dofs[66:99], th)
```
This gives the **4 BASE curves in meters**. simsopt uses `numquadpoints = 10*order = 160`;
N is arbitrary for your own resampling.

### 1.3 Symmetry expansion (4 base -> 16 physical)

The 16 physical coils are `coils_via_symmetries(base_curves, base_currents, nfp=2, stellsym=True)`.
The MAKEGRID-style explicit list is NOT stored — you must expand yourself. Two options:

- **Simplest (recommended here, since simsopt imports):** reproduce the exact construction and
  set `bs.x`:
  ```python
  from simsopt.geo import create_equally_spaced_curves
  from simsopt.field import Current, coils_via_symmetries, BiotSavart
  from simsopt.field.coil import ScaledCurrent
  base = create_equally_spaced_curves(4, 2, stellsym=True, R0=1.1, R1=0.6,
                                       order=16, numquadpoints=160)
  bc = []
  for i in range(4):
      c = Current(1.)
      if i == 0: c.fix_all()
      bc.append(ScaledCurrent(c, 1e5))
  coils = coils_via_symmetries(base, bc, 2, True)   # -> 16 coils
  bs = BiotSavart(coils); bs.x = np.loadtxt("xmin_QA24.txt")
  centerlines = np.array([c.curve.gamma() for c in coils])   # (16,160,3) meters
  currents    = np.array([c.current.get_value() for c in coils])
  ```
- **Pure-numpy:** after reconstructing the 4 base gammas, apply, for p in 0..nfp-1 and
  flip in {+,−}: rotation by angle `2*pi*p/nfp` about z; the stellsym partner maps
  `(x,y,z)->(x,-y,-z)` with a `phi -> -phi` traversal reversal and **current sign flip**.
  (Replicating `coils_via_symmetries` by hand is error-prone; prefer the simsopt path
  unless a hard simsopt-free requirement exists.)

### 1.4 Currents

Physical current = (current DOF) × 1e5 A. Coil-0 DOF is fixed at 1.0 → **1e5 A**.
For **QA[24]** the 16 physical currents (A) are, in coil order:
```
+1.000e5, +9.8287e4, +1.01677e5, +1.03531e5,   (field period 0)
-1.000e5, -9.8287e4, -1.01677e5, -1.03531e5,   (stellsym partners, sign-flipped)
+1.000e5, +9.8287e4, +1.01677e5, +1.03531e5,   (field period 1)
-1.000e5, -9.8287e4, -1.01677e5, -1.03531e5
```
i.e. 4 distinct magnitudes {1.0, 0.98287, 1.01677, 1.03531}×1e5 A, signs flipped by stellsym.

### 1.5 Scale / units

Meters. `create_equally_spaced_curves` uses R0=1.1, R1=0.6 → coil centerlines span
**R ≈ 0.24–2.23 m, Z ≈ ±1.09 m** (experiment/unit scale, matching the
`input.LandremanPaul2021_QA` equilibrium, aspect ratio ≈ 6). No rescaling applied.

### 1.6 License

**None.** The repo has no LICENSE file (GitHub license API returns null; raw LICENSE = 404).
Usage terms are unspecified — treat as "all rights reserved" and cite the PNAS paper /
contact the authors before redistribution. **FLAG for the user.**

---

## DATASET 2 — Precise QH (Wiedman, Buller, Landreman 2023)

- Zenodo: https://doi.org/10.5281/zenodo.10211349 (record 10211349)
- Title: *Data and scripts for "Coil Optimization for Quasi-helically Symmetric
  Stellarator Configurations"*; creators **Wiedman, Alexander V.; Buller, Stefan;
  Landreman, Matt**; published 2023-11-27. (Author attribution matches the task.)
- Download API used:
  `https://zenodo.org/api/records/10211349/files/2023-11-24_QH_coil_optimization_paper.zip/content`

Cached zip: `data/raw/QH_coil_optimization_paper.zip`
(84,194,060 bytes; sha256
`28f03e503a83e92ab41783dc05114f3a67b6548d7c3eb1eda0816012e6106929`),
unzipped to `data/raw/QH_unzipped/zenodo/`.

Single record = one 84 MB zip → `zenodo/configurations/{LandremanPaulQH_coils,
LandremanPaulQH_scaled, LBDQH5, LBDQH5_coils, Mercier, Mercier_coils}` plus
`optimization_scripts/`, `plot_scripts/`, `perturbation_script/`.

### 2.1 Canonical coil file (Landreman–Paul precise-QH VACUUM set)

`data/raw/QH_unzipped/zenodo/configurations/LandremanPaulQH_coils/coils.curves_22_7_21`
(203,947 bytes; sha256
`3897fd5e3460d37118e858e30a1e2dd4cb0db9cbc83cf95bed90ae89484cc95e`)

Companion `extcur.curves_22_7_21` = MAKEGRID `EXTCUR(i)` list (20 values, all equal).
(`LBDQH5_coils/` and `Mercier_coils/` are finite-beta configs; `LandremanPaulQH_coils`
is the vacuum precise-QH set.)

### 2.2 File format: MAKEGRID `coils.` ASCII

NOT a simsopt JSON and NOT Fourier coefficients — it is a **MAKEGRID coils file** listing
every physical coil as an explicit **Cartesian polyline**. Header:
```
periods   4              <- nfp = 4
begin filament
mirror NIL
```
Then rows `X  Y  Z  CURRENT`. Each coil is a run of rows; its **last row is a terminator**
carrying `CURRENT = 0.0` plus two trailing tokens `<group_index> <label>`
(e.g. `... 0.000E+00 20 RotatedCurve50`). File ends with `end`.

- **20 coils**, **106 rows each** (105 distinct points + 1 terminator row that is an
  **exact duplicate of the first point**, closing the loop; verified closing-distance = 0).
- Coil labels: `CurveXYZFourier1..5` (the 5 unique base coils) + `RotatedCurve36..50`
  (15 rotated copies) = **5 unique × nfp(4) = 20 total**.

### 2.3 Reconstruction — trivial, NO Fourier, NO expansion

Points are already Cartesian and **all 20 coils are explicit**, so the loader needs
**no Fourier evaluation and no symmetry expansion**:
```python
lines = open(path).read().splitlines()[3:]     # skip 3-line header
coils, cur = [], []
for ln in lines:
    p = ln.split()
    if p[0] == 'end': break
    cur.append([float(p[0]), float(p[1]), float(p[2])])   # X,Y,Z  (meters)
    if len(p) > 4:                # terminator row (has group+label, current==0)
        coils.append(np.array(cur[:-1]))   # drop the duplicate closing point
        cur = []
# coils -> list of 20 arrays, each (105,3). Stack for (20,105,3).
```
Current per coil = column 4 of any non-terminator row (constant within a coil).
Verified with `simsopt.field.load_coils_from_makegrid_file(path, order=16, ppp=1)`:
20 coils, gamma resampled to (20,16,3), same geometry.

### 2.4 Currents

Column 4. **All 20 coils carry the same current 1.2779754812e7 A** (≈12.78 MA);
`extcur.curves_22_7_21` lists 20 identical `EXTCUR(i) = 1.27797548115612E+07`.

### 2.5 Symmetry / count

`periods 4` → **nfp = 4**. The file already contains all physical coils (5 unique + 15
rotated = 20); **no expansion required**. NOTE: this is **20 coils, not 40** — the file
applies only the 4× field-period rotation, not an additional stellarator-symmetry
doubling. If the downstream model assumes "40 coils / 5 per half-period," reconcile
against this: the vacuum set as distributed is 20 physical coils. **FLAG for the user.**

### 2.6 Scale / units

Meters (MAKEGRID convention). Centerlines span **R ≈ 5.13–22.65 m, Z ≈ −8.86..+4.09 m**
→ reactor-scale (major radius ≈ 14 m). This is a **scaled-up** configuration; a companion
`LandremanPaulQH_scaled/` directory exists (VMEC-only, no coils file). Rescale if a
different device size is needed.

### 2.7 License

**CC-BY-4.0** (Zenodo metadata `license.id = cc-by-4.0`). Attribution required.

---

## Summary table

| | Dataset 1 (QA) | Dataset 2 (QH) |
|---|---|---|
| Canonical file | `qs_xmin/xmin_QA24.txt` (+3 more) | `LandremanPaulQH_coils/coils.curves_22_7_21` |
| Format | simsopt `BiotSavart.x` DOF vector (399 floats, order 16) | MAKEGRID `coils.` Cartesian polylines |
| Fourier? | Yes — CurveXYZFourier `[c0,s1,c1,...,s16,c16]` per x/y/z block | No — explicit points |
| Reconstruction verified | pure-numpy vs simsopt = 2.2e-15 | simsopt makegrid loader ✓ |
| Unique coils stored | 4 base (need symmetry expansion) | 20 (already fully expanded) |
| Symmetry | nfp=2, stellsym → 16 physical | nfp=4, file lists all 20 (no expansion) |
| Currents | dof×1e5 A; coil0=1e5 A fixed; ±{1,0.98287,1.01677,1.03531}×1e5 | all 20 = 1.2779754812e7 A |
| Units / scale | meters, R≈0.24–2.23 (experiment scale) | meters, R≈5.1–22.6 (reactor scale) |
| License | **NONE (no LICENSE in repo)** | **CC-BY-4.0** |
| sha256 (canonical) | QA24 `29c71871…9fd358a` | `3897fd5e…84cc95e` |

## Flags / ambiguities
1. **QA repo has no license** — redistribution rights unspecified.
2. **QH file is 20 coils, not 40** — only field-period (×4) expansion, no stellsym
   doubling in the distributed vacuum set. Reconcile with any "40-coil" assumption.
3. **QA geometry is not stored as an explicit coil list** — only the optimizer state
   vector `xmin.txt`. Loading requires the exact `create_curves` recipe (R0=1.1, R1=0.6,
   order=16, nfp=2, ncoils=4, coil-0 current fixed). A pure-numpy loader must also
   replicate `coils_via_symmetries` for the 4→16 expansion; using simsopt is safer.
4. Both datasets are at their native scales (QA ~1 m, QH ~14 m). Rescale as needed.
5. QA `objective.py` in the repo imports the **old** simsopt API
   (`simsopt._core.graph_optimizable`) and will not import under simsopt 1.10.6; the
   construction was re-implemented with the current API (shown in §1.3).
