# Precise-QA device standup plan

Goal: stand up a **precise quasi-axisymmetric (precise-QA)** stellarator device to sit
alongside the existing precise-QH device, feeding the same pipeline:

```
equilibrium (wout / DESC)  ->  spf_fluxmap_v1  ->  native StellaratorSource
                           \->  ParaStell thickness_matrix build  ->  OpenMC DAGMC transport
coils  ->  coil_adapter.load_qa("QA24")  ->  MAKEGRID .coils  ->  ParaStell magnets
```

All facts below were verified against this checkout on 2026-07-22.

---

## 0. Bottom line (decisions)

| Question | Answer (verified) |
|---|---|
| **Which QA boundary** | **Landreman & Paul 2021 precise-QA** (PRL 128, 035001), nfp=2, aspect=6, iota~0.42. VMEC INDATA already vendored: `shield_opt/data/raw/equilibria/input.LandremanPaul2021_QA` (NFP=2, MPOL=16, NTOR=12, LASYM=F, PHIEDGE=0.0838573, unit scale R0~1.01 m / a~0.168 m / \<B\>~1 T). |
| **Which coils** | **Wechsung_QA24** (`shield_opt/data/processed/Wechsung_QA24.npz`, 16 coils, nfp=2). Per `CoilsForPreciseQS/README.md`, these modular coils were optimized to **approximate this exact Landreman-Paul QA boundary** (QA24 = length-bound-24, order-16). **Matched pair — no mismatch.** They are NOT the QA+Well variant (that's `xmin_QAWell*`). |
| **Reactor scale** | Native QA a=0.1683 m; reactor target a=1.7 m (same as QH) => **scale = 1.7/0.1683 = 10.10**, giving R0~10.2 m, a~1.7 m. `coil_adapter.load_qa` already applies this factor to the coils by default, so the plasma fluxmap must use the same 10.10 to stay in one frame. |
| **Do we already have a fluxmap?** | Yes but unit-scale + DESC-made: `data/precise_QA_fluxmap.npz` (nfp=2, R0~0.98 m, \|B\|~1 T), built via `quasr_fluxmap.build(example="precise_QA")` using `desc.examples.get("precise_QA")`. Good enough for the analytic anarrima benchmark; **not** reactor-scaled for the ParaStell/OpenMC build. |
| **New artifact needed** | `build_qa_fluxmap.py` (in this dir): turns a VMEC `wout_*.nc` into spf_fluxmap_v1 with an optional reactor scale. Smoke-tested on a real wout (gates pass). |

---

## 1. Getting the precise-QA equilibrium

The boundary is **already in the repo** — you do NOT need to fetch anything from SIMSOPT.
(Verified: the pip `simsopt` 1.10.6 package ships only `mhd/input.default`; the
`input.LandremanPaul2021_QA` file lives in SIMSOPT's source `tests/test_files`, and a
copy is already vendored at `shield_opt/data/raw/equilibria/input.LandremanPaul2021_QA`.)

Two independent routes to an equilibrium/fluxmap. **Route A (VMEC) is canonical** because
it matches the QH pipeline (`vmec_fluxmap.py`) and produces an SI wout that ParaStell and
BOOZXFORM can also consume. **Route B (DESC) is already done** and is the analytic-benchmark input.

### Route A — VMEC (canonical; matches QH pipeline)

**A1. Run VMEC on the vendored INDATA.** VMEC is not installed locally (no `vmecpp`, no
`vmec_venv`); this step needs the **Ginsburg cluster** (or any box with VMEC/VMEC++).

Option A1a — VMEC++ (Python, same engine `vmec_fluxmap.py` uses):
```bash
# on a machine/venv with vmecpp:
python - <<'PY'
import vmecpp
inp = vmecpp.VmecInput.from_file(
    "shield_opt/data/raw/equilibria/input.LandremanPaul2021_QA")
out = vmecpp.run(inp, verbose=True)          # fixed-boundary, LFREEB=F already
out.wout.save("wout_LandremanPaul2021_QA.nc")  # or out.save(...) per vmecpp version
PY
```

Option A1b — classic Fortran VMEC (STELLOPT/PARVMEC on Ginsburg):
```bash
# INDATA already has NS_ARRAY up to 201 and FTOL 2e-17; runs as-is.
srun -n 8 xvmec2000 input.LandremanPaul2021_QA
#   -> wout_LandremanPaul2021_QA.nc
```
Either produces `wout_LandremanPaul2021_QA.nc` (SI: meters, Tesla). Sanity after the run:
`nfp==2`, `aspect~6.0`, `Aminor_p~0.168`, `Rmajor_p~1.01`.

**A2. Synthesize the fluxmap (LOCAL — only numpy + a netCDF reader needed):**
```bash
cd shield_opt/qa_setup

# reactor-scaled map for the ParaStell build + OpenMC transport:
python build_qa_fluxmap.py wout_LandremanPaul2021_QA.nc \
       --scale 10.10 --out ../../data/precise_QA_reactor_fluxmap

# native unit-scale map for the analytic anarrima cross-check:
python build_qa_fluxmap.py wout_LandremanPaul2021_QA.nc \
       --scale 1.0   --out ../../data/precise_QA_vmec_fluxmap
```
Emits `<out>.npz` + `<out>.meta` + `<out>.bin` (spf_fluxmap_v1). Same synthesis and the
same two rigor gates as `vmec_fluxmap.py`: **V2** single-signed sqrt(g) (hard assert) and
**V1b** divergence-theorem volume vs VMEC volume (<1.5e-2). Verified end-to-end on the
LandremanPaulQH wout: relV=3.4e-3, gates pass.

### Route B — DESC (already produced the current unit-scale map)

`data/precise_QA_fluxmap.npz` was built by `python/quasr_fluxmap.py build(example="precise_QA")`
(`desc.examples.get("precise_QA")` -> a converged Landreman-Paul QA DESC equilibrium ->
same `write_fluxmap_bin`). DESC is not installed locally either (`import desc` fails), so
regenerating it also needs the cluster / a DESC venv. Keep this as the **analytic-benchmark**
input (it is what `python/ana_precise_qa.py` + anarrima consume in the native frame). For
a reactor build, prefer Route A so plasma and coils share the VMEC-derived geometry.

---

## 2. Coils (already done, just confirm the frame)

```bash
python shield_opt/coil_adapter.py qa --qa_config QA24 --save
#   -> data/processed/Wechsung_QA24.npz + .coils (MAKEGRID), scaled to reactor a=1.7 m
```
`load_qa("QA24")` expands the 4 base `CurveXYZFourier` (order 16) via
`coils_via_symmetries(nfp=2, stellsym=True)` -> 16 physical coils and multiplies by
`scale = 1.7/0.1683 = 10.10`. **This is the same 10.10 passed to `build_qa_fluxmap.py --scale`,
so coils and plasma land in one reactor frame.** The `.coils` MAKEGRID file is what ParaStell's
magnet builder consumes.

---

## 3. ParaStell build (cluster)

Mirror the QH radial build. ParaStell needs the cluster (Ginsburg) — Coreform Cubit +
`parastell` are not local. Inputs:
- LCFS + flux surfaces: from `precise_QA_reactor_fluxmap.npz` (`R`,`Z`,`phi` in reactor meters)
  or directly from the VMEC wout (ParaStell reads VMEC).
- `thickness_matrix` radial build: reuse the QH `thickness_field.py` machinery (in `shield_opt/`).
- Magnets: `data/processed/Wechsung_QA24.coils`.
- **Units:** wout/fluxmap are meters; OpenMC/DAGMC is cm. Apply the m->cm (x100) factor in the
  ParaStell build step, consistently for the plasma surface AND the coils.

Output: a DAGMC `.h5m` for OpenMC, exactly as for QH.

---

## 4. OpenMC transport

Point the native `StellaratorSource` at `precise_QA_reactor_fluxmap.{meta,bin}` and run against
the QA DAGMC geometry, same driver as QH (`python/run_conformal.py` / `run_ginsburg.py`).
Cross-check: the free-streaming NWL directionality (iso vs A vs B modes) from
`python/ana_precise_qa.py` (anarrima, native frame) should match the OpenMC free-streaming run.

---

## 5. What runs where

| Step | Where | Why |
|---|---|---|
| VMEC solve (A1) | **Ginsburg** (or vmecpp venv) | no `vmecpp`/VMEC locally |
| Fluxmap synthesis (A2, `build_qa_fluxmap.py`) | **local** | numpy + netCDF4/scipy only (both present) |
| DESC regen (Route B) | **Ginsburg / DESC venv** | no `desc` locally |
| Coil adapter | **local** | simsopt 1.10.6 present |
| ParaStell build | **Ginsburg** | Cubit + parastell |
| Analytic anarrima NWL | **anarrima venv** (`ana_precise_qa.py`) | jax + anarrima |
| OpenMC DAGMC transport | **Ginsburg** | DAGMC-enabled OpenMC + geometry |

---

## 6. Open items / caveats

- **PHIEDGE sign/magnitude:** the vendored INDATA has PHIEDGE=+0.0838573 at unit scale.
  If VMEC++ complains about flux sign, flip per the QH convention; does not affect b_hat direction.
- **iota check:** after the VMEC solve, confirm iota ~ 0.42 across the profile (Landreman-Paul QA
  signature) as an equilibrium sanity check before trusting the map.
- **Scale consistency is the one thing to not get wrong:** fluxmap `--scale` and
  `coil_adapter` scale must both be 10.10; ParaStell m->cm must be applied to both. A mismatch
  puts the source inside/outside the vessel.
- The existing `data/precise_QA_fluxmap.npz` is **unit-scale DESC** — do not feed it to the
  reactor ParaStell build; use it only for the native-frame anarrima benchmark.
