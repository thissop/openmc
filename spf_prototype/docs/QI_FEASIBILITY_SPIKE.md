# QI feasibility spike — per-device cost of adding a QI reactor to the adjoint pipeline

*Executes Step 4 of `SCIENCE_STRATEGY.md` §5: measure the true per-device cost of a NEW QS-type
(quasi-isodynamic) device BEFORE committing to n≈5. Every device so far is QA/QH; a QI point is the
standoff / QS-type spread the predictive law needs. Run 2026-07-26 on Ginsburg (login node only; no
sbatch jobs launched — this was a scope/de-risk spike, not a build). Nothing was committed; boundaries
staged under `/ginsburg/astro/users/tjk2147/spf_work/qi_spike/`.*

**Bottom line: CONDITIONAL-GO.** ConStellaration is fully reachable and I pulled 4 reactor-relevant QI
boundaries tonight (cheap, done). The equilibrium and DAGMC→adjoint steps are reusable from the QA/QH
pipeline at modest added cost. The **coil-design step is the long pole and is currently un-tooled**
(simsopt is not installed in ANY env; ConStellaration ships boundaries, not coils, and no coil benchmark).
First QI device ≈ **3–5 working days** human-in-the-loop (dominated by first-time simsopt/vmecpp env
setup + coil design + nfp handling); steady-state ≈ **1.5–2.5 days/device**. That is **under the
"> ~1 week ⇒ tight law out of scope" gate**, so n≈2–3 QI is feasible this cycle — **but only if the
coil-design script is de-risked on device 1 first, and only after the cheap kill-shots (Steps 1–2)
run.** Do NOT front-load QI DAGMC ahead of the transport kill-shot.

---

## 1. Access method (verified tonight)

- **Ginsburg login node has working outbound internet** (HTTPS 200 to `huggingface.co` and `pypi.org`).
  Login-node rule enforced: "processes longer than a few seconds / >1 core STRICTLY FORBIDDEN." All
  heavy work (parquet reads, pip installs, VMEC/coil solves) MUST go through `sbatch`/`srun`. Tonight I
  stayed within small curl calls only.
- **Best access path = the HF datasets-server REST API** (`datasets-server.huggingface.co`), NOT a bulk
  download. It exposes the full column schema and a SQL-like `/filter` endpoint, so I could select
  reactor-relevant boundaries server-side and pull only the `boundary.json` strings (≈2 KB each) —
  no parquet download, no package install. This is the cheap, repeatable way to slice the dataset.
  - Schema (per row): `metrics.qi`, `metrics.vacuum_well`, `metrics.aspect_ratio`, `metrics.max_elongation`,
    `metrics.edge_rotational_transform_over_n_field_periods`, `boundary.n_field_periods`,
    `boundary.is_stellarator_symmetric`, `boundary.r_cos`/`boundary.z_sin` (Fourier), `boundary.json`
    (serialized `SurfaceRZFourier`), plus an `omnigenous_field_and_targets.*` block and DESC opt settings.
  - Filter example (URL-encoded `where`): `"boundary.n_field_periods"=1 AND "metrics.vacuum_well">0 AND
    "metrics.qi"<0.02`, `orderby="metrics.qi" ASC`. The DuckDB index rebuilds on first hit ("index is
    loading, ~1 min") — retry after ~60 s.
- **Dataset layout** (`huggingface.co/datasets/proxima-fusion/constellaration`, `/tree/main`):
  - `data/train-00000..2-of-00003.parquet` (~610 MB total) = boundaries + all metrics + omnigenous-field
    targets. This is what the datasets-server indexes.
  - `vmecpp_wout/part.{0..1000+}.parquet` (each ~175–540 MB) = the **ideal-MHD VMEC++ equilibria**, packed,
    indexed by `vmecpp_wout/id_to_file_map.parquet` (2.5 MB). Equilibria are NOT individual small files —
    getting one wout means id_to_file_map → identify part.N → download that ~300–540 MB part → extract row.
    That is an sbatch job, not a login-node curl.
  - `finite_beta_{1..5}pct/` and `vmecpp_wout_finite_beta_*` = finite-β variants (equilibria at 1–5% β).
- **GitHub** `github.com/proximafusion/constellaration`: the `constellaration` PyPI package = geometry
  (`surface_rz_fourier.py`), omnigeneity/QI (`omnigeneity/qi.py`), a **VMEC++ forward model**
  (`forward_model.py`), and boundary-level optimization examples (incl. `launch_alm_simple_to_build_stellarator.py`).
  **No coil-design module and no coil sets anywhere in the repo** — the "simple-to-build" example optimizes
  the *boundary*, it does not produce coils. Two small example boundaries live in-repo
  (`hugging_face_competition/inputs/boundary.json`, `tests/geometry/test_data/test_boundary_*.json`).

## 2. QI boundaries pulled (landed on Ginsburg)

`/ginsburg/astro/users/tjk2147/spf_work/qi_spike/` (all parse cleanly: 5×9 Fourier `r_cos`/`z_sin`,
stellarator-symmetric, R0≈1.0 — see scale caveat below):

| file | nfp | metrics.qi (↓good) | vacuum_well (>0 = stabilizing) | aspect_ratio | max_elongation | edge ι/nfp |
|---|---|---|---|---|---|---|
| `boundary_nfp1_0_D25DWpUofpybrrLbbDtM6dt.json` | 1 | 6.1e-4 | +0.0264 | 8.61 | 6.63 | 0.081 |
| `boundary_nfp1_1_DF3hpfVE8MML46j5tUhZpiP.json` | 1 | 7.5e-4 | +0.0322 | 8.45 | 7.18 | 0.244 |
| `boundary_nfp2_0_DQSM8EDhMGCDxBDKyvKVvYL.json` | 2 | 1.8e-4 | +0.0023 | 8.21 | 7.15 | — |
| `boundary_nfp2_1_DEBFz4pTnvnra3TcNmNXSqs.json` | 2 | 1.9e-4 | +0.0133 | 7.70 | 6.88 | — |

Filter matched **2568 nfp=1** and **1178 nfp=2** boundaries with a positive vacuum well and qi<0.02, so
supply is abundant. All four have excellent QI (qi≈1e-4–1e-3) and a stabilizing magnetic well, i.e.
reactor-attractive per Kappel (small nfp → large plasma–coil separation). Aspect ratios ~8 and
elongation ~7 are on the high side (typical of the ConStellaration QI family) — worth noting for coil
standoff and for comparability against our QA/QH pair.
Also pulled the 2 GitHub example boundaries (`gh_*.json`, nfp=4) as format references.

**Equilibria: NOT pulled tonight** (each wout part is ~300–540 MB → needs sbatch). Two routes exist and
both are cheap-ish as batch jobs: (a) extract the shipped VMEC++ wout from `vmecpp_wout/` via
`id_to_file_map`; or (b) **regenerate** from the boundary with ConStellaration's `forward_model.py`
(fixed-boundary VMEC++), which is the cleaner route because it also rescales to reactor size.

## 3. Coil-design step — the long pole (scoped, NOT run)

**Finding: no coils, no benchmark, no installed tooling.**
- ConStellaration ships **boundaries only** — confirmed against the repo tree and pyproject. There is no
  reusable QI coil benchmark to copy.
- **simsopt is MISSING in every env on Ginsburg** (checked: `pstl`, `parastell`, `parastell_env`,
  `spf-stellarator`, base, `fiqus-env`, `fiqus312`). This matches the standing memo note "QUASR coils
  loaded from JSON, no simsopt" — our QA/QH work never needed a coil optimizer because QUASR *shipped*
  coils. QI has none, so we must design them, and the tool isn't here yet.
- **Install path is known and self-consistent:** `constellaration`'s own dependency set already pins
  `vmecpp==0.5.4` + `simsopt==1.8.1` + `desc-opt==0.17.1` + `booz-xform==0.0.9` (needs `numpy<2.3`,
  `pybind11<3`). So a single `pip install constellaration` into a **fresh** env yields BOTH the
  equilibrium generator (vmecpp) AND the coil-opt library (simsopt). Caveat: `numpy<2.3` conflicts with
  the pstl/parastell/openmc env → it must be a **separate** env; hand the VMEC wout + coil file to the
  pstl env for ParaStell. Wheels exist, but this is a multi-minute install with real dependency-resolution
  risk → sbatch, first-time debug budget.
- **Design workflow for one QI boundary** (must be written; not shipped): define a winding surface
  (uniform offset from the plasma boundary) → REGCOIL / current-potential solve (fast, linear, minutes)
  to get Bnormal-minimizing sheet current → cut discrete filament coils from the current potential (or run
  a simsopt stage-2 FOCUS-style filament optimization) → sanity-check Bnormal on the boundary + coil-coil
  clearance. **nfp=1 is the hard case**: least symmetry, coils wrap the full torus, contour-cutting is
  finicky and often needs manual iteration on #coils / regularization. **Prefer an nfp=2 boundary for the
  first device** to keep the stage-2 problem well-conditioned and comparable to QA.
- **Estimated cost:** did NOT run tonight (simsopt absent + login-node rule + >30 min ⇒ out of scope per
  instructions). First QI coil set: **1–2 days** human-in-the-loop incl. env bring-up and tuning;
  ~0.5 day/device once the script exists.

## 4. ParaStell → DAGMC → adjoint risk delta vs QA/QH

Real working env is **`pstl`** (`parastell 0.1.0`, `pymoab 5.5.1`, `cad_to_dagmc 0.11.9`, `gmsh 4.15.2`,
`openmc 0.15.0`) — the gmsh/cad-to-dagmc path, **no Cubit**. (The `parastell` env, by contrast, lacks
pymoab/cad_to_dagmc — don't use it for builds.)

| dimension | QA/QH (done) | QI (this boundary) | risk delta |
|---|---|---|---|
| **Equilibrium reader** | QUASR-shipped data | need VMEC wout ingested by ParaStell `invessel_build` | **NEW GAP**: pstl env has **no** VMEC reader (`pystell`/`pystell_uw`/`booz_xform`/`vmecpp` all MISSING). Wout ingestion path is untested here → integration + debug. |
| **nfp / symmetry** | QUASR QA/QH exploit nfp≥2 wedge + reflect | **nfp=1 forces full-torus build** | **BIGGEST geometric delta** for nfp=1: no 1/nfp shortcut → larger h5m, more faceting/memory in gmsh, slower transport. **Mitigation: use the nfp=2 QI boundary** (keeps the wedge, comparable to QA). |
| **Scale** | QUASR at physical scale | boundaries **normalized to R0≈1 m** | NEW per-device step: rescale plasma+coils+radial build consistently to reactor size. Easy silent unit bug — flag it. |
| **Material tags** | QA segfault on 9-tag set w/o `Vacuum` → device-parameterized `build_materials` | same machinery | REPEAT, not worse. Reuse the parameterized builder. |
| **Coil-cell ID** | `Cell.fill` NotImplemented → centroid-matching; Wechsung ×100 frame-match saga (221 cm offset) | self-designed coils | **Partly BETTER**: we own the coil frame/units → the QUASR ×100 scale saga disappears. But new coil count/topology → centroid-match thresholds need re-tuning. |
| **Adjoint + reciprocity** | random-ray adjoint (reactivity-independent, re-weight free) + forward per-coil | identical | REUSE as-is; 1 cluster job each + first-time coil-cell ID debug. |

Net: the DAGMC/adjoint core **reuses** cleanly; the genuinely new/risky items are (1) the VMEC-reader
gap in the pstl env, (2) nfp=1 full-torus cost (avoidable by choosing nfp=2), and (3) R0=1 m rescaling.
The QA coil-frame saga actually gets *easier* with self-designed coils.

## 5. Per-device cost estimate + go/no-go

**First QI device (end-to-end, human-in-the-loop):**
| stage | estimate | notes |
|---|---|---|
| Pull boundary | done tonight (~min) | datasets-server filter |
| Env: `pip install constellaration` (vmecpp+simsopt+desc) | 0.5 day | fresh env, dependency-resolution risk, sbatch |
| Equilibrium (VMEC++ from boundary, rescaled) | 0.25–0.5 day | forward_model.py; minutes to solve once env works |
| **Coil design (REGCOIL/stage-2 → discrete coils)** | **1–2 days** | dominant new cost; nfp=1 hardest |
| ParaStell radial build → watertight DAGMC | 0.5–1 day | reuse; nfp=1 full-torus + material-tag + watertightness debug |
| Adjoint + forward reciprocity validation | 0.5 day | reuse pipeline + coil-cell ID re-tune |
| **Total device 1** | **≈ 3–5 working days** | first-time setup heavy |
| **Steady-state per additional QI device** | **≈ 1.5–2.5 days** | coil script + env amortized |

This sits **under the SCIENCE_STRATEGY gate** ("if a device costs > ~1 week, the tight law is out of
scope"). So n≈2–3 QI this cycle is **feasible** (~1–2 focused weeks after device 1), and n≈5 for the
ambitious law is plausible-but-tight — consistent with the strategy's "n≈7 = illustrative trend, not a
law" honesty rule.

**Recommendation — CONDITIONAL-GO:**
1. **GO now (cheap):** boundaries are pulled; abundant supply. This part is retired.
2. **Prefer nfp=2** for the first QI DAGMC device (`boundary_nfp2_1_...`, AR 7.7, well +0.013, or
   `nfp2_0`): keeps the symmetry wedge, avoids the full-torus blow-up, most comparable to QA.
3. **GATE the coil→DAGMC→adjoint investment behind Steps 1–2** (the transport placement kill-shot on
   QH+QA, and the zoo-wide free-streaming coil law). Per the strategy, no expensive QI infra is bought
   until the mechanism survives Step 1. Do NOT front-load QI DAGMC.
4. **De-risk the long pole first:** stand up the fresh `constellaration` env and prove the coil-design
   script on **one** device end-to-end (the real "spike") before promising n=2–3. The coil step is the
   only unvalidated part; everything downstream reuses QA/QH machinery.

## 6. Honest limits — what I could NOT verify tonight

- **No sbatch jobs run** (login-node rule + spike scope): did not actually install simsopt/constellaration,
  did not generate or download a VMEC equilibrium, did not run a coil solve, did not build a DAGMC. All
  cost figures are **estimates** grounded in the QA/QH bring-up evidence (SCIENCE_STRATEGY §3), not measured.
- **Equilibria not obtained** — only boundaries. Confirmed they're available two ways but both need a batch job.
- **ParaStell VMEC-ingestion for these boundaries is unproven** in the pstl env (no VMEC reader installed);
  I inferred the gap from package probes, did not exercise the path.
- Did not confirm whether nfp=1 actually breaks ParaStell's periodicity handling vs merely being slower —
  flagged as risk, not tested.
- simsopt/constellaration install success on Ginsburg's numpy/compiler stack is **assumed from wheels**,
  not verified.
