# RESULTS — gentle published-QA free-streaming: OpenMC ↔ analytic on QUASR QA 59509

The toy QA-low cross-check ([[RESULTS_toy_qalow]]) repeated on a **real, published**
quasi-axisymmetric equilibrium from the QUASR database, so the Monte-Carlo paper's
analytic↔MC validation rests on a literature config, not only an idealized toy.

## Why a QUASR config (and how it was selected)
A DESC survey of the standard published QAs showed **none is natively low-ε_eff**
(precise_QA/reactor_QA 2.3, ESTELL 1.0, NCSX 1.3, ARIES-CS 1.6, WISTELL-A 1.6 at a
tight wall) — real QAs are strongly shaped, which is the MC regime. To find a genuinely
gentle published QA we mined the QUASR database **without the ~13 GB bulk download**:
- `quasr_catalogue.py` pulls the single `database.json.gz` (37 MB, **371,701 devices**)
  → our own CSV; filters QA (`helicity==0`, 71,794 in pool) and shortlists 60 gentle
  candidates by `max_elongation` + a random sample.
- `quasr_eps_eff.py` + `quasr_geom.py` download only the small per-device **VMEC
  namelist** (RBC/ZBS, ~18 kB) and compute the geometric ε_eff exactly in pure numpy.
  The metadata proxy works: **corr(max_elongation, e_geom) = 0.94**.
- Gentlest published QA: **device 59509** (nfp=3, aspect 6.7, R0=1, a≈0.175), with
  intrinsic excursion **e_geom ≈ 0.29** and, at the validation wall, **ε_eff ≈ 0.43** —
  comfortably in the analytic series/quadrature-convergent regime.

`simsopt`'s compiled core SIGILLs on this CPU (both surface `gamma()` and
`BiotSavart`), so geometry comes from the namelist and the field from a pure-numpy
model (below) — no simsopt at run time.

## Setup (`quasr_config.py`, shared by both sides via `SPF_CONFIG=quasr_config`)
- **Source:** device-59509 boundary (221 modes) + ρ=0.6 surface, 14 poloidal loops
  each, ρ² weighted. Geometry R,Z(θ,φ) evaluated directly from RBC/ZBS.
- **Field B̂:** flux-surface-tangent field at the device's transform,
  `B̂ ∝ ∂r/∂φ + ι·∂r/∂θ` (ι=0.1) — the equilibrium-implied field direction from the
  real boundary, NOT a toy field. It uses the VMEC poloidal angle rather than the
  straight-field-line angle (an O(shaping) approximation, small at ε_eff≈0.43). **Both
  the analytic and OpenMC sides call the identical B̂**, so the cross-check is valid
  regardless; B̂ only sets how realistic the pattern is, not whether the two methods agree.
- **Wall:** square-cross-section torus enclosing the boundary + 0.4·a gap (R∈[0.734,
  1.223], |Z|<0.264). Axisymmetric; all non-axisymmetry is in the source.
- **Analytic:** anarrima `free_streaming_quadrature(normalize=True)`, 384 nodes, trig
  order K=34 (≥ max |n|·nfp = 30); loop reconstruction max|err| **6.7e-16**.
- **OpenMC:** 3×10⁵ pre-sampled births/mode, free-streaming into the void wall.

## Result — agreement to MC statistics
Per-wall directionality A/iso, B/iso (the SPF steering):

| wall | A/iso MC | A/iso ana | rel | B/iso MC | B/iso ana | rel |
|---|---|---|---|---|---|---|
| inboard | 1.247 | 1.210 | 3.0% | 0.758 | 0.790 | 4.0% |
| outboard | 0.839 | 0.856 | 2.0% | 1.166 | 1.144 | 1.9% |
| floor | 1.042 | 1.055 | 1.2% | 0.958 | 0.945 | 1.4% |
| ceiling | 1.040 | 1.055 | 1.4% | 0.947 | 0.945 | 0.2% |

**Worst-wall discrepancy 4.0%** (inboard B/iso — lowest-current wall, noisiest); most
walls 1–3%. Consistent with MC statistics. Physics as expected: A (∝sin²θ) enhances
the inboard/floor and suppresses the outboard load; B/C (∝¼+¾cos²θ) is the mirror image
— now on a real QA at ε_eff≈0.43, three× more shaped than the toy.

Figures (`SPF_CONFIG=quasr_config plot_toy_qalow.py` -> `figs/`):
`quasr59509_directionality` (profiles, analytic vs OpenMC, 4 walls),
`quasr59509_agreement` (per-bin scatter on y=x), `quasr59509_geometry` (the rotating
QA cross-section in the wall).

## Reproduce
```
# 1) catalogue + shortlist + eps_eff (downloads ~37 MB once + small namelists)
$HOME/spf_venv/bin/python spf_prototype/python/quasr_catalogue.py
$HOME/spf_venv/bin/python spf_prototype/python/quasr_eps_eff.py
# 2) analytic THEN openmc (run analytic first; openmc reads its npz)
SPF_CONFIG=quasr_config $HOME/ana_venv/bin/python spf_prototype/python/toy_qalow_analytic.py
SPF_CONFIG=quasr_config PATH=$HOME/spf_venv/bin:$PATH OMP_NUM_THREADS=4 \
    $HOME/spf_venv/bin/python spf_prototype/python/toy_qalow_openmc.py
# 3) figures
SPF_CONFIG=quasr_config $HOME/spf_venv/bin/python spf_prototype/python/plot_toy_qalow.py
```

## Scope / caveats
Field is the flux-surface-tangent model above (not a VMEC `wout` or BiotSavart field —
simsopt unusable on this CPU); for a future upgrade, a `wout` would give the exact
field with no change to either driver. Single energy, filamentary/ρ-weighted-ring
source, square axisymmetric wall, free-streaming (no scattering). This is the low-ε_eff
QA analytic↔MC overlap; high-ε_eff exact QA (e.g. precise_QA) remains MC-only.
