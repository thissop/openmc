# Coil+equilibrium databases and relational coil–plasma descriptors

*Research memo for the cross-device "geometry → magnet-load attribution / shieldability" study.
Written 2026-07-25. Scope: (1) which PUBLIC sources ship a plasma boundary/equilibrium AND a
realistic coil set together, ranked for a ParaStell → DAGMC → adjoint filament-coil pipeline;
(2) whether relational (coil × plasma) geometric descriptors are a defensible predictor axis for
magnet neutron loading. Grounds and extends `SCIENCE_STRATEGY.md` (which already found QUASR,
ConStellaration, CIEMAT-QI4X, Stellaris). Every non-obvious claim is cited inline. Honesty rule:
where a source is paper-only, gated, or unverified, it is flagged as such — do not treat a described
coil set as a downloadable one.*

---

## 0. Report-back (read this first)

**Top-3 datasets to pull (with coils):**
1. **QUASR** — the only large public dataset that ships filament coils for *every* device (~370k QA+QH,
   MIT, simsopt JSON, dimensionless/rescalable). This is the spine.
2. **simsopt bundled configs** (`simsopt.configs`: NCSX, HSX, W7-X, Landreman–Paul QA/QH) — a handful of
   named, vetted coil+boundary sets; also the conversion layer that reads QUASR coils into polylines.
3. **Wechsung "CoilsForPreciseQS"** (Zenodo 10.5281/zenodo.5975323) — validated precise-QA filament coils
   + boundary; a clean reference/validation geometry.
   *(For QI variety add **ConStellaration** boundaries — but you must generate coils yourself; and
   **CIEMAT-QI4X** as one reactor-relevant QI point design whose coils are described but not confirmed
   publicly deposited.)*

**Top-3 relational descriptors to compute first:**
1. **Coil–plasma standoff, min/mean** (per-coil, poloidally/toroidally resolved) — strongest physics
   support; sets available shield thickness; magnet load ≈ NWL_geom · exp(−standoff/λ).
2. **Coil–plasma conformity / normal-distance variance** — variance of standoff along each coil; at fixed
   mean standoff, a wavy coil produces an exponential local hotspot. Genuinely under-occupied on the
   neutronics axis.
3. **Inter-coil gap streaming solid angle (as seen from the plasma)** — fraction of outward solid angle
   that reaches a coil without threading the winding pack; the one mechanism that bypasses exp(−t/λ)
   shielding. Highest novelty, clear physics, essentially no quantitative prior art.

**One-line verdict:** *"Coil–plasma relationship" is a stronger and more defensible predictor thesis than
"plasma geometry alone" for the magnet side — standoff-family relational descriptors are the physically
correct predictors of magnet load, and adopting them promotes standoff from a lurking confound to a named,
measured predictor — provided the paper stays honest that the strongest descriptor (min standoff) may
dominate the rest and that conformity/streaming are physically motivated but not yet literature-validated.*

---

## 1. Databases with coils + equilibrium (ranked)

Ranking is by usefulness for a **filament-coil** ParaStell → DAGMC → adjoint pipeline that needs coils +
boundary together and reactor scaling.

| # | Source | Boundary? | Coils? | Reactor-scale | Format | License / access | nfp & QS-type | ~#devices | Download | Pipeline note |
|---|---|---|---|---|---|---|---|---|---|---|
| **1** | **QUASR** | Yes (VMEC + Fourier target surf) | **Yes — all devices** | Dimensionless R₀=1 m → rescalable | **simsopt JSON** + VMEC | **MIT**, open | **QA & QH**, nfp ≈ 1–6 | **~370,000** | Web navigator (per-device); Zenodo `10.5281/zenodo.10050655` (concept) / `10.5281/zenodo.13717741` (~12.7 GB tar); simsopt `get_data` (v1.10.6+) | **Primary.** Only large dataset with filament coils + matching boundary; browsable UI; rescale for reactor studies. Vacuum fields only (no finite-β). |
| **2** | **simsopt bundled configs** | Yes (VMEC inputs) | Yes (few named) | Experiment-scale | simsopt JSON + VMEC | Open (MIT-style) | NCSX nfp=3 (QA), HSX nfp=4 (QH), W7-X nfp=5 (QI-ish), LP-QA nfp=2, LP-QH nfp=4 | ~5–10 | `pip install simsopt`; `simsopt.configs` (`get_ncsx_data`, `get_hsx_data`, `get_w7x_data`) | Vetted small test geometries + the **converter** you'll use to turn QUASR coils into ParaStell polylines. NCSX/HSX/W7-X are one-filament modular-coil approximations. |
| **3** | **Wechsung CoilsForPreciseQS** | Yes (VMEC) | **Yes (QA, QA+Well)** | Dimensionless, rescalable | simsopt | "Other (Open)"; Zenodo `10.5281/zenodo.5975323` (~1.3 MB) | **QA (nfp=2)**; QH **not** confirmed in this deposit | ~1 base + length-bound variants | Zenodo zip; GitHub `florianwechsung/CoilsForPreciseQS` | High-quality precise-QA reference coils; validation-grade, not dataset-scale. |
| **4** | **ConStellaration** (Proxima) | **Yes (VMEC/DESC + metrics)** | **No** | QI reactor-relevant boundaries | parquet/JSON (`wout`, DESC) | **MIT**, open | **QI**, varied nfp | ~158k–182k | HuggingFace `proxima-fusion/constellaration` | Best QI *boundary* variety. **Ships no coils** — coil-complexity *proxy metrics* only; you must run stage-2 coil optimization before DAGMC. |
| **5** | **CIEMAT-QI4X** | Yes (fixed- & free-boundary VMEC; HINT-3D) | **Yes (filamentary set), described** | **Reactor-relevant** QI | (paper) filamentary coils | Paper: IOP `10.1088/1741-4326/ae54ad`; arXiv:2512.08825 | **QI, nfp=4** | 1 point design | **Public deposit UNVERIFIED** — coils described in-paper; request from authors | A concrete reactor-scale QI with a *simple* filament coil set + island divertor. Best single QI reactor point if coils are obtainable. |
| **6** | **Stellaris** (Proxima) | Yes (free-boundary VMEC) | Paper-only | High-field QI power plant | — | *Fusion Eng. Des.* `S0920379625000705` | QI | 1 | Not publicly deposited | Reactor-relevant high-field regime; cite as context, coils not a supply. |
| **7** | **Thea Energy Eos** | Yes (QA, nfp=2, R=2.7 m, A=6) | **Planar-array (dipoles), not filaments** | Reactor-scale | Papers/patents | arXiv:2502.07702, 2503.18960; US patents | QA nfp=2 | 1 | Not deposited as filament set | **Topology mismatch**: 12 planar TF coils + array of ~50 cm-offset planar shaping dipoles → does NOT map onto a filament adjoint source. Discussion/stretch only. |
| **8** | **ARIES-CS / HELIAS-5B / WISTELL-D** | Yes (in studies) | Paper-only (reconstructable) | Reactor studies | — | ARIES-CS IAEA FT/P5-26; HELIAS/Helios arXiv:2512.08027; WISTELL-D via ParaStell (Frontiers 2024) | QA (ARIES-CS), QI (HELIAS) | 1 each | Not clean public coil deposits | Canonical reactor context + neutronics precedent (HELIAS used ADVANTG; ParaStell built WISTELL-D DAGMC). Coils not drop-in. |

**Notes / honesty flags:**
- **QUASR coils are for ALL ~370k devices, not a subset** — confirmed on the Zenodo deposit ("~370,000
  vacuum-field stellarators along with the electromagnetic coils that generate them") and the JPP paper
  (Giuliani et al., arXiv:2409.04826). Engineering constraints are baked in (max curvature 5 m⁻¹, min
  coil–surface & coil–coil distance 0.1 m), which makes the coils reactor-*plausible* after rescaling.
  Exact per-nfp / per-QS-type counts were not enumerated here — verify against the Zenodo metadata if a
  balanced sub-sample matters.
- **ConStellaration ships no coils** and there is (as of this memo) **no verified public release** of the
  matched coil dataset from the Proxima augmented-Lagrangian coil work (arXiv:2507.12681) or the
  coil-nonplanarity study (arXiv:2604.26763). Treat "ConStellaration coils" as not-public until confirmed
  with the authors.
- **"CoilForge"** (named in the tasking) **could not be verified** as an existing public stellarator coil
  dataset/tool — no matching deposit or paper surfaced. Do not cite it until located; the near-axis /
  single-stage coil *generators* (Giuliani, Jorge/Wechsung single-stage arXiv:2010.02033) are tools, not
  datasets.
- **Landreman–Paul precise-QS** ships as VMEC *boundaries only* inside simsopt
  (`input.LandremanPaul2021_QA` / `_QH`); no standalone coil Zenodo DOI distinct from the Wechsung deposit
  was found. Pair the boundary with §3 coils or re-optimize.

---

## 2. Relational (coil × plasma) descriptors

Each is computable from the two things the group already parses — filament polylines + a VMEC boundary.
"Support" = strength of the *published* physics link to **neutron/magnet loading** (not to plasma-surface
geometry). Physical kernel behind most rows: the deterministic 3-D NWL formula (Fischer/Bykov et al.,
IOP `10.1088/1741-4326/ac6a67`) is `Q ∝ Σ (n·l̂)/|l|²` — **inverse-square distance × incidence cosine** —
and magnet load ≈ NWL_geom · exp(−shield/λ) with shield λ ≈ 33–53 mm (Windsor et al., IOP
`10.1088/1741-4326/aabdb0`). Together these make *local minimum standoff* the theoretically dominant driver.

| Descriptor | Definition (from polylines + VMEC) | Cheap? | Physical rationale → target | Literature support |
|---|---|---|---|---|
| **(i) Standoff, min/mean** | For each coil point, distance to LCFS (KD-tree on sampled boundary or Newton on Fourier surface); report min + mean, poloidally/toroidally resolved | **Yes** — O(N·log N) | Sets available shield thickness; load ≈ NWL_geom·exp(−d/λ) with NWL_geom~1/d²; min standoff is the reactor size constraint → **attribution concentration & shieldability** | **Strong** — ARIES-CS (FT/P5-26); Helios (arXiv:2512.08027); ParaStell "greatest heating near min LCFS–coil spacing"; shield λ (10.1088/1741-4326/aabdb0) |
| **(ii) Conformity / normal-distance variance** | Variance of standoff `d(s)` along each coil arclength (and toroidal-position variance of per-coil min); "how parallel to a constant-offset winding surface" | **Yes** — reuse (i)'s field | At fixed mean d, high variance → a stretch dips close → exponential local hotspot & 1/r² NWL spike → **concentration** (peaky vs diffuse) | **Moderate / indirect** — REGCOIL premise (arXiv:1609.04378); ParaStell hotspot; *no* paper names it as a neutronics predictor |
| **(iii) Inter-coil gap streaming solid angle** | From LCFS source samples, fraction of outward solid angle reaching a coil without intersecting the winding-pack footprint (ray-cast; must inflate 1-D polyline to a pack width) | **Moderate** — needs pack model + ray-casting | Streaming through gaps bypasses exp(−t/λ) shield → dumps fast neutrons on coil cases → **concentration & (un)shieldability** | **Physical but unvalidated** — only qualitative "avoid straight gaps" (FNSF, OSTI 1239593); no stellarator magnet-load metric published |
| **(iv) Coil-to-coil spacing / minor radius** | Min coil–coil polyline distance ÷ plasma minor radius a (a from VMEC) | **Yes** | Wider gaps → larger streaming solid angle (worse, via iii); competes with engineering access → **concentration** | **Moderate / indirect** — canonical stage-2 penalty; neutronics role only via (iii) |
| **(v) Relative non-planarity (coil vs local surface)** | Coil excursion in the local LCFS *normal* direction not explained by following the surface (residual after projecting onto local offset surface) | **Yes** — reuse (i) frames | Proxy for conformity (ii); inherits (ii)'s weak indirect link | **Speculative for neutronics** — established result is manufacturability, best-predicted by plasma `pdrot`, **no neutronics** (arXiv:2604.26763) |
| **(vi) Coil vs plasma curvature/torsion** | Frenet κ,τ along coils vs LCFS principal curvatures at nearest point; mismatch \|κ_coil − κ_plasma\| | **Yes** | Correlate of (i)/(ii)/(iii), not an independent neutron driver | **Weak** — engineering/complexity penalty (Wechsung PoP 31,112501; FOCUS), not a neutron metric |
| **(vii) \|B\| ripple / field spectrum at coil** | Spectral content of \|B\| along coil (Biot–Savart on polylines) | **Moderate** | No published neutronics link; ripple is a *field-quality* quantity | **Speculative** — beware ε_eff "effective ripple" is a *neoclassical-transport* metric, unrelated to coil loads |

**Supporting physics scaffolding worth reusing:**
- **Deterministic 3-D NWL kernel** (IOP `10.1088/1741-4326/ac6a67`): ~1 CPU-s, gives the exact
  1/r²·cosθ rationale and can generate a cheap *surrogate label* for descriptor screening. Its known
  weakness — it overcounts at concave ridges because it lacks line-of-sight self-shadowing — is precisely
  the physics that makes a **concavity / streaming-solid-angle descriptor (iii)** matter.
- **Standoff is itself predictable from plasma geometry** via `min(L_∇B)`, the magnetic-gradient scale
  length (Landreman, arXiv:2309.11342; filament-coil follow-up arXiv:2602.18974). This is the physics
  bridge that lets you connect *plasma* geometry → *achieved* standoff → magnet load, i.e. it explains
  *why* standoff is the mediating relational variable rather than an accident.

---

## 3. Recommendation

### 3a. Which datasets to pull for an n > 2 study
1. **QUASR (primary, n large):** pull a **standoff-spanning, QS-balanced sub-sample** (both QA and QH,
   nfp 2–5) via the web navigator or the Zenodo tarball. This is the only way to get a *real n* with
   filament coils + boundary already matched, and it directly enables the SCIENCE_STRATEGY Step-2
   free-streaming coil-concentration law with the standoff decomposition — no new DAGMC builds.
2. **ConStellaration (QI boundaries) + self-generated coils:** for the QI point(s) needed to break the
   all-QS/all-QUASR degeneracy, take 1–3 QI boundaries and run simsopt stage-2 (or the augmented-Lagrangian
   recipe of arXiv:2507.12681) to make coils. This is the SCIENCE_STRATEGY Step-4 "QI feasibility spike."
3. **CIEMAT-QI4X:** one reactor-relevant QI point design with a described simple filament coil set + island
   divertor — request the coils from the authors; it is the most reactor-credible QI in the set.
4. **simsopt configs (NCSX/HSX/W7-X) + Wechsung QA:** as *validation/anchor* geometries — named real
   devices to sanity-check the descriptor pipeline and the DAGMC bring-up against literature.

So: **QUASR for n**, **ConStellaration + CIEMAT-QI4X for QI/standoff spread**, **simsopt/Wechsung as
anchors** — 3 core + 2 anchor sources.

### 3b. Which relational descriptors to compute first
Compute in this order (cheapest + best-supported first): **(i) standoff min/mean → (ii) conformity /
normal-distance variance → (iii) inter-coil streaming solid angle**, plus **(iv) coil–coil spacing/a** as a
near-free add-on. Defer (v)–(vii): (v) is subsumed by (ii); (vi) is a correlate; (vii) has no neutronics
link and a terminology trap.

Run all four against **free-streaming per-coil-flux concentration (Gini/PR) and peak** across the QUASR
sub-sample first (no DAGMC needed), using the committed `concentration.py` /
`descriptor_correlation.chatterjee_xi`. Explicitly test whether (ii)–(iv) add predictive power **beyond**
(i) alone (partial dependence / standoff-stratified ξ) — that is the whole scientific question.

### 3c. Is the relational framing defensible as the paper's predictor axis?
**Yes — with a stated caveat.** Three reasons it is stronger than "plasma geometry alone":
- **Physics is on your side.** The magnet load literally keys off coil-referenced quantities (available
  shield thickness = local standoff; streaming = inter-coil gap geometry). Plasma-shape descriptors
  (elongation, iota) drive the *first-wall* free-streaming projection, but the magnet side is through
  ~1 m of blanket+shield and is dominated by standoff and leak paths — i.e. relational quantities. This is
  the SCIENCE_STRATEGY §1.2 point, now backed by the standoff/shield-λ physics.
- **It dissolves the standoff confound.** The n=2 QA/QH comparison is confounded because QA/QH differ
  mostly in standoff (1.9×). Making standoff an *explicit named predictor* (descriptor i) converts the
  lurking variable into a measured axis; conformity/streaming (ii/iii) then test whether anything survives
  *beyond* standoff. That is exactly the decomposition a reviewer will demand.
- **The niche is open.** The nearest geometric-descriptor paper (arXiv:2604.26763) is on the
  coil-shape ↔ *plasma-surface* axis with **zero neutronics**. Your response variable (magnet neutron
  loading) is genuinely unoccupied, and (ii)/(iii) as *neutronics* predictors have no prior art to be
  scooped by.

**The honest caveat to state in the paper:** the strongest relational descriptor (min standoff) is so
dominant, via the exponential shield law, that it may *absorb* most of the predictive signal — so the
defensible headline is "the coil–plasma **relationship** (led by standoff) predicts magnet load, and
conformity/streaming are physically motivated candidate second-order predictors we test but do not yet
claim." Frame (ii)/(iii) as *hypotheses with clean falsification tests*, not established laws — consistent
with the group's n-honesty rules. Under that framing the relational axis is not just defensible, it is the
*correct* axis for the magnet side.

---

## Sources

**Databases / coils**
- QUASR — Giuliani et al., *J. Plasma Phys.* 2025, arXiv:2409.04826 (https://arxiv.org/abs/2409.04826);
  QA-only method arXiv:2310.19097; navigator https://quasr.flatironinstitute.org/;
  Zenodo `10.5281/zenodo.10050655` (concept) / `10.5281/zenodo.13717741`.
- simsopt — https://github.com/hiddenSymmetries/simsopt; configs docs
  https://simsopt.readthedocs.io/en/stable/simsopt.configs.html (`get_ncsx_data`, `get_hsx_data`,
  `get_w7x_data`); stage-2 example `examples/2_Intermediate/stage_two_optimization.py`.
- Wechsung CoilsForPreciseQS — Zenodo `10.5281/zenodo.5975323`;
  https://github.com/florianwechsung/CoilsForPreciseQS; PNAS 119, e2202084119 (2022).
- ConStellaration — arXiv:2506.19583; https://huggingface.co/datasets/proxima-fusion/constellaration;
  coil follow-ups arXiv:2507.12681, arXiv:2604.26763 (public coil release **unverified**).
- CIEMAT-QI4X — IOP `10.1088/1741-4326/ae54ad`; arXiv:2512.08825 (coils described; deposit unverified).
- Stellaris — *Fusion Eng. Des.* `S0920379625000705`.
- Thea Energy Eos — arXiv:2502.07702, arXiv:2503.18960 (planar dipole array; topology mismatch).
- ARIES-CS — IAEA FEC2006 FT/P5-26 (https://www-pub.iaea.org/MTCD/Meetings/FEC2006/ft_p5-26.pdf);
  HELIAS/Helios arXiv:2512.08027; ParaStell/WISTELL-D — Frontiers Nucl. Eng. 2024
  (10.3389/fnuen.2024.1384788).

**Relational-descriptor physics**
- Deterministic 3-D NWL kernel — IOP `10.1088/1741-4326/ac6a67`.
- Shield exponential attenuation (λ≈33–53 mm) — Windsor et al., IOP `10.1088/1741-4326/aabdb0`;
  IOP `10.1088/1741-4326/aa7e3e`.
- Standoff proxy `min(L_∇B)` — Landreman arXiv:2309.11342 (PPCF 10.1088/1361-6587/ad1a3e);
  filament follow-up arXiv:2602.18974.
- REGCOIL winding-surface conformity — Landreman & Hanson arXiv:1609.04378 (10.1088/1741-4326/aa57d4).
- Coil complexity penalties — Wechsung *Phys. Plasmas* 31, 112501; FOCUS (Zhu et al.).
- Neutron streaming lore — FNSF magnet studies, OSTI 1239593.
- Contrast paper (relational geometry, NO neutronics) — arXiv:2604.26763.
- Parametric stellarator neutronics — Lyytinen et al., IOP `10.1088/1741-4326/ad4f9f`.

*Unverified / could not confirm: "CoilForge" (no matching public dataset found); public deposits of
CIEMAT-QI4X, Stellaris, and ConStellaration-matched coils; a QH set inside the Wechsung deposit; exact
per-nfp/per-QS QUASR counts. Two 2026 arXiv IDs (2604.26763, 2602.18974) should have their final venue
confirmed before formal citation.*
