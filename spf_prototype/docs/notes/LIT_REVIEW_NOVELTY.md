# Novelty & Prior-Art Literature Review — SPF Neutron Source + Stellarator NWL

**Purpose.** Establish what is genuinely novel vs. already-published for the two methods papers, and identify the closest prior works to cite and differentiate against. Written adversarially: the goal is to surface work that would *undermine* a novelty claim, not confirm it.

**Date:** 2026-07-07. Sources are cited inline with URLs / arXiv IDs. All verdicts reflect what could be verified from primary sources (papers fetched and read, not just search snippets).

---

## TL;DR — verdicts up front

| Claimed contribution | Verdict | Who got closest / why |
|---|---|---|
| **1. Native SPF polarized angular neutron source in OpenMC** | **PARTIALLY NOVEL — weakest claim** | **Bae et al. 2025 (Nucl. Fusion 65 086051)** already sampled the polarized DT angular distribution `W(θ)∝(a,b,c)` in OpenMC for a spherical tokamak. They did it as a *Python-precomputed point-source cloud*, not a native C++ compiled source with on-the-fly local-B̂ sampling. That distinction is real but modest. **Do not claim "first polarized fusion source in a MC transport code."** |
| **2. StellaratorSource sampling from real VMEC equilibria (extending Peterson TokamakSource)** | **PARTIALLY NOVEL** | The *concept* of a VMEC/flux-surface stellarator neutron source for MC transport already exists (HELIAS source routine in MCNP/Serpent2; Wu & Warmer et al. 2024). Novelty is: **in OpenMC**, built on Peterson's parametric `TokamakSource` API, **with local B̂ direction carried for polarization**. The polarized + VMEC + OpenMC combination is new; "stellarator neutron source from VMEC" alone is not. |
| **3. Free-streaming polarized NWL on DAGMC conformal walls for real QUASR devices; polarization-optimal a₂ per device** | **NOVEL (in combination)** | No prior work computes *polarized* NWL on *real optimized-stellarator* conformal walls, nor optimizes a₂ per device. Unpolarized stellarator NWL is well-trodden (Lion 2022; HELIAS DAGMC studies). The polarization axis on real stellarator geometry is the new part. |
| **4. anarrima-style analytic *perturbative* polarized NWL for weakly-shaped stellarators (field-period series, singularity subtraction, differentiable polarized quadrature)** | **NOVEL** | Schwartz 2025 (arXiv:2507.11758) does the *axisymmetric* analytic polarized NWL and ships `anarrima`, but explicitly restricts to axisymmetric convex cross sections. No analytic/perturbative field-period expansion of NWL was found in the literature. This is the strongest novelty claim of the four. |

**Bottom line:** Contribution 1 is where a reviewer will push hardest, because Bae et al. beat you to "polarized angular fusion source inside OpenMC." Frame the paper around the *native, differentiable, local-B̂, stellarator-capable* source and the *stellarator NWL / a₂-optimization* results, not around being first-to-OpenMC.

---

## A. Spin-polarized fusion (SPF) neutronics — who has modeled transport / wall loading?

**The central novelty question, and the one with the most dangerous prior art.**

### A.1 The one that hurts: Bae et al. 2025 — polarized angular source already in OpenMC
- **J. W. Bae, K. Borowiec, V. Badalassi, J. Parisi, A. Diallo, J. Menard, A. Khodak, T. Brown, "Neutronics analysis of spin-polarized fuel in spherical tokamaks," Nucl. Fusion 65, 086051 (2025).** DOI 10.1088/1741-4326/adf3c6. Open access PDF: https://www.osti.gov/servlets/purl/2583820
- **What they did (verified from the PDF):** They implement *exactly* the polarized DT differential cross section `dσ_{a,b,c}/dΩ = (σ₀/4π)·W(θ)` with the **same (a,b,c) deuteron/triton polarization parameterization** you use (their Eq. 1–2; `a,b` vector, `c` tensor). They define three schemes — unpolarized `W=1`, "parallel" (anti-aligned) `W=(1+3cos²θ)/2`, and "perpendicular" (spin-aligned) `W=(9/4)sin²θ`. They **sample neutron birth directions from `W(θ)·sin(θ)` vs the pitch angle θ to the magnetic field**, generate 10M samples per scheme via the **OpenMC Python API**, compile them into **discrete point sources**, transport 40M particles through a DAGMC spherical-tokamak (STAR) model, and report first-wall interaction poloidal distributions, TBR, and magnet fluence. Their Fig. 3 is literally the sampled-`cosθ`-vs-analytic-`W(θ)sinθ` histogram check you also do in tier 1.
- **How your work differs (defensible deltas):**
  1. **Native compiled source vs precomputed cloud.** Bae et al. sample in Python and *freeze* the result into a static point-source list exported to the settings file. Yours is a native `CompiledSource` that samples birth position *and* direction on the fly, in the **local B̂-aligned frame per position** (constant or toroidal field), thread-safe under OpenMC's OpenMP. This is an architectural/generality difference, not a physics difference.
  2. **Local, position-dependent B̂ and geometry generality.** They fix a per-scheme global emission pattern in a tokamak; you carry local B̂ direction, which is what makes the extension to **stellarators** (3D, spatially varying B̂) possible at all.
  3. **They study TBR and magnet fluence with full scattering; you study free-streaming NWL patterns** and a₂-optimization. Different observable.
- **Reviewer's likely words:** *"Bae et al. (2025) already put the polarized (a,b,c) DT emission distribution into OpenMC and sampled neutron directions from it. Isn't your source just a compiled-code repackaging of that?"* — You need a crisp answer ready (native/differentiable/local-B̂/stellarator). **This is contribution 1's Achilles heel.**

### A.2 The analytic sibling: Schwartz 2025 (your own Paper-2 basis and closest analytic prior art)
- **J. A. Schwartz, "Analytic neutron wall loading from spin-polarized fusion in axisymmetric geometries," arXiv:2507.11758 (submitted 15 Jul 2025).** https://arxiv.org/abs/2507.11758
- **What it did (verified):** Closed-form, fully **differentiable** formulas for NWL on an **axisymmetric torus with convex poloidal cross section**, from filamentary ring sources, arbitrary fuel polarization and arbitrary local B̂. **Analytic only — explicitly excludes neutron scattering.** Ships the **`anarrima`** Python package (its §5). Provides the (a,b,c)→angular-factor parameterization and the ±43%/±22% inboard/outboard midplane oracle numbers your prototype validates against.
- **Prior-work chain it cites:** Kulsrud et al. (SPF physics foundation); Ciullo et al. (SPF book); Heidbrink (tensor/anti-B "Tensor D" mode); **Lion 2022** (stellarator NWL, numerical); Chapin & Price, Chan et al., Yang et al. (multifilament/ray-trace NWL). Note it frames itself as the *analytic* alternative to numerical multifilament NWL.
- **Differentiation:** This is *the* reference framework you build on. Your Paper 1 differs by being a *Monte-Carlo transport* source (scattering-capable, arbitrary geometry). Your Paper 2 differs by extending the analytic method **off-axisymmetry** via a field-period perturbation series — which Schwartz explicitly does *not* do (he restricts to convex axisymmetric, "with one exception").

### A.3 Note on the arXiv ID in the brief
- **arXiv:2512.09242 is NOT an SPF-source paper.** It is **Parisi, Schwartz, Wurzel, Rutkowski, Harter, "Isotope Production in Fusion Systems"** (submitted 10 Dec 2025). https://arxiv.org/abs/2512.09242 It uses *asymmetric* NWL as a transmutation lever and cites the SPF work, but does not itself model a polarized source. Don't cite it as the SPF-source primary reference.

### A.4 SPF physics-only prior art (cross section & reactivity, NOT transport) — the "distinguish from" pile
These model the polarized *nuclear physics* (cross section, reactivity, spin transport) but **not neutron transport / wall loading**. Cite them to establish physics provenance and to show the transport step is the gap they leave open:
- **Kulsrud, Furth, Valeo, Goldhaber (1982)**, foundational SPF, 50% σ enhancement for S=3/2 alignment.
- **Temporal et al., Nucl. Fusion 52, 103011 (2012)** — ICF ignition with polarized DT.
- **Ab initio polarized DT cross sections:** e.g. https://arxiv.org/abs/1803.11378 (PMC6341121).
- **Spin transport hydrodynamics of polarized DT plasma, arXiv:2204.11523** — computes neutron angular distribution/polarization from spin transport, still no wall-loading transport.
- **Sadler et al.** and the polarized-fuel-lifetime program (Frontiers in Physics 2024, https://www.frontiersin.org/articles/10.3389/fphy.2024.1355212/full) — depolarization/lifetime, not transport.

**Verdict A:** The polarized DT *angular distribution* is standard physics; putting it into a *MC transport code and computing wall effects* has been done exactly once before you (Bae et al. 2025, OpenMC, tokamak, precomputed point sources) and analytically once (Schwartz 2025, axisymmetric). **You are second, not first, to polarized-source-in-OpenMC.** Your defensible novelty is native/differentiable/local-B̂/stellarator.

---

## B. Parametric plasma sources in OpenMC / MC codes

### B.1 Peterson TokamakSource — the thing you extend (and it's brand new/unmerged)
- **OpenMC PR #3999, "native parametric tokamak source," author eepeterson (Ethan Peterson), opened 6 Jul 2026, status OPEN/unmerged.** https://github.com/openmc-dev/openmc/pull/3999 Adds a native `TokamakSource` C++ class: rejection-free conformal sampling from parametric flux surfaces (major/minor radius, triangularity, elongation, Shafranov & vertical shift), Python + XML interface. **Isotropic emission, no B̂, tokamak-only.**
- **Differentiation:** Your `StellaratorSource` generalizes this API to (i) 3D VMEC flux-surface geometry and (ii) anisotropic polarized emission with local B̂. Because PR #3999 is *unmerged and days old*, cite it as concurrent/base work; your extension is a natural but genuinely new capability on top.

### B.2 openmc-plasma-source (fusion-energy org / Shimwell et al.)
- **https://github.com/fusion-energy/openmc-plasma-source**, PyPI `openmc-plasma-source`. Provides `tokamak_source`, `fusion_ring_source`, `fusion_point_source`. **Tokamak/ring/point only, isotropic emission, no VMEC, no stellarator, no polarization.** Predecessor: `parametric-plasma-source` (makeclean / open-radiation-sources), a C++/OpenMC compiled tokamak source. https://github.com/open-radiation-sources/parametric-plasma-source
- **Differentiation:** All tokamak, all isotropic. None sample a stellarator or carry B̂.

### B.3 Physics-informed / experimental-profile tokamak sources
- **"Physics-Informed Plasma Profile Input for OpenMC Fusion Neutronics," Fusion Sci. Technol. 82 (2025).** https://doi.org/10.1080/15361055.2025.2525028 — builds OpenMC fixed sources from simulated/experimental tokamak plasma states. Tokamak, isotropic.
- **MC-kit `parametrized-tokamak-source`** (OpenMC/MCNP tabulated profiles). Tokamak, isotropic.
- **Tokamak DT source models for confinement modes** (ResearchGate 257318050) — source-strength profiles, isotropic.

### B.4 Stellarator source samplers from VMEC — this is the one to differentiate carefully
- **Wu, Warmer, et al., "Proof-of-principle of parametric stellarator neutronics modeling using Serpent2," Nucl. Fusion 64 (2024).** https://iopscience.iop.org/article/10.1088/1741-4326/ad4f9f — builds the stellarator neutron source **directly from VMEC flux surfaces (via the SBGeom package, Fourier-coefficient surface set, interpolated between surfaces)** and transports in **Serpent2** using the **HELIAS user-defined source routine (ported from MCNP)**. **Isotropic emission, no polarization.**
- **Differentiation:** The *idea* "sample a stellarator DT source from VMEC flux surfaces for MC transport" is **already published** (Serpent2/MCNP, HELIAS). Your `StellaratorSource` is new by being **in OpenMC**, on **Peterson's parametric API**, and — the real differentiator — **carrying local B̂ for polarized anisotropic emission**. Do **not** claim "first VMEC-based stellarator neutron source for MC transport." Claim "first *polarized*, and first in OpenMC."

**Verdict B:** Parametric tokamak sources in OpenMC are a crowded space (Peterson, Shimwell, makeclean, MC-kit). VMEC-based *stellarator* sources exist in Serpent2/MCNP (Wu & Warmer). Your combination — OpenMC + Peterson-API + VMEC + polarized B̂ — is new, but each ingredient except polarization has precedent.

---

## C. Stellarator neutron wall loading — the conventional baseline

Establish the standard so you can say what's different. Two method families:

**(i) Deterministic / fast NWL (line-of-sight, inverse-square):**
- **Lion, Warmer, Wang, "A deterministic method for the fast evaluation and optimisation of the 3D neutron wall load for generic stellarator configurations," Nucl. Fusion 62, 076040 (2022).** https://doi.org/10.1088/1741-4326/ac6a67 (verified from PDF). Splits the plasma into point sources and the first wall into elements, sums **isotropic** inverse-square contributions over visible source points — a "matrix method," few CPU-seconds, built for FW-shape optimization to flatten NWL heterogeneity. HELIAS. **Isotropic, unpolarized, numerical (not analytic).** *This is the single closest prior work to your stellarator-NWL results and your a₂-optimization framing — cite prominently and differentiate on polarization.*
- **Chapin & Price; Chan et al.; Yang et al.** — earlier ray-trace / multifilament NWL (as cited by Schwartz).

**(ii) Full Monte Carlo NWL on CAD/DAGMC:**
- **HELIAS neutronics via DAGMC** (Häußler, Warmer, Fischer et al., EUROfusion; e.g. WPS2-CPR 17310) — MCNP/DAGMC particle tracking on CAD, 2D NWL maps, 3D flux, TBR. **Isotropic.**
- **ARIES-CS** (Najmabadi et al.) — compact-stellarator power-plant NWL, MC. Isotropic.
- **Stellaris (Type One Energy), Helios (arXiv:2512.08027)** — recent QI-stellarator FPP neutronics, MC/DAGMC. Isotropic.
- **Newer deterministic:** discontinuous-Galerkin deterministic neutronics for stellarators, Nucl. Fusion (2025), https://iopscience.iop.org/article/10.1088/1741-4326/ae4e48. Isotropic.

**Verdict C:** Stellarator NWL is a mature field, but **uniformly isotropic-emission**. Nobody computes *polarization-dependent* NWL on a stellarator wall, and nobody optimizes an emission-anisotropy parameter (a₂) per real device. That is your genuine contribution-3 novelty; the *geometry/first-wall* machinery is well-precedented (Lion 2022 is the closest).

---

## D. Native spin-polarized sources in any MC transport code

- **Built-in polarized *fusion birth source*: none found in MCNP, Geant4, Serpent, FLUKA, or (before Bae 2025) OpenMC.** General-purpose MC codes do not carry polarized emission angular distributions natively. Confirmed by e.g. the polarized-nuclear-imaging literature explicitly stating spin/polarization angular effects "are not implemented in general Monte Carlo codes such as MCNP, EGSnrc, or Geant4" (S1120179725001917, GATE polarized extension, 2025).
- **Polarized *particle transport* extensions exist but are a different thing:** Geant4 tracks polarized leptons/photons (polarized positron source studies, arXiv:physics/0512192; Compton polarimetry); GATE polarized nuclear-imaging extension (2025). These transport *particle spin/polarization through matter* — they are **not a polarized fusion birth source** and do not model the DT `W(θ)` emission. Distinguish explicitly: your source sets the *birth angular distribution* from nuclear spin state; it does not track neutron spin during transport (nor does Bae).
- **The only true prior "polarized fusion source in a MC code" is Bae et al. 2025 in OpenMC (§A.1)** — via precomputed point sources, not a native sampler.

**Verdict D:** A *native, in-code, on-the-fly* polarized fusion source is essentially unprecedented (Bae's is a precomputed cloud; everyone else is unpolarized or does spin-transport-not-source). This is the most defensible framing of contribution 1 — "native compiled polarized source," not "first polarized source."

---

## E. Perturbative / analytic stellarator NWL (supports Paper 2)

- **Schwartz 2025 (arXiv:2507.11758)** is the only analytic polarized NWL, and it is **axisymmetric-convex only** — explicitly *not* extended off-axisymmetry. It is your launch point, not competition.
- **Lion 2022** is *fast* but **numerical** (matrix/inverse-square summation), not a closed-form perturbation series, and isotropic.
- **No analytic field-period / near-axisymmetric perturbation expansion of NWL was found.** Near-axis / field-period expansions in the stellarator literature (Landreman-Sengupta near-axis expansion, quasisymmetry construction — arXiv:2209.11849, 2204.10234, etc.) are for **MHD equilibrium and field geometry**, not for the neutron-wall-loading integral. No one appears to have Taylor/Fourier-expanded the NWL surface integral in a shaping (field-period) parameter, done analytic singularity subtraction of the inverse-square kernel, or built a differentiable *polarized* hybrid quadrature for it.
- **Adjacent but different:** deterministic DG neutronics (2025, ae4e48) and Lion's matrix method give *fast numerical* NWL; systems codes (stellarator PROCESS, Lion 2021) use *scaling-law* NWL. None are perturbative-analytic.

**Verdict E:** The perturbative analytic polarized NWL for weakly-shaped stellarators (field-period series + singularity subtraction + differentiable polarized quadrature) appears **genuinely novel**. The closest prior art is your own axisymmetric basis (Schwartz) and the numerical fast method (Lion). This is the cleanest novelty story of the set.

---

## Ranked closest prior works (most-dangerous first)

| # | Work | What they did | How your work differs |
|---|---|---|---|
| 1 | **Bae et al. 2025, Nucl. Fusion 65 086051** (OSTI 2583820) | Sampled polarized DT `W(θ;a,b,c)` neutron directions in OpenMC (Python-precomputed point sources), transported in DAGMC spherical tokamak; TBR & magnet fluence. | Native compiled on-the-fly source, local per-position B̂, thread-safe; free-streaming NWL & a₂-optimization; extends to stellarators. Same core physics — **your #1 differentiator target.** |
| 2 | **Schwartz 2025, arXiv:2507.11758 + `anarrima`** | Analytic differentiable polarized NWL, axisymmetric convex torus; ±43/±22% oracle. | You do MC transport (scattering-capable, any geometry) and, for Paper 2, an off-axisymmetric field-period perturbation series (he stays axisymmetric). |
| 3 | **Lion, Warmer, Wang 2022, Nucl. Fusion 62 076040** | Fast deterministic matrix NWL for generic stellarators; FW-shape optimization to flatten NWL. | Isotropic & numerical; you add polarization and MC transport / analytic perturbation, and a₂-per-device optimization. |
| 4 | **Wu & Warmer et al. 2024, Nucl. Fusion 64 (ad4f9f)** | VMEC/SBGeom stellarator neutron source in Serpent2 (HELIAS routine from MCNP), isotropic. | Yours is in OpenMC on Peterson's API and carries B̂ for polarized emission; VMEC-source-for-MC concept itself is precedented here. |
| 5 | **Peterson, OpenMC PR #3999 (2026)** | Native parametric `TokamakSource` (flux-surface conformal sampling), isotropic, tokamak. | You generalize the API to 3D VMEC + polarized anisotropic emission (StellaratorSource). Concurrent base work. |
| 6 | **openmc-plasma-source (Shimwell et al.) / parametric-plasma-source (makeclean)** | Parametric tokamak/ring/point OpenMC sources, isotropic. | Tokamak-only, isotropic, no B̂/VMEC. |
| 7 | **HELIAS DAGMC neutronics (Häußler/Warmer/Fischer, EUROfusion)** | Full-MC NWL/TBR maps on CAD stellarator FW. | Isotropic; you add polarization and real QUASR conformal walls / a₂ optimization. |
| 8 | **ARIES-CS (Najmabadi et al.)** | Compact-stellarator power-plant NWL (MC), FW/blanket design. | Isotropic, single design; you do polarization + multi-device optimization. |
| 9 | **Parisi, Schwartz et al. 2025, arXiv:2512.09242** | Isotope production via asymmetric NWL; cites SPF. | Application paper, not a source/method; not a novelty threat, but a natural citation for "why asymmetric NWL matters." |
| 10 | **Temporal et al. 2012; Kulsrud et al. 1982; ab-initio σ (1803.11378); spin-transport (2204.11523)** | Polarized DT cross section / reactivity / spin-transport physics. | Physics provenance only — none do neutron transport or wall loading. |
| 11 | **GATE/Geant4 polarized extensions; Geant4 polarized lepton/photon transport** | Transport of *particle* polarization through matter. | Different problem: you set a birth angular distribution, not spin transport. Cite to pre-empt "isn't polarized transport already in Geant4?" |
| 12 | **DG deterministic stellarator neutronics 2025 (ae4e48); stellarator PROCESS (Lion 2021)** | Fast/deterministic & systems-code NWL. | Isotropic, non-polarized, numerical/scaling. |

---

## Reviewer's strongest objections (anticipate and disarm)

1. **"Bae et al. (2025) already sampled the polarized (a,b,c) DT emission in OpenMC. What's new?"** — *The strongest objection to Paper 1.* Answer: Bae precomputed a **static point-source cloud in Python for a fixed tokamak scheme**; ours is a **native, thread-safe, on-the-fly `CompiledSource`** that samples in the **local B̂ frame at each birth position**, which is (a) differentiable/composable within OpenMC's source machinery and (b) the *enabling* generalization to spatially-varying-B̂ **stellarator** geometry — something Bae's fixed-scheme tokamak approach cannot do. Concede priority on "polarized-in-OpenMC"; claim novelty on native + local-B̂ + stellarator. **Cite Bae prominently and early.**

2. **"VMEC-based stellarator neutron sources already exist (Wu & Warmer, Serpent2/MCNP HELIAS)."** — Answer: correct for *isotropic* sources; the new element is **polarized emission tied to local B̂ on VMEC surfaces, in OpenMC on the Peterson API**. Don't overclaim the VMEC-sampling itself.

3. **"Stellarator NWL is a solved, mature problem (Lion 2022, HELIAS DAGMC, ARIES-CS)."** — Answer: yes, for **isotropic** emission. No prior stellarator NWL is polarization-dependent, and none optimize an emission-anisotropy parameter per real optimized (QUASR) device. That axis is new.

4. **"Isn't polarized transport already in Geant4?"** — Answer: that is *spin transport of particles through matter* (leptons/photons), a different problem; there is **no polarized fusion birth source** in general MC codes, and neither we nor Bae track neutron spin post-birth.

5. **"Your Paper-2 analytic method is just Schwartz with more terms."** — Answer: Schwartz is axisymmetric-convex and stops there; the **field-period perturbation series + analytic singularity subtraction + differentiable polarized hybrid quadrature** for weakly-shaped *stellarators* is not in his paper or anywhere else found. This is the cleanest novelty — lead Paper 2 with it.

---

## Honest self-assessment

- **Most novel:** Paper 2's perturbative analytic stellarator NWL (contribution 4), then contribution 3 (polarized NWL on real QUASR conformal walls + a₂ optimization).
- **Least novel / most exposed:** contribution 1's "polarized source in OpenMC," because **Bae et al. 2025 got there first**. It survives as "native/differentiable/local-B̂/stellarator-enabling," but the paper must be framed accordingly and must cite Bae up front. Do not use language like "first polarized fusion source in a Monte Carlo transport code" — it is false.
- **Watch item:** Peterson's `TokamakSource` (PR #3999) is *days old and unmerged*; coordinate framing so your StellaratorSource reads as a collaborative extension of concurrent work, not a land-grab.

---

### Source list (primary)
- Bae et al. 2025, Nucl. Fusion 65 086051 — https://www.osti.gov/servlets/purl/2583820 (DOI 10.1088/1741-4326/adf3c6)
- Schwartz 2025, arXiv:2507.11758 — https://arxiv.org/abs/2507.11758
- Parisi, Schwartz et al. 2025, arXiv:2512.09242 — https://arxiv.org/abs/2512.09242
- Lion, Warmer, Wang 2022, Nucl. Fusion 62 076040 — https://doi.org/10.1088/1741-4326/ac6a67
- Wu & Warmer et al. 2024, Nucl. Fusion 64 (Serpent2 parametric stellarator) — https://iopscience.iop.org/article/10.1088/1741-4326/ad4f9f
- Peterson, OpenMC PR #3999 (native parametric tokamak source, 2026) — https://github.com/openmc-dev/openmc/pull/3999
- openmc-plasma-source — https://github.com/fusion-energy/openmc-plasma-source
- DG deterministic stellarator neutronics 2025 — https://iopscience.iop.org/article/10.1088/1741-4326/ae4e48
- GATE/Geant4 polarized nuclear imaging extension 2025 — https://www.sciencedirect.com/science/article/pii/S1120179725001917
- Physics-informed plasma profile OpenMC source, FST 2025 — https://doi.org/10.1080/15361055.2025.2525028
</content>
</invoke>
