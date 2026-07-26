# Commute brief — lit review + answers (2026-07-16)

Reading order: **§0 first** (the Thea papers change the competitive picture and answer several of your
questions with hard numbers). Then the Q&A sections. Items tagged **[LIT PENDING]** are being
checked by background research agents; I'll append their findings.

Saved papers: `spf_prototype/docs/papers/thea_2512.08027.pdf` (Helios overview, Swanson et al. Dec 2025)
and `thea_helios_leadlithium_blanket.pdf` (Pasmann et al., the neutronics paper). The Thea neutronics
*slideshow* is behind a Cloudflare bot-wall; curl can't fetch it — an agent is trying WebFetch.

---

## §0. The Thea papers — this is the most important thing to read

The Helios blanket paper (Pasmann, van Riel, Khera, Romano, Koen, Donovan, Swanson, **Gates**; and the
overview has **E.J. Paul, D.A. Gates, S. Pasmann**) is essentially a parallel effort in *exactly* our
space, using the *same stack* (OpenMC + DAGMC + weight windows). Read it as both a validation and a
competitive clock.

**What they did (so we DON'T claim it):**
- A **stellarator neutron source in OpenMC** (§3.2): a point cloud of discrete `IndependentSource` sites,
  **isotropic** 14.1 MeV, with **realistic Bosch–Hale fusion-reactivity weighting** `s_i = I_i V_i / (I_T V_T)`
  (core-peaked, Fig. 5). → "native stellarator source in OpenMC" is **not** novel. Our differentiators:
  **spin polarization** (angular emission) and a **compiled continuous sampler** vs. their discrete cloud.
- **3D DAGMC forward neutronics** with **global variance reduction**: they use **The Random Ray Method
  (TRRM) + FW-CADIS** for neutron weight windows (§3.3), noting it beats MAGIC and gives uniform
  uncertainty in deep-shield low-flux regions. **Photon transport is not yet in TRRM** → photon WW via MAGIC.
- Forward **NWL map** (Fig. 10, poloidal×toroidal, Boozer coords): mean 0.95 MW/m², **peak 1.86 → peaking ≈ 2×**;
  "comparable to ARIES-CS NWL peaking ~2× and dynamic factor ~10×."
- **First-wall DPA** with their **cell-under-voxel (CUV)** method (`MeshMaterialFilter`): mean 11.7 DPA/FPY,
  min/max 4.6/19.0 → FW life 17.7 FPY (10.5–42.2). Non-uniform (higher on broad faces). Their proposed fix
  is **sector replacement**, not source reshaping.
- **Coil nuclear heating** (§3.6) — *this is our problem, published and unsolved*:
  - Design limit **150 W/m³** (cooling power). Binding constraint is **HEATING, not fluence.**
  - Mean coil heating **80.6 W/m³**, **max 179.3 ± 4.2 W/m³ = 2.2× the mean**, **8 of 324 coils exceed**,
    concentrated **inboard**. Stated remedy: "additional cooling or **shielding** considerations."
  - (Fluence limit is separately met: min shaping-coil life 40.3 FPY.)

**What they explicitly DEFER (= our contribution):**
- Blanket is a **uniform radial offset** (§3.1, "incrementally and uniformly offset").
- §2.3: "**Future work will include additional shielding optimization... non-uniform radial thicknesses.**"

**Their radial build (Table 2)** — near-identical philosophy to our corrected build, thinner total (~100 cm):
FW 2.05 cm (V-4Cr-4Ti) · breeder **50 cm Pb-17Li** (86.2% PbLi/EUROFER/He/SiC), Li6 **65%** · interior shield
**23 cm WC + 7 cm B4C** · plasma vessel ~13 cm (SS316L + vac + borated water) · exterior **5 cm borated poly**.
TBR target **1.3**, met at Li6 65–70% (3D: 1.307). Power multiplication 1.2.
Note: their device is **QA, R=8 m, A=4.5, B0=6 T, Bmax-coil 20 T, min plasma–coil 1.2 m**, PbLi not FLiBe,
12 planar encircling + 324 shaping coils.

**Bottom line for framing:** Thea (Gates/Paul/Pasmann) is sophisticated, active, and publishing in our exact
niche as of Dec 2025. They have a **real, concentrated, inboard coil-heating hotspot (2.2× peaking, 8/324
over limit)** and **no spatial shield solution** — they defer it. Our closed-loop non-uniform shield
optimizer, **driven by coil hotspots**, is precisely the deferred piece. But we must be **surgical about
novelty**: claim polarization, attribution/adjoint, and the closed spatial loop — NOT "stellarator source in
OpenMC" or "3D stellarator neutronics," which they (and ParaStell) already have.

---

## §1. Why is neutronics almost always downstream? What do people actually do?

Short answer, and the Thea papers illustrate it perfectly: **the design pipeline is plasma → coils →
(freeze) → neutronics evaluation → hand-iterate**, and neutronics is kept out of the plasma/coil optimizer
for a few compounding reasons:

1. **Cost/differentiability mismatch.** Plasma+coil optimization (SIMSOPT/DESC/STELLOPT) runs thousands of
   cheap, differentiable evaluations of physics figures of merit (quasisymmetry, MHD, transport). A single
   3D Monte Carlo neutronics evaluation with deep-shield weight windows is minutes-to-hours and **non-
   differentiable** (stochastic). You cannot put it inside a gradient-based shape loop without a surrogate.
2. **Weak *shape* coupling, strong *scenario* coupling.** To leading order the neutron **source** is set by
   the fusion rate (power, profiles), and the wall/coil load is set by geometry+standoff — both of which are
   only loosely sensitive to the fine plasma-shape knobs the optimizer is actually turning. So neutronics is
   treated as separable and checked afterward.
3. **Tooling separation.** Different codes, teams, meshes. ParaStell/Stellarmesh/OpenMC live downstream of
   DESC/SIMSOPT by construction.

**What they do instead of co-optimizing:** they **design to a limit and add bulk shielding / margin.** Two
strategies, both visible in the Thea papers:
- **Buy margin with standoff** (Helios philosophy): make the plasma *less strongly shaped* so coils sit
  further out (1.2 m), leaving room for a **uniform** thick blanket+shield. Overview paper, verbatim: *"The
  Helios equilibrium is less strongly shaped, allowing the coils to be further away. This significantly eases
  the design of the breeding blanket while protecting the coils... ARIES-CS required a highly optimized
  non-uniform blanket that prioritizes shielding at the expense of breeding zone in certain places."*
- **Spatially optimize the blanket by hand** (ARIES-CS philosophy): strongly-shaped/compact → not enough
  uniform room → thin the breeder and insert shield at the hot spots, **by engineering judgment**.

So yes — essentially **"add enough shielding to meet the peak limit, iterate the plasma by hand if you can't."**
Peaking factor **is** effectively the binding constraint (it sets how thick the shield must be to protect the
*worst* coil), but nobody closes it in an automated spatial loop. That gap — **automated, spatially-resolved,
hotspot-driven** shield allocation — is the opening, and it exists precisely *because* MC-in-the-loop is
expensive and non-differentiable (which is why our cheap surrogate + adjoint engine matters). **[LIT PENDING:
confirming nobody co-optimizes neutronics peaking with plasma/coil shape — agent running.]**

---

## §2. Wall vs. coil shielding; why sacrifice breeder; why hasn't anyone done it for stellarators

**Why thin the *breeder* to make room for *shield* (rather than move coils or grow the machine)?**
- **Machine cost scales steeply with size.** Major radius drives everything (magnet stored energy, structure,
  building, cost). Pushing coils out to make room is the most expensive possible fix. The whole point of a
  *compact* stellarator reactor is to NOT do that. So the radial build between plasma and coil is a **fixed,
  scarce budget**, and within it breeder and shield compete. (This is the ARIES-CS vs Helios split above.)
- **Coils must stay close** because the field they must produce at the plasma scales with distance —
  Thea (blanket §2): *"The further away the shaping coils, the higher the magnetic field required on the coil,
  which incentivizes a minimum distance between coils and plasma."* Move coils out → need higher on-coil field
  → bigger/more-expensive HTS, higher stress. So standoff is bounded from both sides.
- **Coils are the irreplaceable, lifetime-limiting component** (40 yr), whereas the **first wall is a
  consumable** (Thea FW life ~17.7 FPY, and they plan **sector replacement**). So you spend breeder margin to
  protect the thing you can't replace.
- **Maintenance drives architecture** (Helios: sector removal between planar coils, 84-day biennial outage,
  88% capacity factor) — but that's about *access*, not about the breeder/shield trade per se.

**Why "wall shielding" isn't a thing but "coil shielding" is** (your instinct, confirmed): you can't put
shield *in front of* the first wall — it *is* the first surface. First-wall load is reduced only by reshaping
the **source** (plasma shape, polarization) or the **wall** (shape/armor) or by **replacing** it. The coil,
sitting behind ~1 m of blanket, is the only place a **shield** lever exists. Hence: **first wall = source/
geometry problem; magnets = shield problem.**

**Why hasn't anyone automated it for stellarators?** Three reasons, all now evidenced:
1. **It's genuinely hard and expensive** (MC-in-loop, non-differentiable, 3D CAD geometry regen per iterate).
2. **The tooling only just matured** — ParaStell (2024) and Stellarmesh (Thea) made parametric 3D
   CAD→DAGMC feasible; before that you couldn't cheaply regenerate a non-uniform stellarator blanket.
3. **The people who could just got here** — ParaStell **defers** the loop; Thea **defers** non-uniform
   thickness. Two independent groups parked at the same "future work." That's the strongest possible signal
   the automated spatial loop is unclaimed *and wanted*. **[LIT PENDING: agent confirming no closed loop.]**

---

## §3. A principled, plasma-power-specific attribution cutoff (not 5%)

Your instinct is right and the Thea numbers make it concrete. The cutoff should be the **engineering limit**,
which *is* plasma-power-specific because the loads scale with fusion power:

- **The load is power-normalized.** NWL and coil heating scale ~linearly with fusion power (W/m² and W/m³ per
  MW of fusion). Thea normalizes everything "per 958 MW." So a device at 2× the power density has ~2× the
  loads for the same geometry → the *same geometry* can be under-limit at low power and over-limit at high
  power. The attribution cutoff = "**where does load exceed the fixed material limit at THIS device's
  power**." That's inherently power- and device-specific, not a fixed percentile.
- **Concretely (Thea):** limit 150 W/m³; at 958 MW, **8/324 coils** exceed. At higher power more coils cross;
  at lower power none do. **The over-limit SET is the attribution target**, and its size is an output, not an
  assumption. This is strictly better than "top 5%."

**"Do they just build thicker walls, and how does that affect magnet distance?"** — Exactly the tension:
thicker shield → either (a) push coils out (expensive, higher on-coil field) or (b) thin the breeder (lose
TBR) at fixed envelope. Thea buys margin by **standoff** (1.2 m, uniform build); ARIES-CS buys it by
**spatial breeder→shield trade**. Our optimizer automates (b) *only where needed* (over the 8 hot coils),
spending the *minimum* breeder — which is the marginal-value idea below.

**Marginal-value weighting (you liked this):** weight each over-limit location by how *binding* it is —
`d(lifetime)/d(load)` or how much relaxing it moves the limiting constraint. For a **sharp-threshold** failure
(insulator dose, quench) the near-threshold coils dominate; for a **wear-out** mode (fluence accumulation) the
integral matters. Thea's coil case is a **hard cooling-power threshold** (150 W/m³) → near-threshold/over-
threshold coils dominate → attribute and shield exactly those 8. This turns "top 5%" into "the set whose
marginal protection is worth the breeder it costs." **[LIT PENDING: whether marginal-value/limit-based
attribution cutoff is used anywhere in fusion neutronics.]**

---

## §4. Device-to-device (QA/QH/QI) variation + concentration index

You want to (a) map how attribution region size / concentration varies across the QS zoo, and (b) formalize a
Gini/participation concentration index, with the thesis *"attribution concentration predicts shield-optimizer
efficacy."* I think this is a real, under-explored science contribution and it connects cleanly to our
saturation result:

- **Why devices differ:** QA (tokamak-like, weaker shaping, Helios) vs QH (stronger helical shaping, our
  device) vs QI have very different |B| spectra, standoff variation, and boundary curvature → different load
  distributions. Thea's QA coil hotspot is **inboard-concentrated** (8/324). Our QH baseline was a
  **localized consecutive cluster** (coils 27–30). Different concentration signatures already, n=2.
- **The concentration index** (Gini / participation ratio of the culprit field `C_v`, or of the coil-load
  distribution) makes "concentrated vs diffuse" a single number. Then the thesis is testable across a QUASR
  batch: does high concentration ⇒ localized shield works well (few patches, low breeder cost), and low
  concentration ⇒ localized shield saturates fast / you need broad coverage or a source-side fix?
- **Ties to our result:** we *measured* that a localized patch saturates on a leak-around floor. Concentration
  should predict *where that floor sits* and *how many patches* you need — the quantitative bridge between
  attribution (diagnosis) and shield optimization (cure). **[LIT PENDING: is QS-zoo load-distribution
  variation studied? Any Gini/participation concentration metric in fusion neutronics? — agent running.]**

Note on your aside "*higher-concentration devices may be the unpolarized ones*": worth testing directly —
polarization redistributes the source, which can either sharpen or spread the coil load depending on B
geometry. Don't assume; it's a measurement (and a nice result either way).

---

## §5. Slant path + your scattering question (this is a sharp question and it refines the idea)

You asked: does scattering randomize the neutron path, and what does that do to the slant-path effective-length
(the `t/cosθ` idea)? Yes — and thinking it through actually *connects the slant-path lever to our saturation
result*:

- **The `t/cosθ` slant-path enhancement is an UNCOLLIDED (first-flight) effect.** It's the extra attenuation a
  neutron gets by entering the shield obliquely and traversing more material *before its first collision*. It
  applies cleanly to the **direct/uncollided** component `exp(−Σ·t/cosθ)`.
- **Once a neutron scatters, it forgets its entry angle.** The direction randomizes; the "I entered at a slant"
  memory is largely lost after ~1 collision. So the scattered component's transmission is governed by the
  **removal cross section + buildup**, not by the original incidence angle.
- **Scattering also *lengthens* the effective path** via random-walk (total track length ≫ straight-line
  thickness) — this is exactly why the *removal* MFP is what matters and why buildup factors exist.
- **The punchline, tied to our data:** our δ-sweep showed the coil flux past ~1 mfp is **dominated by
  scattered / leak-around paths** (that's the saturation floor). Those paths have *forgotten* the entry angle.
  So **the slant-path (channel 2) lever of polarization mostly affects the uncollided/first-flight component,
  which is exactly the component that's already small deep behind thick shield.** Therefore:
  **polarization's coil-dose lever is probably dominated by channel 1 (source redistribution), with channel 2
  (slant path) a smaller, shallow-shield / near-surface effect.** That's a non-obvious, honest refinement — and
  it's *testable*: decompose the coil-dose change from polarization into uncollided (angle-sensitive) vs
  scattered (angle-washed) components in the MC. If channel 2 is real it lives in the uncollided tally.

This also means the **first-wall** slant-path story is cleaner than the coil one: the first wall sees mostly
uncollided source neutrons (transparent plasma), so incidence-angle effects survive there; deep at the coil,
scattering has washed them out.

---

## §6. Random Ray, the angular idea, and the adjoint plan

**Partial answer already from Thea's paper:** Random Ray **is** in OpenMC, **is** coupled to **FW-CADIS**
(an adjoint-importance method), and Thea *uses* it (TRRM) for weight windows — so the adjoint importance is
*already being computed* in this stack, but for **variance reduction**, not **source attribution**. Two gaps
Thea names that map onto your PI's-student idea:
- **CORRECTION (agent found this — do not overclaim):** *anisotropic / angular random ray is ALREADY
  PUBLISHED* — **Neame & Cosgrove, "Linear Sources and Anisotropic Scattering in the Random Ray Method," NSE
  2024, DOI 10.1080/00295639.2024.2394729** — flat→linear source and isotropic→**anisotropic scattering up to
  P3**, implemented in **SCONE** (not OpenMC). So "anisotropic random ray" as a method is **not novel.** What
  *is* still open, and is the defensible engine claim:
  1. It does **not exist in OpenMC's** random ray (OpenMC RR is MG-only, **isotropic source, scalar/angle-
     integrated flux, no anisotropic scattering**) — porting it in is engineering-novel;
  2. retaining an **angularly-resolved adjoint importance field** `ψ†(r,Ω,E)` as the *solution output* (not
     merely an anisotropic scattering source that still collapses to scalar flux) is under-explored even in
     SCONE, and is exactly what a **streaming-sensitive coil** attribution needs (precedent: FW/CADIS-Ω,
     arXiv:1612.00793, shows adjoint *angular* flux matters for deep-penetration/streaming);
  3. the **application** — angular adjoint importance → stellarator **coil source-attribution** — is unclaimed.
  So the methods-paper engine claim is (ii)+(iii), NOT "we invented anisotropic random ray."
- **OpenMC has no general continuous-energy adjoint solver**; the accepted routes are (a) **multigroup
  adjoint** (what random-ray/FW-CADIS effectively do) and (b) **reciprocity**: forward-transport a pseudo-
  source placed *at the coil* into the shield and tally where it reaches the plasma = the contributon/importance
  map. Both are implementable in-stack.

**Your "do both to validate" plan is exactly right** and is standard practice: forward MC with **source-region
tagging** on one reference config = ground truth; adjoint/contributon importance = the cheap engine; show they
agree; then run the adjoint for the sweep. **[LIT PENDING: full Random Ray doc read + adjoint-in-stellarator
novelty + whether angular random ray is being pursued — agent running; I'll append the doc-level details.]**

**[LIT PENDING: "first adjoint/contributon coil-load source-attribution in a 3D stellarator" novelty check.]**

---

## §7. Real emissivity — visualization, pluggable infra, and the factorization

- **Thea already uses real emissivity** (Bosch–Hale, core-peaked, Fig. 5) for their *forward* NWL. So our
  uniform-emissivity culprit map is a *simplification relative to Thea's forward source*. That's fine — but it
  means we should (a) build the pluggable-emissivity infra and (b) USE realistic profiles for the *science*,
  keeping uniform only for the clean decomposition. This makes our tool a superset of theirs (adds attribution
  + polarization + pluggable scenario) rather than a subset.
- **Pluggable weight function** `r(ρ,θ,ζ)`: keep the geometry+polarization kernel fixed; users drop in a
  profile (or a kinetic-profile file). Clean software boundary; lets others study their scenarios.
- **The 3D translucent-voxel visualization you sketched:** yes — render the plasma as a volume of voxels whose
  **opacity ∝ emissivity·contribution** and **color ∝ contribution-to-peak-coil-load**, so low-culprit plasma
  is transparent and culprit regions "light up" as solid colored blobs inside the boundary; add **rotating
  views + a couple of poloidal/toroidal cross-section cut-planes** so you can see interior structure (core vs
  edge culprits). This is a genuinely better figure than the surface paint, and it's the natural way to show
  the *volumetric* culprit field with realistic (core-weighted) emissivity. Cheap to build (matplotlib 3D
  scatter with alpha, or a proper volume render).
- **The factorization** `C_v ≈ [geometry kernel] × [polarization g] × [emissivity r]` is the reason the
  analytic tool is the decomposition workhorse: the three factors enter **multiplicatively** in the free-
  streaming first-wall limit, so a **factorial study {uniform, realistic} × {unpol, A, B/C} × {device zoo}**
  cleanly separates how much of the culprit *location* and *concentration* is set by geometry vs shifted by
  polarization vs shifted by profile. Full transport breaks strict separability (scattering mixes terms), so:
  **analytic factorization = the clean decomposition; MC = validation + the deep-coil (scattered) regime.**
  Best-case finding: a stable trend like "geometry sets location; polarization moves it by X%; profile changes
  concentration not location" — a transferable result, not just a method.

---

## §8. Two-paper split — I strongly agree, and the timing makes it urgent

Your instinct to split is right, and the Thea Dec-2025 paper makes the **methods paper time-critical** (stake
the tooling claims before the active groups — Thea/CNERG — extend into them).

- **Paper 1 (Methods / tools):**
  1. compiled **polarized** `StellaratorSource` (+ SPF in `TokamakSource`) — differentiate vs Thea's isotropic
     point-cloud and vs Bae's polarized tokamak source;
  2. **automatic non-uniform deformation** of an initially-constant-offset blanket (the thing ParaStell &
     Thea both defer);
  3. the **adjoint/contributon attribution engine** (+ angular random-ray if your PI's student's idea lands),
     validated against forward source-tagged MC;
  4. the **analytic free-streaming factorization** tool + pluggable emissivity.
- **Paper 2 (Science / across the QS zoo):** your 5-part spine — attribution (culprit maps, limit-based
  cutoff, concentration index) → concentration-as-predictor → interventions (source-side polarization+shape
  for the wall; coverage-not-thickness shield for the coils, with the leak-around floor) → geometry/
  polarization/scenario decomposition → adjoint-validated engine.

Honest note: Paper 1 must be **surgical on novelty** given Thea/ParaStell. The claims that survive are
polarization, attribution/adjoint, closed-loop non-uniform shield, and (maybe) angular random-ray — NOT
"stellarator neutronics in OpenMC."

---

## §9. Literature-search results (3 agents, cited) — the novelty verdicts

Tags: **[EST]** read from primary/authoritative source · **[INF]** inferred from adjacent sources ·
**[NULL]** searched hard, apparently unpublished · **[CNV]** could-not-verify (paywalled/under-indexed).
Numbers from abstract-level snippets are flagged — spot-check primaries before quoting in a paper.

### 9.1 Why neutronics is downstream / co-optimization (agent A)
- **[EST]** The stated reason is overwhelmingly **MC is ~5 orders of magnitude too slow for an optimization
  loop** + it needs a **fully-built CAD device**. Lion et al. 2022 (Nucl. Fusion 62 076040, DOI
  10.1088/1741-4326/ac6a67) say resolving neutronics needs 1e8–1e9 histories, "tractable only by HPC," and
  their deterministic surrogate is "~5 orders faster." Bogaarts & Warmer 2024 (arXiv:2411.16369) echo it.
- **[EST]** Practice is **limit-driven radial-build sizing**: set a peak-NWL / magnet-dose limit, then tailor
  (often non-uniform) blanket+shield to meet it. Peaking managed by **first-wall shaping** (Lion: peaking
  **1.69→1.23**) or **non-uniform shield** (ARIES-CS). STELLOPT/SIMSOPT/DESC objectives contain **no
  neutronics**.
- **[EST/strong NULL]** Neutronics enters the *inner* plasma/coil optimization loop **only via a geometric
  standoff proxy** — min plasma–coil distance, or the **magnetic-gradient-scale-length** differentiable proxy
  (DOI 10.1088/1741-4326/ae69fa). **NWL-peaking / magnet-dose as a transport-derived objective on plasma
  boundary or coil shape is NOT demonstrated** and is **repeatedly named as future work** (ParaStell: "coupling
  ParaStell to machine-driven optimization is planned"; Lion: "can then... be used within an optimization
  framework"). → confirms §1/§2 above.

### 9.2 Adjoint / contributon attribution + Random Ray (agent B) — the strongest novelty
- **[NULL]** **Adjoint/contributon backward attribution of first-wall OR coil load to plasma SOURCE REGIONS in a
  3-D stellarator = unpublished.** Searches over CNERG/ParaStell, W7-X, HELIAS, Thea, ORNL/Denovo/Shift,
  contributon theory returned only (i) *forward* NWL maps (Lion 2022, ParaStell), (ii) adjoint for *tokamak
  diagnostics* (reciprocity, axisymmetric — arXiv:1504.03073), (iii) *geometry/MHD* adjoints (Paul/Antonsen).
  **Defensible novelty** if framed as *backward source-region attribution* (contributon map) vs the abundant
  forward maps. Formal basis: contributon theory, Williams UWFDM-1338. Hedge with "to our knowledge" (proceedings
  under-indexed).
- **[EST]** **OpenMC has NO continuous-energy adjoint.** Only **multigroup Random-Ray adjoint** (transpose
  scattering matrix + swap νΣ_f↔χ), used for FW-CADIS weight windows. Workarounds: MG RR adjoint; **reciprocity**
  (detector-side forward source); CE perturbation/sensitivity (not an importance field); external S_N (Denovo/
  ADVANTG). CADIS/FW-CADIS = the mature hybrid VR standard (Wagner; MS-CADIS for ITER shutdown dose, CNERG).
- **[EST]** **Random Ray = stochastic Method of Characteristics, multigroup only, isotropic source, angular
  dependence integrated out (scalar flux), no anisotropic scattering in OpenMC, flat/linear source; adjoint IS
  supported** (`settings.random_ray['adjoint']=True`). Tramm et al. NSE 2023.
- **[EST] correction to §6:** anisotropic random ray is **already published** (Neame & Cosgrove, NSE 2024,
  SCONE, P3). Open + defensible: angular-resolved adjoint *output* field, porting anisotropy into OpenMC, and
  the stellarator-coil-attribution application.

### 9.3 QS-zoo variation, concentration metrics, wall-vs-coil trade (agent C)
- **[NULL]** **No systematic cross-configuration (QA/QH/QI) or QUASR-wide NWL-distribution/peaking study.** Lion
  2022 computes peaking per-config (HELIAS-3/4/5 + a QA case) but it's a *method* paper, not a controlled zoo
  comparison. **QUASR has only QA+QH** (~370k configs; Giuliani, arXiv:2409.04826) — a QI-inclusive comparison
  needs QI configs from elsewhere (W7-X line / Stellaris). Novelty = "peaking computed device-by-device but never
  **mapped across the quasisymmetry landscape**." Do NOT claim to invent NWL peaking.
- **[CNV→opportunity]** **Gini / Lorenz / participation-ratio concentration metrics are NOT used in fusion
  neutronics** — only crude peak/mean. Importing them is **novel in application** (not in invention). Participation
  ratio (PR = (Σvᵢ²)²/Σvᵢ⁴ = "effective number of hot spots") is the most natural for "how many patches carry the
  load." Frame as borrowed from economics (Gini) / condensed-matter (PR).
- **[EST]** **Wall-vs-coil / breeder-for-shield trade is well documented.** ARIES-CS swapped blanket→WC shield at
  min-standoff spots, cutting local plasma-to-mid-coil build **~1.79 m → ~1.31 m (~20–30%)** while preserving
  **TBR ≈ 1.1** (El-Guebaly UWFDM-1336). Reason: machine size/cost scale steeply with standoff; coils must stay
  close (field-on-coil rises with distance). Corroborated by Helios (1.2 m budget) and ARC (FLiBe does triple duty).
- **[NULL, but "reviewer may call it obvious"]** **"Attribution concentration predicts localized-intervention
  effectiveness" appears unpublished.** Defend it by (a) quantifying with a real concentration metric, (b) making
  it a *falsifiable* cross-device prediction, (c) showing failure cases. **Honesty flag from the agent, which I
  second:** this is the same shape of claim as our earlier cheap-predictor that did NOT validate at n=5 (benefit
  tracked headroom, not shape) — hold it to the same skeptical bar before asserting.

### 9.4 Net novelty scorecard (what survives adversarial review)
| Claim | Verdict |
|---|---|
| First **native, reusable, general** stellarator `Source` class in OpenMC (+ **polarization**) | **Novel** — Thea has *a* stellarator source *in* OpenMC but it is an internal workflow emitting a **discrete point-cloud of IndependentSource sites** (isotropic, scenario-specific), NOT a drop-in reusable Source class; Bae's polarized source is tokamak. **Wording:** claim "first native/reusable/general," NOT "first stellarator source in OpenMC." **To substantiate "others can use," upstream it or release a documented package.** |
| Automated **non-uniform** blanket/shield spatial optimization (closed loop) | **Novel** (ParaStell + Thea both defer it) |
| **Adjoint/contributon coil-load → plasma-source attribution in 3D stellarator** | **Novel [NULL]** — strongest claim |
| Coverage-not-thickness + leak-around floor finding | **Novel** (our measurement) |
| Concentration index (Gini/PR) for load/attribution + "predicts efficacy" | **Novel in application**; the *predictor* thesis needs skeptical validation |
| Cross-QS-zoo peaking/concentration map | **Novel** (never mapped across the landscape) |
| "Stellarator source in OpenMC" / "3D stellarator neutronics" | **NOT novel** (Thea, ParaStell) — don't claim |
| "Anisotropic random ray" as a method | **NOT novel** (Neame & Cosgrove 2024) — don't claim |
| NWL peaking / magnet dose as plasma/coil optimization objective | **[strong NULL]** — future-work everywhere; a separate big claim if you pursue it |

### 9.5 Still open
- Thea neutronics **slideshow** stayed Cloudflare-blocked (agent couldn't fetch it either). Not critical — the
  two papers cover the substance. Grab it from your own browser if you want the slides.
- Several ARIES-CS / Lion exact numbers are abstract-level; pull primary PDFs (Lion TU/e copy; El-Guebaly
  UWFDM-1336) before quoting figures in a manuscript.
