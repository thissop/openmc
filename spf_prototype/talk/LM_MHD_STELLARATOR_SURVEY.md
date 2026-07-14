# Liquid-Metal MHD in Stellarator Geometry — State-of-the-Art Survey

**Prepared for:** T. Kiker (Columbia), for the spin-polarized-fusion (SPF) stellarator **FLiBe** blanket study.
**Question being answered:** How much does liquid-metal MHD constrain a *static-neutronics + non-uniform-shield-thickness* FLiBe stellarator blanket study, and is it worth coupling to an existing COMSOL/OFT liquid-metal MHD effort?
**Date:** 2026-07-13. Survey compiled from open literature via web search; claims flagged as *established*, *inferred*, or *sparse/absent* where the record is thin.

---

## TL;DR (bottom line up front)

1. **MHD in fusion blankets is fundamentally a *conductivity* problem.** The MHD pressure drop scales as `σ·B²·L·U` (roughly `∝ Ha²` in the conducting-wall limit). Liquid metals (PbLi, Li) have σ ≈ 0.7–1×10⁶ S/m; molten fluoride salt **FLiBe has σ ≈ 1–2×10² S/m — about 4 orders of magnitude lower**. This is *the* reason FLiBe blankets are routinely described as having "practically eliminated" MHD pressure-drop problems. (*Established.*)
2. **For a FLiBe blanket, MHD is a second-order effect, not a design driver.** The dominant blanket-design constraints for a neutronics/shielding study are tritium breeding ratio (TBR), shield thickness / coil dose, and thermal-hydraulics — **not** MHD pressure drop. MHD does *not* materially constrain a static neutronics + shield-thickness optimization for FLiBe. (*Established physics + inference for this specific study.*)
3. **There is essentially no validated 3D liquid-metal-MHD capability for a fully non-axisymmetric stellarator field.** All mature LM-MHD blanket codes assume near-straight ducts and near-axisymmetric (tokamak) fields. The one active stellarator LM-blanket line (EUROfusion/CIEMAT **HELIAS DCLL**) explicitly treats the 3D field as an open challenge and works around it with geometry (quasi-toroidal segmentation), not with a validated 3D solver. (*Established, and the gap is real.*)
4. **The novel opening is real but is a *liquid-metal* opening, not a FLiBe one.** 3D quasisymmetric-field × liquid-metal-MHD (gradient-B driven currents, non-transferable flow-channel-insert mitigation, field direction varying *along* a channel) is genuinely under-served. But it is only load-bearing if the blanket uses PbLi/Li. If your SPF study stays FLiBe, MHD coupling is **low priority** — worth a paragraph of justification, not a coupled simulation.
5. **On tooling:** COMSOL is a legitimate, community-validated LM-MHD workhorse (via coupled physics interfaces, since COMSOL has no dedicated MHD module). **OFT's "MUG" solver is a *plasma* extended-MHD code, not a liquid-metal duct-flow solver** — flag this, because "OFT MUG LM MHD" conflates two different meanings of "MHD." (*Established from OFT docs.*)

---

## 1. Core physics: MHD pressure drop & flow redistribution in LM ducts

A conducting liquid moving across **B** induces currents `J = σ(E + U×B)`; the resulting `J×B` Lorentz force opposes the motion and reshapes the velocity profile. The governing dimensionless groups (*established*, standard blanket-MHD framework — Smolentsev 2021 review [1]; Smolentsev 2010 [2]):

- **Hartmann number** `Ha = B·L·√(σ/(ρν))` — ratio of Lorentz to viscous forces; fusion blankets run `Ha ~ 10³–10⁴`.
- **Interaction parameter / Stuart number** `N = Ha²/Re = σB²L/(ρU)` — ratio of Lorentz to inertial forces; typically `N ≫ 1` (MHD-dominated, laminarized flow).
- **Wall conductance ratio** `c = σ_wall·t_wall/(σ_fluid·L)` — how much current the walls short-circuit; a first-order driver of pressure drop.

**Boundary-layer structure** (*established*):
- **Hartmann layers** (walls ⟂ B), thickness `∝ 1/Ha`, carry the return current.
- **Side / Shercliff layers** (walls ∥ B), thickness `∝ 1/√Ha`, much thicker.
- In conducting-wall ducts (Hunt flow) the side layers develop **"M-shaped" velocity profiles with high-velocity side jets** that can carry a large fraction of the flow — these jets matter for heat transfer, corrosion, and MHD instability. (Smolentsev 2021 [1]; classic duct-flow numerics [3].)

**Pressure-drop scaling** (*established*): for a conducting duct the fully-developed MHD pressure gradient `∝ σ·U·B²·(c/(1+c))`. In fusion conditions the MHD pressure drop exceeds the ordinary hydrodynamic drop by **5–7 orders of magnitude** [1]. Additional **3D pressure drop** arises wherever `B` varies along the flow or the duct bends/expands (manifolds, bends, fringing field): local 3D currents brake the core and push fluid to the walls, adding `∝ N`-scaled losses [4,5].

**Why σ is the whole game.** Because every MHD force term carries a factor of σ, the ~10⁴× conductivity gap between liquid metal (~10⁶ S/m) and FLiBe (~10² S/m) maps almost directly onto a ~10⁴× reduction in MHD pressure drop and current-driven flow distortion. Fluoride salts are repeatedly described in the blanket literature as having "extremely low electrical conductivity which practically eliminates MHD problems" [6, review of advanced blankets; 7, molten-salt MHD heat-transfer analysis]. This is why FHR/FLiBe and molten-salt blanket concepts largely sidestep the MHD program that dominates PbLi blanket R&D. (*Established.*)

---

## 2. State of the art in *modeling* LM-MHD for blankets

**Main groups & concepts** (*established*):
- **UCLA (Smolentsev, Abdou, Kirillov-lineage methods)** — the reference group for blanket LM-MHD. Anchors the **DCLL (Dual-Coolant Lead-Lithium)** concept, where PbLi flows slowly (~10 cm/s) for breeding/power conversion and the MHD pressure drop is mitigated by a **SiC flow-channel insert (FCI)** acting as an electrical (and thermal) insulator that decouples the PbLi from the conducting steel wall [8, DCLL status & R&D needs; 9, FNSF DCLL MHD thermohydraulics].
- **KIT (Karlsruhe)** — WCLL/HCLL and PbLi loop experiments, MHD benchmarking, insulating-coating R&D (e.g., Al₂O₃ coatings for PbLi [10]).
- **EUROfusion breeding-blanket program** — MHD design guidelines & numerical analysis for LM blankets [5,11].

**What current codes can actually do** (*established, with an important caveat*):
- Fully-developed and quasi-3D MHD **duct** flows at fusion `Ha`, including FCIs, fringing fields, manifolds, and mixed convection (buoyancy + MHD). Tools: **COMSOL** (coupled physics interfaces — see §4), UCLA's in-house **HIMAG/other codes**, OpenFOAM-based MHD solvers [12], and specialized fully-developed solvers.
- **Geometric/field assumptions are the key limitation:** the validated toolset assumes **straight or gently varying ducts** and **near-axisymmetric tokamak fields** where **B** is dominantly toroidal, roughly `1/R`, and (locally) unidirectional across a channel. High-`Ha` LM-MHD is numerically brutal (thin `1/Ha` Hartmann layers demand extreme resolution), so 3D full-geometry simulations remain expensive and are validated only against tokamak-relevant benchmarks (fully-developed duct, fringed-field duct, magneto-convection) [11,13]. There is **no community-validated, general 3D LM-MHD capability for arbitrary non-axisymmetric fields**. (*Established as a gap.*)

---

## 3. The stellarator-specific gap

**What breaks going from tokamak duct flow to a stellarator** (*inferred from established MHD physics + confirmed as open by the HELIAS work*):
- **|B| and field *direction* vary along a single channel.** In a tokamak a poloidal duct sees a roughly constant-direction toroidal field; in a stellarator the field rotates and its magnitude ripples along the flow path. The angle between **U** and **B** (which sets the Lorentz braking) is no longer quasi-constant.
- **Gradient-B and 3D currents everywhere.** 3D pressure-drop mechanisms that are *localized* to bends/manifolds in a tokamak become *distributed* along stellarator channels, and FCI mitigation tuned for a fixed field orientation is **not transferable** when the field direction sweeps.
- **No axisymmetry to exploit** for reducing the modeling problem — the whole point of stellarator symmetry (quasisymmetry is a symmetry of `|B|`, *not* of the vector field the fluid feels) does not simplify the LM-MHD boundary-value problem.

**Does any validated 3D stellarator LM-MHD capability exist? Essentially no.** (*Sparse/absent — this is an honest negative result.*)
- The only sustained stellarator *liquid-metal breeding-blanket* effort is the **EUROfusion/CIEMAT HELIAS DCLL** line (Palermo et al.) [14, Energy 2024; 15, Energies 2023 neutronics]. It explicitly acknowledges that the "un-even magnetic field" makes mixed-convection PbLi flow segment-dependent, notes that the **tangential field component suppresses flow instabilities but does *not* reduce total pressure drop**, and — tellingly — proposes a **"quasi-toroidal" blanket segmentation to minimize MHD pressure drop while avoiding complex insulating components**. That is a *geometry* workaround, not a validated 3D MHD solution; they flag magnetoconvective analysis in true HELIAS conditions as an open challenge. (*Established that they treat it as open.*)
- **Stellarator liquid-metal PFC/divertor** work exists but is a *different* problem (free-surface heat exhaust, not breeding-duct pressure drop): capillary-porous-system (CPS) targets with Li/Sn (OLMAT facility [16]), the **ASTER** "liquid stellarator" reactor concept with Li-on-molten-salt PFCs [17], and FFHR-d1 LM limiter/divertor concepts. W7-X uses solid divertors; no operating stellarator runs a liquid-metal breeding blanket.
- I found **no paper presenting a validated, fully-3D non-axisymmetric-field liquid-metal duct-MHD simulation for stellarator breeding channels.** The literature here is genuinely thin. (*Sparse — flagged honestly.*)

---

## 4. Tooling: COMSOL and the Open FUSION Toolkit (OFT / "MUG")

**COMSOL for LM-MHD** (*established*): COMSOL ships **no dedicated MHD module**; blanket MHD is done by coupling its Fluid-Flow (Navier–Stokes), Electric-Currents, and (optionally) Heat-Transfer interfaces. This coupled approach has been **verified & validated against the fusion MHD benchmark suite** (fully-developed duct, fringed-field duct, magneto-convection) — see the V&V studies [18, Energies 14:5413, 2021; 19, Fusion Eng. Des. validation]. COMSOL is a legitimate, widely-used LM-MHD workhorse for blanket ducts; its practical ceiling is high-`Ha` 3D resolution (mesh cost of `1/Ha`-thin layers).

**Open FUSION Toolkit (OFT)** (Chris Hansen et al., Columbia/U. Washington lineage) [20, OFT docs; 21, Columbia FRC]: a finite-element suite for plasma/fusion problems in 2D/3D geometry. Modules:
- **TokaMaker** — free-boundary time-dependent Grad–Shafranov equilibrium [22, arXiv:2311.07719].
- **ThinCurr** — 3D thin-wall eddy-current / electromagnetics.
- **Marklin** — 3D force-free (uniform-λ) equilibrium solver.
- **MUG** — **"2D/3D linear/nonlinear *extended-MHD* simulation package."**

**Critical clarification / sanity-check:** **OFT's MUG is a *plasma* extended-MHD solver (the fusion-plasma fluid, Hall/two-fluid terms, instabilities), NOT a liquid-metal engineering duct-flow / heat-transfer solver.** The phrase "OFT MUG LM MHD" conflates two distinct meanings of "MHD." I found **no evidence of a liquid-metal blanket-flow module in OFT.** In a COMSOL+OFT effort the natural—and defensible—division of labor is: **OFT provides the accurate 3D vacuum/equilibrium field** (`B(x)` on the stellarator geometry), and **COMSOL consumes that field to solve the liquid-metal Navier–Stokes + current problem.** If the collaborator is genuinely running liquid-metal duct MHD *inside* MUG, that would be a novel/off-label use worth verifying directly with them before relying on it. (*Established re: MUG's stated purpose; the coupling interpretation is inferred.*)

---

## 5. Bottom line for the SPF FLiBe stellarator blanket study

**How much does LM-MHD constrain your design? For FLiBe: very little.** (*Established physics, applied to your case.*)
- FLiBe's ~10⁴× lower conductivity makes MHD pressure drop and current-driven flow redistribution **second-order**. Your stated study — **static neutronics + non-uniform shield-thickness optimization** — is governed by TBR, shielding/coil dose, and geometry, none of which MHD moves at leading order. You can justify *decoupling* MHD from the neutronics optimization with a short scaling argument (σ_FLiBe/σ_PbLi ~ 10⁻⁴ ⟹ Ha² and pressure drop down ~10⁴), citing the molten-salt-blanket MHD literature [6,7]. **A coupled COMSOL/OFT LM-MHD run is not needed to make the neutronics/shielding result credible.**
- The one place FLiBe MHD is *not* fully negligible is **magneto-convective heat transfer** (even weak `J×B` can damp turbulent heat transfer and alter thermal stratification). If your study makes thermal-hydraulic or first-wall-temperature claims, cite this as a caveat rather than modeling it. (*Inferred.*)

**Where is the genuine research opening?** At the intersection of **3D quasisymmetric fields and *liquid-metal* MHD** — distributed gradient-B currents, field direction sweeping *along* a breeding channel, and FCI/insulation strategies that don't transfer from tokamak geometry. This is real, under-served, and only lightly touched by the HELIAS-DCLL workarounds. **But it is a PbLi/Li opening, not a FLiBe one.** For your SPF study the strategic choices are:
  1. **Stay FLiBe →** treat MHD as a cited-and-dismissed second-order effect; **do not** invest in the COMSOL/OFT coupling for this paper. (Recommended if the SPF neutronics is the contribution.)
  2. **Pivot/extend to a PbLi or dual-coolant stellarator blanket →** *then* the 3D LM-MHD × quasisymmetry coupling becomes a legitimately novel, publishable thrust, and coupling OFT's 3D field to COMSOL's LM solver is the right architecture — with the caveat that you'd be near the frontier of *validated* capability (expect to do your own benchmarking).

**Recommended framing for the collaborator conversation:** confirm (a) whether their "OFT MUG" work is plasma-MHD field generation vs. actual liquid-metal duct flow (they are different codes/physics), and (b) whether the blanket is FLiBe or PbLi — that single choice determines whether LM-MHD is a footnote or a co-equal research thrust.

---

## References

1. S. Smolentsev, *Physical Background, Computations and Practical Issues of the MHD Pressure Drop in a Fusion Liquid Metal Blanket*, **Fluids** 6(3):110 (2021). https://www.mdpi.com/2311-5521/6/3/110 (DOI: 10.3390/fluids6030110)
2. S. Smolentsev, R. Moreau, L. Bühler, C. Mistrangelo, *MHD thermofluid issues of liquid-metal blankets: Phenomena and advances*, **Fusion Eng. Des.** 85 (2010) 1196–1205. https://www.sciencedirect.com/science/article/abs/pii/S0920379610000645
3. *Numerical simulation of liquid-metal MHD flows in rectangular ducts*, **J. Fluid Mech.** https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/numerical-simulation-of-liquidmetal-mhd-flows-in-rectangular-ducts/C414DA36B3E66D7B4D509159C4D282F2
4. *MHD pressure drop and flow balancing of liquid metal flow in a prototypic fusion blanket manifold*, **Phys. Fluids** 30, 057101 (2018). https://pubs.aip.org/aip/pof/article/30/5/057101/911826
5. C. Mistrangelo, L. Bühler et al., *MHD flow in liquid metal blankets: Major design issues, MHD guidelines and numerical analysis*, **Fusion Eng. Des.** (2021). https://www.sciencedirect.com/science/article/abs/pii/S0920379621005718
6. *Review of blanket designs for advanced fusion reactors*, **Fusion Eng. Des.** (2008). https://www.sciencedirect.com/science/article/abs/pii/S0920379608002391 (molten salt: "extremely low electrical conductivity... practically eliminates MHD problems")
7. S. Smolentsev et al., *MHD Effects on Heat Transfer in a Molten Salt Blanket* (TOFE), UCLA/UW-Madison FTI. https://fti.neep.wisc.edu/fti.neep.wisc.edu/tofeprogram/pdf/SergeyATfusion.ucla.edu1083344123.pdf
8. S. Smolentsev et al., *Dual-coolant lead–lithium (DCLL) blanket status and R&D needs*, **Fusion Eng. Des.** 100 (2015) 44–54. https://www.fusion.ucla.edu/files/2019/08/FED-v100-Smolentsev-Dual_Coolant_Lead_Lithium_Blanket_Status2015.pdf
9. S. Smolentsev et al., *MHD thermohydraulics analysis and supporting R&D for DCLL blanket in the FNSF*, OSTI 1610064. https://www.osti.gov/pages/servlets/purl/1610064
10. *Development and validation of electrical-insulating Al₂O₃ coatings for high-temperature liquid PbLi applications*, arXiv:2107.03244. https://arxiv.org/pdf/2107.03244
11. *MHD R&D Activities for Liquid Metal Blankets*, **Energies** 14:6640 (2021). https://doi.org/10.3390/en14206640
12. *A multi-region and multiphase MHD OpenFOAM solver for fusion reactor analysis*, **Fusion Eng. Des.** (2024). https://www.sciencedirect.com/science/article/pii/S0920379624000693
13. *Characterization of the MHD flow and pressure drop in the access ducts of a liquid metal fusion blanket*, **Fusion Eng. Des.** (2024). https://www.sciencedirect.com/science/article/pii/S0920379624001157
14. I. Palermo et al., *Challenges towards an acceleration in stellarator reactors engineering: the DCLL breeding blanket helical-axis advanced stellarator (HELIAS) case*, **Energy** (2024). https://www.sciencedirect.com/science/article/pii/S0360544223033649
15. I. Palermo et al., *Neutronic Assessments towards a Novel First Wall Design for a Stellarator Fusion Reactor with DCLL Breeding Blanket*, **Energies** 16(11):4430 (2023). https://www.mdpi.com/1996-1073/16/11/4430 (DOI: 10.3390/en16114430)
16. *Physics and Technology Research for Liquid-Metal Divertor Development (Tin-CPS) at the OLMAT High-Heat-Flux Facility*, **J. Fusion Energy** (2023). https://link.springer.com/article/10.1007/s10894-023-00373-9
17. *Lithium Divertor Targets and Walls for the ASTER Liquid Stellarator Reactor*, **J. Fusion Energy** (2025). https://link.springer.com/article/10.1007/s10894-025-00481-8
18. *Verification and Validation of COMSOL Magnetohydrodynamic Models for Liquid Metal Breeding Blankets Technologies*, **Energies** 14(17):5413 (2021). https://www.mdpi.com/1996-1073/14/17/5413
19. *Validation of COMSOL code for analyzing liquid metal magnetohydrodynamic flow*, **Fusion Eng. Des.** (2018). https://www.sciencedirect.com/science/article/abs/pii/S0920379618300115
20. Open FUSION Toolkit — documentation & module list (TokaMaker, ThinCurr, Marklin, MUG). https://openfusiontoolkit.github.io/OpenFUSIONToolkit/ ; GitHub: https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit
21. Open FUSION Toolkit, Columbia Fusion Research Center. https://fusion.columbia.edu/research-projects/open-fusion-toolkit-oft
22. C. Hansen et al., *TokaMaker: An open-source time-dependent Grad–Shafranov tool*, arXiv:2311.07719. https://arxiv.org/pdf/2311.07719

*Confidence notes:* §1, §2, §4 (COMSOL and MUG's stated purpose) rest on multiple converging sources and are **established**. The **absence** of validated 3D stellarator LM-MHD (§3) is an honest negative from a targeted search — the HELIAS-DCLL group's own "open challenge / geometry-workaround" framing is the strongest positive evidence that the gap is real. The FLiBe-specific conclusion in §5 is **inference from established σ-scaling**, not a paper that ran your exact case.
