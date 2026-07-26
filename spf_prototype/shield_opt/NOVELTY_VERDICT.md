# Novelty verdict — closed-loop breeder↔shield magnet-protection optimizer

*Deep adversarial lit search, 2026-07-14 (three parallel web-research passes). Verdict below;
"[FETCHED]" = read from primary text, "[SEARCH]" = search snippets/abstract only.*

## VERDICT

- **The CONCEPT is NOT novel.** Locally thinning the tritium breeder and substituting a
  high-efficiency shield at the minimum plasma–coil standoff, to fit + protect the
  superconducting magnet inside a fixed radial envelope, was done by **ARIES-CS (2008)** —
  by hand. Do **not** claim the idea; cite it.
- **The AUTOMATED, SPATIAL-HOTSPOT-DRIVEN CLOSED LOOP appears UNPUBLISHED.** No single work
  occupies all three of: (i) *spatially-resolved* (poloidal+toroidal) breeder↔shield
  allocation, (ii) at a *fixed outer envelope*, (iii) in a *closed automated loop driven by
  magnet dose/hotspot feedback*. That intersection is the defensible contribution.
- **3D stellarator source-attribution (the "culprit map") also appears UNPUBLISHED.** No paper
  inverts wall/coil load back to plasma source regions in 3D (closest is a *forward* NWL map).
- **Novelty is real but NARROW and TIME-SENSITIVE.** Frame it as *method/automation*, not
  concept. The biggest risk is an imminent Moreno/Wilson (UW-Madison CNERG) publication.

## The three anchor prior-arts (cite all three, differentiate carefully)

1. **ARIES-CS — the CONCEPT anchor.** El-Guebaly et al., UWFDM-1336 / *Fusion Sci. Technol.*
   54 (2008); Najmabadi et al., *FST* 54(3):655, DOI 10.13182/FST54-655. [FETCHED] Thinned the
   blanket and inserted WC shield at the 4 Δmin hotspots per field period (~24% of FW area) to
   fit breeder+shield and protect the SC magnet at fixed envelope; NWL 3.3 vs 5.3 MW/m². Done by
   **engineering judgment + iterative 1-D/3-D MCNP checks** — NOT automated, NOT hotspot-driven
   in a loop. *Differs from you only in automation.* This is your #1 prior art.
2. **ParaStell — the TOOLING anchor.** Moreno, Bader, Wilson, *Front. Nucl. Eng.* 3, 1384788
   (2024), DOI 10.3389/fnuen.2024.1384788. [FETCHED] 3D t(φ,θ) variable build, fixed envelope,
   breeder↔shield tradeoff, magnet-heating+TBR in OpenMC — but a **manual 99-config sweep** with
   **global** magnet totals, and it **explicitly names the closed loop as future work**
   ("coupling ParaStell to machine-driven optimization is planned"). *Differs: no automation, no
   spatial-hotspot driver.*
3. **FERMI/Tracer — the AUTOMATION anchor.** Borowiec, Bae, Badalassi, *Fusion Eng. Des.* 200,
   114159 (2024), DOI 10.1016/j.fusengdes.2024.114159. [SEARCH] A genuine automated closed-loop
   neutronics optimizer — but **ARC tokamak**, perturbs **uniform layer thicknesses**, objective
   is **TBR/cost**, not magnet-dose-driven, not a spatial field. *Differs: not spatial, not
   stellarator, not magnet-protection.* (Note: J.W. Bae = also the polarized-source-in-OpenMC
   author — this group is active in your exact tooling space.)

## Other close hits (all differ materially)
- **ARC / Kuang et al., *Fus. Eng. Des.* 137:221 (2018)** — tokamak analog: dedicated inboard
  shield behind FLiBe to protect the TF leg. Hand-designed, binary inboard/outboard, not a loop.
- **HELIAS (Häußler/Warmer 2017; Palermo *Nucl. Fusion* 61:076019 2021)** — thickness varies only
  to space-fill to the coil envelope; not a deliberate breeder↔shield trade. Distinguishable.
- **Lion, Warmer, Wang, *Nucl. Fusion* 62:076040 (2022)** — fast 3D NWL, but *forward*
  (source→wall) and optimizes *wall shape*, not source attribution and not shield. Nearest to the
  culprit-map idea but different direction.
- **Scalar/GA/Bayesian blanket optimizers** (Humphrey/SLEDO *Front. Nucl. Eng.* 2025; GA shield
  studies; TBR ANN surrogates; FUSE arXiv:2409.05894) — all optimize **global/scalar** metrics
  (TBR, 1-D stack), none a spatial magnet-dose-driven shield field.
- **Bae et al., *Nucl. Fusion* 65:086051 (2025)** — closest "polarized source in OpenMC" (forward
  TBR sensitivity to polarization); no shield opt, no attribution. (Frame SPF vs this carefully.)
- **EU-DEMO / CFETR inboard-shield-vs-breeder scans** — establish the *tension* but are manual
  1-D/3-D parameter scans, per-side uniform, no loop.

## Priority-defeaters (ranked) — check before any submission/claim
1. **Connor Moreno PhD thesis (UW-Madison / Wilson CNERG)** — his dissertation is literally this
   loop ("numerical optimization... modification of breeding blanket and shielding systems").
   In-progress, no paper/preprint/thesis found yet. **Highest risk. Recheck ProQuest/ANS/ISFNT
   2025–26 right before you claim priority.**
2. A **ParaStell/CNERG follow-on** paper closing the loop (TOFE/ISFNT proceedings).
3. **FERMI group (Bae/Badalassi/ORNL)** extending Tracer to stellarators or magnet-driven objectives.

## Recommended framing (honest + defensible)
> "ARIES-CS established the concept of a spatial breeder-for-shield trade to protect magnets at
> fixed envelope by hand; ParaStell built the 3D-variable-thickness geometry engine but evaluates
> it with manual sweeps and defers automation to future work; automated fusion-neutronics
> optimizers (FERMI, SLEDO) target global scalar metrics. We contribute the **automated closed
> loop that maps a spatially-resolved magnet-hotspot tally to an optimal t_shield(θ,φ) field**,
> plus the **source-region attribution** that identifies which plasma regions drive the load —
> neither of which is in the literature."

## Due diligence to do before a formal claim
Pull full text on: ARIES-CS *FST* journal version, Kuang 2018, Borowiec 2024; and re-search for a
2025–26 Moreno thesis / CNERG proceedings. Those are the items search snippets couldn't fully verify.
