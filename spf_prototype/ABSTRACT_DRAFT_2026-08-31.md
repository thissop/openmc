# Abstract/title draft for Peterson discussion — 2026-08-31

## Recommended title

**Harder, Better, Faster, Stronger: Physics-Informed Fusion Neutron Sources and Response-Aware Magnet-Shield Optimization in OpenMC**

## Final-evidence abstract

Protecting superconducting magnets from 14 MeV fusion neutrons requires source models that represent the plasma physics and shield designs judged by transport, mass, and engineering feasibility. We extend OpenMC with two composable capabilities that pose a **harder** problem than axisymmetric source models: a native `StellaratorSource` that samples fusion births from three-dimensional VMEC/DESC flux maps with volumetric weighting `S(ρ)|√g|` and local magnetic-field direction, and a source-agnostic deuterium–tritium spin-polarization kernel implementing the Schwartz angular distribution. Together they provide a **better** account of where neutrons are born and directed within a common Monte Carlo transport stack. We then build a **faster** response-aware shield-design workflow on periodic QUASR/ParaStell/DAGMC geometry; direct-forward transport resolves the peak coil-load objective without surface-source replay or production weight windows. The resulting design picture is **stronger**, but not simply thicker. At fixed outer envelope, attribution-guided placement reduced peak coil fast-neutron load by 4.2× relative to uniform, although the limiting coil migrated. In a controlled equal-mass comparison at 2.163×10⁶ kg, however, clearance- and single-response-guided layouts were statistically indistinguishable from uniform, while an extrapolating linear-minimax candidate drove the peak to 9.91±0.36 times uniform, falsifying its surrogate. Independent 40-million-sample volume calculations show that identical +20 cm patches add 163–400 tonnes depending on location, a 2.45× spread, while available growth falls to 4.7 cm inboard—where shielding is most needed. These results separate attribution, intervention worth, and feasible controllability, and motivate bounded transport-in-the-loop optimization.

## Claim boundary for the meeting

- The 4.2× comparison is **fixed envelope**, not equal mass.
- The equal-mass comparison is complete: clearance and response are unresolved versus uniform; the
  extrapolating linear-minimax stress test is 9.91±0.36× uniform and decisively worse.
- The stress-test result falsifies the one-shot linear surrogate; it is not evidence that a physically
  constrained, sequential transport-in-the-loop optimizer cannot improve the design.
- The conformal spin-polarization smoke study supports the **direction** of the coil-flux change, not a precise device-scale percentage.
- The stellarator birth density is `S(ρ)|√g|`; describing it as merely “proportional to `√g`” drops the emissivity profile.

## Equal-mass result sentence (gate passed)

> At a common shield mass of 2.163×10⁶ kg, neither clearance- nor single-response-guided placement changed the transported peak coil fast-neutron load significantly relative to uniform bulk, while an extrapolating linear-minimax candidate increased it by 891±36%; the limiting coil remained c23.
