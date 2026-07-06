# PLAN — next phase: scattering, a realistic source, and a simple blanket

Status: proposal for review. Grounded in an OpenMC-API + data-coverage scoping
pass (all required nuclides, tally scores, and anarríma angled-field kernels
confirmed present — see "Feasibility", below).

## The organizing problem: verification after the analytic anchor fades

Everything so far was trusted because OpenMC could be checked against an **exact,
free-streaming analytic answer** (Schwartz / anarríma). The moment we add
scattering, **that closed-form ground truth disappears**. So the central design
rule of this phase is:

> **Every new physics step must ship with a replacement verification** — a
> limiting case that recovers something already validated, a conservation law, a
> code-to-code or benchmark comparison, or a physical-sanity signature — *not*
> "the plot looks reasonable."

Two axes of realism, deliberately separated by how well we can still verify:

| Axis | Loses analytic anchor? | Replacement verification |
|---|---|---|
| **Source realism** (energy spectrum, pitched/real field) | **No** — anarríma's angled kernels keep a free-streaming analytic answer | direct recovery test (as in Tier 2b) |
| **Transport realism** (scattering wall → heating/damage → blanket/TBR) | **Yes** | limiting cases + conservation + benchmarks + spectrum sanity |

So we do the **well-anchored source work first** (high confidence, and it's the
bridge to Peterson's source), then the **transport work** with new checks.

## Feasibility (confirmed this pass)

- **Data:** all of Fe/Cr/Ni/W/C/Si (structure), Li6/Li7/Be9/F19/Pb/O16 (breeder/
  multiplier) are in our NNDC library with elastic-scattering data. **No new data
  download needed.**
- **Tally scores valid in our OpenMC:** `H3-production` (tritium), `heating` /
  `heating-local`, `damage-energy`, plus `flux` with an `EnergyFilter`.
- **Realistic energy:** `openmc.stats.fusion_neutron_spectrum(ion_temp,
  reactants='DT')` gives a Ballabio-broadened ~14.1 MeV Gaussian.
- **Angled field analytic anchor:** anarríma exposes `g_HAa, g_Hca, g_HBa, g_VAa,
  g_Vca, g_VBa` (args include pitch angles α, β). Our sampler's Gram-Schmidt
  rotation already emits correctly for **any** B̂ — so a pitched field keeps an
  exact free-streaming cross-check.

---

## Tier 4 — A realistic source (still free-streaming; analytic anchor preserved)

Cheap, high-confidence, and the bridge to a physics-informed plasma source.

**4a. Birth energy spectrum.** Replace the monoenergetic 14.1 MeV with a
Ballabio-broadened Gaussian sampled inside the compiled source (Box–Muller on
`openmc::prn`), parameters from `fusion_neutron_spectrum(ion_temp,'DT')`.
- *Verify:* (i) C++↔Python parity of the energy sampler on fixed seeds; (ii) the
  sampled spectrum matches `openmc.stats` mean/σ; (iii) the **geometric NWL
  pattern is unchanged** (energy is decoupled from direction) — a clean
  regression against Tier 2b.

**4b. Angled / spatially-varying field (the Peterson bridge).** Feed a pitched
field B̂(position) (toroidal + poloidal, i.e. α,β) into the source; extend
`analytic_nwl.py` to call anarríma's `g_*a` angled kernels.
- *Verify:* the **same recovery test as Tier 2b**, but for the pitched field —
  OpenMC (free-streaming) must reproduce the angled-field analytic. Because the
  anchor survives, we keep full confidence.
- *Payoff:* this is exactly the interface a real equilibrium field plugs into →
  the natural coupling point to **Peterson's physics-informed tokamak source**.

**4c. (Optional, larger) Realistic plasma shape.** Sample birth positions from a
D-shaped / equilibrium last-closed-flux-surface instead of `1−ρ²` rings, with the
local B̂ from the equilibrium. Still anchored by 4b's angled kernels.

*Effort:* 4a small; 4b small–moderate (sampler ready, kernels exist); 4c moderate.

---

## Tier 5 — Scattering on: a real first wall (the first OpenMC-only physics)

The first results the analytic theory **structurally cannot** produce — and the
first honest step toward an actual wall dose.

**5a. Real wall material.** Swap the void first wall for tungsten or a steel
(`set_density` + `add_nuclide`/`add_element`, `fill=material`). Keep the plasma
and exterior as void.
- *Verify (the anchor):* the **"sandwich test"** — run (1) void (current Tier
  2b), (2) the same wall at density ×10⁻⁴ (mean free path ≫ device), (3)
  full density. Step (2) **must recover the free-streaming NWL within MC
  statistics**; step (3) then departs monotonically. This is how we keep a
  verification anchor with scattering on.

**5b. The new, scattering-only observables.** Add tallies that were meaningless
in the void: wall `heating` (→ MW/m²), `damage-energy` (→ dpa proxy), and an
energy-resolved `flux` spectrum at the wall.
- *Verify (no analytic):* (i) the **down-scattered tail** below the 14.1 MeV peak
  must appear (and vanish in the thin-wall limit) — the direct signature that
  scattering is active; (ii) neutron balance / leakage conservation; (iii)
  optional **code-to-code** (OpenMC vs Serpent/MCNP) on the same model.
- *Science question:* does the spin-polarization steering **survive scattering**
  and still reduce real center-stack *heating/damage* (not just geometric
  current)?

*Effort:* moderate. Mechanically easy; the work is the verification strategy.

---

## Tier 6 — A simple breeder blanket and TBR (the headline reactor metric)

A new physics tier — highest payoff, hardest to verify.

**6a. Blanket geometry + material.** Add an annular breeder layer behind the
first wall: Li-ceramic or Li/FLiBe, Li-6 enriched
(`add_element('Li',1.0,enrichment=60,enrichment_target='Li6')`), optionally a
Be/Pb neutron multiplier.

**6b. TBR.** Tally `H3-production` over the breeder (`CellFilter`/
`MaterialFilter`); in fixed-source mode the result is already **tritium per
source neutron = TBR**. Report TBR per polarization mode.
- *Science question:* does steering neutrons outboard (B/C) **change the TBR** —
  turning our geometric "inboard/outboard fraction" precursor into a real number?
- *Verify (no analytic):* (i) reproduce a **textbook/benchmark TBR** for a simple
  Li sphere or slab (code-to-code or against literature) before trusting the
  torus number; (ii) neutron **conservation** (leakage + absorption + (n,2n)
  balance); (iii) **sensitivity** behaves correctly (TBR rises with Li-6
  enrichment and with a multiplier); (iv) cross-check that the breeding hot-spots
  line up with the Tier-3 geometric current map.

*Effort:* new physics tier (most work). Scoping TBR only — no depletion, single
zone, no tritium-extraction/coolant detail.

---

## Cross-cutting scope (state in every result, as before)

Even after this phase: single-zone scoping blanket; no depletion/burnup; photon
transport only if `heating-local` is used; idealized geometry (square or simple
D-shape, no divertor); polarization still a uniform input; neutron channel only.

## Effort / value summary

| Tier | What | Effort | Analytic anchor? | Headline value |
|---|---|---|---|---|
| 4a | Ballabio energy | small | yes (regression) | realism, cheap |
| 4b | angled/real field | small–moderate | **yes (g_\*a)** | **Peterson bridge** |
| 4c | D-shape plasma | moderate | yes | realistic source |
| 5 | scattering + wall heating/dpa/spectrum | moderate | limiting case | **first dose-relevant result** |
| 6 | blanket + TBR | large (new tier) | benchmarks only | **the reactor metric** |

## Recommended order
**4b → 4a → 5 → 6.** Lead with the angled-field bridge (cheap, keeps the anchor,
connects to the PI's source), fold in energy broadening, then scattering (first
new physics), then the blanket/TBR (the big one).

## Open questions for you
1. **How far this phase?** Stop at scattering + wall heating (Tier 5), or commit
   to the blanket/TBR (Tier 6)?
2. **Peterson's source now or later?** Should Tier 4b target a generic pitched
   field, or do you already have access to his equilibrium/field so we couple
   directly?
3. **Wall material of interest** (tungsten vs a steel vs FLiBe-facing) — drives
   the Tier-5 numbers.
