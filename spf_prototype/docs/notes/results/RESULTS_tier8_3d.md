# RESULTS — Tier 8 / Part B (3D-field-driven stellarator SPF: end-to-end)

> **Machinery validation with PUBLIC placeholder geometry — NOT a Helios physics claim.** Helios-class radial-build rig (R0=8 m, a=1.8 m, ≥1.2 m standoff; arXiv:2512.08027 / ARIES-CS); every device-specific value is `INJECT(helios)` (see DATA_NEEDED.md). The geometry is **axisymmetric** (radial-build rig); the 3-D physics is in the SOURCE (2-field-period field map B̂(x) + flux-surface birth). True 3-D shaping is DEFERRED.

Radial build (cavity→out, cm): W 0.1, steel 1.9, Be 5.0, FLiBe 50.0, shield 40.0, coil 15.0. FLiBe 65% ⁶Li. 2,000,000 histories/config. Modes (Bae naming): unpolarized=iso, perpendicular=A(∝sin²θ), parallel=B/C(∝1+3cos²θ).

## Anchor B — emission moments about the local 3-D field B̂(x)

⟨(u·B̂_local)²⟩ over 200,000 flux-surface births (B̂ from the QA map):

| mode | ⟨(u·B̂)²⟩ | analytic | Δ |
|---|---|---|---|
| unpolarized | 0.3328 | 0.3333 | -0.0005 |
| perpendicular | 0.2007 | 0.2000 | +0.0007 |
| parallel | 0.4668 | 0.4667 | +0.0001 |

→ the 3-D-field-driven emission has the correct angular distribution about the **local, spatially varying** B̂. **Anchor C**: B̂ genuinely varies over the source (std of B̂_z across births = 0.092 ≠ 0 ⇒ poloidal tilt + φ-dependence present; a pure-toroidal field would give 0).

## Anchor A — reduce-to-axisymmetric (thin-wall free-streaming)

Field-map TOROIDAL map vs trusted bmode=toroidal, inboard/outboard directionality (mode/iso−1):

| mode | inboard fmap | inboard bmode | outboard fmap | outboard bmode |
|---|---|---|---|---|
| perpendicular | +19.5% | +20.3% | -16.1% | -16.1% |
| parallel | -19.8% | -20.4% | +16.3% | +16.2% |

→ max |field-map − bmode| = **0.8 pts** (limited by the field-map grid resolution + MC stats; the field-map interpolation reproduces the analytic toroidal steering).

## End-to-end scattering run (per source neutron; rate held fixed)

| metric | unpolarized | perpendicular (A) | parallel (B/C) |
|---|---|---|---|
| global TBR | 1.368±0.001 | 1.366±0.001 | 1.373±0.001 |
| inboard first-wall load (rel iso) | 0.0% | +2.3% | -2.2% |
| outboard first-wall load (rel iso) | 0.0% | -7.0% | +7.2% |
| coil fast flux >0.1 MeV [/cm²/src] | 5.35e-04 (±24%) | 5.83e-04 (±23%) | 4.68e-04 (±30%) |

**Directional efficiency η (resolved, first-wall):** parallel reduces the inboard load by **2.2%** with the 3-D QA field vs **2.4%** with the axisymmetric toroidal field ⇒ **η = 0.89** (fraction of the axisymmetric steering retained under the pitched 3-D field). Coil fast-flux directional-η is left to a variance-reduced run (the deep-coil tally is 30% stat. error here; DEFERRED).

## Anchors / honesty

- **Reduce-to-axisymmetric** (Anchor A) and **emission moments** (Anchor B) pass → the field-map source + 3-D B̂ plumbing is correct.

- **Geometry is axisymmetric** (radial-build rig); no field-period/QA *geometry* effect is claimed. **Coil precision** needs variance reduction (MAGIC, as Bae 2025). **Single energy** (14.1 MeV). See DEFERRED.md.

**Tier-8/Part-B gate: the pipeline runs end-to-end (unpolarized/parallel/perpendicular) on a 3-D-field-driven flux-surface source through a Helios-class radial build, anchored by reduce-to-axisymmetric + emission moments, and produces TBR, first-wall steering, coil flux, and η.**
