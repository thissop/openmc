# JOINT_VALIDATION.md — analytic ↔ Monte-Carlo cross-check protocol

How the companion analytic paper (`docs/Analytic_Neutron_Wall_Loading_for_Spin_
Polarized_Plasmas_in_Quasisymmetric_Stellarators_by_Field_Period_Perturbation`)
and this OpenMC SPF pipeline validate each other, per the paper's §9/§11:

> "Because the analytic kernel is exact in the free-streaming limit, it is the
> quantity a Monte Carlo neutronics calculation with a polarized source must
> reproduce, per wall patch and per polarization mode, in its uncollided
> (first-flight) tally."

So the joint result is a three-step ladder:
1. **Free-streaming agreement** — OpenMC (near-void) reproduces the analytic
   free-streaming NWL, per patch, per mode (iso / A / B-C). Validates the source
   sampling, the polarized angular distribution, and the line-of-sight geometry.
2. **Scattering departure** — turn materials on; quantify how much the collided
   wall load departs from free-streaming (the part the analytic cannot give).
3. **Quasisymmetry scan** — sweep QA→QH equilibria; correlate the steering benefit
   / power / loading with the (non-)axisymmetry, the headline science question.

## What makes it a *joint* validation (shared inputs)

Both codes consume the **same equilibrium** — `precise_QA` (Landreman–Paul 2022;
NFP=2, aspect 6, R0≈1.03 m, a≈0.17 m). This is deliberate: precise_QA geometrically
matches the paper's §11 model QA test case (NFP=2, R0≈1, a≈0.17), so the analytic
series and the MC run are the *same configuration*, not two lookalikes.

- **Analytic side**: reads precise_QA's boundary harmonics + field angles (α,β) +
  field-direction ripple directly (the paper's Eqs. 6–7, B0, btilde).
- **MC side**: `desc_to_fieldmap.py` exports, from the *same* precise_QA:
  - `equil_precise_qa.{meta,bin}` — the field map B̂(x) (the source's steering axis),
  - `equil_precise_qa_surface.npz` — the LCFS surface grid (the conformal wall),
  - `equil_precise_qa_boundary.txt` — the boundary Fourier modes (reference).

## Conventions — audited against anarrima, confirmed to match (June 2026)

The analytic-side author audited the anarrima repo; both conventions match this
repo's implementation, verified in code (no reconciliation needed).

**Geometry / φ.** φ is the ordinary geometric cylindrical toroidal angle,
right-handed about +ẑ: `r_s(φ) = (R(φ)cosφ, R(φ)sinφ, Z(φ))`. Boundary uses the
simsopt/VMEC sign `R ~ cos(mθ − n·Nfp·φ)`, `Z ~ sin(mθ − n·Nfp·φ)`, so **φ=0 is a
stellarator-symmetry plane** and **θ=0 is the outboard midplane** (max R, Z=0).
- MC side matches: the source and `analytic_nwl.py` both use
  `r_s = (R cosφ, R sinφ, Z)` (polarized_fusion_source.cpp:138-154,
  analytic_nwl.py:87-88). My DESC precise_QA surface is up-down symmetric at φ=0
  (verified earlier), i.e. φ=0 is a symmetry plane — consistent with the sign above.
- If φ ever disagrees: left-handed φ → set `n → −n` (≡ sample φ → −φ); origin offset
  Δφ₀ → sample at `φ − Δφ₀`. Inherited automatically if both build geometry from the
  *same* precise_QA `SurfaceRZFourier`.

**Polarization.** K(cosθ) with `cosθ = B̂·Δ̂` — B̂ the LOCAL unit field at the
source, Δ̂ the unit line of sight source→wall (θ = angle between the outgoing
neutron and the local field). Modes: isotropic, A = sin²θ, B/C = ¼+¾cos²θ; physical
mix `G_phys = ¾·a·G_A + (⅔·b + ⅓·c)·G_B`, a+b+c=1.
- MC side matches exactly (verified in code): `analytic_nwl.FACTORS` A=`1−cos²θ`,
  BC=`¼+¾cos²θ`, `bracket = ¾a·G_A + (⅔b+⅓c)·G_BC` (analytic_nwl.py:33-35,289-294);
  the emission kernel I(θ)=¾a sin²θ+(⅔b+⅓c)(¼+¾cos²θ) (spf_sampler.hpp:15-20); and
  the source emits `u` about the **local** B̂ = `field_->bhat({px,py,pz})`, so
  `cosθ = u·B̂` with u the line of sight to the wall
  (polarized_fusion_source.cpp:159-162). The mean-field (α,β) is also identical to
  anarrima's angled kernels, pinned to ~2×10⁻¹⁴ (`RESULTS_tier8_angled.md`).

## CAVEAT — the numbers-level cross-check is gated on the analytic side

The conventions agree, but a **numbers-level** precise_QA comparison is **not yet
possible**: the analytic Tier-3 validation currently runs on a 6-mode toy spectrum +
an analytic *model* field, NOT precise_QA. anarrima itself speaks the real-data
convention; only its driver substitutes toy producers. A numbers cross-check needs
BOTH sides loading the **same real precise_QA** (same R_mn/Z_mn) and, on the analytic
side, sampling the **real vacuum field** at `r_s(φ)`. Until that swap is done, do NOT
treat the analytic Tier-3 outputs as precise_QA ground truth.

Same-source requirement: the MC field map + geometry here are built from **DESC**'s
precise_QA; the analytic side speaks **simsopt/VMEC**. For the numbers comparison,
both should load the *same* boundary (ideally simsopt `input.LandremanPaul2021_QA`)
— DESC's re-solved precise_QA may differ slightly from simsopt's. (Adding a
simsopt-boundary producer path is the clean way to guarantee a shared source.)

## Where the check runs (when both sides are on precise_QA)

`run_conformal.py` / `run_ginsburg.py` produce the free-streaming statepoints +
φ-resolved tallies; the per-patch, per-mode comparison to the analytic NWL is the
postprocessing step (not yet automated — needs the analytic side on precise_QA).
Practical matching:
- **Key the comparison on each patch's (R, Z, n̂), not on index.**
- The analytic driver currently sweeps **poloidal θ_w in [0,2π) at a single toroidal
  cut** (no toroidal φ_w sweep yet), so compare at that one cut for now; the MC
  φ-resolved tally already has the full toroidal×poloidal map for when the analytic
  side adds a φ_w sweep.
