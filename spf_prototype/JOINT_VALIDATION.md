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

## Field convention is already identical (verified)

The paper's mean-field parameterization (Eq. B0)

    B̂₀ = (cosα sinβ cosφ − cosβ sinφ,  cosβ cosφ + cosα sinβ sinφ,  −sinα sinβ)

is **algebraically identical** to this repo's `AngledField`
(`src/spf_field.hpp`), i.e. B̂ = cosβ·φ̂ + sinβ(cosα·R̂ − sinα·ẑ), and that
convention was independently pinned to anarrima's angled kernels to ~2×10⁻¹⁴
(`RESULTS_tier8_angled.md`). So the (α,β) the analytic code uses and the (α,β) the
MC source uses are the same axis — no reconciliation needed for the *mean* field.

## The one thing to confirm with the analytic author

The field map carries the **full** B̂(x) including the field-period **ripple** b̃
(the n≠0 content), which is exactly what the source samples. The analytic series
treats the ripple perturbatively (its α₁,β₁ excursions, Eq. btilde / btilde2). For
the free-streaming cross-check to agree to the analytic's stated ~1% at ε_eff≈0.15,
both codes must use the **same ripple definition along the loop**. Open question:
is the analytic code's per-loop (α(φ),β(φ)) sampled from precise_QA in the **same
toroidal-angle convention** (geometric φ, same sign/origin) as the field map grid?
If the free-streaming residual shows a coherent phase shift, an α-sign or a φ-origin
offset between the two samplings is the first suspect (cf. `RESULTS_tier8_angled.md`
§7 ambiguity note).

## Where the check runs

Step 1 (free-streaming agreement) is the **smoke test** in `RUN_ON_GINSBURG.md`
§1(iii): a low-history near-void conformal run whose per-patch, per-mode wall load
is compared to the analytic NWL on precise_QA. Passing it gates the scattering and
scan runs. φ-resolved tallies (`run_conformal.py`, 32 toroidal bins) expose the
toroidal/poloidal structure the comparison needs.
