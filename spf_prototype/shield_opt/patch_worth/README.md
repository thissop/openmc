# Patch-worth experiment: attribution vs shielding actionability (QH)

**The decisive experiment** (per the ChatGPT red-team §15 + our prior-art convergence). It tests the
one claim that turns Paper B from "a technique" into a result: **does the adjoint *attribution* (which
plasma region drives the coil) actually predict the *shielding worth* (where adding shield reduces the
coil response per unit mass)?** If they differ, "attribution &ne; actionability" is the novel finding.
If they agree, attribution is a validated cheap proxy. Either outcome is publishable.

## The clean idea
The response-contribution density is the **contributon** `C(x) = phi(x) * psi_dagger(x)` (forward flux ×
adjoint flux), and it is defined *everywhere*, not just in the plasma. So:
- **Attribution** `A(theta,phi)` = the contributon in the **plasma** (`q * psi_dagger`), projected to the
  shield surface = "where were the responsible neutrons born." (This is exactly our current `placed`
  logic via `adjoint_placement.placement_priority`.)
- **Worth** `W(theta,phi)` = the contributon in the **shield** (`phi * psi_dagger` summed over the shield
  radial band per angular column) = "where does blocking a neutron actually help." Times the added-material
  removal cross section, this is the first-order sensitivity `-dR/dm` (adjoint shielding worth).

They live on the same `(theta,phi)` grid but come from different regions; transport moves the contributon
from plasma to shield along streaming paths, so `A` and `W` need not agree. **The experiment compares
their rankings on `(theta,phi)` patches, and validates `W` against finite-difference truth.**

CAVEAT (red-team §5.1): the scalar product `phi*psi_dagger` is the *leading-order* worth; for anisotropic
scattering the exact derivative needs angular fluxes/moments (this is where Jack's angular tally becomes
relevant to Paper B too). So `W` here is a scalar heuristic that we **validate against finite-difference**
patch worth — we do not claim it is exact.

## Geometry / response (reuse what we have)
- QH corrected 8-layer DAGMC (`dagmc_corr_w35f` / the step1 build), coil-20 adjoint-peak.
- Response `R` = **peak-coil (cell 20) fast flux**, same as the kill-shot. (Extend to heating/dpa later.)
- `(theta,phi)` grid = the build grid (13 toroidal × 19 poloidal) coarsened to **~16 patches** (4×4).
- Adjoint field `psi_dagger`: `~/adjoint_test/recip_w35f_coil20_flat/adjoint_importance_flat_P0.npz`
  (key `importance`, mesh 30×30×20 over the full bounding box) — HAVE.
- Plasma source `q`: `~/qh_freestream_fluxmap.npz` (uniform emissivity, to match the forward source) — HAVE.

## Protocol (tiered — cheap first, expensive only to validate)

### Tier 0 — one new forward run (prerequisite, cheap)
Run the forward source on the baseline QH DAGMC with a **regular mesh flux tally matching the adjoint mesh**
(30×30×20, same bounds) so `phi(x)` and `psi_dagger(x)` are co-located. -> `qh_fwd_meshflux.npz`.
(`slurm_patch_worth.sh` step 1.) One forward run with the FW-CADIS weight windows -> fast.

### Tier 1 — the cheap comparison (no new geometry)
`patch_worth.py`:
1. Define ~16 `(theta,phi)` patches about the magnetic axis R0.
2. `A(theta,phi)` = `placement_priority(contributon(coil20_map, fluxmap, uniform))` binned to patches
   (attribution).
3. `W(theta,phi)` = sum of `phi*psi_dagger` over the **shield radial band** (radial distance from LCFS in
   [breeder_out, shield_out] &asymp; [0.59, 0.99] m) per patch (worth, scalar leading-order).
4. Rank-compare `A` vs `W`: **Spearman rho + Chatterjee xi**, plus a sign/top-k overlap. Report the
   **attribution-actionability gap** (how different the two rankings are). -> `patch_worth_tier1.csv` + figure.

**This alone answers the headline** (do attribution and worth rank patches the same?) with zero new DAGMC.

### Tier 2 — finite-difference validation (expensive, small subset)
For ~6-8 patches spanning {high-A/low-W, low-A/high-W, high-both, low-both} (chosen from Tier 1):
1. `build_patch.py` builds a QH DAGMC with a fixed extra shield `delta` at that single patch (breeder
   thinned to hold the envelope, exactly the kill-shot trade, but localized to one patch).
2. `slurm_patch_worth.sh` runs the peak-coil dose (with FW-CADIS weight windows) -> `R_j`.
3. Finite-difference worth `W_FD_j = -(R_j - R_baseline) / delta_mass_j`.
4. `patch_worth_analyze.py`: correlate `W_FD` vs `W` (validate the adjoint worth) and vs `A` (test
   attribution&ne;worth) — Spearman + sign accuracy + relative-magnitude error. -> the decisive table + figure.

**Statistical rigor** (red-team §13): correlated sampling (same seed across patches), independent-seed
verification of the baseline + the top design, per-run relerr on the peak coil reported. `delta` chosen so
the expected `Delta R` exceeds the per-bin relerr (check against the kill-shot relerr first).

## Outcomes (all publishable)
- `W_FD` &asymp; `W` and `W` &ne; `A`  -> **attribution is NOT a good placement proxy; worth-guided is the lever** (strong novel result).
- `W_FD` &asymp; `W` &asymp; `A`  -> attribution is a validated cheap proxy (also useful).
- `W_FD` &ne; `W`  -> the scalar worth is inadequate; need angular (Jack) or finite-difference ranking (a real finding about the method).

## Files
- `patch_worth.py` — Tier 1 (attribution vs scalar worth on patches).
- `build_patch.py` — per-patch single-patch delta DAGMC (Tier 2).
- `slurm_patch_worth.sh` — Tier-0 forward mesh-flux run + Tier-2 per-patch coil doses.
- `patch_worth_analyze.py` — final W_FD vs W vs A correlations + the decisive figure.

## Cost
Tier 0: 1 forward run. Tier 1: seconds (fields we have). Tier 2: ~6-8 DAGMC builds + weight-windowed coil
doses. That is the whole expense — small vs a continuous shape derivative, and it settles whether the
"attribution -> actionability" thesis is real before any zoo-scale work.
