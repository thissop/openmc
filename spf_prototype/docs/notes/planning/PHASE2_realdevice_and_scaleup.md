# Phase 2 — real-device √g ([C1] fix) + intelligent device scale-up

Status log for: (a) landing the QUASR-boundary→DESC `[C1]` fix so real *device* equilibria solve, and
(b) scaling the geometry→peaking correlation study to ~10× more devices spanning the S_φ range (plus a
QI arm). Deferred for now: (d) the scattering arm A(τ), and paper writing.

Newest status at the top of each section.

## (a) [C1] fix — real per-device √g

**Diagnosis (2026-07-07, `ginsburg_jobs/c1_diag.sbatch`).** The QUASR boundaries are NOT degenerate
(803097 axisym cross-section area 2.07; 886079 area 20.7; sensible R/Z). The failure is a **θ-winding
convention**: QUASR's LCFS parameterization winds θ *clockwise* (the m=1 R mode is negative, e.g.
803097 R₁=−0.085; the poloidal Jacobian is consistently negative), so DESC's
`ensure_positive_jacobian` hits a zero building the axisymmetric seed inside
`solve_continuation_automatic`.

**Fix.** `quasr_equilibrium_field._quasr_boundary_modes(..., flip_theta=True)` now evaluates the
boundary at **−θ** so the fitted Fourier modes wind counter-clockwise (positive Jacobian). This is a
pure reparameterization of the *same* physical boundary. Verification is the existing rigor gates in
`quasr_fluxmap.build`:
- **V4** (rho=1 reconstructs the fitted boundary, in the same flipped convention) — the key
  self-consistency check that the flip didn't distort the surface;
- **V1** (own-quadrature volume = DESC V), **V1b** (independent divergence-theorem volume), **V2**
  (single-signed √g). If all pass on 803097, the fix is correct and general.

**Status (2026-07-07): the theta-flip did NOT fix it — DROPPED as a rabbit hole.** Re-tested 803097
with `flip_theta=True`: identical `sqrt(g) == 0` failure in `ensure_positive_jacobian` during DESC's
axisymmetric-seed construction. So the cause is NOT theta-winding — DESC's default initial guess for
803097's boundary **self-intersects** (√g crosses zero near the axis, from the m≤8 boundary detail).
This is a genuine DESC fixed-boundary convergence problem (needs a supplied magnetic axis / reduced
boundary resolution / a different seed), with ~34-min queue-bound iterations. **Not worth grinding:**
the StellaratorSource is already validated on a REAL equilibrium (`precise_QA` — all rigor gates + the
transport ±40.6%/−21.1% capstone + the sampler on real √g), so a specific QUASR device would add
capability we've shown, not physics we haven't. If revisited later: pass a good axis guess to
`Equilibrium(...)`, or lower `mpol/ntor` in `_quasr_boundary_modes`, or use VMEC instead of DESC.
Meanwhile the real-equilibrium demonstration runs on `precise_QA` (`python/realdev_stell_demo.py`):
native StellaratorSource in transport on the real equilibrium (real √g births + real b̂), 3 modes.

## (b) Intelligent device scale-up (~10× more, spanning S_φ)

Goal: populate the `PF_perp/unpol vs S_φ` (and vs nfp / class) correlation graphs with ~100+ devices
that **uniformly span the S_φ range** and the QA/QH(/QI) classes — NOT all clustered at high-S_φ QA.

Pipeline (each step a script + a SLURM array where heavy):
1. **Candidate pool** from `data/quasr/catalogue.csv.gz` (has qs_error, nfp, aspect_ratio, helicity,
   mean_iota, elongation for the whole database). Stratify candidates by (helicity-class × nfp ×
   qs_error-bin × aspect-bin) so the pool already spans design space.
2. **Download serials** for the pool (`quasr_loader` fetches `serial<ID>.json` from
   quasr.flatironinstitute.org; we currently have only 12 local). Cache under `data/quasr/serials/`.
3. **Compute S_φ** for every candidate (cheap: coil Biot-Savart b̂ on the (1−ρ²) source volume →
   `coherence_metrics`; no transport). Parallel.
4. **Stratified select** ~100–120 devices to make the S_φ histogram roughly FLAT across [min,max]
   within each class (equal-count S_φ bins), avoiding S_φ clustering. Record the selection + why.
5. **Run `geometry_peaking`** on the selected set (analytic free-streaming, already validated) →
   refreshed correlation figures + stats with n≈100.
6. **QI arm:** QUASR has no QI (it's quasisymmetric-only). Generate a few QI equilibria in DESC
   (omnigenity/QI examples + published QI boundaries), run them through `quasr_fluxmap`→peaking with
   b̂ from the equilibrium, and drop QI points on the same S_φ plot → 3-class QA/QH/QI landscape.

Selection-bias guard: report the S_φ histogram of the final set (must be ~flat, not peaked); report
per-class counts; keep the reversal-device audit (`reversal_frac`) so 1190023-like configs are flagged.

**Status (2026-07-07): DONE for QA/QH (QI still pending DESC).** Pool of 354 stratified candidates →
S_phi computed for all (0 failures), spanning **S_φ ∈ [−0.264, 0.999]** → flat down-select to **83
devices** (26 QA, 57 QH; ~flat 8–16/bin) → parallel peaking run.
**Result — the scale-up REVISED the T3 finding (honest):** Spearman(S_φ, PF_perp/unpol) = **−0.30
(p=0.006)** at n=83, vs −0.83 at n=12. The strong small-n value was the QA↔QH cluster contrast. Class
separation is robust (Mann-Whitney p=0.001), but S_φ as a *continuous* predictor is weak — the peak
factor is geometry-dominated (the **G[M]** point). See `FINDINGS_geometry_peaking.md` revision box.
**Next:** (i) test **G[M]** (the geometry-folded contraction) as the predictor — it includes the wall
geometry S_φ lacks, so it should track PF far better; (ii) the **QI arm** once `[C1]`/DESC is back;
(iii) optionally expand the pool (larger N_PER_STRATUM) to fill the sparse low-S_φ-QA / high-S_φ-QH
bins toward n~120.
