# Rung 2 (exact) — PREFLIGHT before projecting adjoint importance onto the shield surface

Per the standing integrity rule (2026-08-24): no torus/cylinder/nearest-coordinate proxy. Use the exact
ParaStell/VMEC geometry that built the transport model, or STOP and report the missing metadata.

## 1. Source quantity + physical meaning
Coil-20 ADJOINT IMPORTANCE psi_dagger(r), flat (energy-independent) response. Physical meaning: the expected
contribution of a neutron at position r to the coil-20 fast-flux response. NOT contributon (no S(r) weight),
NOT material worth (no scattering/GPT term), NOT forward-adjoint overlap. It is one scalar field for ONE response.

## 2. Source coords / units / shape / normalization
File adjoint_coil20_flat_P0.npz: importance shape (30,30,20) on a CARTESIAN mesh, lower_left=(-1800,-1800,-550),
upper_right=(1800,1800,550) cm -> voxels 120x120x55 cm in (x,y,z). Units: adjoint flux (relative). coil_cells=[20].
MUST VERIFY: which geometry/build this adjoint was computed on (step1_uniform full-torus?), and coil-20's period.

## 3. Target quantity + physical meaning
A per-patch scalar that should PREDICT the empirical coil-20 worth: the mean coil-20 importance over the exact
3D SHIELD SHELL VOLUME that each patch (i_tor, j_pol) converts breeder->shield (+20cm). High importance in a
patched shell => shielding there should reduce coil-20 (worth), if attribution predicts placement.

## 4. Target coords / units / shape / normalization
Patches: ParaStell (13 toroidal x 19 poloidal) grid, coarsened 4x4; nfp=4 replicated. Each patch = +20cm shield
over a (theta_pol, phi_tor) sector, in the SHIELD radial layer (between back_wall and gap). The shield shell in
the REAL twisted VMEC geometry -- points (theta_pol, phi_tor, radial) -> (x,y,z) via the ParaStell build off
wout_qh.nc at wall_s=1.08 + radial-build offset. ONLY the coil-20-period image of the (nfp-replicated) patch
contributes to coil-20's (single-location) importance.

## 5. Exact transform / correspondence
Need (theta_pol, phi_tor) -> (x,y,z) on the SHIELD mid-surface, in the SAME convention ParaStell used to build
the DAGMC (so coordinates are consistent with the adjoint mesh). Preferred: run ParaStell's own
construct_invessel_build (same call as build_step1.py) and extract the shield-layer surface point cloud per
(theta,phi). This guarantees coordinate consistency (same code path as the transport geometry). Alternative
(only if ParaStell exposes it cleanly): read wout_qh.nc R_mnc/Z_mns, evaluate the wall_s surface + ParaStell's
exact normal-offset to the shield -- but replicating ParaStell's offset convention is itself an approximation
to be avoided unless verified against ParaStell output.

## 6. Checks that MUST pass
- GEOMETRY-OVERLAY: each projected (x,y,z) must fall inside a SHIELD DAGMC cell (openmc.lib.find_cell on the
  actual dagmc_qh_step1_uniform.h5m). If points land in breeder/gap/void, the map is WRONG -> BLOCK.
- COORDINATE CONSISTENCY: ParaStell toroidal_angles (0..90, one period) and poloidal_angles (0..360) map to the
  adjoint mesh's (x,y,z) with the correct period offset for coil 20.
- RECONSTRUCTION: the set of patched (theta,phi) cells must reproduce the delta-field footprint in delta_qh_patch_*.
- DIMENSIONAL: importance averaged over shell voxels -> scalar per patch (consistent).

## 7. Approximation introduced
Target: NONE. Use ParaStell's exact surface. If the exact surface cannot be extracted from the repo, STOP and
report exactly what is missing (do NOT fall back to the torus proxy that produced the anti-correlated Rung-2 result).

## 8. Validation
The geometry-overlay find_cell check (point-in-shield-cell) is the pass/fail gate. Only if >~90% of projected
points land in shield cells is the projection trusted; then compute correlations vs empirical worth with the
same metrics as before, and report with the overlay diagnostic figure.
