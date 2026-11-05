# RUNG 4 — physical metric + intervention definitions (locked before any computation)

Per the standing integrity rule: define exact quantities now; each computation below gets a preflight and uses
real geometry/volumes, not proxies. This is critical-path (engineering ladder), done locally while the cluster
transport queue clears; the compute steps flagged [CLUSTER] wait for allocation but do not block the definitions.

## A. Response metrics (winding-pack / magnet) — distinguish integral vs density
The current `coil_fast_flux` (CellFilter flux) is a CELL-INTEGRATED track length, NOT a density. Definitions:

1. **Cell-integrated fast flux** F_j = integral over coil-j cell of phi(E in [0.1,25]MeV) dV  [units cm; what we have].
2. **Coil-average fast flux density** <phi>_j = F_j / V_j  [1/cm^2]. Requires the ACTUAL coil cell volume V_j
   from an openmc.VolumeCalculation on the DAGMC (stochastic volume), NOT a formula. [CLUSTER]
3. **Peak local fast flux density** phi_peak,j = max over a fine mesh covering coil-j of the flux density —
   the binding magnet constraint (insulation/HTS damage at the hot spot). Requires a mesh tally over each coil
   region + mesh element volumes. [CLUSTER]
4. **Nuclear heating density** q'''_j (W/cm^3) = (heating-local integral)/V or mesh-resolved peak — the cryo /
   quench constraint. `coil_heating` (heating-local, CellFilter) is integrated -> divide by V or mesh-resolve. [CLUSTER]
5. **TBR** (global) — design output/cost, already tallied per config.
6. **Blanket/shield mass & max local build** — engineering cost (see C).

Reporting rule: never call F_j a "density". Report integral and density separately; the paper's magnet metric is
the peak local density (3)/(4), with the cell-integrated flux as a fast proxy whose relationship to density is
verified (F_j vs <phi>_j across coils once V_j is computed).

## B. Distinct intervention definitions — NEVER combine as one "intervention"
Each is a different physical design and must be compared separately (equal cost computed from real geometry):

- **(I1) shield<->breeder at FIXED envelope**: t_shield+t_breeder const; swap material only. (the coarse patches, DOE)
- **(I2) shield added into available VOID/clearance**: total build grows into the plasma-coil gap; breeder unchanged.
- **(I3) total build growth, coils fixed**: envelope grows inward toward the (fixed) coils; reduces clearance.
- **(I4) coil displacement / machine growth**: coils move out; changes device size.
- **(I5) breeder-only thickness control**: TBR lever, ~no shielding change.

The coarse patch study is I1 (fixed-envelope swap). The placement-vs-bulk and Pareto questions require I2/I3 with
mass and clearance computed exactly. Do NOT reuse an I1 result to claim an I2/I3 conclusion.

## C. "Equal cost" definitions — computed, not assumed (retracts the earlier cell.cm claim)
The earlier "t3p3 (500 cell.cm) vs uniform (630 cell.cm) at equal mass" was WRONG (cell.cm != mass). Exact costs:

- **Added shield mass** = sum over changed cells of (rho_shield - rho_replaced) * V_cell, with rho_shield(WC)=15.6,
  rho_breeder(FLiBe)=1.94 g/cc, and V_cell = the ACTUAL 3D shell volume of that (theta,phi) cell in the twisted
  geometry (from ParaStell surface area x radial thickness, or a VolumeCalculation). [CLUSTER/PARASTELL]
- **Max local build** = max over (theta,phi) of the total plasma->coil-facing build thickness.
- **Machine-size penalty** = change in the outer envelope / coil standoff (for I3/I4).
- **TBR** = tallied.

An "equal-mass" or "equal-cost" comparison is only valid once these are computed from the real geometries+densities.
Until then, placement-vs-bulk is reported as PRELIMINARY (a shield ADDITION, not equal-mass).

## D. What is blocked on the cluster vs doable now
- DEFINITIONS (this file): done, local.
- Coil/mesh VOLUMES (metrics 2-4), added-mass (C): need openmc.VolumeCalculation on the DAGMC. [CLUSTER, critical-path]
- Exact clearance field C(theta,phi) (Rung 6): min distance blanket-surface -> coil exclusion; needs ParaStell
  blanket surface + coils_qh. [CLUSTER/PARASTELL, critical-path]
- Rung-2 exact importance + projection: SECONDARY, queued opportunistically (job 9489689, backfill walltime).

Next critical-path compute when allocation frees: a single VolumeCalculation run (coil + shield-sublayer + breeder
cell volumes on dagmc_qh_step1_uniform) unblocks metrics (2-4) and the exact added-mass -> enables the real
placement-vs-bulk (I1) and the equal-cost baselines. This is cheap (~minutes) and high-value; prioritize it over the
secondary importance job when a slot opens.
