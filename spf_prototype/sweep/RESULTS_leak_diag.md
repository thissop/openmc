# DAGMC leaky-yield diagnosis (native, no transport)

Diagnosed 2026-07-26. Tools: `leak_diag.py` (probe), `leak_gate.py` (staged fix).
Ground truth = `trimesh` watertightness + a shell-crossing metric. No OpenMC / no VM.

## The leak level: NOT faceting, NOT boundary-watertightness

For every leaky device, each conformal boundary surface is **individually perfect**:
`is_watertight=True`, `is_winding_consistent=True`, 0 naked edges, 0 non-manifold
edges, 0 degenerate faces. The sweep's own gate agrees (manifest: watertight=True,
simple_cross_section=True, badphi=0 on all layers) -- yet these devices leak.

The real defect is **adjacent nested SHELLS crossing each other**. Metric =
`worst_shell_crossing` = max over adjacent boundary pairs of the fraction of the
inner boundary's vertices lying OUTSIDE the next boundary out (trimesh `contains`).

| device   | status (Ginsburg) | worst_shell_crossing @64x96 |
|----------|-------------------|-----------------------------|
| 59509    | done              | **0.000** |
| 115887   | done              | **0.000** |
| 134412   | done              | **0.000** |
| 157482   | done              | **0.000** |
| 1505944  | error (leaky)     | **0.376** |
| 2190977  | error (leaky)     | **0.324** |

The metric separates good from leaky perfectly (0.0 vs 0.15-0.38). Cause: the naive
poloidal normal-offset (`stellarator_geometry.poloidal_outward_normals` +
`offset_surface`) inverts direction on the concave-inboard indentation of a bean
cross-section, so inner shells poke through outer shells. `_poly_is_simple` (the
existing gate) only tests single-polygon self-crossing, never shell-vs-shell, so it
is blind to this.

## Neither resolution nor scale recovers them

- **Scale** (native, `leak_diag` scale sweep): inflating the build scale 4x (device
  minor radius ~250 -> ~1000 cm, blanket fixed 116 cm) drops worst_crossing only
  0.376 -> 0.169 (1505944) and 0.324 -> 0.159 (2190977). Monotone but nowhere near
  sealed -> it's an offset-DIRECTION defect, not a magnitude one.
- **Resolution**: higher poloidal/toroidal density does the OPPOSITE of sealing. The
  coarse 64x96 mesh SMOOTHS the fold (individually watertight boundaries); at 96x144 /
  128x192 the original per-polygon self-intersection gate trips instead (leak -> hard
  fold). Corroborated by the parallel batch agent (1880092: perpendicular-mode lost
  particles persist past an honest lost-particle cap of 100).

**Recovery yield ~= 0.** Tight concave-inboard shapes are intrinsically leaky for
this builder.

## The fix (staged, not wired into the running sweep): clean reject gate

`leak_gate.worst_shell_crossing(geom_dir)` -- cheap, native, no transport. Drop it in
after `dagmc_writer.build_from_stls` and before `run_transport`:
if crossing > tol (1e-3) set `status="dagmc_leaky"` (a clean reject, NOT a transport
"error") and skip transport. See the module docstring for the exact snippet. This
converts the ~38% silent transport failures into an explicit pre-transport verdict and
stops wasting Lima/Ginsburg time on non-convergent meshes. An optional bounded
scale-retry (`build_geometry_watertight`) is included but recovers ~0.

## Yield impact on the eta(C) law

Errors are NOT recoverable geometry -- they are legitimately leaky and must be
rejected, not transported. So:
- Committed law n = 64 stays **n = 64** (no reclaimed devices).
- The honest sweep success fraction stays ~36% (Ginsburg: 68 done / 69 error / 34
  audit_failed / 12 folded = 68/183 ~ 37%).
- Real gain is operational, not statistical: the gate removes 69+ wasted transport
  attempts and reclassifies them cleanly, so future batches spend cores only on
  meshes that can converge. To grow n, the lever is a better conformal offset
  (curvature-limited / signed-distance morphological offset) or a larger QUASR draw,
  NOT resolution.
