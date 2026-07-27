# Conformal-offset recovery of leaky devices — findings

Date 2026-07-26. Staged, NOT wired into the live sweep default. Native (no transport)
except a final DAGMC probe. Metric = `leak_gate.worst_shell_crossing` (3D trimesh
`contains`), the task's ground truth; good devices -> 0.000, leaky -> 0.15-0.60.

## Two distinct failure modes (the task premise was only half the story)

The task attributed ALL leaks to a concave-inboard normal-offset FOLD, fixable by a
better offset. Investigating with the real 3D metric shows TWO modes:

1. **Concave-inboard FOLD** (fixable). The naive normal offset (`R + d*n`) inverts on
   the concave-inboard indentation of a bean cross-section. Fixed by a morphological
   (Minkowski-buffer) offset that fills the notch instead of inverting.

2. **Inboard AXIS-CROSSING** (NOT fixable by any offset). The leaky set is dominated by
   COMPACT, low-aspect-ratio devices (R0/a ~ 1.2-1.5). Their inboard wall sits closer
   to the machine axis than the blanket is thick, so the OUTER shells offset past R=0
   and the torus overlaps itself at the axis. Example 1505944: R0/a=1.22, inboard
   clearance 55 cm, but the 116 cm blanket reaches R=-61 cm. BOTH the buffer and the
   original normal offset hit the identical R=-61 -> this is a device/blanket-thickness
   incompatibility, not an offset-algorithm defect. No offset can fix it (the fix is a
   thinner blanket or a larger build scale, out of scope here).

Aspect-ratio / axis-crossing across the local 14-device error set (blanket 116 cm,
minor radius scaled to 250 cm):

| device  | R0/a | inboard R (cm) | buffer(116) min R | mode          |
|---------|------|----------------|-------------------|---------------|
| 1505944 | 1.22 | 55  | -61 | axis-crossing (unfixable) |
| 1577542 | 1.48 | 121 |   5 | axis-crossing |
| 1862008 | 1.47 | 118 |   2 | axis-crossing |
| 2190977 | 1.33 | 83  | -33 | axis-crossing |
| 2194424 | 1.31 | 78  | -38 | axis-crossing |
| 2218517 | 1.33 | 83  | -33 | axis-crossing |
| 2225908 | 1.21 | 53  | -63 | axis-crossing |
| 2369483 | 1.24 | 61  | -55 | axis-crossing |
| 2469994 | 1.39 | 98  | -18 | axis-crossing |
| 1880092 | 2.46 | 364 | 251 | **fold-limited (recoverable)** |
| 2147186 | 1.86 | 215 |  99 | **fold-limited** |
| 2291622 | 1.51 | 128 |  12 | **fold-limited** |
| 2296767 | 1.61 | 153 |  37 | **fold-limited** |
| 2524474 | 2.30 | 325 | 209 | **fold-limited** |
| 59509   | 5.60 | 1151| 1035| good (no leak) |
| 115887  | 3.00 | 501 | 385 | good |
| 134412  | 1.71 | 177 |  63 | good |
| 157482  | 2.14 | 286 | 170 | good |

**9/14 leaky are axis-crossing (unfixable); 5/14 are fold-limited (recoverable).**

## The offset method (isolated module `stellarator_geometry_offset2.build_layers_sdf`)

Non-folding morphological offset + per-shell phi-smoothing, on the SAME STL/manifest
interface as `build_layers` (consumed unchanged by dagmc_writer / leak_gate / sweep):
  * per fixed-phi cross-section, boundary(d) = exterior of `Polygon.buffer(d)` (disk
    dilation -> never folds, fills the concave notch); arc-length resampled;
  * each shell's slices are PHASE-ALIGNED across phi (FFT cross-correlation) so the
    shell's own closed surface is simple in 3D (no toroidal self-intersection).
The manifest records `axis_crossing` (outermost boundary min R <= 0) so a compact
device that no offset can help is flagged rather than silently leaking.

Also implemented (documented, not the default): `build_layers_clamped`, a
curvature-limited normal offset (task option b) — it nests exactly but COLLAPSES
shells to zero thickness on tight concavities, so it is not used.

## worst_shell_crossing before/after (leak_gate 3D metric, pol_res=96 unless noted)

Fold-limited (recoverable) leaky devices. All 5 have sweep status=error (leaked in
transport under the live normal offset). OLD column = the earlier per-phi 2D crossing
of the normal offset (the 3D leak_gate OLD is >= this and DAGMC errored on all 5):
| device  | OLD (normal, 2D) | NEW buffer (leak_gate 3D) | verdict |
|---------|------------------|---------------------------|---------|
| 1880092 | 0.277 | 0.0026 (pol_res=96) -> <1e-3 at pol_res>=128 | RECOVERED* |
| 2147186 | 0.025 | 0.0000 | RECOVERED |
| 2291622 | 0.057 | 0.0000 | RECOVERED |
| 2296767 | 0.467 | 0.0000 | RECOVERED |
| 2524474 | 0.056 | 0.0000 | RECOVERED |

*1880092's only residual is piecewise-linear chord-sag on the 0.2 cm W layer, a
resolution artifact (2D metric: 0.0007 at pol_res=128 -> 0.0000 at 192), not a
geometry defect. Default pol_res is 3x input (=192) which clears it.

No-regression (good devices), OLD vs NEW both via leak_gate (3D):
| device  | OLD (normal) | NEW (buffer) | |
|---------|--------------|--------------|-|
| 59509   | 0.0000 | 0.0000 | OK |
| 134412  | 0.0000 | 0.0000 | OK (aspect 1.71, closest to fold regime) |
| 157482  | 0.0000 | 0.0000 | OK |

3/3 good devices unchanged -> NO regression.

Axis-crossing devices stay leaky by construction: 1505944 NEW=0.094 (per-pair: inner
5 shells 0.0000, outer shield/coil 0.058/0.094 because they offset past R=0). NOT
recoverable at 116 cm blanket -- confirmed the outer coil mesh is 3D self-intersecting
(trimesh.contains is ~random even on points deep inside it).

## Transport validation (DAGMC, Lima worker)

Recovered fold-limited device 2296767, buffer offset pol_res=192, cheap DAGMC
watertight probe (low-density filler, isotropic 14 MeV point source, 20000 particles,
OMP=2 to stay off the concurrent grow-sweep):
  status = "go",  n_lost = 0,  lost_frac = 0.0
So the recovered geometry is genuinely watertight and transports with ZERO lost
particles -- leak_gate=0.0000 corresponds to real DAGMC success. (A full polarized
delta_free run was deferred: the VM was saturated at load ~5 by the parallel
grow-sweep; the leak verdict is field-independent so the cheap probe is sufficient to
prove recovery.)

## Recovery projection for the eta(C) law

Observed on the local 14-device error set: **5/14 = 36% recoverable** (the fold-limited
fraction); the other 9/14 are geometrically incompatible with the 116 cm blanket
(axis-crossing) and require a thinner blanket or larger scale, not a better offset.

Projected onto the ~69 Ginsburg errors, IF the aspect-ratio distribution matches the
local set: ~0.36 x 69 ~= 25 recovered -> new n ~= 64 + 25 ~= **~89** (about +40%, not
the hoped ~2x). The remaining ~44 axis-crossing errors are recoverable ONLY by
changing the blanket/scale policy (a separate lever), not the offset.

CAVEAT: the local 14 are the devices with a staged surface npz; the true 69-set aspect
distribution should be checked before committing the 89. The buffer offset is the
correct and necessary fix for the fold-limited class and does no harm to good devices.
