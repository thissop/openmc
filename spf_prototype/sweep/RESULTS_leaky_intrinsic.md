# Conformal-shell leaks are INTRINSIC, not a numerical artifact (defines the studiable-device boundary)

**Question.** Can we recover the ~38% of conformal-sweep devices that fail with
`status="error"` (OpenMC "Maximum lost particles / No intersection with DAGMC cell N")?

**Answer: no — the leaks are physical, tied to plasma shape, and higher resolution makes them worse.**
This is a real result: it defines which devices the free-streaming conformal method *can* study.

## Evidence (local worker, native diagnosis)
1. **The leak is real and mode-specific.** Device 1880092: `free_unpolarized` loses 0 particles, but
   `free_perpendicular` leaks ~thousands (~0.2% of histories) — neutrons emitted perpendicular to B̂ are
   aimed at the **tangled concave-inboard** region of the conformal shell. An honest lost-particle cap of
   100 (= 2.5×10⁻⁵ of 4M, negligible vs the ~0.2% coil-flux stat error) does NOT recover it; all 14 sampled
   error devices still abort → correctly stay `error`.
2. **Higher resolution folds, it doesn't seal.** Rebuilding at 96×144 and 128×192 makes the geometry
   **self-intersect (fold)** at every scale. The coarse 64×96 mesh was *smoothing over* a concave-inboard
   near-self-intersection; higher resolution resolves it as a hard fold. **leak → fold, not leak → watertight.**
3. **Root cause of the `error` (vs a clean gate):** `run_conformal.build_model` never set
   `max_lost_particles`, so OpenMC's default cap of 10 hard-aborts a leaky device before the sweep's own
   `lost_particle_max=50` warn logic runs. Fixed to be env-gated (`SPF_MAX_LOST` etc.), defaults preserve the
   exact prior behavior — but raising the cap does not honestly recover the devices (point 1).

## Consequence
- The **~36% usable yield is a physics boundary**, not a bug: tight concave-inboard plasma cross-sections
  cannot be given a watertight conformal offset shell. These should be a **clean `dagmc_leaky` reject**
  (documented device class), not an error.
- **Growing n past 64 requires FRESH devices** (new engineering-relevant QUASR configs outside the
  178-device `production_scaleup.json` manifest), NOT re-running the manifest (which is exactly
  {64 successes} + {114 known failures}). The local worker reproduces the failures deterministically.
- The audit-failers (separate ~19%) are marginal QUASR coil-field quality — the **own-coil design** path
  (design our own good coils) would remove that failure class; the leaky class is plasma-shape-intrinsic and
  own-coils won't help it (it's the wall shape, not the coils).

## Artifacts
- `sweep/recover_leaky.py` — the native leak-probe (per-mode lost-particle count + resolution sweep).
- `python/run_conformal.py` — env-gated lost-particle cap override (behavior-preserving default).
- `sweep/sweep_out_local/batch_classification_summary.csv` — the 34+ re-run classifications (0 new done, as expected).
- Combined-law harvest (`shield_opt/eta_C_zoo_harvest_local.py`) independently reproduced the committed
  law: **n=64, ρ=+0.464 — unchanged** (a validation of the CSV/analysis).
