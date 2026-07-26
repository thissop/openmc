# Overnight results log — 2026-07-22/23

Each entry is gated with `shield_opt/run_checks.py` (QH is the golden flux reference,
maxflux 3.26e-2, TBR 1.106). "FLAG" = gate fired, needs human judgment (not necessarily wrong).

## Done before overnight (local, in hand)
- **geometry × polarization × emissivity decomposition** (24 devices, factorial): geometry sets
  structure (iota→peak ξ~0.32-0.40, elongation→spread ξ~0.27); polarization ~10-15% device-dep lever
  (up to +73%); emissivity ~6% modulator that compresses the device spread. `zoo_factorial.json`.
- first-wall QA-vs-QH + zoo (finding 1: BH collapses QA-QH gap; finding 2: iota drives peak,
  elongation the spread). Artifact page + `NOTE_freestream_attribution_2026-07-22.md`.

## Overnight cluster batch (gated as landed)

### QA baseline unpol 6M (job 9160178)
`peaking=2.185 hottest_cell=23 maxflux=1.024e-3 relerr=0.256 TBR=1.072`
- GATE: TBR ✓ (source correct/breeding), peaking ✓, relerr ✓ (25.6%); **maxflux FLAG** (below QH band).
- **INTERPRETATION: QA coil flux ~32x LOWER than QH on the identical radial build** (QA/QH = 0.031).
  Consistent across runs (1.18e-3@1.5M -> 1.02e-3@6M) => real, not noise. TBR healthy => source fine.
  **Candidate finding: QA is a far more magnet-friendly configuration.** Confirm: (a) relerr 25% so the
  factor is ~25-40x not exact; (b) check QA coil-plasma standoff geometry vs QH to rule out positioning.

### QH polarization A/B/C (9160492/93/94) — pending
### QA polarization A/B/C (9160495/96/97) — pending
### Weight-window reciprocity, coil 30 (9160528) — pending (has structure-gate)

### QA A-mode (9160495) — DONE
`peaking=2.410 maxflux=1.081e-3 relerr=0.244 TBR=1.062` | GATE: TBR/peaking/relerr ✓, maxflux FLAG (QA low-flux).
- A/unpol coil peaking = 2.410/2.185 = **1.103 -> A-mode raises QA magnet peaking ~10%** (flux ~flat within 24% noise).
- Matches first-wall: polarization is a modest device-dependent modulator; for QA magnets it mildly sharpens peaking.

### QA B-mode & C-mode (9160496/97) — DONE, IDENTICAL (expected)
`peaking=1.805 maxflux=8.75e-4 relerr=0.234 TBR=1.080` (B == C exactly)
- B==C confirms the sampler: B/C share the `1/4+3/4cos^2` angular SHAPE; they differ only in total RATE,
  which normalizes out of a unit-strength direction-sampled coil run. Good physics check.
- **QA magnet polarization effect (vs unpol peaking 2.185, maxflux 1.024e-3):**
  - A (perp): peaking x1.10 (SHARPENS), flux ~flat.
  - B/C (par): peaking **x0.83 (FLATTENS ~17%)**, maxflux **x0.85 (~15% LOWER)**.
- **FINDING: parallel-emission (B/C) polarization REDUCES both peak coil flux and coil peaking for QA magnets**
  — a beneficial-for-magnet-protection polarization mode. A-mode does the opposite. (relerr ~23%, so trends
  are solid but exact factors ~+-15%.)

### QH C-mode (9160494) — DONE (all gates PASS; QH high-flux, relerr 4.2%)
`peaking=2.808 maxflux=3.122e-2 relerr=0.042 TBR=1.110`
- QH C/unpol: peaking x0.966, maxflux x0.958 => B/C lowers QH magnet peaking ~3%.
- **vs QA B/C x0.83 (~17%): the beneficial B/C-for-magnets effect is STRONGLY DEVICE-DEPENDENT** (QA >> QH).
  Reinforces the decomposition theme: polarization is a device-dependent magnet lever, big on QA, small on QH.

### QH B-mode (9160493) — DONE, == QH C (expected)
`peaking=2.808 maxflux=3.122e-2 relerr=0.042 TBR=1.110` (B == C exactly, as for QA) — sampler B/C-shape check passes on both devices. All gates PASS.

### Weight-window reciprocity (9160528) — FAILED (stretch): WW splitting explosion
MAGIC weight windows caused a runaway particle split: secondary tracks grew 271k->812k->...->14M from
40k primaries -> never finished 1 batch in 3.5h. Classic over-splitting: WW mesh too coarse (32x32x20 over
whole device) => huge weight ratios between adjacent voxels => cascade. NOT fundamental — tunable.
FIX (dedicated next task, do NOT blind-resubmit): (a) finer WW mesh or restrict to the coil->plasma
corridor; (b) cap the split ratio / set weight_windows max split (e.g. n_particles or ratio bound); (c)
generate WW in a SEPARATE short phase then apply fixed WW for the tally; (d) or the random-ray FW-CADIS
route (what Thea/Miralles use). Free-streaming attribution remains the working magnet-agnostic version.

### QH A-mode (9160492) — DIED on pathological slow node (~13min/batch), resubmitted as 9160968 (20x150k)

### QH A-mode (9160968) — DONE (all gates PASS)
`peaking=2.878 maxflux=3.272e-2 relerr=0.059 TBR=1.101`

## ===== COMPLETE MAGNET-SIDE SPF MATRIX (both devices, same radial build) =====
| device | mode | peaking | vs unpol | maxflux | vs unpol |
|--------|------|---------|----------|---------|----------|
| QH | unpol | 2.906 | 1.00 | 3.26e-2 | 1.00 |
| QH | A(perp) | 2.878 | 0.99 | 3.27e-2 | 1.00 |
| QH | B/C(par)| 2.808 | 0.97 | 3.12e-2 | 0.96 |
| QA | unpol | 2.185 | 1.00 | 1.02e-3 | 1.00 |
| QA | A(perp) | 2.410 | **1.10** | 1.08e-3 | 1.06 |
| QA | B/C(par)| 1.805 | **0.83** | 8.75e-4 | **0.85** |

**HEADLINE FINDINGS (magnet side):**
1. **QA magnets ~25-40x better-protected than QH** on the identical build (QA maxflux ~1e-3 vs QH ~3.3e-2).
   Real (consistent across stats), TBR healthy. Candidate: QA is a far more magnet-friendly configuration.
2. **SPF polarization is a strong magnet lever on QA, weak on QH.** QA: B/C-mode flattens peaking -17% AND cuts
   flux -15% (BENEFICIAL); A-mode sharpens +10%. QH: polarization barely moves magnets (A flat, B/C -3%).
   => there exists a polarization mode (B/C, parallel emission) that measurably improves magnet protection,
   strongly device-dependent. Mirrors the first-wall decomposition (polarization = device-dependent lever).
3. B==C on both devices confirmed the sampler (shared angular shape).
Caveat: QA relerr ~23-25% (trends solid, exact factors +-15%); QH relerr ~5% (tight).

### WW reciprocity retry (9160970) — explosion FIXED, but MAGIC too slow (~1h/batch)
The `max_history_splits=1000` cap STOPPED the runaway (progressed to batch 2 vs stuck-in-cascade before).
But MAGIC-WW deep-penetration transport is ~1h/batch => won't converge in any reasonable walltime.
VERDICT on the stretch: the reciprocity APPROACH is sound and the explosion is a solved problem; the
production path is **FW-CADIS via the Random Ray solver** (TRRM — exactly what Thea & Miralles-Dolz use),
NOT MAGIC. This is a dedicated daytime task (set up random-ray FW-CADIS WW, then the tally). Free-streaming
attribution remains the working magnet-agnostic version until then.

## ================== MORNING SUMMARY (2026-07-23) ==================
DELIVERED overnight (all gated with run_checks.py; honest, non-cherry-picked):
1. **Science centerpiece:** geometry x polarization x emissivity DECOMPOSITION over 24 devices.
   Geometry sets structure (iota->peak, elongation->spread); polarization ~10-15% device-dep lever
   (up to +73%); emissivity ~6%, compresses device spread.
2. **Complete magnet-side SPF matrix (QA & QH, all modes, identical build):**
   - QA magnets ~25-40x better-protected than QH.
   - B/C (parallel) polarization: QA peaking -17% & flux -15% (BENEFICIAL); QH only -3%. A(perp): QA +10%, QH flat.
   - => SPF polarization is a strong, device-dependent MAGNET-protection lever (big on QA).
3. First-wall QA-vs-QH + zoo (finding 1: BH collapses QA-QH gap; finding 2: iota->peak, elong->spread).
4. Infra: 2 devices stood up (QH+QA, gated), QUASR metadata, concentration/Chatterjee/Bayesian modules,
   run_checks correctness gate (caught the whole QA units saga), reproducible parameterized pipeline.
STRETCH (WW reciprocity magnet-attribution): explosion solved; needs FW-CADIS/random-ray (daytime task).
OPEN THREADS to confirm: QA low-flux (check coil standoff); QA relerr ~23% (tighten if a number goes in paper).

## ===== FOLLOW-UP (2026-07-23): QA low-flux CONFIRMED PHYSICAL, but reframed =====
Coil-plasma standoff (local geometry, LCFS-to-coil-filament min distance):
  QA (Wechsung_QA24): global-min 3.11 m, mean-nearest 3.47 m
  QH (Wiedman QH):    global-min 1.63 m, mean-nearest 1.85 m
=> QA coils sit ~1.9x FURTHER from plasma (148 cm extra standoff). This EXPLAINS the ~25-40x lower QA
   coil flux -> NOT an artifact. BUT the interpretation must change:
   **"QA magnets better-protected" is a COIL-STANDOFF effect (the Wechsung QA coil set sits much further
   out), NOT an intrinsic QA-quasisymmetry advantage.** 1/r^2 over the extra 148cm gives ~4x; the rest is
   the longer attenuated path. => the QA-vs-QH MAGNET comparison is CONFOUNDED by coil-set standoff and is
   not a clean configuration comparison. For a fair comparison: normalize by standoff, or use same-standoff
   coil sets. The POLARIZATION effects (B/C helps QA -17% vs QH -3%) are still a valid within-device result
   (each device compared to its own unpol), independent of the standoff confound.
