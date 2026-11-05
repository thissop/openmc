# MEETING CONTEXT FOR AI ASSISTANT — Stellarator Neutronics Weekly Check-In (2026-08-19)

You are assisting me (Thaddaeus Kiker, PhD student) during a weekly research check-in with my
advisor, Ethan Peterson. Below is everything you need to answer his questions accurately. If he asks
something not covered here, say what is known and flag the uncertainty honestly — do NOT invent
numbers or overclaim. Prefer "that's an open question, my current read is…" over false confidence.

===============================================================================
## 1. ONE-PARAGRAPH SUMMARY
===============================================================================
I'm building a tool that uses Monte Carlo neutron transport to protect the superconducting magnets
of a stellarator fusion reactor. The idea: at a fixed machine size, redistribute the blanket — locally
trade tritium-breeder thickness for neutron-shield thickness — to cut the peak neutron load on the
coils, driven automatically by transport (not hand-tuning). The bottleneck is that a full
plasma→magnet Monte Carlo run is ~9 hours, too slow to repeat for every candidate design. This week's
work is about a method (surface-source decoupling) to make each design evaluation cheaper, a bug I
found in it, the fix, and an honest accounting of what the method actually saves.

===============================================================================
## 2. THE DEVICE
===============================================================================
- Quasi-helically symmetric (QH) stellarator, Landreman–Paul configuration.
- 4 field periods (nfp = 4). 20 modular coils. Reactor-scale, major radius R0 ≈ 13.7 m.
- It is a STELLARATOR, not a tokamak (no toroidal symmetry; the plasma cross-section twists as it
  goes around — that's why the neutronics must be fully 3-D).
- Geometry built with ParaStell (Moreno, Bader, Wilson 2024) as a DAGMC model.
- Radial build, 8 layers from plasma outward: first wall → Be multiplier → FLiBe breeder →
  back wall → tungsten-carbide (WC) shield → gap → vacuum vessel → thermal shield → magnets.
- In the DAGMC model: plasma void = cell 2; first wall = cell 3; … shield = cell 7; magnets = cells
  11–30. The adjoint importance used here targets coil 20.
- Neutron source: OpenMC StellaratorSource on the QH flux surfaces, monoenergetic 14.06 MeV,
  unpolarized. (Monoenergetic matters — see the fix below.)

===============================================================================
## 3. THE BIG GOAL (the thesis)
===============================================================================
At a FIXED outer envelope (coils don't move, machine size fixed), redistribute blanket material —
trade breeder thickness for shield thickness, LOCALLY, on a (toroidal φ, poloidal θ) grid — to
minimize the PEAK neutron/heat load on the magnets, inside an automated optimization loop
(change geometry → re-tally with transport → update → converge). Primary metric: peak coil fast-neutron
flux / nuclear heating. This is exactly the "machine-driven FWBS optimization" that the ParaStell
paper names as future work — I'm building that layer.

===============================================================================
## 4. PRIOR CONTEXT (how we got here — in case he references it)
===============================================================================
- Earlier I ran a "patch-worth" experiment: does a CHEAP adjoint proxy predict WHERE adding shield
  actually helps, patch by patch? Ground truth = full finite-difference transport per patch (~17 h each,
  12 patches). VERDICT: the cheap proxy does NOT rank the patches. It only reliably flags the one
  patch where adding shield BACKFIRES (raises coil dose, because removing breeder removes neutron
  moderation and releases fast flux — a multigroup effect a single-group proxy can't see).
- Because the cheap proxy can't rank, real transport must stay in the optimization loop. That is WHY
  we pivoted to the surface-source decoupling (SSR) — to keep Monte Carlo in the loop but make each
  iteration cheaper.

===============================================================================
## 5. THE METHOD: SURFACE-SOURCE DECOUPLING (SSR)
===============================================================================
- SSR = "surface source read/write," an OpenMC feature. WRITE: while a run transports particles, bank
  every particle that crosses a chosen surface (its position, direction, energy, weight) to a file.
  READ: start a later run by launching particles FROM that banked file instead of a physics source.
- Plan: run the expensive plasma transport ONCE, freeze the neutron field at the FIRST WALL, then for
  each candidate shield/breeder design, replay those frozen neutrons through only the (varying)
  blanket+shield out to the magnets. Skip re-sampling the plasma each time.
- Why freeze at the FIRST WALL specifically: the thing we vary (breeder AND shield) sits BEHIND the
  first wall, so the freeze surface must be in front of both to keep them free to change. (In OpenMC
  this is surf_source_write with cellfrom = 2, the plasma-void/first-wall boundary.)
- The validation test ("Gate 2"): run SSR on the SAME uniform build the full run used. If the frozen
  source is faithful, the per-coil ratio (full/SSR) is CONSTANT across coils (the constant is a benign
  source-normalization factor) and the hottest coil matches. Target CV < 5%.

===============================================================================
## 6. THIS WEEK: THE BUG, THE DIAGNOSIS, THE FIX
===============================================================================
THE BUG (double-counting):
- surf_source_write records a neutron EVERY time it crosses the first wall.
- The first wall wraps a torus → it is NON-CONVEX (the inboard wall faces itself across the donut
  hole). So a neutron can leave the wall, scatter in the blanket, cross back through the hole, and be
  recorded AGAIN somewhere else.
- Measured: 56.7 million banked crossings from 24 million source neutrons = 2.36 crossings per neutron.
- Composition: 42% are "fresh" uncollided plasma neutrons (still at exactly 14.06 MeV; 0.99 per birth =
  exactly one per neutron, as it should be); 58% are COLLIDED backscatter re-crossings (degraded below
  14.06 MeV).
- Effect: replaying all crossings as independent particles DOUBLE-COUNTS the backscatter population and
  distorts the coil-load pattern (worst on the far coils). Gate-2 CV = 48%.

FAILED FIRST FIX (worth mentioning if asked, shows the reasoning):
- Tried filling the plasma cavity with a B-10 absorber to kill re-entering neutrons ("void-kill").
- WRONG TOOL: B-10 is a THERMAL absorber, nearly transparent to the FAST neutrons that hit the coils,
  and it mostly SCATTERS them rather than absorbing. Result: CV only 48% → 47%. Essentially no help.

THE REAL FIX (pristine source):
- Keep ONLY the 42% fresh (uncollided, 14.06 MeV) crossings — the true plasma→wall current, which is
  build-independent — and leave the plasma void as a REAL VACUUM so re-crossing/void-transit happens
  naturally during the replay.
- This is rigorous: no double-count (backscatter never banked), no black-material hack, no energy bias.
- Reused the existing 4.5 GB master run (filtered the fresh neutrons out of it) — no re-run needed.
- I found the split cleanly BECAUSE the source is monoenergetic: anything at 14.06 MeV is fresh,
  anything below is a collided re-crossing.

===============================================================================
## 7. THE RESULT (be precise about these)
===============================================================================
Partial success — it works where it matters, with one open residual.
- INBOARD PEAK-LOAD COILS (cells 16–21): full/SSR ratio collapsed to k = 1.26, CV = 4.1% — FAITHFUL
  (passes the <5% gate for these coils). These carry the peak magnet load; they are what the optimizer
  targets.
- HOTTEST COIL: full run peaks at coil 20, SSR peaks at coil 19 — but they're within 0.6% of each other
  (a tie, not a real mismatch).
- Spearman rank correlation (full vs SSR across all coils) = +0.91.
- OPEN RESIDUAL: far/outboard coils (cells 26–30) are still under-predicted by 3.5–5×. This is REAL,
  not statistical noise (per-coil relative errors are 5–9%, but the deficit is 3–5×).
- GLOBAL CV (all 20 coils) = 47.4% — still "fails" the gate, but that number is dominated by the far
  coils. Mean k over all coils = 2.38.
- The deficit grows with "downstream-ness": coils directly behind the high-current inboard wall are
  faithful; coils reached only by indirect, scattered, void-crossing paths are progressively worse.

===============================================================================
## 8. THE "WHY IS IT STILL LOW ON THE FAR COILS?" ANSWER (say this)
===============================================================================
"In principle the replay should be exact — every fusion neutron streams straight through the empty
plasma to the first wall, so freezing that current and leaving the void empty should reproduce
everything downstream. It does on the inboard peak-load coils. The far/outboard coils are fed mostly
by neutrons that scatter in the inboard blanket and cross the plasma void to the far side — an indirect
channel — and that channel comes out systematically low. My leading explanation is that dropping the
collided re-crossings removed part of what feeds the far side, and the fresh neutrons don't fully
regenerate it. Two concrete next moves: tally the void-crossing current in the full run vs the replay
to localize exactly where the flux is lost, and test freezing at the neutron-exit layer (right before
the magnets), where the coil load is captured directly and should be exact for every coil."
- If asked "is it just statistics?": No. Per-coil errors are 5–9%; the deficit is 3–5×. Systematic.
- If asked "did you check for a re-crossing kill in the code?": Yes — OpenMC's surf_source_read does
  NOT terminate particles that re-cross the surface (confirmed the master banked 2.36 crossings/birth,
  so particles keep transporting). So that clean explanation is ruled out.

===============================================================================
## 9. THE SPEED-UP REALITY (important, and honest)
===============================================================================
- At EQUAL particle count (24 M histories), one SSR replay ran ~7h48m vs the full run's ~8h48m —
  only ~11–13% faster.
- Why so small: the replay still redoes ALL the breeder+shield scattering (the expensive part). It only
  skips the cheap plasma-source sampling and the empty-void streaming.
- So raw per-run speed is NOT the win. The real value of SSR here is:
  (a) ANY-RESOLUTION MESHING for free: the frozen surface source stores exact neutron positions, so you
      can bin the neutron field to any (φ,θ) grid — e.g. on the neutron-exit layer coming out of the
      shield onto the magnets — statistics-limited, not voxel-limited. (This is DEMONSTRATED.) The old
      3-D transport mesh was 30×30×20 with ~120 cm voxels — coarser than the entire ~130 cm radial
      build — so it couldn't resolve the ~40 cm shield band and capped the worth grid near 4×4. The
      frozen source removes that wall.
  (b) CORRELATED SAMPLING on differences: if every design replays the SAME frozen source, the common
      source noise cancels in the DIFFERENCE between designs, so each design needs far fewer histories
      to tell which is better. (Mechanism in hand; NOT yet quantified — don't claim a specific speedup.)

===============================================================================
## 10. OPEN DECISIONS I WANT PETERSON'S STEER ON
===============================================================================
1. Is "peak-load-faithful" enough to start the breeder-for-shield optimizer now (the objective IS the
   peak magnet load, and SSR nails the peak at CV 4%), or hold until the far-side residual is closed?
2. FREEZE-SURFACE CHOICE (the core tradeoff):
   - First wall (what we did): keeps the breeder trade live, but only ~13% raw saving + the far-side
     residual.
   - Blanket exit: cheaper real transport, but bakes the breeder in (can't vary it) — kills the trade.
   - Neutron-exit layer (right before magnets): magnet load captured directly, exact for every coil,
     ideal for the fine meshing/attribution — but bakes in the entire radial build.
3. Run the cross-build / correlated-sampling test now (replay through modified builds vs their
   finite-difference truth) to quantify the speedup and prove the optimizer premise?

===============================================================================
## 11. ANTICIPATED QUESTIONS + ANSWERS
===============================================================================
Q: Why freeze at the first wall and not deeper?
A: Because both breeder and shield vary, and they're behind the first wall; freezing deeper would fix
   the breeder and kill the trade we're studying. That's also why the raw time saving is small.

Q: Is the far-coil residual a showstopper?
A: Not for the peak-load objective — SSR is faithful there. It matters if the optimizer needs the full
   spatial pattern. It's localized and I have a concrete diagnostic + an exit-layer freeze that would
   sidestep it.

Q: Why is even the peak coil ratio 1.26 and not exactly 1?
A: There's an overall benign source-normalization factor (crossings vs births); the peak coils are the
   part closest to faithful. The residual above 1 and the far-side growth are what I'm chasing.

Q: What's the "58%" number?
A: 58% of the banked first-wall crossings were collided backscatter re-crossings (below 14.06 MeV),
   only 42% were fresh plasma neutrons. Dropping the 58% is the fix.

Q: Did the B-10 void-kill do anything?
A: Almost nothing (CV 48→47%). B-10 is a thermal absorber, transparent to fast neutrons — wrong tool.

Q: How confident are you in the geometry / is this a real reactor?
A: It's a clean, reactor-scale QH Landreman–Paul equilibrium built with ParaStell; materials are
   homogenized (fine for transport, not for activation). The method is device-agnostic; this is the
   test case.

Q: What resolution can you mesh the shield/magnet load at now?
A: Any resolution, from the frozen surface source (statistics-limited). The old 3-D transport mesh
   capped us near 4×4; the frozen source goes to 13×19, 60×90, whatever we want, for free.

Q: Is spin-polarized fusion part of this?
A: Separate thread (a directional-source idea that composes with the same adjoint magnet-importance
   machinery). Not part of this week's work.

Q: What's next concretely?
A: (1) Localize the residual with a void-crossing current tally. (2) Test the exit-layer freeze.
   (3) Quantify correlated sampling via a cross-build test. (4) Settle the freeze-surface choice.

===============================================================================
## 12. GLOSSARY (for the assistant)
===============================================================================
- SSR / surface source read-write: OpenMC feature to bank particles crossing a surface (write) and
  later launch from that bank (read).
- cellfrom = 2: bank neutrons leaving DAGMC cell 2 (the plasma void), i.e. crossing the first wall.
- Coil fast flux: neutron flux in the 0.1–25 MeV group inside a magnet cell — the magnet-damage proxy.
- Contributon (C = φ·ψ†): forward flux × adjoint coil-importance; how much coil load passes through a
  point. Attribution reads it in the plasma (where neutrons are born); worth reads it in the shield
  (where blocking helps).
- Adjoint / ψ†: importance function; solves the transport problem "backwards" from a coil response.
- TBR: tritium breeding ratio (breeder performance). NWL: neutron wall loading.
- Breeder (FLiBe): makes tritium AND moderates fast neutrons. Shield (WC): stops neutrons reaching the
  coils. The trade removes moderator when it adds shield — that coupling is why the cheap proxy failed.
- LCFS / flux surface / nfp / QH: last closed flux surface / nested magnetic surfaces / number of field
  periods (4) / quasi-helically symmetric (the optimization class of this stellarator).
- ParaStell: parametric tool that turns a VMEC plasma equilibrium + coil filaments into a DAGMC
  neutronics model. DAGMC: lets OpenMC transport directly on CAD surfaces.

===============================================================================
## 13. TONE / HOW TO HELP ME LIVE
===============================================================================
- Be concise and specific; give me the number or the one-line answer I can say out loud.
- If Peterson challenges a claim, help me concede honestly where the method is unproven (correlated
  sampling not yet quantified; far-side residual open; only peak coils validated) rather than bluff.
- The strong, honest framing of the week: "I diagnosed the failure to a measured mechanism, fixed the
  part that carries the peak load, learned the raw speed saving is small but the meshing + correlated-
  sampling value is real, and here are the design decisions I need your steer on."
