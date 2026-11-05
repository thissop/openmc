# Peterson check-in — talking points (2026-08-31)

Figure: `fig_peterson_summary.png` (regenerate: `python fig_peterson_summary.py`). All numbers are from
gated artifacts; `verify_closeout.py` PASSES every gate.

Meeting deck: `../../../talk/peterson_checkin_2026-08-31.pptx`.

## Six-slide run of show

1. **Harder, Better, Faster, Stronger neutronics** — use the four words as an organizing device, not
   four independent novelty claims.
2. **One transport stack connects the four claims** — harder/better are source fidelity;
   faster/stronger are the direct-forward design loop and placement result.
3. **The shield-design thesis is quantified** — walk left-to-right through mass, response, clearance.
4. **Three ready results, three caveats** — say the caveat immediately after each number.
5. **Equal mass overturns the heuristics** — neither clearance nor response improves uniform; the
   extrapolating minimax stress test backfires by 891±36%.
6. **Advisor decisions** — choose the abstract/paper lead and the next proof target.

Preferred summer-project title:
**Harder, Better, Faster, Stronger: Physics-Informed Fusion Neutron Sources and Response-Aware
Magnet-Shield Optimization in OpenMC.** The 228-word final-evidence abstract and the passed
equal-mass result sentence are in `../../../ABSTRACT_DRAFT_2026-08-31.md`.

## One-line thesis
**Response-aware stellarator shield design, quantified:** *attribution ≠ intervention worth ≠ feasible
controllability.* Where you most need shield (inboard, coil-20 hotspot) is where you have the least mass
budget and the least clearance to grow.

## Three results ready to present (verified, no pending compute)

**A. Equal thickness ≠ equal mass** — 40M-sample MC VolumeCalculation on the real DAGMC.
- The *same* +20 cm shield patch costs **163–400 t depending on sector — a 2.45× spread**.
- Every patch is **lighter than the uniform +2.55 cm control** (529 t); t3p3 is 129 t lighter at 8σ.
- Refutes the earlier "t3p3 ≈ uniform at equal mass" guess. Use MC for absolute mass; ParaStell area×t is
  a shape-only proxy, +6–19% biased.

**B. Placement-dependent control** — 12-patch finite difference on the 4×4 sector grid.
- Inboard sectors **reduce** the peak coil load; wrong sectors **backfire** above baseline. **The sign
  pattern is robust.**
- Honest caveat (say it out loud): the absolute magnitude (best patch ~5× down) is **baseline-confounded**
  by the shield-outside offset (patch builds at 40 cm vs the 42.55 cm baseline), and the **hotspot
  migrates** off coil-20 (to c25/c30/c15) even for the "reducers." So: real steering exists; the headline
  number is not yet clean.

**C. Need collides with feasibility** — clearance field to the coil filaments.
- Available growth toward the coils dips to **4.7 cm inboard (poloidal ≈173°)** — exactly where coil-20
  peaks. Median room is ~290 cm; the bind is entirely inboard.
- This is the concrete engineering tension in one number.

## Equal-mass transport complete — the decisive negative result
- All four ParaStell/DAGMC families were evaluated at **M*=2.163×10⁶ kg** and 4M histories with the same
  direct-forward objective. Peak-value and hotspot-label gates both pass.
- Clearance-driven is **+0.2% ± 4.8% vs uniform (z=+0.05)** and response-aware is **+0.6% ± 4.9%
  (z=+0.12)**: both are statistically indistinguishable from uniform, with the resolved peak remaining at c23.
- The **linear-minimax stress test catastrophically backfires: +891% ± 36% (z=+88)**, a 9.91× peak, again at
  c23. Its linear surrogate had predicted a nonphysical negative response; transport falsifies that extrapolation.
- Meeting sentence: *"At fixed shield mass, neither clearance nor single-patch response guidance improves on
  uniform bulk, while unconstrained linear-minimax extrapolation worsens the transported peak by nearly 10×."*
- Interpretation: the data reject these one-shot placement heuristics, not transport-informed optimization in
  principle. A production optimizer needs bounded steps plus direct-forward relinearization/nonlinear training data.
- Integrity note: the aggressive minimax breeder layer split into seven CAD solids. The guarded exporter reconciled
  25 tags to 31 volumes, and the tally used H5M-derived magnet cells 25–32 with stable physical labels 19–26.

## Method note (de-risks the story)
- Coils are **not** variance-reduction-starved: the nominal high-flux 2M-history analog benchmark resolves all
  8 coils/period at 2–4% in ~6 min, so the optimizer uses **direct forward transport** — no VR/SSR needed for
  the peak coil metric. In the heavily shielded 4M `famU` run, the peak is still well resolved (3.5%), although
  only 5/8 low-to-high coils are below 10%; do not generalize the all-eight precision claim to every family.
  (The earlier "zero coils" was an SSR black-inside artifact, since fixed.)
