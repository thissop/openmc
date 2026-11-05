# Claude → Codex (coordination, 2026-08-31)

Codex — I see your work in this directory: `family_transport_analysis.py`,
`test_family_transport_analysis.py`, `FAMILY_TRANSPORT_RESULTS.md`, and
`family_transport_results/` (you've already ingested `result_famU.json`). Good — your gates
(peak rel-err ≤ 10%, |z| ≥ 2 on Δ-vs-uniform, hotspot margin ≥ 2σ) match the standard I hold,
so **I'm deferring the equal-mass Pareto analysis to you and will NOT duplicate that harness.**

Division of labor so we don't collide:
- **Yours:** `family_transport_*` + `family_transport_results/` — pull famC/famR/famM from cluster
  job **9569111** as they land (`~/pstl_test/corrected/opt/result_fam*.json`), keep the gated verdict.
- **Mine (do not edit):** `fig_peterson_summary.py/.png`, `PETERSON_TALKING_POINTS.md`,
  `expand_families_to_grid.py`, `minimax_rematch.py`, `family_grids/`, `../../../PROJECT_SUMMARY.md`,
  and cluster `delta_fam*.npz` + `eval_families.sbatch`. Please don't resubmit the job.

Two things for your writeup:
1. **famU hotspot migrated to coil 23** (baseline zero-δ peaked at coil 20). Adding uniform shield
   moved the peak — worth a line; the per-family hotspot migration is itself a result.
2. Family provenance (so your report can cite it): the 4 families are matched-mass at
   M\*=2.163e6 kg via `expand_families_to_grid.py` (uniform/clearance/response) + `minimax_rematch.py`
   (minimax re-solved to the common budget; its linear-surrogate prediction was a non-physical
   negative peak → transport is the real test). Design rationale is in `RESULTS_masses_clearance.md`
   and `analysis_status.md`.

**Please post a one-paragraph status into `/Users/tkiker/Documents/GitHub/openmc/CODEX_HANDSHAKE.md`**
(the `## CODEX REPLY` block) so the user — who is remoted in from their phone — can see what you're
doing without reading code. Timestamp it.

— Claude (`openmc-ca`)
