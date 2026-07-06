# Legacy: the precise_QA single-equilibrium Ginsburg workflow (SUPERSEDED)

These are the ORIGINAL Ginsburg operational docs for the single-equilibrium **precise_QA**
conformal run (Tier-8), **superseded by the QUASR quasisymmetry sweep**. Kept for
reference, not deleted.

- `SETUP_GINSBURG.md` — offline SLURM runbook (env + storage + login-node rule + the
  precise_QA `ginsburg_job.sh` submit + the QA→QH scan). Its still-current ENV / STORAGE /
  login-node FACTS have been folded into the sweep runbook `sweep/docs/README.md`.
- `RUN_ON_GINSBURG.md` — the conceptual precise_QA pipeline (desc_to_fieldmap →
  stellarator_geometry → build_dagmc (pymoab) → run_conformal).
- `GINSBURG_CLAUDE_HANDOFF.md` — the (stale) mission brief for a Claude-on-Ginsburg session
  to smoke-test the precise_QA package.
- `ginsburg_job.sh` — the precise_QA sbatch (`run_ginsburg.py --stem equil_precise_qa`).

## What to use instead (the QUASR sweep — the current work)

- **Runbook:** `spf_prototype/sweep/docs/README.md` (self-contained).
- **SLURM:** `spf_prototype/sweep/sweep.sbatch`.
- **DAGMC:** `python/dagmc_writer.py` (pymoab-free), NOT `build_dagmc`/`stl_to_h5m`.
- **Env preflight (still current):** `spf_prototype/ginsburg_preflight.sh` (at repo root).

> NOTE: the CODE these docs drive — `python/run_ginsburg.py`, `python/run_conformal.py`,
> `python/build_dagmc.py` — is NOT archived; it stays in `python/` because the sweep
> reuses `run_ginsburg._mat_score` and `run_conformal.build_model`.
