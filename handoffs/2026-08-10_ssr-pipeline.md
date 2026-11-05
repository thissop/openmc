# Session Handoff — SSR Decoupled Shielding Pipeline

Written at the end of the 2026-08-10 session (after the Peterson meeting). Read this first next
session to pick up cleanly. Companion memory: `spf-ssr-decoupled-pipeline`, `spf-patch-worth-n12-verdict`,
`spf-killshot-placed-vs-uniform`.

---

## 1. Where we are in one paragraph

We finished the patch-worth experiment (does a cheap adjoint proxy predict where shielding helps),
and the honest answer at full resolution is **no** — neither attribution, scalar worth, nor signed
worth ranks the 12 measured patches; the signed worth only reliably catches the one patch where
shielding *backfires*. After the meeting we pivoted to Peterson's architecture: **decouple plasma
transport from shield optimization** so we can keep Monte Carlo in the loop but run it cheaply at high
resolution. That pipeline is now built, tested, and queued on Ginsburg.

## 2. Why we pivoted (the logic, not just the what)

- The proxies fail because they are single-group and ignore that the breeder **moderates** fast
  neutrons. Trading breeder for shield removes moderator and can *release* fast flux to the coil —
  a sign flip no total-cross-section proxy can see. (This is the n=12 verdict.)
- If the cheap proxy can't rank, you must evaluate placements with real transport. But a full
  plasma→magnet run per optimizer iteration is far too slow for high resolution.
- Peterson's fix: run **one** expensive transport, freeze the neutron phase space at a boundary
  (`surf_source_write`), then re-launch from that frozen source through only the part you're changing
  (`surf_source_read`). Now each optimizer iteration is a cheap shielding-only run, and you can re-bin
  the frozen source into any patch resolution offline for free.

## 3. Teaching — the four ideas you need fluent

**Contributon / worth.** `C = φ·ψ†` (forward flux × adjoint coil-importance) is how much of the coil
load passes through a point. Attribution reads it in the plasma (where neutrons are born); worth reads
it in the shield (where blocking helps). Measured worth `W_FD = −Δ(coil dose)/Δ(material)` is the
finite-difference truth. The experiment tested whether the free proxies predict W_FD. They don't rank.

**Why the proxy fails (multigroup moderation).** First-order theory says adding attenuator always
lowers dose. But the coil sees *fast* flux, and the breeder downscatters fast neutrons. Remove breeder
→ less moderation → more fast flux out. On one patch that beats the added shield and dose rises +41%.
Getting the sign right needs fast-vs-thermal treatment, i.e. real transport, not a scalar proxy.

**The SSR decoupling.** Freeze the inward neutron current at the **first wall** (not the blanket
exit). The first-wall source is set by the plasma and is independent of the radial build, so SSR from
there through the *variable* blanket+shield captures the breeder moderation and skips the expensive
plasma sampling. Freezing at the blanket exit instead would bake the breeder in and lose exactly the
physics that matters. And do **not** use an analytic `exp(−Σx)` shortcut — it ignores scattering, the
same error that sank the proxy.

**Collective, not per-coil.** The coil load is dominated by a shared inboard region every coil passes
through, not by 20 independent footprints. So the real optimization is low-dimensional (shape the
inboard band), which is what makes an automatic optimizer tractable.

## 4. What is queued right now (check these first)

**UPDATE 2026-08-11:** smoke chain PASSED as a plumbing test overnight (ran end-to-end, files written,
hottest coil matched cell 18, Spearman +0.72). Its CV=112% "FAIL" is **statistics only** — 40k
histories can't populate 20 coils, so far coils get 0 counts → NaN ratios. Not the backscatter
double-count; ignore the smoke Gate-2 verdict.

**Burst preemption bit us:** the real_master (`9358863`, 48M histories) ran ~10h to batch 7/12 then got
**preempted and restarted from zero — twice** (`Restarts=2`). Burst is opportunistic/preemptible and the
run needs ~17h uninterrupted, which it never got. **Fix applied:** moved the real chain to the
non-preemptible `short` queue (12h cap) and cut to **2M/batch = 24M total** each (still 600× smoke) so
master (~8-9h) + dose (~7h) fit 12h. Both `slurm_ssw_master.sh` and `slurm_ssr.sh` were edited
`-p burst`→`-p short`, `-t 1-00:00:00`→`-t 12:00:00`. New submit script: `ssr/submit_ssr_real_short.sh`.

Current real chain (on `short`, non-preemptible):
- `9372052` real_master → `9372053` real_gate2 (afterok). Submitted 2026-08-11 ~16:45 EDT.
- When both finish, run Gate 2: `python3 ssr/gate2_analyze.py coil_real_master.npz coil_real_gate2.npz`.
- If Gate 2 is statistics-marginal (not a physics fail), scale batches/particles back up on `short`
  via a second chained job or accept the 24M verdict.

Superseded (do not use): `9358861/62/63/64` (old burst chain, cancelled).

**UPDATE 2026-08-12 — real Gate 2 ran (24M) and FAILED as anticipated; void-kill fix implemented + rerunning.**
Full chain finished clean on `short`, no preemption (master 8h48m, dose 6h56m). Verdict (`GATE2_REAL_VERDICT.txt`):
CV=**48%** (need <5%), Spearman **+0.90**, hottest coil full=20 vs ssr=19 (0.6% apart = tie). NOT statistics:
600× more particles than smoke only moved CV 112%→48%, so ~48% is a systematic floor. Distortion is
geometric — full/ssr ratio ~2 at hot inboard coils (16–22) but 4–9 at cold coils. Cause (confirmed in
code): `surf_source_write cellfrom=2` banks EVERY plasma-cavity crossing incl. blanket-backscatter
re-exits; the SSR replay lets those stream back across the empty cavity and land decorrelated on the far
wall — the exact "plasma-void re-entry" the handoff pre-flagged.
FIX (in `ssr/ssr_dose.py`, `VOID_KILL=1` default): fill DAGMC cell 2 (plasma cavity) with a B-10 absorber
(2.2 g/cc) so replayed neutrons heading back in are killed; their real far-wall contribution is still carried
by the separately-banked re-exit particle. Master run + its 4.7GB `surf_real_master.h5` are VALID and reused;
only the ~7h dose reruns.

**Two false starts before it landed (important API gotchas, all caught by a hard-assert so NO compute wasted):**
1. First tried `dag._material_overrides={2:[void_abs]}` — validated on the LOGIN node's conda openmc (0.15.3),
   but the ACTUAL run uses the custom transport build at `$SRC` (`PYTHONPATH`+`LD_PRELOAD`, also self-reports
   0.15.3 but different code). `model.run()`'s export path there dropped the dict: `model.xml` had an empty
   `<dagmc_universe/>` and the run was a byte-for-byte no-op vs the pre-fix result. LESSON: probe the `$SRC`
   openmc, not conda.
2. Then tried `dag.replace_material_assignment("Vacuum", void_abs)` — that method DOES NOT EXIST in the `$SRC`
   build (only conda has it). AttributeError, instant fail.
CORRECT FIX (works): after building the model, `model.init_lib(output=False); model.sync_dagmc_universes()`
(populates `dag.cells`, 30 cells), then `dag.add_material_override(2, void_abs)` — in this build that method
sets `cells[2].fill` directly, and `create_xml_subelement` writes every synced cell's fill to model.xml (no
clobbering loop), so it survives `model.run()`. Then `finalize_lib()`, then assert `cells[2].fill` is
`void_absorber`. KEY GEOMETRY FACT (sync-debug `ssr/debug_sync.slurm`): cell 2 (plasma cavity) is a TRUE VOID
`fill=None`, NOT a 'Vacuum'-named material — target it by ID, not by material name. Cell 31 = outer graveyard
void (leave it, so leakage still works). Cells 3–30 = the radial build. Validated by an ON-vs-OFF smoke pair
(`vk4_on`/`vk4_off`, same source+seed): assert passed + coil fluxes DIFFER = physics genuinely changed.

Current real chain on `short` (non-preemptible): `9387596` real_vk (24M, VOID_KILL=1) → `9387597` verdict_vk
(writes `GATE2_REAL_VK_VERDICT.txt`, compares vs still-present pre-fix `coil_real_gate2.npz` CV=48%).
Superseded/cancelled VK attempts: `9379468/69/70` (no-op dict), `9387265`/`9387490` (API/void-name fails).
If CV collapses to single digits → decoupling validated → build the breeder optimizer. If it drops but not
enough, residual is fast-neutron scatter-back in the absorber (bump B-10 density) or reconsider freeze surface.

Earlier jobs (from the previous direction, still worth having):
- `9354539` mcgreedy — concentrate-multicoil dose (~9.5 h of ~17 h). Tests whether concentrating
  shield on the collective inboard beats single-coil concentration (peak 8.28e-3).
- `9354536` fwdfast — fast-group forward flux (superseded by the pivot, but nearly done).

## 5. First thing next session — checklist

Connect (login node, via the live control socket; the transfer node has no SLURM):
```
ssh ginsburg            # user does Duo; agent reuses ~/.ssh/cm-* control socket
export PATH=/cm/shared/apps/slurm/current/bin:$PATH
squeue -u tjk2147
```

Then, in order:
1. **Did smoke pass?** Look for `surf_smoke_master.h5` + `coil_smoke_master.npz` + `coil_smoke_gate2.npz`
   in `~/pstl_test/corrected`. If smoke failed, read `ssr/slurm_smoke_*.log` — most likely a
   `surf_source_read`/`surf_source_write` API detail; fix in `ssr/ssw_master.py` or `ssr/ssr_dose.py`
   and resubmit `bash ssr/submit_ssr_pipeline.sh`.
2. **When real_master + real_gate2 finish, run Gate 2** (the make-or-break test):
   ```
   cd ~/pstl_test/corrected && python3 ssr/gate2_analyze.py coil_real_master.npz coil_real_gate2.npz
   ```
   PASS = per-coil ratio full/ssr is constant (CV < 5%) and the hottest coil matches. That validates
   the decoupling. If CV is high, it's the interior-surface backscatter double-count — refine by
   killing plasma-void re-entry before trusting anything downstream.
3. **If Gate 2 passes**, build the high-resolution breeder optimizer: re-bin the frozen source
   (`ssr/patch_rebin.py`) to a fine mesh, cost = max-in-patch, keep the TBR floor, each iteration a
   cheap SSR run. This is the payoff.
4. **The NWL-vs-leakage figure** (Peterson asked for it, good for the Aug 26 talk): needs a shield-exit
   surface source too (add a second SSR that writes `cellfrom:7`), then
   `python3 ssr/flux_maps.py surf_real_master.h5 <shieldexit>.h5 maps.png --wout <wout_qh.nc>`.
5. **Earlier jobs**: pull `coil_step1b_multicoil_greedy.npz` and compare its peak coil flux to the
   single-coil placed (8.28e-3) and spread multicoil (2.33e-2); grab `qh_fwd_meshflux_fast.npz` and
   rerun `shield_opt/patch_worth/multigroup_signed_worth.py --fwd-fast` for completeness.

## 6. The pipeline files (`spf_prototype/shield_opt/ssr/`)

- `ssw_master.py` — full run; freezes inward first-wall source (`cellfrom=2`) + full-run coil ref.
- `ssr_dose.py` — reads the frozen source, transports the variable build to the magnets.
- `patch_rebin.py` — re-bin the source to any (tor,pol) resolution; VMEC-surface **area
  normalization** + **max-in-patch** (both Peterson flags). 6 local tests green.
- `flux_maps.py` — three-panel NWL / shield-exit leakage / transmission-ratio maps.
- `gate2_analyze.py` — the pass/fail verdict.
- `slurm_ssw_master.sh`, `slurm_ssr.sh`, `submit_ssr_pipeline.sh` — the SLURM layer.

QH DAGMC cell map (`dagmc_qh_step1_uniform.h5m`): plasma void = 2, first_wall = 3, multiplier = 4,
breeder = 5, back_wall = 6, shield = 7, gap = 8, vac_vessel = 9, thermal_shield = 10, magnets = 11–30.

## 7. Open decisions to settle with Peterson

- **Freeze surface**: first wall (keeps the breeder trade in the loop) vs blanket-exit (cleaner but
  fixes the breeder). We chose first wall; confirm.
- **Trade definition**: strict fixed-envelope breeder↔shield vs shield-thickness-with-fixed-breeder.
- **Terminology**: publish "ENL" as transmitted wall loading (TWL) / magnet leakage flux / magnet
  fluence map — the field's standard terms.
- **Patch resolution / adaptive mesh**: once Gate 2 passes, decide fine-uniform vs physics-adaptive
  (concentrate patches at the inboard pinch).

## 8. Deadlines

- **Aug 19** — present a relevant paper (group intro). Candidate: the Gil et al. augmented-Lagrangian
  coil paper (our coil source) or a shielding/SSR-decoupling methods paper.
- **Aug 26** — 30-minute research talk. Spine: honest patch-worth result → why proxies fail
  (moderation) → SSR decoupling → NWL-vs-leakage maps → high-res breeder optimization.

## 9. Access / environment notes

- SSH: `ssh ginsburg` reuses `~/.ssh/cm-tjk2147@...:22` control socket (ControlPersist 24h). The
  agent multiplexes over it; user handles Duo. **Transfer nodes (`motion`) have no SLURM** — need a
  login node (`ruth` etc.) to submit. `module`/`squeue` need `PATH=/cm/shared/apps/slurm/current/bin`.
- Work dir: `/burg-archive/home/tjk2147/pstl_test/corrected`; SSR code in `ssr/`.
- Envs: transport `conda activate spf-stellarator` (+ `PYTHONPATH=$SRC`, `LD_PRELOAD=$SRC/build/lib/libopenmc.so`,
  `SRC=/ginsburg/astro/users/tjk2147/spf_work/openmc_src`); geometry build `.../envs/pstl`; simsopt in `~/vmec_venv`.
- Local repo mirrors the analysis code; `spf_prototype/talk/` has the deck, prep notes, and the real
  QH 3D patch/coil figures (`gen_patch_fig.py`, `view_patch_geometry.py`).
