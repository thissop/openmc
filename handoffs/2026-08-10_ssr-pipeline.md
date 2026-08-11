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

SSR pipeline (smoke chain gates the real overnight chain via `afterok`):
- `9358861` smoke_master → `9358862` smoke_gate2 → `9358863` real_master → `9358864` real_gate2
- Last seen: smoke_master PENDING (burst nodes down/reserved). It is tiny; runs in minutes once a
  node frees, then real runs auto-fire.

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
