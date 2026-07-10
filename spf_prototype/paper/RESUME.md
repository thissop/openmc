# RESUME — pick up the capstone paper from here

*State as of 2026-07-09. Read this first; it is the single entry point.*

## Where everything is

```
spf_prototype/paper/
  PAPER_DIRECTION.md      # thesis, structure, IN/OUT, novelty framing, fidelity ladder
  RESUME.md               # this file
  manuscript/
    spf_openmc.tex        # THE paper (14 pp, compiles clean). Built on your mc_spf_nwl.tex voice.
    spf_openmc.pdf        # compiled
  data/
    PROVENANCE.md         # every number -> its source run
    *.csv                 # frozen tables (steering, rate, bare-cylinder, validation, predictor_n15)
    conformal_maps/*.npz  # 15-device (theta,phi) NWL maps (unpol/perp/par + dA)
  src/
    _style.py             # smplotlib minimalist style
    fig_*.py              # one script per figure; read only data/, write figures/
  figures/*.pdf,*.png     # generated
```

Compile the paper: `cd manuscript && pdflatex spf_openmc && pdflatex spf_openmc`.
Regenerate all figures: `cd src && for f in fig_*.py; do python $f; done`.

## What is DONE

- Native SPF source + StellaratorSource, verified (Tiers 1-2, sampler self-test, a2-linearity).
- Free-streaming analytic<->MC cross-validation on real QUASR QA 59509; point-patch-limit proof.
- Transport results, all frozen in `data/`: **3x scattering dilution** (Tier 5), **rate-vs-steering /
  B-sweet-spot** (Tier 7), **3-D field end-to-end / eta=0.89** (Tier 8b), and **15-device conformal
  NWL maps** with the peaking-optimal-a2 scan.
- **n=15 predictor result (NEW, this session):** benefit is weak (0-16%, median 4%) and NOT
  forecastable from plasma shape (nematic rho=+0.40 n.s.); best predictor is the transport-free
  headroom proxy (PF_geo rho=+0.53, borderline). The strong n=5 correlation (rho=+0.90) regressed to
  noise -- a properly-powered null that reinforces the "must run transport" thesis. Frozen in
  `data/predictor_n15.csv`, figure `fig_predictor`.
- Manuscript §1-§11 written in your dense voice; abstract/intro/summary reframed; Bae (a,b,c)
  symbol-collision footnote, precise_QA quadrature-breakdown, and data-availability paragraph added.

## Punch list to finish (ranked)

1. **Co-author sign-off.** Ethan Peterson is on the byline (from your draft). Confirm he is
   co-authoring and has seen §7-§10 before submission. His call, not assumed.
2. **Firm two citations.** `schwartz2025` = arXiv:2507.11758 (set). `lion2022` ("A deterministic
   method for stellarator NWL") and `parisi2024` ("Spin-polarized fusion in magnetic confinement")
   are plausible but need exact journal/volume/DOI -- confirm before submission. Bae is
   Nucl. Fusion 65, 086051 (2025); Lyytinen is Nucl. Fusion 64, 076042 (2024) (both verified).
3. **Conformal + scattering bridge run (the one worthwhile new run).** §7-8 are box-torus, §9 is
   conformal; Limitations states this. A single scattering-on run on a conformal wall closes the
   fidelity gap and shows the 3x dilution survives on real geometry. **The script already exists:**
   `python/run_conformal.py` runs BOTH free-streaming (density_scale=1e-4) and scattering
   (density_scale=1.0) per mode on a conformal DAGMC wall with the W/steel/Be/FLiBe/shield/coil
   blanket, tallying phi-resolved first-wall current + per-material heating. **Not yet run on a QUASR
   device** -- its `main()` is wired for the `precise_QA` demo. To run on e.g. 886079: (a) build the
   multi-layer mesh, `python build_dagmc.py quasr886079_wall quasr886079_blanket.h5m` (uses the
   existing `quasr886079_wall_surface.npz`); (b) call `build_model(abc, h5m, fluxmap_stem, R0_cm,
   a_cm, density_scale=dscale)` with the device's `quasr886079_cm_fluxmap` stem and its R0/a, looping
   dscale in (1e-4, 1.0) and abc over the three modes; (c) tally inboard/outboard current per
   stream/mode and save to `data/quasr886079_conformal_scatter.npz`. Expect the free stream to match
   the frozen `conformal_maps/quasr886079_conformalmap.npz` steering (an anchor) and the scatter
   stream to show the ~3x dilution on real geometry.
4. **Target venue decision.** Nucl. Fusion (fits Bae/Lyytinen), Comp. Phys. Comm. (tool angle), or
   Fusion Eng. & Design. Sets length/format; the draft is article-class and venue-agnostic for now.
5. **Prose density pass on §7-§10** against your §1-§6 bar -- you are the calibration.

## Ginsburg (how to get more data)

Connection: you authenticate once (Duo), then `ssh -o BatchMode=yes ginsburg` reuses the socket.
See `docs/GINSBURG_VIA_CLAUDE.md`. Working tree: `/ginsburg/astro/users/tjk2147/spf_work/spf_pp`.
Env: conda `spf-stellarator`, `LD_PRELOAD=$SPF/openmc_src/build/lib/libopenmc.so` (built with
`-DOPENMC_USE_DAGMC=ON`), `OPENMC_CROSS_SECTIONS=$SPF/xs/endfb-viii.0-hdf5/cross_sections.xml`.

Batches that produced the current data:
- `vmec_batch30.sbatch` (job 8885457): 30 QUASR devices -> 22 converged VMEC fluxmaps. QH devices
  self-intersect in fixed-boundary VMEC and fail (expected).
- `conformal_batch30.sbatch` (job 8885458, afterany): conformal free-streaming NWL maps -> 15 devices.

To grow the device set: pick more IDs from `sweep/sweep_out/candidate_pool.csv`, edit the `IDS=(...)`
arrays in those two sbatch files, resubmit chained (`sbatch --parsable vmec... ; sbatch --dependency
=afterany:<id> conformal...`), then pull with `tar czf` on the remote + scp (streaming tar-over-ssh
corrupts; build the tarball remotely and scp the whole file). Re-run `python/predictor_test.py` after
updating its `DEV` list.

## Bridge run (conformal + scattering) recipe

`ginsburg_jobs/conformal_scatter.sbatch` builds a layered conformal wall (LCFS + firstwall/blanket
layers via `stellarator_geometry.build_layers`) with materials and runs scattering, tallying wall
current + heating. If it lands cleanly, results go to `data/quasr<ID>_conformal_scatter.npz`; scp
them into `paper/data/` and add a figure comparing free-streaming vs scattered conformal steering
(mirror of `fig_dilution` but on the real wall). If it fails, the box-torus Tier-5 dilution already
carries the result and the fidelity-ladder framing in Limitations covers the gap.
