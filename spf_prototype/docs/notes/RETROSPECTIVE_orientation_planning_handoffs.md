# Retrospective — orientation, planning, and handoffs

*A condensed, in-retrospect narrative of the project's own scaffolding: the §0
API/environment ground-truth (`NOTES_orientation.md`, now in `notes/orientation/`),
the compaction-surviving state dumps (`HANDOFF.md`, `CONTEXT_REPORT.md`, in
`notes/handoffs/`), the roadmaps (`NEXT_STEPS.md`, `PLAN_next_phase.md`, in
`notes/planning/`), and the still-live scope/data docs (`DATA_NEEDED.md`,
`DEFERRED.md`, filed in `notes/planning/` and `reference/`). Substance preserved;
blockquotes verbatim and attributed. Positioning material draws on `PRIOR_WORK_BAE2025.md`
(kept intact in `reference/`).*

---

## Orientation — ground truth for this checkout (before any code)

**At the time** the OpenMC source API had refactored around 0.14/0.15 and could not be
trusted from training data, *and* the environment turned out to be bare. `NOTES_orientation.md`
pinned both, line-cited, cross-checked by a 6-agent search.

**API anchor (verified, `include/openmc/...`):**
- `openmc::Source::sample(uint64_t* seed) const` returns a `SourceSite` **by value**
  (source.h:88). The `.so` must export `extern "C" unique_ptr<openmc::Source>
  openmc_create_source(std::string)`; `dlsym` needs exactly that symbol.
- `SourceSite` fields: `r` (Position), `u` (Direction, unit), `E` (double eV, **no
  default — must set**), `time`, `wgt {1.0}`, `delayed_group {0}`, `surf_id`,
  `particle`. `Position`/`Direction` are the same struct with `.x/.y/.z`, `.dot`,
  `.cross`, `.norm` — everything Gram-Schmidt needs is built in.
- **The API differs from the spec** in one place that mattered:

  > "ParticleType::neutron() is a function call, not a plain enum value ::neutron" —
  > NOTES_orientation.md

- `double prn(uint64_t* seed)` returns [0,1) (PCG-RXS-M-XS, stream `STREAM_SOURCE=1`);
  `sample()` is called concurrently per-particle ⇒ the hard rule:

  > "no mutable members in our source; everything set in ctor, all methods const, only
  > prn(seed) for randomness" — NOTES_orientation.md

**Environment finding (the big one):**

> "This checkout is not built and has no Python environment" — NOTES_orientation.md

No `build/`, no `libopenmc.so`, no `openmc` on PATH, `import openmc` fails, system
python has no pip, no cross-section data. **We achieved** a fully-verified bring-up
recipe (§0.9): system HDF5, a dedicated `$HOME/spf_venv` with editable OpenMC + deps,
a `make install` into the venv prefix (a build-tree config is *not* consumable), NNDC
cross sections at `$HOME/nndc_hdf5/cross_sections.xml`, and a `find_package(OpenMC) +
target_link_libraries(... OpenMC::libopenmc)` custom-source build — proven end-to-end
by an in-tree ring-source smoke test (Leakage=1.0, surface current 1.0±0.0/src). Also
re-derived the physics oracle from scratch (corner cross sections 2/3, 1, 2/3, 1/3;
**C is not isotropic**; the sin²θ inverse-CDF cubic is the u→1−u mirror of the spec's
form — equivalent, endpoint-verified form used).

---

## Planning — verification-first roadmaps

**At the time** only the free-streaming tiers existed. Two planning docs charted the
path past them, both organised around the same rule.

`NEXT_STEPS.md` ranked three directions by readiness × value:
1. **angled / spatially-varying B̂** — cheapest; *"The C++ sampler already emits
   correctly for any B̂"* (only a field model needs adding); keeps the analytic anchor;
   the integration bridge to Peterson's physics-informed source.
2. **scattering ON with a real wall** — cheap machinery, but *"The hard part is losing
   the analytic ground truth: once scattering is on there is no closed-form reference"*.
3. **breeder blanket + TBR** — *"#3 is the headline reactor metric and a genuine new
   tier of work."*

`PLAN_next_phase.md` turned that into Tiers 4–6 with an explicit organising rule:

> "Every new physics step must ship with a replacement verification — a limiting case
> that recovers something already validated, a conservation law, a code-to-code or
> benchmark comparison, or a physical-sanity signature — not 'the plot looks
> reasonable.'" — PLAN_next_phase.md

It separated *source realism* (energy spectrum, pitched field — **keeps** an analytic
anchor via anarrima's angled kernels) from *transport realism* (scattering → heating →
TBR — **loses** it, needs new checks), confirmed data/tally feasibility (all nuclides
present; `H3-production`, `heating`, `damage-energy` valid), and recommended order
**4b → 4a → 5 → 6**. The delivered work (Tiers 5/6/7/8) followed this framing.

---

## Handoffs — the compaction-surviving state dumps

Two documents were written to survive context compaction and to hand the project to a
collaborator.

`HANDOFF.md` (Bae-comparison session) is the "read this first to continue" dump:
environment recipe, commit history, file map, don't-re-derive physics/API facts, all
per-tier key numbers, and the novelty recalibration. Its headline verdict:

> "our work = independent validation + reusable verified tool, not discovery" — HANDOFF.md

It records the hard workflow constraints:

> "commit to branch spf-prototype only, NEVER push, NEVER to main/develop, do not commit
> CLAUDE.md" — HANDOFF.md

`CONTEXT_REPORT.md` is a read-only, `file:line`-cited snapshot for a collaborator
scoping the 3-D stellarator phase — repo map (~3,170 LOC), sampler internals, exactly
what is and isn't validated, and a candid gaps/risks list. Its scope statement:

> "This is a research prototype, honest about scope … it is axisymmetric (tokamak-limit)
> only. There is no stellarator / 3D / general-field code yet" — CONTEXT_REPORT.md

It also caught the subtle honesty flags later folded into LIMITATIONS: the Tier-2b
residual σ≈0.68 framed loosely as "stdev≈1", and the B≡C "exact" agreement being a
bit-identical-RNG artifact rather than statistics. And the recurring clarification:

> "'isotropic' in the figures/results = unpolarized fuel baseline, not a modeling error"
> — CONTEXT_REPORT.md

Its top gaps to a 3-D stellarator MC study — real B̂(x) (only toroidal/constant exist,
but *the header rotation already handles any B̂*), true 3-D geometry (all CSG is
axisymmetric; `analytic_nwl.py` won't transfer), no 3-D analytic ground truth, and no
coil model — are exactly what the later Tier-8 field-map work began to address.

---

## Positioning — the Bae 2025 recalibration

**At the time** we found the closest prior work and recalibrated novelty downward.
`PRIOR_WORK_BAE2025.md` (kept in `reference/`) maps our modes to Bae's airtight
(perpendicular = our A ∝ sin²θ, +50% rate; parallel = our B/C ∝ 1+3cos²θ; unpolarized
= iso), warns of a **symbol collision** (Bae's (a,b,c) are Hupin–Navrátil vector/tensor
polarizations, *not* Schwartz's collision-mode fractions that we use), and states the
contribution honestly: *"No new physics."* — we independently reproduce Schwartz
(analytic + MC) and Bae's TBR-polarization sign in a different code, breeder (FLiBe vs
Pb-Li), and geometry (square torus vs spherical tokamak); the smaller magnitude is
explained by Bae's inboard-starved ST blanket. The one possible thin novelty — the
explicit free-streaming-vs-transport "~3× dilution" quantification — is *complementary*
to Bae, not a correction; do not oversell.

---

## Still-live scope and data docs

Two docs remain forward-looking rather than purely historical.

`DATA_NEEDED.md` is the ranked map of real **Helios** inputs that would replace the
public placeholders in Tier-8, each tagged to its `INJECT(helios)` code seam
(equilibrium/field #1, radial build, breeder spec, coil spec, source rate, profiles,
config scalars). The collaborator had offered boundary Fourier fits + plasma params.

> "Nothing here is a Helios physics claim until these are filled." — DATA_NEEDED.md

`DEFERRED.md` is the honest list of what the 3-D-field tier does and does *not* support
— the decisive caveat being that the **3-D physics is only in the source**:

> "Transport geometry is axisymmetric (a layered box-torus radial-build rig), NOT a 3D
> stellarator first wall. The 3D physics is only in the SOURCE … No field-period / QA
> geometry effect is claimed." — DEFERRED.md

Both are filed with the historical notes for tidiness, but a future agent should treat
them as **live** inputs to the next phase, not closed history.
