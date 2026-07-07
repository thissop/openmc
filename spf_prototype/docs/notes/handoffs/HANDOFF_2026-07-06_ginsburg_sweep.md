# HANDOFF — Ginsburg SPF stellarator sweep (2026-07-06)

**Audience:** Thaddaeus Kiker (Columbia) + Ethan Peterson (MIT).
**Purpose:** resume the spin-polarized-fusion (SPF) neutron-steering-in-stellarators work
from a clean, reproducible state. This is a *null with a diagnosed cause and a concrete fix*,
not a dead end — read §3 and §6 carefully.

---

## 1. TL;DR

We stood up the entire SPF stellarator pipeline end-to-end on Columbia's Ginsburg cluster
(conda OpenMC 0.15.3 + DAGMC/MOAB, ENDF/B-VIII.0, compiled SPF `openmc::Source` `.so` built
against the conda OpenMC) and ran the first-ever multi-device QUASR transport sweep: QUASR
simsopt-serial coils → exact Biot-Savart `b̂` → hard-gated field audit → coherence `C` &
nematic order `S_φ` → conformal blanket → pymoab-free DAGMC `.h5m` → OpenMC → directional
efficiency `η`. A 6-config free-streaming sweep produced a **null**: `η_source` did **not**
track `C` or `S_φ` (fit R² ≈ 0.1). The smoking gun is device **886079**, which has high
coherence (`S_φ = 0.92`) yet `η_source ≈ 0.08` — a coherent field physically cannot lose its
steering, so that number proves the **observable is broken, not the hypothesis**. Diagnosis:
`η` was read from the free-streaming **coil fast flux**, a per-device-geometry-confounded ruler
(exactly the confound an adversarial pre-analysis flagged at project start); `S_φ` itself is a
field-only quantity, unchanged and uncompromised. Fix: build a purpose-built inboard/outboard
partial-current NWL wall observable, validate it against anarrima at the axisymmetric limit,
then recover the low-`C` QH range.

---

## 2. What we stood up on Ginsburg (concrete, reproducible setup)

**Compute space.** Warm astro-group space at `/ginsburg/astro/users/<user>/spf_work`
(~15 TB), not scratch — persists between jobs. Repo pushed there and built in place.

**Environment.** Conda env `spf-stellarator`:
- **OpenMC 0.15.3**, built **WITH DAGMC + MOAB** (needed for the `.h5m` CAD geometry path).
- Cross sections: **ENDF/B-VIII.0** HDF5 library, `OPENMC_CROSS_SECTIONS` pointed at it.
- Compiler toolchain: **GNU 14** (this is what the `.so` was linked against — keep the
  compiler consistent with the conda OpenMC ABI or the source will fail to load).

**The compiled SPF source.** The `openmc::Source` C++ sampler was built into a shared object
**against the conda OpenMC** (historically the single most fragile step in this project) — it
**built on the first try** under GNU 14. **36 architecture-independent tests pass on x86**
(the tier-1 direction-sampler suite: `(a,b,c)↔(w_perp,w_par)` identities, χ² on sampled
directions, C++/Python parity on fixed seeds, inverse-CDF vs rejection KS).

**Pipeline stages (all executed end to end for the first time on multiple devices):**
1. Load QUASR simsopt-serial coil JSON (`simsopt_serials`) — parsed **without importing
   simsopt** (numpy only).
2. Exact **Biot-Savart** evaluation of `b̂(x)` from the coil filaments.
3. **Hard-gate field audit** on every device: unit-norm `|b̂|`, divergence-free check, and
   `B·n̂` on the LCFS. Devices failing the audit are dropped, not warned-through.
4. Compute the field-only order parameters: coherence `C = |⟨b̂⟩_s|` and nematic order
   `S_φ = (3λ_φ − 1)/2`, with `b̂` in the **local cylindrical basis** (purely toroidal
   tokamak field ⇒ `C = 1`).
5. **Conformal radial-build blanket** geometry (offset first wall around the LCFS).
6. **pymoab-free DAGMC `.h5m`** via the in-house `dagmc_writer` (numpy + h5py only — no
   pymoab dependency).
7. **OpenMC transport** with the compiled SPF source → `η`.

**Login-node discipline.** All compute runs via `sbatch`; the login node is used only for
quick queries (parsing, plotting, status). Nothing heavy on the login node.

---

## 3. Results — the 6-config free-streaming sweep (the null + the smoking gun)

Free-streaming (anisotropic SPF source, near-void materials) `η_source`, read out from the
**coil fast-flux observable**. 5 QA devices + 1 QH:

| device  |  C    | S_φ   | η_source |
|---------|-------|-------|----------|
| 803097  | 0.999 | 0.996 | 1.00  (anchor) |
| 923602  | 0.998 | 0.995 | 1.01  |
| 932746  | 0.998 | 0.994 | 0.53  |
| 59509   | 0.997 | 0.991 | 0.77  |
| 886079  | 0.974 | 0.921 | 0.08  ← **SMOKING GUN** (high coherence, ~0 steering ⇒ observable broken) |
| 1960314 | 0.760 | 0.365 | 0.49  ← the only QH; the single low-C point |

**Fits:** `η_source(C)` R² = **0.09**, `η_source(S_φ)` R² = **0.10** — both **null**.

See `../../decks/figs/eta_vs_Sphi.png` (η_source vs S_φ, with 886079 highlighted as the
confound) and `../../decks/figs/eta_vs_C.png` (η_source vs C).

**Honest interpretation.** The null is **not** evidence against `S_φ`. Device 886079 has a
highly coherent field (`S_φ = 0.92`), and a coherent field cannot physically lose its
directional steering — so `η_source ≈ 0.08` is a measurement artifact, not physics. The
`η` values here were derived from the free-streaming **coil fast flux**, which is dominated by
per-device coil geometry (standoff, coil shape) rather than by the wall steering signal — a
**confounded ruler**. This exact confound was predicted by an adversarial pre-analysis at the
start of the project. `S_φ` is a field-only quantity, fully decoupled from the neutronics, and
is **unchanged and unaffected** by this problem.

A second attempt to re-derive `η` from a **saved spatially-resolved wall-current mesh**
(without re-running transport) **also washed out**: at `C = 0.999`, where Schwartz predicts
**+43%** inboard, the saved tally shows **+0.1%** — the mesh is too coarse and net current
cancels. Conclusion: the fix requires a **new, purpose-built observable**, not reprocessing of
existing tallies.

---

## 4. Why S_φ is still the correct predictor (parity / order-matching) + the λ_φ ≈ C² identity

**Order matching (the crux).** The SPF emission kernel is `w(θ_B) ∼ 1 + a₂ P₂(cos θ_B)` —
**quadrupolar (P₂, second Legendre) and EVEN** under `b̂ → −b̂` (no dipole/P₁ term). Using the
identity

```
P₂(n̂·b̂) = (3/2)(n̂·b̂)² − 1/2 = (3/2) n̂ᵀ (b̂ b̂ᵀ) n̂ − 1/2,
```

**any** wall observable is a linear functional of the **second moment** `⟨b̂ b̂ᵀ⟩` (the
direction tensor `T`) **only**, and is **independent of the first moment `C`**. A rank-2
(spin-2) kernel couples only to the rank-2 field moment — multipole matching. `C` is the
**first** moment (rank-1, l=1, dipole) and is **off-resonance**: it predicts `η` only by
correlation, never causally. The right predictor is the **nematic order**
`S_φ = (3λ_φ − 1)/2` with `λ_φ = ⟨b_φ²⟩_s = T_φφ` — the liquid-crystal / Woodcock orientation-
tensor scalar for a **headless (apolar) director**, which is exactly what an even-in-`b̂` kernel
sees.

**Executable parity check:** flip `b̂ → −b̂` on half the source points. `S_φ` / the tensor is
**invariant**; `C` **collapses**. That is the operational proof that `C` is off-resonance and
`S_φ` is on-resonance.

**Why C nonetheless *looks* predictive — the λ_φ ≈ C² identity.** On real QUASR coils
`λ_φ ≈ C²` (correlation **0.998** over 12 devices, RMS(`λ_φ − C²`) = **0.0075**), because these
fields are **cap-shaped / co-toroidal** (the director points the same toroidal sense
everywhere). So `C` works as a **proxy** only. It **decouples** from `S_φ` for reversal devices —
e.g. **1190023** (high-iota QH, ~4% toroidal-sense reversal) — and those are precisely the
devices the coil audit happens to filter out. See `../../decks/figs/lambda_vs_C2.png` (near-
identity line; the sole reversal device 1190023 sits off the line).

**The predictor ladder: `C → S_φ → G[M]`.** The rungs are two *different kinds* of object. `C` and
`S_φ` are field-only statistics — they screen the field's steering potential and know nothing about
where the walls are. **`G[M]`** folds the field tensor against the *actual* detector geometry:
`G[r_w] = ∫ s(x) [ (3/2) n̂ᵀ(b̂ b̂ᵀ) n̂ − ½ ] dV / |r_w − x|²` — the tensor `b̂ b̂ᵀ` contracted with the
line-of-sight kernel `n̂ n̂ᵀ / r²`. Loosely `G[M] ≈ (field coherence ~ S_φ) ⊗ geometry`. So `G[M]`
does **not** replace `S_φ`: `S_φ` stays the field-physics variable; `G[M]` is the *analytic η itself*
(the free-streaming NWL — exactly what `anarrima` computes), i.e. a validator / cheap surrogate, not
a competing predictor. This is why today's η was confounded (a geometry-laden η read against a
geometry-blind predictor). `G[M]` is not far-future — it *is* what phase B validates against
anarrima; on the ladder it's the immediate next rung ("next to test").

---

## 5. The QA device family + why QA can't show the collapse + the QH geometry-fold problem

The 3D coil wireframes of the QA devices ("the QA forms" from QUASR) are in
`../../decks/figs/wireframes_qa.png`.

**Why QA alone can't test S_φ.** Quasi-axisymmetry (QA) *is* a coherent (co-toroidal) field,
so a QA-only family clusters near `C ≈ 1` / `S_φ ≈ 1` and **cannot exhibit the low-order
collapse** — there's no low-`C` lever arm inside QA. The only low-`C` point we have is the
single QH device (**1960314**, `C = 0.760`, `S_φ = 0.365`). To actually stress-test
`η_source(S_φ)` we need **more QH devices** spanning the low-`S_φ` range.

**The QH geometry-fold problem.** The current **conformal offset** first wall **self-
intersects** for aggressive high-nfp QH devices — the offset surface folds through itself.
Confirmed folds: **2630694** (nfp7), **2632243** (nfp8), **2588494** (nfp6). These are exactly
the low-`C` devices we most want, and they are currently unrunnable with the conformal build.
Fixing this (thinner build or ParaStell) is on the critical path to recovering the low-`C`
range.

---

## 6. NEXT STEPS (phase B — the fix; numbered and actionable)

1. **Build a proper first-wall NWL observable.** Replace the confounded coil-flux ruler with a
   clean **inboard/outboard partial-current** tally at the actual first wall (or a DAGMC first-
   wall-surface partial-current tally). It must be poloidally resolved (≥50 bins/wall) and use
   **partial** (directional) current so inboard/outboard net-cancellation cannot wash the
   signal out. This is the `n̂ᵀ (b̂ b̂ᵀ) n̂` wall contraction — the true steering quantity and the
   analytic observable. Design the reduction so the answer is unambiguous (the previous saved-
   mesh attempt failed because it was too coarse and net-current-cancelling).
2. **Validate per-mode against anarrima** (Schwartz's own analytic NWL package). Use
   `examples/square_torus.py`, which exposes `irrad_{inboard,outboard,floor}_{A,B}`. At the
   `C ≈ 1` axisymmetric limit the new observable must reproduce **+43% inboard / −22% outboard**
   for mode A (and the B/C mirror). This is the Paper-1 (methods) validation gate — do not
   proceed until it passes.
3. **Recompute `η_source` from the NWL observable** (S_φ is unchanged — no field recompute
   needed) and **re-test** `η_source(S_φ)` and `η_source(C)` on the 6-config set. Expect 886079
   to move to ~1.0 once the confound is removed; that is the falsifiable prediction.
4. **Recover the low-C range.** Re-run the folded QH devices with a **much thinner conformal
   build** (`configs/qh_queue.json`) or switch to **ParaStell** geometry to avoid the self-
   intersection. This is where the `S_φ` collapse (if real) actually appears — without it the
   family stays clustered at high `C`/`S_φ` and neither predictor can be distinguished.
5. **Then the scattering arm.** With a validated field-geometry `η_source(S_φ)`, add the
   **blanket scattering washout** `A(τ_scatter)` and test the factorization hypothesis
   `η(C, blanket) ≈ η_source(C) · A(τ)` (field geometry in `η_source`, scattering in a scalar
   `A(τ)`).

---

*Status: pipeline live on Ginsburg; first QUASR sweep complete; observable confound diagnosed;
`S_φ` intact. Resume at step 1 above.*
