# Check-in 2026-07-06 — Speaker notes

*Spoken script for the check-in deck. One block per slide; bullets are what to say out loud, not slide text. The slides carry minimal text — the reasoning lives here. Q&A for Ethan Peterson at the end. Presenter: Thaddaeus Kiker (Columbia); collaborator: Ethan Peterson (MIT).*

---

## Slide 1: Title

- One opener: "This is a working check-in on spin-polarized fusion as a neutron-steering tool — I'll show the OpenMC source we built, the full transport pipeline now live on Ginsburg, and the field-direction order parameter we use to predict whether the steering survives in a stellarator."

---

## Slide 2: The Spin-Polarized DT Source in OpenMC

- Start with the physics. When you align the deuteron and triton nuclear spins, the DT reaction stops emitting neutrons isotropically and instead emits *anisotropically about the local magnetic field direction* b̂(x). That directional bias is the entire lever this project pulls on.
- The emission kernel is w(θ_B) ∝ 1 + a₂·P₂(cos θ_B), where θ_B is the angle between the neutron direction and b̂, and P₂ is the second Legendre polynomial. Two structural facts to plant now, because the predictor argument later hangs on them: the kernel is **quadrupolar** (P₂, second order) and it is **even** under b̂ → −b̂. There's no dipole term — the emission cares about the *axis* b̂ lies on, not which way it points.
- The three modes: unpolarized is isotropic (a₂ = 0). Perpendicular polarization emits *across* B — peaks perpendicular to b̂ — and carries a roughly **50% reactivity enhancement**. Parallel polarization emits *along* B — peaks parallel to b̂. Same kernel, different sign of a₂.
- Implementation: it's a compiled C++ `openmc::Source`. Neutrons are born on the flux surfaces, the source reads b̂ from a field map, samples the birth direction in the *field-aligned local frame*, and rotates that direction into the lab frame.
- It's verified at the unit level: C++/Python parity on identical seed streams, χ² tests on the sampled directions, and recovery of the analytic free-streaming neutron wall load. The sampler itself is solid — it's the trustworthy actuator that everything downstream leans on.
- Attribution: the polarized-emission physics is Kulsrud 1982 for the reactivity enhancement and Schwartz 2025 for the wall-loading treatment we validate against.

---

## Slide 3: Building the Stack on Ginsburg

- Lead with the environment: the full transport stack is standing up on Columbia's Ginsburg cluster — a conda environment with OpenMC 0.15.3 built *with* DAGMC, plus MOAB and ENDF/B-VIII.0 cross sections, living in the group's warm storage.
- The headline for anyone who's built custom OpenMC sources before: the compiled SPF source `.so` linked against the conda OpenMC **on the first try**. That's historically the single most fragile step in this kind of work, and it's behind us.
- All architecture-independent tests pass on the cluster — the sampler behaves identically on Ginsburg's hardware as it does locally, so we can trust results generated in batch.
- Net message: the compute home is set, the environment is reproducible, and the risky build step is done. From here it's science, not toolchain.

---

## Slide 4: The Device-to-Transport Pipeline

- Walk the pipeline as a chain, because every link is a hard gate. It starts from a QUASR simsopt-serial coil description, which we parse **without importing simsopt** — we read the serialization directly, so simsopt is not on the critical path.
- From the coils we compute the *exact* Biot–Savart field — not an interpolant — to get b̂ everywhere in the plasma volume.
- Then the field passes a hard-gate audit: b̂ must be unit, the field must be divergence-free, and B·n must be tangent (≈ 0) on the plasma boundary so the flux surface is genuinely a flux surface. Devices that fail don't proceed.
- The audited field feeds the field-coherence statistics — the C and S_φ numbers this talk is about — and, in parallel, a conformal blanket geometry built around the plasma.
- That geometry is written to a **pymoab-free** DAGMC `.h5m` file using only numpy and h5py, which is a second brittle dependency removed from the critical path. That mesh plus the SPF source drives OpenMC transport, which produces the directional wall load.
- Operational note worth saying out loud: everything compute-heavy runs as a batch job. The login node is only for quick queries — we're being good citizens on a shared cluster.

---

## Slide 5: Scientific Goals

- Frame the payoff first: spin-polarized fusion lets us *steer where neutrons deposit their load*. That's a real engineering knob — you can aim flux away from the magnets, or toward the breeding blanket where you want it.
- Why the geometry matters: in a tokamak the field is nearly uniform in the local frame, so every source point steers the same way and the effect is coherent and strong. In a stellarator the field twists through the volume, so the steering from different points can partly cancel.
- So the two questions that define the project. **Q1, the physics question:** does the steering benefit survive the fully 3D stellarator field? **Q2, the screening question:** can a cheap field-only number predict *how much* survives — so we can rank a whole device family without running transport on every one?
- Keep those two threads distinct through the talk: one is "does the effect exist off-axisymmetry," the other is "can we predict it from the field alone." The order parameter we build is the answer to Q2.

---

## Slide 6: Directional Statistics: What We Measure

- The central object: the quantity that governs the steering is the *dispersion of the field direction* across the plasma — not |B|, and not quasisymmetry error. Those metrics constrain the field *magnitude*; the SPF kernel couples to the *direction* b̂.
- So we study the **source-weighted distribution of the unit field direction** b̂(x), weighted by the birth-rate density over the plasma volume. It's a distribution on the unit sphere, and its shape is what predicts the steering.
- We compute it in the **local cylindrical basis** — this matters. In that frame a purely toroidal (tokamak-like) field is a single point on the sphere: perfectly coherent. Compute it in lab Cartesian instead and even a perfect tokamak's toroidal winding would look dispersed, which is a frame artifact, not physics.
- The intuition to carry forward: a *tight cluster* of directions means a coherent field and strong steering; a *spread-out* distribution means the emission lobes point every which way and the benefit cancels.
- The next few slides show what that field shaping actually looks like on real QUASR devices — that's the geometry that sets b̂ and therefore the dispersion.

---

## Slide 7: Two Quasisymmetric Forms from QUASR

- These are two quasi-axisymmetric configurations pulled from the QUASR database — coils and plasma boundary rendered together in one view.
- Point at the coil set: note the complexity of the winding, the non-planar shaping. Those coils are what produce the field, and the shaped plasma boundary is what the field wraps around.
- The takeaway to say: this is the geometry that determines b̂(x) throughout the volume. Everything about the steering starts from how these coils shape the field, so it's worth grounding the audience in what the actual hardware and plasma shape look like before we go to cross-sections.
- These are quasi-axisymmetric, so they're relatively *gentle* — close to a tokamak in field structure. That's the near-coherent end of the family; the strongly-shaped devices come later.

---

## Slide 8: A Quasi-Axisymmetric Form: 803097

- On the left is the device itself — the coil filaments together with the plasma boundary. On the right is the poloidal cross-section of the plasma taken at three toroidal angles, drawn solid, dashed, and dotted.
- The key point to make on all three of these device slides: as you move around toroidally, the cross-section *rotates and breathes* — the shape changes from one toroidal cut to the next. That toroidal variation of the plasma shape is exactly what disperses the field direction across the volume.
- And that dispersion is precisely what our statistic captures. The three overlaid cross-sections are a direct visual proxy for how much b̂ wanders — more variation between the cuts means more spread in the direction distribution.
- 803097 is quasi-axisymmetric, so the variation between the three cuts is modest — the shape is close to constant around the torus. That's the near-coherent regime, high C, where we expect the steering to survive well.

---

## Slide 9: A Quasi-Axisymmetric Form: 923602

- Same layout: coils plus boundary on the left, three toroidal poloidal cross-sections — solid, dashed, dotted — on the right.
- Look again at how the cross-section evolves toroidally. The rotating-and-breathing of the shape between the three cuts is the fingerprint of the field-direction dispersion we're quantifying — it's the geometric origin of the whole effect.
- Because this is another QA form, the toroidal variation is again gentle. The three curves stay close to one another, which is the visual signature of a near-axisymmetric, high-coherence field.
- Worth contrasting device to device: even within the QA family the amount of breathing differs, and that's exactly the kind of variation the statistic turns into a single number we can rank on.

---

## Slide 10: A Quasi-Axisymmetric Form: 59509

- Same story, third device: filaments and boundary on the left, three toroidal cross-sections on the right in solid, dashed, dotted.
- Emphasize the mechanism one more time so it lands: the shape changing as you sweep toroidally is *what non-axisymmetry is* in real space, and it's what pushes b̂ off a single direction. The statistic reads that off the field.
- 59509 is again quasi-axisymmetric and therefore gentle — modest toroidal variation, so we'd expect it to sit on the coherent end of the family alongside the other two.
- Setup for the next slides: we've now seen the *geometry* that produces field-direction dispersion. The question is what *scalar* summarizes that dispersion in a way that predicts the wall steering — that's where the order parameter comes in.

---

## Slide 11: Coherence: The Mean Resultant Length

- Introduce the natural first statistic from directional statistics: the **mean resultant length**, C = |⟨b̂⟩| = |∫ s(x) b̂(x) dV| / ∫ s(x) dV, computed in the local cylindrical basis. It's the length of the average unit-direction vector, and it lives in [0,1].
- Physical reading: C is "how aligned are the field directions, on average." All b̂ parallel means the unit vectors add up and C = 1 — that's the perfectly coherent tokamak-like limit. As the directions disperse, they start cancelling and C → 0.
- It's the natural first guess, it's cheap, it's field-only, and on the real devices it *tracks* the shaping — coherent devices score high. So it earns its place as a diagnostic and we compute it in the pipeline anyway.
- But — and this is the setup for the next slide — C is the **first moment** of the direction distribution. It carries the mean, signed direction. It knows the difference between b̂ and −b̂. And that's the wrong *order* of statistic for a kernel that doesn't. Hold that tension for the next slide.

---

## Slide 12: Order Matching: The Kernel Is Quadrupolar

- This is the theoretical heart — go slowly. The emission kernel is P₂: second order in the field, and *even* under b̂ → −b̂. Write the identity: P₂(n̂·b̂) = (3/2)(n̂·b̂)² − 1/2 = (3/2) n̂ᵀ (b̂ b̂ᵀ) n̂ − 1/2, where n̂ is the emission (wall) direction and b̂ b̂ᵀ is the outer product.
- Read what that means: the kernel's angular dependence is a *quadratic form* in b̂. When you integrate emission over all source points to get any wall observable, the only field quantity that survives is ⟨b̂ b̂ᵀ⟩ — the **second moment**, the direction tensor. Every wall observable is a linear functional of that tensor.
- The consequence: any wall observable is **independent of the first moment C**. C literally doesn't appear in the expression for what the wall sees. It can only *predict* the steering by being correlated with the second moment, never by causing it.
- Say it as a resonance / selection-rule argument, because it's clean: a rank-2 kernel couples only to the rank-2 moment of the field-direction distribution — multipole matching. C is rank-1, a dipole, so it's off-resonance with the kernel. The right predictor has to be built from the second moment.

---

## Slide 13: The Nematic Order Parameter

- Define the matched predictor: the **nematic order parameter** S_φ = (3λ_φ − 1)/2, where λ_φ = ⟨b_φ²⟩ is the second moment of the field direction along the toroidal axis — the toroidal-toroidal component of the direction tensor.
- Name what it is: this is exactly the liquid-crystal orientation-tensor order parameter for a *headless* (apolar) director — the standard quantity for describing alignment of objects that don't distinguish head from tail. That's precisely what an even-in-b̂ kernel sees, so S_φ is on-resonance with the physics by construction.
- Be candid about the empirical near-degeneracy, up front so it's not a surprise later: on the real coils, λ_φ ≈ C² almost exactly — correlation about 0.998. So C *looks* perfectly adequate on this family, but only because it's acting as a proxy for S_φ through that near-identity.
- The two separate only for fields that **reverse toroidal sense** somewhere in the source region. Under such a reversal the first moment C partly cancels across the sign flip while the second moment — and therefore S_φ — is invariant. That's the regime that distinguishes the causal predictor from the proxy.
- The right panel shows λ_φ against C² across the family — the tight line is the near-degeneracy, and the physics predicts it should break precisely for sign-reversing fields.

---

## Slide 14: Fractionalizing the Directional Efficiency

- State the factorization: the directional efficiency splits as η(C, blanket) ≈ η_source(C) · A(τ). The slide is deliberately just the equation — let it breathe.
- η_source is the **field-geometry factor**: how much of the ideal steering survives the field's direction dispersion, measured in vacuum / free-streaming.
- A(τ) is the **blanket factor**: how much of the surviving steering makes it through a real scattering blanket, as a function of the scattering optical depth τ.
- The point of writing it this way: the field problem and the blanket problem separate. Screen devices on the cheap field-only factor, then apply one blanket curve. The next slide says why we expect the split to be clean.

---

## Slide 15: Two Ways the Steering Erodes

- There are two physically independent mechanisms that erode the steering, and they live in the two factors from the previous slide.
- Mechanism (i), in vacuum: the field-direction dispersion makes the coherent sum of emission lobes *partially cancel*. It's deterministic and purely geometric — set entirely by how b̂ is arranged across the volume. This is the stellarator-specific erosion, and it lives in η_source.
- Mechanism (ii), in the blanket: scattering randomizes the neutron direction roughly in proportion to the optical depth τ, washing out the birth-direction bias. It's stochastic, it's a material property, and it lives in A(τ).
- Because the two act on *different variables at different stages* — geometry at birth versus scattering en route — we expect them to factorize cleanly. That's the hypothesis the phase-B measurement will test directly, and the reason the screening strategy is worth building.

---

## Slide 16: Status and Next Steps

- Status: the pipeline is validated end to end on Ginsburg — from QUASR coils through the field audit, the directional statistics, the DAGMC geometry, and OpenMC transport. And C and S_φ are computed on the real QUASR coils across the device family.
- The directional-efficiency measurement is **implemented**. The immediate next run compares η against S_φ using the inboard/outboard first-wall load as the observable — the geometry-aware quantity the P₂ kernel actually produces.
- That observable is validated per mode against Schwartz's analytic neutron wall load through anarrima, his open-source package: at the axisymmetric limit our tally must reproduce the known oracles — **+43% inboard / −22% outboard** for the perpendicular mode, mirror image for parallel. That's the methods-validation gate before we trust any stellarator number.
- Then we add the scattering factor A(τ) by thickening the blanket and test whether η_source · A(τ) reproduces full transport — the factorization hypothesis.
- Finally we extend to the strongly-shaped, low-coherence devices, where b̂ disperses the most and the S_φ prediction is the most stringent. That's where the predictor either proves itself or gets falsified, and it's the destination the whole pipeline was built to reach.

---

## Anticipated Questions (Ethan Peterson)

**Q1. How exactly is η measured, and how is the observable validated?**
η comes from the inboard/outboard first-wall load — a directional (partial-current) neutron wall-loading tally on the actual first-wall surface, which is exactly the n̂ᵀ(b̂ b̂ᵀ)n̂ contraction the kernel produces. We validate it per mode against anarrima, Schwartz's analytic NWL package, at the C ≈ 1 axisymmetric limit: the tally has to reproduce the paper oracles of +43% inboard and −22% outboard for the perpendicular mode, and the mirror for parallel, before we trust it on any 3D field. That analytic cross-check is the pass/fail gate for the observable itself, independent of the predictor question.

**Q2. How does S_φ relate to ε_eff or the standard quasisymmetry metrics — is it just a repackaging of QS error?**
No — they live on different axes. ε_eff and QS error constrain |B| as a function of Boozer coordinates; they're magnitude-symmetry statements built for neoclassical transport. S_φ is a real-space second moment of the field *direction* b̂, built to match the P₂ emission kernel. Two devices can share a QS error and differ in S_φ, and vice versa. They may correlate on a given family because coherent fields tend to score well on both, but that's a family accident, not a derivation.

**Q3. You're using the vacuum coil field. Is that adequate versus the equilibrium flux-function field, especially for aggressive devices?**
For the near-axisymmetric QA devices the vacuum coil field is an excellent approximation to b̂, and the structural argument — even kernel → second moment → S_φ — holds for *either* field. Where it matters is the strongly-shaped, high-β, or sign-reversing devices, where the finite-β equilibrium field can differ from the coil field precisely in the region that decouples C from S_φ. So the plan is vacuum field for the QA screening now, and equilibrium b̂ for the aggressive cases where we make quantitative decoupling claims. It's a known refinement, not a blocker.

**Q4. QUASR has no quasi-isodynamic configurations — where's the QI arm?**
Correct, QUASR is QA and QH only, so we can't test QI from this catalog. The S_φ argument is configuration-agnostic — it's about b̂-direction dispersion, not the symmetry class — so it should carry over to QI, and QI is where dispersion is largest and the predictor would be most stressed, which makes it the most interesting missing arm. Getting a QI equilibrium means pulling from a different catalogue or generating one, so it's a longer-term extension rather than a near-term item.

**Q5. How will you build the strongly-shaped high-nfp devices whose conformal blanket self-intersects?**
That's a geometry-construction issue, not a physics one — the conformal offset wall self-intersects at high curvature. Two routes: shrink the radial build so the offset surface stays non-self-intersecting, which handles the moderately shaped cases; or drop the conformal-offset construction and use ParaStell, which is designed for tight high-nfp stellarator shaping. We only need a handful of low-coherence points to populate the regime where S_φ makes its sharp prediction, so even a few successful aggressive builds are enough.

**Q6. Does free-streaming η_source actually capture the benefit that matters once there's a real scattering blanket?**
On its own it captures only mechanism (i), the field-geometry dephasing — that's exactly what η_source is scoped to. The scattering-relevant part is the other factor, A(τ), which we measure separately by thickening the blanket. The whole point of the factorization is that geometry can be screened transport-free while the blanket contributes a single attenuation curve — and we verify the product reproduces full transport before trusting it, which is an explicit phase-B deliverable.

**Q7. Why compute the directional statistics in the local cylindrical basis rather than lab Cartesian?**
Because the kernel couples to the field direction relative to the *local* geometry, and the cylindrical basis is what makes a purely toroidal field a single coherent point. In lab Cartesian, the toroidal winding alone would smear even a perfect tokamak's b̂ around the sphere and register as dispersion — a pure frame artifact. The cylindrical frame removes that trivial rotation so the statistic measures genuine, physical departure from a coherent field.

**Q8. What's the concrete confirm-or-falsify criterion for S_φ as the predictor?**
Confirm: across a range that includes low-coherence QH devices, η collapses onto a monotone curve in S_φ with high R², and the b̂ → −b̂ flip test leaves both S_φ and η invariant while C moves — that's the in-code test that separates the causal predictor from the proxy. Falsify: a low-coherence, high-dispersion device that retains full steering while S_φ is low would break the tensor picture; so would a sign-reversing device whose η tracks C rather than S_φ once you use the equilibrium field. The upcoming low-coherence QH runs are what make that test possible.

**Q9. Why trust that the emission kernel is exactly even in b̂ — is there a sub-leading dipole you're dropping?**
To the order we work in — Kulsrud 1982 and Schwartz 2025 — the polarized DT emission is purely quadrupolar in the field axis: 1 + a₂P₂, no P₁ term. It's even because the physics is the alignment of the deuteron/triton spin *axis* with B, and the emission is symmetric about that axis with no handedness to seed a dipole. If a real depolarization or spin-precession effect introduced a small P₁, C would re-enter at that order — but that would itself be a physics finding worth isolating, and for the leading directional benefit even-in-b̂ is exact.

**Q10. How coarse is the field-map grid, and could it bias S_φ?**
S_φ is a volume-averaged second moment weighted by the birth density — a smooth integral, so it converges faster than pointwise field values and is forgiving of grid resolution. It also has to clear the hard field audit (unit b̂, divergence-free, B·n ≈ 0 on the boundary) before we use it. We bound any residual sensitivity directly and cheaply by refining the field map and confirming S_φ is stationary, and that's a check worth showing alongside the first η-vs-S_φ results.
