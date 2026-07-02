# Speaker notes — Does spin-polarized-fuel neutron steering survive a real stellarator wall?

## Slide 1 — Title

- Quick line: I'm going to walk through whether spin-polarized-fuel neutron steering still shows up once you put it through a real, strongly-shaped stellarator wall — and I want to be clear up front this is machinery validation on a public equilibrium, not a physics claim about any reactor.

## Slide 2 — The question

- The starting idea is simple: spin-polarized fuel biases the direction DT neutrons are born, relative to the local magnetic field.
- If you can bias emission, maybe you can steer the fast-neutron load relative to where the magnets sit.
- Two things could realistically kill that steering: a strongly-shaped 3D wall that scrambles the geometry, and scattering in a thick breeding blanket that washes out the birth direction.
- So the question I'm actually testing is whether the directional signal survives both of those.
- One equation to anchor it: emission goes like sin-squared theta for the perpendicular mode and one-plus-three-cos-squared theta for parallel, where theta is measured from B.
- And again — this is validating the pipeline on a public equilibrium, not making a reactor claim.

## Slide 3 — The device

- On the left is the plasma boundary — this is precise_QA, a quasi-axisymmetric stellarator with two field periods.
- On the right you can see the conformal multi-layer breeding blanket wrapping that plasma shape.
- Conformal means the blanket layers are constant-offset surfaces that follow the plasma, not a simple torus.
- One line on scale: the exported geometry is at the base equilibrium scale, so the numbers on the axes are small — don't read anything into that, the shape is what matters.

## Slide 4 — Blanket cross-section

- This is the point of the whole geometry story: the blanket is non-axisymmetric.
- At the phi-equals-zero cut the plasma is a crescent pushed toward the outboard side, sitting inside a nearly circular blanket, so you get thick inboard FLiBe.
- At phi-equals-ninety it's rounder and more centered.
- The reason is that the offset is large compared to the plasma minor radius, so the outer shells relax toward circles and de-center relative to the bean-shaped plasma.
- One honest caveat if it comes up: the shells stay area-nested, but at phi-equals-zero the thick outer offsets slightly self-fold at the inboard midplane, where the offset overshoots the machine axis — a known artifact of a naive constant-offset wall that a physics-grade run would fix with a proper offset.

## Slide 5 — Radial stack

- Here's the radial layer stack, drawn to scale — 116 cm total from the plasma surface out to the coil case.
- Going outward: a 5 cm vacuum scrape-off, a thin tungsten first wall, steel, a beryllium multiplier, then the 50 cm FLiBe breeder, then shield, then the coil case.
- The FLiBe is the workhorse — it's the breeder and coolant and it's where most of the action happens.
- Keep this stack in mind, because the coil tallies I'll show are sitting behind all 116 cm of this.

## Slide 6 — The pipeline

- This slide is really about the fact that the whole pipeline now runs end to end.
- It starts with a polarized compiled source — the angular sampler was verified against the DT differential cross-section in a separate tier before this run.
- One line on geometry: the watertight DAGMC geometry was hand-built to get around a mesher tooling blocker.
- From there it's OpenMC transport, run as two streams — free-streaming, which is source directions only with no interactions, and scatter, which is the full blanket physics.
- And then the tallies: TBR, coil fast flux, per-layer heating, and a first-wall current mesh.
- The free stream isolates the pure geometric steering; the scatter stream tells me what's left after real physics.

## Slide 7 — TBR

- TBR here is tritium produced per source neutron in the FLiBe.
- With scattering on, it's about 0.68 to 0.69 across all three modes — unpolarized 0.680, perpendicular 0.687, parallel 0.693.
- Polarization barely moves it — the spread is within about one sigma — so I am not claiming polarization changes breeding. Breeding is a bulk volumetric effect and it's nearly blind to birth direction.
- The free-streaming TBR is essentially zero, and that's by construction: with no scattering the 14 MeV neutrons never slow down to where lithium-6 capture is large, and there's no capture channel, so they just free-stream out.

## Slide 8 — Coil fast flux

- This is where the directional signal shows up: fast flux above 0.1 MeV at the coil, which is a magnet-lifetime proxy.
- The ordering is clear and consistent — parallel raises the coil load, perpendicular lowers it, relative to unpolarized.
- Being honest about the sign: parallel is the one that puts more fast flux on the magnets; if you wanted a magnet-protection benefit in this geometry it would come from the perpendicular mode.
- The full parallel-to-perpendicular swing is about 34% in the free stream and about 22% with scattering.
- The key thing is the ordering survives scattering — the signal is degraded but it doesn't disappear.

## Slide 9 — Directional efficiency η (headline)

- This is the headline concept. Take the parallel signal: it's plus 13.5% over unpolarized in the free stream, and plus 8.5% with scattering.
- The ratio of those, 8.5 over 13.5, is about 0.63 — that's eta: the fraction of the ideal free-streaming steering signal that survives the real conformal wall plus scattering.
- So roughly two-thirds of the geometric steering makes it through — that's the story I want people to remember.
- But I have to be straight about the error: eta is 0.63 plus or minus 0.22. This is a 10,000-history smoke run and the coil tallies are buried behind 116 cm of blanket, so the MC error is large.
- The robust result is the trend, not the precise number — a coherent directional ordering appears in both streams. The exact eta is indicative and not yet converged.

## Slide 10 — Heating

- Per-layer heating with scattering on: the FLiBe absorbs most of the energy, on the order of 80% of what's shown — no surprise, it's the thick breeder.
- What matters here is the coil: parallel deposits more than unpolarized, which deposits more than perpendicular.
- That's plus 10.5% for parallel over unpolarized, and plus 23.1% over perpendicular.
- The important point is this is the same directional sign as the coil flux — a completely separate tally telling the same story.
- So the steering signal corroborates itself: it's visible in both coil flux and coil heating.

## Slide 11 — Caveats + what's next

- Let me be honest about the limits: this is a 10k-history smoke run, so the MC error is large, especially at depth behind the blanket.
- This is not a reactor claim — it's a public equilibrium and the goal was to validate the machinery.
- The wall itself is a naive constant-offset construction that self-folds a little on the inboard side at the most-shaped angle; a physics-grade run needs a proper non-self-intersecting offset.
- Eta needs variance reduction — weight windows or FW-CADIS — before I'd quote it to more than one significant figure.
- There's also an analytic wall-load cross-check still pending, against a stellarator NWL kernel on the same first-wall patches.
- After that: higher statistics, and then moving to physics-relevant equilibria rather than this demonstration case.

## Slide 12 — Backup: raw numbers table

- This is the backup with all the raw numbers and error bars — TBR, coil fast flux for both streams, coil heating, and eta — if anyone wants to dig into a specific value.
