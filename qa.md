# Anticipated Q&A — SPF conformal-stellarator smoke study

Prepared for discussion with E. Peterson (MIT PSFC). Framing throughout: this is
machinery validation on a public equilibrium (`precise_QA`, nfp = 2) at 10,000
histories per config — a smoke run, not a reactor claim.

## Q: Why is the free-streaming TBR essentially zero?

A: In the "free" stream we disable all material interactions and transport only
the source birth directions, so there is no moderation. The 14 MeV neutrons never
slow down to the energies where the ⁶Li(n,t) cross-section is large, and with no
scattering there is also no capture channel, so they simply free-stream out of the
blanket. The free-streaming TBR ≈ 0 is therefore a construction artifact that
isolates the pure geometric steering; all real breeding shows up only in the
"scatter" stream (TBR ≈ 0.68–0.69 per source neutron).

## Q: How is the polarized source angular distribution verified against the DT differential cross-section?

A: The angular sampler was validated in a separate tier before this transport run.
We histogram sampled directions in cosθ (θ measured from the local magnetic field B)
and run a chi-square test against the analytic dσ/dΩ for each mode — sin²θ for the
perpendicular ("A") mode and 1 + 3cos²θ for the parallel ("B/C") mode. We confirm
the isotropic limit is recovered for the unpolarized reference, and we cross-check
the inverse-CDF sampler against an independent rejection sampler. Only after those
pass did we wire the CompiledSource into the transport model.

## Q: Is the DAGMC geometry watertight, and how was that validated?

A: The `.h5m` geometry was hand-built — offset surfaces faceted into closed,
nested watertight shells — to bypass a mesher/tooling blocker in the normal CAD
pipeline. Watertightness is a hard requirement for DAGMC ray tracing: a leaky
surface manifests as lost particles during transport, and this pipeline now runs
end-to-end across all configs, which is the practical smoke-level evidence that the
shells close. A rigorous `make_watertight` audit and a formal lost-particle count
are part of hardening this beyond the smoke run.

## Q: How is the implicit complement / plasma void handled?

A: The blanket is a stack of nested closed shells, so every point is inside exactly
one defined cell — the innermost region is the 5 cm vacuum scrape-off (SOL) that
stands in for the plasma void, and the layers build outward (W · steel · Be · FLiBe
· shield · coil case). The DAGMC implicit complement is then only the exterior of
the outermost coil-case shell, which is treated as vacuum where escaping particles
terminate. Because the shells are nested rather than abutting arbitrary CAD solids,
the complement never has to fill an interior gap, which keeps the geometry simple
and avoids overlap ambiguity.

## Q: Is the TBR difference between modes statistically significant right now?

A: No. The scatter-stream TBRs are 0.680 ± 0.008 (unpolarized), 0.687 ± 0.007
(perpendicular), and 0.693 ± 0.015 (parallel) — the parallel-vs-unpolarized gap is
about 2% and the error bars overlap within roughly 1σ. Breeding is a bulk
volumetric effect and is nearly blind to birth direction, so we do not claim
polarization changes the TBR. The real directional signal lives in the coil fast
flux and coil heating, not in the TBR.

## Q: How will the coil fast flux be converged?

A: The coil tallies sit behind 116 cm of blanket + shield + coil case, so at 10k
histories the MC error is large (e.g. free-stream fluxes of 31.68 ± 0.48,
26.84 ± 0.19, 35.96 ± 0.54 /cm²/src). The plan is global variance reduction —
weight windows generated via FW-CADIS from an adjoint/importance solve tied to the
coil tally — plus higher raw statistics, so that the deep-penetration flux and the
derived η tighten up. That is the main thing standing between the current indicative
trend and a converged number.

## Q: What is the analytic cross-check plan for the wall load?

A: We intend to cross-check the free-streaming first-wall load against an independent
analytic neutron-wall-loading kernel evaluated on the same first-wall patches. The
approach is a stellarator NWL kernel (e.g. an anarrima-style line-of-sight / solid-
angle integral, extended by a field-period perturbation to handle the nfp = 2
non-axisymmetry) computed on the identical patch geometry the MC tallies use. Because
the free stream is pure geometric steering with no material physics, it is exactly the
quantity an analytic kernel can reproduce, which makes it the cleanest validation
target before trusting the scattered results.

## Q: Does the conformal wall self-intersect, and how did you check?

A: Both — the outer shells de-center and, at the most strongly-shaped angle, lightly
self-fold. Because the constant offset (up to 116 cm) is large compared with the
plasma's inboard clearance, the shells become nearly circular and de-center from the
bean-shaped plasma (crescent pushed outboard at φ=0, rounder and centered at φ=90°). I
checked with a segment-crossing test on each boundary's 64-point poloidal contour at
every exported φ, plus a polygon-area monotonicity check: enclosed areas increase
strictly outward at all φ (the shells stay nested in that sense), but at φ=0 the shield
and coil-case shells develop a small self-intersection right at the inboard midplane,
where the offset overshoots the machine axis (coil_outer reaches R ≈ −2.4 cm). It is
localized to a couple of poloidal points and invisible at device scale, but it is a
genuine artifact of the naive constant-normal-offset construction — a physics-grade run
needs a proper non-self-intersecting offset (signed-distance / level-set, or inboard
clipping), which is also part of why the DAGMC had to be hand-built.

## Q: What ⁶Li enrichment did you use, and why does it matter?

A: The bundle specifies the breeder only as "FLiBe" and does not pin an enrichment
level, so I will not quote a value we did not set. Natural lithium is about 7.5% ⁶Li;
enrichment is a design knob that raises the TBR and hardens the breeding profile by
favoring the ⁶Li(n,t)α channel. Because this run is a directional-steering smoke test
with the reaction rate held fixed per source neutron, the enrichment was not a variable
we scanned — fixing and scanning it is a natural follow-up once the geometry and
variance reduction are hardened.

## Q: How does this scale to a real HELIAS / Helios-class geometry?

A: The pipeline is deliberately equilibrium-agnostic: the polarized CompiledSource,
the conformal offset-shell blanket builder, and the tally stack all take the
equilibrium as input, so swapping `precise_QA` (nfp = 2) for a higher-field-period
HELIAS-class device is a matter of feeding a new boundary and re-facetting the shells.
The physics that would change is the stronger 3D shaping and the different B-field
geometry, which is precisely what would test whether the directional steering survives
a more realistic wall. This is why the current run is framed as machinery validation
on a public equilibrium rather than a claim about any specific reactor.

## Q: Why is η quoted on the parallel mode, and what is the sign convention?

A: The convention is θ from the local B field: the parallel (B/C, 1 + 3cos²θ) mode
RAISES the coil fast flux (+13.5% free → +8.5% scatter vs unpolarized) and the
perpendicular (A, sin²θ) mode LOWERS it (−15.3% free → −11.3% scatter). We quote
η ≡ (scatter signal)/(free signal) = 8.5/13.5 = 0.63 ± 0.22 on the parallel mode
because it is the larger, cleaner signal, but the story is symmetric — a magnet-
protection benefit in this geometry would actually come from the perpendicular mode,
which is the one that lowers the coil load. η = 0.63 is indicative given the large MC
error; the robust result is the coherent directional ordering appearing in both coil
flux and coil heating across both streams.
