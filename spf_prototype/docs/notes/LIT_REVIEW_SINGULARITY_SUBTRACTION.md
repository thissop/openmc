# Prior-art / novelty review: closed-form singularity subtraction for near-singular trig-polynomial line integrals

**Scope.** Numerical-analysis and computational-physics literature only. The fusion/plasma/NWL literature is covered separately and is deliberately excluded here. The object under review is the *numerical method*, not the application.

**The method (restated).** We evaluate near-singular 1-D line integrals

    I = ∫ f(φ) / D(φ)^{3/2} dφ,     D(φ) = |Δ(φ)|²,

arising from free-streaming / line-of-sight / inverse-square integrals over a *parametric source curve* seen from a wall point. The distinguishing structural fact is that **D(φ) is a trigonometric polynomial** (from a Fourier/harmonic description of the geometry). At grazing incidence min|Δ| → 0 and the integrand grows a peak of height ~D_min^{-3/2}·(width) with width ~√D_min; fixed-node Gauss–Legendre (GL) fails (>6000 nodes or diverges). We claim four novel ingredients:

1. Because D is a trig polynomial, the peak location φ* = argmin D, depth D_min, and curvature D''(φ*) are available in **closed form** (critical points/roots of a trig polynomial) — no root search.
2. Subtract an **analytically integrable** peak model m(φ) = (c0 + c1(φ−φ*) + c2(φ−φ*)²)/(D_min + a(φ−φ*)²)^{3/2} — with the **linear c1 term essential** near grazing — then integrate the smooth remainder f/D^{3/2} − m with few-node GL.
3. Keep the whole construction **differentiable / autodiff (JAX)**: c0,c1,c2,a,D_min are smooth closed-form functions of geometry; φ* is `stop_gradient`'d (exact, because the subtraction identity holds for any expansion center).
4. Result: 12 nodes beats 6144 GL; gradients match FD to 1e-8.

---

## A. Singularity subtraction, generally (the textbook baseline)

**Verdict: textbook. Claim nothing here.** "Subtraction of the singularity" is classical. The canonical origin is **Kantorovich & Krylov, *Approximate Methods of Higher Analysis* (1958)**: rearrange so the kernel singularity is (partly) cancelled by adding/subtracting an auxiliary function reproducing the leading singular behavior, then quadrature the smooth remainder. It is standard in **Kress, *Linear Integral Equations*** and Atkinson's numerical-integral-equations texts, and is a workhorse of the boundary-element/boundary-integral (BEM/BIE) community.

- Singularity subtraction as a unified BEM tool across singularity orders: overview and elementary analytic integrals — Järvenpää/Ylä-Oijala-style treatments; "Elementary analytical integrals required in subtraction of singularity method for evaluation of weakly singular boundary integrals," *Eng. Anal. Bound. Elem.* — https://www.sciencedirect.com/science/article/abs/pii/S095579970600097X
- Singularity subtraction in the numerical solution of integral equations, *ANZIAM J.* — https://www.cambridge.org/core/journals/anziam-journal/article/singularity-subtraction-in-the-numerical-solution-of-integral-equations/C42E161C2F215CFC2278C29E8BA88C9D
- Kress, *Linear Integral Equations* (Springer) — https://link.springer.com/book/10.1007/978-1-4614-9593-2

**True vs near singularity — the important distinction for us.** Classical singularity subtraction addresses **true** singularities: the field point lies *on* the source (integrable-but-singular kernels, e.g. weakly singular log or 1/r in BIE). Our problem is a **near-singular** integral: the field/wall point is *close to but not on* the source curve, so the integrand is smooth but has a sharp, tall, narrow peak whose scale is set by the standoff D_min. The peak is bounded (not a true singularity), so the classical "subtract the exact singular kernel" recipe does not directly apply — one must subtract a *finite-standoff regularized model* (the D_min in the denominator). That is a real and well-recognized distinction, and it moves us out of the plain-textbook regime into the near-singular quadrature literature (§B), which is where the genuine prior art lives.

---

## B. Near-singular quadrature (the field that actually contains our competitors)

**Verdict: large, active field. Our approach is a member of the "subtract/factor the near-singularity" family, closest to singularity-swap. It is not subsumed by transform-based methods, but it must be positioned carefully against singularity-swap.** The main techniques:

- **Quadrature by Expansion (QBX)** — Klöckner, Barnett, Greengard, O'Neil, *J. Comput. Phys.* 252 (2013) 332 — https://www.sciencedirect.com/science/article/abs/pii/S0021999113004579 ; arXiv:1207.4461 — https://arxiv.org/abs/1207.4461. Forms a local (Taylor/multipole) expansion of the layer potential about a center offset from the boundary and evaluates the near/on-surface value from the smooth far part. Dimension-independent, FMM-compatible. **Differs from us:** QBX expands the *potential* (the whole integral as a function of target position), not the integrand's peak; it never needs the distance function's critical point in closed form. It is more general but heavier, and it is designed for layer potentials of elliptic PDEs, not a bare inverse-cube geometric kernel over a 1-D curve.
- **Singularity swap / singularity-cancelling line-integral quadrature** — af Klinteberg & Barnett, "Accurate quadrature of nearly singular line integrals in two and three dimensions by singularity swapping," *BIT Numer. Math.* 61 (2021) — https://link.springer.com/article/10.1007/s10543-020-00820-5 ; and the follow-up "Singularity swap quadrature for nearly singular line integrals on closed curves in two dimensions," *BIT* (2024) — https://link.springer.com/article/10.1007/s10543-024-01013-0. **This is the closest prior art.** It targets exactly our object: nearly singular *line integrals over a source curve*. The near-singularity is a **root of the (complexified) distance function |Δ(φ)|²**; they locate the complex preimage φ* by **Newton iteration**, factor out the singular factor, expand the remainder in monomials, and integrate the singular×monomial products analytically. See §D for the head-to-head.
- **Sinh transformation** — Johnston & Elliott (2005) and many extensions; "The sinh transformation for evaluating nearly singular boundary element integrals over high-order geometry elements," *Eng. Anal. Bound. Elem.* — https://www.sciencedirect.com/science/article/abs/pii/S0955799712002238 ; distance-sinh — https://link.springer.com/article/10.1007/s00466-013-0913-0. A change of variable that analytically clusters nodes into the peak. **Differs from us:** it *reshapes* the quadrature abscissae rather than *subtracting* the peak; it needs the foot-of-perpendicular parameter and standoff but not the peak curvature in closed form, and it does not produce an analytic peak value. It is a strong, simple competitor for accuracy but is a transform, not a subtraction, and is less obviously differentiable-friendly.
- **Telles transformation** (cubic polynomial coordinate transform with vanishing Jacobian at the projection of the singular point) — and **Duffy transformation** (vertex singularity cancellation), generalized Duffy / augmented Duffy / distance-Duffy — https://www.sciencedirect.com/science/article/abs/pii/S0955799719305405 . Same "reshape, don't subtract" character as sinh; transform-family, not subtraction-family.
- **Continuation / preimage-Newton + singularity subtraction + transplanted Gauss** for curved elements — Montanelli et al., "Computing weakly singular and near-singular integrals over curved boundary elements," *SIAM J. Sci. Comput.* — https://epubs.siam.org/doi/10.1137/21M1462027 ; arXiv:2111.13151 — https://arxiv.org/abs/2111.13151 ; strongly-singular sequel arXiv:2309.17039. Uses Newton to find the singularity preimage, then subtracts. Again root-finding, not closed-form.
- **Adaptive / graded panels, generalized Gaussian quadrature, error-estimate-driven adaptivity** — e.g. axisymmetric-surface adaptive layer-potential quadrature arXiv:2412.19575. General but node-hungry; this is precisely what we beat (the "6144 GL" baseline).

**Summary for B:** the transform methods (sinh/Telles/Duffy) are a *different family* and do not subsume us. Within the *subtraction/factoring family*, singularity-swap (af Klinteberg–Barnett) is functionally the same idea specialized to line integrals, and is the work we must distinguish ourselves from.

---

## C. View-factor / radiative / rendering illumination integrals

**Verdict: the inverse-square grazing singularity is well known here, and analytic/edge-integral solutions exist — but for *polygonal/flat* geometry, not for a curved parametric source with a trig-polynomial distance. No direct precedent for our specific subtraction.**

- Singularity-subtraction is explicitly used in radiative transfer: "Application of the singularity-subtraction technique to isotropic scattering in a planar layer," *JQSRT* (1985) — https://www.sciencedirect.com/science/article/abs/pii/0022407385901414 . Confirms subtraction is standard practice at large optical thickness / grazing incidence, but for the scattering integral equation, not a geometric line-of-sight kernel.
- View-factor community knows grazing/near-contact singularities and mostly answers them with **closed-form analytic results for polygons** rather than peak subtraction: Narayanaswamy, "An analytic expression for radiation view factor between two arbitrarily oriented planar polygons," *Int. J. Heat Mass Transfer* (2015) — https://www.sciencedirect.com/science/article/abs/pii/S0017931015008492 ; contour-integral technique — https://www.sciencedirect.com/science/article/abs/pii/S0306454909003016 .
- Rendering: **analytic irradiance from polygonal area lights** via Arvo's irradiance tensors and edge (boundary) integrals — "Geometric Derivation of the Irradiance of Polygonal Lights," https://hal.science/hal-01458129/document ; projected solid angle via Gauss–Bonnet, arXiv:2205.12241 ; linearly-transformed cosines. These convert the 2-D inverse-square area integral into a **1-D edge integral solved in closed form for straight edges** — conceptually parallel to our reduction to a 1-D curve integral, but they exploit *flatness/straightness*, so no near-singular quadrature is needed. When the emitter is a *curve* rather than a polygon, that closed form is lost and our trig-polynomial-distance subtraction becomes relevant.

**Net:** the *phenomenon* (inverse-square grazing peak) and the *tool* (singularity subtraction) both appear in this literature, but the *specific move* — closed-form peak parameters from a trig-polynomial distance of a curved source, subtracted with a linear term — does not appear.

---

## D. Exploiting the closed-form structure of D(φ) (the narrowest, strongest claim)

**Verdict: this is the real novelty locus, but it is a *refinement*, not a category, over singularity-swap. Partially novel.**

The entire near-singular-line-integral family already *knows* the near-singularity is a **root/critical point of the distance function** D(φ) = |Δ(φ)|² continued to the complex plane. That is the founding observation of singularity-swap (af Klinteberg–Barnett) and of the preimage-Newton subtraction schemes (Montanelli et al.). The literature's standing practice is to find φ* by **Newton iteration** because the source curve is a general spline/panel and D has no special algebraic form.

Our added structural fact — **D is a trigonometric polynomial**, so φ* (argmin), D_min, and D''(φ*) come from roots/critical points of a trig polynomial in closed form (or via a companion-matrix eigenproblem), *no iterative search* — is:
- **Genuinely useful and, in this exact packaging, not something I found published.** Search for "trigonometric polynomial distance ... critical point ... quadrature" surfaces the algebraic machinery (critical points of squared-distance functions reduce to univariate polynomial roots; Keplerian-distance critical-point work, arXiv:2305.13900) and the fact that "roots of the displacement/distance function limit polynomial-approximation convergence" (arXiv:2111.13151), but **not** anyone using a trig-polynomial D to obtain closed-form singularity-subtraction parameters.
- **But narrow.** A referee will note that Newton on D converges in 2–3 iterations from an obvious initial guess, so "closed-form vs Newton" is a convenience, not a complexity change — *unless* the differentiability argument (§E) is invoked, which is where closed-form actually earns its keep (a `stop_gradient`'d exact φ* avoids differentiating through a Newton solve).

**Also note:** the **linear c1 term** in the peak model — argued essential because c0 ~ √D_min → 0 while c1 = O(1) at grazing — is a specific, correct, and not-obviously-published detail. Most subtraction models subtract an even model (c0 + c2 ξ²)/(D_min + a ξ²)^{3/2}; keeping the odd term to capture the asymmetry of f and of higher D-derivatives about φ* is a real refinement. This strengthens the D-claim.

---

## E. Differentiable near-singular quadrature (the second novelty hook)

**Verdict: differentiable BEM/near-singular quadrature now exists, but "differentiable + near-singular-line-integral + closed-form-peak subtraction with stop_gradient'd exact expansion center" as a package appears novel.**

- **JAX-BEM: Gradient-Based Acoustic Shape Optimisation via a Differentiable Boundary Element Method** (Hipperson, Hargreaves & Cox, Salford/Funktion-One), arXiv:2604.21431 — https://arxiv.org/abs/2604.21431 . **[Full text read directly 2026-07 — caveat resolved.]** Differentiable BEM in JAX for acoustic shape optimization. Confirmed details: (1) kernel is the **1/r Helmholtz Green's function** `e^{ik|r−r′|}/(4π|r−r′|)` — a *true on-surface* self-term plus adjacent-element near-singular terms, a weaker singularity than our D^{−3/2} inverse-cube streaming kernel. (2) It does **not** disclose or exploit any closed-form-peak / distance-function structure — §3.1 only states *"handling of singular integration routines… special handling is required for element self-interaction (singular) and adjacent element interaction (nearly-singular)"* and effectively defers to bempp-style reference handling. (3) **Crucially, its differentiability strategy is the *opposite* of ours:** it pushes the non-differentiable control flow *out* of the autodiff path (*"adjacency and selection of integration types is re-computed once per iteration before the differentiable sequence"*) and applies `JAX.custom_vjp` only to bypass the **GMRES linear solve** — **not** to the singular integration. So JAX-BEM makes the *solver* differentiable while keeping the singular quadrature *outside* the differentiated region. anarrima instead differentiates *through* the near-singular quadrature itself (subtraction-identity center-invariance + stop_gradient'd closed-form φ*). This makes JAX-BEM a **contrast citation, not prior art**, for hook E: the state-of-the-art differentiable-BEM paper explicitly sidesteps the very thing anarrima does.
- **Fully Differentiable Boundary Element Solver for Hydrodynamic Sensitivity** (arXiv:2501.06988), **JAX-FEM** (arXiv:2212.00964), **JAX-SSO** (arXiv:2407.20026), **quadax** (JAX adaptive tanh-sinh for singular/near-singular integrands, https://github.com/f0uriest/quadax) — establish that differentiable solvers and even differentiable singular-quadrature *tooling* are mainstream. quadax notably provides tanh-sinh (a transform method) for near-singular integrands under autodiff, i.e. a differentiable competitor exists at the tooling level.
- Classical **shape-derivative of boundary integral operators** theory (Costabel–Le Louër, arXiv:1105.2474; Boundary Integral Evaluation of Surface Derivatives, SIAM SISC) shows analytic differentiability of BIE operators, but by hand-derived shape calculus, not autodiff-through-the-quadrature.
- The `stop_gradient`-on-nodes idea is generic in the JAX world ("quadrature nodes/weights are non-trainable buffers ⇒ autodiff treats them as stop-gradient"). Our specific and *correct* twist is that stop_gradient on **φ*** is not an approximation: the add-subtract identity ∫f/D^{3/2} = ∫(f/D^{3/2} − m) + ∫m holds for *any* expansion center, so freezing φ* leaves the value (and its exact gradient) unchanged. That justification is sound and I did not find it stated in the near-singular-quadrature literature.

**Net for E:** differentiable near-singular quadrature is no longer novel per se; the *specific, provably-exact-gradient closed-form subtraction* for a curved source is the defensible new piece. The one paper that could have undercut this (JAX-BEM, arXiv:2604.21431) was read in full and does the *opposite* — it makes the linear solve differentiable while hoisting the singular integration *out* of the autodiff path — so it strengthens the contrast rather than pre-empting the hook.

---

## Four-way verdict

**(ii) A novel combination of known ingredients — bordering on (iii) a specialization of singularity-swap.**

- Singularity subtraction (A): textbook (Kantorovich–Krylov, Kress).
- Near-singular line-integral factoring at a root of the complex distance function (B, D): known — **singularity-swap, af Klinteberg & Barnett**.
- Autodiff-through-a-solver / stop_gradient'd quadrature nodes (E): known/mainstream (JAX-BEM, quadax).
- **The combination that is new:** (closed-form φ*/D_min/D'' from a *trigonometric-polynomial* distance) × (an *odd-term-inclusive* rational peak model with elementary antiderivative) × (provably-exact autodiff gradients via a stop_gradient'd, expansion-center-invariant subtraction identity), applied to the **D^{-3/2} inverse-cube 1-D kernel**. No single located work assembles these four.

It is **not** (i) genuinely-new-from-scratch — every ingredient has ancestry — and it is **not** cleanly (iv) subsumed, because the closest general method (singularity-swap) uses Newton root-finding and monomial expansion, not a closed-form differentiable subtraction model, and does not make the autodiff-exactness argument.

---

## Closest prior work, and precisely how we differ

**Closest single work:** L. af Klinteberg & A. H. Barnett, *"Accurate quadrature of nearly singular line integrals in two and three dimensions by singularity swapping,"* BIT Numer. Math. 61 (2021), https://link.springer.com/article/10.1007/s10543-020-00820-5 (with the 2024 closed-curve sequel, https://link.springer.com/article/10.1007/s10543-024-01013-0).

Same problem class (nearly singular integrals over a source *curve*), same founding observation (the near-singularity is a root of the complexified distance function), same high-level cure (isolate that factor, integrate the singular part in closed form, quadrature the smooth remainder with few nodes).

**How we differ, precisely:**
1. **Location:** they find the complex preimage φ* by **Newton iteration**; we get φ*, D_min, D''(φ*) in **closed form** because D is a trig polynomial (companion-matrix roots / critical points), no iteration.
2. **What is isolated:** they *swap* the singularity (map to the complexified preimage and expand the remainder in monomials); we *subtract* an explicit rational peak model m with an **elementary antiderivative**, keeping an **odd (linear) term** tuned to the √D_min → 0 grazing asymmetry.
3. **Differentiability:** they do not target autodiff/exact geometry-gradients; we make the construction **JAX-differentiable with provably exact gradients**, exploiting that closed-form (stop_gradient'd) φ* + the expansion-center-invariance of the subtraction identity yields machine-precision derivatives w.r.t. geometry through the near-singular quadrature.
4. **Kernel:** our specific target is the D^{-3/2} inverse-cube view/streaming kernel; theirs is developed for Laplace/Stokes-type 1/|Δ| line-integral kernels.

---

## The numerical-analysis reviewer's strongest objection

*"This is singularity-swap (af Klinteberg–Barnett) / preimage-subtraction (Montanelli et al.) re-skinned. Those methods already locate the near-singularity as a root of the complex distance function and integrate the singular factor in closed form; your only structural addition is that a trig-polynomial D makes the root closed-form instead of a 2–3-iteration Newton solve — a convenience, not a new capability. The autodiff angle is likewise near-standard (JAX-BEM, quadax already do differentiable near-singular quadrature), and freezing quadrature-related quantities with stop_gradient is routine. Show me a case where closed-form φ* changes the asymptotic cost or the achievable accuracy versus Newton-based singularity-swap at the same node count — and prove the exact-gradient claim is more than 'we differentiated through a converged Newton solve, which is also exact by the implicit function theorem.' Absent that, this is an application-tailored recombination, publishable as a methods note in a domain venue, not as a new NA technique in SISC/SINUM."*

**Rebuttals available (and their limits):**
- The **exact-gradient-via-invariance** argument is genuinely cleaner than implicit-differentiation of Newton and is worth foregrounding — but a strong referee will counter that differentiating a converged root-find is *also* exact, so this is elegance, not a new result.
- The **linear-term-at-grazing** analysis (c0 ~ √D_min while c1 = O(1)) is a concrete accuracy claim; if backed by a convergence figure showing the even-only model stalls where the odd model does not, it materially strengthens the paper.
- The **closed-form D_min/D''** genuinely removes a failure mode of Newton-based swap (bad initial guess / near-degenerate double roots at deep grazing where two complex preimages coalesce). If demonstrated at extreme grazing (min|Δ| → 0), this is the most defensible differentiator.

**Bottom line:** honest framing is a *novel, well-motivated combination and specialization* — publishable as such if positioned explicitly against singularity-swap and if the closed-form/differentiable advantages are demonstrated (not merely asserted) at deep grazing. Overclaiming it as a brand-new quadrature primitive would not survive review.
