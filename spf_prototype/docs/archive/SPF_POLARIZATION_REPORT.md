# Spin-polarized fusion emission for OpenMC plasma sources

Author: Thaddaeus Kiker (tjk2147@columbia.edu)
Branch: `spf-polarization` (one commit on top of `stellarator_source` @ `cf616f432`)
Scope: source, tests, and docs. The representation study in Section 5 is research context; its
study code and validation harness are not in this PR, only summarized here.

---

## 1. What this adds

A spin-polarized D-T birth-direction sampler (Schwartz 2025), wired into both plasma source
classes through one shared sampler. With no polarization input, the direction sampling and RNG
draw order are identical to before.

- `spin_fractions_to_abc(d_plus, d_zero, d_minus, t_plus, t_minus)`: maps deuteron and triton
  spin-projection fractions to the three mode weights `(a, b, c)`.
- `TokamakSource`: optional `polarization` with `field_model='toroidal'` or `'pitched'`
  (Miller-surface poloidal tilt from a supplied `safety_factor` profile).
- `StellaratorSource`: optional `polarization` with the local field direction reconstructed on
  the fly from the rotational transform `iota` and the `lambda` stream function. `from_vmec`
  reads these from the wout; `from_desc` extracts them from the DESC Fourier-Zernike basis and
  iota profile.
- C++ sampler: closed-form inverse-CDF sampling of the two P2 angular shapes, with a
  least-aligned-axis Gram-Schmidt frame rotation. One helper serves both classes.

### Physics

Emission about the local field `B̂` follows Schwartz Eq. 2:

```
dσ/dΩ ∝ (3/4)·a·sin²θ  +  ((2/3)·b + (1/3)·c)·(1/4 + (3/4)·cos²θ)
```

The birth-direction pdf then needs only two angular weights, `w_perp = (3/4)a` (the `sin²θ`
shape) and `w_par = (2/3)b + (1/3)c` (the `1/4 + 3/4 cos²θ` shape). B and C modes share the
parallel shape and differ only in total rate. The total-rate factor `η = a + (2/3)b + (1/3)c`
is kept separate from the direction sampling; it scales source strength, not the pdf.

---

## 2. Tests (56 unit tests, all passing)

`tests/unit_tests/test_source_tokamak.py` and `test_source_stellarator.py`:

- Mode algebra: `spin_fractions_to_abc` corners; unpolarized gives `(1/3, 1/3, 1/3)`.
- Input validation: negative/NaN/out-of-range polarization; bad safety-factor grids; polarized
  stellarator requires `field_rho/iota/lmns`; pitched requires `safety_factor`.
- XML round-trips for every new field, and confirmation that an unpolarized source writes no
  SPF elements (backward compatibility).
- Compiled angular-moment recovery (real `sample_external_source`): about the local `B̂`,
  `⟨cos²θ⟩ = 1/5` (A), `1/3` (unpolarized), `7/15` (B/C), for the tokamak (toroidal and
  q-pitched) and the stellarator (with a 3-D lambda harmonic so a toroidal-only implementation
  cannot pass).
- VMEC half-mesh handling (classic wout lambda/iota row-0 padding) and the DESC HDF5
  iota-profile reader.

### Field-direction validation against an independent oracle

The reconstructed `B̂` is checked against each code's own computed field, not against the
formula that produced it:

| Equilibrium set | Oracle | Field-direction error |
|---|---|---|
| VMEC QA / QH | exported contravariant `bsupumnc/bsupvmnc` | p99 = 1.3e-4° (QA), 3.8e-4° (QH) |
| DESC precise_QA | `eq.compute("B")` | median 0.031° to 0.010° at n_rho = 16 to 49 |
| DESC precise_QH | `eq.compute("B")` | median 0.049° to 0.015° at n_rho = 16 to 49 |

The DESC error converges with radial resolution (a sign or convention bug would sit near 90°),
which confirms the convention.

### End-to-end compiled sampling from a DESC equilibrium

Full chain `from_desc(polarization=…) -> XML -> C++ constructor -> C++ field_direction -> C++
sampler`. Each Cartesian birth is inverted back to flux coordinates (Newton, residual ~1e-13
cm) to recover the local `B̂` and measure `⟨cos²θ⟩`:

| Device | perp | par | unpolarized |
|---|---|---|---|
| precise_QA | 0.1998 | 0.4680 | 0.3330 |
| precise_QH | 0.1994 | 0.4651 | 0.3327 |
| analytic | 0.2000 | 0.4667 | 0.3333 |

---

## 3. VMEC vs DESC integration

Both feed the same compiled spectral field path; only the extraction differs.

- VMEC (`from_vmec`): `iota` and the `lambda` coefficients (`lmns`/`lmnc`) are read from the
  wout. The classic half-mesh (row-0 padding on `lmns`/`iotas`) and the VMEC++ full-mesh form
  are both handled, on the native radial mesh.
- DESC (`from_desc`): `lambda` (`L_lmn` on `L_basis`) is Zernike-collapsed onto the rho grid
  and Ptolemy-folded into the same combined VMEC-style mode table as the geometry (the mode
  union is extended to include lambda modes). `iota` is evaluated on rho (live equilibrium via
  DESC `compute`; HDF5 via the stored power-series profile, reading the even-power basis from
  the profile rather than assuming it). Validated as in Section 2 on QA and QH.

The field-line relation is `B^θ ∝ iota − ∂λ/∂ζ`, `B^ζ ∝ 1 + ∂λ/∂θ`, evaluated in the
coordinate basis to avoid dividing by a small `1 + ∂λ/∂θ`.

---

## 4. Spin-polarized emission across source types

The kernel is source-agnostic; it steers about whatever `B̂` the source provides.

| Source | Field direction | Status |
|---|---|---|
| Tokamak | toroidal, or q-pitched Miller poloidal tilt | in this PR, validated |
| Stellarator, analytic (spectral) | on-the-fly from VMEC/DESC `iota` + `lambda` | in this PR, validated |
| Stellarator, interpolated (fluxmap grid) | precomputed `(B_R, B_φ, B_Z)` on a grid | same kernel; used in my transport work (Section 5), not part of this PR |

---

## 5. Representation study: interpolated vs analytic vs adaptive (research context)

Three representations of the flux-surface data the source needs at each birth (position `R, Z`,
birth weight `√g`, field direction `B̂`):

- Analytic (spectral): keep the VMEC/DESC Fourier coefficients and evaluate `R, Z, √g, B̂` on
  the fly. Exact to the equilibrium's own truncation, small memory, one Fourier sum per
  evaluation.
- Interpolated (grid): precompute `R, Z, √g, B̂` on a `(ρ, θ, ζ)` grid and trilinearly
  interpolate. O(1) per evaluation, with interpolation error and larger memory. This is the
  representation my transport work uses.
- Adaptive interpolation: the same grid, with node placement driven by the spectral truth
  (radial nodes equidistributed on `√g` curvature; angular nodes on a field-aware monitor that
  clusters where `B̂` rotates with the field-line pitch).

### 5.1 Single device (886079), head-to-head on one machine

| Representation | Speed (M evals/s) | Memory | Field-dir error (source-wtd p99) |
|---|---|---|---|
| Spectral (analytic) | 0.66 | 0.09 MiB | 0° (exact) |
| Grid 8×32×96 | 2.25 | 1.1 MiB | 15.6° |
| Grid 24×96×288 (production) | 1.68 | 30.4 MiB | 0.92° |
| Adaptive 24×48×96 | 2.08 | 5.1 MiB | 0.18° |

The grid is 2-3× faster per evaluation than the Fourier sum but uses 12-330× more memory, and
the coarse fast grid is 15.6° off in field direction. The adaptive grid matches the grid's
speed, uses 6× less memory than the production grid, and is close to spectral accuracy. (886079
is a strong steerer, so its coarse-grid field error sits at the high end of the ensemble in
5.2.)

### 5.2 20-device QUASR ensemble, median [IQR]

Grid-vs-spectral over 20 QA/QH devices. Spectral baseline: 407 k evals/s, 0.093 MiB, exact.

| Grid (ρ×θ×ζ) | Memory | Speed | √g RMS error | Field-dir error p99 |
|---|---|---|---|---|
| 8×32×96 | 1.1 MiB | 1.28 M/s | 0.211 [0.196, 0.218] | 5.36° [1.28, 11.46] |
| 12×48×144 | 3.8 MiB | 1.20 M/s | 0.168 [0.157, 0.174] | 0.98° [0.50, 1.37] |
| 16×64×192 | 9.0 MiB | 1.15 M/s | 0.145 [0.135, 0.149] | 0.62° [0.40, 0.83] |
| 24×96×288 | 30.4 MiB | 1.06 M/s | 0.119 [0.110, 0.122] | 0.44° [0.28, 0.58] |

### 5.3 Jacobian error vs grid resolution

`√g` sets the birth density. Its interpolation error is axis-dominated, so refining the grid
past production resolution helps little while memory grows quickly:

| Grid | Memory | √g RMS | Field-dir error p99 |
|---|---|---|---|
| 24×96×288 | 30.4 MiB | 0.119 | 0.44° |
| 32×128×384 | 72.0 MiB | 0.102 | 0.33° |
| 40×160×480 | 140.6 MiB | 0.092 | 0.26° |
| 48×192×576 | 243.0 MiB | 0.088 | 0.20° |

Going from 30 to 243 MiB (8×) drops the `√g` error only from 0.119 to 0.088. Grid refinement
does not reach spectral accuracy.

### 5.4 End-to-end tally impact

The polarized source drives a free-streaming first-wall loading tally. Gridded and spectral
wall patterns are compared bin by bin with a per-bin z-score, so a real bias shows up as a
coherent mean that grows like √N while noise stays flat.

- The coarse grid has a coherent under-steering bias. On device 886079 the 8³ grid's inboard
  bias grows −0.57σ, −0.90σ, −1.60σ at 2 M, 8 M, 20 M histories (max |z| = 5.3 at 20 M), which
  is √N growth rather than noise. The 30 MiB 24³ grid stays near the noise floor (|z| ≈ 0.3,
  max ≈ 2.6).
- The bias tracks steering strength. Across the 20 devices, the coarse-grid wall bias
  correlates with how strongly each device steers: Pearson r = +0.83 for the 8³ grid, r = +0.54
  at 24³. The grid error is largest in the high-steering cases where polarization matters most.
- Reproduced by the transport engine. An OpenMC `FileSource` run of the same births shows the
  8³ grid under-predicting inboard load by −0.82 % and over-predicting outboard by +0.79 %; the
  24³ grid is −0.004 % / +0.22 %.

### 5.5 Adaptive interpolation

Node placement driven by the spectral truth:

- Radial (`√g`) adaptivity alone, at equal node count, cuts `√g` RMS about 3.2× and the
  source-weighted field error about 6× (0.92° to 0.15° at 30 MiB) versus the uniform grid.
- Angular adaptivity needs a field-aware monitor. Geometry curvature is nearly flat in θ/ζ for
  these near-QS surfaces, but the field direction rotates in ζ; adding `B̂` curvature to the
  monitor lets a 24×48×96 grid (5.1 MiB) reach 0.18° field error, matching the uniform 30 MiB
  production grid at 6× less memory.
- End-to-end: at equal 30 MiB the adaptive grid reduces the coherent wall bias from 0.32σ/0.45σ
  (inboard/outboard) to 0.04σ/0.17σ, near the noise floor, and it beats the uniform 61 MiB
  (48³) grid at half the memory.

### 5.6 Summary

The spectral source is exact and compact but costs a Fourier sum per draw. The interpolated
grid gives O(1) evaluation with a fixed accuracy floor: 12-330× more memory, a `√g` error that
does not refine away, and a coherent ~1 % wall-load bias that scales with steering strength.
The spectral-driven adaptive grid reaches spectral accuracy at grid speed and about 6× less
memory. Numbers and figures are reproducible; scripts available on request.

---

## 6. Backward compatibility

With no `polarization` argument, both sources are byte-identical to the current behavior: the
isotropic direction branch and its random-number draw order are unchanged. Verified by the
unpolarized round-trip test and the gated C++ sampler.
