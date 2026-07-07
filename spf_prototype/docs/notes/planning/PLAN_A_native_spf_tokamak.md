# PLAN A — Native SPF emission for `openmc::TokamakSource` (PR #3999)

**Goal.** Add spin-polarized-fusion (SPF) anisotropic neutron emission to Ethan
Peterson's native `TokamakSource` (PR #3999) by replacing its single isotropic
direction draw with our tier-1-validated Schwartz-P2 mode-mixture sampler about a
local field direction `b̂`. Backward-compatible: with no polarization input the
source is bit-for-bit the current isotropic behavior.

**Non-goal / hard constraint.** This document does **not** edit OpenMC core. It is
a design + concrete draft code that Ethan (or we, in a follow-up branch) can drop
into `openmc/source.py`, `include/openmc/source.h`, and `src/source.cpp`. All code
below is written against the *actual* signatures in the extracted PR sources.

Reference sources read for this plan:
- PR Python `TokamakSource(SourceBase)` — `scratchpad/pr3999/openmc_source.py`
- PR C++ header `class TokamakSource : public Source` — `scratchpad/pr3999/include_openmc_source.h`
- PR C++ `TokamakSource::sample()` — `scratchpad/pr3999/src_source.cpp` (lines 1148–1188; the isotropic step is line 1174)
- Validated sampler — `spf_prototype/src/spf_sampler.hpp`, `spf_field.hpp`, `polarized_fusion_source.cpp`

---

## 0. The one-line summary of the change

In `TokamakSource::sample()` the current birth-direction step is:

```cpp
// 5. Sample isotropic direction
site.u = angle_->sample(seed).first;
```

SPF replaces this single line with a branch: if the source is unpolarized, keep
the isotropic draw unchanged; otherwise sample the neutron direction from the
Schwartz P2 mode-mixture distribution about the **local magnetic-field direction
`b̂`** at the birth point. Everything else in `sample()` (position, energy, time,
weight) is untouched. That locality is the whole reason this is easy and mergeable.

The key enabling fact: at the point of the isotropic draw, the flux coordinates
`r` (minor radius, cm), `alpha` (poloidal angle, rad), and `phi` (toroidal angle,
rad) are **all still in local scope** (declared at lines 1156, 1160, 1163). So the
field direction — which depends on the flux-surface tangent — is computable with
zero extra plumbing.

---

## 1. The `b̂` model for a tokamak

The neutron angular distribution is azimuthally symmetric about the **local
magnetic-field direction** `b̂` (Schwartz 2025, Eq. 2; sampled by
`spf::sample_local_direction` in the `b̂`-aligned frame, then rotated to global by
the Gram–Schmidt `build_frame`). So SPF needs exactly one new physical input at
each birth point: the unit vector `b̂(r, α, φ)`.

In a tokamak the field is **tangent to the flux surface**: a dominant toroidal
part `B_φ φ̂` plus a poloidal part `B_θ p̂` that follows the flux-surface contour,
with the ratio set by the safety factor `q`.

### 1(a) Ship-first model: pure toroidal, `b̂ = φ̂`

```
φ̂(φ) = (−sin φ, cos φ, 0)
b̂    = φ̂
```

This is the `β = 0` limit and is **exactly** our existing `spf::ToroidalField`
(`spf_field.hpp`, lines 58–65), whose formula `{−y/r_xy, x/r_xy, 0}` is the
toroidal regression anchor validated to ~2e-14 against anarrima. Two decisive
advantages:

1. **It is the axisymmetric case Schwartz analyzed.** Schwartz's ±43% inboard /
   ∓22% outboard midplane NWL numbers assume `B̂ = φ̂`. So the ship-first model is
   the one the analytic oracle (anarrima) directly validates — no model mismatch
   in the headline test.
2. **Parameter-free.** No `q`-profile, no handedness, nothing to get wrong. It
   depends only on position, so it needs neither the flux coordinates nor any new
   plasma input.

### 1(b) Better model: toroidal + poloidal pitch from a `q`-profile

The physical field line has a poloidal tilt. Using the PR's own Miller
parameterization

```
R(r,α) = R0 + r·cos(ψ) + Δ·(1 − (r/a)²),   ψ = α + δ·sin α
Z(r,α) = κ·r·sin α
```

the **exact unit tangent to the flux-surface contour** (poloidal direction) is the
`α`-derivative at fixed `r`:

```
tR = ∂R/∂α = −r·sin(ψ)·(1 + δ·cos α)
tZ = ∂Z/∂α =  κ·r·cos α
tnorm = sqrt(tR² + tZ²)
p̂ = (tR/tnorm)·R̂ + (tZ/tnorm)·ẑ          # unit poloidal tangent, in the (R̂,ẑ) plane
```

with the poloidal-plane basis at toroidal angle `φ`:

```
R̂(φ) = (cos φ, sin φ, 0)
ẑ     = (0, 0, 1)
φ̂(φ) = (−sin φ, cos φ, 0)
```

The pitch magnitude comes from the safety factor. Using the standard
large-aspect-ratio (cylindrical) relation `q ≈ r·B_φ / (R·B_θ)`, evaluated with
the **local** major radius `R` at the birth point:

```
λ(r,α) = |B_θ| / |B_φ| = r / (q(r)·R(r,α))          # → 0 as q → ∞ (recovers toroidal)
```

The field direction (unnormalized, then normalized) is:

```
B = φ̂ + s·λ·p̂                                        # s = ±1 handedness (co-/counter-current)
b̂ = B / |B|
```

In explicit Cartesian components (this is what goes in the code):

```
pR = tR / tnorm
pZ = tZ / tnorm
Bx = −sin φ + s·λ·pR·cos φ
By =  cos φ + s·λ·pR·sin φ
Bz =  s·λ·pZ
b̂  = (Bx, By, Bz) / sqrt(Bx² + By² + Bz²)
```

Notes on correctness and scope:
- This is **exact geometry** for the poloidal tangent (from the parametric
  derivatives) and a **leading-order** pitch magnitude (`q_cyl`). A fully
  consistent `B_θ(α)` from `∇ψ` needs a real equilibrium and is out of scope for
  the prototype; the `q`-profile parameterization is the pragmatic middle ground
  and reduces cleanly to the toroidal limit.
- `λ ~ ε/q ~ 0.1–0.3` in practice, so the poloidal tilt is a **subleading**
  correction that slightly rotates the emission lobe in the poloidal plane and
  breaks exact up-down symmetry. It is physically real but small.
- `s` (handedness) only matters for the pitched model and only tilts the lobe; it
  has no effect in the toroidal limit.

### Recommendation

**Ship 1(a) (pure toroidal) first.** It is the validated, parameter-free,
oracle-matching case and captures the dominant SPF effect. Implement 1(b) behind a
`field_model='pitched'` switch as a fast follow, once the toroidal case is green
against anarrima. The code below wires both so 1(b) is a small, contained addition.

---

## 2. The direction hook — exact spot + concrete C++

### 2.1 Where

`src/source.cpp`, inside `TokamakSource::sample()`, replacing the single line at
1174. All of `r`, `alpha`, `phi`, and `site.r` are in scope.

### 2.2 New private state on `TokamakSource` (in `include/openmc/source.h`)

Add to the private data members of `class TokamakSource` (near the existing
`angle_` member). All set in the constructor, all const at sample time — no
mutable state, thread-safe by construction, matching the SPF discipline.

```cpp
  // --- SPF (spin-polarized fusion) emission -------------------------------
  bool polarized_ {false};        //!< false => isotropic (backward compatible)
  double w_perp_ {0.0};           //!< (3/4)·a          — sin²θ weight
  double w_par_  {0.0};           //!< (2/3)·b + (1/3)·c — (1/4+3/4cos²θ) weight
  double p_perp_ {0.0};           //!< a/η — stratified selection prob of sin²θ shape
  int    field_model_ {0};        //!< 0 = toroidal (b̂=φ̂), 1 = pitched (q-profile)
  double field_sign_ {1.0};       //!< s = ±1 handedness (pitched model only)
  vector<double> q_r_over_a_;     //!< r/a grid for q-profile (pitched model)
  vector<double> q_values_;       //!< q(r/a) values          (pitched model)
```

`η = a + (2/3)b + (1/3)c` is the total-rate factor. For a single tokamak source
with a user-set `strength`, η rescales absolute yield but **not** the normalized
NWL shape, so it need not touch the per-particle weight; see §5 open question on
whether to fold η into `strength`.

### 2.3 New private helpers (declarations in the header, near the other samplers)

```cpp
  //! Local magnetic-field DIRECTION b̂ at a birth point.
  //! \param r     minor radius [cm]
  //! \param alpha poloidal angle [rad]
  //! \param phi   toroidal angle [rad]
  //! \return unit field direction in global Cartesian
  Direction field_direction(double r, double alpha, double phi) const;

  //! Sample a birth direction from the Schwartz P2 mode-mixture about b̂.
  //! \param bhat local field direction (unit)
  //! \param seed PRN seed pointer
  //! \return unit neutron direction in global Cartesian
  Direction sample_polarized_direction(Direction bhat, uint64_t* seed) const;
```

### 2.4 The hook itself (replace line 1174)

```cpp
  // 5. Sample birth direction.
  //    Unpolarized  -> isotropic (bit-for-bit the previous behavior).
  //    Polarized    -> Schwartz P2 mode-mixture about the local field b̂.
  if (polarized_) {
    Direction bhat = field_direction(r, alpha, phi);
    site.u = sample_polarized_direction(bhat, seed);
  } else {
    site.u = angle_->sample(seed).first;
  }
```

### 2.5 `field_direction()` — implementation sketch (`src/source.cpp`)

```cpp
Direction TokamakSource::field_direction(
  double r, double alpha, double phi) const
{
  const double sphi = std::sin(phi), cphi = std::cos(phi);

  if (field_model_ == 0) {
    // Pure toroidal b̂ = φ̂ = (−sin φ, cos φ, 0). Already unit.
    return {-sphi, cphi, 0.0};
  }

  // Pitched model: toroidal + poloidal tangent scaled by local pitch λ = r/(qR).
  const double delta = triangularity_;
  const double psi = alpha + delta * std::sin(alpha);
  const double R = major_radius_ + r * std::cos(psi) +
                   shafranov_shift_ * (1.0 - (r * r) /
                     (minor_radius_ * minor_radius_));

  // Exact unit poloidal tangent from the parametric α-derivatives.
  const double tR = -r * std::sin(psi) * (1.0 + delta * std::cos(alpha));
  const double tZ = elongation_ * r * std::cos(alpha);
  const double tnorm = std::sqrt(tR * tR + tZ * tZ);
  // At r == 0 the tangent is degenerate; there the field is essentially toroidal.
  if (!(tnorm > 0.0))
    return {-sphi, cphi, 0.0};
  const double pR = tR / tnorm, pZ = tZ / tnorm;

  // Local pitch magnitude from the q-profile (q_cyl relation), s = handedness.
  const double q = interp_q(r / minor_radius_);   // linear interp on q_r_over_a_
  const double lambda = (q > 0.0) ? (r / (q * R)) : 0.0;
  const double sl = field_sign_ * lambda;

  // B = φ̂ + s·λ·p̂, expressed in Cartesian, then normalized.
  double Bx = -sphi + sl * pR * cphi;
  double By =  cphi + sl * pR * sphi;
  double Bz =  sl * pZ;
  const double Bnorm = std::sqrt(Bx * Bx + By * By + Bz * Bz);
  return {Bx / Bnorm, By / Bnorm, Bz / Bnorm};
}
```

(`interp_q` is a trivial linear interpolation over `q_r_over_a_`/`q_values_`,
clamped at the ends — one small private helper, or inline `std::lower_bound`.)

### 2.6 `sample_polarized_direction()` — native port of the validated sampler

This is a **line-for-line transcription** of `spf_sampler.hpp`
(`sample_local_direction` + `build_frame` + `to_global`), rewritten in
`openmc::Direction` and `openmc::prn`. The math is already tier-1 validated
(χ²/KS, inverse-CDF-vs-rejection, off-axis rotation isotropy, degenerate `b̂=±ẑ`),
so the port only needs a "port vs `spf_sampler.hpp` on identical seeds" parity test
(§4), mirroring the existing C++↔Python parity discipline.

```cpp
namespace {
// sin²θ shape: p(x)∝(1−x²). Inverse CDF root in [−1,1] of x³−3x+(4u−2)=0.
inline double spf_costheta_perp(double u)
{
  double x = 2.0 * std::cos(std::acos(1.0 - 2.0 * u) / 3.0 - 2.0 * PI / 3.0);
  return std::clamp(x, -1.0, 1.0);
}
// 1/4+3/4cos²θ shape: p(x)=1/4+3/4x². Single real root of x³+x+(2−4u)=0 (Cardano).
inline double spf_costheta_par(double u)
{
  const double q = 2.0 - 4.0 * u;
  const double sq = std::sqrt(q * q / 4.0 + 1.0 / 27.0);
  double x = std::cbrt(-q / 2.0 + sq) + std::cbrt(-q / 2.0 - sq);
  return std::clamp(x, -1.0, 1.0);
}
} // namespace

Direction TokamakSource::sample_polarized_direction(
  Direction bhat, uint64_t* seed) const
{
  // 1. cosθ relative to b̂ via stratified selection of the two P2 shapes.
  //    RNG draw order: selector, shape-u, phi  (matches spf_sampler.hpp).
  double x;
  if (prn(seed) < p_perp_) {
    x = spf_costheta_perp(prn(seed));
  } else {
    x = spf_costheta_par(prn(seed));
  }
  const double az = 2.0 * PI * prn(seed);
  const double st = std::sqrt(std::max(0.0, 1.0 - x * x));
  const Direction local {st * std::cos(az), st * std::sin(az), x};

  // 2. Rotate from the b̂-aligned local frame to global via Gram–Schmidt with a
  //    least-aligned reference axis (avoids the b̂ ∥ ẑ degeneracy).
  Direction b = bhat / bhat.norm();
  const double ax = std::abs(b.x), ay = std::abs(b.y), az2 = std::abs(b.z);
  Direction ref = (ax <= ay && ax <= az2) ? Direction{1, 0, 0}
                : (ay <= az2)             ? Direction{0, 1, 0}
                                          : Direction{0, 0, 1};
  Direction ex = ref - b.dot(ref) * b;
  ex /= ex.norm();
  Direction ey = b.cross(ex);

  Direction u = local.x * ex + local.y * ey + local.z * b;
  return u / u.norm();  // guard fp drift; |u|≈1
}
```

(`openmc::Direction` provides `.x/.y/.z`, `.dot()`, `.cross()`, `.norm()`, and the
arithmetic operators used above; `PI` and `prn` are already in scope in
`source.cpp`. If any operator is missing on the exact `Position`/`Direction` in the
target checkout, fall back to explicit component math exactly as in
`spf_sampler.hpp` — the two are interchangeable.)

### 2.7 Constructor wiring (`TokamakSource::TokamakSource`, `src/source.cpp`)

After the existing validation block, read the optional polarization inputs and
precompute the two angular weights + the stratified selection probability. Derive
constants — do not paste magic numbers (identities documented in `spf_sampler.hpp`
lines 22–31):

```cpp
  // --- SPF: optional polarization (absent => isotropic, unchanged behavior) ---
  if (check_for_node(node, "polarization")) {
    auto abc = get_node_array<double>(node, "polarization");  // [a, b, c]
    if (abc.size() != 3)
      fatal_error("TokamakSource: polarization must be 3 numbers a b c.");
    double a = abc[0], b = abc[1], c = abc[2];
    if (a < 0.0 || b < 0.0 || c < 0.0)
      fatal_error("TokamakSource: polarization fractions must be >= 0.");
    double s = a + b + c;
    if (s <= 0.0)
      fatal_error("TokamakSource: polarization a+b+c must be > 0.");
    a /= s; b /= s; c /= s;                       // renormalize (Python already warns)

    w_perp_ = 0.75 * a;                           // (3/4) a
    w_par_  = (2.0 / 3.0) * b + (1.0 / 3.0) * c;  // (2/3)b + (1/3)c
    const double eta = a + (2.0 / 3.0) * b + (1.0 / 3.0) * c;
    // P_perp = W_perp/(W_perp+W_par) with W_perp=w_perp·(8π/3)=2πa, W_par=w_par·2π.
    p_perp_ = a / eta;                            // == 2πa / (2πa + 2π w_par)
    // Invariant (not user input): P_perp is a valid probability. (η = a + w_par.)
    assert(p_perp_ >= -1e-15 && p_perp_ <= 1.0 + 1e-15 &&
           "TokamakSource SPF: P_perp must lie in [0,1]");
    polarized_ = true;

    // Field model (default toroidal). "pitched" additionally reads q-profile+sign.
    field_model_ = 0;
    if (check_for_node(node, "field_model")) {
      std::string fm = get_node_value(node, "field_model");
      if (fm == "toroidal") field_model_ = 0;
      else if (fm == "pitched") field_model_ = 1;
      else fatal_error("TokamakSource: field_model must be 'toroidal' or 'pitched'.");
    }
    if (field_model_ == 1) {
      q_r_over_a_ = get_node_array<double>(node, "q_r_over_a");
      q_values_   = get_node_array<double>(node, "q_values");
      if (check_for_node(node, "field_sign"))
        field_sign_ = std::stod(get_node_value(node, "field_sign"));
    }
  }
```

The existing `angle_ = UPtrAngle{new Isotropic()};` line stays as the
unpolarized/fallback path.

---

## 3. Python API — `openmc.TokamakSource`

Additive, backward-compatible. Default `polarization=None` ⇒ nothing written to
XML ⇒ C++ `polarized_ = false` ⇒ identical to today. Existing tests
(`test_source_tokamak.py`) pass unchanged.

### 3.1 New constructor parameters

```python
def __init__(
    self,
    major_radius, minor_radius, elongation, triangularity, shafranov_shift,
    r_over_a, emission_density, energy, time=None,
    phi_start=0.0, phi_extent=2.0*np.pi, n_alpha=101, vertical_shift=0.0,
    strength=1.0, constraints=None,
    # --- SPF additions -------------------------------------------------------
    polarization=None,           # None | (a, b, c) | dict of spin fractions
    field_model='toroidal',      # 'toroidal' (b̂=φ̂) | 'pitched' (q-profile)
    safety_factor=None,          # (r_over_a, q) tuple; required if field_model='pitched'
    field_sign=1,                # +1 / -1 handedness (pitched only)
):
    ...
    self.polarization = polarization   # setter normalizes to an (a,b,c) tuple or None
    self.field_model = field_model
    self.safety_factor = safety_factor
    self.field_sign = field_sign
```

### 3.2 Polarization setter — accept `(a,b,c)` OR spin fractions

Canonical stored form is the `(a, b, c)` triple (matches the C++ XML contract and
`abc_modes.py`). A dict of deuteron/triton spin fractions is converted with the
Schwartz Eq. 1 map, so users can specify the physical knobs directly:

```python
@property
def polarization(self):
    return self._polarization  # None or (a, b, c)

@polarization.setter
def polarization(self, value):
    if value is None:
        self._polarization = None
        return
    if isinstance(value, dict):
        # Schwartz 2025 Eq. 1: deuteron d+,d0,d- ; triton t+,t- (each set sums to 1)
        dp, d0, dm = value['d_plus'], value['d_zero'], value['d_minus']
        tp, tm = value['t_plus'], value['t_minus']
        a = dp*tp + dm*tm
        b = d0                    # = d0*(tp+tm)
        c = dp*tm + dm*tp
    else:
        a, b, c = value           # (a, b, c) triple
    for name, x in (('a', a), ('b', b), ('c', c)):
        cv.check_greater_than(f'polarization {name}', x, 0.0, equality=True)
    s = a + b + c
    if abs(s - 1.0) > 1e-6:       # user input => warn (not assert); C++ renormalizes too
        warnings.warn(f'polarization (a,b,c) sums to {s}, renormalizing to 1.')
        a, b, c = a/s, b/s, c/s
    self._polarization = (a, b, c)
```

Convenience module helper (mirrors `spf_prototype/python/abc_modes.py`, sympy-verified):

```python
def spin_fractions_to_abc(d_plus, d_zero, d_minus, t_plus, t_minus):
    """Schwartz 2025 Eq. 1 mode fractions from deuteron/triton spin fractions."""
    a = d_plus*t_plus + d_minus*t_minus
    b = d_zero
    c = d_plus*t_minus + d_minus*t_plus
    return a, b, c
```

### 3.3 XML round-trip

`populate_xml_element()` — append only when polarized (keeps existing files clean):

```python
    if self.polarization is not None:
        a, b, c = self.polarization
        ET.SubElement(element, "polarization").text = f"{a} {b} {c}"
        ET.SubElement(element, "field_model").text = self.field_model
        if self.field_model == 'pitched':
            r_q, q = self.safety_factor
            ET.SubElement(element, "q_r_over_a").text = ' '.join(map(str, r_q))
            ET.SubElement(element, "q_values").text = ' '.join(map(str, q))
            ET.SubElement(element, "field_sign").text = str(self.field_sign)
```

`from_xml_element()` — read them back (all optional; absence ⇒ unpolarized):

```python
    pol_text = get_text(elem, 'polarization')
    polarization = tuple(float(x) for x in pol_text.split()) if pol_text else None
    field_model = get_text(elem, 'field_model') or 'toroidal'
    # ... q_r_over_a / q_values / field_sign analogous, only if present ...
```

### 3.4 User-facing example

```python
import openmc
# Unpolarized (unchanged): polarization omitted -> isotropic.
src = openmc.TokamakSource(
    major_radius=620.0, minor_radius=200.0, elongation=1.8, triangularity=0.45,
    shafranov_shift=10.0, r_over_a=r, emission_density=1 - r**2,
    energy=openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=2.0e4))

# Fully perpendicular-polarized (A mode), toroidal field: peaks the NWL inboard.
src_A = openmc.TokamakSource(..., polarization=(1.0, 0.0, 0.0))

# From physical spin fractions, with a pitched field from a q-profile.
a, b, c = openmc.spin_fractions_to_abc(d_plus=0.9, d_zero=0.05, d_minus=0.05,
                                       t_plus=0.9, t_minus=0.1)
src_phys = openmc.TokamakSource(..., polarization=(a, b, c),
                                field_model='pitched',
                                safety_factor=(r_grid, q_grid), field_sign=+1)
```

---

## 4. Validation plan

Gate order mirrors the existing prototype tiers: unit parity first, then the
analytic NWL headline.

### 4.1 Sampler parity (cheap, first gate)

The native `sample_polarized_direction` port must reproduce `spf_sampler.hpp`
`sample_global_direction` **bit-for-bit on a fixed seed sequence** (same RNG draw
order: selector → shape-u → phi). Reuse `tests/test_field_parity.py` /
`tests/test_sampler_stats.py` machinery: draw ~1e6 directions for {non-pol, A, B,
C, mixed (0.5,0.3,0.2)} with `b̂=ẑ`, χ² against `spf::pdf_costheta`, and confirm
A vanishes at cosθ=±1, B/C peak there, non-pol is flat. This re-uses the already
green tier-1 oracles — it only proves the transcription into `openmc::Direction`
is faithful.

### 4.2 Headline: axisymmetric NWL vs Schwartz / anarrima

Configuration that makes the tokamak source reproduce Schwartz's axisymmetric
example:
- Geometry knobs for pure axisymmetry: `triangularity=0`, `shafranov_shift=0`,
  `phi_extent=2π`; circular plasma at aspect ratio `R0/a = 2.5` (Schwartz §2),
  `elongation=1`. Square-cross-section vacuum vessel around it (reuse the tier-2
  vessel from `spf_prototype/python/build_and_run.py`).
- **Field model = toroidal** (ship-first 1(a)) — this is exactly `B̂=φ̂`, the case
  the oracle assumes.
- Scattering off: near-void interior / current tally so neutrons free-stream ring
  → wall (reuse the tier-2 near-void material; confirm mean free path ≫ device).
- Poloidal surface-current tally on the first wall, ≥50 bins/wall, each pattern
  normalized to its own integral (compare shapes).

Configs: {isotropic, pure-A `(1,0,0)`, pure-B `(0,1,0)`, pure-C `(0,0,1)`, mixed
`(0.5,0.3,0.2)`}.

**Oracle = anarrima** (Schwartz's own package), already wired as the primary
analytic reference in `spf_prototype/python/analytic_nwl.py` /
`verify_nwl_analytic.py`. The TokamakSource samples the plasma *volumetrically*, so
compare against anarrima's ring-superposition weighted by the same
`emission_density × flux-Jacobian` the source uses (the volume source = weighted
sum of filamentary rings). Do it two ways and require agreement: (a) our own φ
quadrature, (b) anarrima.

**Pass/fail oracles (must reproduce from the OpenMC run, not just analytics):**
- Pure-A: **+43% inboard** midplane, **−22% outboard** midplane vs isotropic.
- Pure-B / pure-C: the mirror image (−43% inboard, +22% outboard).
- **Isotropic ≡ pure-C normalized shape** (B and C share directionality; differ
  only in total rate η) — a direct test of the §1.2 correction.
- **Linearity**: normalized `(0.5,0.3,0.2)` pattern = `0.5·A + 0.3·B + 0.2·C` of
  the pure-mode normalized patterns, within statistics.
- Residual normality: `(Q_MC − Q_analytic)/σ_MC` over wall bins ~ N(0,1); a
  coherent nonzero mean flags a real bug (wrong constant / rotation), not noise.

### 4.3 Pitched model (follow-up)

Once toroidal is green, run `field_model='pitched'` with a flat `q` profile and
confirm: (i) as `q → ∞` the pattern converges to the toroidal result; (ii) a
finite `q` produces the expected small poloidal rotation / up-down asymmetry of the
lobe, consistent with anarrima's angled-field kernels (`g_HAa/g_Hca/g_HBa`, …), the
same ones `spf_field.hpp`'s `AngledField` is pinned to. There is no external oracle
for a full q-profile equilibrium field — that is out of scope; the check is
internal consistency + the toroidal limit.

---

## 5. Why this is easy/mergeable + open questions for Ethan

### Why it is low-risk and mergeable

- **Surgical.** One branch at the direction step; everything else in `sample()`
  untouched. No change to position/energy/time sampling, CDFs, or XML for existing
  users.
- **Zero behavior change by default.** `polarization=None` ⇒ no XML ⇒
  `polarized_=false` ⇒ the exact current isotropic draw. Existing
  `test_source_tokamak.py` passes unmodified.
- **Reuses validated math.** The angular sampler, the stratified selection, the
  Gram–Schmidt rotation, and the toroidal field are already tier-1/tier-2 green in
  `spf_prototype`; the native version is a transcription guarded by a parity test.
- **No new dependencies, no mutable state.** All SPF members set in the
  constructor, const at sample time; `openmc::prn` only — same thread-safety
  contract as the rest of `TokamakSource`.
- **Additive API.** New Python kwargs default to the unpolarized path; new XML
  elements are optional and round-trip cleanly.
- **The physics knob is tiny.** Three numbers `(a,b,c)` (or five spin fractions)
  plus a field-model switch. Everything else is derived.

### Open questions for Ethan

1. **Scope: this PR or a follow-up?** SPF is a clean additive layer on #3999.
   Land it in #3999 (one reviewer, one geometry) or as an immediate follow-up PR
   on top? Recommendation: follow-up PR so #3999 merges on its own schedule, but
   the hook (§2.4) can go in now as a no-op branch to minimize later churn.
2. **`b̂` model to ship.** Recommend toroidal-only first (validated, parameter-free,
   oracle-matching). Is a `q`-profile input acceptable for the pitched follow-up,
   or would he prefer to defer any poloidal pitch until a real equilibrium field
   map is available?
3. **Where should the P2 sampler live?** As private `TokamakSource` methods
   (simplest, shown here), or factored into a small shared util
   (`math_functions`/a new `polarized_source.h`) so a future general/compiled
   polarized source can reuse it? The latter is more work but avoids a second copy.
4. **Total-rate factor `η`.** Should `η = a + (2/3)b + (1/3)c` scale the source
   `strength` (so A-mode vs C-mode fuel emit physically different absolute yields),
   or stay purely a shape sampler (η ignored)? Irrelevant to normalized NWL shape;
   matters for absolute yields and multi-source normalization. Recommend exposing
   it as an opt-in (`scale_strength_by_eta=False` default) to avoid surprising
   existing yield normalizations.
5. **Spin fractions vs `(a,b,c)` in the public API.** Store `(a,b,c)` canonically
   (done), but is the `spin_fractions_to_abc` helper + dict form worth carrying, or
   keep the API to the three mode fractions only?

### Explicitly out of scope (state in `LIMITATIONS`)

- **Depolarization** in transport/along field lines — birth distribution only.
- **Energy–angle correlation** — energy stays decoupled (monoenergetic / Ballabio
  as today).
- **Alpha channel, non-axisymmetric / real-equilibrium `b̂`** — pitched model is a
  `q_cyl` parameterization, not a solved equilibrium.
- **Convex-wall / free-stream NWL validation only** — scattering is off in the
  headline comparison, matching Schwartz's own caveats.

---

## Appendix — identity crib (all sympy-verified in `spf_prototype`)

```
Schwartz Eq. 2:  dσ/dΩ = (σ0/2π)[ (3/4)a·sin²θ + ((2/3)b+(1/3)c)(1/4+(3/4)cos²θ) ]
Angular weights: w_perp = (3/4)a ;  w_par = (2/3)b + (1/3)c
Solid-angle integrals: ∫sin²θ dΩ = 8π/3 ;  ∫(1/4+3/4cos²θ) dΩ = 2π
Total-rate factor:     η = σ_tot/σ0 = a + (2/3)b + (1/3)c
  corners:  nonpol(1/3,1/3,1/3)→2/3 ;  A(1,0,0)→1 ;  B(0,1,0)→2/3 ;  C(0,0,1)→1/3
Stratified selection:  W_perp = w_perp·(8π/3) = 2πa ;  W_par = w_par·2π
                       P_perp = W_perp/(W_perp+W_par) = a/η
B and C share the (1/4+3/4cos²θ) angular shape; they differ ONLY in η (rate),
so isotropic and pure-C have identical normalized wall patterns.
```
```
Miller flux surface (PR #3999):
  R(r,α) = R0 + r·cos(ψ) + Δ(1−(r/a)²),  ψ = α + δ·sin α ;  Z(r,α) = κ·r·sin α
Poloidal tangent (∂/∂α):  tR = −r·sin ψ·(1+δ cos α) ;  tZ = κ·r·cos α
Local pitch (q_cyl):       λ = |Bθ|/|Bφ| = r/(q·R)
Field direction:           b̂ ∝ φ̂ + s·λ·(pR·R̂ + pZ·ẑ),  p̂ = (tR,tZ)/‖·‖
```
```
```
