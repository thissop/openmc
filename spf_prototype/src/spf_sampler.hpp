// spf_sampler.hpp -- Pure-math core of the spin-polarized DT neutron birth-
// direction sampler (Schwartz 2025). Header-only; depends only on <cmath>,
// <cassert>, <stdexcept>. It is deliberately FREE of any OpenMC type so it can
// be (a) compiled into a standalone tier-1 driver with no libopenmc, and
// (b) mirrored 1:1 in Python (python/spf_mirror.py). The OpenMC `Source`
// wrapper (tier 2) includes this header and adapts openmc::Position/prn to the
// Vec3/Rng interface below.
//
// Randomness is injected via a callable `rng()` returning a double in [0,1)
// (OpenMC's prn). The sampler is otherwise pure and const: no mutable state.
//
// ----------------------------------------------------------------------------
// Physics (Schwartz 2025, Eqs. 1-2), with all constants DERIVED, not pasted:
//
//   dσ/dΩ = (σ0/2π) [ (3/4) a sin²θ + ((2/3)b + (1/3)c)(1/4 + (3/4)cos²θ) ]
//
//   θ = angle to the local magnetic field B̂. Define the two angular weights
//       w_perp = (3/4) a                 (the sin²θ piece)
//       w_par  = (2/3) b + (1/3) c       (the 1/4 + 3/4 cos²θ piece)
//   Unnormalized birth intensity on the sphere: I(θ) = w_perp sin²θ + w_par(1/4+3/4cos²θ).
//
//   Solid-angle integrals (verify in tests/test_abc_modes.py with sympy):
//       ∫ sin²θ dΩ          = 8π/3
//       ∫ (1/4 + 3/4 cos²θ) dΩ = 2π
//   ⇒ total cross section / σ0 = a + (2/3)b + (1/3)c ≡ η  (the total-RATE factor)
//     Corners: nonpol(1/3,1/3,1/3)→2/3,  A(1,0,0)→1,  B(0,1,0)→2/3,  C(0,0,1)→1/3.
//
//   Stratified angular sampling: the direction pdf depends ONLY on (w_perp,w_par).
//   Solid-angle-integrated weights  W_perp = w_perp·(8π/3) = 2π a,
//                                   W_par  = w_par ·(2π)    = 2π·((2/3)b+(1/3)c).
//   ⇒ P_perp = W_perp/(W_perp+W_par) = a/η.  (rate η is carried SEPARATELY via strength.)
// ----------------------------------------------------------------------------
#ifndef SPF_SAMPLER_HPP
#define SPF_SAMPLER_HPP

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace spf {

constexpr double PI = 3.14159265358979323846;

// Closed-form solid-angle integrals of the two angular shapes (documented above).
constexpr double INT_SIN2 = 8.0 * PI / 3.0;  // ∫ sin²θ dΩ
constexpr double INT_PAR = 2.0 * PI;         // ∫ (1/4 + 3/4 cos²θ) dΩ

// --------------------------------------------------------------------------
// Minimal 3-vector (kept separate from openmc::Position so the header stands
// alone; the tier-2 wrapper converts at the boundary).
// --------------------------------------------------------------------------
struct Vec3 {
  double x {0.0}, y {0.0}, z {0.0};
};

inline double dot(const Vec3& a, const Vec3& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a)
{
  return std::sqrt(dot(a, a));
}
inline Vec3 operator*(double s, const Vec3& a)
{
  return {s * a.x, s * a.y, s * a.z};
}
inline Vec3 operator+(const Vec3& a, const Vec3& b)
{
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 operator-(const Vec3& a, const Vec3& b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 normalized(const Vec3& a)
{
  double n = norm(a);
  if (!(n > 0.0))
    throw std::invalid_argument("spf: cannot normalize a zero/NaN vector");
  return (1.0 / n) * a;
}

// --------------------------------------------------------------------------
// Mode weights: (a,b,c) -> (w_perp, w_par, eta, P_perp). Input is renormalized
// to sum 1 (silently here; the OpenMC wrapper does the user-facing warn-on-
// large-deviation check). Invariants are hard-asserted, not warned.
// --------------------------------------------------------------------------
struct ModeWeights {
  double a, b, c;   // renormalized mode fractions (sum to 1)
  double w_perp;    // (3/4) a
  double w_par;     // (2/3) b + (1/3) c
  double eta;       // a + (2/3) b + (1/3) c  (= σ_tot/σ0, the total-rate factor)
  double P_perp;    // a/eta  (stratified selection prob of the sin²θ shape)
};

inline ModeWeights make_mode_weights(double a, double b, double c)
{
  // Input validation (user data) -> throw. NaNs fail the >= 0 test.
  if (!(a >= 0.0) || !(b >= 0.0) || !(c >= 0.0))
    throw std::invalid_argument("spf: a,b,c must each be >= 0");
  const double s = a + b + c;
  if (!(s > 0.0))
    throw std::invalid_argument("spf: a+b+c must be > 0");
  a /= s;
  b /= s;
  c /= s;

  ModeWeights m;
  m.a = a;
  m.b = b;
  m.c = c;
  m.w_perp = 0.75 * a;
  m.w_par = (2.0 / 3.0) * b + (1.0 / 3.0) * c;
  m.eta = a + (2.0 / 3.0) * b + (1.0 / 3.0) * c;

  // P_perp from solid-angle-integrated weights (closed form == a/eta).
  const double W_perp = m.w_perp * INT_SIN2;  // = 2π a
  const double W_par = m.w_par * INT_PAR;      // = 2π w_par
  const double W_tot = W_perp + W_par;
  assert(W_tot > 0.0 && "eta must be positive");
  m.P_perp = W_perp / W_tot;

  // Invariant: selection probabilities sum to 1 (P_par = 1 - P_perp).
  assert(std::fabs(m.P_perp + (W_par / W_tot) - 1.0) < 1e-12 &&
         "selection probabilities must sum to 1");
  assert(m.P_perp >= -1e-15 && m.P_perp <= 1.0 + 1e-15 && "P_perp in [0,1]");
  return m;
}

// --------------------------------------------------------------------------
// Per-shape marginal samplers in x = cosθ. Both inverse-CDF forms were derived
// and endpoint-verified (u=0 -> x=-1, u=1 -> x=+1); see tests for KS cross-
// check against the rejection variants.
// --------------------------------------------------------------------------

// sin²θ shape: p(x) ∝ (1 - x²). Inverse CDF solves x³ - 3x + (4u - 2) = 0.
// Three real roots (|q|<=2); the one in [-1,1] is the k=1 trig root:
//   x = 2 cos( (1/3) arccos(1 - 2u) - 2π/3 ).
inline double sample_costheta_perp_invcdf(double u)
{
  double x = 2.0 * std::cos(std::acos(1.0 - 2.0 * u) / 3.0 - 2.0 * PI / 3.0);
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;  // guard fp drift exactly at the endpoints
  return x;
}

// 1/4 + 3/4 cos²θ shape: p(x) = 1/4 + 3/4 x² (already normalized on [-1,1]).
// Inverse CDF solves x³ + x + (2 - 4u) = 0 (p=1>0 -> single real root, Cardano).
inline double sample_costheta_par_invcdf(double u)
{
  const double q = 2.0 - 4.0 * u;          // depressed cubic x³ + 1·x + q = 0
  const double disc = q * q / 4.0 + 1.0 / 27.0;  // > 0 always
  const double sq = std::sqrt(disc);
  // std::cbrt gives the real cube root (handles negative args correctly).
  double x = std::cbrt(-q / 2.0 + sq) + std::cbrt(-q / 2.0 - sq);
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;
  return x;
}

// Rejection variants (envelope = uniform x on [-1,1]); used only to cross-check
// the inverse-CDF forms statistically (two-sample KS).
template<class Rng>
double sample_costheta_perp_reject(Rng& rng)
{
  for (;;) {
    double x = 2.0 * rng() - 1.0;
    if (rng() <= (1.0 - x * x))  // accept ∝ (1-x²), max 1 at x=0  (eff 2/3)
      return x;
  }
}
template<class Rng>
double sample_costheta_par_reject(Rng& rng)
{
  for (;;) {
    double x = 2.0 * rng() - 1.0;
    if (rng() <= (0.25 + 0.75 * x * x))  // accept ∝ 1/4+3/4x², max 1 at |x|=1
      return x;
  }
}

// --------------------------------------------------------------------------
// Local-frame direction sample (ẑ_local ∥ B̂). RNG draw order (invcdf path):
//   1) selector, 2) shape u, 3) phi  -- this fixed order is what the Python
//   mirror replicates for the bit-for-bit parity test.
// --------------------------------------------------------------------------
template<class Rng>
Vec3 sample_local_direction(const ModeWeights& m, Rng& rng, bool use_rejection = false)
{
  double x;  // cosθ relative to B̂
  if (rng() < m.P_perp) {
    x = use_rejection ? sample_costheta_perp_reject(rng)
                      : sample_costheta_perp_invcdf(rng());
  } else {
    x = use_rejection ? sample_costheta_par_reject(rng)
                      : sample_costheta_par_invcdf(rng());
  }
  const double phi = 2.0 * PI * rng();
  const double st = std::sqrt(std::fmax(0.0, 1.0 - x * x));
  return {st * std::cos(phi), st * std::sin(phi), x};
}

// --------------------------------------------------------------------------
// Local -> global rotation via least-aligned reference axis + Gram-Schmidt.
// Avoids the B̂ × ẑ_world degeneracy when B̂ ∥ ẑ.
// --------------------------------------------------------------------------
struct Frame {
  Vec3 ex, ey, ez;  // ez = B̂
};

inline Frame build_frame(const Vec3& Bhat)
{
  const Vec3 b = normalized(Bhat);
  const double ax = std::fabs(b.x), ay = std::fabs(b.y), az = std::fabs(b.z);
  Vec3 ref;
  if (ax <= ay && ax <= az)
    ref = {1.0, 0.0, 0.0};
  else if (ay <= az)
    ref = {0.0, 1.0, 0.0};
  else
    ref = {0.0, 0.0, 1.0};
  const Vec3 ex = normalized(ref - dot(ref, b) * b);
  const Vec3 ey = cross(b, ex);
  return {ex, ey, b};
}

inline Vec3 to_global(const Vec3& local, const Frame& f)
{
  return local.x * f.ex + local.y * f.ey + local.z * f.ez;
}

// Full birth-direction sample in global Cartesian. Asserts the result is unit.
template<class Rng>
Vec3 sample_global_direction(
  const ModeWeights& m, const Vec3& Bhat, Rng& rng, bool use_rejection = false)
{
  const Vec3 local = sample_local_direction(m, rng, use_rejection);
  const Vec3 u = to_global(local, build_frame(Bhat));
  assert(std::fabs(norm(u) - 1.0) < 1e-12 && "birth direction must be unit");
  return u;
}

// Analytic normalized marginal pdf in x = cosθ (oracle for the chi-square test):
//   p(x) = I(x) / ∫ I dΩ', with the φ integral folded in, so
//   p(x) = [w_perp(1-x²) + w_par(1/4+3/4x²)] / [w_perp·(4/3) + w_par·1].
inline double pdf_costheta(const ModeWeights& m, double x)
{
  const double I = m.w_perp * (1.0 - x * x) + m.w_par * (0.25 + 0.75 * x * x);
  const double Z = m.w_perp * (4.0 / 3.0) + m.w_par * 1.0;  // ∫ I dx over [-1,1]
  return I / Z;
}

}  // namespace spf

#endif  // SPF_SAMPLER_HPP
