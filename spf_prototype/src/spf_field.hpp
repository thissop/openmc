// spf_field.hpp -- Magnetic-field DIRECTION models B-hat(x) that feed the verified
// spin-polarized birth-direction sampler. Header-only and OpenMC-free (depends
// only on spf_sampler.hpp's Vec3 + <cmath>).
//
// SCOPE: this file ONLY supplies the unit field direction B-hat at a position.
// The sampler (spf_sampler.hpp) rotates emission about whatever B-hat it is given
// (Gram-Schmidt frame, tested off-axis and degenerate); that math is NOT touched.
// Every field here is mirrored 1:1 in python/spf_mirror.py and checked by
// tests/test_field_parity.py (C++ vs Python bhat to <= 1e-12).
//
// Conventions (lengths in cm; global Cartesian). At a position r=(x,y,z) define
// phi_p = atan2(y,x) and the local cylindrical basis
//     phi_hat = (-sin phi_p, cos phi_p, 0),  R_hat = (cos phi_p, sin phi_p, 0),  z_hat = (0,0,1).
//
//   ToroidalField:  B-hat = phi_hat                              (the beta=0 limit)
//   AngledField:    B-hat = cos(beta)*phi_hat + sin(beta)*( cos(alpha)*R_hat - sin(alpha)*z_hat )
//       beta = polar tilt of B-hat away from the toroidal direction phi_hat
//              (beta=0  => purely toroidal; beta=pi/2 => B-hat in the poloidal plane),
//       alpha = azimuth of that tilt in the poloidal (R_hat, z_hat) plane.
//   This (alpha,beta) convention is pinned NUMERICALLY to anarrima's angled-field
//   kernels (g_HAa/g_Hca/g_HBa, g_VAa/g_Vca/g_VBa) to ~2e-14: beta=0 reduces to the
//   toroidal kernel exactly, and the generalized quad cos(theta)=Delta_hat . B-hat
//   reproduces the angled kernels across walls/modes/pitches. See
//   RESULTS_tier8_angled.md and python/analytic_nwl.py (angled engine).
//
// INJECT(helios): a general 3D field B-hat(x) from a real equilibrium plugs in as
// another MagneticField subclass (FieldMapField, added for the 3D tier) without
// changing the sampler or the plugin call site.
#ifndef SPF_FIELD_HPP
#define SPF_FIELD_HPP

#include <cmath>

#include "spf_sampler.hpp"

namespace spf {

// Abstract field-direction model. bhat() is const + thread-safe (no mutable state).
struct MagneticField {
  virtual ~MagneticField() = default;
  // Field direction at global position r (cm). Need not be exactly unit: the
  // sampler's build_frame() normalizes; field-map subclasses renormalize.
  virtual Vec3 bhat(const Vec3& r) const = 0;
};

// Fixed direction everywhere (normalized once). Matches the old bmode=constant.
class ConstantField : public MagneticField {
public:
  explicit ConstantField(const Vec3& b0) : b_(normalized(b0)) {}
  Vec3 bhat(const Vec3&) const override { return b_; }

private:
  Vec3 b_;
};

// Purely toroidal B-hat = phi_hat. Bit-identical to the original plugin formula
// {-y/rxy, x/rxy, 0} (the toroidal regression anchor).
class ToroidalField : public MagneticField {
public:
  Vec3 bhat(const Vec3& r) const override
  {
    const double rxy = std::sqrt(r.x * r.x + r.y * r.y);
    return {-r.y / rxy, r.x / rxy, 0.0}; // phi_hat = (-sin phi, cos phi, 0)
  }
};

// Axisymmetric pitched field: a toroidal component plus a poloidal tilt.
//   B-hat = cos(beta)*phi_hat + sin(beta)*( cos(alpha)*R_hat - sin(alpha)*z_hat ).
// Unit by construction (|B-hat|^2 = cos^2 b + sin^2 b = 1); beta=0 => toroidal.
class AngledField : public MagneticField {
public:
  AngledField(double alpha, double beta)
    : ca_(std::cos(alpha)), sa_(std::sin(alpha)), cb_(std::cos(beta)),
      sb_(std::sin(beta))
  {}
  Vec3 bhat(const Vec3& r) const override
  {
    const double rxy = std::sqrt(r.x * r.x + r.y * r.y);
    const double cphi = r.x / rxy, sphi = r.y / rxy;
    // cb*phi_hat + sb*ca*R_hat - sb*sa*z_hat
    return {cb_ * (-sphi) + sb_ * ca_ * cphi,
            cb_ * (cphi) + sb_ * ca_ * sphi,
            -sb_ * sa_};
  }

private:
  double ca_, sa_, cb_, sb_;
};

} // namespace spf

#endif // SPF_FIELD_HPP
