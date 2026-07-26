//! \file stellarator_source.h
//! \brief Conformal spin-polarized stellarator plasma neutron source
//!
//! Native OpenMC-core Source subclass, patterned on Ethan Peterson's
//! TokamakSource (PR #3999). It samples birth POSITIONS from a REAL-equilibrium
//! flux map (the spf_fluxmap_v1 .meta/.bin written by
//! spf_prototype/python/quasr_fluxmap.py) via a rejection-free
//! marginal->conditional CDF cascade, and birth DIRECTIONS from the Schwartz P2
//! spin-polarized-fusion mixture about the local field b-hat read from the SAME
//! grid. It is a line-for-line C++ port of the validated Python reference
//! ConformalStellaratorSampler in spf_prototype/python/stellarator_source.py.
//!
//! To register in an OpenMC build:
//!   * add this header under include/openmc/,
//!   * add stellarator_source.cpp under src/ and to CMakeLists.txt next to
//!     src/source.cpp,
//!   * add a "stellarator" branch to Source::create (src/source.cpp) and to
//!     SourceBase.from_xml_element (openmc/source.py).
//! See native_spf/README_stellarator.md for the exact snippets.

#ifndef OPENMC_STELLARATOR_SOURCE_H
#define OPENMC_STELLARATOR_SOURCE_H

#include <cstddef> // for size_t
#include <string>

#include "pugixml.hpp"

#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/memory.h" // for unique_ptr, UPtrDist, UPtrAngle
#include "openmc/position.h"
#include "openmc/source.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Conformal spin-polarized stellarator plasma neutron source
//!
//! Birth positions are drawn from  p(rho,theta,zeta) ~ S(rho) * sqrt(g) via a
//! REJECTION-FREE marginal->conditional CDF cascade
//!   p(rho) * p(zeta|rho) * p(theta|rho,zeta)
//! so every birth weight is exactly 1. Geometry (R,Z,phi), the metric Jacobian
//! sqrt(g), and the field (BR,Bphi,BZ) all come from ONE (rho,theta,zeta) grid
//! -- the spf_fluxmap_v1 file -- so position and field are self-consistent on a
//! single equilibrium.
//!
//! Sampling algorithm (matches ConformalStellaratorSampler line-for-line):
//!  1. Cascade cell selection: searchsorted(side="right") on cdf_rho, then on
//!     cdf_zeta[i], then on cdf_theta[i,:,k] -> discrete cell (i,j,k).
//!  2. Jitter uniformly inside the selected cell in (rho,theta,zeta).
//!  3. Trilinear interp of R,Z,phi (rho non-periodic; theta,zeta 2*pi-periodic)
//!     -> cylindrical -> Cartesian position.
//!  4. b-hat from trilinear interp of BR,Bphi,BZ -> Cartesian -> normalize
//!     (or the pure-toroidal phi-hat override).
//!  5. Direction: isotropic if unpolarized, else the Schwartz P2 mode-mixture
//!     about b-hat (reused verbatim from the native SPF TokamakSource work).
//!
//! All sampler state is set in the constructor and const at sample() time: no
//! mutable members, thread-safe by construction, openmc::prn only.
//==============================================================================

class StellaratorSource : public Source {
public:
  // Constructors
  explicit StellaratorSource(pugi::xml_node node);

  //! Sample from the stellarator source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site (weight exactly 1)
  SourceSite sample(uint64_t* seed) const override;

private:
  //==========================================================================
  // Construction helpers

  //! Read the spf_fluxmap_v1 <stem>.meta + <stem>.bin into the grid members.
  void read_fluxmap(const std::string& stem);

  //! Build the marginal/conditional CDF cascade from S(rho)*sqrt(g). Mirrors
  //! ConformalStellaratorSampler.__init__ exactly (same cumsum axis order and
  //! per-slice normalization), so the C++ cascade selects the same cells as the
  //! validated Python reference given the same uniform draws.
  void build_cdf_cascade();

  //==========================================================================
  // Sampling helpers

  //! Flat C-order index into an (nr,nt,nz) grid array: [i_rho][j_theta][k_zeta].
  std::size_t idx3(int i, int j, int k) const
  {
    return (static_cast<std::size_t>(i) * nt_ + j) * nz_ + k;
  }
  //! Flat index into cdf_zeta_ laid out [nr][nz+1].
  std::size_t idxZ(int i, int k) const
  {
    return static_cast<std::size_t>(i) * (nz_ + 1) + k;
  }
  //! Flat index into cdf_theta_ laid out [nr][nt+1][nz] (matches the Python
  //! numpy layout _cdf_theta[i, j, k]).
  std::size_t idxT(int i, int j, int k) const
  {
    return (static_cast<std::size_t>(i) * (nt_ + 1) + j) * nz_ + k;
  }

  //! searchsorted(a, u, side="right") over a strided ascending array, i.e. the
  //! count of elements a[m*stride] <= u (== std::upper_bound index). The
  //! selected cell is this minus 1. Exactly reproduces numpy's convention used
  //! by ConformalStellaratorSampler._sample_cells.
  static int upper_bound_strided(
    const double* a, int n, int stride, double u);

  //! Trilinear interpolation of a grid array A at the bracketed cell. Term order
  //! is theta (inner) -> zeta -> rho (outer), identical to the reference
  //! interp() closure. rho is clamped (non-periodic); theta/zeta wrap.
  double interp(const vector<double>& A, int ir, int jt0, int jt1, int kz0,
    int kz1, double fr, double ft, double fz) const;

  //! Sample a birth direction from the Schwartz P2 mode-mixture about b-hat.
  //! Reused verbatim from the native SPF TokamakSource (native_spf/source.cpp):
  //! stratified P2 shape selection (RNG order selector -> shape-u -> phi) then a
  //! Gram-Schmidt rotation with a least-aligned reference axis.
  //! \param bhat Local field direction (renormalized here)
  //! \param seed PRN seed pointer
  //! \return Unit neutron direction in global Cartesian
  Direction sample_polarized_direction(Direction bhat, uint64_t* seed) const;

  //==========================================================================
  // Data members (all set in the constructor; const at sample time)

  // Fluxmap grid dimensions and provenance
  int nr_ {0};              //!< number of rho grid points
  int nt_ {0};              //!< number of theta grid points (full torus)
  int nz_ {0};              //!< number of zeta grid points (full torus)
  int nfp_ {1};             //!< number of field periods (informational)
  double sign_sqrtg_ {1.0}; //!< sign of raw sqrt(g) before |.| (provenance)

  // 1-D flux-coordinate grids
  vector<double> rho_;   //!< rho grid  [nr]  (drho..1]
  vector<double> theta_; //!< theta grid [nt] [0,2*pi)
  vector<double> zeta_;  //!< zeta grid  [nz] [0,2*pi)
  double drho_ {0.0};    //!< uniform rho spacing
  double dtheta_ {0.0};  //!< uniform theta spacing (2*pi/nt)
  double dzeta_ {0.0};   //!< uniform zeta spacing (2*pi/nz)

  // 3-D geometry + field, flat C-order [i_rho][j_theta][k_zeta]
  vector<double> sqrtg_; //!< |sqrt(g)| metric Jacobian (>= 0)
  vector<double> R_;     //!< cylindrical major radius R
  vector<double> Z_;     //!< cylindrical height Z
  vector<double> phi_;   //!< cylindrical toroidal angle phi
  vector<double> BR_;    //!< field component B_R
  vector<double> Bphi_;  //!< field component B_phi
  vector<double> BZ_;    //!< field component B_Z

  // Radial emission weight S(rho) on the rho grid (default parabolic 1-rho^2)
  vector<double> S_;

  // The CDF cascade (rejection-free marginal -> conditional sampling)
  vector<double> cdf_rho_;   //!< marginal p(rho)   [nr+1], 0..1
  vector<double> cdf_zeta_;  //!< p(zeta|rho)       [nr][nz+1], per-i 0..1
  vector<double> cdf_theta_; //!< p(theta|rho,zeta) [nr][nt+1][nz], per-(i,k) 0..1

  // Energy / time / fallback-angle distributions
  vector<unique_ptr<Distribution>> energy_dists_; //!< energy (single for all r)
  UPtrDist time_;  //!< time distribution (default delta at t=0)
  UPtrAngle angle_; //!< isotropic fallback (unpolarized path)

  // --- SPF (spin-polarized fusion) emission ------------------------------
  // Derived from (a,b,c) in the constructor from Schwartz Eq. 2 coefficients;
  // no magic constants stored. Same discipline as the native SPF TokamakSource.
  bool polarized_ {false}; //!< false => isotropic (backward compatible)
  double w_perp_ {0.0};    //!< (3/4)*a          -- sin^2 theta weight
  double w_par_ {0.0};     //!< (2/3)*b + (1/3)*c -- (1/4+3/4cos^2 theta) weight
  double p_perp_ {0.0};    //!< a/eta -- stratified selection prob of sin^2 shape
  int field_model_ {0};    //!< 0 = fluxmap b-hat (grid B), 1 = pure toroidal phi-hat
};

//==============================================================================
// Functions
//==============================================================================

//! Map deuteron/triton spin-projection fractions to the (a,b,c) collision-mode
//! fractions (Schwartz 2025, Eq. 1). Provided for parity with the native SPF
//! TokamakSource; the Python StellaratorSource also exposes it.
//!   a = d_plus*t_plus + d_minus*t_minus
//!   b = d_zero
//!   c = d_plus*t_minus + d_minus*t_plus

} // namespace openmc

#endif // OPENMC_STELLARATOR_SOURCE_H
