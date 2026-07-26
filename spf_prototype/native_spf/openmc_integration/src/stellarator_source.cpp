//! \file stellarator_source.cpp
//! \brief Implementation of the conformal spin-polarized StellaratorSource.
//!
//! Line-for-line C++ port of the validated Python reference
//! ConformalStellaratorSampler (spf_prototype/python/stellarator_source.py).
//! The CDF cascade, jitter, and trilinear interpolation reproduce that
//! reference's numpy operations exactly (same cumsum axis order, same
//! searchsorted side conventions, same interp term order, same period wraps);
//! the SPF birth-direction sampler is reused verbatim from the native SPF
//! TokamakSource work (spf_prototype/native_spf/source.cpp).

#include "openmc/stellarator_source.h"

#include <algorithm> // for clamp, lower_bound, max
#include <cassert>   // for assert (hard-gated invariants)
#include <cmath>     // for sin, cos, sqrt, floor, acos, cbrt, abs
#include <cstdint>   // for uint64_t
#include <fstream>   // for ifstream (portable, HDF5-free fluxmap reader)
#include <string>
#include <unordered_map>
#include <utility> // for pair, move

#include <fmt/core.h>

#include "openmc/constants.h"        // for PI
#include "openmc/distribution.h"     // for Discrete, distribution_from_xml
#include "openmc/distribution_multi.h" // for Isotropic
#include "openmc/error.h"            // for fatal_error, warning
#include "openmc/particle_type.h"    // for ParticleType
#include "openmc/random_lcg.h"       // for prn
#include "openmc/xml_interface.h"    // for check_for_node, get_node_*

namespace openmc {

namespace {

// --------------------------------------------------------------------------
// SPF per-shape inverse-CDF marginal samplers in x = cos(theta). These are a
// verbatim copy of the native SPF TokamakSource helpers (native_spf/source.cpp,
// themselves ports of spf::sample_costheta_{perp,par}_invcdf in
// spf_prototype/src/spf_sampler.hpp). The roots were derived and
// endpoint-verified there and cross-checked against rejection samplers via
// two-sample KS. No magic constants: every coefficient traces to Schwartz
// Eq. 2's two angular shapes.
// --------------------------------------------------------------------------

// sin^2 theta shape: p(x) ~ (1 - x^2). Inverse CDF solves x^3 - 3x + (4u-2) = 0.
// The root in [-1,1] is the k=1 trig root:
//   x = 2 cos( (1/3) arccos(1 - 2u) - 2*pi/3 ).
inline double spf_costheta_perp(double u)
{
  double x = 2.0 * std::cos(std::acos(1.0 - 2.0 * u) / 3.0 - 2.0 * PI / 3.0);
  return std::clamp(x, -1.0, 1.0); // guard fp drift exactly at the endpoints
}

// 1/4 + 3/4 cos^2 theta shape: p(x) = 1/4 + 3/4 x^2 (normalized on [-1,1]).
// Inverse CDF solves x^3 + x + (2 - 4u) = 0 (single real root, Cardano).
inline double spf_costheta_par(double u)
{
  const double q = 2.0 - 4.0 * u;                        // x^3 + x + q = 0
  const double sq = std::sqrt(q * q / 4.0 + 1.0 / 27.0); // disc > 0 always
  double x = std::cbrt(-q / 2.0 + sq) + std::cbrt(-q / 2.0 - sq);
  return std::clamp(x, -1.0, 1.0);
}

// Periodic-linear interpolation index helper, a line-for-line port of
// _wrap_index() in stellarator_source.py: for a uniform `grid` covering one
// `period`, return the lower/upper bracket cells (i0,i1) and the fractional
// weight `frac` in [0,1). Correct for negative t (jitter just below grid[0]).
inline void wrap_index(double x, const vector<double>& grid, double period,
  int& i0, int& i1, double& frac)
{
  const int n = static_cast<int>(grid.size());
  const double d = period / n;
  const double t = (x - grid[0]) / d;
  const double ft = std::floor(t);
  int idx = static_cast<int>(ft);
  i0 = ((idx % n) + n) % n; // non-negative modulo
  frac = t - ft;            // always in [0,1)
  i1 = (i0 + 1) % n;
}

} // namespace

//==============================================================================
// StellaratorSource implementation
//==============================================================================

int StellaratorSource::upper_bound_strided(
  const double* a, int n, int stride, double u)
{
  // First index m in [0,n] with a[m*stride] > u == std::upper_bound offset ==
  // np.searchsorted(a, u, side="right"). The selected cell is this minus 1.
  int lo = 0;
  int hi = n;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (a[static_cast<std::size_t>(mid) * stride] <= u) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

double StellaratorSource::interp(const vector<double>& A, int ir, int jt0,
  int jt1, int kz0, int kz1, double fr, double ft, double fz) const
{
  // Term order theta (inner) -> zeta -> rho (outer), matching the reference
  // interp() closure exactly. rho bracket is (ir, ir+1) [non-periodic];
  // theta/zeta brackets wrap.
  const double c00 = A[idx3(ir, jt0, kz0)] * (1 - ft) + A[idx3(ir, jt1, kz0)] * ft;
  const double c01 = A[idx3(ir, jt0, kz1)] * (1 - ft) + A[idx3(ir, jt1, kz1)] * ft;
  const double c0 = c00 * (1 - fz) + c01 * fz;
  const double d00 =
    A[idx3(ir + 1, jt0, kz0)] * (1 - ft) + A[idx3(ir + 1, jt1, kz0)] * ft;
  const double d01 =
    A[idx3(ir + 1, jt0, kz1)] * (1 - ft) + A[idx3(ir + 1, jt1, kz1)] * ft;
  const double d0 = d00 * (1 - fz) + d01 * fz;
  return c0 * (1 - fr) + d0 * fr;
}

void StellaratorSource::read_fluxmap(const std::string& stem_in)
{
  // Accept either a bare stem or a path ending in .meta/.bin.
  std::string stem = stem_in;
  auto strip = [&](const std::string& ext) {
    if (stem.size() > ext.size() &&
        stem.compare(stem.size() - ext.size(), ext.size(), ext) == 0) {
      stem.erase(stem.size() - ext.size());
    }
  };
  strip(".meta");
  strip(".bin");

  // ---- parse <stem>.meta (text 'key value' per line) ----
  std::ifstream m(stem + ".meta");
  if (!m) {
    fatal_error(fmt::format(
      "StellaratorSource: cannot open fluxmap header '{}.meta'.", stem));
  }
  std::unordered_map<std::string, std::string> kv;
  std::string k, v;
  while (m >> k >> v) {
    kv[k] = v;
  }
  if (kv["magic"] != "spf_fluxmap_v1") {
    fatal_error(fmt::format(
      "StellaratorSource: bad/missing magic in '{}.meta' (expected "
      "spf_fluxmap_v1).",
      stem));
  }
  nr_ = std::stoi(kv.at("nR"));
  nt_ = std::stoi(kv.at("ntheta"));
  nz_ = std::stoi(kv.at("nzeta"));
  nfp_ = kv.count("nfp") ? std::stoi(kv.at("nfp")) : 1;
  sign_sqrtg_ = kv.count("sign_sqrtg") ? std::stod(kv.at("sign_sqrtg")) : 1.0;
  if (nr_ < 2 || nt_ < 2 || nz_ < 2) {
    fatal_error("StellaratorSource: fluxmap needs nR,ntheta,nzeta >= 2.");
  }

  // ---- read <stem>.bin (little-endian f8, C-order, fixed array order) ----
  std::ifstream b(stem + ".bin", std::ios::binary);
  if (!b) {
    fatal_error(fmt::format(
      "StellaratorSource: cannot open fluxmap payload '{}.bin'.", stem));
  }
  auto read_vec = [&](vector<double>& out, std::size_t n) {
    out.resize(n);
    b.read(reinterpret_cast<char*>(out.data()), n * sizeof(double));
  };
  const std::size_t n3 =
    static_cast<std::size_t>(nr_) * static_cast<std::size_t>(nt_) *
    static_cast<std::size_t>(nz_);
  read_vec(rho_, nr_);
  read_vec(theta_, nt_);
  read_vec(zeta_, nz_);
  read_vec(sqrtg_, n3);
  read_vec(R_, n3);
  read_vec(Z_, n3);
  read_vec(phi_, n3);
  read_vec(BR_, n3);
  read_vec(Bphi_, n3);
  read_vec(BZ_, n3);
  if (!b) {
    fatal_error(fmt::format(
      "StellaratorSource: short read on fluxmap payload '{}.bin'.", stem));
  }

  // Uniform spacings (the producer emits uniform grids).
  drho_ = rho_[1] - rho_[0];
  dtheta_ = theta_[1] - theta_[0];
  dzeta_ = zeta_[1] - zeta_[0];
  if (!(drho_ > 0.0) || !(dtheta_ > 0.0) || !(dzeta_ > 0.0)) {
    fatal_error("StellaratorSource: fluxmap grids must be strictly increasing.");
  }

  // Invariant (not user input): sqrt(g) is stored as a non-negative magnitude.
  double sg_min = sqrtg_[0];
  for (std::size_t n = 1; n < sqrtg_.size(); ++n) {
    sg_min = std::min(sg_min, sqrtg_[n]);
  }
  assert(sg_min >= -1e-12 && "StellaratorSource: sqrt(g) must be non-negative");
  (void)sg_min;
}

void StellaratorSource::build_cdf_cascade()
{
  // Per-cell emission weight w[i][j][k] = S(rho_i) * |sqrt(g)|. Uniform
  // (rho,theta,zeta) spacing means the constant cell volume cancels in every
  // normalized CDF, so it is omitted (exactly as in the reference).
  const std::size_t n3 = sqrtg_.size();
  vector<double> w(n3);
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0; j < nt_; ++j) {
      for (int kk = 0; kk < nz_; ++kk) {
        w[idx3(i, j, kk)] = S_[i] * sqrtg_[idx3(i, j, kk)];
      }
    }
  }

  // ---- marginal p(rho): cdf_rho_ [nr+1], cumulative then normalized ----
  cdf_rho_.assign(nr_ + 1, 0.0);
  for (int i = 0; i < nr_; ++i) {
    double m = 0.0;
    for (int j = 0; j < nt_; ++j) {
      for (int kk = 0; kk < nz_; ++kk) {
        m += w[idx3(i, j, kk)];
      }
    }
    cdf_rho_[i + 1] = cdf_rho_[i] + m;
  }
  const double tot = cdf_rho_[nr_];
  if (!(tot > 0.0)) {
    fatal_error("StellaratorSource: integrated emission S(rho)*sqrt(g) is zero "
                "or negative. Check the fluxmap / emission profile.");
  }
  const double inv_tot = 1.0 / tot;
  for (int i = 0; i <= nr_; ++i) {
    cdf_rho_[i] *= inv_tot;
  }

  // ---- conditional p(zeta|rho): cdf_zeta_ [nr][nz+1], per-rho normalized ----
  cdf_zeta_.assign(static_cast<std::size_t>(nr_) * (nz_ + 1), 0.0);
  for (int i = 0; i < nr_; ++i) {
    for (int kk = 0; kk < nz_; ++kk) {
      double s = 0.0;
      for (int j = 0; j < nt_; ++j) {
        s += w[idx3(i, j, kk)]; // cz[i][k] = sum_theta w
      }
      cdf_zeta_[idxZ(i, kk + 1)] = cdf_zeta_[idxZ(i, kk)] + s;
    }
    const double dz = cdf_zeta_[idxZ(i, nz_)];
    const double inv = (dz > 0.0) ? 1.0 / dz : 1.0; // where(dz>0, dz, 1.0)
    for (int kk = 0; kk <= nz_; ++kk) {
      cdf_zeta_[idxZ(i, kk)] *= inv;
    }
  }

  // ---- conditional p(theta|rho,zeta): cdf_theta_ [nr][nt+1][nz] ----
  cdf_theta_.assign(
    static_cast<std::size_t>(nr_) * (nt_ + 1) * nz_, 0.0);
  for (int i = 0; i < nr_; ++i) {
    for (int kk = 0; kk < nz_; ++kk) {
      for (int j = 0; j < nt_; ++j) { // cumsum over theta
        cdf_theta_[idxT(i, j + 1, kk)] =
          cdf_theta_[idxT(i, j, kk)] + w[idx3(i, j, kk)];
      }
      const double dt = cdf_theta_[idxT(i, nt_, kk)];
      const double inv = (dt > 0.0) ? 1.0 / dt : 1.0; // where(dt>0, dt, 1.0)
      for (int j = 0; j <= nt_; ++j) {
        cdf_theta_[idxT(i, j, kk)] *= inv;
      }
    }
  }

  // Invariants (not user input): the marginal CDF is monotone and ends at 1.
#ifndef NDEBUG
  for (int i = 0; i < nr_; ++i) {
    assert(cdf_rho_[i] <= cdf_rho_[i + 1] + 1e-12 &&
           "StellaratorSource: cdf_rho must be monotone non-decreasing");
  }
  assert(std::abs(cdf_rho_[nr_] - 1.0) < 1e-9 &&
         "StellaratorSource: cdf_rho must end at 1");
#endif
}

StellaratorSource::StellaratorSource(pugi::xml_node node) : Source(node)
{
  // ---- fluxmap path (required) ----
  if (!check_for_node(node, "fluxmap")) {
    fatal_error("StellaratorSource: a 'fluxmap' path is required.");
  }
  read_fluxmap(get_node_value(node, "fluxmap", false, true));

  // ---- radial emission weight S(rho) ----
  // Optional (emission_rho, emission_density) override, else parabolic 1-rho^2
  // (the reference default). Only the shape matters; normalization cancels.
  S_.resize(nr_);
  if (check_for_node(node, "emission_rho") &&
      check_for_node(node, "emission_density")) {
    auto er = get_node_array<double>(node, "emission_rho");
    auto ed = get_node_array<double>(node, "emission_density");
    if (er.size() != ed.size() || er.size() < 2) {
      fatal_error("StellaratorSource: emission_rho and emission_density must "
                  "have equal length >= 2.");
    }
    for (std::size_t p = 1; p < er.size(); ++p) {
      if (er[p] <= er[p - 1]) {
        fatal_error("StellaratorSource: emission_rho must be increasing.");
      }
    }
    for (double s : ed) {
      if (s < 0.0) {
        fatal_error("StellaratorSource: emission_density must be >= 0.");
      }
    }
    for (int i = 0; i < nr_; ++i) {
      const double r = rho_[i];
      double s;
      if (r <= er.front()) {
        s = ed.front();
      } else if (r >= er.back()) {
        s = ed.back();
      } else {
        auto it = std::lower_bound(er.begin(), er.end(), r);
        const std::size_t p = static_cast<std::size_t>(it - er.begin());
        const double t = (r - er[p - 1]) / (er[p] - er[p - 1]);
        s = ed[p - 1] + t * (ed[p] - ed[p - 1]);
      }
      S_[i] = std::max(0.0, s);
    }
  } else {
    for (int i = 0; i < nr_; ++i) {
      const double r = rho_[i];
      S_[i] = std::max(0.0, 1.0 - r * r); // parabolic flux function
    }
  }

  build_cdf_cascade();

  // ---- energy distribution(s): single distribution used for all rho ----
  for (auto energy_node : node.children("energy")) {
    energy_dists_.push_back(distribution_from_xml(energy_node));
  }
  if (energy_dists_.empty()) {
    // Default monoenergetic 14.06 MeV (DT birth energy; energy is decoupled
    // from the polarized angular distribution to leading order).
    double E[] {14.06e6};
    double p[] {1.0};
    energy_dists_.push_back(UPtrDist {new Discrete {E, p, 1}});
  }

  // ---- time distribution: default delta at t=0 ----
  if (check_for_node(node, "time")) {
    time_ = distribution_from_xml(node.child("time"));
  } else {
    double T[] {0.0};
    double p[] {1.0};
    time_ = UPtrDist {new Discrete {T, p, 1}};
  }

  // ---- isotropic fallback (unpolarized path) ----
  angle_ = UPtrAngle {new Isotropic()};

  // ---- SPF polarization (optional; absent => isotropic, unchanged behavior) --
  // Physics (Schwartz 2025, Eq. 2), constants DERIVED not pasted:
  //   dsigma/dOmega = (s0/2pi)[ (3/4)a sin^2 t + ((2/3)b+(1/3)c)(1/4+3/4cos^2 t) ]
  //   w_perp = (3/4)a ; w_par = (2/3)b + (1/3)c ; eta = a + (2/3)b + (1/3)c.
  //   int sin^2 dOmega = 8pi/3, int (1/4+3/4cos^2) dOmega = 2pi
  //   => W_perp = 2pi a, W_par = 2pi w_par => P_perp = a/eta. Same as A's derivation.
  if (check_for_node(node, "polarization")) {
    auto abc = get_node_array<double>(node, "polarization"); // [a, b, c]
    if (abc.size() != 3) {
      fatal_error("StellaratorSource: polarization must be 3 numbers 'a b c'.");
    }
    double a = abc[0], b = abc[1], c = abc[2];
    // User input -> fatal_error on invalid; NaN fails the >= 0 test.
    if (!(a >= 0.0) || !(b >= 0.0) || !(c >= 0.0)) {
      fatal_error("StellaratorSource: polarization fractions a,b,c must be >= 0.");
    }
    const double s = a + b + c;
    if (!(s > 0.0)) {
      fatal_error("StellaratorSource: polarization a+b+c must be > 0.");
    }
    a /= s; // renormalize to sum 1 (Python side already warns on deviation)
    b /= s;
    c /= s;

    w_perp_ = 0.75 * a;                         // (3/4) a
    w_par_ = (2.0 / 3.0) * b + (1.0 / 3.0) * c; // (2/3)b + (1/3)c
    const double eta = a + (2.0 / 3.0) * b + (1.0 / 3.0) * c;
    p_perp_ = a / eta; // == 2pi a / (2pi a + 2pi w_par)
    // Invariant (not user input): P_perp is a valid probability. Hard-assert.
    assert(p_perp_ >= -1e-15 && p_perp_ <= 1.0 + 1e-15 &&
           "StellaratorSource SPF: P_perp must lie in [0,1]");
    polarized_ = true;

    // Field model: default reads b-hat from the fluxmap grid. "toroidal" is an
    // override that uses the pure phi-hat field (axisymmetric-limit gate vs the
    // TokamakSource with field_model='toroidal').
    field_model_ = 0;
    if (check_for_node(node, "field_model")) {
      std::string fm = get_node_value(node, "field_model");
      if (fm == "fluxmap") {
        field_model_ = 0;
      } else if (fm == "toroidal") {
        field_model_ = 1;
      } else {
        fatal_error(
          "StellaratorSource: field_model must be 'fluxmap' or 'toroidal'.");
      }
    }
  }
}

Direction StellaratorSource::sample_polarized_direction(
  Direction bhat, uint64_t* seed) const
{
  // 1. cos(theta) relative to b-hat via stratified selection of the two P2
  //    shapes. RNG draw order: selector -> shape-u -> phi (matches
  //    spf_sampler.hpp so the native port is bit-for-bit identical on a fixed
  //    seed sequence).
  double x; // cos(theta) relative to b-hat
  if (prn(seed) < p_perp_) {
    x = spf_costheta_perp(prn(seed));
  } else {
    x = spf_costheta_par(prn(seed));
  }
  const double az = 2.0 * PI * prn(seed);
  const double st = std::sqrt(std::max(0.0, 1.0 - x * x));
  const Direction local {st * std::cos(az), st * std::sin(az), x};

  // 2. Rotate from the b-hat-aligned local frame to global via Gram-Schmidt
  //    with a least-aligned reference axis (avoids the b-hat || z_world
  //    degeneracy).
  Direction b = bhat / bhat.norm();
  const double ax = std::abs(b.x), ay = std::abs(b.y), az2 = std::abs(b.z);
  Direction ref = (ax <= ay && ax <= az2) ? Direction {1.0, 0.0, 0.0}
                  : (ay <= az2)           ? Direction {0.0, 1.0, 0.0}
                                          : Direction {0.0, 0.0, 1.0};
  Direction ex = ref - b.dot(ref) * b;
  ex /= ex.norm();
  Direction ey = b.cross(ex);

  Direction u = local.x * ex + local.y * ey + local.z * b;
  return u / u.norm(); // guard fp drift; |u| ~ 1
}

SourceSite StellaratorSource::sample(uint64_t* seed) const
{
  SourceSite site;
  site.particle = ParticleType::neutron();
  site.delayed_group = 0;

  // ---- 1. cascade cell selection (rho -> zeta -> theta), side="right" ----
  const double u_rho = prn(seed);
  int i = upper_bound_strided(cdf_rho_.data(), nr_ + 1, 1, u_rho) - 1;
  i = std::clamp(i, 0, nr_ - 1);

  const double u_zeta = prn(seed);
  int kk = upper_bound_strided(&cdf_zeta_[idxZ(i, 0)], nz_ + 1, 1, u_zeta) - 1;
  kk = std::clamp(kk, 0, nz_ - 1);

  const double u_theta = prn(seed);
  int j =
    upper_bound_strided(&cdf_theta_[idxT(i, 0, kk)], nt_ + 1, nz_, u_theta) - 1;
  j = std::clamp(j, 0, nt_ - 1);

  // ---- 2. jitter uniformly inside the selected cell (order rho,theta,zeta) --
  const double rr = rho_[i] + (prn(seed) - 0.5) * drho_;
  const double tt = theta_[j] + (prn(seed) - 0.5) * dtheta_;
  const double zz = zeta_[kk] + (prn(seed) - 0.5) * dzeta_;

  // ---- 3. trilinear-interp brackets: rho clamped, theta/zeta 2*pi-periodic --
  int ir = static_cast<int>(
             std::lower_bound(rho_.begin(), rho_.end(), rr) - rho_.begin()) -
           1;
  ir = std::clamp(ir, 0, nr_ - 2);
  const double fr = (rr - rho_[ir]) / drho_;
  int jt0, jt1;
  double ft;
  wrap_index(tt, theta_, 2.0 * PI, jt0, jt1, ft);
  int kz0, kz1;
  double fz;
  wrap_index(zz, zeta_, 2.0 * PI, kz0, kz1, fz);

  // ---- 4. interp geometry -> Cartesian position ----
  const double R = interp(R_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
  const double Z = interp(Z_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
  const double PH = interp(phi_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
  const double cph = std::cos(PH), sph = std::sin(PH);
  site.r = {R * cph, R * sph, Z};

  // ---- 5. local field b-hat ----
  Direction bhat;
  if (field_model_ == 1) {
    // Pure toroidal phi-hat = (-sin PH, cos PH, 0) (already unit).
    bhat = {-sph, cph, 0.0};
  } else {
    // b-hat from the fluxmap: interp (BR,Bphi,BZ), rotate rpz->Cartesian, norm.
    const double br = interp(BR_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
    const double bp = interp(Bphi_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
    const double bz = interp(BZ_, ir, jt0, jt1, kz0, kz1, fr, ft, fz);
    const double bx = br * cph - bp * sph;
    const double by = br * sph + bp * cph;
    double bn = std::sqrt(bx * bx + by * by + bz * bz);
    if (!(bn > 0.0)) {
      bn = 1.0; // degenerate field guard (matches reference +1e-300)
    }
    bhat = {bx / bn, by / bn, bz / bn};
  }

  // ---- 6. birth direction ----
  if (polarized_) {
    site.u = sample_polarized_direction(bhat, seed);
  } else {
    site.u = angle_->sample(seed).first;
  }

  // ---- 7. energy + time ----
  auto [E, E_wgt] = energy_dists_[0]->sample(seed);
  site.E = E;
  auto [t, t_wgt] = time_->sample(seed);
  site.time = t;

  // The rejection-free cascade gives an exactly-1 spatial birth weight; the
  // energy/time importance weights are 1 for the default monoenergetic/delta
  // distributions and are carried through for biased distributions.
  site.wgt = E_wgt * t_wgt;

  return site;
}

} // namespace openmc
