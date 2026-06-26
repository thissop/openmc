// spf_fieldmap.hpp -- FieldMapField: a general B-hat(x) read from a gridded field
// map (the 3D-tier field). Reads the portable raw format written by
// python/make_standin_field.py (<stem>.meta + <stem>.bin) and interpolates with a
// hand-coded multilinear scheme that is BIT-IDENTICAL to python/spf_mirror.py's
// FieldMapField (clamped R/Z, PERIODIC phi, fixed term order R->phi->Z, then
// renormalize). OpenMC-free (std::ifstream only) so the field-parity driver
// builds with plain g++ and the C++<->Python parity test covers the field map.
//
// INJECT(helios): a real DESC/VMEC QA equilibrium writes the same format; this
// reader is producer-agnostic.
// DEFERRED: raw little-endian f8 assumes same-endian platforms. HDF5 (already
// linked by libopenmc) is the documented portable upgrade (see DEFERRED.md).
#ifndef SPF_FIELDMAP_HPP
#define SPF_FIELDMAP_HPP

#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "spf_field.hpp"

namespace spf {

class FieldMapField : public MagneticField {
public:
  explicit FieldMapField(const std::string& stem)
  {
    // ---- parse <stem>.meta (text key value per line) ----
    std::ifstream m(stem + ".meta");
    if (!m)
      throw std::runtime_error("FieldMapField: cannot open " + stem + ".meta");
    std::unordered_map<std::string, std::string> kv;
    std::string k, v;
    while (m >> k >> v)
      kv[k] = v;
    if (kv["magic"] != "spf_fieldmap_v1")
      throw std::runtime_error("FieldMapField: bad/missing magic in " + stem);
    nR_ = std::stoi(kv.at("nR"));
    nphi_ = std::stoi(kv.at("nphi"));
    nZ_ = std::stoi(kv.at("nZ"));
    Rmin_ = std::stod(kv.at("R_min"));
    Zmin_ = std::stod(kv.at("Z_min"));
    phimin_ = std::stod(kv.at("phi_min"));
    period_ = std::stod(kv.at("period"));
    dR_ = (std::stod(kv.at("R_max")) - Rmin_) / (nR_ - 1);
    dZ_ = (std::stod(kv.at("Z_max")) - Zmin_) / (nZ_ - 1);
    dphi_ = period_ / nphi_;

    // ---- read <stem>.bin : little-endian f8, Bx then By then Bz, C-order ----
    const std::size_t n = std::size_t(nR_) * nphi_ * nZ_;
    std::ifstream b(stem + ".bin", std::ios::binary);
    if (!b)
      throw std::runtime_error("FieldMapField: cannot open " + stem + ".bin");
    bx_.resize(n);
    by_.resize(n);
    bz_.resize(n);
    b.read(reinterpret_cast<char*>(bx_.data()), n * sizeof(double));
    b.read(reinterpret_cast<char*>(by_.data()), n * sizeof(double));
    b.read(reinterpret_cast<char*>(bz_.data()), n * sizeof(double));
    if (!b)
      throw std::runtime_error("FieldMapField: short read on " + stem + ".bin");
  }

  Vec3 bhat(const Vec3& r) const override
  {
    const double R = std::sqrt(r.x * r.x + r.y * r.y);
    const double phi = std::atan2(r.y, r.x);
    const double Z = r.z;

    // R, Z: clamped uniform index + fractional weight
    const double tR = (R - Rmin_) / dR_;
    const int iR = clamp_idx(int(std::floor(tR)), nR_);
    const double fR = tR - iR;
    const double tZ = (Z - Zmin_) / dZ_;
    const int iZ = clamp_idx(int(std::floor(tZ)), nZ_);
    const double fZ = tZ - iZ;
    // phi: periodic wrap into [phimin, phimin+period)
    double pp = std::fmod(phi - phimin_, period_);
    if (pp < 0.0)
      pp += period_;
    const double tphi = pp / dphi_;
    int ip = int(std::floor(tphi));
    const double fphi = tphi - ip;
    ip = ip % nphi_;
    if (ip < 0)
      ip += nphi_;
    const int ip1 = (ip + 1) % nphi_;
    const int iR1 = iR + 1, iZ1 = iZ + 1;

    Vec3 out;
    out.x = tri(bx_, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ);
    out.y = tri(by_, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ);
    out.z = tri(bz_, iR, iR1, ip, ip1, iZ, iZ1, fR, fphi, fZ);
    return normalized(out); // interp of unit vectors is sub-unit -> renormalize
  }

private:
  static int clamp_idx(int i, int n)
  {
    if (i < 0)
      return 0;
    if (i > n - 2)
      return n - 2;
    return i;
  }
  std::size_t idx(int iR, int ip, int iZ) const
  {
    return (std::size_t(iR) * nphi_ + ip) * nZ_ + iZ;
  }
  // trilinear with fixed term order: R outer, phi mid, Z inner.
  double tri(const std::vector<double>& a, int iR, int iR1, int ip, int ip1,
    int iZ, int iZ1, double fR, double fp, double fZ) const
  {
    const double c00 = a[idx(iR, ip, iZ)] * (1 - fZ) + a[idx(iR, ip, iZ1)] * fZ;
    const double c01 = a[idx(iR, ip1, iZ)] * (1 - fZ) + a[idx(iR, ip1, iZ1)] * fZ;
    const double c10 = a[idx(iR1, ip, iZ)] * (1 - fZ) + a[idx(iR1, ip, iZ1)] * fZ;
    const double c11 = a[idx(iR1, ip1, iZ)] * (1 - fZ) + a[idx(iR1, ip1, iZ1)] * fZ;
    const double c0 = c00 * (1 - fp) + c01 * fp;
    const double c1 = c10 * (1 - fp) + c11 * fp;
    return c0 * (1 - fR) + c1 * fR;
  }

  int nR_, nphi_, nZ_;
  double Rmin_, Zmin_, phimin_, period_, dR_, dZ_, dphi_;
  std::vector<double> bx_, by_, bz_;
};

} // namespace spf

#endif // SPF_FIELDMAP_HPP
