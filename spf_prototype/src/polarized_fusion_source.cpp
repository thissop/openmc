// polarized_fusion_source.cpp -- OpenMC CompiledSource that emits spin-polarized
// DT neutrons with the Schwartz (2025) birth-direction distribution. It is a
// thin openmc::Source wrapper around the tier-1-verified pure-math header
// spf_sampler.hpp: position + B̂ here, the angular physics there.
//
// Thread-safety: all state is set in the constructor and is const thereafter;
// sample() uses only openmc::prn(seed). No mutable members.
//
// Parameter string (comma-separated key=value, lengths in cm):
//   a,b,c          polarization mode fractions (renormalized; warn if |sum-1|>1e-6)
//   bmode          field model: "toroidal" (B̂=φ̂), "constant", or "angled" (pitched)
//   bx,by,bz       field direction for bmode=constant (default +z)
//   alpha,beta     pitch angles (rad) for bmode=angled (see spf_field.hpp; beta=0=toroidal)
//   shape          "ring" (filamentary) or "plasma" (parabolic circular)
//   r0,z0          ring major radius & height            (shape=ring)
//   R0,aminor      plasma major & minor radius           (shape=plasma)
//   energy         birth energy in eV (default 14.1e6, monoenergetic)
//                  INJECT(helios): Ballabio/D-D broadening hook -- NOT implemented.
// Example: "a=1,b=0,c=0,bmode=angled,alpha=0.0,beta=0.4,shape=plasma,R0=100,aminor=50"
#include <cmath>
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/source.h"

#include "spf_field.hpp"
#include "spf_sampler.hpp"

namespace {

std::unordered_map<std::string, std::string> parse_params(const std::string& s)
{
  std::unordered_map<std::string, std::string> kv;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    auto eq = tok.find('=');
    if (eq == std::string::npos)
      continue;
    std::string key = tok.substr(0, tok.find_first_of('=', 0));
    // strip surrounding whitespace from key/value
    auto strip = [](std::string x) {
      const char* ws = " \t\n\r";
      auto a = x.find_first_not_of(ws);
      if (a == std::string::npos)
        return std::string {};
      auto b = x.find_last_not_of(ws);
      return x.substr(a, b - a + 1);
    };
    kv[strip(tok.substr(0, eq))] = strip(tok.substr(eq + 1));
  }
  return kv;
}

double get_d(const std::unordered_map<std::string, std::string>& kv,
  const std::string& k, double def)
{
  auto it = kv.find(k);
  return it == kv.end() ? def : std::stod(it->second);
}

std::string get_s(const std::unordered_map<std::string, std::string>& kv,
  const std::string& k, const std::string& def)
{
  auto it = kv.find(k);
  return it == kv.end() ? def : it->second;
}

} // namespace

class PolarizedFusionSource : public openmc::Source {
public:
  explicit PolarizedFusionSource(const std::string& params)
  {
    auto kv = parse_params(params);

    // --- polarization weights: warn on user input, hard-assert invariants ---
    double a = get_d(kv, "a", 1.0 / 3.0);
    double b = get_d(kv, "b", 1.0 / 3.0);
    double c = get_d(kv, "c", 1.0 / 3.0);
    double sum = a + b + c;
    if (std::fabs(sum - 1.0) > 1e-6) {
      std::fprintf(stderr,
        "[PolarizedFusionSource] warning: a+b+c=%.8g != 1; renormalizing.\n", sum);
    }
    mw_ = spf::make_mode_weights(a, b, c); // renormalizes + asserts invariants

    // --- field model: build a MagneticField (B-hat(x); see spf_field.hpp) ---
    const std::string bmode = get_s(kv, "bmode", "toroidal");
    if (bmode == "toroidal") {
      field_ = std::make_unique<spf::ToroidalField>();
    } else if (bmode == "constant") {
      spf::Vec3 b0 {get_d(kv, "bx", 0.0), get_d(kv, "by", 0.0), get_d(kv, "bz", 1.0)};
      field_ = std::make_unique<spf::ConstantField>(b0);
    } else if (bmode == "angled") {
      field_ = std::make_unique<spf::AngledField>(
        get_d(kv, "alpha", 0.0), get_d(kv, "beta", 0.0));
    } else {
      throw std::invalid_argument(
        "bmode must be 'toroidal', 'constant', or 'angled'");
    }

    // --- spatial shape ---
    shape_ = get_s(kv, "shape", "ring");
    if (shape_ == "ring") {
      r0_ = get_d(kv, "r0", 100.0);
      z0_ = get_d(kv, "z0", 0.0);
    } else if (shape_ == "plasma") {
      R0_ = get_d(kv, "R0", 100.0);
      aminor_ = get_d(kv, "aminor", 50.0);
      if (!(aminor_ > 0.0) || !(R0_ - aminor_ > 0.0))
        throw std::invalid_argument("plasma needs aminor>0 and R0>aminor");
      pmax_ = R0_ + aminor_;
    } else {
      throw std::invalid_argument("shape must be 'ring' or 'plasma'");
    }

    energy_ = get_d(kv, "energy", 14.1e6);
  }

  openmc::SourceSite sample(uint64_t* seed) const override
  {
    auto rng = [seed]() { return openmc::prn(seed); };

    // 1. birth position (cm)
    double px, py, pz;
    if (shape_ == "ring") {
      double phi = 2.0 * spf::PI * rng();
      px = r0_ * std::cos(phi);
      py = r0_ * std::sin(phi);
      pz = z0_;
    } else {
      // parabolic plasma: sample (p,z0) ∝ (1-ρ²)·p over the disk (volume-weighted),
      // φ uniform. ρ² = ((p-R0)² + z0²)/aminor².
      double p, z0w, rho2;
      for (;;) {
        p = (R0_ - aminor_) + 2.0 * aminor_ * rng();
        z0w = -aminor_ + 2.0 * aminor_ * rng();
        rho2 = ((p - R0_) * (p - R0_) + z0w * z0w) / (aminor_ * aminor_);
        if (rho2 <= 1.0 && rng() <= (1.0 - rho2) * p / pmax_)
          break;
      }
      double phi = 2.0 * spf::PI * rng();
      px = p * std::cos(phi);
      py = p * std::sin(phi);
      pz = z0w;
    }

    // 2. local magnetic field direction B-hat(x) (see spf_field.hpp)
    const spf::Vec3 Bhat = field_->bhat({px, py, pz});

    // 3. birth direction from the verified angular sampler
    spf::Vec3 u = spf::sample_global_direction(mw_, Bhat, rng);

    // 4. assemble site
    openmc::SourceSite site;
    site.particle = openmc::ParticleType::neutron();
    site.r.x = px;
    site.r.y = py;
    site.r.z = pz;
    site.u.x = u.x;
    site.u.y = u.y;
    site.u.z = u.z;
    site.E = energy_;
    site.wgt = 1.0;
    site.delayed_group = 0;
    return site;
  }

private:
  // All set in the constructor and only read in sample() -> thread-safe.
  // (field_->bhat() is const and has no mutable state.)
  spf::ModeWeights mw_;
  std::unique_ptr<const spf::MagneticField> field_;
  std::string shape_;
  double r0_ {100.0}, z0_ {0.0};
  double R0_ {100.0}, aminor_ {50.0}, pmax_ {150.0};
  double energy_ {14.1e6};
};

// dlopen factory (external C linkage required for dlsym to find it).
extern "C" std::unique_ptr<openmc::Source> openmc_create_source(
  std::string parameters)
{
  return std::make_unique<PolarizedFusionSource>(parameters);
}
