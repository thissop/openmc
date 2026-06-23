// Tier-1 standalone driver for spf_sampler.hpp. Links OpenMC's real prn
// (src/random_lcg.cpp) but NOT libopenmc, so it tests the sampler MATH in
// isolation. Emits sampled birth directions for the Python statistical tests
// and the C++/Python bit-for-bit parity check.
//
// Build (see python/build_and_run.py / tests/conftest.py for the canonical call):
//   g++ -O2 -std=c++17 -I<openmc>/include -I<openmc>/spf_prototype/src \
//       standalone_driver.cpp <openmc>/src/random_lcg.cpp -o spf_driver
//
// Usage:  spf_driver a b c Bx By Bz N seed [invcdf|reject]
//   prints N lines:  cz ux uy uz
//     cz       = u·B̂  (= cosθ, the angle to the field) -- for per-mode chi-square
//     ux uy uz = global unit birth direction          -- for isotropy/steering
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "openmc/random_lcg.h"  // openmc::prn(uint64_t*)

#include "spf_sampler.hpp"

int main(int argc, char** argv)
{
  if (argc < 9) {
    std::fprintf(stderr,
      "usage: %s a b c Bx By Bz N seed [invcdf|reject]\n"
      "  prints N lines: cz ux uy uz  (cz=u.Bhat=cos(theta); uxyz=global unit dir)\n",
      argv[0]);
    return 2;
  }

  const double a = std::atof(argv[1]);
  const double b = std::atof(argv[2]);
  const double c = std::atof(argv[3]);
  spf::Vec3 Bhat {std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  const long N = std::atol(argv[7]);
  uint64_t seed = std::strtoull(argv[8], nullptr, 10);
  const bool reject = (argc >= 10 && std::string(argv[9]) == "reject");

  const spf::ModeWeights m = spf::make_mode_weights(a, b, c);
  const spf::Vec3 Bn = spf::normalized(Bhat);

  // RNG closure: OpenMC's real PCG generator over a single seed stream.
  auto rng = [&seed]() { return openmc::prn(&seed); };

  for (long i = 0; i < N; ++i) {
    const spf::Vec3 u = spf::sample_global_direction(m, Bn, rng, reject);
    const double cz = spf::dot(u, Bn);
    std::printf("%.17g %.17g %.17g %.17g\n", cz, u.x, u.y, u.z);
  }
  return 0;
}
