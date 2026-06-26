// Field-parity driver for spf_field.hpp. Prints B-hat(x) for a chosen field
// model at positions read from stdin, so tests/test_field_parity.py can compare
// the C++ field to its Python mirror (python/spf_mirror.py) to <= 1e-12.
//
// Links neither libopenmc nor random_lcg (spf_field.hpp is OpenMC-free).
// Build (see tests/conftest.py for the canonical call):
//   g++ -O2 -std=c++17 -I<openmc>/spf_prototype/src spf_field_driver.cpp -o spf_field_driver
//
// Usage:  spf_field_driver <bmode> [params...]   ; positions on stdin: "x y z" per line
//   toroidal
//   constant bx by bz
//   angled   alpha beta        (radians)
//   fieldmap <stem>            (reads <stem>.meta + <stem>.bin)
// Output: one line "bx by bz" per input position (17 sig figs).
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "spf_field.hpp"
#include "spf_fieldmap.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <toroidal|constant bx by bz|angled alpha beta>\n",
      argv[0]);
    return 2;
  }
  const std::string bmode = argv[1];
  std::unique_ptr<const spf::MagneticField> field;
  if (bmode == "toroidal") {
    field = std::make_unique<spf::ToroidalField>();
  } else if (bmode == "constant") {
    if (argc < 5) { std::fprintf(stderr, "constant needs bx by bz\n"); return 2; }
    field = std::make_unique<spf::ConstantField>(
      spf::Vec3 {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])});
  } else if (bmode == "angled") {
    if (argc < 4) { std::fprintf(stderr, "angled needs alpha beta\n"); return 2; }
    field = std::make_unique<spf::AngledField>(std::atof(argv[2]), std::atof(argv[3]));
  } else if (bmode == "fieldmap") {
    if (argc < 3) { std::fprintf(stderr, "fieldmap needs <stem>\n"); return 2; }
    field = std::make_unique<spf::FieldMapField>(std::string(argv[2]));
  } else {
    std::fprintf(stderr, "unknown bmode '%s'\n", bmode.c_str());
    return 2;
  }

  double x, y, z;
  while (std::cin >> x >> y >> z) {
    const spf::Vec3 b = field->bhat({x, y, z});
    std::printf("%.17g %.17g %.17g\n", b.x, b.y, b.z);
  }
  return 0;
}
