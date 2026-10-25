// core/npcompat np::sum / np::mean beyond numpy's 8192-element reduction buffer
// (tests/vectors/analog_vtx/gen_core_sum.py): bit-identical for float32, float64, complex64.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/npcompat.hpp"

using namespace si;

TEST_CASE("npcompat: sum/mean in NPY_BUFSIZE blocks, bit-identical") {
  std::ifstream f(std::string(SI_TESTVEC) + "/analog_vtx/core_sum.bin", std::ios::binary);
  REQUIRE_MESSAGE(f.good(), "missing vector (run tests/vectors/analog_vtx/gen_core_sum.py)");
  for (size_t n : {7, 100, 8192, 8193, 30000}) {
    CAPTURE(n);
    fvec x32(n);
    dvec x(n);
    cfvec z(n);
    float r32[2];
    double r64[2];
    cf rz;
    f.read(reinterpret_cast<char*>(x32.data()), n * 4);
    f.read(reinterpret_cast<char*>(x.data()), n * 8);
    f.read(reinterpret_cast<char*>(z.data()), n * 8);
    f.read(reinterpret_cast<char*>(r32), 8);
    f.read(reinterpret_cast<char*>(r64), 16);
    f.read(reinterpret_cast<char*>(&rz), 8);
    REQUIRE(f.good());
    CHECK(np::sum(x32) == r32[0]);
    CHECK(np::mean(x32) == r32[1]);
    CHECK(np::sum(x) == r64[0]);
    CHECK(np::mean(x) == r64[1]);
    CHECK(np::pairwise_sum(z.data(), n) == rz);
  }
}
