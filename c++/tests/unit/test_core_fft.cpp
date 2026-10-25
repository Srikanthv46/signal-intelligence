// core/fft against numpy 2.x / scipy (tests/vectors/core/gen_fft.py): bit-identical.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/fft.hpp"

using namespace si;

namespace {

template <class T>
std::vector<T> load(const std::string& name, size_t n, size_t count = 0) {
  std::vector<T> v(count ? count : n);
  std::ifstream f(std::string(SI_TESTVEC) + "/core/fft/" + name + "_" + std::to_string(n) + ".bin", std::ios::binary);
  REQUIRE_MESSAGE(f.good(), "missing vector " << name << " (run tests/vectors/core/gen_fft.py)");
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
  REQUIRE(f.gcount() == static_cast<std::streamsize>(v.size() * sizeof(T)));
  return v;
}

template <class T>
size_t mismatches(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) return std::max(a.size(), b.size());
  size_t bad = 0;
  for (size_t i = 0; i < a.size(); ++i) bad += !(a[i] == b[i]);
  return bad;
}

}  // namespace

TEST_CASE("fft: scipy and numpy semantics, bit-identical") {
  for (size_t n : {1000, 1967}) {
    CAPTURE(n);
    const auto x = load<cf>("in_c64", n);
    const auto r = load<float>("in_f32", n);
    const auto d = load<cd>("in_c128", n);
    const size_t h = n / 2 + 1;
    CHECK(mismatches(fft::fft(x), load<cf>("sp_fft", n)) == 0);
    CHECK(mismatches(fft::ifft(x), load<cf>("sp_ifft", n)) == 0);
    CHECK(mismatches(fft::rfft(r), load<cf>("sp_rfft", n, h)) == 0);
    CHECK(mismatches(fft::np_rfft(r), load<cf>("np_rfft", n, h)) == 0);
    CHECK(mismatches(fft::np_fft(x), load<cf>("np_fft", n)) == 0);
    CHECK(mismatches(fft::np_ifft(x), load<cf>("np_ifft", n)) == 0);
    CHECK(mismatches(fft::np_fft(x, 2 * n), load<cf>("np_fft_pad", n, 2 * n)) == 0);
    CHECK(mismatches(fft::fft(d), load<cd>("fft_c128", n)) == 0);
    CHECK(mismatches(fft::ifft(d), load<cd>("ifft_c128", n)) == 0);
  }
}
