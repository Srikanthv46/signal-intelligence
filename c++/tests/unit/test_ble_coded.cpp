// BLE LE Coded PHY (BT5 Long Range Remote ID, rc_link/ble_coded): research capability vs python on synthetic vectors
// (tests/vectors/ble_coded/gen_ble_coded.py), incl. the whole rc_link pipeline on a 15.36 MSps capture.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/pyfmt.hpp"
#include "receivers/rc_link/ble_coded.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/fsk.hpp"
#include "receivers/rc_link/gfsk_coherent.hpp"

using namespace si;
namespace C = si::rx::rc_link::coded;

namespace {
Json ref() {
  std::ifstream f(std::string(SI_TESTVEC) + "/ble_coded/cases.json");
  REQUIRE(f);
  return Json::parse(f);
}
cfvec load_cf(const std::string& n) {
  std::ifstream f(std::string(SI_TESTVEC) + "/ble_coded/" + n, std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t m = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec y(m);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(y.data()), static_cast<std::streamsize>(m * sizeof(cf)));
  return y;
}
std::string bstr(const std::vector<uint8_t>& b) {
  std::string s;
  for (uint8_t v : b) s += v ? '1' : '0';
  return s;
}
std::vector<uint8_t> aa_bits() {
  std::vector<uint8_t> a(32);
  for (int i = 0; i < 32; ++i) a[i] = (0x8E89BED6u >> i) & 1;
  return a;
}
}  // namespace

TEST_CASE("ble coded: viterbi, known symbols, extended header == python") {
  const Json r = ref();
  CHECK(bstr(C::viterbi(r["viterbi"]["llr"].get<dvec>())) == r["viterbi"]["out"].get<std::string>());
  CHECK(bstr(C::known_symbols(aa_bits(), 0)) == r["known0"].get<std::string>());
  const std::string h = r["parse_extended"]["hex"];
  std::vector<uint8_t> by;
  for (size_t i = 0; i + 1 < h.size(); i += 2) by.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  std::vector<uint8_t> data;
  Json got = C::parse_extended(by, data);
  std::string dh;
  for (uint8_t v : data) dh += fmt("%02x", v);
  got["adv_data"] = dh;
  CHECK(got.dump() == r["parse_extended"]["out"].dump());
}

TEST_CASE("ble coded: preamble search, decode and ble_coded_decode == python") {
  const Json r = ref();
  const cfvec all = load_cf("bursts.cf32");
  for (const Json& c : r["cases"]) {
    const size_t o = c["offset"], n = c["n"];
    const cfvec y(all.begin() + static_cast<ptrdiff_t>(o), all.begin() + static_cast<ptrdiff_t>(o + n));
    const auto f = C::find_coded(rx::rc_link::gfsk::diff_products(y), aa_bits());
    REQUIRE(f.size() == c["find"].size());
    CHECK(f[0].first == c["find"][0][0].get<int64_t>());
    const auto d = C::decode(y, f[0].first, aa_bits());
    REQUIRE(d);
    CHECK(d->ci == c["ci"].get<int>());
    CHECK(d->aa_bit_errors == c["aa_err"].get<int64_t>());
    CHECK(d->h == doctest::Approx(c["h"].get<double>()).epsilon(1e-9));
    CHECK(bstr(d->bits) == c["bits"].get<std::string>());
    CHECK(rx::rc_link::detail::ble_coded_decode(y, 8e6, c["ch"].get<int>()).dump() == c["out"].dump());
  }
}

TEST_CASE("ble coded: whole rc_link pipeline on a BT5 Long Range capture == python") {
  const Json r = ref();
  const cfvec x = load_cf("capture.cf32");
  const auto em = rx::rc_link::detect(x, r["capture"]["fs"], r["capture"]["fc"]);
  Json got = Json::array();
  for (const auto& e : em)
    if (e.decoded.is_object() && e.decoded.contains("phy")) got.push_back(e.decoded);
  CHECK(got.dump() == r["capture"]["decoded"].dump());
}
