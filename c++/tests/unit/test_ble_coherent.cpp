// Coherent BLE GFSK detector + CRC Chase (rc_link/gfsk_coherent): research capability vs python on synthetic vectors
// (tests/vectors/ble_coherent/gen_ble_coherent.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/fsk.hpp"
#include "receivers/rc_link/gfsk_coherent.hpp"

using namespace si;
namespace G = si::rx::rc_link::gfsk;

namespace {
Json load() {
  std::ifstream f(std::string(SI_TESTVEC) + "/ble_coherent/cases.json");
  REQUIRE(f);
  return Json::parse(f);
}
cfvec bursts() {
  std::ifstream f(std::string(SI_TESTVEC) + "/ble_coherent/bursts.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec y(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(y.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  return y;
}
std::vector<uint8_t> bitvec(const std::string& s) {
  std::vector<uint8_t> b(s.size());
  for (size_t i = 0; i < s.size(); ++i) b[i] = s[i] == '1';
  return b;
}
}  // namespace

TEST_CASE("ble coherent: access address, demodulation and decode == python") {
  const Json ref = load();
  const cfvec all = bursts();
  std::vector<uint8_t> aa(32);
  for (int i = 0; i < 32; ++i) aa[i] = (0x8E89BED6u >> i) & 1;
  for (const Json& c : ref["cases"]) {
    const size_t o = c["offset"].get<size_t>(), n = c["n"].get<size_t>();
    const cfvec y(all.begin() + static_cast<ptrdiff_t>(o), all.begin() + static_cast<ptrdiff_t>(o + n));
    const auto got = G::find_access_address(G::diff_products(y), aa);
    REQUIRE(got.size() == c["aa"].size());
    for (size_t i = 0; i < got.size(); ++i) {
      CHECK(got[i].first == c["aa"][i][0].get<int64_t>());
      CHECK(got[i].second == doctest::Approx(c["aa"][i][1].get<double>()).epsilon(1e-9));
    }
    const int64_t s0 = c["s0"].get<int64_t>();
    const int64_t nb = (got[0].first - s0) / 8;
    const auto dm = G::demodulate(y, s0, nb + 32 + 16 + 8 * 255 + 24);
    REQUIRE(dm);
    const auto want = bitvec(c["bits"].get<std::string>());
    REQUIRE(dm->bits.size() == want.size());
    size_t mism = 0;
    for (size_t i = 0; i < want.size(); ++i) mism += dm->bits[i] != want[i];
    CHECK(mism == 0);
    CHECK(dm->h == doctest::Approx(c["h"].get<double>()).epsilon(1e-9));
    CHECK(dm->cfo_hz == doctest::Approx(c["cfo_hz"].get<double>()).epsilon(1e-6));
    CHECK(dm->dt == c["dt"].get<int>());
    CHECK(rx::rc_link::detail::ble_coherent_decode(y, 8e6, 37).dump() == c["coherent"].dump());
    CHECK(rx::rc_link::detail::ble_try_decode(y, 8e6, 2402e6).dump() == c["try_decode"].dump());
  }
}

TEST_CASE("ble coherent: CRC syndrome and Chase == python") {
  const Json ref = load();
  const size_t n_pdu = ref["n_pdu"].get<size_t>();
  for (const Json& c : ref["chase"]) {
    const auto b = bitvec(c["bits"].get<std::string>());
    const auto [syn, contrib] = G::crc_syndrome_parts(b, n_pdu);
    CHECK(syn == c["syndrome"].get<int64_t>());
    CHECK(Json(contrib).dump() == ref["contrib"].dump());
    const auto got = G::chase(b, n_pdu, c["rel"].get<dvec>(), 40);
    if (c["out"].is_null()) {
      CHECK(!got);
    } else {
      REQUIRE(got);
      std::string s;
      for (size_t i = 0; i < n_pdu + 24; ++i) s += got->first[i] ? '1' : '0';
      CHECK(s == c["out"][0].get<std::string>());
      CHECK(got->second == c["out"][1].get<int>());
    }
  }
}
