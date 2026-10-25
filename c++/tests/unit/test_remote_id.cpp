// Remote ID over Bluetooth LE (research capability): core/remote_id parse_ble_ad and rc_link ble_try_decode vs the
// Python reference on synthetic vectors (tests/vectors/remote_id/gen_remote_id.py; no capture).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/remote_id.hpp"
#include "receivers/rc_link/fsk.hpp"

using namespace si;

namespace {
Json load(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/remote_id/" + name);
  REQUIRE(f);
  return Json::parse(f);
}
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> v;
  for (size_t i = 0; i + 1 < h.size(); i += 2) v.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  return v;
}
}  // namespace

TEST_CASE("remote_id: BLE advertising data parse == python") {
  const Json ref = load("ad.json");
  const Json got = remote_id::parse_ble_ad(unhex(ref["adv_data_hex"].get<std::string>()));
  CHECK(got.dump() == ref["out"].dump());
}

TEST_CASE("remote_id: synthetic BLE advert demodulated and decoded == python") {
  const Json ref = load("ble.json");
  std::ifstream f(std::string(SI_TESTVEC) + "/remote_id/ble_rid.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec y(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(y.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  const Json got = rx::rc_link::detail::ble_try_decode(y, ref["fs"], ref["f_hz"]);
  REQUIRE(got.is_object());
  CHECK(got["crc_ok"] == true);
  CHECK(got["ad"].dump() == ref["out"]["ad"].dump());
  CHECK(got["adv_addr"] == ref["out"]["adv_addr"]);
}
