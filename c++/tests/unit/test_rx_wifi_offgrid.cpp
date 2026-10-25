// receivers/wifi ofdm_rx scan_plan / offgrid_plan / merge_packets (python tests/test_core.py
// test_wifi_offgrid_link_found / _weak_low_duty; the live links are checked by parity on LIVE_MIMOSA_4880MHz and
// LIVE_MIMOSA_4880MHz_weak).
#include <doctest/doctest.h>

#include <random>

#include "receivers/wifi/ofdm_rx.hpp"

using namespace si;
using namespace si::rx::wifi;

TEST_CASE("wifi scan_plan: 1 MHz centres, standard ones named after their channel") {
  const auto plan = ofdm::scan_plan(5190e6, 61.44e6, 60e6);
  bool ch36 = false, off = false;
  for (const auto& c : plan) {
    if (c.freq_hz == 5180e6) ch36 = c.standard && c.name == "5g36";
    if (c.freq_hz == 5183e6) off = !c.standard && c.name == "5183.0MHz";
  }
  CHECK(ch36);
  CHECK(off);
}

TEST_CASE("wifi offgrid_plan: noise plans nothing") {
  std::mt19937_64 g(1);
  std::normal_distribution<float> n(0.f, 0.70710678f);
  cfvec x(1 << 19);
  for (auto& v : x) v = cf(n(g), n(g));
  CHECK(ofdm::offgrid_plan(x, 61.44e6, 4877e6, {}).empty());
}

TEST_CASE("wifi offgrid_plan: 2^21 noise plans nothing (windowed fallback gated, PLAN G13)") {
  std::mt19937_64 g(2);
  std::normal_distribution<float> n(0.f, 0.70710678f);
  cfvec x(1 << 21);
  for (auto& v : x) v = cf(n(g), n(g));
  CHECK(ofdm::offgrid_plan(x, 61.44e6, 4877e6, {}).empty());
}

TEST_CASE("wifi merge_packets: known kept in order, duplicates dropped, the best duplicate kept") {
  auto pk = [](double t_us, double f, double stf) {
    ofdm::Packet p;
    p.start_us = t_us;
    p.start_sample = static_cast<int64_t>(t_us * 61.44);
    p.f_lo_mhz = f - 8.3;
    p.f_hi_mhz = f + 8.3;
    p.ch_freq_mhz = f;
    p.stf_metric = stf;
    return p;
  };
  const auto m = ofdm::merge_packets({pk(500, 5180, 0.9), pk(100, 5180, 0.9)},
                                     {pk(101, 5181, 0.95), pk(300, 4880, 0.8), pk(301, 4881, 0.92), pk(200, 4880, 0.9)});
  REQUIRE(m.size() == 4);
  CHECK(m[0].start_us == 500);
  CHECK(m[1].start_us == 100);
  CHECK(m[2].start_us == 200);
  CHECK(m[3].start_us == 301);  // 300 and 301 are one PPDU: the higher STF metric stays
}
