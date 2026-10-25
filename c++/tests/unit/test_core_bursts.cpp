// core/bursts (python core/bursts.py) and the blind hop cadence (receivers/digital_vtx/blind.py
// hop_period): rectangle fit through in-band dips, seeds returned unchanged on noise.
#include <doctest/doctest.h>

#include <random>

#include "core/bursts.hpp"
#include "receivers/digital_vtx/blind.hpp"

using namespace si;

TEST_CASE("bursts: rect_fit spans in-band dips") {
  dvec d(40, 0.0);
  for (int k = 0; k < 12; ++k)
    for (double v : {0.0, 0.0, -2.5, 0.0, 0.0}) d.push_back(3 + v);
  d.insert(d.end(), 40, 0.0);
  const auto r = bursts::rect_fit(d);
  CHECK(r.lo == 40);
  CHECK(r.hi == 100);
}

TEST_CASE("bursts: energy_bursts on noise returns the seeds") {
  std::mt19937_64 g(3);
  std::normal_distribution<float> n(0.f, 1.f);
  cfvec x(1 << 18);
  for (auto& v : x) v = cf(n(g), n(g));
  const bursts::TFPower tf(x, 30.72e6);
  CHECK(tf.L == 512);
  const bursts::Runs seeds = {{1e-3, 2e-3}, {3.5e-3, 4.4e-3}};
  const auto [runs, info] = bursts::energy_bursts(tf, seeds, -1.5e6, 7.5e6, 13e6);
  CHECK(!info.used);
  CHECK(runs == seeds);
}

TEST_CASE("blind: hop_period of an O4 hop train re-timed mid-capture") {
  // 1.37, 3.87, 6.38, 8.88 ms on one phase, 17.08, 24.57 on another (D_20260609-151042 ch0)
  const auto T = rx::digital_vtx::blind::hop_period({1.37e-3, 3.87e-3, 6.38e-3, 8.88e-3, 17.08e-3, 24.57e-3});
  REQUIRE(T);
  CHECK(std::fabs(*T - 2.505e-3) < 10e-6);
  CHECK(!rx::digital_vtx::blind::hop_period({1.37e-3, 8.87e-3, 17.08e-3, 24.58e-3}));  // two pairs: ambiguous
}
