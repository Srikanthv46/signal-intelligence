// receivers/jammer: synthetic swept carriers and look-alikes (as python tests/receivers/test_det_jammer.py), and
// the bench capture of a real jammer (testdata/live/LIVE_JAMMER_SWEEP_2460MHz). Python/C++ equality of the report
// is checked by tools/parity.py.
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <random>
#include <string>

#include "io/sigmf.hpp"
#include "receivers/jammer/detector.hpp"

using namespace si;
namespace jm = si::rx::jammer;

namespace {

constexpr double kFs = 61.44e6, kFc = 2.45e9;
constexpr size_t kN = 2000000;

// carrier swept between f_lo and f_hi (Hz from the LO) `rate` times a second, snr_db over the noise in 1 MHz,
// blanked beyond +-0.45 fs (the receiver's anti-alias filter); unit-variance complex noise
cfvec sweep(double rate, double f_lo, double f_hi, bool sawtooth, double snr_db, unsigned seed) {
  std::mt19937 g(seed);
  std::normal_distribution<double> n01(0.0, std::sqrt(0.5));
  const double a = std::pow(10.0, snr_db / 20) * std::sqrt(1e6 / kFs);
  cfvec x(kN);
  double ph = 0;
  for (size_t i = 0; i < kN; ++i) {
    const double u = std::fmod(i / kFs * rate, 1.0);
    const double f = f_lo + (f_hi - f_lo) * (sawtooth ? u : 1 - std::fabs(2 * u - 1));
    ph += 2 * kPi * f / kFs;
    const double amp = std::fabs(f) < 0.45 * kFs ? a : 0.0;
    x[i] = cf(static_cast<float>(amp * std::cos(ph) + n01(g)), static_cast<float>(amp * std::sin(ph) + n01(g)));
  }
  return x;
}

cfvec noise_cw(double cw_amp, unsigned seed) {
  std::mt19937 g(seed);
  std::normal_distribution<double> n01(0.0, std::sqrt(0.5));
  cfvec x(kN);
  for (size_t i = 0; i < kN; ++i) {
    const double ph = 2 * kPi * 3e6 * i / kFs;
    x[i] = cf(static_cast<float>(cw_amp * std::cos(ph) + n01(g)), static_cast<float>(cw_amp * std::sin(ph) + n01(g)));
  }
  return x;
}

}  // namespace

TEST_CASE("jammer: triangle sweep 10 kHz over 45 MHz") {
  const auto em = jm::detect(sweep(10e3, -25e6, 20e6, false, 10, 1), kFs, kFc);
  REQUIRE(em.size() == 1);
  const Json& F = em[0].features;
  CHECK(std::fabs(F["sweep_rate_hz"].get<double>() / 10e3 - 1) < 50e-6);
  CHECK(std::fabs(F["sweep_low_mhz"].get<double>() - (kFc - 25e6) / 1e6) < 1.5);
  CHECK(std::fabs(F["sweep_high_mhz"].get<double>() - (kFc + 20e6) / 1e6) < 1.5);
  CHECK(F["shape"] == "up-down");
  CHECK(F["beyond_capture"].empty());
  CHECK(em[0].label == "swept_carrier");
}

TEST_CASE("jammer: sawtooth sweep 3 kHz over 30 MHz") {
  const auto em = jm::detect(sweep(3e3, -15e6, 15e6, true, 10, 2), kFs, kFc);
  REQUIRE(em.size() == 1);
  CHECK(std::fabs(em[0].features["sweep_rate_hz"].get<double>() / 3e3 - 1) < 50e-6);
  CHECK(em[0].features["shape"] == "sawtooth");
}

TEST_CASE("jammer: sweep leaving the band below") {
  const auto em = jm::detect(sweep(10e3, -35e6, 15e6, false, 20, 3), kFs, kFc);
  REQUIRE(em.size() == 1);
  CHECK(em[0].features["beyond_capture"] == Json::array({"below"}));
  CHECK(em[0].f_lo_mhz == doctest::Approx((kFc - kFs / 2) / 1e6));
}

TEST_CASE("jammer: nothing in noise, a CW tone or a narrow sweep") {
  CHECK(jm::detect(noise_cw(0.0, 4), kFs, kFc).empty());
  CHECK(jm::detect(noise_cw(30.0, 5), kFs, kFc).empty());
  CHECK(jm::detect(sweep(10e3, -2e6, 2e6, false, 10, 6), kFs, kFc).empty());
}

TEST_CASE("jammer: bench captures LIVE_JAMMER_SWEEP_2460MHz(_jitter)") {
  for (const char* name : {"LIVE_JAMMER_SWEEP_2460MHz", "LIVE_JAMMER_SWEEP_2460MHz_jitter"}) {
  CAPTURE(name);
  const std::string p = std::string(SI_TESTDATA) + "/live/" + name + "/" + name + ".sigmf-meta";
  if (!std::ifstream(p).good()) continue;
  const auto cap = load_sigmf(p);
  const auto em = jm::detect(cap.iq(), cap.fs, cap.fc);
  REQUIRE(em.size() == 1);
  const Json& F = em[0].features;
  CHECK(std::fabs(F["sweep_rate_hz"].get<double>() - 10000.15) < 0.5);   // long-baseline reference 10000.09-.21
  CHECK(F["beyond_capture"] == Json::array({"below"}));
  CHECK(F["sweep_high_mhz"].get<double>() > 2468);
  CHECK(F["sweep_high_mhz"].get<double>() < 2480);
  CHECK(F["shape"] == "up-down");
  CHECK(F["flybacks"] == 0);
  }
}
