// detect/merge merge_pulse_trains and the features/measure OFDM gate for short pulses (python tests/test_core.py
// test_pulse_train_pieces_merge, test_short_pulses_are_not_ofdm; live LIVE_PULSED_3040MHz).
#include <doctest/doctest.h>

#include "detect/merge.hpp"
#include "features/measure.hpp"

using namespace si;

namespace {
dvec train(double shift) {
  dvec t;
  for (int k = 0; k < 9; ++k) t.push_back(1.548e-3 + 4e-3 * k + shift);
  return t;
}
Signal sig(double f0, double f1, const dvec& pulses, bool host) {
  Signal s;
  s.f0 = f0;
  s.f1 = f1;
  s.snr_db = 3.0;
  s.detection["description"] = "x";
  for (double p : pulses) s.bursts.push_back(Box{p, p + 10e-6, f0, f1, 5, 3, "fine"});
  if (host) {
    s.features["pri_s"] = 4e-3;
    s.features["_pulse_t_s"] = pulses;
  }
  return s;
}
}  // namespace

TEST_CASE("merge_pulse_trains: pieces of one train and its slivers become one signal") {
  const dvec t = train(0);
  Signal a = sig(3036.8e6, 3039.5e6, t, true), b = sig(3040.1e6, 3041.8e6, train(2e-6), true);
  Signal sliver = sig(3036.7e6, 3037.1e6, {}, false);
  sliver.bursts = {Box{t[3] - 5e-6, t[3] + 12e-6, 3036.7e6, 3037.1e6, 4, 2, "fine"}};
  Signal off = sig(3038.4e6, 3038.7e6, {}, false);
  off.bursts = {Box{0.0, 83e-6, 3038.4e6, 3038.7e6, 4, 2, "fine"}};
  Signal other = sig(3029e6, 3030e6, train(1e-3), true);
  auto [out, merged] = detect::merge_pulse_trains({a, off, b, sliver, other});
  REQUIRE(out.size() == 3);
  REQUIRE(merged.size() == 1);
  const Signal& m = out[merged[0]];
  CHECK(m.f0 == 3036.7e6);
  CHECK(m.f1 == 3041.8e6);
  CHECK(m.bursts.size() == 19);
  CHECK(m.detection["description"].get<std::string>().find("3 energy detections with one pulse train") !=
        std::string::npos);
  CHECK(out[0].f0 == 3038.4e6);  // off-pulse burst kept, in input order
  CHECK(out[1].f0 == 3029e6);    // other PRI phase kept
}

TEST_CASE("infer_modulation: an autocorrelation lag longer than a quarter of the pulses is not OFDM") {
  Signal s;
  s.f0 = 3036.8e6;
  s.f1 = 3041.8e6;
  s.snr_db = 8.0;
  Json f = {{"ofdm_scs_hz", 180410.0}, {"acf_peak", 0.08}, {"acf_peak_over_median", 17.1}, {"acf_peak_lag_s", 5.54e-6}};
  CHECK(features::infer_modulation(s, f, 1e3).first == "OFDM");
  f["pulse_n"] = 9;
  f["pulse_width_s_median"] = 200e-6;
  CHECK(features::infer_modulation(s, f, 1e3).first == "OFDM");
  f["pulse_width_s_median"] = 11e-6;
  CHECK(features::infer_modulation(s, f, 1e3).first != "OFDM");
}
