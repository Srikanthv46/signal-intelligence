// Pipeline plugin: FMCW radio altimeters, 4200-4400 MHz (research; python: receivers/radalt/plugin.py).
#include <algorithm>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/radalt/detector.hpp"

namespace si::rx::radalt {

namespace {

class RadaltDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "radalt"; }
  std::string package() const override { return "radalt"; }
  std::string protocol() const override { return "FMCW radio altimeter"; }
  std::string key() const override { return "radalt"; }
  std::set<std::string> bands() const override { return {"radalt"}; }
  std::string method() const override { return "receivers.radalt: STFT ridge crossings of the FMCW sweep (research)"; }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"gate: the capture overlaps 4200-4400 MHz, fs >= 20 MSps",
            "ridge: 128-point Hann FFT every 64 samples, the peak bin (parabolic) >= 10 dB over the frame median",
            "crossings: ridge frames joined into tracks (<= 4 bins per elapsed frame, gaps <= 8 frames), >= 5 frames "
            "on in >= half the frames spanned, linear fit R^2 >= 0.99 (a V split at its extreme), same-direction "
            "pieces within 20 us merged, >= 15 MHz covered",
            "sweep: rate = median |slope| (2-1500 MHz/ms); same-direction spacing in whole periods -> period; "
            "up / down pairs -> turnarounds, sweep bandwidth and the triangle period 2 B / r",
            "confidence 0.7, 0.8 when the bandwidth is 95-190 MHz or the period matches an ITU-R M.2059 type "
            "(20, 6.67, 7, 1 ms)"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    if (std::min(ctx.cap.fs, ctx.cap.bandwidth_hz.value_or(ctx.cap.fs)) < 20e6) return {};
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc);
  }
  std::string group_key(const Emission& e) const override { return e.label; }

  void finish(Signal& sig, const std::vector<Emission>& em, DecodeContext&) const override {
    const Json& F = em[0].features;
    sig.detection["description"] = "defined by the radalt receiver: " + py_str(F["n_crossings"]) +
                                   " straight frequency ramp(s) at " + py_str(F["rate_mhz_per_ms"]) +
                                   " MHz/ms through the capture, not by energy alone";
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission& e = em[0];
    const Json& F = e.features;
    std::string what = F["shape"].get<std::string>() + ", " + py_str(F["rate_mhz_per_ms"]) + " MHz/ms";
    if (F.contains("sweep_bw_mhz")) what += ", " + py_str(F["sweep_bw_mhz"]) + " MHz sweep";
    if (F.contains("period_ms")) what += ", period " + py_str(F["period_ms"]) + " ms";
    std::string ev = py_str(F["n_crossings"]) + " linear crossing(s) (R^2 >= " + py_str(F["r2_min"]) + ") at " +
                     py_str(F["rate_mhz_per_ms"]) + " MHz/ms, seen " + py_str(F["seen_mhz"][0]) + "-" +
                     py_str(F["seen_mhz"][1]) + " MHz";
    if (F.contains("f_top_mhz")) ev += "; top turnaround " + py_str(F["f_top_mhz"]) + " MHz";
    if (F.contains("f_bottom_mhz")) ev += "; bottom turnaround " + py_str(F["f_bottom_mhz"]) + " MHz";
    return result("partial", method(), steps(), F,
                  "sweep measured (an FMCW altimeter carries no payload); the platform is not identified",
                  "Radio altimeter (FMCW)-like (" + what + ")", e.confidence, "FMCW", ev,
                  std::make_pair(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<RadaltDecoder>()); }

}  // namespace si::rx::radalt
