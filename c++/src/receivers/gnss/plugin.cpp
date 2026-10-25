// Pipeline plugin: GPS L1 C/A satellites, GNSS jamming / spoofing indicators (research; python: receivers/gnss/plugin.py).
#include <cmath>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/plugin.hpp"
#include "receivers/gnss/acquire.hpp"

namespace si::rx::gnss {

namespace {

class GnssDecoder : public Decoder {
 public:
  std::string name() const override { return "gnss"; }
  std::string package() const override { return "gnss"; }
  std::string protocol() const override { return "GPS L1 C/A"; }
  bool capture_level() const override { return true; }
  std::vector<std::string> after() const override { return {"lte"}; }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    const std::vector<std::string> steps = {
        "gate: 1575.42 MHz +-1.1 MHz in the usable band, fs >= 2.1 MSps, >= 16 ms",
        "block FFT cropped to +-1.024 MHz around L1 / (1 + ppm) -> 2.048 MSps; gate at 1 sample / chip",
        "32 PRNs (IS-GPS-200 Gold codes, ICD self-test octals), Doppler 125 Hz bins over +-8.2 kHz with the "
        "receiver ppm (+-25 kHz without), 4 ms coherent folds x 4 non-coherent, one FFT per fold shared by "
        "all PRNs; threshold chi2(2M) at P_fa 1e-3 per PRN per block",
        "each hit refined at 2 samples / chip over the whole block: code phase, Doppler, C/N0 (parabolic "
        "straddle correction + measured implementation loss)",
        "indicators: in-band J/N over the capture's floor outside 1559-1610 MHz, narrow lines near L1 (receiver "
        "lines excluded); spoofing S1 C/N0 > 55, S2 uniform C/N0 (>= 4 PRNs > 45 dB-Hz within 1.5 dB), S3 J/N "
        ">= 1.5 dB with >= 4 strong PRNs, S4 second correlation peaks in >= 2 PRNs",
        "swept jammer (personal privacy device): 32-point FFT ridge over 400 us, sweep period from its "
        "autocorrelation (0.5-50 us), fold coherence >= 0.6 and span >= 1 MHz (fs >= 30 MSps)"};
    const Json* cal = ctx.cache.get<Json>("ppm");
    if (!cal || cal->is_null() || cal->empty()) cal = ctx.cache.get<Json>("sweep_ppm");
    std::optional<double> ppm;
    if (cal && cal->is_object() && cal->contains("ppm") && (*cal)["ppm"].is_number()) ppm = (*cal)["ppm"].get<double>();
    dvec excl;
    if (const Json* ln = ctx.cache.get<Json>("receiver_lines"))
      for (const auto& l : *ln) excl.push_back(l["freq_hz"].get<double>());
    const auto r = analyse(ctx.ch, ppm, excl);
    if (!r) return {};
    const Json& sats = (*r)["satellites"];
    const Json& ind = (*r)["indicators"];
    const size_t n = sats.size();
    Json names = Json::array();
    for (const Json& s : sats) names.push_back(s["prn"]);
    const std::pair<double, double> span = {kFL1 - 1.023e6, kFL1 + 1.023e6};
    std::string label, ev;
    double conf;
    if (ind["spoofing"]["verdict"].get<bool>()) {
      const int64_t k = ind["spoofing"]["n"].get<int64_t>();
      label = fmt("GNSS spoofing-like (GPS L1 C/A, %zu satellites)", n);
      conf = k >= 3 ? 0.75 : 0.6;
      Json fired = Json::array();
      for (auto it = ind["spoofing"].begin(); it != ind["spoofing"].end(); ++it)
        if (it.value().is_boolean() && it.value().get<bool>() && it.key() != "verdict") fired.push_back(it.key());
      ev = fmt("%zu GPS satellites ", n) + py_repr(names) + "; spoofing indicators " + py_repr(fired);
    } else if (ind["jamming"]["verdict"].get<bool>()) {
      const Json& itf = (*r)["interference"];
      const bool sw = itf.contains("swept");
      label = sw ? "GNSS jamming-like (GPS L1, swept PPD)" : "GNSS jamming-like (GPS L1)";
      conf = sw ? 0.85 : 0.6;
      ev = "in-band J/N " + py_str(itf["jn_db"]) + fmt(" dB, %zu narrow line(s) near L1; %zu satellite(s) acquired",
                                                         itf.contains("lines") ? itf["lines"].size() : size_t(0), n);
      if (sw) {
        const Json& w = itf["swept"];
        const double c = w["centre_offset_mhz"].get<double>();
        ev += "; swept carrier (" + w["shape"].get<std::string>() + ") every " + py_str(w["period_us"]) + " us over " +
              py_str(w["span_mhz"]) + " MHz, centre " + (std::signbit(c) ? "" : "+") + py_str(w["centre_offset_mhz"]) +
              " MHz from L1" + (w["covers_l1"].get<bool>() ? "" : " (does not cross L1)");
      }
    } else if (n >= 2 || (n == 1 && sats[0]["cn0_dbhz"].get<double>() >= 40)) {
      label = fmt("GNSS: GPS L1 C/A (%zu satellite%s)", n, n > 1 ? "s" : "");
      conf = 0.9;
      ev = fmt("%zu PRN(s) above the chi2 threshold (P_fa 1e-3 per PRN): ", n);
      for (size_t i = 0; i < n; ++i)
        ev += (i ? ", " : "") + std::string("PRN ") + py_str(sats[i]["prn"]) + " " + py_str(sats[i]["cn0_dbhz"]) + " dB-Hz";
    } else {
      ctx.cache.set<std::string>("gnss_note", fmt("L1 searched: %zu satellite(s) above threshold, no jamming / "
                                                  "spoofing indicator", n));
      return {};
    }
    DecodeResult res = result("partial", "receivers.gnss: GPS L1 C/A acquisition (research)", steps, *r,
                              "satellites acquired by code correlation; the navigation message is not decoded", label,
                              conf, std::string("BPSK(1) DSSS 1.023 Mcps (GPS C/A)"), ev, span);
    const double dur = static_cast<double>(ctx.ch.L()) / ctx.ch.fs();
    return {receiver_signal({{0.0, dur}}, span.first, span.second, res, std::nullopt, name())};
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<GnssDecoder>()); }

}  // namespace si::rx::gnss
