// Pipeline plugin: swept-carrier jammers (python: receivers/jammer/plugin.py).
#include <algorithm>
#include <memory>

#include "core/bands_india.hpp"
#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/jammer/barrage.hpp"
#include "receivers/jammer/detector.hpp"

namespace si::rx::jammer {

namespace {

class JammerDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "jammer"; }
  std::string package() const override { return "jammer"; }
  std::string protocol() const override { return "swept-carrier jammer"; }
  std::string key() const override { return "jammer"; }
  std::string method() const override { return "receivers.jammer (swept-carrier detector)"; }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"occupancy: Hann STFT power, fs/64 bins, ~1 us rows; cyclic spectrum of per-bin log occupancy, "
            "whitened per bin and summed -> up to 3 sweep-rate candidates (lines >= 10 dB over the local median)",
            "period: the shortest multiple (<= 8) of the line period whose fold has >= 80 % of the best phase "
            "structure, refused if 2-5x that period has more (pulsed); rate corrected from the fold drift "
            "across 2-4 block segments",
            "fold on the period: one dominant carrier per phase moving continuously (<= 1 flyback, <= 2 band "
            "exits), in band >= 30 % of the period, spanning >= 10 MHz -> sweep rate, extent, shape, trajectory"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    // python: schema/requirements.json jammer.min_sample_rate_hz (20 MHz): a >= 10 MHz sweep
    const double rate = std::min(ctx.cap.fs, ctx.cap.bandwidth_hz.value_or(ctx.cap.fs));
    if (rate < 20e6) return {};
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc);
  }
  bool in_range(double, double) const override { return true; }   // jammers are not tied to an allocation
  std::string group_key(const Emission& e) const override { return py_float_repr(e.center_mhz); }

  void finish(Signal& sig, const std::vector<Emission>& em, DecodeContext&) const override {
    const Json& F = em[0].features;
    sig.detection["description"] =
        fmt("defined by the jammer receiver: a cyclic line in the per-bin STFT occupancy at the sweep rate, "
            "confirmed by folding the STFT on the %.2f us period (one carrier moving continuously through "
            "%.1f MHz), not by energy alone",
            F["sweep_period_us"].get<double>(), F["sweep_span_mhz"].get<double>());
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission& e = em[0];
    const Json& F = e.features;
    static const char* keys[] = {"sweep_rate_hz", "sweep_period_us", "sweep_low_mhz", "sweep_high_mhz",
                                 "sweep_span_mhz", "beyond_capture", "shape", "in_band_fraction", "flybacks",
                                 "band_exits", "dominance_db", "fold_peak_db", "trajectory_us_mhz"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k)) fields[k] = F[k];
    auto num = [&](const char* k) { return F[k].get<double>(); };
    std::string beyond;
    for (const auto& b : F["beyond_capture"]) beyond += (beyond.empty() ? "" : " and ") + b.get<std::string>();
    const std::string extent = fmt("%.1f-%.1f MHz seen", num("sweep_low_mhz"), num("sweep_high_mhz")) +
                               (beyond.empty() ? "" : ", continuing " + beyond + " the captured band");
    // a continuous swept carrier in the radio-altimeter band is the band's own FMCW radar, not a jammer
    const auto bk = bands::band_keys(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6);
    const bool fmcw = std::find(bk.begin(), bk.end(), "radalt") != bk.end();
    const std::string label = fmcw ? "Swept carrier (FMCW-like)" : "Jammer (swept carrier)";
    const double conf = fmcw ? 0.6 : e.confidence;
    const std::string ev =
        fmt("one carrier swept %s every %.2f us (%.1f Hz), ", F["shape"].get<std::string>().c_str(),
            num("sweep_period_us"), num("sweep_rate_hz")) +
        extent + fmt("; in band %.0f%% of the period, %.1f dB over the swept span", num("in_band_fraction") * 100,
                     num("dominance_db")) +
        (fmcw ? "; in the radio-altimeter band: FMCW radar, not called a jammer" : "");
    return result("success", method(), steps(), fields,
                  "sweep measured (a jammer carries no payload; the repeating sweep is the verification)", label,
                  conf, F.value("modulation", std::string("swept CW")), ev,
                  std::make_pair(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6));
  }
};

// python BarrageDecoder: barrage-noise / noise-FM jammers on unexplained wide continuous energy boxes
class BarrageDecoder : public Decoder {
 public:
  std::string name() const override { return "jammer_barrage"; }
  std::string package() const override { return "jammer"; }
  std::string protocol() const override { return "barrage-noise / noise-FM jammer"; }
  std::string method() const {
    return "receivers.jammer.barrage: plateau, Gaussianity, no CP / cyclic structure, stationarity (research)";
  }
  std::vector<std::string> steps() const {
    return {"gate: an unexplained energy box >= 1 MHz wide, continuous over the block, >= 6 dB",
            "plateau (block FFT, 64 sub-bands of the inner 90 %): >= 90 % within +-3 dB, outside strips >= 6 dB "
            "lower, no narrow line",
            "statistics at 1.25 x the width: envelope k = E|x|^4 / E^2|x|^2, I / Q kurtosis",
            "no structure: Hann-tapered autocorrelation (cyclic prefix) < 0.025 over 1-1000 us, no cyclic line in "
            "|x|^2, x^2 or the instantaneous frequency, power per 1 ms within +-1 dB",
            "barrage: Gaussian + no structure + flat plateau, >= 10 MHz or edges on an allocation (covers_band); "
            "noise-FM: constant envelope, wide Gaussian-like instantaneous frequency, no structure"};
  }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return sig.bw() >= 1e6 && sig.snr_db >= 6.0 &&
           std::any_of(sig.bursts.begin(), sig.bursts.end(), [](const Box& b) { return b.kind == "continuous"; });
  }

  DecodeResult decode(const Signal& sig, DecodeContext& ctx) const override {
    const auto r = barrage::analyse(ctx.ch, sig.f0, sig.f1);
    if (!r || (*r)["class"].is_null()) {
      std::string why = "not a noise jammer";
      if (r)
        why += " (envelope k " + py_str((*r)["envelope_k"]) + ", cyclic-prefix correlation " +
               py_str((*r)["cp_corr"]) + ", plateau " + ((*r)["plateau_ok"].get<bool>() ? "flat" : "not flat") + ")";
      return result("not_possible", method(), steps(), Json::object(), why);
    }
    const Json& R = *r;
    std::string cbs;
    for (const auto& b : R["covers_band"]) cbs += (cbs.empty() ? "" : ", ") + b.get<std::string>();
    std::string label, ev, mod;
    double conf;
    if (R["class"] == "barrage") {
      label = "Jammer (barrage noise)-like";
      conf = cbs.empty() ? 0.6 : 0.75;
      ev = "flat Gaussian plateau " + py_str(R["plateau_bw_mhz"]) + " MHz (envelope k " + py_str(R["envelope_k"]) +
           ", I / Q kurtosis " + py_str(R["kurtosis_i"]) + " / " + py_str(R["kurtosis_q"]) + "), edges " +
           py_repr(R["edge_db"]) + " dB; no cyclic prefix (max " + py_str(R["cp_corr"]) +
           "), no cyclic line, stationary within " + py_str(R["stationarity_db"]) + " dB";
      mod = "noise";
    } else {
      label = "Jammer (noise-FM)-like";
      conf = 0.7;
      ev = "constant envelope (k " + py_str(R["envelope_k"]) + ") over " + py_str(R["plateau_bw_mhz"]) +
           " MHz, instantaneous frequency spread " + py_str(R["if_std_mhz"]) + " MHz with kurtosis " +
           py_str(R["if_kurtosis"]) + "; no cyclic line, stationary";
      mod = "noise-FM";
    }
    if (!cbs.empty()) ev += "; edges on the " + cbs + " allocation";
    return result("partial", method(), steps(), R, "a noise jammer carries no payload: the statistics are the evidence",
                  label, conf, mod, ev, std::make_pair(sig.f0, sig.f1));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) {
  out.push_back(std::make_unique<JammerDecoder>());
  out.push_back(std::make_unique<BarrageDecoder>());
}

}  // namespace si::rx::jammer
