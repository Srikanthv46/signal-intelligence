// Pipeline plugin: ILS glide slope / localizer and VOR (research; python: receivers/navaids/plugin.py).
#include <algorithm>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/navaids/detector.hpp"

namespace si::rx::navaids {

namespace {

class NavaidsDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "navaids"; }
  std::string package() const override { return "navaids"; }
  std::string protocol() const override { return "ILS / VOR navigation aids"; }
  std::string key() const override { return "navaids"; }
  std::set<std::string> bands() const override { return {"ils_glidepath", "vor_ils_loc", "vor"}; }
  std::string method() const override {
    return "receivers.navaids: carrier on the ICAO channel grid + tone fits (research)";
  }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"channels inside the capture: ILS glide path 329.15-335.00 MHz and localizer 108.10-111.95 MHz "
            "(ICAO pairing), VOR 108-117.95 MHz; a carrier line >= 15 dB over the +-35 kHz median, receiver ppm "
            "applied; a second line 4-32 kHz away = two-frequency system",
            "ILS: envelope of the strongest carrier (+-1.5 kHz, 6 kSps) over whole 1/30 s windows, least squares "
            "on [1, 90 s, 150 s Hz] with the common tone scale s searched over 0.97-1.03 -> m90, m150, DDM, SDM",
            "VOR: envelope (+-12.5 kHz, 32 kSps) over whole 30 Hz periods: 30 Hz AM by least squares; 9960 Hz "
            "subcarrier (+-700 Hz) from the envelope FFT, its instantaneous frequency fitted at 30 Hz -> "
            "deviation; radial = AM phase - FM phase",
            "validity: ILS tones >= 5 sigma and SDM of the kind (glide path 0.6-1.0, localizer 0.25-0.65); VOR "
            "subcarrier >= 10 dB, AM depth 0.15-0.45, deviation 300-700 Hz; else '-like'"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    const Json* cal = ctx.cache.get<Json>("ppm");
    if (!cal || cal->is_null() || cal->empty()) cal = ctx.cache.get<Json>("sweep_ppm");
    double ppm = 0.0;
    if (cal && cal->is_object() && cal->contains("ppm") && (*cal)["ppm"].is_number()) ppm = (*cal)["ppm"].get<double>();
    return detect(ctx.ch, ppm);
  }
  std::string group_key(const Emission& e) const override {
    return e.label + "|" + py_float_repr(e.features["channel_mhz"].get<double>());
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    size_t bi = 0;
    for (size_t i = 1; i < em.size(); ++i)
      if (em[i].confidence > em[bi].confidence) bi = i;
    const Emission& e = em[bi];
    const Json& F = e.features;
    const Json fields = F;
    const bool ok = F["verified_structure"].get<bool>();
    std::string kind, what, ev, mod;
    if (e.label == "vor") {
      kind = "VOR beacon";
      what = fmt("%.2f MHz", F["channel_mhz"].get<double>());
      ev = "carrier " + py_str(F["line_db"]) + " dB over its surroundings; 30 Hz AM " +
           fmt("%.0f", F["am_depth"].get<double>() * 100) + " %, 9960 Hz subcarrier +" + py_str(F["subcarrier_db"]) +
           " dB with 30 Hz FM +-" + py_str(F["deviation_hz"]) + " Hz; radial " + py_str(F["radial_deg"]) +
           " deg from the station";
      mod = "AM 30 Hz + 9960 Hz subcarrier FM (VOR)";
    } else {
      kind = e.label == "ils_gs" ? "ILS glide slope" : "ILS localizer";
      what = fmt("%.2f MHz, paired %s %.2f MHz", F["channel_mhz"].get<double>(), e.label == "ils_gs" ? "LOC" : "GS",
                 F["paired_mhz"].get<double>());
      ev = "carrier " + py_str(F["line_db"]) + " dB over its surroundings; 90 Hz " +
           fmt("%.1f", F["m90"].get<double>() * 100) + " % and 150 Hz " + fmt("%.1f", F["m150"].get<double>() * 100) +
           " % (tone scale " + py_str(F["tone_scale"]) + "); DDM " + fmt("%+.4f", F["ddm"].get<double>()) + " (" +
           F["side"].get<std::string>() + "), SDM " + fmt("%.3f", F["sdm"].get<double>());
      mod = "AM 90 / 150 Hz (ILS)";
    }
    if (F.contains("two_frequency_hz")) {
      what += ", two-frequency";
      ev += "; second carrier " + py_str(F["two_frequency_hz"]) + " Hz away (two-frequency system)";
    }
    const std::string label = kind + (ok ? "" : "-like") + " (" + what + ")";
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    return result("partial", method(), steps(), fields,
                  "navigation tones measured (an ILS / VOR carries no payload); ident (Morse 1020 Hz) needs a "
                  "dwell of several seconds",
                  label, e.confidence, mod, ev, std::make_pair(lo * 1e6, hi * 1e6));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<NavaidsDecoder>()); }

}  // namespace si::rx::navaids
