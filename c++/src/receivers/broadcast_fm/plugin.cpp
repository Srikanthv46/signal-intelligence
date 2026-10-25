// Pipeline plugin: FM broadcast stations and RDS (research; python: receivers/broadcast_fm/plugin.py).
#include <algorithm>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/broadcast_fm/detector.hpp"

namespace si::rx::broadcast_fm {

namespace {

class BroadcastFmDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "broadcast_fm"; }
  std::string package() const override { return "broadcast_fm"; }
  std::string protocol() const override { return "FM broadcast + RDS"; }
  std::string key() const override { return "broadcast_fm"; }
  std::set<std::string> bands() const override { return {"fm_bcast"}; }
  std::string method() const override {
    return "receivers.broadcast_fm: FM discriminator, stereo pilot, RDS block code (research)";
  }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"100 kHz raster 87.6-107.9 MHz: channel power over +-100 kHz >= 10 dB over the band floor (local "
            "maxima on the raster)",
            "+-120 kHz at 240 kSps, FM discriminator -> MPX: peak deviation 15-110 kHz; stereo when the 19 kHz "
            "pilot line is >= 10 dB over its surroundings",
            "RDS (blocks >= 45 ms): MPX shifted by 3 x the pilot, +-2.4 kHz, BPSK phase from the squared signal, "
            "biphase matched filter at 1187.5 bit/s, differential decoding, block sync on the offset words "
            "(error-free (26,16) blocks only) -> PI, PTY, TP, programme name, RadioText, clock time"};
  }
  std::vector<Emission> run_detector(DecodeContext& ctx) const override { return detect(ctx.ch); }
  std::string group_key(const Emission& e) const override {
    return py_float_repr(e.features["channel_mhz"].get<double>());
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    size_t bi = 0;
    for (size_t i = 1; i < em.size(); ++i)
      if (em[i].confidence > em[bi].confidence) bi = i;
    const Emission& e = em[bi];
    const Json& F = e.features;
    const Json rds = F.contains("rds") ? F["rds"] : Json::object();
    std::string parts = fmt("%.1f MHz", F["channel_mhz"].get<double>());
    if (F["stereo"].get<bool>()) parts += ", stereo";
    if (rds.contains("ps")) {
      std::string ps = rds["ps"].get<std::string>();
      std::replace(ps.begin(), ps.end(), '_', ' ');
      const auto a = ps.find_first_not_of(" \t\n\r\x0b\x0c"), b = ps.find_last_not_of(" \t\n\r\x0b\x0c");
      parts += ", RDS '" + (a == std::string::npos ? std::string() : ps.substr(a, b - a + 1)) + "'";
    } else if (rds.contains("pi")) {
      parts += ", RDS PI " + rds["pi"].get<std::string>();
    }
    const bool ok = rds.contains("pi");
    std::string ev = "FM carrier " + py_str(F["level_db"]) + " dB over the band floor, deviation " +
                     py_str(F["deviation_khz"]) + " kHz" +
                     (F["stereo"].get<bool>() ? ", 19 kHz pilot +" + py_str(F["pilot_db"]) + " dB" : ", no stereo pilot");
    if (ok) {
      ev += "; RDS: " + py_str(rds["blocks_ok"]) + " error-free block(s), PI " + rds["pi"].get<std::string>();
      if (rds.contains("pty_name")) ev += ", PTY " + rds["pty_name"].get<std::string>();
    }
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    return result(ok ? "success" : "partial", method(), steps(), F,
                  ok ? "RDS blocks verified by their (26,16) checkwords"
                     : "FM broadcast measured; no RDS in this block (RDS needs >= 45 ms, many stations send none)",
                  "FM broadcast (" + parts + ")", e.confidence,
                  F["stereo"].get<bool>() ? "wideband FM (+-75 kHz), stereo MPX" : "wideband FM (+-75 kHz)", ev,
                  std::make_pair(lo * 1e6, hi * 1e6));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<BroadcastFmDecoder>()); }

}  // namespace si::rx::broadcast_fm
