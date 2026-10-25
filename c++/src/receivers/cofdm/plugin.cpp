// Pipeline plugin: DVB-T-framed COFDM links (python: receivers/cofdm/plugin.py; research capability).
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/cofdm/detector.hpp"

namespace si::rx::cofdm {

namespace {

class CofdmDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "cofdm"; }
  std::string package() const override { return "cofdm"; }
  std::string protocol() const override { return "DVB-T COFDM video link"; }
  std::string key() const override { return "cofdm"; }
  std::string method() const override { return "receivers.cofdm (research receiver: ETSI EN 300 744, synthetic validation only)"; }
  std::vector<std::string> steps() const override {
    return {"candidates: Welch spectrum runs 0.951 x B wide (B = 8 / 7 / 6 / 5 / 4 / 2.5 / 2 MHz) above floor + 2 dB",
            "numerology: slice at 8 B / 7 MHz, CP metric folded on N (1 + g) for 2K / 8K x GI 1/4..1/32",
            "carrier: fractional offset from the CP phase, integer offset from the 4/3-boosted continual pilots",
            "TPS: DBPSK across symbols on the 17 / 68 TPS carriers, sync word, BCH(67,53) -> mode, GI, "
            "constellation, hierarchy, code rates, frame number, cell-ID byte"};
  }
  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    const double rate = std::min(ctx.cap.fs, ctx.cap.bandwidth_hz.value_or(ctx.cap.fs));
    if (rate < 2.5e6) return {};
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc);
  }
  bool in_range(double, double) const override { return true; }
  std::string group_key(const Emission& e) const override { return py_float_repr(e.center_mhz); }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission& e = em[0];
    const Json& F = e.features;
    if (F.contains("p1") && F["p1"].is_object()) {
      // DVB-T2 P1 preamble (python summarise; research, PLAN D5)
      const Json& p1 = F["p1"];
      Json fields = Json::object();
      for (const char* k : {"standard", "channel_bw_mhz", "occupied_mhz", "elementary_rate_mhz", "carrier_mhz", "p1_count",
                            "p1_metric", "p1_start_us", "p1"})
        if (F.contains(k)) fields[k] = F[k];
      const std::string B = fmt("%g MHz", F["channel_bw_mhz"].get<double>());
      const bool like = e.label != "dvbt2";
      const std::string label = std::string("COFDM video link") + (like ? "-like" : "") + " (DVB-T2 " +
                                p1["preamble"].get<std::string>() + ", " + B + ", " + p1["fft_gi"].get<std::string>() + ")";
      const std::string ev = "P1 C-A-B correlation " + py_str(F["p1_metric"]) + " (" + py_str(F["p1_count"]) +
                             " in the block); S1 " + p1["s1"].get<std::string>() + " / S2 " + p1["s2"].get<std::string>() +
                             " (pattern margins " + py_str(p1["s1_margin"]) + " / " + py_str(p1["s2_margin"]) + ")";
      return result(like ? "partial" : "success", method(), steps(), fields,
                    "P1 signalling decoded (S1 / S2 by pattern correlation); L1 / payload not decoded", label, e.confidence,
                    "COFDM (DVB-T2)", ev, std::make_pair(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6));
    }
    const Json tps = F.contains("tps") && F["tps"].is_object() ? F["tps"] : Json::object();
    static const char* keys[] = {"standard", "channel_bw_mhz", "occupied_mhz", "mode", "guard_interval", "fft_size",
                                 "elementary_rate_mhz", "scs_hz", "tu_us", "symbol_us", "carrier_mhz",
                                 "integer_offset_carriers", "cp_peak", "cp_contrast", "symbols", "continual_pilot_db",
                                 "continual_pilot_margin", "tps", "tps_note"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k)) fields[k] = F[k];
    const std::string B = fmt("%g MHz", F["channel_bw_mhz"].get<double>());
    const std::string mode = F["mode"].get<std::string>(), gi = F["guard_interval"].get<std::string>();
    std::string label, status, reason;
    const bool ok = tps.contains("bch_ok") && tps["bch_ok"].get<bool>();
    if (ok) {
      label = "COFDM video link (DVB-T " + mode + ", " + B + ", " + tps["constellation"].get<std::string>() + ", HP " +
              tps["code_rate_hp"].get<std::string>() + ", GI " + tps["guard_interval"].get<std::string>() + ")";
      status = "success";
      reason = "TPS decoded (BCH(67,53) passed); the video payload (MPEG-TS, often AES) is not decoded";
    } else {
      label = "COFDM video link-like (DVB-T " + mode + " numerology, " + B + ", GI " + gi + ")";
      status = "partial";
      reason = F.contains("tps_note") && F["tps_note"].is_string() ? F["tps_note"].get<std::string>() : "TPS not BCH-verified";
    }
    std::string ev = "CP fold peak " + py_str(F["cp_peak"]) + " (" + py_str(F["cp_contrast"]) + "x median) at " + mode + " GI " +
                     gi + ", " + py_str(F["symbols"]) + " symbols; continual pilots +" + py_str(F["continual_pilot_db"]) +
                     " dB at offset " + py_str(F["integer_offset_carriers"]);
    if (!tps.empty()) ev += "; TPS frame " + py_str(tps["frame_number"]) + ", BCH " + (ok ? "ok" : "failed");
    const std::string con = tps.contains("constellation") ? tps["constellation"].get<std::string>() : "QPSK/16/64-QAM";
    return result(status, method(), steps(), fields, reason, label, e.confidence, "COFDM " + mode + " (" + con + ")", ev,
                  std::make_pair(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<CofdmDecoder>()); }

}  // namespace si::rx::cofdm
