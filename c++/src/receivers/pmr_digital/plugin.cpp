// Pipeline plugin: digital land-mobile radio (TETRA, DMR, NXDN) (python: receivers/pmr_digital/plugin.py).
#include <algorithm>
#include <memory>
#include <string>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/pmr_digital/detector.hpp"

namespace si::rx::pmr_digital {

namespace {

const char* kTetraKeys[] = {"carrier_hz", "mcc", "mnc", "colour_code", "system_code", "sharing_mode",
                            "ts_reserved_frames", "u_plane_dtx", "frame18_extension", "cell_service_level",
                            "late_entry", "neighbour_cell_broadcast", "n_bsch_crc_ok", "n_sync_bursts",
                            "n_normal_training", "slot_consistent_hits", "bsch", "bsch_identities",
                            "symbol_rate_line", "phase_error_rad", "snr_db"};
const char* kDmrKeys[] = {"source", "colour_code", "colour_code_votes", "n_bursts", "sync_types", "voice_bursts",
                          "data_types", "crc_checked", "crc_ok", "csbk", "symbol_rate", "symbol_rate_line",
                          "deviation_hz_p99", "snr_db", "bursts"};
const char* kNxdnKeys[] = {"symbol_rate", "channel_bw_khz", "n_frame_sync", "frame_consistent_hits",
                           "symbol_rate_line", "deviation_hz_p99", "snr_db"};

template <size_t N>
Json pick(const Json& F, const char* const (&keys)[N]) {
  Json out = Json::object();
  for (const char* k : keys)
    if (F.contains(k) && !F[k].is_null()) out[k] = F[k];
  return out;
}

std::string decoder_of(const Signal& s) {
  return s.features.contains("decoder") && s.features["decoder"].is_string() ? s.features["decoder"].get<std::string>()
                                                                             : std::string();
}

class PmrDigitalDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "pmr_digital"; }
  std::string package() const override { return "pmr_digital"; }
  std::string protocol() const override { return "digital land-mobile radio (TETRA, DMR, NXDN)"; }
  std::string key() const override { return "pmr_digital"; }
  std::string method() const override {
    return "receivers.pmr_digital: TETRA pi/4-DQPSK + BSCH; 4FSK discriminator + DMR sync/Golay/BPTC, NXDN sync";
  }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"candidates: narrow carriers (<= 30 kHz, >= 8 dB) in 380-400/410-430/450-470/870-876/915-921 MHz "
            "(TETRA) and 136-174/300-520/806-870 MHz (4FSK), 1 kHz bins",
            "TETRA: 25 kHz channel at 144 kSps, 18 ksym/s envelope line; pi/4-DQPSK (RRC 0.35) with the "
            "4th-power carrier estimate and its quarter-turn ambiguity resolved by the training sequences; "
            "synchronisation (38-bit) and normal (22-bit) training sequences on the 510-bit slot grid; "
            "BSCH: descramble (e=0), de-interleave (120, 11), RCPC 2/3 -> rate-1/4 K=5 soft Viterbi, CRC-CCITT "
            "-> SYNC PDU (MCC, MNC, colour code, TN/FN/MN)",
            "4FSK: 12.5 kHz channel at 48 kSps, FM discriminator, 2400/4800 sym/s line; DMR: SYNC (9 patterns) "
            "-> per-burst RRC timing/levels from the sync -> slot type Golay(20,8) (colour code, data type) "
            "-> BPTC(196,96) Hamming(13,9)/(15,11) -> CRC-CCITT with the data-type mask (CSBK opcode when valid); "
            "NXDN: 20-bit frame sync on the 192-symbol frame"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    dvec excl;
    if (const Json* ln = ctx.cache.get<Json>("receiver_lines"))
      for (const auto& l : *ln) excl.push_back(l["freq_hz"].get<double>());
    std::optional<std::pair<double, double>> band;
    if (ctx.det) band = ctx.det->usable_hz;
    return detect(ctx.ch, band, excl);
  }

  bool in_range(double lo, double hi) const override {
    std::vector<std::pair<double, double>> r = tetra_ranges();
    r.insert(r.end(), fsk_ranges().begin(), fsk_ranges().end());
    return in_ranges(lo, hi, r);
  }
  std::string group_key(const Emission& e) const override {
    return e.family + "|" + py_float_repr(py_round(e.center_mhz, 3));
  }
  bool keep(const Emission&) const override { return true; }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    size_t bi = 0;
    for (size_t i = 1; i < em.size(); ++i)
      if (em[i].confidence > em[bi].confidence) bi = i;
    const Emission& e = em[bi];
    const Json& F = e.features;
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    const auto span = std::make_pair(lo * 1e6, hi * 1e6);
    const std::string where = fmt("%.4f MHz", e.center_mhz);
    if (e.family == "tetra") {
      const Json fields = pick(F, kTetraKeys);
      const std::string line = py_str(F["symbol_rate_line"]);
      if (e.label == "tetra_downlink_bsch") {
        return result("success", method(), steps(), fields,
                      "BSCH (SYNC PDU) CRC-verified; MAC/MLE broadcast (SYSINFO) not decoded",
                      "TETRA downlink MCC " + py_str(F["mcc"]) + " MNC " + py_str(F["mnc"]) + " colour code " +
                          py_str(F["colour_code"]) + " (" + where + ")",
                      e.confidence, std::string("pi/4-DQPSK 18 ksym/s (TETRA)"),
                      py_str(F["n_bsch_crc_ok"]) + " of " + py_str(F["n_sync_bursts"]) +
                          " synchronisation bursts' BSCH pass the CRC; training sequences on one slot phase (" +
                          py_str(F["slot_consistent_hits"]) + "); 18 ksym/s line x" + line,
                      span);
      }
      if (e.label == "tetra_downlink") {
        return result("partial", method(), steps(), fields,
                      "TETRA burst structure (training sequences) found; no synchronisation burst decoded",
                      "TETRA downlink (" + where + ")", e.confidence, std::string("pi/4-DQPSK 18 ksym/s (TETRA)"),
                      py_str(F["slot_consistent_hits"]) + " training sequences on one 510-bit slot phase; " +
                          "18 ksym/s line x" + line,
                      span);
      }
      return result("partial", method(), steps(), fields,
                    "18 ksym/s pi/4-DQPSK-like carrier, no training sequence in this block",
                    "TETRA-like carrier (" + where + ")", e.confidence, std::string("pi/4-DQPSK-like 18 ksym/s"),
                    "18 ksym/s line of the matched-filtered envelope x" + line, span);
    }
    if (e.family == "dmr") {
      const Json fields = pick(F, kDmrKeys);
      const bool has_cc = F.contains("colour_code") && !F["colour_code"].is_null();
      const std::string cc = has_cc ? ", colour code " + py_str(F["colour_code"]) : "";
      std::string what;
      if (F.contains("voice_bursts") && F["voice_bursts"].get<int>() != 0) {
        what = "voice call";
      } else if (F.contains("data_types") && !F["data_types"].empty()) {
        std::string top;
        int topn = -1;
        for (const auto& [k, v] : F["data_types"].items())  // python max(dts, key=dts.get): first maximum
          if (v.get<int>() > topn) {
            topn = v.get<int>();
            top = k;
          }
        if (top == "CSBK") what = "control channel (CSBK)";
      }
      const std::string lab =
          "DMR " + F["source"].get<std::string>() + cc + (what.empty() ? "" : ", " + what) + " (" + where + ")";
      const bool ok = F.contains("csbk") && !F["csbk"].empty();
      std::string reason;
      if (ok) {
        reason = "CSBK CRC-verified";
      } else {
        reason = "SYNC, slot type (Golay) and BPTC (Hamming) verified; ";
        if (F.contains("crc_checked") && F["crc_checked"].get<int64_t>() != 0) {
          const int64_t n = F["crc_checked"].get<int64_t>();
          const int64_t okn = F.contains("crc_ok") ? F["crc_ok"].get<int64_t>() : 0;
          reason += "CSBK/header CRC failed " + std::to_string(n - okn) + " of " + std::to_string(n) +
                    " (non-standard or vendor control channel?)";
        } else {
          reason += "no CRC-protected block in this capture";
        }
      }
      std::string ev = py_str(F["n_bursts"]) + " DMR SYNC word(s) " + py_repr(F["sync_types"]);
      if (has_cc)
        ev += "; colour code " + py_str(F["colour_code"]) + " in " + py_str(F["colour_code_votes"]) +
              " slot types (Golay)";
      return result(ok ? "success" : "partial", method(), steps(), fields, reason, lab, e.confidence,
                    std::string("4FSK 4800 sym/s (DMR)"), ev, span);
    }
    const Json fields = pick(F, kNxdnKeys);
    const double sr = F["symbol_rate"].get<double>();
    return result("partial", method(), steps(), fields, "NXDN frame sync verified; LICH/SACCH not decoded",
                  fmt("NXDN %g sym/s (", sr) + where + ")", e.confidence, fmt("4FSK %g sym/s (NXDN)", sr),
                  py_str(F["n_frame_sync"]) + " NXDN frame sync word(s), " + py_str(F["frame_consistent_hits"]) +
                      " on the 192-symbol frame grid; symbol-rate line x" + py_str(F["symbol_rate_line"]),
                  span);
  }

  // python arbitrate: a carrier verified here is digital; drop voice_fm's analog reading of it
  std::vector<Signal> arbitrate(std::vector<Signal> signals, DecodeContext&) const override {
    std::vector<double> mine;
    for (const Signal& s : signals)
      if (decoder_of(s) == name()) mine.push_back(s.fc());
    if (mine.empty()) return signals;
    std::vector<Signal> out;
    for (auto& s : signals) {
      bool drop = false;
      if (decoder_of(s) == "voice_fm")
        for (double m : mine) drop = drop || std::fabs(s.fc() - m) < 4e3;
      if (!drop) out.push_back(std::move(s));
    }
    return out;
  }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<PmrDigitalDecoder>()); }

}  // namespace si::rx::pmr_digital
