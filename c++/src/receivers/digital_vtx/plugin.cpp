// Pipeline plugin: Herelink and unknown-OFDM digital video links (python: receivers/digital_vtx/plugin.py).
#include <map>
#include <memory>
#include <regex>

#include "core/bands_india.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/digital_vtx/catalog.hpp"
#include "receivers/digital_vtx/detector.hpp"

namespace si::rx::digital_vtx {
namespace {

constexpr double kHerelinkMinHz = 2.400e9;  // below: LTE band 40 (same PHY), left to the LTE decoder

// python NUMEROLOGY_NAMES
std::string numerology_name(const std::string& n) {
  static const std::map<std::string, std::string> m = {
      {"lte_like", "LTE-like: DJI OcuSync / O3, CubePilot, Autel"}, {"dji_like", "DJI OcuSync-like"},
      {"autel_skylink_like", "Autel SkyLink-like"},                 {"herelink", "Herelink-like"},
      {"dji_o4_like", "DJI O4-like"},                               {"lightbridge_like", "DJI Lightbridge / WiMAX-like"},
      {"walksnail_like", "Walksnail / WiMAX-like"},                 {"pmp450_like", "Cambium PMP 450-like"},
      {"hdzero_like", "HDZero-like"}};
  const auto it = m.find(n);
  return it == m.end() ? n : it->second;
}

bool starts_with(const std::string& s, const char* p) { return s.rfind(p, 0) == 0; }

bool truthy_num(const Json& v) { return v.is_number() && v.get<double>() != 0; }

class DigitalVtxDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "digital_vtx"; }
  std::string package() const override { return "digital_vtx"; }
  std::string protocol() const override { return "Herelink / HDZero / unknown OFDM video link"; }
  std::string key() const override { return "digital_vtx"; }
  std::string method() const override { return "receivers.digital_vtx (research receiver)"; }
  std::vector<std::string> steps() const override {
    return {"Herelink: LTE PSS (ZC 25/29/34) fold on the 5 ms half-frame over the band, carrier search, "
            "SSS decode -> PCI, CP-phase CFO, subframe activity -> control / video / uplink",
            "HDZero (27 MHz mode): on each HDZero channel in the capture, 1024-FFT at 27 MSps slice, CP-fold symbol "
            "timing and fractional CFO, continual-pilot coherence vs the data carriers beside them over the integer "
            "carrier offsets -> channel, ppm, per-frame preamble -> video frame rate",
            "any other CP-OFDM link (Tu >= 10 us): coherent lag-Tu autocorrelation scan (time-AGC'd), symbol-rate "
            "line / per-window fold confirmation, CP cross-spectrum band, bursts on local symbol grids, TDD frame "
            "completion, band and burst extent from in-band power seeded by the CP evidence; numerology table -> "
            "drone family / infrastructure / unknown",
            "hopping narrowband CP-OFDM (RC uplinks): windowed CP-fold scan at the capture's and the drone families' "
            "numerologies, envelope edges, on-vs-off band fit, CP-fold verification, frame-cadence tracking"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, false);
  }

  bool in_range(double lo, double hi) const override {
    return in_ranges(lo, hi, bands::digital_link_ranges());
  }

  bool keep(const Emission& e) const override {
    if (e.family == "herelink") return e.center_mhz * 1e6 >= kHerelinkMinHz && e.confidence >= 0.5;
    if (e.family == "hdzero") return true;
    if (e.family != "digital_vtx") return false;
    // inside a licensed cellular band 15/30 kHz OFDM is the operator's carrier
    static const std::regex cellular("^[bn]\\d");
    for (const auto& k : bands::band_keys(e.center_mhz * 1e6, e.center_mhz * 1e6))
      if (std::regex_search(k, cellular)) return false;
    // 802.11-family numerology at a scaled clock is the Wi-Fi receivers'
    return !(e.features.contains("role") && e.features["role"] == "wifi");
  }

  // python: herelink (family, "", round(center, 0)); blind links (family, label, carrier_mhz)
  std::string group_key(const Emission& e) const override {
    if (e.family == "herelink" || e.family == "hdzero") return e.family + "||" + py_float_repr(np::round(e.center_mhz, 0));
    if (starts_with(e.label, "ofdm_hop"))   // one signal per hop channel
      return e.family + "|" + e.label + "|" + py_repr(e.features.contains("channel_mhz") ? e.features["channel_mhz"] : Json(nullptr));
    return e.family + "|" + e.label + "|" + py_repr(e.features.contains("carrier_mhz") ? e.features["carrier_mhz"] : Json(nullptr));
  }

  // python arbitrate: a blind OFDM link that a DJI video link or Herelink of the same subcarrier
  // spacing covers over more than half of its band is that receiver's signal (recorded on the host)
  std::vector<Signal> arbitrate(std::vector<Signal> signals, DecodeContext&) const override {
    auto is_gen = [&](const Signal& s) {
      return s.features.value("decoder", std::string()) == name() && !starts_with(s.label, "Herelink") &&
             !starts_with(s.label, "HDZero");
    };
    bool any = false;
    for (const auto& s : signals) any = any || is_gen(s);
    if (!any) return signals;
    std::vector<std::pair<size_t, dvec>> hosts;
    for (size_t i = 0; i < signals.size(); ++i) {
      const Signal& s = signals[i];
      const std::string dec = s.features.value("decoder", std::string());
      dvec scs;
      if (dec == "dji") {
        for (const auto& e : s.emissions)
          if (e.label != "uplink" && e.features.contains("scs_khz") && e.features["scs_khz"].is_number())
            scs.push_back(e.features["scs_khz"].get<double>());
      } else if (dec == name() && starts_with(s.label, "Herelink")) {
        scs.push_back(15.0);
      } else if (dec == name() && starts_with(s.label, "HDZero")) {
        scs.push_back(26.367);
      } else {
        continue;
      }
      hosts.emplace_back(i, scs);
    }
    // a blind hop that a DJI RC uplink hop covers (same time, same band) is the DJI receiver's; hops
    // that mostly sit inside a DJI video link's bursts and band are slices of that video
    std::vector<size_t> ups, vids;
    for (size_t i = 0; i < signals.size(); ++i) {
      if (signals[i].features.value("decoder", std::string()) != "dji") continue;
      if (signals[i].label.find("uplink") != std::string::npos) ups.push_back(i);
      if (signals[i].label.find("video") != std::string::npos) vids.push_back(i);
    }
    auto record = [&](Signal& h, const Signal& s, const char* note) {
      Json rec = Json::object();
      rec["decoder"] = name();
      rec["label"] = s.label;
      rec["f_low_hz"] = py_round_int(s.f0);
      rec["f_high_hz"] = py_round_int(s.f1);
      rec["note"] = note;
      if (!h.detection.contains("explained_receiver_detections")) h.detection["explained_receiver_detections"] = Json::array();
      h.detection["explained_receiver_detections"].push_back(rec);
    };
    std::vector<uint8_t> drop(signals.size(), 0);
    for (size_t i = 0; i < signals.size(); ++i) {
      const Signal& s = signals[i];
      if (!is_gen(s)) continue;
      const Json& F = s.decoded.fields;
      if (F.contains("hop_us") && !F["hop_us"].is_null()) {
        bool done = false;
        for (size_t u : ups) {
          const Signal& h = signals[u];
          if (!(std::min(h.f1, s.f1) - std::max(h.f0, s.f0) > 0.5 * std::min(h.bw(), s.bw()))) continue;
          bool tov = false;
          for (const auto& a : h.emissions)
            for (const auto& b : s.emissions) tov = tov || std::min(a.end_us, b.end_us) - std::max(a.start_us, b.start_us) > 0;
          if (!tov) continue;
          record(signals[u], s, "the same hops found blind (numerology, band, timing): corroborates it");
          drop[i] = done = true;
          break;
        }
        if (done) continue;
        for (size_t v : vids) {
          const Signal& h = signals[v];
          if (!(std::min(h.f1, s.f1) - std::max(h.f0, s.f0) > 0.5 * s.bw())) continue;
          int64_t inside = 0;
          for (const auto& b : s.emissions) {
            bool any = false;
            for (const auto& a : h.emissions)
              any = any || std::min(b.end_us, a.end_us) - std::max(b.start_us, a.start_us) > 0.5 * (b.end_us - b.start_us);
            inside += any;
          }
          if (!(2 * inside > static_cast<int64_t>(s.emissions.size()))) continue;
          record(signals[v], s, "narrowband slices of this link's bursts, found blind: part of it");
          drop[i] = true;
          break;
        }
        continue;
      }
      const double scs = F.contains("scs_khz") && F["scs_khz"].is_number() ? F["scs_khz"].get<double>() : 0.0;
      for (const auto& [hi, hs] : hosts) {
        Signal& h = signals[hi];
        bool same = false;
        for (double v : hs) same = same || std::fabs(scs / v - 1) < 0.03;
        if (!(std::min(h.f1, s.f1) - std::max(h.f0, s.f0) > 0.5 * std::min(h.bw(), s.bw()) && same)) continue;
        Json rec = Json::object();
        rec["decoder"] = name();
        rec["label"] = s.label;
        rec["f_low_hz"] = py_round_int(s.f0);
        rec["f_high_hz"] = py_round_int(s.f1);
        rec["note"] = "the same OFDM link found blind (numerology, band, bursts): corroborates it";
        if (!h.detection.contains("explained_receiver_detections"))
          h.detection["explained_receiver_detections"] = Json::array();
        h.detection["explained_receiver_detections"].push_back(rec);
        drop[i] = 1;
        break;
      }
    }
    std::vector<Signal> keep;
    for (size_t i = 0; i < signals.size(); ++i)
      if (!drop[i]) keep.push_back(std::move(signals[i]));
    return keep;
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission* best = &em[0];
    for (const auto& x : em)
      if (x.confidence > best->confidence) best = &x;
    const Emission& e = *best;
    const Json& F = e.features;
    Json labels = Json::object();
    for (const auto& x : em) labels[x.label] = labels.value(x.label, 0) + 1;
    static const char* keys[] = {"pci",          "nid1",       "carrier_mhz",    "cfo_hz",          "cfo_ppm",
                                 "bw_mhz",       "scs_khz",    "tu_us",          "cp_us",           "symbol_us",
                                 "numerology",   "frame_period_ms", "burst_us",  "clock_ppm",       "pss_instances",
                                 "role",         "matches",    "cp_fraction",    "duty_cycle",      "symbol_grids",
                                 "line_z",       "strong_windows", "link_snr_db", "channel_mhz",   "hop_bw_mhz",
                                 "hop_us",       "n_hops",     "n_channels",     "hop_period_ms",   "paired_link_mhz", "tx_antennas",
                                 "system_candidates", "channel",  "bandwidth_mode", "fft_size", "sample_rate_msps",
                                 "used_subcarriers", "continual_pilots", "pilot_coherence", "data_coherence",
                                 "cp_fold_contrast", "n_symbols", "preamble_ms", "frame_symbols", "frame_rate_fps",
                                 "band_context"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k) && !F[k].is_null()) fields[k] = F[k];
    fields["n_bursts"] = em.size();
    fields["labels"] = labels;
    if (e.label.rfind("ofdm_hop", 0) == 0) {   // model level from the hop timing (research, PLAN M6)
      auto opt = [&](const char* k) {
        return F.contains(k) && F[k].is_number() ? std::optional<double>(F[k].get<double>()) : std::nullopt;
      };
      Json mb = catalog::hop_model(opt("hop_bw_mhz"), opt("hop_us"), opt("hop_period_ms"));
      if (!mb.is_null()) fields["model"] = std::move(mb);
    }
    if (F.contains("system_candidates") && F["system_candidates"].is_array() && !F["system_candidates"].empty()) {
      // model level (research, PLAN M3 / M4)
      auto opt = [&](const char* k) {
        return F.contains(k) && F[k].is_number() ? std::optional<double>(F[k].get<double>()) : std::nullopt;
      };
      std::optional<int> na;
      if (F.contains("tx_antennas") && F["tx_antennas"].is_object() && F["tx_antennas"].contains("tx_antennas") &&
          F["tx_antennas"]["tx_antennas"].is_number())
        na = F["tx_antennas"]["tx_antennas"].get<int>();
      Json mb = catalog::model_candidates(F["system_candidates"], opt("carrier_mhz"), opt("bw_mhz"), na);
      if (!mb.is_null()) fields["model"] = std::move(mb);
    }
    auto get = [&](const char* k) { return F.contains(k) ? F[k] : Json(nullptr); };
    std::string status = "partial", reason, label, mod, ev;
    if (e.family == "hdzero") {
      fields["model"] = catalog::hdzero_model(py_str(get("bandwidth_mode")));
      reason = "HDZero PHY synchronised (symbol timing, carrier, pilots, video frame); the video samples are "
               "mapped straight onto the subcarriers (JSCC) with no published format: not decoded";
      label = "HDZero digital FPV video link";
      mod = "OFDM " + py_str(get("scs_khz")) + " kHz SCS (1024 FFT at 27 MSps), CP " + py_str(get("cp_us")) + " us, " +
            py_str(get("used_subcarriers")) + " subcarriers";
      const Json fps = get("frame_rate_fps"), pre = get("preamble_ms");
      ev = "HDZero continual pilots coherent on channel " + py_str(get("channel")) + " (" + py_str(get("carrier_mhz")) +
           " MHz, " + py_str(get("cfo_ppm")) + " ppm): pilot coherence " + py_str(get("pilot_coherence")) + " vs " +
           py_str(get("data_coherence")) + " on the data carriers beside them over " + py_str(get("n_symbols")) +
           " symbols, CP fold contrast " + py_str(get("cp_fold_contrast")) + ", " + py_str(get("bw_mhz")) + " MHz, SNR " +
           py_str(get("link_snr_db")) + " dB";
      if (truthy_num(fps))
        ev += "; per-frame preamble every " + py_str(get("frame_period_ms")) + " ms (" + py_str(fps) + " fps)";
      else if (pre.is_array() && !pre.empty())
        ev += "; preamble at " + py_str(pre) + " ms";
    } else if (e.family == "herelink") {
      const bool has_pci = F.contains("pci") && !F["pci"].is_null();
      reason = std::string("LTE-like PHY synchronised (PSS") + (has_pci ? "+SSS -> PCI" : "") +
               "); video payload is proprietary/encrypted";
      label = "Herelink (CubePilot) drone video/control link";
      mod = "OFDM 15 kHz SCS (LTE Rel-8 DL PHY)";
      const Json cmj = get("carrier_mhz");
      const double cm = truthy_num(cmj) ? cmj.get<double>() : e.center_mhz;
      const std::string ism = cm < 3000 ? "2.4 GHz" : cm >= 5650 ? "5.8 GHz" : "5 GHz";
      ev = "LTE PSS at " + py_str(cmj) + " MHz in the " + ism + " band with Herelink TDD subframe rhythm (" +
           py_repr(labels) + ")" + (has_pci ? ", PCI " + py_str(F["pci"]) : "");
    } else if (starts_with(e.label, "ofdm_hop")) {
      const Json sj = get("scs_khz");
      const double scs = sj.is_number() ? sj.get<double>() : 0.0;
      std::string sc = fmt("%.2f", scs);
      while (!sc.empty() && sc.back() == '0') sc.pop_back();
      if (!sc.empty() && sc.back() == '.') sc.pop_back();
      const std::string role = F.value("role", std::string()), num = F.value("numerology", std::string());
      const Json nch = get("n_channels");
      const bool hopping = nch.is_number() && nch.get<double>() >= 2;
      if (role == "drone")
        label = std::string("OFDM ") + (hopping ? "hopping " : "") + "RC uplink-like (" + numerology_name(num) +
                " numerology, " + sc + " kHz SCS)";
      else if (role == "infrastructure")
        label = "OFDM fixed-wireless uplink-like (" + numerology_name(num) + " numerology, " + sc + " kHz SCS)";
      else
        label = std::string("unknown OFDM ") + (hopping ? "hopping" : "narrowband burst") + " link (" + sc + " kHz SCS)";
      reason = std::string("narrowband CP-OFDM bursts (") + (hopping ? "hopping between channels" : "on one channel") +
               "), found and measured blindly; no protocol receiver identified them";
      mod = "OFDM " + sc + " kHz SCS, CP " + py_str(get("cp_us")) + " us, ~" + py_str(get("hop_bw_mhz")) + " MHz hops";
      const Json pl = get("paired_link_mhz");
      ev = "blind CP-OFDM hops: Tu " + py_str(get("tu_us")) + " us, CP " + py_str(get("cp_us")) + " us, " +
           std::to_string(em.size()) + " hop(s) on this channel (" + py_str(get("channel_mhz")) + " MHz), " +
           py_str(get("n_hops")) + " hops on " + py_str(get("n_channels")) + " channel(s) in the capture, " +
           py_str(get("hop_us")) + " us x " + py_str(get("hop_bw_mhz")) + " MHz each, hop period " +
           py_str(get("hop_period_ms")) + " ms" +
           (pl.is_array() && !pl.empty() ? "; same numerology as the OFDM link at " + py_str(pl) + " MHz" : "");
    } else {
      const Json s = get("scs_khz");
      const double scs = s.is_number() ? s.get<double>() : 0.0;
      std::string sc = fmt("%.2f", scs);  // f"{scs:.2f}".rstrip("0").rstrip(".")
      while (!sc.empty() && sc.back() == '0') sc.pop_back();
      if (!sc.empty() && sc.back() == '.') sc.pop_back();
      const std::string role = F.value("role", std::string()), num = F.value("numerology", std::string());
      if (role == "drone")
        label = "OFDM drone video link-like (" + numerology_name(num) + " numerology, " + sc + " kHz SCS)";
      else if (role == "infrastructure")
        label = "OFDM fixed-wireless link-like (" + numerology_name(num) + " numerology, " + sc + " kHz SCS)";
      else
        label = std::string("unknown OFDM ") + (truthy_num(get("frame_period_ms")) ? "TDD " : "") + "link (" + sc +
                " kHz SCS)";
      const std::string band_ctx = F.contains("band_context") && F["band_context"].is_string()
                                       ? F["band_context"].get<std::string>() : std::string();
      if (!band_ctx.empty()) label += " in the " + band_ctx.substr(0, band_ctx.find(" ("));  // PLAN G8
      reason = role == "drone" || role == "infrastructure"
                   ? "CP-OFDM link found and measured blindly; the numerology resembles a known system but no "
                     "protocol receiver identified it"
                   : "CP-OFDM link found and measured blindly; numerology unknown, no protocol to decode";
      mod = "OFDM " + sc + " kHz SCS, CP " + py_str(get("cp_us")) + " us";
      const Json tx = F.contains("tx_antennas") && F["tx_antennas"].is_object() ? F["tx_antennas"] : Json::object();
      if (tx.contains("tx_antennas") && tx["tx_antennas"].is_number() && tx["tx_antennas"].get<int>() == 2) {
        const std::string scheme = tx.value("verdict", std::string()) == "2tx_stbc_like" ? "STBC" : "SFBC";
        if (!label.empty() && label.back() == ')')  // inferred from the data statistics: "-like"
          label = label.substr(0, label.size() - 1) + ", 2-TX " + scheme + "-like)";
        else
          label += " (2-TX " + scheme + "-like)";
        mod += ", 2 transmit antennas (" + scheme + "-like)";
      }
      ev = "blind CP-OFDM: Tu " + py_str(get("tu_us")) + " us, CP " + py_str(get("cp_us")) + " us (" +
           py_str(get("cp_fraction")) + " Tu), " + py_str(get("bw_mhz")) + " MHz at " + py_str(get("carrier_mhz")) +
           " MHz, " + std::to_string(em.size()) + " burst(s), duty " + py_str(get("duty_cycle")) + ", frame period " +
           py_str(get("frame_period_ms")) + " ms, symbol-rate line z " + py_str(get("line_z")) + ", " +
           py_str(get("strong_windows")) + " strong 1 ms windows, SNR " + py_str(get("link_snr_db")) + " dB";
      if (!band_ctx.empty()) ev += "; band: " + band_ctx;
      if (tx.contains("verdict") && tx["verdict"].is_string() && !tx["verdict"].get<std::string>().empty()) {
        ev += "; transmit streams: " + tx["verdict"].get<std::string>();
        if (tx.contains("note") && tx["note"].is_string() && !tx["note"].get<std::string>().empty())
          ev += " (" + tx["note"].get<std::string>() + ")";
      }
      const Json cands = F.contains("system_candidates") && F["system_candidates"].is_array() ? F["system_candidates"] : Json::array();
      if (!cands.empty()) {
        auto join = [](const Json& a) {
          std::string r;
          for (const auto& x : a) r += (r.empty() ? "" : ", ") + x.get<std::string>();
          return r;
        };
        ev += "; most like " + cands[0]["system"].get<std::string>() + " (research catalog, matched " +
              join(cands[0]["matched"]) +
              (cands[0]["conflicts"].empty() ? std::string() : "; conflicts " + join(cands[0]["conflicts"])) + ")";
      }
    }
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    return result(status, method(), steps(), fields, reason, label, e.confidence, mod, ev,
                  std::make_pair(lo * 1e6, hi * 1e6));
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<DigitalVtxDecoder>()); }

}  // namespace si::rx::digital_vtx
