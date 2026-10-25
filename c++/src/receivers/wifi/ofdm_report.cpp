// Port of the report half of python receivers/wifi/plugin_ofdm.py (see ofdm_report.hpp).
#include "receivers/wifi/ofdm_report.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include "core/drone_models.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "receivers/wifi/field.hpp"
#include "receivers/wifi/links.hpp"
#include "receivers/wifi/payload.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::ofdm_report {

namespace {

using detail::jget;
using detail::json_key;
using detail::pstr;
using detail::truthy;

// collections.Counter -> dict (first-seen order)
void count(Json& d, const Json& key, double add = 1) {
  const std::string k = json_key(key);
  if (d.contains(k)) d[k] = d[k].get<int64_t>() + static_cast<int64_t>(add);
  else d[k] = static_cast<int64_t>(add);
}

std::string clip120(const std::string& s) { return detail::head_codepoints(s, 120); }

struct Payload {
  bool done = false;
  std::optional<payload::NonHt> nonht;
  std::optional<payload::Ht> ht;
  std::optional<std::string> error;
  bool decodable = true;
  std::string note;
};

std::string vht_note(const Json& bw, std::optional<int> nss, int ldpc, double fs) {
  const double b = bw.get<double>();
  if (b != 0 && b * 1e6 > fs) return pstr(bw) + " MHz payload is wider than the " + fmt("%g", fs / 1e6) + " MHz capture";
  if (nss && *nss > 1) return std::to_string(*nss) + " spatial streams (needs as many antennas)";
  if (ldpc) return "LDPC-coded VHT payload (only 20 MHz HT LDPC implemented)";
  return pstr(bw) + " MHz VHT payload decoding not implemented";
}

std::string he_note(const Json& bw, std::optional<int> nss, double fs) {
  const double b = bw.get<double>();
  if (b != 0 && b * 1e6 > fs) return pstr(bw) + " MHz payload is wider than the " + fmt("%g", fs / 1e6) + " MHz capture";
  if (nss && *nss > 1) return std::to_string(*nss) + " spatial streams (needs as many antennas)";
  return "HE (802.11ax) payload decoding not implemented";
}

// python _he_fields(mine, all_recs): airtime by direction, HE formats, BSS colours, colour -> BSSID (unique beacon)
void he_fields(const std::vector<Json*>& mine, const std::vector<Json>& all, Json& fields) {
  std::vector<const Json*> he;
  for (const Json* p : mine)
    if (jget(*p, "fmt") == "HE" && !jget(*p, "bss_color").is_null()) he.push_back(p);
  if (he.empty()) return;
  std::vector<std::pair<std::string, double>> air;  // Counter, first-seen order
  for (const Json* p : he) {
    const std::string d = (*p)["direction"].get<std::string>();
    const Json& du = jget(*p, "dur_us");
    const double v = truthy(du) ? du.get<double>() : 0.0;
    auto it = std::find_if(air.begin(), air.end(), [&](const auto& a) { return a.first == d; });
    if (it == air.end()) air.push_back({d, 0.0 + v});
    else it->second += v;
  }
  std::stable_sort(air.begin(), air.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
  Json a = Json::object();
  for (const auto& [k, v] : air) a[k] = py_round(v / 1e3, 3);
  fields["he_airtime_ms_by_direction"] = a;
  Json fm = Json::object();
  for (const Json* p : he) count(fm, (*p)["he_format"]);
  fields["he_formats"] = fm;
  std::map<int, int64_t> cols;
  for (const Json* p : he) ++cols[(*p)["bss_color"].get<int>()];
  Json c = Json::object();
  for (const auto& [k, v] : cols) c[std::to_string(k)] = v;
  fields["he_bss_colors"] = c;
  std::map<int, std::set<std::string>> bss;
  for (const auto& r : all) {
    if (!r.contains("macs")) continue;
    for (const auto& m : r["macs"]) {
      const Json& op = jget(m, "he_operation");
      if (!op.is_object()) continue;
      const Json& cc = jget(op, "bss_color");
      if (truthy(cc) && truthy(jget(m, "bssid"))) bss[cc.get<int>()].insert(m["bssid"].get<std::string>());
    }
  }
  Json hits = Json::object();
  for (const auto& [k, v] : cols) {
    auto it = bss.find(k);
    if (it != bss.end() && it->second.size() == 1) hits[std::to_string(k)] = *it->second.begin();
  }
  if (!hits.empty()) fields["he_bss_from_color"] = hits;
}

const std::vector<std::string>& steps_common() {
  static const std::vector<std::string> s = {
      "non-HT PPDUs: OFDM demod, deinterleave, Viterbi, descramble, CRC-32 FCS (receivers.wifi.payload)",
      "HT (802.11n) MCS 0-7, 20 MHz: channel per space-time stream from the HT-LTFs (smoothed when the HT-SIG "
      "allows), Alamouti combining for STBC, pilots p_3.. with the rotating HT pattern, per-symbol phase and "
      "timing-drift tracking; BCC: N_COL=13 deinterleave + Viterbi; LDPC: 802.11 codeword layout "
      "(shortening/puncturing/repetition), log-MAP LLRs, min-sum decoding, data-aided channel re-estimation "
      "from codewords that pass their parity checks; descramble; an A-MPDU is split on its delimiters and "
      "every MPDU checked on its own FCS",
      "MAC dissection (mac80211.py) incl. Remote ID (ASTM F3411 / DJI) when a drone broadcasts it",
      "transmitter attribution (links.py): MAC header, control responses (ACK/BA/CTS RA), sequential "
      "scrambler-seed counters read from the SERVICE field; received power as a weak hint"};
  return s;
}

constexpr int64_t kCtrlMaxBytes = 40;

// python _field_frames(recs): PPDU records -> field.py frames
std::vector<Json> field_frames(const std::vector<Json*>& recs) {
  std::vector<Json> out;
  for (const Json* p : recs) {
    Json ms = Json::array();
    if (truthy(jget(*p, "macs"))) ms = (*p)["macs"];
    else if (truthy(jget(*p, "mac"))) ms.push_back((*p)["mac"]);
    Json b = nullptr;
    if ((*p)["fmt"] == "NON_HT") b = (*p)["length"];
    else if ((*p)["fmt"] == "HT") b = jget(*p, "psdu_bytes");
    bool ctrl = false, video = false;
    for (const auto& m : ms) {
      ctrl = ctrl || jget(m, "type") == "ctrl";
      video = video || truthy(jget(jget(m, "ipv4"), "video"));
    }
    ctrl = ctrl || (ms.empty() && !b.is_null() && b.get<int64_t>() <= kCtrlMaxBytes);
    const Json& du = jget(*p, "dur_us");
    out.push_back({{"t0", (*p)["t0"]}, {"t1", (*p)["t1"]}, {"dur_us", truthy(du) ? du.get<double>() : 0.0},
                   {"bytes", b}, {"ta", jget(*p, "ta")}, {"macs", ms}, {"ctrl", ctrl}, {"video", video}});
  }
  return out;
}

// python _fields(mine, sifs_us, block_s) -> (fields, fcs records)
Json fields_of(const std::vector<Json*>& mine, double sifs_us, std::vector<Json*>& fcs, double block_s) {
  Json fmts = Json::object(), chans = Json::object(), bws = Json::object();
  for (const Json* p : mine) {
    count(fmts, (*p)["fmt"]);
    count(chans, (*p)["channel"]);
    count(bws, (*p)["bw_mhz"]);
  }
  fcs.clear();
  for (Json* p : mine)
    if (truthy(jget(*p, "fcs_ok"))) fcs.push_back(p);
  std::vector<Json> macs;
  for (const Json* p : fcs) {
    if (p->contains("macs")) {
      for (const auto& m : (*p)["macs"]) macs.push_back(m);
    } else if (p->contains("mac")) {
      macs.push_back((*p)["mac"]);
    }
  }
  const Json tx = links::attribute(mine, sifs_us);
  std::set<std::string> bss, ssid;
  for (const auto& m : macs) {
    if (truthy(jget(m, "bssid"))) bss.insert(m["bssid"].get<std::string>());
    if (truthy(jget(m, "ssid"))) ssid.insert(m["ssid"].get<std::string>());
  }
  Json f = Json::object();
  f["n_ppdus"] = mine.size();
  f["formats"] = fmts;
  f["channels"] = chans;
  f["bandwidths_mhz"] = bws;
  f["n_fcs_ok"] = fcs.size();
  f["n_mpdus"] = macs.size();
  f["bssids"] = Json(std::vector<std::string>(bss.begin(), bss.end()));
  f["ssids"] = Json(std::vector<std::string>(ssid.begin(), ssid.end()));
  f["transmitters"] = tx;
  int64_t n_att = 0;
  for (const Json* p : mine) n_att += truthy(jget(*p, "ta"));
  f["n_attributed"] = n_att;
  Json pp = Json::array();
  for (size_t i = 0; i < mine.size() && i < 200; ++i) {
    Json r = Json::object();
    for (const auto& [k, v] : mine[i]->items())
      if (k != "f0" && k != "f1" && k != "t0" && k != "t1") r[k] = v;
    pp.push_back(std::move(r));
  }
  f["ppdus"] = pp;
  int64_t n_ht = 0, n_ht_ok = 0;
  for (const Json* p : mine)
    if ((*p)["fmt"] == "HT") {
      ++n_ht;
      n_ht_ok += truthy(jget(*p, "fcs_ok"));
    }
  if (n_ht) {
    f["n_ht_ppdus"] = n_ht;
    f["n_ht_decoded"] = n_ht_ok;
  }
  {  // research, PLAN M2: manufacturer from the Remote ID serial (python: the last frame's remote_id)
    const Json* rid = nullptr;
    for (const auto& m : macs)
      if (truthy(jget(m, "remote_id"))) rid = &m["remote_id"];
    if (rid) f["model"] = drone_models::remote_id_model(*rid);
  }
  Json notes = Json::object();
  for (const Json* p : mine)
    if (truthy(jget(*p, "payload_note"))) count(notes, (*p)["payload_note"]);
  if (!notes.empty()) f["payload_not_decodable"] = notes;
  if (block_s > 0) {  // research, PLAN W1
    const Json fs = field_summary(field_frames(mine), block_s);
    for (auto it = fs.begin(); it != fs.end(); ++it) f[it.key()] = it.value();
  }
  return f;
}

std::pair<double, double> span_of(const std::vector<Json*>& recs) {
  double lo = 0, hi = 0;
  for (size_t i = 0; i < recs.size(); ++i) {
    const double a = (*recs[i])["f0"].get<double>(), b = (*recs[i])["f1"].get<double>();
    if (i == 0 || a < lo) lo = a;
    if (i == 0 || b > hi) hi = b;
  }
  return {lo, hi};
}

}  // namespace

std::vector<Json> records(const std::vector<ofdm::Packet>& pkts, double fs, double clock_hz, bool scaled, int nthreads) {
  const double scale = clock_hz / 20e6;
  // payload decoding, every PPDU independent
  std::vector<Payload> pay(pkts.size());
  parallel_for(
      pkts.size(),
      [&](size_t i) {
        const ofdm::Packet& p = pkts[i];
        Payload& out = pay[i];
        if (!p.rx) return;
        try {
          if (p.fmt == "NON_HT") {
            out.nonht = payload::decode_nonht_payload(*p.rx, p.lsig_rate_mbps, p.lsig_length);
          } else if (p.fmt == "HT") {
            const double bw = clock_hz != 20e6 ? 20.0 : p.bandwidth_mhz;
            auto dec = payload::ht_decodable(p.mcs, p.ldpc, p.stbc, p.n_ess, bw);
            if (!dec.first) {
              out.decodable = false;
              out.note = dec.second;
              return;
            }
            out.ht = payload::decode_ht_payload(*p.rx, *p.mcs, p.ht_length, p.short_gi, p.aggregation, p.n_ess,
                                                p.chan_offset_hz, p.chan_usable_half_hz, p.stbc, p.ldpc, p.smoothing);
          }
        } catch (const std::exception& e) {
          out.error = clip120(e.what());
        }
      },
      nthreads);
  std::vector<Json> out;
  for (size_t i = 0; i < pkts.size(); ++i) {
    const ofdm::Packet& p = pkts[i];
    Json rec = Json::object();
    rec["t0"] = static_cast<double>(p.start_sample) / fs;
    rec["t1"] = static_cast<double>(p.end_sample) / fs;
    rec["f0"] = p.f_lo_mhz * 1e6;
    rec["f1"] = p.f_hi_mhz * 1e6;
    rec["channel"] = p.channel;
    rec["center_mhz"] = p.center_mhz;
    rec["fmt"] = p.fmt;
    const Json bwj = scaled ? Json(p.bandwidth_mhz) : Json(static_cast<int64_t>(p.bandwidth_mhz));
    rec["bw_mhz"] = bwj;
    rec["rate_mbps"] = py_round(p.lsig_rate_mbps * scale, 3);
    rec["length"] = p.lsig_length;
    rec["mcs"] = p.mcs ? Json(*p.mcs) : Json(nullptr);
    rec["sig_crc_ok"] = p.sig_crc_ok ? Json(*p.sig_crc_ok) : Json(nullptr);
    rec["snr_db"] = p.snr_db;
    rec["cfo_khz"] = p.cfo_khz;
    rec["start_us"] = p.start_us;
    rec["dur_us"] = p.duration_us;
    if (p.rx) {
      const cdvec& r = *p.rx->r;
      const int64_t a = std::min<int64_t>(p.rx->start, static_cast<int64_t>(r.size()));
      const int64_t b = std::min<int64_t>(p.rx->start + 320, static_cast<int64_t>(r.size()));
      if (b > a) {
        dvec e(b - a);
        for (int64_t k = a; k < b; ++k) e[k - a] = np::cabs2(r[k]);
        rec["power_db"] = py_round(10 * std::log10(np::mean(e) + 1e-30), 1);
      }
    }
    const Payload& pl = pay[i];
    if (p.fmt == "NON_HT" && p.rx) {
      if (pl.error) {
        rec["payload_error"] = *pl.error;
      } else {
        const auto& d = *pl.nonht;
        rec["fcs_ok"] = d.fcs_ok;
        int64_t seed = 0;
        for (size_t k = 0; k < d.seed.size(); ++k) seed += static_cast<int64_t>(d.seed[k]) << k;
        rec["scrambler_seed"] = seed;
        rec["seed_ok"] = d.fcs_ok || d.service_ok;
        if (d.fcs_ok) {
          Json mac = Json::object();
          for (const auto& [k, v] : d.info.items())
            if (k.empty() || k[0] != '_') mac[k] = v;
          rec["mac"] = mac;
          rec["macs"] = Json::array({mac});
        }
      }
    } else if (p.fmt == "HT" && p.rx) {
      rec["psdu_bytes"] = p.ht_length;
      rec["ht_rate_mbps"] = payload::ht_mcs_known(p.mcs)
                                ? Json(py_round(payload::ht_rate_mbps(*p.mcs, p.short_gi) * scale, 3))
                                : Json(nullptr);
      rec["short_gi"] = p.short_gi != 0;
      rec["aggregated"] = p.aggregation != 0;
      rec["coding"] = p.ldpc ? "LDPC" : "BCC";
      rec["stbc"] = p.stbc;
      if (!pl.decodable) {
        rec["payload_note"] = pl.note;
      } else if (pl.error) {
        rec["payload_error"] = *pl.error;
      } else {
        const auto& r = *pl.ht;
        rec["fcs_ok"] = r.n_fcs_ok > 0;
        rec["n_subframes"] = r.subframes.size();
        rec["n_fcs_ok"] = r.n_fcs_ok;
        rec["scrambler_seed"] = r.scrambler_seed;
        rec["seed_ok"] = r.service_ok || r.n_fcs_ok > 0;
        if (r.ldpc) rec["ldpc_codewords_ok"] = std::to_string(r.n_codewords_ok) + "/" + std::to_string(r.n_codewords);
        Json macs = Json::array();
        for (const auto& s : r.subframes)
          if (s.mac && truthy(*s.mac)) macs.push_back(*s.mac);
        if (!macs.empty()) {
          rec["mac"] = macs[0];
          rec["macs"] = macs;
        }
      }
    } else if (p.fmt == "VHT" && p.sig_crc_ok && *p.sig_crc_ok) {
      rec["nss"] = p.nsts_su;
      rec["ldpc"] = p.ldpc != 0;
      if (p.group_id) {
        const int g = *p.group_id;
        rec["direction"] = g == 0 ? "uplink (station to AP)" : g == 63 ? "downlink (AP to station)" : "multi-user MIMO";
      } else {
        rec["direction"] = nullptr;
      }
      if (p.partial_aid) rec["partial_aid"] = *p.partial_aid;
      rec["payload_note"] = vht_note(bwj, p.nsts_su, p.ldpc, fs);
    } else if (p.fmt == "HE" && p.sig_crc_ok && *p.sig_crc_ok) {
      rec["he_format"] = p.he_format;
      rec["bss_color"] = p.bss_color;
      rec["txop_us"] = p.txop_us ? Json(*p.txop_us) : Json(nullptr);
      rec["direction"] = p.he_uplink ? "uplink (station to AP)" : "downlink (AP to station)";
      if (p.he_nsts) {
        rec["nss"] = *p.he_nsts;
        rec["ldpc"] = p.ldpc != 0;
        rec["dcm"] = p.dcm != 0;
      }
      rec["payload_note"] = he_note(bwj, p.he_nsts, fs);
    }
    out.push_back(std::move(rec));
  }
  return out;
}

void attribute_vht(std::vector<Json>& recs) {
  std::set<std::string> bss;
  for (const auto& r : recs) {
    if (!r.contains("macs")) continue;
    for (const auto& m : r["macs"]) {
      if (!truthy(jget(m, "bssid"))) continue;
      const std::string b = m["bssid"].get<std::string>();
      if (!(std::stoi(b.substr(0, 2), nullptr, 16) & 1)) bss.insert(b);
    }
  }
  std::map<int, std::vector<std::string>> paid;
  for (const auto& b : bss) paid[payload::vht_uplink_paid(b)].push_back(b);
  for (auto& r : recs) {
    const Json& d = jget(r, "direction");
    if (!(d.is_string() && d.get<std::string>().rfind("uplink", 0) == 0)) continue;
    if (!truthy(jget(r, "partial_aid"))) continue;
    auto it = paid.find(r["partial_aid"].get<int>());
    if (it != paid.end() && it->second.size() == 1) r["bssid_from_partial_aid"] = it->second[0];
  }
}

Summary summarise_ofdm(const std::vector<Json*>& mine, size_t n_all, size_t n_rej, const std::vector<Json>* all,
                       double block_s) {
  Summary s;
  s.steps = {"capture-wide 802.11 legacy-preamble scan (receivers.wifi.ofdm_rx.WifiDetector.run): every 20 MHz "
             "channel resampled to 20 MSps, L-STF autocorrelation plateau, L-LTF channel estimate, "
             "L-SIG Viterbi + parity + bit-agreement check, HT-SIG/VHT-SIG-A CRC-8, HE-SIG-A CRC-4 (BSS colour, UL/DL, MCS, "
             "streams, TXOP)",
             "off the standard grid / outside the Wi-Fi bands: 20 MHz slices every 5 MHz with the decoded PPDUs "
             "blanked; >= 3 adjacent slices with a 3.2 us CP correlation -> blind 1 MHz centres around them "
             "(ofdm_rx.offgrid_plan)",
             std::to_string(n_all) + " PPDUs in capture, " + std::to_string(n_rej) + " preambles rejected; " +
                 std::to_string(mine.size()) + " attributed to this signal"};
  for (const auto& t : steps_common()) s.steps.push_back(t);
  if (mine.empty()) {
    s.status = "not_possible";
    s.method = "receivers.wifi ofdm_rx L-STF/L-LTF/L-SIG receiver";
    s.reason = "no 802.11 OFDM preamble found in this signal's time-frequency span";
    return s;
  }
  std::vector<Json*> fcs;
  Json fields = fields_of(mine, 16.0, fcs, block_s);
  const Json fmts = fields["formats"];
  std::vector<const Json*> vht;
  for (const Json* p : mine)
    if (truthy(jget(*p, "direction")) && jget(*p, "fmt") == "VHT") vht.push_back(p);
  if (!vht.empty()) {
    std::vector<std::pair<std::string, double>> air;  // Counter, first-seen order
    for (const Json* p : vht) {
      const std::string d = (*p)["direction"].get<std::string>();
      const Json& du = jget(*p, "dur_us");
      const double v = truthy(du) ? du.get<double>() : 0.0;
      auto it = std::find_if(air.begin(), air.end(), [&](const auto& a) { return a.first == d; });
      if (it == air.end()) air.push_back({d, 0.0 + v});
      else it->second += v;
    }
    std::stable_sort(air.begin(), air.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    Json a = Json::object();
    for (const auto& [k, v] : air) a[k] = py_round(v / 1e3, 3);
    fields["vht_airtime_ms_by_direction"] = a;
    Json st = Json::object();
    for (const Json* p : vht) count(st, jget(*p, "nss"));
    fields["vht_streams"] = st;
    Json hits = Json::object();
    for (const Json* p : vht)
      if (truthy(jget(*p, "bssid_from_partial_aid"))) count(hits, (*p)["bssid_from_partial_aid"]);
    if (!hits.empty()) fields["vht_uplink_bss_from_partial_aid"] = hits;
  }
  if (all && !all->empty()) {
    he_fields(mine, *all, fields);
  } else {
    std::vector<Json> own;
    for (const Json* p : mine) own.push_back(*p);
    he_fields(mine, own, fields);
  }
  // Counter(channel).most_common(1)
  std::string ch;
  int64_t best = -1;
  for (const auto& [k, v] : fields["channels"].items())
    if (v.get<int64_t>() > best) {
      best = v.get<int64_t>();
      ch = k;
    }
  std::vector<std::string> fl;
  for (const auto& [k, v] : fmts.items()) fl.push_back(k);
  std::sort(fl.begin(), fl.end());
  std::string join;
  for (size_t i = 0; i < fl.size(); ++i) join += (i ? "/" : "") + fl[i];
  bool any_crc = false;
  for (const Json* p : mine) any_crc |= truthy(jget(*p, "sig_crc_ok"));
  s.status = fcs.empty() ? "partial" : "success";
  s.method = "receivers.wifi ofdm_rx + payload (research receiver)";
  s.reason = fcs.empty() ? "PHY headers (L-SIG/HT-SIG) decoded; payloads need multiple antennas "
                           "(2+ spatial streams) or a wider channel, or they failed FCS"
                         : "";
  const auto [suf, dev] = drone_label(fields);
  s.label = "Wi-Fi 802.11 OFDM (" + join + ") " + ch + suf;
  s.confidence = fcs.empty() ? 0.9 : 0.97;
  s.modulation = "OFDM 64-pt, 312.5 kHz SCS (BPSK..1024-QAM)";
  s.evidence = std::to_string(mine.size()) + " 802.11 PPDUs decoded (L-SIG parity/agreement" +
               (any_crc ? ", HT/VHT/HE-SIG CRC" : "") + ")" +
               (fcs.empty() ? "" : "; " + std::to_string(fcs.size()) + " MPDUs pass CRC-32 FCS") + dev;
  s.span = span_of(mine);
  s.fields = std::move(fields);
  return s;
}

Summary summarise_scaled(const scaled::Link& ln, const std::vector<Json*>& recs, double block_s) {
  const scaled::Numerology& nm = ln.numerology;
  const double clk = ln.clock_hz / 1e6;
  const std::string g = fmt("%g", clk);
  Summary s;
  s.steps = {
      "bursty wideband bands with repeatable edges from a ~60 kHz spectrogram (ofdm_scaled.find_bands)",
      fmt("numerology of the band: CP autocorrelation peak at T_u = %.3f us (%.2f, %sx median), symbol T_s = %.3f us, ",
          nm.tu_s * 1e6, nm.cp_corr, py_float_repr(nm.cp_corr_over_median).c_str(), nm.ts_s * 1e6) +
          fmt("guard %s of T_u (802.11: 0.25) -> 64-subcarrier clock %.4f MHz, taken as %s MHz; occupied %.2f MHz = "
              "%s x clock",
              py_float_repr(nm.guard).c_str(), *nm.clock_measured_hz / 1e6, g.c_str(), nm.occupied_hz / 1e6,
              py_float_repr(*nm.occupied_over_clock).c_str()),
      "capture time-scaled by 20/" + g + " into a standard 20 MHz channel" +
          (ln.inverted ? " after mirroring the spectrum (complex conjugate): only the mirrored orientation gives "
                         "parity/CRC-valid SIG fields"
                       : "") +
          "; ofdm_rx: L-STF plateau, L-LTF, L-SIG parity + bit agreement, HT-SIG CRC",
      std::to_string(recs.size()) + " PPDUs, " + std::to_string(ln.n_rejects) + " preambles rejected"};
  for (const auto& t : steps_common()) s.steps.push_back(t);
  std::vector<Json*> fcs;
  Json fields = fields_of(recs, 16.0 * 20 / clk, fcs, block_s);
  const Json fmts = fields["formats"];
  fields["center_mhz"] = py_round(ln.center_hz / 1e6, 4);
  fields["clock_mhz"] = clk;
  fields["channel_width_mhz"] = clk;
  fields["subcarrier_spacing_khz"] = py_round(312.5 * clk / 20, 3);
  fields["spectrum_inverted"] = ln.inverted;
  fields["numerology"] = nm.json();
  fields["standard_equivalent"] = clk == 5    ? "802.11 quarter-clocked (5 MHz) channel"
                                  : clk == 10 ? "802.11 half-clocked (10 MHz) channel"
                                  : clk == 20 ? "802.11 at 20 MHz, off the channel grid"
                                              : "802.11 run at a " + g + " MHz clock (non-standard)";
  std::set<std::string> vendors;
  for (const auto& t : fields["transmitters"]) {
    const Json& v = jget(t, "vendor");
    if (!truthy(v)) continue;
    const std::string vs = v.get<std::string>();
    if (detail::ascii_lower(vs).find("unknown") != std::string::npos) continue;
    if (!truthy(jget(jget(t, "evidence"), "mac header"))) continue;
    vendors.insert(vs);
  }
  const bool ht = fmts.contains("HT");
  std::string lab = "Wi-Fi " + std::string(ht ? "802.11n" : "802.11a/g") + "-type link, " + g + " MHz channel at " +
                    fmt("%.2f", ln.center_hz / 1e6) + " MHz" + (ln.inverted ? ", spectrum-inverted" : "");
  if (!vendors.empty()) {
    std::string j;
    for (const auto& v : vendors) j += (j.empty() ? "" : ", ") + v;
    lab += " (" + j + ")";
  }
  std::string ev = std::to_string(recs.size()) + " PPDUs with parity/CRC-checked SIG fields at a " + g + " MHz OFDM clock";
  if (!fcs.empty()) ev += "; " + std::to_string(fcs.size()) + " PPDUs carry CRC-32-valid MPDUs";
  if (!fields["transmitters"].empty()) {
    std::string j;
    for (const auto& t : fields["transmitters"]) j += (j.empty() ? "" : ", ") + pstr(t["ta"]);
    ev += "; transmitters " + j;
  }
  {
    const auto [suf, dev] = drone_label(fields);
    lab += suf;
    ev += dev;
  }
  s.status = fcs.empty() ? "partial" : "success";
  s.method = "receivers.wifi ofdm_scaled + ofdm_rx + payload (time-scaled 802.11)";
  s.reason = fcs.empty() ? "PHY headers decoded; no MPDU passed its FCS" : "";
  s.label = lab;
  s.confidence = fcs.empty() ? 0.9 : 0.97;
  s.modulation = fmt("OFDM 64-pt, %.2f kHz SCS (802.11 at a %s MHz clock)", 312.5 * clk / 20, g.c_str());
  s.evidence = ev;
  s.span = span_of(recs);
  s.fields = std::move(fields);
  return s;
}

Json field_summary(const std::vector<Json>& frames, double block_s) {
  Json out = Json::object();
  const Json lk = field::links(frames, block_s);
  if (lk.empty()) return out;
  Json top = Json::array(), dl = Json::array();
  for (size_t i = 0; i < lk.size() && i < 10; ++i) top.push_back(lk[i]);
  out["wifi_links"] = top;
  for (const auto& r : lk) {
    if (!r.contains("drone_link")) continue;
    Json d = {{"bssid", r["bssid"]}, {"ssid", jget(r, "ssid")}, {"clients", r["clients"]}};
    for (auto it = r["drone_link"].begin(); it != r["drone_link"].end(); ++it) d[it.key()] = it.value();
    dl.push_back(d);
  }
  if (!dl.empty()) out["drone_links"] = dl;
  return out;
}

std::pair<std::string, std::string> drone_label(const Json& fields) {
  const Json& dl = jget(fields, "drone_links");
  if (!truthy(dl)) return {"", ""};
  auto conf = [](const Json& x) { return x["confidence"].is_number() ? x["confidence"].get<double>() : 0.0; };
  const Json* d = &dl[0];
  for (const auto& x : dl)
    if (conf(x) > conf(*d)) d = &x;
  const std::string what = (*d)["platform"].is_null() ? "by behaviour" : pstr((*d)["platform"]);
  return {", drone link-like (" + what + ")", "; drone link-like on BSS " + pstr((*d)["bssid"]) + " (" + what + ", " +
                                                  pstr((*d)["basis"]) + ", confidence " + pstr((*d)["confidence"]) + ")"};
}

}  // namespace si::rx::wifi::ofdm_report
