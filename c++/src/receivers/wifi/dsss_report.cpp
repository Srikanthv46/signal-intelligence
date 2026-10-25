// Port of the summaries in python receivers/wifi/plugin_dsss.py (see dsss_report.hpp).
#include "receivers/wifi/dsss_report.hpp"
#include "receivers/wifi/ofdm_report.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include "core/bands_india.hpp"
#include "core/drone_models.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::dsss {

using detail::jget;
using detail::pstr;
using detail::truthy;

namespace {

const char* kMethod = "802.11b DSSS/CCK receiver (this package)";
const char* kMethodBarker = "802.11b DSSS receiver (this package)";

double wifi24(int ch) {
  for (const auto& [c, f] : bands::wifi_24())
    if (c == ch) return f;
  return 0.0;
}

void count(Json& d, const Json& key) {
  const std::string k = detail::json_key(key);
  if (d.contains(k)) d[k] = d[k].get<int64_t>() + 1;
  else d[k] = 1;
}

bool group_addr(const std::string& mac) { return std::stoi(mac.substr(0, 2), nullptr, 16) & 1; }

Json sorted_strings(const std::set<std::string>& s) {
  Json out = Json::array();
  for (const auto& x : s) out.push_back(x);
  return out;
}

const char* kBssKeys[] = {"ssid", "ssid_hidden", "security", "phy", "ap_uptime", "beacon_interval_ms",
                          "capabilities", "rates_mbps", "basic_rates_mbps", "country", "bss_load", "ht",
                          "ht_operation", "vht_operation", "wps", "p2p", "mobility_domain",
                          "radio_measurement_802_11k", "interworking_802_11u", "vendor_elements", "erp", "tim"};

}  // namespace

const std::vector<std::string>& steps() {
  static const std::vector<std::string> s = {
      "channelise each 2.4 GHz channel: 22 MSps (detection) and 44 MSps (4 samples/chip, decoding)",
      "detect: |Barker-11 matched filter|^2 summed over 16 symbols at each of the 22 sample phases; DSSS where "
      "the best phase stands out (peak/mean > 2.8, noise < 2.6); runs split where the best phase moves",
      "acquire SFD-first: Barker symbols on a 1 us grid at every sample phase, DBPSK, descramble x^7+x^4+1, "
      "SFD (long 0xF3A0 / short 0x05CF, <= 1 bit error) after a SYNC run; sidelobe hits merged",
      "track each PPDU from its own SYNC (early/late loop gated on its amplitude); rebuild the scrambled "
      "SYNC+SFD exactly (best of 128 scrambler states) -> data-aided CFO, phase, SNR",
      "PLCP header coherent (phase-tracked), differential fallback; CRC-16",
      "PSDU: 1/2 Mb/s coherent Barker symbols; 5.5/11 Mb/s CCK on a CFO-locked chip grid through a 9-tap "
      "LS equaliser trained on preamble + header chips, ML codeword + phi1 decisions",
      "descramble, CRC-32 FCS; on failure FCS-guided correction of <= 4 of the 20 least reliable symbol "
      "decisions, then STFT interference excision and a re-decode",
      "MAC dissection (mac80211.py): addresses and roles, IEs, security, capabilities, WPS / P2P "
      "identity, EAPOL / ARP / IP headers when unencrypted; vendor from the IEEE registry",
      "dedup: a PPDU seen on neighbouring channels is kept where its Barker amplitude is highest; "
      "a beacon's DS parameter names the channel"};
  return s;
}

Json frame_line(const Json& p) {
  const Json& m = jget(p, "mac");
  Json o = Json::object();
  o["t_ms"] = py_round(p["start_s"].get<double>() * 1e3, 3);
  o["rate_mbps"] = jget(p, "rate_mbps");
  o["preamble"] = jget(p, "preamble");
  o["length_us"] = jget(p, "length_us");
  o["snr_db"] = jget(p, "snr_db");
  o["cfo_khz"] = py_round(p["cfo_hz"].get<double>() / 1e3, 1);
  o["fcs_ok"] = truthy(jget(p, "fcs_ok"));
  o["symbols_corrected"] = p.contains("symbols_corrected") ? p["symbols_corrected"] : Json(0);
  const Json& s = jget(m, "summary");
  o["summary"] = truthy(s) ? s : Json(truthy(jget(p, "psdu_bytes")) ? "FCS failed" : "header only");
  return o;
}

Json bss_table(const std::vector<const Json*>& frames) {
  std::vector<Json> rows;
  std::vector<std::set<std::string>> stations;
  std::vector<dvec> snrs;
  for (const Json* pp : frames) {
    const Json& p = *pp;
    const Json& m = p["mac"];
    const Json& bj = jget(m, "bssid");
    if (!truthy(bj) || bj == "ff:ff:ff:ff:ff:ff") continue;
    const std::string b = bj.get<std::string>();
    size_t k = 0;
    while (k < rows.size() && rows[k]["bssid"] != b) ++k;
    if (k == rows.size()) {
      Json r = Json::object();
      r["bssid"] = b;
      r["vendor"] = jget(m, "bssid_vendor");
      r["n_frames"] = 0;
      r["frame_types"] = Json::object();
      r["stations"] = Json::array();
      r["snr_db"] = Json::array();
      rows.push_back(r);
      stations.emplace_back();
      snrs.emplace_back();
    }
    Json& r = rows[k];
    r["n_frames"] = r["n_frames"].get<int64_t>() + 1;
    count(r["frame_types"], jget(m, "subtype"));
    if (!jget(p, "snr_db").is_null() && jget(m, "ta") == b) snrs[k].push_back(p["snr_db"].get<double>());
    const Json& st = jget(m, "subtype");
    if (st == "beacon" || st == "probe-resp") {
      r["channel"] = truthy(jget(m, "ds_channel")) ? m["ds_channel"] : p["channel"];
      r["cfo_hz"] = p["cfo_hz"];
      for (const char* key : kBssKeys)
        if (m.contains(key)) r[key] = m[key];
    }
    for (const char* key : {"sa", "da", "ta", "ra"}) {
      const Json& v = jget(m, key);
      if (truthy(v) && v != b && !group_addr(v.get<std::string>())) stations[k].insert(v.get<std::string>());
    }
  }
  for (size_t k = 0; k < rows.size(); ++k) {
    rows[k]["stations"] = sorted_strings(stations[k]);
    rows[k]["snr_db"] = snrs[k].empty() ? Json(nullptr) : Json(py_round(np::median(snrs[k]), 1));
  }
  std::stable_sort(rows.begin(), rows.end(), [](const Json& a, const Json& b) {
    const bool na = jget(a, "ssid").is_null(), nb = jget(b, "ssid").is_null();
    if (na != nb) return !na;
    return a["n_frames"].get<int64_t>() > b["n_frames"].get<int64_t>();
  });
  Json out = Json::array();
  for (auto& r : rows) out.push_back(std::move(r));
  return out;
}

Json station_table(const std::vector<const Json*>& frames) {
  std::vector<Json> aps;
  for (const Json* p : frames) {
    const Json& m = (*p)["mac"];
    const Json& st = jget(m, "subtype");
    if (st == "beacon" || st == "probe-resp") aps.push_back(jget(m, "bssid"));
  }
  std::vector<Json> rows;
  std::vector<std::set<std::string>> probed;
  for (const Json* pp : frames) {
    const Json& m = (*pp)["mac"];
    const Json& txj = truthy(jget(m, "sa")) ? m["sa"] : jget(m, "ta");
    if (!truthy(txj) || std::find(aps.begin(), aps.end(), txj) != aps.end()) continue;
    const std::string tx = txj.get<std::string>();
    if (group_addr(tx)) continue;
    size_t k = 0;
    while (k < rows.size() && rows[k]["mac"] != tx) ++k;
    if (k == rows.size()) {
      Json r = Json::object();
      r["mac"] = tx;
      r["vendor"] = truthy(jget(m, "sa_vendor")) ? m["sa_vendor"] : jget(m, "ta_vendor");
      r["randomized"] = static_cast<bool>(std::stoi(tx.substr(0, 2), nullptr, 16) & 2);
      r["n_frames"] = 0;
      r["frame_types"] = Json::object();
      r["probed_ssids"] = Json::array();
      r["bssid"] = nullptr;
      rows.push_back(r);
      probed.emplace_back();
    }
    Json& r = rows[k];
    r["n_frames"] = r["n_frames"].get<int64_t>() + 1;
    count(r["frame_types"], jget(m, "subtype"));
    const Json& ssid = jget(m, "ssid");
    if (jget(m, "subtype") == "probe-req" && !ssid.is_null() && ssid != "<broadcast>")
      probed[k].insert(ssid.get<std::string>());
    const Json& bs = jget(m, "bssid");
    if (truthy(bs) && bs != "ff:ff:ff:ff:ff:ff" && jget(m, "type") == "data") r["bssid"] = bs;
    const Json& n1 = jget(jget(m, "p2p"), "device_name");
    const Json& name = truthy(n1) ? n1 : jget(jget(m, "wps"), "device_name");
    if (truthy(name)) r["device_name"] = name;
    if (truthy(jget(m, "phy"))) r["phy"] = m["phy"];
  }
  for (size_t k = 0; k < rows.size(); ++k) rows[k]["probed_ssids"] = sorted_strings(probed[k]);
  std::stable_sort(rows.begin(), rows.end(), [](const Json& a, const Json& b) {
    return a["n_frames"].get<int64_t>() > b["n_frames"].get<int64_t>();
  });
  Json out = Json::array();
  for (auto& r : rows) out.push_back(std::move(r));
  return out;
}

Summary summarise(const std::vector<const Json*>& mine, size_t n_all, double block_s) {
  Summary res;
  res.steps = steps();
  res.steps.push_back(fmt("%zu PLCP headers with valid CRC-16 in capture; %zu attributed to this signal", n_all,
                          mine.size()));
  res.method = kMethod;
  if (mine.empty()) {
    res.method = kMethodBarker;
    res.status = "not_possible";
    res.reason = "no Barker-spread PLCP header with valid CRC-16 in this signal's span";
    return res;
  }
  std::vector<const Json*> fcs;
  for (const Json* p : mine)
    if (truthy(jget(*p, "fcs_ok"))) fcs.push_back(p);
  // Counter(channel).most_common(1): the first channel with the highest count
  std::vector<std::pair<int, int>> cc;
  for (const Json* p : mine) {
    const int c = (*p)["channel"].get<int>();
    auto it = std::find_if(cc.begin(), cc.end(), [&](const auto& e) { return e.first == c; });
    if (it == cc.end()) cc.push_back({c, 1});
    else ++it->second;
  }
  int ch = cc[0].first, best = cc[0].second;
  for (const auto& [c, n] : cc)
    if (n > best) ch = c, best = n;
  const Json bss = bss_table(fcs), stations = station_table(fcs);
  dvec snr, cfo;
  int64_t ncorr = 0;
  Json rates = Json::object(), pres = Json::object(), chans = Json::object(), ftypes = Json::object();
  std::set<std::string> bssids, ssids;
  for (const Json* p : mine) {
    if (!jget(*p, "snr_db").is_null()) snr.push_back((*p)["snr_db"].get<double>());
    cfo.push_back((*p)["cfo_hz"].get<double>());
    count(rates, jget(*p, "rate_mbps"));
    count(pres, jget(*p, "preamble"));
    count(chans, jget(*p, "channel"));
  }
  for (const Json* p : fcs) {
    if (truthy(jget(*p, "symbols_corrected"))) ++ncorr;
    const Json& m = (*p)["mac"];
    if (truthy(jget(m, "bssid"))) bssids.insert(m["bssid"].get<std::string>());
    if (truthy(jget(m, "ssid"))) ssids.insert(m["ssid"].get<std::string>());
    count(ftypes, jget(m, "subtype"));
  }
  const double fch = wifi24(ch);
  Json f = Json::object();
  f["channel"] = ch;
  f["center_mhz"] = fch / 1e6;
  f["n_ppdus"] = mine.size();
  f["n_fcs_ok"] = fcs.size();
  f["n_fcs_corrected"] = ncorr;
  f["rates_mbps"] = rates;
  f["preambles"] = pres;
  f["channels"] = chans;
  f["snr_db_inband_median"] = snr.empty() ? Json(nullptr) : Json(py_round(np::median(snr), 1));
  f["cfo_hz_median"] = static_cast<int64_t>(np::median(cfo));
  f["bssids"] = sorted_strings(bssids);
  f["ssids"] = sorted_strings(ssids);
  f["frame_types"] = ftypes;
  f["n_frames"] = fcs.size();
  f["bss"] = bss;
  f["stations"] = stations;
  Json frames = Json::array(), ppdus = Json::array();
  for (size_t i = 0; i < mine.size() && i < 300; ++i) frames.push_back(frame_line(*mine[i]));
  for (size_t i = 0; i < mine.size() && i < 200; ++i) ppdus.push_back(*mine[i]);
  f["frames"] = frames;
  f["ppdus"] = ppdus;
  {  // research, PLAN M2: manufacturer from the Remote ID serial (python: the last frame's remote_id)
    const Json* rid = nullptr;
    for (const Json* p : fcs)
      if (truthy(jget((*p)["mac"], "remote_id"))) rid = &(*p)["mac"]["remote_id"];
    if (rid) f["model"] = drone_models::remote_id_model(*rid);
  }
  if (block_s > 0) {  // research, PLAN W1 (python plugin_dsss: field.py frames)
    std::vector<Json> frames;
    for (const Json* p : mine) {
      const Json& m = jget(*p, "mac");
      const bool ok = truthy(jget(*p, "fcs_ok"));
      const Json& pb = jget(*p, "psdu_bytes");
      const int64_t pbn = truthy(pb) ? pb.get<int64_t>() : 99;
      Json ta = jget(m, "ta");
      if (!truthy(ta)) ta = jget(m, "sa");
      const double t0 = (*p)["start_s"].get<double>(), t1 = (*p)["end_s"].get<double>();
      frames.push_back({{"t0", t0}, {"t1", t1}, {"dur_us", (t1 - t0) * 1e6}, {"bytes", pb}, {"ta", ta},
                        {"macs", ok && truthy(m) ? Json::array({m}) : Json::array()},
                        {"ctrl", jget(m, "type") == "ctrl" || (!ok && pbn <= 40)},
                        {"video", truthy(jget(jget(m, "ipv4"), "video"))}});
    }
    const Json fs = ofdm_report::field_summary(frames, block_s);
    for (auto it = fs.begin(); it != fs.end(); ++it) f[it.key()] = it.value();
  }
  std::string names;
  for (const auto& b : bss)
    if (truthy(jget(b, "ssid"))) names += (names.empty() ? "" : ", ") + ("'" + pstr(b["ssid"]) + "'");
  names = detail::head_codepoints(names, 200);
  const bool ok = !fcs.empty();
  res.status = ok ? "success" : "partial";
  res.fields = f;
  res.reason = ok ? "" : "PLCP headers decoded (CRC-16 ok) but no PSDU passed FCS";
  const auto [suf, dev] = ofdm_report::drone_label(f);
  res.label = "Wi-Fi 802.11b DSSS ch" + std::to_string(ch) + suf;
  res.confidence = ok ? 0.97 : 0.9;
  res.modulation = "DSSS Barker-11 DBPSK/DQPSK, CCK";
  res.span = std::make_pair(fch - 11e6, fch + 11e6);
  std::string ev = std::to_string(mine.size()) + " 802.11b PLCP headers with valid CRC-16";
  if (ok) ev += "; " + std::to_string(fcs.size()) + " MPDUs pass FCS";
  if (!bss.empty()) ev += "; " + std::to_string(bss.size()) + " BSS: " + names;
  res.evidence = ev + dev;
  return res;
}

Summary barker_only(const std::vector<const Json*>& regs) {
  const Json* r = regs[0];
  double t0min = (*regs[0])["t0"].get<double>(), dur = 0;
  for (const Json* x : regs) {
    if ((*x)["coherence"].get<double>() > (*r)["coherence"].get<double>()) r = x;
    t0min = std::min(t0min, (*x)["t0"].get<double>());
    dur += (*x)["t1"].get<double>() - (*x)["t0"].get<double>();
  }
  const bool early = t0min < 60e-6;
  Json f = Json::object();
  f["channel"] = (*r)["channel"];
  f["center_mhz"] = (*r)["ch_freq_mhz"];
  f["n_regions"] = regs.size();
  f["barker_coherence"] = (*r)["coherence"];
  f["duration_us"] = py_round_int(dur * 1e6);
  Json rg = Json::array();
  for (size_t i = 0; i < regs.size() && i < 50; ++i) {
    const Json& x = *regs[i];
    Json o = Json::object();
    o["t_start_s"] = py_round(x["t0"].get<double>(), 6);
    o["t_end_s"] = py_round(x["t1"].get<double>(), 6);
    o["channel"] = x["channel"];
    o["coherence"] = x["coherence"];
    rg.push_back(o);
  }
  f["regions"] = rg;
  Summary res;
  std::string why;
  if (early) {
    res.label = "Wi-Fi 802.11b DSSS-like (preamble before capture start)";
    res.confidence = 0.7;
    why = "Barker-11 spreading at 1 Msym/s from the start of the capture: the PPDU began before it";
  } else {
    res.label = "Barker-11 DSSS (no 802.11b PLCP header)-like";
    res.confidence = 0.55;
    why = "Barker-11 spreading at 1 Msym/s but no SFD / valid PLCP header: the preamble was lost "
          "(collision, weak) or the framing is not 802.11b";
  }
  const double fr = (*r)["ch_freq_mhz"].get<double>() * 1e6;
  res.status = "partial";
  res.method = kMethodBarker;
  res.steps = steps();
  res.fields = f;
  res.reason = why;
  res.modulation = "DSSS Barker-11 1 Msym/s";
  res.span = std::make_pair(fr - 11e6, fr + 11e6);
  res.evidence = fmt("Barker-11 matched-filter coherence %.2f (noise ~0.05) over %.0f us on ch%s; ",
                     (*r)["coherence"].get<double>(), dur * 1e6, pstr((*r)["channel"]).c_str()) +
                 why;
  return res;
}

Summary no_dsss() {
  Summary res;
  res.status = "not_possible";
  res.method = kMethodBarker;
  res.steps = steps();
  res.reason = "no Barker-11 spreading (coherence < 0.16) and no PLCP header in this signal's span";
  return res;
}

}  // namespace si::rx::wifi::dsss
