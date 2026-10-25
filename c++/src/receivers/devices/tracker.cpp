// python: receivers/devices/tracker.py
#include "receivers/devices/tracker.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <set>

#include "receivers/devices/pyjson.hpp"

namespace si::rx::devices {
namespace {

using namespace detail;

const std::set<std::string> kApSubtypes = {"beacon", "probe-resp", "assoc-resp", "reassoc-resp"};

bool in_ap_subtypes(const Json& st) { return st.is_string() && kApSubtypes.count(st.get<std::string>()) > 0; }

int first_octet(const std::string& mac) { return std::stoi(mac.substr(0, 2), nullptr, 16); }

bool is_group(const Json& mac) { return truthy(mac) && (first_octet(mac.get<std::string>()) & 0x01); }

bool is_random(const std::string& mac) {
  if (mac.empty()) return false;
  const int b = first_octet(mac);
  return !(b & 0x01) && (b & 0x02);
}

// python _probe_fingerprint
std::optional<std::string> probe_fingerprint(const Json& m) {
  if (get(m, "subtype") != Json("probe-req")) return std::nullopt;
  const Json ids = py_or(get(m, "element_ids"), Json::array());
  if (!truthy(ids)) return std::nullopt;
  const Json ext = py_or(get(m, "ext_element_ids"), Json::array());
  const Json vend = py_or(get(m, "vendor_elements"), Json::array());
  auto join = [](const Json& a) {
    std::string s;
    for (size_t i = 0; i < a.size(); ++i) s += (i ? "," : "") + py_str(a[i]);
    return s;
  };
  std::string out = "ie:" + join(ids);
  if (truthy(ext)) out += "|ext:" + join(ext);
  if (truthy(vend)) {
    std::vector<std::string> v;
    for (const auto& x : vend) v.push_back(py_str(x));
    out += "|vnd:" + join(str_array(sorted_unique(v)));
  }
  const Json ht = py_or(get(m, "ht"), Json::object());
  if (truthy(ht))
    out += "|ht:" + py_str(get(ht, "spatial_streams")) + "ss" + (truthy(get(ht, "40mhz")) ? "/40" : "");
  if (truthy(get(m, "rates_mbps"))) out += "|rates:" + join(m["rates_mbps"]);
  return out;
}

// python _median: sorted non-None values, middle one or the mean of the two middle ones
Json median(const std::vector<Json>& in) {
  std::vector<Json> xs;
  for (const auto& x : in)
    if (!x.is_null()) xs.push_back(x);
  if (xs.empty()) return nullptr;
  std::stable_sort(xs.begin(), xs.end(), [](const Json& a, const Json& b) { return num(a) < num(b); });
  const size_t n = xs.size();
  if (n % 2) return xs[n / 2];
  return (num(xs[n / 2 - 1]) + num(xs[n / 2])) / 2;
}

// python _round: round(v, n) keeps an int an int
Json round_json(const Json& v, int n) {
  if (v.is_null() || !v.is_number_float()) return v;
  return py_round(v.get<double>(), n);
}

// a + b with Python int/float promotion
Json py_add(const Json& a, const Json& b) {
  if (a.is_number_integer() && b.is_number_integer()) return a.get<int64_t>() + b.get<int64_t>();
  return num(a) + num(b);
}

struct Device {
  std::string mac;
  int64_t n_frames = 0;
  std::optional<double> t0, t1;
  std::vector<Json> snr, cfo;
  Counter subtypes, channels, phys, peers, fingerprints, rates, sources;
  std::set<std::string> probed;
  double airtime = 0.0;
  Json bytes = 0;
  Json vendor = nullptr;
  bool is_ap = false;
  Json ssid = nullptr;
  bool ssid_hidden = false;
  Json bssid = nullptr;
  Json security = nullptr, phy = nullptr, uptime = nullptr, country = nullptr, bss_load = nullptr, wps = nullptr;
  Json device_name = nullptr, drone = nullptr, remote_id = nullptr;
  Json pn_max = nullptr;
  int64_t retries = 0, pwr_save = 0;
  std::optional<std::pair<double, Json>> seq_first, seq_last;  // (time, sequence number) of the first / last frame

  void see(const Frame& f, const Json& m, bool role_ap) {
    ++n_frames;
    const double t = f.t0;
    t0 = t0 ? std::min(*t0, t) : t;
    const double a = (t1 && *t1 != 0.0) ? *t1 : t;  // self.t1 or t0
    const double b = t + f.dur_s;
    t1 = a >= b ? a : b;
    airtime += f.dur_s;
    if (truthy(f.length_bytes)) bytes = py_add(bytes, f.length_bytes);
    if (!f.snr_db.is_null()) snr.push_back(f.snr_db);
    if (!f.cfo_hz.is_null()) cfo.push_back(f.cfo_hz);
    if (!f.channel.is_null()) channels.add(f.channel);
    subtypes.add(m.contains("subtype") ? m["subtype"] : Json("?"));
    if (!f.source.empty()) sources.add(f.source);
    if (truthy(f.rate_mbps)) rates.add(f.rate_mbps);
    if (const Json sq = get(m, "seq"); !sq.is_null()) {
      if (!seq_first || t < seq_first->first) seq_first = std::make_pair(t, sq);
      if (!seq_last || t >= seq_last->first) seq_last = std::make_pair(t, sq);
    }
    if (truthy(get(m, "retry"))) ++retries;
    if (truthy(get(m, "pwr_mgmt"))) ++pwr_save;
    is_ap = is_ap || role_ap;
    if (!truthy(vendor)) vendor = py_or(py_or(get(m, "ta_vendor"), get(m, "sa_vendor")), get(m, "bssid_vendor"));
    if (auto fp = probe_fingerprint(m)) fingerprints.add(*fp);
    const Json ss = get(m, "ssid");
    if (get(m, "subtype") == Json("probe-req") && !ss.is_null() && ss != Json("<broadcast>"))
      probed.insert(py_str(ss));
    const Json name = py_or(get(py_or(get(m, "p2p"), Json::object()), "device_name"),
                            get(py_or(get(m, "wps"), Json::object()), "device_name"));
    if (truthy(name)) device_name = name;
    const Json dr = get(m, "drone");
    if (truthy(dr) && (drone.is_null() || num(dr["confidence"]) > num(drone["confidence"]))) drone = dr;
    if (truthy(get(m, "remote_id"))) remote_id = m["remote_id"];
    if (truthy(get(m, "phy"))) phys.add(m["phy"]);
    const Json enc = py_or(get(m, "encryption"), Json::object());
    const Json pn = get(enc, "pn");
    if (!pn.is_null()) {
      const Json cur = truthy(pn_max) ? pn_max : Json(0);
      pn_max = num(cur) >= num(pn) ? cur : pn;
    }
    if (role_ap) {
      const std::pair<const char*, Json*> keys[] = {{"ssid", &ssid},         {"security", &security},
                                                    {"phy", &phy},           {"ap_uptime", &uptime},
                                                    {"country", &country},   {"bss_load", &bss_load},
                                                    {"wps", &wps}};
      for (const auto& [k, dst] : keys) {
        const Json v = get(m, k);
        if (!v.is_null()) *dst = v;
      }
      if (truthy(get(m, "ssid_hidden"))) ssid_hidden = true;
    }
  }

  Json to_dict() const {
    Json d = Json::object();
    d["mac"] = mac;
    d["vendor"] = vendor;
    d["address"] = is_random(mac) ? "randomized (locally administered)" : "global (OUI)";
    d["role"] = is_ap ? "access point" : "station";
    d["n_frames"] = n_frames;
    d["first_seen_s"] = t0 ? Json(py_round(*t0, 6)) : Json(nullptr);
    d["last_seen_s"] = t1 ? Json(py_round(*t1, 6)) : Json(nullptr);
    d["airtime_ms"] = py_round(airtime * 1e3, 3);
    d["bytes_seen"] = truthy(bytes) ? bytes : Json(nullptr);
    d["snr_db_median"] = round_json(median(snr), 1);
    Json mx = nullptr;
    for (const auto& x : snr)
      if (mx.is_null() || num(x) > num(mx)) mx = x;
    d["snr_db_max"] = round_json(mx, 1);
    d["cfo_hz_median"] = round_json(median(cfo), 0);
    d["channels"] = channels.to_json();
    d["frame_types"] = subtypes.to_json();
    d["phy"] = truthy(phy) ? phy : (phys.empty() ? Json(nullptr) : phys.most_common(1)[0].first);
    if (is_ap) {
      d["ssid"] = (ssid_hidden && !truthy(ssid)) ? Json("<hidden>") : ssid;
      d["ssid_hidden"] = ssid_hidden ? Json(true) : Json(nullptr);
      d["security"] = security;
      d["ap_uptime"] = uptime;
      d["country"] = get(py_or(country, Json::object()), "code");
      d["bss_load"] = bss_load;
      d["wps"] = wps;
    } else {
      d["associated_bssid"] = bssid;
      d["probed_ssids"] = probed.empty() ? Json(nullptr) : str_array({probed.begin(), probed.end()});
    }
    if (!rates.empty()) d["rates_mbps"] = rates.to_json();
    if (retries) d["retry_fraction"] = py_round(double(retries) / double(n_frames), 3);
    if (pwr_save) d["power_save_fraction"] = py_round(double(pwr_save) / double(n_frames), 3);
    if (!pn_max.is_null()) d["encryption_packet_counter_max"] = pn_max;
    if (truthy(device_name)) d["device_name"] = device_name;
    if (!peers.empty()) {
      Json t = Json::array();
      for (const auto& e : peers.most_common(8)) t.push_back(e.first);
      d["talks_to"] = t;
    }
    if (!fingerprints.empty()) d["probe_fingerprint"] = fingerprints.most_common(1)[0].first;
    if (seq_first) {
      d["seq_first"] = seq_first->second;
      d["seq_last"] = seq_last->second;
    }
    if (truthy(drone)) d["drone"] = drone;
    if (truthy(remote_id)) d["remote_id"] = remote_id;
    if (!sources.empty()) {
      std::vector<std::string> s;
      for (const auto& e : sources.items()) s.push_back(e.first.get<std::string>());
      d["decoded_by"] = str_array(sorted_unique(s));
    }
    Json out = Json::object();
    for (auto it = d.begin(); it != d.end(); ++it)
      if (!it.value().is_null()) out[it.key()] = it.value();
    return out;
  }
};

}  // namespace

Frame frame(const Json& mac, double t0, double dur_s, Json snr_db, Json cfo_hz, Json channel, Json phy,
            Json rate_mbps, std::string source, Json length_bytes) {
  Frame f;
  f.mac = mac;
  f.t0 = t0;
  f.dur_s = dur_s;
  f.snr_db = std::move(snr_db);
  f.cfo_hz = std::move(cfo_hz);
  f.channel = std::move(channel);
  f.phy = std::move(phy);
  f.rate_mbps = std::move(rate_mbps);
  f.length_bytes = !length_bytes.is_null() ? std::move(length_bytes) : get(mac, "length_bytes");
  f.source = std::move(source);
  return f;
}

Json tracks(const std::vector<Frame>& frames) {
  std::vector<Device> devs;
  std::map<std::string, size_t> idx;
  for (const auto& f : frames) {
    const Json m = py_or(f.mac, Json::object());
    const Json st = get(m, "subtype");
    Json tx = py_or(py_or(get(m, "ta"), get(m, "sa")), get(m, "addr2"));
    const Json bssid = get(m, "bssid");
    if (in_ap_subtypes(st) && truthy(bssid) && !is_group(bssid)) tx = py_or(tx, bssid);
    if (!truthy(tx) || is_group(tx)) continue;
    const std::string txs = tx.get<std::string>();
    auto [it, fresh] = idx.emplace(txs, devs.size());
    if (fresh) {
      devs.emplace_back();
      devs.back().mac = txs;
    }
    Device& d = devs[it->second];
    d.see(f, m, in_ap_subtypes(st) && bssid == tx);
    const Json rx = py_or(py_or(get(m, "da"), get(m, "ra")), get(m, "addr1"));
    if (truthy(rx) && !is_group(rx)) d.peers.add(rx);
    if (truthy(bssid) && !is_group(bssid) && bssid != tx && get(m, "type") == Json("data")) d.bssid = bssid;
  }

  std::set<std::string> aps;
  for (const auto& d : devs)
    if (d.is_ap) aps.insert(d.mac);
  for (auto& d : devs) {
    if (d.bssid.is_null() && !d.is_ap)
      for (const auto& [p, n] : d.peers.most_common())
        if (p.is_string() && aps.count(p.get<std::string>())) {
          d.bssid = p;
          break;
        }
  }

  std::vector<Json> out;
  for (const auto& d : devs) out.push_back(d.to_dict());
  auto snr_key = [](const Json& d) {
    const Json v = get(d, "snr_db_median");
    return -(v.is_null() ? -99.0 : num(v));
  };
  std::stable_sort(out.begin(), out.end(), [&](const Json& a, const Json& b) {
    const double ka = snr_key(a), kb = snr_key(b);
    if (ka != kb) return ka < kb;
    return -a["n_frames"].get<int64_t>() < -b["n_frames"].get<int64_t>();
  });
  return Json(out);
}

constexpr int64_t kSeqGapMax = 64;  // python SEQ_GAP_MAX

Json link_randomized(const Json& devices, double cfo_tol_hz) {
  struct Group {
    std::vector<std::string> macs;
    Json fp, cfo;
    int64_t n_frames = 0;
  };
  std::vector<Group> groups;
  for (const auto& d : devices) {
    if (!str_startswith(get(d, "address"), "randomized") || !truthy(get(d, "probe_fingerprint"))) continue;
    bool placed = false;
    const Json dc = get(d, "cfo_hz_median");
    for (auto& g : groups) {
      if (g.fp != d["probe_fingerprint"]) continue;
      if (!dc.is_null() && !g.cfo.is_null() && std::fabs(num(dc) - num(g.cfo)) > cfo_tol_hz) continue;
      g.macs.push_back(d["mac"].get<std::string>());
      g.n_frames += d["n_frames"].get<int64_t>();
      placed = true;
      break;
    }
    if (!placed) groups.push_back({{d["mac"].get<std::string>()}, d["probe_fingerprint"], dc,
                                   d["n_frames"].get<int64_t>()});
  }
  std::map<std::string, const Json*> by_mac;
  for (const auto& d : devices)
    if (str_startswith(get(d, "address"), "randomized") && truthy(get(d, "probe_fingerprint")))
      by_mac[d["mac"].get<std::string>()] = &d;
  std::vector<Json> out;
  for (const auto& g : groups) {
    if (g.macs.size() < 2) continue;
    const bool same_cfo = !g.cfo.is_null();
    // research, PLAN W7: the next address continuing the previous one's sequence counter
    std::vector<std::string> order = g.macs;
    auto first = [&](const std::string& a) {
      const Json v = get(*by_mac[a], "first_seen_s");
      return truthy(v) ? num(v) : 0.0;
    };
    std::sort(order.begin(), order.end(), [&](const std::string& a, const std::string& b) {
      return first(a) != first(b) ? first(a) < first(b) : a < b;
    });
    bool seq = false;
    for (size_t i = 0; i + 1 < order.size(); ++i) {
      const Json sl = get(*by_mac[order[i]], "seq_last"), sf = get(*by_mac[order[i + 1]], "seq_first");
      if (sl.is_null() || sf.is_null()) continue;
      const int64_t gap = ((sf.get<int64_t>() - sl.get<int64_t>()) % 4096 + 4096) % 4096;
      seq = seq || (gap > 0 && gap <= kSeqGapMax);
    }
    std::vector<std::string> macs = g.macs;
    std::sort(macs.begin(), macs.end());
    Json o = Json::object();
    o["macs"] = str_array(macs);
    o["n_addresses"] = g.macs.size();
    o["n_frames"] = g.n_frames;
    o["probe_fingerprint"] = g.fp;
    o["cfo_hz_median"] = g.cfo;
    o["sequence_continues"] = seq;
    o["confidence"] = seq ? 0.9 : same_cfo ? 0.8 : 0.6;
    o["evidence"] = "identical probe-request element signature" +
                    (same_cfo ? fmt(" and transmitter clock offset within %.0f Hz", cfo_tol_hz) : std::string()) +
                    (seq ? " and one sequence-number counter continuing across the address change" : "") +
                    " -- one device rotating its MAC address";
    out.push_back(o);
  }
  std::stable_sort(out.begin(), out.end(), [](const Json& a, const Json& b) {
    return -a["n_frames"].get<int64_t>() < -b["n_frames"].get<int64_t>();
  });
  return Json(out);
}

}  // namespace si::rx::devices
