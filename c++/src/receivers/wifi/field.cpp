// python receivers/wifi/field.py (see there for the method and the sources).
#include "receivers/wifi/field.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <string>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::field {

namespace {

using detail::jget;
using detail::truthy;

constexpr int64_t kBigBytes = 600;

bool group(const Json& a) { return !a.is_string() || (std::stoi(a.get<std::string>().substr(0, 2), nullptr, 16) & 1); }

struct Side {
  int64_t ppdus = 0, data_ppdus = 0, bytes = 0, big = 0;
  double airtime_us = 0, data_airtime_us = 0;
  std::vector<double> sizes;
};

struct Info {
  Side down, up;
  int64_t beacons = 0, video = 0;
  std::vector<Json> drone;
  std::vector<double> t;
};

Json side_json(const Side& s) {
  Json o = Json::object();
  o["ppdus"] = s.ppdus;
  o["data_ppdus"] = s.data_ppdus;
  o["airtime_ms"] = py_round(s.airtime_us / 1e3, 3);
  o["data_airtime_ms"] = py_round(s.data_airtime_us / 1e3, 3);
  o["bytes"] = s.bytes;
  o["big_ppdus"] = s.big;
  o["median_bytes"] = s.sizes.empty() ? Json(nullptr) : Json(static_cast<int64_t>(np::median(s.sizes)));
  return o;
}

// python verdict(row, drone, block_s)
Json verdict(const Json& row, const std::vector<Json>& drone, double block_s) {
  const Json &d = row["downlink"], &u = row["uplink"];
  Json beh = Json::object();
  beh["one_client"] = row["n_clients"].get<int64_t>() == 1;
  beh["downlink_dominant"] = d["data_ppdus"].get<int64_t>() >= 5 &&
                             d["data_airtime_ms"].get<double>() >= 3 * u["data_airtime_ms"].get<double>();
  beh["sustained_downlink"] = d["data_airtime_ms"].get<double>() >= 0.05 * block_s * 1e3 ||
                              row["downlink_mbps"].get<double>() >= 1.0;
  beh["uplink_small"] = u["median_bytes"].is_null() || u["median_bytes"].get<int64_t>() < 200;
  beh["video_payload"] = row["video_frames"].get<int64_t>() > 0;
  if (!drone.empty()) {
    const Json* best = &drone[0];
    auto conf = [](const Json& x) { return x.contains("confidence") && x["confidence"].is_number() ? x["confidence"].get<double>() : 0.0; };
    for (const auto& x : drone)
      if (conf(x) > conf(*best)) best = &x;
    Json v = Json::object();
    v["verdict"] = "drone link";
    v["platform"] = (*best)["platform"];
    v["basis"] = "identity";
    v["source"] = jget(*best, "source");
    v["confidence"] = jget(*best, "confidence");
    v["evidence"] = jget(*best, "evidence");
    v["behaviour"] = beh;
    return v;
  }
  if (beh["one_client"].get<bool>() && beh["downlink_dominant"].get<bool>() && beh["sustained_downlink"].get<bool>() &&
      beh["uplink_small"].get<bool>()) {
    Json v = Json::object();
    v["verdict"] = "drone / FPV video link-like";
    v["platform"] = nullptr;
    v["basis"] = "behaviour";
    v["confidence"] = beh["video_payload"].get<bool>() ? 0.65 : 0.55;
    v["behaviour"] = beh;
    v["caveat"] = "one block shows a single client receiving a steady heavy downlink; an IP camera or a "
                  "phone hotspot streaming video looks the same: confirm over ~1 s of blocks or by identity";
    return v;
  }
  return nullptr;
}

}  // namespace

Json links(const std::vector<Json>& frames, double block_s) {
  std::vector<std::string> order;                       // python dict insertion order of aps
  std::map<std::string, Json> aps;
  std::map<std::string, std::set<std::string>> clients;
  std::vector<std::string> client_order;                // insertion order of clients' BSSs (python dict)
  for (const auto& f : frames) {
    for (const auto& m : f["macs"]) {
      const Json& bj = jget(m, "bssid");
      if (bj.is_null() || group(bj)) continue;
      const std::string b = bj.get<std::string>();
      if (!aps.count(b)) {
        aps[b] = Json{{"bssid", b}};
        order.push_back(b);
      }
      const Json& st = jget(m, "subtype");
      if (st == "beacon" || st == "probe-resp")
        for (const char* k : {"ssid", "security", "phy"})
          if (!jget(m, k).is_null() && !aps[b].contains(k)) aps[b][k] = m[k];
      std::vector<Json> peers;
      const bool to_ds = truthy(jget(m, "to_ds")), from_ds = truthy(jget(m, "from_ds"));
      if (jget(m, "type") == "data" && !(to_ds && from_ds)) {
        if (to_ds) peers = {jget(m, "ta")};
        else if (from_ds) peers = {jget(m, "ra")};
        else peers = {jget(m, "ta"), jget(m, "ra")};
      } else if (st == "assoc-req" || st == "reassoc-req" || st == "assoc-resp" || st == "reassoc-resp" || st == "auth") {
        peers = {jget(m, "ta"), jget(m, "ra")};
      }
      for (const auto& a : peers) {
        if (!truthy(a) || a == b || group(a)) continue;
        if (!clients.count(b)) client_order.push_back(b);
        clients[b].insert(a.get<std::string>());
      }
    }
  }
  if (aps.empty()) return Json::array();
  std::map<std::string, std::string> owner;            // "" = more than one BSS
  for (const auto& b : client_order)
    for (const auto& c : clients[b]) owner[c] = owner.count(c) ? std::string() : b;
  std::map<std::string, Info> info;
  for (const auto& b : order) info[b];
  for (const auto& f : frames) {
    const Json& tj = jget(f, "ta");
    std::string b;
    bool down = false;
    if (tj.is_string() && aps.count(tj.get<std::string>())) {
      b = tj.get<std::string>();
      down = true;
    } else if (tj.is_string()) {
      auto it = owner.find(tj.get<std::string>());
      if (it != owner.end() && !it->second.empty()) b = it->second;
    }
    if (b.empty()) continue;
    Info& r = info[b];
    bool beacon = false;
    for (const auto& m : f["macs"]) beacon = beacon || jget(m, "subtype") == "beacon";
    if (beacon) {
      ++r.beacons;
      continue;
    }
    Side& s = down ? r.down : r.up;
    const double dur = f["dur_us"].get<double>();
    ++s.ppdus;
    s.airtime_us += dur;
    if (!truthy(jget(f, "ctrl"))) {
      ++s.data_ppdus;
      s.data_airtime_us += dur;
      if (!jget(f, "bytes").is_null()) {
        const int64_t n = f["bytes"].get<int64_t>();
        s.bytes += n;
        s.sizes.push_back(static_cast<double>(n));
        s.big += n >= kBigBytes;
      }
      if (down) r.t.push_back(f["t0"].get<double>());
    }
    r.video += truthy(jget(f, "video"));
    for (const auto& m : f["macs"])
      if (truthy(jget(m, "drone"))) r.drone.push_back(m["drone"]);
  }
  std::vector<Json> out;
  for (const auto& b : order) {
    const Info& r = info[b];
    if (!(r.down.ppdus || r.up.ppdus)) continue;
    Json row = aps[b];
    const auto& cs = clients[b];
    row["clients"] = Json(std::vector<std::string>(cs.begin(), cs.end()));
    row["n_clients"] = static_cast<int64_t>(cs.size());
    row["beacons"] = r.beacons;
    row["downlink"] = side_json(r.down);
    row["uplink"] = side_json(r.up);
    row["data_airtime_pct"] = py_round(100 * (r.down.data_airtime_us + r.up.data_airtime_us) / (block_s * 1e6), 2);
    row["downlink_mbps"] = py_round(static_cast<double>(r.down.bytes) * 8 / block_s / 1e6, 3);
    if (r.t.size() >= 3) {
      std::vector<double> t = r.t;
      std::sort(t.begin(), t.end());
      std::vector<double> g(t.size() - 1);
      for (size_t i = 1; i < t.size(); ++i) g[i - 1] = (t[i] - t[i - 1]) * 1e3;
      row["downlink_gap_ms"] = {{"median", py_round(np::median(g), 3)},
                                {"max", py_round(*std::max_element(g.begin(), g.end()), 3)}};
    }
    row["video_frames"] = r.video;
    const Json v = verdict(row, r.drone, block_s);
    if (!v.is_null()) row["drone_link"] = v;
    out.push_back(std::move(row));
  }
  auto key = [](const Json& x) { return x["downlink"]["airtime_ms"].get<double>() + x["uplink"]["airtime_ms"].get<double>(); };
  std::stable_sort(out.begin(), out.end(), [&](const Json& a, const Json& c) { return key(a) > key(c); });
  return Json(out);
}

}  // namespace si::rx::wifi::field
