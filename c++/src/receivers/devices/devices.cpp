// python: receivers/devices/plugin.py (_frames_of, DeviceTracking.postprocess / finalize_report)
#include "receivers/devices/devices.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

#include "receivers/devices/pyjson.hpp"

namespace si::rx::devices {

using namespace detail;

std::vector<Frame> frames_of(const Json& fields, const std::string& source) {
  std::vector<Frame> out;
  const Json ppdus = py_or(get(fields, "ppdus"), Json::array());   // PPDU_KEYS = ("ppdus",)
  for (const auto& p : ppdus) {
    // only FCS-verified frames become devices
    const Json ok = get(p, "fcs_ok");
    if (ok.is_boolean() && !ok.get<bool>()) continue;
    Json macs = get(p, "macs");
    if (!truthy(macs)) macs = truthy(get(p, "mac")) ? Json::array({p["mac"]}) : Json::array();
    if (macs.empty()) continue;
    double t0, dur;
    Json cfo, phy;
    if (p.contains("start_s")) {                                    // 802.11b DSSS records
      t0 = num(p["start_s"]);
      const double end = truthy(get(p, "end_s")) ? num(p["end_s"]) : t0;
      dur = (end - t0) > 0.0 ? end - t0 : 0.0;
      cfo = get(p, "cfo_hz");
      phy = "802.11b";
    } else {                                                        // OFDM records
      const Json su = get(p, "start_us"), du = get(p, "dur_us"), ck = get(p, "cfo_khz");
      t0 = (truthy(su) ? num(su) : 0.0) / 1e6;
      dur = (truthy(du) ? num(du) : 0.0) / 1e6;
      cfo = (truthy(ck) ? num(ck) : 0.0) * 1e3;
      phy = get(p, "fmt");
    }
    const double share = dur / double(macs.size());
    const Json rate = py_or(get(p, "rate_mbps"), get(p, "ht_rate_mbps"));
    for (const auto& m : macs)
      out.push_back(frame(m, t0, share, get(p, "snr_db"), cfo, get(p, "channel"), phy, rate, source));
  }
  return out;
}

namespace {

// python: (dec.get("method") or "").split()[0] or "wifi"
std::string source_name(const Decoded& dec) {
  const std::string m = dec.method.value_or("");
  const auto ws = " \t\n\r\v\f";
  const size_t a = m.find_first_not_of(ws);
  if (a == std::string::npos) throw std::runtime_error("IndexError: list index out of range");
  const std::string tok = m.substr(a, m.find_first_of(ws, a) - a);
  return tok.empty() ? "wifi" : tok;
}

std::string tx_of(const Json& mac) {
  const Json tx = py_or(py_or(get(mac, "ta"), get(mac, "sa")), get(mac, "addr2"));
  return tx.is_string() ? tx.get<std::string>() : std::string();
}

}  // namespace

Json device_report(std::vector<Signal>& signals) {
  std::vector<Frame> frames;
  std::vector<std::vector<Frame>> per_signal;
  for (const auto& s : signals) {
    per_signal.push_back(frames_of(s.decoded.fields, source_name(s.decoded)));
    frames.insert(frames.end(), per_signal.back().begin(), per_signal.back().end());
  }
  if (frames.empty()) return Json::object();

  const Json devices = tracks(frames);
  const Json groups = link_randomized(devices);
  Json dlinks = Json::array();  // python: the signals' fields.drone_links (PLAN W1)
  for (const auto& s : signals) {
    const Json dl = get(s.decoded.fields, "drone_links");
    if (dl.is_array())
      for (const auto& d : dl) dlinks.push_back(d);
  }
  const Json found = findings(devices, frames, groups, dlinks);

  std::set<std::string> by_mac;
  for (const auto& d : devices) by_mac.insert(d["mac"].get<std::string>());
  for (size_t i = 0; i < signals.size(); ++i) {
    std::set<std::string> macs;
    for (const auto& f : per_signal[i]) {
      const std::string tx = tx_of(f.mac);
      if (!tx.empty() && by_mac.count(tx)) macs.insert(tx);
    }
    if (macs.empty()) continue;
    Json& fields = signals[i].decoded.fields;
    if (!fields.is_object()) fields = Json::object();
    fields["device_macs"] = str_array({macs.begin(), macs.end()});
  }

  int64_t n_ap = 0, n_rnd = 0;
  for (const auto& d : devices) {
    n_ap += get(d, "role") == Json("access point");
    n_rnd += str_startswith(get(d, "address"), "randomized");
  }
  Json sev = Json::object();
  for (const char* k : {"critical", "high", "medium", "low", "info"}) {
    int64_t n = 0;
    for (const auto& f : found) n += f["severity"] == Json(k);
    if (n) sev[k] = n;
  }
  Json out = Json::object();
  out["devices"] = devices;
  out["findings"] = found;
  out["device_summary"] = {{"n_devices", devices.size()},
                           {"n_access_points", n_ap},
                           {"n_stations", int64_t(devices.size()) - n_ap},
                           {"n_randomized", n_rnd},
                           {"n_frames", frames.size()},
                           {"method", "device tracking over decoded 802.11 frames (receivers.devices)"},
                           {"severity_counts", sev}};
  if (!groups.empty()) out["device_groups"] = groups;
  return out;
}

void finalize(Json& rep) {
  std::map<std::string, Json> where;   // mac -> signal ids, in signal order
  for (const auto& s : py_or(get(rep, "signals"), Json::array())) {
    const Json fields = py_or(get(py_or(get(s, "decoded"), Json::object()), "fields"), Json::object());
    for (const auto& mac : py_or(get(fields, "device_macs"), Json::array())) {
      Json& v = where[mac.get<std::string>()];
      if (v.is_null()) v = Json::array();
      v.push_back(s["id"]);
    }
  }
  if (rep.contains("devices"))
    for (auto& d : rep["devices"]) {
      auto it = where.find(d["mac"].get<std::string>());
      if (it != where.end()) d["signal_ids"] = it->second;
    }
  if (rep.contains("findings"))
    for (auto& f : rep["findings"]) {
      std::set<int64_t> ids;
      for (const auto& m : py_or(get(f, "devices"), Json::array())) {
        auto it = where.find(m.get<std::string>());
        if (it != where.end())
          for (const auto& i : it->second) ids.insert(i.get<int64_t>());
      }
      if (!ids.empty()) {
        Json a = Json::array();
        for (auto i : ids) a.push_back(i);
        f["signal_ids"] = a;
      }
    }
}

}  // namespace si::rx::devices
