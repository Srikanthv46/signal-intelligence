// python: receivers/devices/threats.py
#include "receivers/devices/threats.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "receivers/devices/pyjson.hpp"

namespace si::rx::devices {
namespace {

using namespace detail;
using Extras = std::vector<std::pair<std::string, Json>>;

constexpr int64_t kDeauthSuspect = 8;
constexpr int64_t kDeauthFlood = 25;
const char* const kWeakCiphers[] = {"WEP", "TKIP", "WPA1"};

int severity_order(const Json& s) {
  static const std::map<std::string, int> o = {{"critical", 0}, {"high", 1}, {"medium", 2}, {"low", 3}, {"info", 4}};
  if (!s.is_string()) return 9;
  auto it = o.find(s.get<std::string>());
  return it == o.end() ? 9 : it->second;
}

bool is_ap(const Json& d) { return get(d, "role") == Json("access point"); }

std::string join(const std::vector<std::string>& v, const std::string& sep) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) s += (i ? sep : "") + v[i];
  return s;
}

// python _f
Json finding(const char* kind, const char* severity, const std::string& title, const std::string& detail,
             const Json& confidence, const Json& evidence, const std::vector<std::string>& devices = {},
             const Extras& extra = {}) {
  Json d = Json::object();
  d["kind"] = kind;
  d["severity"] = severity;
  d["title"] = title;
  d["detail"] = detail;
  d["confidence"] = confidence;
  d["evidence"] = evidence;
  if (!devices.empty()) d["devices"] = str_array(sorted_unique(devices));
  for (const auto& [k, v] : extra) d[k] = v;
  return d;
}

// python _norm: MAC with the locally-administered bit cleared
std::string norm(const std::string& mac) {
  return fmt("%02x", std::stoi(mac.substr(0, 2), nullptr, 16) & ~0x02) + mac.substr(2);
}

std::vector<std::string> split(const std::string& s, char c) {
  std::vector<std::string> out;
  size_t a = 0;
  for (size_t b; (b = s.find(c, a)) != std::string::npos; a = b + 1) out.push_back(s.substr(a, b - a));
  out.push_back(s.substr(a));
  return out;
}

// python _same_radio
bool same_radio(const std::vector<std::string>& macs) {
  std::vector<std::string> n;
  std::set<std::string> oui;
  for (const auto& m : macs) {
    n.push_back(norm(m));
    oui.insert(n.back().substr(0, 8));
  }
  if (oui.size() != 1) return false;
  std::vector<std::vector<std::string>> tails;
  for (const auto& x : n) tails.push_back(split(x.size() > 9 ? x.substr(9) : std::string(), ':'));
  for (size_t i = 1; i < tails.size(); ++i) {
    int diff = 0;
    for (size_t j = 0; j < std::min(tails[0].size(), tails[i].size()); ++j) diff += tails[0][j] != tails[i][j];
    if (diff > 1) return false;
  }
  return true;
}

Json drones(const Json& devices) {
  Json out = Json::array();
  for (const auto& d : devices) {
    const Json dr = get(d, "drone");
    if (!truthy(dr)) continue;
    const Json rid = py_or(get(d, "remote_id"), Json::object());
    std::string where, op;
    if (!get(rid, "lat").is_null()) {
      where = fmt(" at %.6f,%.6f", num(rid["lat"]), num(rid["lon"]));
      if (!get(rid, "height_m").is_null()) where += fmt(", %.0f m above takeoff", num(rid["height_m"]));
    }
    if (!get(rid, "operator_lat").is_null())
      op = fmt(" Operator at %.6f,%.6f.", num(rid["operator_lat"]), num(rid["operator_lon"]));
    const Json evidence = dr.contains("evidence") ? dr["evidence"] : Json("");
    const std::vector<std::string> dev = {d["mac"].get<std::string>()};
    if (truthy(rid)) {
      const Json oploc = get(rid, "operator_lat").is_null() ? Json(nullptr)
                                                            : Json::array({rid["operator_lat"], rid["operator_lon"]});
      out.push_back(finding("drone_remote_id", "high", "Drone broadcasting Remote ID: " + py_str(dr["platform"]),
                            "An unmanned aircraft is identifying itself" + where + "." + op +
                                " Remote ID is transmitted in clear by regulation, so these are the "
                                "aircraft's own declared values.",
                            0.99, evidence, dev,
                            {{"remote_id", rid}, {"uas_id", get(rid, "uas_id")}, {"operator_location", oploc}}));
    } else {
      out.push_back(finding("drone_signature", "medium", "Possible drone link: " + py_str(dr["platform"]),
                            "A transmitter matches a drone platform by manufacturer and/or network name. "
                            "No Remote ID was broadcast, so this identifies the equipment, not a flight.",
                            dr.contains("confidence") ? dr["confidence"] : Json(0.7), evidence, dev));
    }
  }
  return out;
}

Json deauth(const std::vector<Frame>& frames) {
  Json out = Json::array();
  Counter by_tx;
  std::map<std::string, std::set<std::string>> targets;
  for (const auto& fr : frames) {
    const Json m = py_or(fr.mac, Json::object());
    const Json st = get(m, "subtype");
    if (st != Json("deauth") && st != Json("disassoc")) continue;
    const Json tx = py_or(py_or(get(m, "ta"), get(m, "sa")), get(m, "addr2"));
    if (!truthy(tx)) continue;
    by_tx.add(tx);
    const Json rx = py_or(py_or(get(m, "da"), get(m, "ra")), get(m, "addr1"));
    if (truthy(rx)) targets[py_str(tx)].insert(py_str(rx));
  }
  for (const auto& [txj, n] : by_tx.items()) {
    if (n < kDeauthSuspect) continue;
    const std::string tx = py_str(txj);
    const auto& tg = targets[tx];
    bool broadcast = false;
    for (std::string t : tg) {
      std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return std::tolower(c); });
      broadcast = broadcast || t == "ff:ff:ff:ff:ff:ff";
    }
    const bool critical = n >= kDeauthFlood || broadcast;
    std::vector<std::string> first16(tg.begin(), tg.end());
    if (first16.size() > 16) first16.resize(16);
    out.push_back(finding(
        "deauth_flood", critical ? "critical" : "high", "Deauthentication burst from " + tx,
        fmt("%lld deauthenticate/disassociate frames in this capture, aimed at %zu address(es)", (long long)n,
            tg.size()) +
            (broadcast ? ", including the broadcast address, which disconnects every client at once" : "") +
            ". This is how a client is forced off a network, either to deny it service or "
            "to make it reconnect so the handshake can be captured. A busy AP shedding "
            "idle clients can also produce these, so check whether one address is the "
            "source of all of them.",
        critical ? 0.85 : 0.75, fmt("%lld deauth/disassoc frames from one transmitter in a single capture", (long long)n),
        {tx}, {{"n_frames", n}, {"targets", str_array(first16)}}));
  }
  return out;
}

Json evil_twin(const Json& devices) {
  Json out = Json::array();
  std::vector<std::pair<std::string, std::vector<Json>>> by_ssid;
  for (const auto& d : devices) {
    const Json ss = get(d, "ssid");
    if (!is_ap(d) || !truthy(ss) || str_startswith(ss, "<")) continue;
    const std::string s = ss.get<std::string>();
    auto it = std::find_if(by_ssid.begin(), by_ssid.end(), [&](const auto& e) { return e.first == s; });
    if (it == by_ssid.end()) {
      by_ssid.push_back({s, {}});
      it = by_ssid.end() - 1;
    }
    it->second.push_back(d);
  }
  for (const auto& [ssid, aps] : by_ssid) {
    if (aps.size() < 2) continue;
    std::vector<Json> secs;
    std::set<std::string> vendors;
    std::vector<std::string> macs;
    bool open_now = false;
    for (const auto& a : aps) {
      const Json s = get(a, "security");
      if (std::none_of(secs.begin(), secs.end(), [&](const Json& x) { return py_eq(x, s); })) secs.push_back(s);
      vendors.insert(split(py_str(py_or(get(a, "vendor"), Json("?"))), ',')[0]);
      macs.push_back(a["mac"].get<std::string>());
      open_now = open_now || s == Json("open") || s == Json("WEP");
    }
    const bool differs = secs.size() > 1;
    if (same_radio(macs)) {
      out.push_back(finding("multi_bssid", "info",
                            "'" + ssid + "' is advertised on " + std::to_string(aps.size()) + " BSSIDs of one radio",
                            "The addresses differ only in the locally-administered bit and a byte, so "
                            "this is one access point running several virtual BSSs (guest network, "
                            "band steering), not a second device.",
                            0.85,
                            "BSSIDs " + join(macs, ", ") +
                                " share an OUI once the locally-administered bit is cleared",
                            macs, {{"ssid", ssid}}));
      continue;
    }
    if (!differs && vendors.size() <= 1) continue;
    const bool high = differs && open_now;
    std::vector<std::string> sec_str;
    for (const auto& s : secs) sec_str.push_back(py_str(s));
    std::sort(sec_str.begin(), sec_str.end());
    Json per = Json::object();
    for (const auto& a : aps) per[a["mac"].get<std::string>()] = get(a, "security");
    out.push_back(finding(
        "duplicate_ssid", high ? "high" : "medium",
        "Network name '" + ssid + "' advertised by " + std::to_string(aps.size()) + " different radios",
        std::string("The same network name is coming from more than one BSSID with ") +
            (differs ? "different security settings" : "different hardware vendors") +
            ". That is what an impersonating access point looks like; it is also what a "
            "multi-AP home mesh or an extender looks like. The telling sign is a copy "
            "offering weaker security than the original." +
            (open_now ? " One copy is open or WEP." : ""),
        high ? 0.7 : 0.55,
        "BSSIDs " + join(macs, ", ") + "; security " + py_repr(str_array(sec_str)) + "; vendors " +
            py_repr(str_array({vendors.begin(), vendors.end()})),
        macs, {{"ssid", ssid}, {"security_per_bssid", per}}));
  }
  return out;
}

Json weak_security(const Json& devices) {
  Json out = Json::array();
  for (const auto& d : devices) {
    if (!is_ap(d)) continue;
    const Json sec = get(d, "security");
    if (sec.is_null()) continue;
    const std::string mac = d["mac"].get<std::string>();
    const Json ssid = get(d, "ssid");
    const std::string ss = py_str(ssid);
    const std::string s = py_str(sec);
    if (sec == Json("open")) {
      out.push_back(finding("open_network", "low", "Open network '" + ss + "'",
                            "No link-layer encryption: everything this network carries is readable "
                            "by anyone in range, and any device may associate.",
                            0.95, "no RSN/WPA element, privacy bit clear (" + mac + ")", {mac}, {{"ssid", ssid}}));
    } else if (sec.is_string() && std::any_of(std::begin(kWeakCiphers), std::end(kWeakCiphers), [&](const char* w) {
                 return s.find(w) != std::string::npos;
               })) {
      const bool strong = s.find("CCMP") != std::string::npos || s.find("GCMP") != std::string::npos;
      if (strong)
        out.push_back(finding("legacy_cipher_offered", "low", "'" + ss + "' still accepts TKIP",
                              "Advertises " + s +
                                  ". Modern clients will use CCMP, but leaving TKIP "
                                  "enabled keeps a weaker option available and blocks 802.11n rates.",
                              0.85, "security = " + s + " (" + mac + ")", {mac}, {{"ssid", ssid}}));
      else
        out.push_back(finding("weak_encryption", "medium", "Obsolete encryption on '" + ss + "'",
                              "Advertises " + s +
                                  " with no modern cipher. WEP and TKIP are broken and "
                                  "WPA1 is withdrawn; this traffic should not be treated as protected.",
                              0.9, "security = " + s + " (" + mac + ")", {mac}, {{"ssid", ssid}}));
    }
    const Json w = py_or(get(d, "wps"), Json::object());
    if (truthy(w) && get(w, "state") == Json("configured") && truthy(get(w, "selected_registrar")))
      out.push_back(finding("wps_active", "medium", "WPS registrar active on '" + ss + "'",
                            "The access point is in WPS enrolment mode, which accepts new devices with "
                            "a PIN or button press. It should not normally be left on.",
                            0.8, "WPS selected registrar set (" + mac + ")", {mac}, {{"ssid", ssid}}));
  }
  return out;
}

Json surveillance_and_identity(const Json& devices) {
  Json out = Json::array();
  for (const auto& d : devices) {
    const Json probed = py_or(get(d, "probed_ssids"), Json::array());
    const std::string mac = d["mac"].get<std::string>();
    if (probed.size() >= 1 && !is_ap(d)) {
      std::vector<std::string> r;
      for (size_t i = 0; i < probed.size() && i < 8; ++i) r.push_back(py_repr(probed[i]));
      out.push_back(finding("client_history", "info", "Client " + mac + " is searching for known networks",
                            "A station lists networks it has joined before when it probes for them. "
                            "That history says where the device has been.",
                            0.9, "probe requests for " + join(r, ", "), {mac},
                            {{"probed_ssids", probed},
                             {"device_name", get(d, "device_name")},
                             {"vendor", get(d, "vendor")}}));
    }
    if (truthy(get(d, "device_name")))
      out.push_back(finding("named_device", "info",
                            "Device identifies itself as '" + py_str(d["device_name"]) + "'",
                            "The device broadcasts its own name (WPS or Wi-Fi Direct), which usually "
                            "gives the model and often the owner.",
                            0.95, mac + " [" + py_str(get(d, "vendor")) + "]", {mac},
                            {{"device_name", d["device_name"]}}));
  }
  return out;
}

Json eapol(const std::vector<Frame>& frames) {
  Json out = Json::array();
  Counter n;
  std::set<std::string> who;
  for (const auto& fr : frames) {
    const Json m = py_or(fr.mac, Json::object());
    const Json ev = get(m, "eapol");
    if (!truthy(ev)) continue;
    const Json h = py_or(py_or(get(ev, "handshake"), get(ev, "eap")), Json("EAPOL"));
    n.add(h);
    for (const char* k : {"ta", "sa", "da", "ra"})
      if (truthy(get(m, k))) who.insert(py_str(m[k]));
  }
  if (!n.empty()) {
    std::vector<std::string> ev;
    for (const auto& [k, v] : n.most_common()) ev.push_back(py_str(k) + " x" + std::to_string(v));
    std::vector<std::string> dev(who.begin(), who.end());
    if (dev.size() > 8) dev.resize(8);
    out.push_back(finding("key_exchange", "info", "WPA key exchange observed",
                          "A client is joining or rekeying. The handshake itself does not expose the "
                          "password, but its presence marks exactly when a device associated.",
                          0.95, join(ev, "; "), dev, {{"messages", n.to_json()}}));
  }
  return out;
}

// python _beacon_flood (research, PLAN W6)
Json beacon_flood(const Json& devices) {
  Json out = Json::array();
  std::vector<std::string> fake;
  std::set<std::string> prefixes;
  for (const auto& d : devices) {
    if (!is_ap(d) || !truthy(get(get(d, "frame_types"), "beacon"))) continue;
    const Json v = get(d, "vendor");
    const std::string vs = truthy(v) ? py_str(v) : "unknown OUI";
    if (vs != "unknown OUI" && vs != "locally-administered (randomized/private)") continue;
    fake.push_back(d["mac"].get<std::string>());
    prefixes.insert(norm(fake.back()).substr(0, 8));
  }
  if (fake.size() < 15 || prefixes.size() < 8) return out;
  std::vector<std::string> shown(fake.begin(), fake.begin() + static_cast<ptrdiff_t>(std::min<size_t>(16, fake.size())));
  out.push_back(finding("beacon_flood", "high", "Beacon flood: " + std::to_string(fake.size()) + " unregistered access points",
                        "Access points with addresses no manufacturer registered, from unrelated address prefixes, all "
                        "beaconing in one capture. That is what a beacon-flooding tool produces; it hides the real "
                        "networks and confuses clients scanning for them.",
                        0.75,
                        std::to_string(fake.size()) + " beaconing BSSIDs without an IEEE OUI, " +
                            std::to_string(prefixes.size()) + " different prefixes",
                        shown, {{"n_access_points", static_cast<int64_t>(fake.size())}}));
  return out;
}

// python _drone_behaviour (research, PLAN W1)
Json drone_behaviour(const Json& links) {
  Json out = Json::array();
  for (const auto& dl : links) {
    if (get(dl, "basis") != Json("behaviour")) continue;
    const std::string b = dl["bssid"].get<std::string>();
    std::vector<std::string> devs = {b}, ev;
    for (const auto& c : get(dl, "clients")) devs.push_back(c.get<std::string>());
    const Json beh = get(dl, "behaviour");
    if (beh.is_object())
      for (auto it = beh.begin(); it != beh.end(); ++it)
        if (truthy(it.value())) ev.push_back(it.key());
    const Json ssid = get(dl, "ssid");
    out.push_back(finding("drone_link_behaviour", "medium",
                          "Possible drone video link: " + b + (truthy(ssid) ? " '" + py_str(ssid) + "'" : ""),
                          "One client receives a steady, heavy downlink from this access point and sends back only "
                          "small frames: how a Wi-Fi drone streams video to the pilot's phone. Nothing identifies the "
                          "equipment, and an IP camera or a phone hotspot streaming video looks the same in one block.",
                          dl.contains("confidence") ? dl["confidence"] : Json(0.55), join(ev, ", "), devs,
                          {{"bssid", b}}));
  }
  return out;
}

Json randomization(const Json& links) {
  Json out = Json::array();
  for (const auto& g : links) {
    std::vector<std::string> macs;
    for (const auto& m : g["macs"]) macs.push_back(m.get<std::string>());
    out.push_back(finding("mac_rotation", "low", py_str(g["n_addresses"]) + " addresses appear to be one device",
                          "The addresses differ but the probe-request signature (and transmitter clock, where "
                          "measured) match, which is what MAC randomisation on one handset looks like. Counting "
                          "these as separate devices would overstate how many are present.",
                          g["confidence"], g["evidence"], macs, {{"cfo_hz_median", get(g, "cfo_hz_median")}}));
  }
  return out;
}

Json inventory(const Json& devices) {
  Json out = Json::array();
  if (devices.empty()) return out;
  size_t n_ap = 0, n_sta = 0;
  std::vector<std::string> named;
  for (const auto& d : devices) {
    if (is_ap(d)) {
      ++n_ap;
      const Json s = get(d, "ssid");
      if (truthy(s) && !str_startswith(s, "<")) named.push_back(s.get<std::string>());
    } else {
      ++n_sta;
    }
  }
  named = sorted_unique(named);
  std::vector<std::string> q;
  for (size_t i = 0; i < named.size() && i < 8; ++i) q.push_back("'" + named[i] + "'");
  std::string ev = q.empty() ? "no named networks" : join(q, ", ");
  if (named.size() > 8) ev += " and " + std::to_string(named.size() - 8) + " more";
  out.push_back(finding("inventory", "info",
                        std::to_string(n_ap) + " access points and " + std::to_string(n_sta) + " client devices heard",
                        "Every transmitter decoded in this capture, with its vendor, signal level and activity.", 1.0,
                        ev, {}, {{"n_access_points", n_ap}, {"n_stations", n_sta}}));
  return out;
}

}  // namespace

Json findings(const Json& devices, const std::vector<Frame>& frames, const Json& links, const Json& drone_links) {
  std::vector<Json> out;
  for (const Json& part : {drones(devices), drone_behaviour(drone_links.is_array() ? drone_links : Json::array()),
                           beacon_flood(devices), deauth(frames), evil_twin(devices), weak_security(devices),
                           eapol(frames), randomization(links.is_array() ? links : Json::array()),
                           surveillance_and_identity(devices), inventory(devices)})
    for (const auto& f : part) out.push_back(f);
  auto conf = [](const Json& f) { return f.contains("confidence") ? num(f["confidence"]) : 0.0; };
  std::stable_sort(out.begin(), out.end(), [&](const Json& a, const Json& b) {
    const int sa = severity_order(a["severity"]), sb = severity_order(b["severity"]);
    if (sa != sb) return sa < sb;
    return -conf(a) < -conf(b);
  });
  return Json(out);
}

}  // namespace si::rx::devices
