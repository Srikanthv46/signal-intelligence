// python core/ble_trackers.py (see there for the layouts and the sources).
#include "core/ble_trackers.hpp"

#include <map>
#include <set>
#include <string>
#include <vector>

#include "core/pyfmt.hpp"
#include "core/pybytes.hpp"

namespace si::ble_trackers {

namespace {

using B = std::vector<uint8_t>;
using pybytes::Bytes;
namespace detail = si::pybytes;

const char* batt(int v) {
  static const char* t[4] = {"critical", "low", "medium", "full"};
  return t[v & 3];
}

const std::map<int, std::string>& apple_models() {
  static const std::map<int, std::string> m = {{0x0E20, "AirPods Pro"}, {0x0220, "AirPods"},
                                               {0x0F20, "AirPods Pro (2nd gen)"}, {0x0A20, "AirPods Max"},
                                               {0x0E30, "AirPods (3rd gen)"}, {0x1420, "AirPods Pro (2nd gen, USB-C)"}};
  return m;
}

std::optional<B> unhex(const std::string& h) {
  if (h.size() % 2) return std::nullopt;
  B b;
  for (size_t i = 0; i < h.size(); i += 2) {
    try {
      b.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    } catch (...) {
      return std::nullopt;
    }
  }
  return b;
}

std::string hex(const uint8_t* p, size_t n) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (size_t i = 0; i < n; ++i) s += d[p[i] >> 4], s += d[p[i] & 15];
  return s;
}

std::optional<B> svc(const Json& ad, int uuid) {
  if (!ad.contains("service_data")) return std::nullopt;
  for (const auto& s : ad["service_data"])
    if (std::stoi(s["uuid16"].get<std::string>(), nullptr, 16) == uuid) return unhex(s["data_hex"].get<std::string>());
  return std::nullopt;
}

Json obj(const char* vendor, const std::string& kind, const std::string& label) {
  Json o = Json::object();
  o["vendor"] = vendor;
  o["kind"] = kind;
  o["label"] = label;
  return o;
}

std::optional<Json> apple(const B& m) {
  std::vector<std::pair<int, B>> tlvs;
  size_t i = 0;
  while (i + 2 <= m.size()) {
    const int t = m[i], ln = m[i + 1];
    if (i + 2 + static_cast<size_t>(ln) > m.size()) break;
    tlvs.emplace_back(t, B(m.begin() + static_cast<long>(i) + 2, m.begin() + static_cast<long>(i) + 2 + ln));
    i += 2 + static_cast<size_t>(ln);
  }
  for (const auto& [t, v] : tlvs)
    if (t == 0x02 && v.size() == 0x15) {
      Json o = obj("Apple", "iBeacon", "iBeacon");
      o["proximity_uuid"] = hex(v.data(), 16);
      o["major"] = (v[16] << 8) | v[17];
      o["minor"] = (v[18] << 8) | v[19];
      o["measured_power_dbm"] = v[20] >= 128 ? v[20] - 256 : v[20];
      return o;
    }
  for (const auto& [t, v] : tlvs) {
    if (t == 0x12 && !v.empty()) {
      const int st = v[0];
      const int dt = (st & 0x30) >> 4;
      const std::string dev = dt == 0 ? "Apple device (Find My)" : dt == 1 ? "AirTag" : dt == 3 ? "AirPods" : "Apple device";
      Json o = obj("Apple", "Find My", "Apple Find My / " + dev);
      o["device_type"] = dev;
      o["offline_finding"] = v.size() >= 0x18;
      o["maintained"] = static_cast<bool>(st & 0x04);
      if (st & 0x04) o["battery"] = batt((st >> 6) & 3);
      return o;
    }
    if (t == 0x07 && v.size() >= 3) {
      const int model = (v[1] << 8) | v[2];
      const auto it = apple_models().find(model);
      Json o = obj("Apple", "proximity pairing", "Apple " + (it != apple_models().end() ? it->second : std::string("device")));
      o["model_id"] = fmt("0x%04X", model);
      o["model"] = it != apple_models().end() ? Json(it->second) : Json(nullptr);
      return o;
    }
  }
  if (!tlvs.empty()) {
    Json o = obj("Apple", "Continuity", "Apple device (Continuity)");
    o["continuity_types"] = Json::array();
    for (const auto& [t, _] : tlvs) o["continuity_types"].push_back(fmt("0x%02X", t));
    return o;
  }
  return std::nullopt;
}

std::optional<Json> recognize_raw(const Json& ad) {
  if (ad.contains("manufacturer_data_hex") && ad.contains("company_id")) {
    const auto m = unhex(ad["manufacturer_data_hex"].get<std::string>());
    if (m) {
      const int cid = std::stoi(ad["company_id"].get<std::string>(), nullptr, 16);
      if (cid == 0x004C && !m->empty())
        if (auto r = apple(*m)) return r;
      if (cid == 0x0006 && m->size() >= 3 && (*m)[2] == 0x80 && (*m)[1] <= 2) {
        static const char* sc[3] = {"LE", "BR/EDR", "LE + BR/EDR"};
        Json o = obj("Microsoft", "Swift Pair", "Microsoft Swift Pair");
        o["pairing"] = sc[(*m)[1]];
        std::string name = m->size() > 3 ? detail::decode_utf8(Bytes(m->data() + 3, m->size() - 3)) : std::string();
        o["display_name"] = name.empty() ? Json(nullptr) : Json(name);
        return o;
      }
    }
  }
  for (const auto& [uuid, reg] : std::vector<std::pair<int, bool>>{{0xFD5A, true}, {0xFD59, false}}) {
    const auto d = svc(ad, uuid);
    if (d && d->size() >= 20) {
      const int state = ((*d)[0] >> 5) & 7;
      if (1 <= state && state <= 6) {
        static const char* sn[7] = {"", "premature-lost", "lost", "overmature-lost", "paired", "connected", "connected (two)"};
        Json o = obj("Samsung", "SmartTag", "Samsung SmartTag");
        o["registered"] = reg;
        o["state"] = sn[state];
        o["battery"] = batt(((*d)[12] >> 6) & 3);
        o["uwb"] = static_cast<bool>((*d)[12] & 0x20);
        return o;
      }
    }
  }
  if (const auto d = svc(ad, 0xFEAA); d && !d->empty()) {
    const int fr = (*d)[0];
    const size_t n = d->size();
    if ((fr == 0x40 || fr == 0x41) && (n == 21 || n == 22 || n == 33 || n == 34)) {
      Json o = obj("Google", "Find My Device (FMDN)", "Google Find My (FMDN)");
      o["unwanted_tracking_protection"] = fr == 0x41;
      o["eid_bytes"] = n < 33 ? 20 : 32;
      return o;
    }
    static const std::map<int, std::string> names = {{0x00, "UID"}, {0x10, "URL"}, {0x20, "TLM"}, {0x30, "EID"}};
    if (const auto it = names.find(fr); it != names.end()) return obj("Eddystone", it->second, "Eddystone-" + it->second);
  }
  if (const auto d = svc(ad, 0xFE2C)) {
    if (d->size() == 3) {
      Json o = obj("Google", "Fast Pair", "Google Fast Pair");
      std::string h = hex(d->data(), 3);
      for (auto& c : h) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      o["model_id"] = h;
      return o;
    }
    if (!d->empty() && (*d)[0] == 0x00) return obj("Google", "Fast Pair", "Google Fast Pair (paired)");
  }
  for (const auto& [uuid, name] : std::vector<std::pair<int, std::string>>{
           {0xFEED, "Tile"}, {0xFEEC, "Tile"}, {0xFD84, "Tile"}, {0xFE33, "Chipolo"}, {0xFD69, "Samsung Find My Mobile"}})
    if (svc(ad, uuid)) {
      const std::string vendor = name.substr(0, name.find(' '));
      Json o = Json::object();
      o["vendor"] = vendor;
      o["kind"] = name;
      o["label"] = name + "-like";
      return o;
    }
  if (const auto d = svc(ad, 0xFCB2); d && d->size() >= 2 && ((*d)[1] & 0xFE) == 0) {
    const bool near = (*d)[1] & 0x01;
    Json o = obj("DULT", "tracker", near ? "Tracker (DULT, near-owner)" : "Tracker (DULT, separated)");
    o["near_owner"] = near;
    return o;
  }
  return std::nullopt;
}

}  // namespace

std::optional<Json> recognize(const Json& ad) {
  static const std::set<std::string> trackers = {"Find My", "SmartTag", "Find My Device (FMDN)", "Tile", "Chipolo",
                                                 "tracker"};
  auto r = recognize_raw(ad);
  if (r) (*r)["tracker"] = trackers.count((*r)["kind"].get<std::string>()) > 0;
  return r;
}

}  // namespace si::ble_trackers
