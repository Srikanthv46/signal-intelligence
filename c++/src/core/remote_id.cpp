// ASTM F3411 / ASD-STAN prEN 4709-002 Open Drone ID message parsing shared by every transport, and the Bluetooth LE
// advertising-data parser (python core/remote_id.py). Moved from receivers/wifi/droneid.cpp (2026-10-21).
#include "core/remote_id.hpp"

#include "core/ble_trackers.hpp"

#include <cmath>
#include <cstring>
#include <map>

#include "core/drone_models.hpp"
#include "core/pybytes.hpp"
#include "core/pyfmt.hpp"

namespace si::remote_id {

using pybytes::Bytes;
using pybytes::jget;
using pybytes::truthy;
namespace detail = si::pybytes;

namespace {

constexpr size_t kMsgLen = 25;

const char* lookup(const std::map<int, const char*>& m, int k) {
  auto it = m.find(k);
  return it == m.end() ? nullptr : it->second;
}
// TABLE.get(k, k)
Json get_or_int(const std::map<int, const char*>& m, int k) {
  const char* s = lookup(m, k);
  return s ? Json(s) : Json(k);
}

const std::map<int, const char*> kMsgType = {{0, "basic_id"}, {1, "location"}, {2, "authentication"},
                                             {3, "self_id"},  {4, "system"},   {5, "operator_id"},
                                             {0xF, "message_pack"}};
const std::map<int, const char*> kIdType = {{0, "none"},
                                            {1, "serial number (CTA-2063-A)"},
                                            {2, "CAA registration ID"},
                                            {3, "UTM assigned UUID"},
                                            {4, "specific session ID"}};
const std::map<int, const char*> kUaType = {
    {0, "undeclared"},  {1, "aeroplane"},     {2, "multirotor"},     {3, "gyroplane"},
    {4, "hybrid lift"}, {5, "ornithopter"},   {6, "glider"},         {7, "kite"},
    {8, "free balloon"}, {9, "captive balloon"}, {10, "airship"},    {11, "free fall / parachute"},
    {12, "rocket"},     {13, "tethered powered aircraft"}, {14, "ground obstacle"}, {15, "other"}};
const std::map<int, const char*> kOpStatus = {
    {0, "undeclared"}, {1, "on ground"}, {2, "airborne"}, {3, "emergency"}, {4, "remote ID system failure"}};
const std::map<int, const char*> kOpLocType = {
    {0, "takeoff point"}, {1, "live GNSS (operator moving)"}, {2, "fixed location"}};
const std::map<int, const char*> kDescType = {{0, "text"}, {1, "emergency"}, {2, "extended status"}};
const std::map<int, const char*> kHorizAcc = {{1, "<10 NM"},  {2, "<4 NM"},   {3, "<2 NM"},   {4, "<1 NM"},
                                              {5, "<0.5 NM"}, {6, "<0.3 NM"}, {7, "<0.1 NM"}, {8, "<0.05 NM"},
                                              {9, "<30 m"},   {10, "<10 m"},  {11, "<3 m"},   {12, "<1 m"}};
const std::map<int, const char*> kVertAcc = {{1, "<150 m"}, {2, "<45 m"}, {3, "<25 m"},
                                             {4, "<10 m"},  {5, "<3 m"},  {6, "<1 m"}};
const std::map<int, const char*> kSpeedAcc = {{1, "<10 m/s"}, {2, "<3 m/s"}, {3, "<1 m/s"}, {4, "<0.3 m/s"}};
const std::map<int, const char*> kUaClassEu = {{0, "undeclared"}, {1, "open"}, {2, "specific"}, {3, "certified"}};
int32_t s32(const Bytes& b, size_t i) {
  uint32_t v = static_cast<uint32_t>(b[i]) | (static_cast<uint32_t>(b[i + 1]) << 8) |
               (static_cast<uint32_t>(b[i + 2]) << 16) | (static_cast<uint32_t>(b[i + 3]) << 24);
  return static_cast<int32_t>(v);
}
uint32_t u32(const Bytes& b, size_t i) { return static_cast<uint32_t>(s32(b, i)); }
int u16(const Bytes& b, size_t i) { return detail::u16le(b, i); }
int16_t s16(const Bytes& b, size_t i) { return static_cast<int16_t>(u16(b, i)); }

// b.split(b"\x00", 1)[0].decode("utf-8", "replace").strip()
std::string text(const Bytes& b) {
  size_t n = 0;
  while (n < b.size() && b.data()[n] != 0) ++n;
  return detail::py_strip(detail::decode_utf8(b.sub(0, n)));
}

// metres (0.5 m steps, -1000 m offset); 0 = not reported
std::optional<double> alt(int raw) {
  if (raw == 0) return std::nullopt;
  return py_round(raw * 0.5 - 1000, 1);
}
std::optional<double> latlon(int32_t v) {
  if (v == 0) return std::nullopt;
  return py_round(v * 1e-7, 7);
}

Json basic_id(const Bytes& b) {
  Json out = Json::object();
  out["id_type"] = get_or_int(kIdType, b[0] >> 4);
  out["ua_type"] = get_or_int(kUaType, b[0] & 0xF);
  std::string uas = text(b.sub(1, 21));
  if (!uas.empty()) {
    out["uas_id"] = uas;
    const Json mf = (b[0] >> 4) == 1 ? drone_models::serial_manufacturer(uas) : Json(nullptr);
    if (mf.is_object()) out["manufacturer"] = mf;   // CTA-2063-A manufacturer code (research, PLAN M2)
  }
  return out;
}

Json location(const Bytes& b) {
  const int st = b[0];
  const int mult = st & 1, ew = (st >> 1) & 1;
  Json out = Json::object();
  out["status"] = get_or_int(kOpStatus, st >> 4);
  out["height_type"] = ((st >> 2) & 1) ? "above ground" : "above takeoff";
  if (b[1] != 181) out["track_deg"] = b[1] + 180 * ew;
  if (b[2] != 255) out["speed_mps"] = py_round(mult ? b[2] * 0.75 + 63.75 : b[2] * 0.25, 2);
  const int vs = static_cast<int8_t>(b[3]);
  if (vs != 63) out["vertical_speed_mps"] = py_round(vs * 0.5, 2);
  auto lat = latlon(s32(b, 4)), lon = latlon(s32(b, 8));
  if (lat && lon && std::fabs(*lat) <= 90 && std::fabs(*lon) <= 180) {
    out["lat"] = *lat;
    out["lon"] = *lon;
  }
  const std::pair<const char*, size_t> alts[] = {
      {"pressure_altitude_m", 12}, {"geodetic_altitude_m", 14}, {"height_m", 16}};
  for (const auto& [k, i] : alts) {
    auto v = alt(u16(b, i));
    if (v) out[k] = *v;
  }
  Json acc = Json::object();
  const std::pair<const char*, const char*> a[] = {{"horizontal", lookup(kHorizAcc, b[18] & 0xF)},
                                                   {"vertical", lookup(kVertAcc, b[18] >> 4)},
                                                   {"speed", lookup(kSpeedAcc, b[19] & 0xF)},
                                                   {"barometric", lookup(kVertAcc, b[19] >> 4)}};
  for (const auto& [k, v] : a)
    if (v) acc[k] = v;
  if (!acc.empty()) out["accuracy"] = acc;
  const int ts = u16(b, 20);
  if (ts != 0xFFFF) out["timestamp_s_into_hour"] = py_round(ts * 0.1, 1);
  return out;
}

Json auth(const Bytes& b) {
  Json out = Json::object();
  out["auth_type"] = b[0] >> 4;
  out["page"] = b[0] & 0xF;
  if ((b[0] & 0xF) == 0 && b.size() >= 7) {
    out["page_count"] = b[1];
    out["length"] = b[2];
    out["timestamp"] = u32(b, 3);
  }
  return out;
}

Json self_id(const Bytes& b) {
  Json out = Json::object();
  out["description_type"] = get_or_int(kDescType, b[0]);
  std::string d = text(b.sub(1, 24));
  if (!d.empty()) out["description"] = d;
  return out;
}

Json system_msg(const Bytes& b) {
  Json out = Json::object();
  out["operator_location_type"] = get_or_int(kOpLocType, b[0] & 3);
  const int cls = (b[0] >> 2) & 7;
  out["classification_type"] = cls == 1 ? Json("EU") : Json(cls);
  auto lat = latlon(s32(b, 1)), lon = latlon(s32(b, 5));
  if (lat && lon && std::fabs(*lat) <= 90 && std::fabs(*lon) <= 180) {
    out["operator_lat"] = *lat;
    out["operator_lon"] = *lon;
  }
  const int n = u16(b, 9);
  if (n > 1) {
    out["area_count"] = n;
    out["area_radius_m"] = b[11] * 10;
  }
  if (cls == 1) {
    out["ua_class_eu"] = get_or_int(kUaClassEu, b[16] >> 4);
    out["ua_category_eu"] = b[16] & 0xF;
  }
  auto v = alt(u16(b, 17));
  if (v) out["operator_altitude_m"] = *v;
  const uint32_t ts = u32(b, 19);
  if (ts) out["timestamp_s_since_2019"] = ts;
  return out;
}

Json operator_id(const Bytes& b) {
  Json out = Json::object();
  out["operator_id_type"] = b[0] == 0 ? Json("CAA-issued") : Json(b[0]);
  std::string v = text(b.sub(1, 21));
  if (!v.empty()) out["operator_id"] = v;
  return out;
}

// one 25-byte message -> {fields..., message_type, protocol_version}
std::optional<Json> message(const Bytes& m) {
  if (m.size() < kMsgLen) return std::nullopt;
  const int t = m[0] >> 4, ver = m[0] & 0xF;
  if (t > 5 || ver > 3) return std::nullopt;
  Json out;
  try {
    const Bytes b = m.sub(1, kMsgLen);
    switch (t) {
      case 0: out = basic_id(b); break;
      case 1: out = location(b); break;
      case 2: out = auth(b); break;
      case 3: out = self_id(b); break;
      case 4: out = system_msg(b); break;
      default: out = operator_id(b); break;
    }
  } catch (const detail::IndexError&) {
    return std::nullopt;
  }
  out["message_type"] = lookup(kMsgType, t);
  out["protocol_version"] = ver;
  return out;
}

bool informative(const std::vector<Json>& msgs) {
  static const char* kContent[] = {"uas_id", "lat", "operator_lat", "operator_id", "description",
                                   "speed_mps", "height_m", "geodetic_altitude_m", "track_deg"};
  for (const auto& m : msgs) {
    for (const char* k : kContent)
      if (!jget(m, k).is_null()) return true;
    const Json& st = jget(m, "status");
    const Json& ua = jget(m, "ua_type");
    if (!(st.is_null() || st == "undeclared") || !(ua.is_null() || ua == "undeclared")) return true;
  }
  return false;
}
}  // namespace

std::optional<Json> parse_gb46750(const std::vector<uint8_t>& bv, const std::string& transport) {
  const Bytes b(bv);
  if (b.size() < 6 || b[0] != 0xFF || (b[1] >> 5) != 1) return std::nullopt;
  const int ln = b[2];
  std::vector<int> flags;
  size_t i = 3;
  while (i < b.size()) {
    flags.push_back(b[i]);
    ++i;
    if (!(b[i - 1] & 1)) break;
  }
  struct Item { int id, n; const char* key; };
  static const Item kItems[] = {
      {1, 20, "uas_id"}, {2, 8, "registration_mark"}, {3, 1, "operation_category"}, {4, 1, "ua_class"},
      {5, 1, "operator_location_type"}, {6, 8, "operator_position"}, {7, 2, "operator_altitude_m"},
      {8, 8, "position"}, {9, 2, "track_deg"}, {10, 2, "speed_mps"}, {11, 2, "height_m"}, {12, 1, "vertical_speed_mps"},
      {13, 2, "geodetic_altitude_m"}, {14, 2, "pressure_altitude_m"}, {15, 1, "status"}, {16, 1, "coordinate_system"},
      {17, 1, "horizontal_accuracy"}, {18, 1, "vertical_accuracy"}, {19, 1, "speed_accuracy"}, {20, 6, "timestamp_ms"},
      {21, 1, "timestamp_accuracy"}};
  auto present = [&](int k) {
    if (k - 1 >= 7 * static_cast<int>(flags.size())) return false;
    return (flags[static_cast<size_t>((k - 1) / 7)] & (0x80 >> ((k - 1) % 7))) != 0;
  };
  int total = 0;
  for (const auto& it : kItems)
    if (present(it.id)) total += it.n;
  if (total != ln || b.size() < i + static_cast<size_t>(ln) || !present(1)) return std::nullopt;
  static const char* kCat[] = {"undefined", "open", "specific", "certified"};
  static const char* kClass[] = {"micro", "light", "small", "medium", "large"};
  static const char* kStatus[] = {"not reported", "ground", "airborne", "emergency", "RID failure", "RID failure + emergency"};
  Json out = Json::object();
  out["standard"] = "China GB 46750-2025";
  out["transport"] = transport;
  out["version"] = std::to_string(b[1] >> 5) + "." + std::to_string(b[1] & 0x1F);
  out["messages"] = Json::array({"GB 46750"});
  for (const auto& it : kItems) {
    if (!present(it.id)) continue;
    const Bytes v = b.sub(i, i + static_cast<size_t>(it.n));
    i += static_cast<size_t>(it.n);
    uint64_t u = 0;
    for (int k = it.n - 1; k >= 0; --k) u = (u << 8) | v[static_cast<size_t>(k)];
    const std::string key = it.key;
    if (key == "uas_id" || key == "registration_mark") {
      const std::string t = text(v);
      if (!t.empty() && t.find_first_not_of('0') != std::string::npos) out[key] = t;
    } else if (key == "operation_category") {
      out[key] = u < 4 ? Json(kCat[u]) : Json(u);
    } else if (key == "ua_class") {
      out[key] = u < 5 ? Json(kClass[u]) : Json(u);
    } else if (key == "operator_location_type") {
      out[key] = u == 0 ? "take-off point" : "remote-station position";
    } else if (key == "operator_position" || key == "position") {
      const int32_t lo = s32(v, 0), la = s32(v, 4);
      if (lo == -1 || la == -1) continue;
      const double lon = lo * 1e-7, lat = la * 1e-7;
      if (std::fabs(lat) > 90 || std::fabs(lon) > 180) continue;
      if (key == "position") {
        out["lat"] = py_round(lat, 7);
        out["lon"] = py_round(lon, 7);
      } else {
        out["operator_lat"] = py_round(lat, 7);
        out["operator_lon"] = py_round(lon, 7);
      }
    } else if (key == "operator_altitude_m" || key == "geodetic_altitude_m" || key == "pressure_altitude_m") {
      if (u) out[key] = py_round(static_cast<double>(u) / 2 - 1000, 1);
    } else if (key == "height_m") {
      if (u) out[key] = py_round(static_cast<double>(u) / 2 - 9000, 1);
    } else if (key == "track_deg" || key == "speed_mps") {
      if (u != 0xFFFF) out[key] = static_cast<double>(u) / 10;
    } else if (key == "vertical_speed_mps") {
      if (u != 0xFF) out[key] = py_round((u & 0x7F) * 0.5 * ((u & 0x80) ? -1 : 1), 1);
    } else if (key == "status") {
      out[key] = u < 6 ? Json(kStatus[u]) : Json(u);
    } else if (key == "coordinate_system") {
      out[key] = u == 1 ? "CGCS2000" : "WGS-84";
    } else {
      out[key] = u;
    }
  }
  if (out.contains("uas_id")) {
    const Json mf = drone_models::serial_manufacturer(out["uas_id"].get<std::string>());
    if (!mf.is_null()) out["manufacturer"] = mf;
  }
  return out;
}

std::optional<Json> parse_odid(const std::vector<uint8_t>& bodyv, const std::string& transport) {
  for (size_t start : {size_t(1), size_t(0)}) {   // China GB 46750-2025 in the same element (0xFF first)
    if (bodyv.size() <= start) continue;
    auto g = parse_gb46750(std::vector<uint8_t>(bodyv.begin() + static_cast<std::ptrdiff_t>(start), bodyv.end()), transport);
    if (g) {
      if (start) (*g)["message_counter"] = bodyv[0];
      return g;
    }
  }
  const Bytes body(bodyv);
  std::vector<Json> msgs;
  std::optional<int> counter;
  for (int start : {1, 0}) {
    const Bytes b = body.from(start);
    if (b.size() < kMsgLen) continue;
    if ((b[0] >> 4) == 0xF && b.size() >= 3 && b[1] == static_cast<int>(kMsgLen)) {  // message pack
      const int n = b[2];
      std::vector<Json> got;
      for (int k = 0; k < n; ++k) {
        auto g = message(b.sub(3 + kMsgLen * k, 3 + kMsgLen * (k + 1)));
        if (g) got.push_back(std::move(*g));
      }
      if (!got.empty() && informative(got)) {
        msgs = std::move(got);
        if (start) counter = body[0];
        break;
      }
    }
    std::vector<Json> got;
    for (size_t k = 0; k < b.size() / kMsgLen; ++k) {
      auto g = message(b.sub(kMsgLen * k, kMsgLen * (k + 1)));
      if (g) got.push_back(std::move(*g));
    }
    if (!got.empty() && informative(got)) {
      msgs = std::move(got);
      if (start) counter = body[0];
      break;
    }
  }
  if (msgs.empty()) return std::nullopt;
  Json out = Json::object();
  out["standard"] = "ASTM F3411 / ASD-STAN prEN 4709-002 (Open Drone ID)";
  out["transport"] = transport;
  out["messages"] = Json::array();
  if (counter) out["message_counter"] = *counter;
  for (const auto& m : msgs) {
    out["messages"].push_back(m["message_type"]);
    for (auto it = m.begin(); it != m.end(); ++it) {
      if (it.key() == "message_type" || it.key() == "protocol_version") continue;
      if (!out.contains(it.key())) out[it.key()] = it.value();
    }
    out["protocol_version"] = m["protocol_version"];
  }
  return out;
}

namespace {
const std::map<int, const char*> kCompanies = {{0x004C, "Apple"}, {0x0006, "Microsoft"}, {0x0075, "Samsung"},
                                               {0x00E0, "Google"}, {0x0059, "Nordic Semiconductor"},
                                               {0x02E5, "Espressif"}, {0x000D, "Texas Instruments"},
                                               {0x0131, "Cypress"}, {0x0157, "Huami"}, {0x038F, "Xiaomi"}, {0x0822, "DJI"}};
std::string hex_bytes(const uint8_t* p, size_t n) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (size_t i = 0; i < n; ++i) s += d[p[i] >> 4], s += d[p[i] & 15];
  return s;
}
}  // namespace

Json parse_ble_ad(const std::vector<uint8_t>& b) {
  Json out = Json::object();
  size_t i = 0;
  while (i < b.size()) {
    const size_t ln = b[i];
    if (ln == 0 || i + 1 + ln > b.size()) break;
    const int t = b[i + 1];
    const uint8_t* v = b.data() + i + 2;
    const size_t n = ln - 1;
    i += 1 + ln;
    if (t == 0x01 && n) {
      out["flags"] = v[0];
    } else if ((t == 0x02 || t == 0x03) && n >= 2) {
      if (!out.contains("service_uuids16")) out["service_uuids16"] = Json::array();
      for (size_t k = 0; k + 1 < n; k += 2) out["service_uuids16"].push_back(fmt("%04X", v[k] | (v[k + 1] << 8)));
    } else if (t == 0x08 || t == 0x09) {
      out["local_name"] = detail::decode_utf8(Bytes(v, n));
    } else if (t == 0x0A && n) {
      out["tx_power_dbm"] = static_cast<int8_t>(v[0]);
    } else if (t == 0x19 && n >= 2) {
      out["appearance"] = v[0] | (v[1] << 8);
    } else if (t == 0xFF && n >= 2) {
      const int cid = v[0] | (v[1] << 8);
      out["company_id"] = fmt("0x%04X", cid);
      auto it = kCompanies.find(cid);
      if (it != kCompanies.end()) out["company"] = it->second;
      out["manufacturer_data_hex"] = hex_bytes(v + 2, n - 2);
    } else if (t == 0x16 && n >= 2) {
      const int uuid = v[0] | (v[1] << 8);
      if (!out.contains("service_data")) out["service_data"] = Json::array();
      out["service_data"].push_back({{"uuid16", fmt("%04X", uuid)}, {"data_hex", hex_bytes(v + 2, n - 2)}});
      if (uuid == 0xFFFA && n >= 3 && v[2] == 0x0D) {
        auto rid = parse_odid(std::vector<uint8_t>(v + 3, v + n), "Bluetooth LE advertising (legacy)");
        if (rid) out["remote_id"] = *rid;
      }
    }
  }
  // item trackers / notable devices (research, signal gaps S7)
  if (auto dev = ble_trackers::recognize(out)) out["device"] = *dev;
  return out;
}

}  // namespace si::remote_id
