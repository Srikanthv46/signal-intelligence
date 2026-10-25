// python core/mavlink.py (same table, order, rounding).
#include "core/mavlink.hpp"

#include <cstring>
#include <map>
#include <string>

#include "core/pybytes.hpp"
#include "core/pyfmt.hpp"

namespace si::mavlink {

namespace {

enum class T { U8, I8, U16, I16, U32, I32, U64, I64, F32, STR16, STR50 };
struct Field { const char* name; T t; size_t off; };
struct Msg { const char* name; int extra; size_t full; std::vector<Field> f; };

size_t tsize(T t) {
  switch (t) {
    case T::U8: case T::I8: return 1;
    case T::U16: case T::I16: return 2;
    case T::U32: case T::I32: case T::F32: return 4;
    case T::U64: case T::I64: return 8;
    case T::STR16: return 16;
    case T::STR50: return 50;
  }
  return 0;
}

const std::map<int, Msg>& msgs() {
  static const std::map<int, Msg> m = {
      {0, {"HEARTBEAT", 50, 9, {{"custom_mode", T::U32, 0}, {"type", T::U8, 4}, {"autopilot", T::U8, 5}, {"base_mode", T::U8, 6},
                                 {"system_status", T::U8, 7}, {"mavlink_version", T::U8, 8}}}},
      {1, {"SYS_STATUS", 124, 43, {{"load", T::U16, 12}, {"voltage_battery", T::U16, 14}, {"current_battery", T::I16, 16},
                                   {"drop_rate_comm", T::U16, 18}, {"errors_comm", T::U16, 20}, {"battery_remaining", T::I8, 30}}}},
      {2, {"SYSTEM_TIME", 137, 12, {{"time_unix_usec", T::U64, 0}, {"time_boot_ms", T::U32, 8}}}},
      {22, {"PARAM_VALUE", 220, 25, {{"param_value", T::F32, 0}, {"param_count", T::U16, 4}, {"param_index", T::U16, 6},
                                     {"param_id", T::STR16, 8}, {"param_type", T::U8, 24}}}},
      {24, {"GPS_RAW_INT", 24, 52, {{"time_usec", T::U64, 0}, {"lat", T::I32, 8}, {"lon", T::I32, 12}, {"alt", T::I32, 16},
                                    {"eph", T::U16, 20}, {"epv", T::U16, 22}, {"vel", T::U16, 24}, {"cog", T::U16, 26},
                                    {"fix_type", T::U8, 28}, {"satellites_visible", T::U8, 29}}}},
      {27, {"RAW_IMU", 144, 29, {{"time_usec", T::U64, 0}}}},
      {29, {"SCALED_PRESSURE", 115, 16, {{"time_boot_ms", T::U32, 0}, {"press_abs", T::F32, 4}, {"temperature", T::I16, 12}}}},
      {30, {"ATTITUDE", 39, 28, {{"time_boot_ms", T::U32, 0}, {"roll", T::F32, 4}, {"pitch", T::F32, 8}, {"yaw", T::F32, 12}}}},
      {33, {"GLOBAL_POSITION_INT", 104, 28, {{"time_boot_ms", T::U32, 0}, {"lat", T::I32, 4}, {"lon", T::I32, 8}, {"alt", T::I32, 12},
                                             {"relative_alt", T::I32, 16}, {"vx", T::I16, 20}, {"vy", T::I16, 22}, {"vz", T::I16, 24},
                                             {"hdg", T::U16, 26}}}},
      {35, {"RC_CHANNELS_RAW", 244, 22, {{"time_boot_ms", T::U32, 0}, {"rssi", T::U8, 21}}}},
      {36, {"SERVO_OUTPUT_RAW", 222, 37, {{"time_usec", T::U32, 0}}}},
      {42, {"MISSION_CURRENT", 28, 18, {{"seq", T::U16, 0}}}},
      {62, {"NAV_CONTROLLER_OUTPUT", 183, 26, {{"nav_bearing", T::I16, 20}, {"target_bearing", T::I16, 22}, {"wp_dist", T::U16, 24}}}},
      {65, {"RC_CHANNELS", 118, 42, {{"time_boot_ms", T::U32, 0}, {"chancount", T::U8, 40}, {"rssi", T::U8, 41}}}},
      {74, {"VFR_HUD", 20, 20, {{"airspeed", T::F32, 0}, {"groundspeed", T::F32, 4}, {"alt", T::F32, 8}, {"climb", T::F32, 12},
                                {"heading", T::I16, 16}, {"throttle", T::U16, 18}}}},
      {76, {"COMMAND_LONG", 152, 33, {{"command", T::U16, 28}, {"target_system", T::U8, 30}, {"target_component", T::U8, 31}}}},
      {77, {"COMMAND_ACK", 143, 10, {{"command", T::U16, 0}, {"result", T::U8, 2}}}},
      {109, {"RADIO_STATUS", 185, 9, {{"rxerrors", T::U16, 0}, {"fixed", T::U16, 2}, {"rssi", T::U8, 4}, {"remrssi", T::U8, 5},
                                      {"txbuf", T::U8, 6}, {"noise", T::U8, 7}, {"remnoise", T::U8, 8}}}},
      {111, {"TIMESYNC", 34, 18, {{"tc1", T::I64, 0}, {"ts1", T::I64, 8}}}},
      {147, {"BATTERY_STATUS", 154, 54, {{"current_consumed", T::I32, 0}, {"energy_consumed", T::I32, 4}, {"temperature", T::I16, 8},
                                         {"voltage_cell1", T::U16, 10}, {"current_battery", T::I16, 30}, {"battery_remaining", T::I8, 35}}}},
      {148, {"AUTOPILOT_VERSION", 178, 78, {{"capabilities", T::U64, 0}, {"uid", T::U64, 8}, {"flight_sw_version", T::U32, 16},
                                            {"board_version", T::U32, 28}, {"vendor_id", T::U16, 32}, {"product_id", T::U16, 34}}}},
      {166, {"RADIO", 21, 9, {{"rxerrors", T::U16, 0}, {"fixed", T::U16, 2}, {"rssi", T::U8, 4}, {"remrssi", T::U8, 5},
                              {"txbuf", T::U8, 6}, {"noise", T::U8, 7}, {"remnoise", T::U8, 8}}}},
      {242, {"HOME_POSITION", 104, 60, {{"latitude", T::I32, 0}, {"longitude", T::I32, 4}, {"altitude", T::I32, 8}}}},
      {245, {"EXTENDED_SYS_STATE", 130, 2, {{"vtol_state", T::U8, 0}, {"landed_state", T::U8, 1}}}},
      {253, {"STATUSTEXT", 83, 54, {{"severity", T::U8, 0}, {"text", T::STR50, 1}}}},
  };
  return m;
}

const std::map<int, const char*> kMavType = {
    {0, "generic"}, {1, "fixed wing"}, {2, "quadrotor"}, {3, "coaxial"}, {4, "helicopter"}, {5, "antenna tracker"}, {6, "GCS"},
    {7, "airship"}, {8, "free balloon"}, {9, "rocket"}, {10, "ground rover"}, {11, "surface boat"}, {12, "submarine"},
    {13, "hexarotor"}, {14, "octorotor"}, {15, "tricopter"}, {16, "flapping wing"}, {17, "kite"}, {18, "onboard controller"},
    {19, "VTOL tailsitter duorotor"}, {20, "VTOL tailsitter quadrotor"}, {21, "VTOL tiltrotor"}, {22, "VTOL fixedrotor"},
    {23, "VTOL tailsitter"}, {24, "VTOL tiltwing"}, {26, "gimbal"}, {27, "ADS-B"}, {28, "parafoil"}, {29, "dodecarotor"},
    {30, "camera"}, {32, "FLARM"}, {34, "ODID"}, {35, "decarotor"}};
const std::map<int, const char*> kAutopilot = {{0, "generic"}, {3, "ArduPilot"}, {4, "OpenPilot"}, {8, "invalid (GCS / peripheral)"},
                                               {9, "Paparazzi"}, {12, "PX4"}, {13, "SMACCMPilot"}, {14, "AutoQuad"},
                                               {18, "SmartAP"}, {20, "Reflex"}, {21, "Flix"}};
const std::map<int, const char*> kState = {{0, "uninit"}, {1, "boot"}, {2, "calibrating"}, {3, "standby"}, {4, "active"},
                                           {5, "critical"}, {6, "emergency"}, {7, "poweroff"}, {8, "flight termination"}};
const std::map<int, const char*> kCopter = {
    {0, "STABILIZE"}, {1, "ACRO"}, {2, "ALT_HOLD"}, {3, "AUTO"}, {4, "GUIDED"}, {5, "LOITER"}, {6, "RTL"}, {7, "CIRCLE"}, {9, "LAND"},
    {11, "DRIFT"}, {13, "SPORT"}, {14, "FLIP"}, {15, "AUTOTUNE"}, {16, "POSHOLD"}, {17, "BRAKE"}, {18, "THROW"}, {19, "AVOID_ADSB"},
    {20, "GUIDED_NOGPS"}, {21, "SMART_RTL"}, {22, "FLOWHOLD"}, {23, "FOLLOW"}, {24, "ZIGZAG"}, {25, "SYSTEMID"},
    {26, "AUTOROTATE"}, {27, "AUTO_RTL"}, {28, "TURTLE"}};
const std::map<int, const char*> kPlane = {
    {0, "MANUAL"}, {1, "CIRCLE"}, {2, "STABILIZE"}, {3, "TRAINING"}, {4, "ACRO"}, {5, "FBWA"}, {6, "FBWB"}, {7, "CRUISE"},
    {8, "AUTOTUNE"}, {10, "AUTO"}, {11, "RTL"}, {12, "LOITER"}, {13, "TAKEOFF"}, {14, "AVOID_ADSB"}, {15, "GUIDED"},
    {16, "INITIALIZING"}, {17, "QSTABILIZE"}, {18, "QHOVER"}, {19, "QLOITER"}, {20, "QLAND"}, {21, "QRTL"}, {22, "QAUTOTUNE"},
    {23, "QACRO"}, {24, "THERMAL"}, {25, "LOITER_ALT_QLAND"}, {26, "AUTOLAND"}};
const std::map<int, const char*> kRover = {{0, "MANUAL"}, {1, "ACRO"}, {3, "STEERING"}, {4, "HOLD"}, {5, "LOITER"}, {6, "FOLLOW"},
                                           {7, "SIMPLE"}, {8, "DOCK"}, {9, "CIRCLE"}, {10, "AUTO"}, {11, "RTL"}, {12, "SMART_RTL"},
                                           {15, "GUIDED"}, {16, "INITIALIZING"}};
const std::map<int, const char*> kPx4Main = {{1, "MANUAL"}, {2, "ALTCTL"}, {3, "POSCTL"}, {4, "AUTO"}, {5, "ACRO"}, {6, "OFFBOARD"},
                                             {7, "STABILIZED"}, {8, "RATTITUDE_LEGACY"}, {9, "SIMPLE"}, {10, "TERMINATION"},
                                             {11, "ALTITUDE_CRUISE"}, {12, "MANUAL_PARKING"}};
const std::map<int, const char*> kPx4Auto = {{1, "READY"}, {2, "TAKEOFF"}, {3, "LOITER"}, {4, "MISSION"}, {5, "RTL"}, {6, "LAND"},
                                             {8, "FOLLOW_TARGET"}, {9, "PRECLAND"}, {10, "VTOL_TAKEOFF"}, {19, "GUIDED_COURSE"},
                                             {20, "AUTO_DESCEND"}};

template <class V>
V rd(const std::vector<uint8_t>& p, size_t o) {
  V v;
  std::memcpy(&v, p.data() + o, sizeof(V));  // little-endian host
  return v;
}

Json field_value(const std::vector<uint8_t>& p, const Field& f) {
  switch (f.t) {
    case T::U8: return p[f.off];
    case T::I8: return static_cast<int8_t>(p[f.off]);
    case T::U16: return rd<uint16_t>(p, f.off);
    case T::I16: return rd<int16_t>(p, f.off);
    case T::U32: return rd<uint32_t>(p, f.off);
    case T::I32: return rd<int32_t>(p, f.off);
    case T::U64: return rd<uint64_t>(p, f.off);
    case T::I64: return rd<int64_t>(p, f.off);
    case T::F32: return py_round(static_cast<double>(rd<float>(p, f.off)), 6);
    case T::STR16:
    case T::STR50: {
      size_t n = tsize(f.t), k = 0;
      while (k < n && p[f.off + k] != 0) ++k;
      return pybytes::decode_utf8(pybytes::Bytes(p.data() + f.off, k));
    }
  }
  return nullptr;
}

std::string lookup(const std::map<int, const char*>& m, int k, const std::string& dflt) {
  auto it = m.find(k);
  return it == m.end() ? dflt : std::string(it->second);
}

Json flight_mode(int autopilot, int type, uint32_t cm) {
  if (autopilot == 3) {
    const bool plane = type == 1 || type == 16 || (type >= 19 && type <= 24);
    const auto& tab = plane ? kPlane : (type == 10 || type == 11) ? kRover : kCopter;
    return lookup(tab, static_cast<int>(cm), "mode " + std::to_string(cm));
  }
  if (autopilot == 12) {
    const int main = (cm >> 16) & 0xFF, sub = (cm >> 24) & 0xFF;
    const std::string name = lookup(kPx4Main, main, "main " + std::to_string(main));
    if (main == 4) {
      auto it = kPx4Auto.find(sub);
      return name + "/" + (it == kPx4Auto.end() ? std::to_string(sub) : std::string(it->second));
    }
    return name;
  }
  return nullptr;
}

Json decode(int msgid, std::vector<uint8_t> p) {
  const Msg& m = msgs().at(msgid);
  if (p.size() < m.full) p.resize(m.full, 0);  // v2 zero-trimmed payloads
  Json out = Json::object();
  out["msg"] = m.name;
  for (const auto& f : m.f)
    if (f.off + tsize(f.t) <= p.size()) out[f.name] = field_value(p, f);
  if (std::string(m.name) == "HEARTBEAT") {
    const int type = out["type"].get<int>(), ap = out["autopilot"].get<int>(), bm = out["base_mode"].get<int>();
    out["vehicle_type"] = lookup(kMavType, type, "type " + std::to_string(type));
    out["autopilot_name"] = lookup(kAutopilot, ap, "autopilot " + std::to_string(ap));
    out["armed"] = static_cast<bool>(bm & 0x80);
    {
      auto it = kState.find(out["system_status"].get<int>());
      out["state"] = it == kState.end() ? Json(out["system_status"]) : Json(it->second);
    }
    if (bm & 0x01) {
      const Json fm = flight_mode(ap, type, out["custom_mode"].get<uint32_t>());
      if (!fm.is_null() && !(fm.is_string() && fm.get<std::string>().empty())) out["flight_mode"] = fm;
    }
  }
  for (const char* k : {"lat", "lon", "latitude", "longitude"})
    if (out.contains(k)) out[std::string(k) + "_deg"] = py_round(out[k].get<double>() * 1e-7, 7);
  return out;
}

}  // namespace

uint16_t x25(const uint8_t* p, size_t n, uint16_t crc) {
  for (size_t i = 0; i < n; ++i) {
    uint8_t t = static_cast<uint8_t>(p[i] ^ (crc & 0xFF));
    t = static_cast<uint8_t>(t ^ (t << 4));
    crc = static_cast<uint16_t>((crc >> 8) ^ (t << 8) ^ (t << 3) ^ (t >> 4));
  }
  return crc;
}

Json scan(const std::vector<uint8_t>& b) {
  Json out = Json::array();
  size_t i = 0;
  while (i < b.size()) {
    const int stx = b[i];
    if ((stx != 0xFE && stx != 0xFD) || i + 2 > b.size()) {
      ++i;
      continue;
    }
    const size_t ln = b[i + 1];
    size_t hdr, end, sig = 0;
    int seq, sysid, compid, msgid;
    if (stx == 0xFE) {
      hdr = 6;
      end = i + 6 + ln + 2;
      if (end > b.size()) { ++i; continue; }
      seq = b[i + 2]; sysid = b[i + 3]; compid = b[i + 4]; msgid = b[i + 5];
    } else {
      hdr = 10;
      if (i + hdr > b.size()) { ++i; continue; }
      const int incompat = b[i + 2];
      if (incompat & ~0x01) { ++i; continue; }
      sig = (incompat & 0x01) ? 13 : 0;
      end = i + hdr + ln + 2;
      if (end + sig > b.size()) { ++i; continue; }
      seq = b[i + 4]; sysid = b[i + 5]; compid = b[i + 6];
      msgid = b[i + 7] | (b[i + 8] << 8) | (b[i + 9] << 16);
    }
    const std::vector<uint8_t> payload(b.begin() + static_cast<ptrdiff_t>(i + hdr), b.begin() + static_cast<ptrdiff_t>(i + hdr + ln));
    const int crc_rx = b[end - 2] | (b[end - 1] << 8);
    auto it = msgs().find(msgid);
    Json ok = nullptr;
    if (it != msgs().end()) {
      const uint8_t ex = static_cast<uint8_t>(it->second.extra);
      ok = x25(&ex, 1, x25(b.data() + i + 1, hdr - 1 + ln)) == crc_rx;
    }
    Json fr = {{"version", stx == 0xFE ? 1 : 2}, {"seq", seq}, {"sysid", sysid}, {"compid", compid},
               {"msgid", msgid}, {"len", ln}, {"crc_ok", ok}};
    if (ok.is_boolean() && ok.get<bool>()) {
      const Json d = decode(msgid, payload);
      for (auto jt = d.begin(); jt != d.end(); ++jt) fr[jt.key()] = jt.value();
      out.push_back(fr);
      i = end + sig;
      continue;
    }
    if (ok.is_null() && stx == 0xFD) {
      out.push_back(fr);
      i = end + sig;
      continue;
    }
    ++i;
  }
  return out;
}

Json summary(const Json& frames) {
  Json s = Json::object();
  for (const auto& f : frames) {
    if (!(f.contains("crc_ok") && f["crc_ok"].is_boolean() && f["crc_ok"].get<bool>())) continue;
    if (!s.contains("sysid")) s["sysid"] = f["sysid"];
    if (!s.contains("messages")) s["messages"] = Json::array();
    const std::string m = f["msg"].get<std::string>();
    bool have = false;
    for (const auto& x : s["messages"]) have = have || x == m;
    if (!have) s["messages"].push_back(m);
    if (m == "HEARTBEAT" && f["sysid"].get<int>() != 255 && f["type"].get<int>() != 6) {
      s["vehicle_type"] = f["vehicle_type"];
      s["autopilot"] = f["autopilot_name"];
      s["armed"] = f["armed"];
      s["state"] = f["state"];
      if (f.contains("flight_mode")) s["flight_mode"] = f["flight_mode"];
    } else if (m == "GLOBAL_POSITION_INT") {
      s["lat"] = f["lat_deg"];
      s["lon"] = f["lon_deg"];
      s["alt_msl_m"] = f["alt"].get<double>() / 1000;
      s["alt_rel_m"] = f["relative_alt"].get<double>() / 1000;
      s["heading_deg"] = f["hdg"].get<int>() == 65535 ? Json(nullptr) : Json(f["hdg"].get<double>() / 100);
    } else if (m == "GPS_RAW_INT") {
      s["gps_fix_type"] = f["fix_type"];
      s["satellites"] = f["satellites_visible"];
      if (!s.contains("lat")) s["lat"] = f["lat_deg"];
      if (!s.contains("lon")) s["lon"] = f["lon_deg"];
    } else if (m == "HOME_POSITION") {
      s["home_lat"] = f["latitude_deg"];
      s["home_lon"] = f["longitude_deg"];
      s["home_alt_m"] = f["altitude"].get<double>() / 1000;
    } else if (m == "VFR_HUD") {
      s["groundspeed_mps"] = f["groundspeed"];
      s["airspeed_mps"] = f["airspeed"];
      s["climb_mps"] = f["climb"];
      s["throttle_pct"] = f["throttle"];
    } else if (m == "SYS_STATUS") {
      s["battery_v"] = f["voltage_battery"].get<double>() / 1000;
      if (f.value("battery_remaining", -1) >= 0) s["battery_pct"] = f["battery_remaining"];
    } else if (m == "ATTITUDE") {
      s["roll_rad"] = f["roll"];
      s["pitch_rad"] = f["pitch"];
      s["yaw_rad"] = f["yaw"];
    } else if (m == "STATUSTEXT") {
      if (!s.contains("statustext")) s["statustext"] = Json::array();
      s["statustext"].push_back(f["text"]);
    } else if (m == "RADIO" || m == "RADIO_STATUS") {
      s["radio_rssi"] = f["rssi"];
      s["radio_remrssi"] = f["remrssi"];
      s["radio_noise"] = f["noise"];
    }
  }
  return s;
}

}  // namespace si::mavlink
