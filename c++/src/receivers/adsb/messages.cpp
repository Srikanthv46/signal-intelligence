// python: receivers/adsb/messages.py. The encodings pyModeS decodes for the Python package are
// implemented from the standards in codec.cpp.
#include "receivers/adsb/messages.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/pyfmt.hpp"
#include "receivers/adsb/codec.hpp"
#include "receivers/adsb/tables.hpp"

namespace si::rx::adsb {

namespace {

// Comm-B register -> (field, register field function), python COMMB
const std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>>& commb_table() {
  static const std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>> t = {
      {"BDS20", {{"callsign", "cs20"}}},
      {"BDS40", {{"selected_alt_mcp_ft", "selalt40mcp"}, {"selected_alt_fms_ft", "selalt40fms"},
                 {"baro_setting_hpa", "p40baro"}}},
      {"BDS44", {{"wind", "wind44"}, {"temperature_c", "temp44"}, {"pressure_hpa", "p44"}, {"humidity_pct", "hum44"},
                 {"turbulence", "turb44"}}},
      {"BDS45", {{"turbulence", "turb45"}, {"wind_shear", "ws45"}, {"microburst", "mb45"}, {"icing", "ic45"},
                 {"wake_vortex", "wv45"}, {"temperature_c", "temp45"}, {"pressure_hpa", "p45"},
                 {"radio_height_ft", "rh45"}}},
      {"BDS50", {{"roll_deg", "roll50"}, {"true_track_deg", "trk50"}, {"ground_speed_kt", "gs50"},
                 {"track_rate_dps", "rtrk50"}, {"tas_kt", "tas50"}}},
      {"BDS60", {{"magnetic_heading_deg", "hdg60"}, {"ias_kt", "ias60"}, {"mach", "mach60"},
                 {"baro_vrate_fpm", "vr60baro"}, {"inertial_vrate_fpm", "vr60ins"}}},
  };
  return t;
}

// python _r: round(v, 3) for floats, anything else unchanged
Json r3(const Json& v) { return v.is_number_float() ? Json(py_round(v.get<double>(), 3)) : v; }

int64_t u(const std::string& b, size_t i, size_t j) { return ubits(b, i, j); }
bool on(const std::string& b, size_t i) { return b[i] == '1'; }

std::string strip_us(const std::string& s) {   // str.strip("_ ")
  size_t i = 0, j = s.size();
  while (i < j && (s[i] == '_' || s[i] == ' ')) ++i;
  while (j > i && (s[j - 1] == '_' || s[j - 1] == ' ')) --j;
  return s.substr(i, j - i);
}

Json cap_text(int ca) {
  auto c = capability(ca);
  return c ? Json(*c) : Json(fmt("CA %d", ca));
}

Json commb(const std::string& b) {
  const Json inferred = commb_infer(b);
  Json out = Json::object();
  out["bds"] = inferred;
  out["candidates"] = Json::object();
  const std::string s = inferred.is_string() ? inferred.get<std::string>() : "";
  size_t p = 0;
  while (true) {
    const size_t q = s.find(',', p);
    const std::string reg = s.substr(p, q == std::string::npos ? std::string::npos : q - p);
    if (reg == "BDS17") {
      out["candidates"][reg] = Json::object({{"supported_registers", commb_cap17(b)}});
    } else {
      const auto& t = commb_table();
      auto it = std::find_if(t.begin(), t.end(), [&](const auto& e) { return e.first == reg; });
      if (it != t.end()) {
        Json dec = Json::object();
        for (const auto& [k, fn] : it->second) {
          const Json v = commb_field(fn, b);
          if (v.is_null()) continue;
          if (v.is_array()) {
            Json l = Json::array();
            for (const auto& x : v) l.push_back(r3(x));
            dec[k] = l;
          } else {
            dec[k] = r3(v);
          }
        }
        out["candidates"][reg] = dec;
      } else if (!reg.empty()) {
        out["candidates"][reg] = Json::object();
      }
    }
    if (q == std::string::npos) break;
    p = q + 1;
  }
  return out;
}

void me_position(const std::string& b, Json& m, int tc) {
  const std::string me = b.substr(32, 56);
  m["odd"] = static_cast<int>(on(me, 21));
  m["lat_cpr"] = u(me, 22, 39);
  m["lon_cpr"] = u(me, 39, 56);
  m["time_sync"] = static_cast<int>(on(me, 20));
  if (tc >= 9) {
    m["surveillance_status"] = surveillance_status(static_cast<int>(u(me, 5, 7)));
    m["nic_b"] = static_cast<int>(on(me, 7));
  }
}

// TC29 subtype 1 (ADS-B v2), DO-260B 2.2.3.2.7.1
void me_target_state(const std::string& b, Json& m) {
  const std::string me = b.substr(32, 56);
  m["tss_subtype"] = u(me, 5, 7);
  if (u(me, 5, 7) != 1) return;
  const int64_t alt = u(me, 9, 20);
  if (alt) {
    m["selected_altitude_ft"] = (alt - 1) * 32;
    m["selected_altitude_source"] = on(me, 8) ? "FMS" : "MCP/FCU";
  }
  const int64_t baro = u(me, 20, 29);
  if (baro) m["baro_setting_hpa"] = py_round(static_cast<double>(baro - 1) * 0.8 + 800, 1);
  if (on(me, 29)) m["selected_heading_deg"] = py_round(static_cast<double>(u(me, 30, 39) * 180) / 256, 1);
  m["nac_p"] = u(me, 39, 43);
  m["nic_baro"] = static_cast<int>(on(me, 43));
  m["sil"] = u(me, 44, 46);
  if (on(me, 46)) {   // mode bits valid
    Json modes = Json::array();
    const std::pair<const char*, size_t> mb[] = {{"autopilot", 47}, {"VNAV", 48}, {"altitude hold", 49},
                                                 {"approach", 51}, {"LNAV", 53}};
    for (const auto& [n, i] : mb)
      if (on(me, i)) modes.push_back(n);
    m["modes"] = modes;
  }
  m["tcas_operational"] = on(me, 52);
}

// TC31 aircraft operational status, DO-260B 2.2.3.2.7.2
void me_status(const std::string& b, Json& m) {
  const std::string me = b.substr(32, 56);
  const int64_t st = u(me, 5, 8);
  m["status_subtype"] = st == 0 ? Json("airborne") : st == 1 ? Json("surface") : Json(st);
  const int64_t ver = u(me, 40, 43);
  m["adsb_version"] = ver;
  m["nic_a"] = static_cast<int>(on(me, 43));
  m["nac_p"] = u(me, 44, 48);
  m["sil"] = u(me, 50, 52);
  m["heading_reference"] = on(me, 53) ? "magnetic north" : "true north";
  m["sil_per"] = on(me, 54) ? "sample" : "hour";
  if (ver < 2) return;
  if (st == 0) {
    m["gva"] = u(me, 48, 50);
    m["nic_baro"] = static_cast<int>(on(me, 52));
    Json c = Json::object();
    c["tcas_operational"] = on(me, 10);
    c["1090es_in"] = on(me, 11);
    c["air_referenced_velocity"] = on(me, 14);
    c["target_state_reports"] = on(me, 15);
    c["trajectory_change_reports"] = u(me, 16, 18);
    c["uat_in"] = on(me, 18);
    m["capabilities"] = c;
  } else {
    m["length_width_code"] = u(me, 20, 24);
    Json c = Json::object();
    c["1090es_in"] = on(me, 11);
    c["uat_in"] = on(me, 15);
    m["capabilities"] = c;
  }
  m["tcas_ra_active"] = on(me, 26);
  m["ident"] = on(me, 27);
  m["single_antenna"] = on(me, 29);
  m["sda"] = u(me, 30, 32);
}

void extended_squitter(const std::string& h, const std::string& b, Json& m) {
  const int df = m["df"].get<int>();
  m["icao"] = h.substr(2, 6);
  const int ca = static_cast<int>(u(b, 5, 8));
  if (df == 17) m["capability"] = cap_text(ca);
  else m["control_field"] = control_field(ca);
  const int tc = es_typecode(b);
  m["typecode"] = tc;
  const std::string me = b.substr(32, 56);
  if (1 <= tc && tc <= 4) {
    m["callsign"] = strip_us(es_callsign(b));
    auto c = category(tc, static_cast<int>(u(me, 5, 8)));
    m["category"] = c ? *c : "reserved";
  } else if (5 <= tc && tc <= 8) {
    const auto mov = es_surface_movement(b);
    m["ground_speed_kt"] = r3(mov.first);
    m["track_deg"] = r3(mov.second);
    me_position(b, m, tc);
  } else if ((9 <= tc && tc <= 18) || (20 <= tc && tc <= 22)) {
    m[tc <= 18 ? "altitude_ft" : "gnss_altitude_ft"] = es_altitude(b, tc);
    me_position(b, m, tc);
  } else if (tc == 19) {
    const Velocity v = es_velocity(b);
    m["velocity_subtype"] = u(me, 5, 8);
    m["nac_v"] = u(me, 10, 13);
    m["intent_change"] = on(me, 8);
    if (v.kind == "GS") {
      m["ground_speed_kt"] = r3(v.spd);
      m["track_deg"] = r3(v.ang);
    } else {
      m[v.kind == "IAS" ? "ias_kt" : "tas_kt"] = v.spd;
      if (!v.ang.is_null()) m["heading_deg"] = r3(v.ang);
      m["heading_source"] = v.ang_src;
    }
    m["vrate_fpm"] = v.vr;
    m["vrate_source"] = v.vr_src;
    m["gnss_minus_baro_ft"] = es_gnss_baro_diff(b);
  } else if (tc == 28) {
    const int64_t sub = u(me, 5, 8);
    m["status_subtype"] = sub;
    if (sub == 1) {
      m["emergency"] = emergency(static_cast<int>(u(me, 8, 11)));
      m["squawk"] = es_emergency_squawk(b);
    } else if (sub == 2) {
      m["tcas_ra_broadcast"] = me.substr(8, 48);
    }
  } else if (tc == 29) {
    me_target_state(b, m);
  } else if (tc == 31) {
    me_status(b, m);
  }
}

void decode(Json& m) {
  const std::string h = m["hex"].get<std::string>();
  const int df = m["df"].get<int>();
  const std::string b = hex_bits(h);
  if (df == 17 || df == 18) {
    extended_squitter(h, b, m);
    return;
  }
  m["icao"] = m["aa"];
  auto ac = [&]() -> Json {
    const auto a = altitude_ac13(b.substr(19, 13));
    return a ? Json(*a) : Json();
  };
  if (df == 11) {
    m["capability"] = cap_text(static_cast<int>(u(b, 5, 8)));
    const int64_t res = std::stoll(m["parity_residual"].get<std::string>(), nullptr, 16);
    const int64_t cl = res >> 4, ic = res & 0xF;   // code label + interrogator code
    if (cl == 0)
      m["interrogator"] = fmt("II %lld", static_cast<long long>(ic)) + (ic == 0 ? " (or acquisition squitter)" : "");
    else
      m["interrogator"] = fmt("SI %lld", static_cast<long long>(ic + 16 * (cl - 1)));
    m["interrogator_code"] = res;
  } else if (df == 0 || df == 16) {   // ACAS (TCAS) air-air replies
    m["vertical_status"] = on(b, 5) ? "on ground" : "airborne";
    m["acas_sensitivity"] = acas_sl(static_cast<int>(u(b, 8, 11)));
    const int ri = static_cast<int>(u(b, 13, 17));
    auto r = acas_ri(ri);
    m["acas_reply_info"] = r ? *r : fmt("RI %d", ri);
    if (df == 0) m["crosslink_capable"] = on(b, 6);
    m["altitude_ft"] = ac();
    if (df == 16) {
      m["mv"] = h.substr(8, 14);
      if (h.compare(8, 2, "30") == 0) m["acas_ra_report"] = true;
    }
  } else if (df == 4 || df == 5 || df == 20 || df == 21) {   // replies to Mode S ground radar
    m["flight_status"] = flight_status(static_cast<int>(u(b, 5, 8)));
    m["downlink_request"] = u(b, 8, 13);
    const int64_t um = u(b, 13, 19);
    m["iis"] = um >> 2;
    m["ids"] = um & 3;
    if (df == 4 || df == 20) m["altitude_ft"] = ac();
    else m["squawk"] = squawk_id13(b.substr(19, 13));
    if (df == 20 || df == 21) m["commb"] = commb(b);
  }
}

}  // namespace

void fields(std::vector<Json>& msgs) {
  for (Json& m : msgs) {
    try {
      decode(m);
    } catch (const std::exception& e) {   // keep the raw frame
      m["field_error"] = std::string(e.what()).substr(0, 80);
    }
  }
}

}  // namespace si::rx::adsb
