// python: receivers/adsb/tracks.py, line by line (dict key order, rounding and float expressions as
// Python evaluates them).
#include "receivers/adsb/tracks.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/pyfmt.hpp"
#include "receivers/adsb/codec.hpp"
#include "receivers/adsb/tables.hpp"

namespace si::rx::adsb {

namespace {

constexpr double kPairS = 10.0;    // max even/odd spacing for a global CPR decode
constexpr double kLocalS = 30.0;   // a fix is a valid local-decode reference for this long
constexpr double kMaxKt = 1000.0;  // faster than this between fixes -> implausible position
constexpr double kPi = 3.141592653589793;
constexpr double kDegToRad = kPi / 180.0, kRadToDeg = 180.0 / kPi;

double radians(double x) { return x * kDegToRad; }
double degrees(double x) { return x * kRadToDeg; }
double py_fmod(double a, double b) {
  double r = std::fmod(a, b);
  if (r != 0.0) {
    if ((b < 0) != (r < 0)) r += b;
  } else {
    r = std::copysign(0.0, b);
  }
  return r;
}

// ---- python dict semantics
bool truthy(const Json& v) {
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number_integer() || v.is_number_unsigned()) return v.get<int64_t>() != 0;
  if (v.is_number_float()) return v.get<double>() != 0.0;
  if (v.is_string()) return !v.get_ref<const std::string&>().empty();
  return !v.empty();
}
Json get(const Json& d, const char* k) {
  auto it = d.find(k);
  return it == d.end() ? Json() : *it;
}
double num(const Json& v) { return v.is_boolean() ? (v.get<bool>() ? 1.0 : 0.0) : v.get<double>(); }
bool py_less(const Json& a, const Json& b) { return num(a) < num(b); }

// dict(Counter(...)) of int keys -> json object (insertion order)
template <class K>
struct Counter {
  std::vector<std::pair<K, int64_t>> v;
  void add(const K& k) {
    for (auto& e : v)
      if (e.first == k) {
        ++e.second;
        return;
      }
    v.emplace_back(k, 1);
  }
};
Json counter_json(const Counter<int64_t>& c) {
  Json o = Json::object();
  for (const auto& [k, n] : c.v) o[std::to_string(k)] = n;
  return o;
}
Json counter_json(const Counter<std::string>& c) {
  Json o = Json::object();
  for (const auto& [k, n] : c.v) o[k] = n;
  return o;
}

Json str_list(const std::set<std::string>& s) {   // sorted(set)
  Json a = Json::array();
  for (const auto& x : s) a.push_back(x);
  return a;
}

using Pos = std::pair<double, double>;

double haversine_km(double lat1, double lon1, double lat2, double lon2) {
  const double la1 = radians(lat1), lo1 = radians(lon1), la2 = radians(lat2), lo2 = radians(lon2);
  const double h = std::pow(std::sin((la2 - la1) / 2), 2.0) +
                   std::cos(la1) * std::cos(la2) * std::pow(std::sin((lo2 - lo1) / 2), 2.0);
  return 2 * 6371.0 * std::asin(std::sqrt(h));
}

double bearing(double lat1, double lon1, double lat2, double lon2) {
  const double la1 = radians(lat1), lo1 = radians(lon1), la2 = radians(lat2), lo2 = radians(lon2);
  const double y = std::sin(lo2 - lo1) * std::cos(la2);
  const double x = std::cos(la1) * std::sin(la2) - std::sin(la1) * std::cos(la2) * std::cos(lo2 - lo1);
  return py_fmod(degrees(std::atan2(y, x)) + 360, 360);
}

struct Ref {
  bool set = false;
  double lat = 0, lon = 0, alt = 0;
  bool has_alt = false;
};

// global CPR from an even/odd pair, then local decoding against the last fix (or the receiver)
std::vector<Json> positions(const std::vector<Json*>& frames, const Ref& ref) {
  std::vector<Json> fixes;
  std::optional<std::pair<double, Pos>> last;
  std::optional<std::pair<std::string, double>> cpr[2];
  for (const Json* mp : frames) {
    const Json& m = *mp;
    const Json tcj = get(m, "typecode");
    const int64_t tc = truthy(tcj) ? tcj.get<int64_t>() : 0;
    if (!m.contains("lat_cpr") ||
        !((9 <= tc && tc <= 18) || (20 <= tc && tc <= 22) || (5 <= tc && tc <= 8 && ref.set)))
      continue;
    const double t = m["t"].get<double>();
    const std::string h = m["hex"].get<std::string>();
    cpr[m["odd"].get<int>()] = std::make_pair(h, t);
    std::optional<Pos> pos;
    std::string method;
    if (last && t - last->first <= kLocalS) {
      pos = cpr_local(h, last->second.first, last->second.second, tc <= 8);
      method = "local (previous fix)";
    } else if (cpr[0] && cpr[1] && std::abs(cpr[0]->second - cpr[1]->second) <= kPairS && tc >= 9) {
      pos = cpr_airborne_global(cpr[0]->first, cpr[1]->first, cpr[0]->second, cpr[1]->second);
      method = "global (even/odd)";
    } else if (ref.set) {
      pos = cpr_local(h, ref.lat, ref.lon, tc <= 8);
      method = "local (receiver)";
    }
    if (!pos) continue;
    Json fix = Json::object();
    fix["t"] = py_round(t, 3);
    fix["lat"] = py_round(pos->first, 5);
    fix["lon"] = py_round(pos->second, 5);
    fix["method"] = method;
    fix["alt_ft"] = m.contains("altitude_ft") ? m["altitude_ft"] : get(m, "gnss_altitude_ft");
    if (last && t > last->first) {
      const double kt = haversine_km(last->second.first, last->second.second, pos->first, pos->second) / 1.852 /
                        ((t - last->first) / 3600);
      if (kt > kMaxKt) {
        fix["implausible_speed_kt"] = py_round_int(kt);
        fixes.push_back(fix);
        continue;   // do not use it as the next reference
      }
    }
    fixes.push_back(fix);
    last = std::make_pair(t, *pos);
  }
  return fixes;
}

// _nearest(items, t, max_dt): the first item closest in time, if within max_dt
template <class V>
const V* nearest(const std::vector<std::pair<double, V>>& items, double t, double max_dt) {
  const std::pair<double, V>* best = nullptr;
  for (const auto& x : items)
    if (!best || std::abs(x.first - t) < std::abs(best->first - t)) best = &x;
  return best && std::abs(best->first - t) <= max_dt ? &best->second : nullptr;
}

using Vel = std::vector<std::pair<double, std::pair<Json, Json>>>;   // (t, (gs, track))
using Alts = std::vector<std::pair<double, Json>>;                   // (t, altitude ft)
struct Cb {
  double t;
  std::string reg;
  Json d;
};

// one register for ambiguous DF20/21 replies (BDS50 vs BDS60) from the aircraft's ADS-B GS/track/alt
std::vector<Cb> resolve_commb(const std::vector<Json*>& frames, const Vel& vel, const Alts& alts) {
  std::vector<Cb> out;
  for (Json* mp : frames) {
    Json& m = *mp;
    const Json cb = get(m, "commb");
    if (!truthy(cb) || !truthy(cb["candidates"])) continue;
    std::vector<std::string> regs;
    for (auto it = cb["candidates"].begin(); it != cb["candidates"].end(); ++it) regs.push_back(it.key());
    Json reg = regs.size() == 1 ? Json(regs[0]) : Json();
    if (reg.is_null() && regs.size() == 2 && std::set<std::string>(regs.begin(), regs.end()) ==
                                                 std::set<std::string>{"BDS50", "BDS60"}) {
      const double t = m["t"].get<double>();
      const auto* v = nearest(vel, t, 10);
      const auto* a = nearest(alts, t, 10);
      if (v && a && !a->is_null())
        reg = commb_is50or60(hex_bits(m["hex"].get<std::string>()), num(v->first), num(v->second), num(*a));
    }
    m["commb"]["register"] = truthy(reg) ? reg : Json("ambiguous");
    if (truthy(reg)) {
      const std::string r = reg.get<std::string>();
      if (!cb["candidates"].contains(r)) throw std::runtime_error("KeyError: '" + r + "'");
      out.push_back({m["t"].get<double>(), r, cb["candidates"][r]});
    }
  }
  return out;
}

// wind from (TAS, heading) vs (GS, track); static air temperature from TAS and Mach
Json met(const std::vector<Cb>& commb, const Vel& vel, const Alts& alts) {
  std::vector<std::pair<double, Json>> tas, hdg, mach;
  for (const auto& c : commb) {
    if (c.reg == "BDS50" && truthy(get(c.d, "tas_kt"))) tas.emplace_back(c.t, c.d["tas_kt"]);
    if (c.reg == "BDS60" && !get(c.d, "magnetic_heading_deg").is_null())
      hdg.emplace_back(c.t, c.d["magnetic_heading_deg"]);
    if (c.reg == "BDS60" && truthy(get(c.d, "mach"))) mach.emplace_back(c.t, c.d["mach"]);
  }
  Json out = Json::object();
  for (const auto& [t, hj] : hdg) {
    const Json* a = nearest(tas, t, 5);
    const auto* g = nearest(vel, t, 5);
    if (a && truthy(*a) && g) {
      const double h = num(hj), av = num(*a);
      const double ax = av * std::sin(radians(h)), ay = av * std::cos(radians(h));
      const double gx = num(g->first) * std::sin(radians(num(g->second)));
      const double gy = num(g->first) * std::cos(radians(num(g->second)));
      const double wx = gx - ax, wy = gy - ay;
      out["wind_kt"] = py_round_int(std::hypot(wx, wy));
      out["wind_from_deg"] = py_round_int(py_fmod(degrees(std::atan2(wx, wy)) + 180, 360));
      out["wind_method"] =
          "ground vector (ADS-B GS/track) - air vector (BDS50 TAS, BDS60 magnetic heading; "
          "declination not applied, ~0-2 deg over India)";
    }
  }
  for (const auto& [t, mcj] : mach) {
    const Json* a = nearest(tas, t, 5);
    const Json* alt = nearest(alts, t, 10);
    const double mc = num(mcj);
    if (a && truthy(*a) && mc > 0.1) {
      const double av = num(*a);
      const double tk = std::pow(av / (mc * 38.967854), 2.0);
      const double sat = tk - 273.15;
      // half-LSB quantisation of TAS (2 kt) and Mach (0.004), doubled by the square
      const double err = tk * 2 * std::hypot(1 / av, 0.002 / mc);
      out["static_air_temp_c"] = py_round(sat, 1);
      out["temp_uncertainty_c"] = py_round(err, 1);
      out["temp_method"] = "SAT = (TAS / (Mach * 38.97))^2, BDS50 TAS + BDS60 Mach";
      if (alt && !alt->is_null() && num(*alt) < 36089)
        out["isa_deviation_c"] = py_round(sat - (15 - 1.98 * num(*alt) / 1000), 1);
    }
  }
  for (const auto& c : commb) {   // MRAR (BDS44): reported directly
    if (c.reg != "BDS44") continue;
    if (truthy(get(c.d, "wind"))) {
      out["wind_kt"] = c.d["wind"][0];
      out["wind_from_deg"] = c.d["wind"][1];
      out["wind_method"] = "BDS44 (reported)";
    }
    if (!get(c.d, "temperature_c").is_null()) {
      out["static_air_temp_c"] = c.d["temperature_c"];
      out["temp_method"] = "BDS44 (reported)";
    }
  }
  return out;
}

Json last_of(const std::vector<Json*>& frames, const char* key) {
  for (auto it = frames.rbegin(); it != frames.rend(); ++it) {
    const Json v = get(**it, key);
    if (!v.is_null()) return v;
  }
  return Json();
}

bool is_digits(const std::string& s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; });
}

// python repr of a sorted list of strings
std::string repr_list(const std::set<std::string>& s) { return py_repr(str_list(s)); }

Json one(const std::string& icao, std::vector<Json*> frames, const Ref& ref) {
  std::stable_sort(frames.begin(), frames.end(),
                   [](const Json* a, const Json* b) { return (*a)["t"].get<double>() < (*b)["t"].get<double>(); });
  Vel vel;
  Alts alts;
  for (const Json* m : frames) {
    if (get(*m, "typecode") == 19 && !get(*m, "ground_speed_kt").is_null())
      vel.emplace_back((*m)["t"].get<double>(), std::make_pair((*m)["ground_speed_kt"], (*m)["track_deg"]));
    if (!get(*m, "altitude_ft").is_null()) alts.emplace_back((*m)["t"].get<double>(), (*m)["altitude_ft"]);
  }
  Json a = Json::object();
  a["icao"] = icao;
  const auto ctry = country(icao);
  a["country"] = ctry ? Json(*ctry) : Json();
  a["first_t"] = py_round(frames.front()->at("t").get<double>(), 3);
  a["last_t"] = py_round(frames.back()->at("t").get<double>(), 3);
  a["n_frames"] = frames.size();
  Counter<int64_t> dfs, tcs;
  for (const Json* m : frames) {
    dfs.add((*m)["df"].get<int64_t>());
    if (truthy(get(*m, "typecode"))) tcs.add((*m)["typecode"].get<int64_t>());
  }
  a["dfs"] = counter_json(dfs);
  a["typecodes"] = counter_json(tcs);
  std::vector<std::string> flags;

  std::vector<std::string> cs;
  for (const Json* m : frames)
    if (truthy(get(*m, "callsign"))) cs.push_back((*m)["callsign"].get<std::string>());
  if (!cs.empty()) {
    const std::string& c = cs.back();
    a["callsign"] = c;
    const auto op = operator_of(c.substr(0, 3));
    if (op && is_digits(c.size() > 3 ? c.substr(3) : "")) {
      a["operator"] = op->first;
      a["flight"] = op->second + c.substr(3);
    }
    const std::set<std::string> u(cs.begin(), cs.end());
    if (u.size() > 1) flags.push_back("callsign changed: " + repr_list(u));
  }
  for (const char* k : {"category", "capability", "control_field", "adsb_version", "squawk", "emergency"}) {
    const Json v = last_of(frames, k);
    if (!v.is_null()) a[k] = v;
  }
  std::set<std::string> sq;
  for (const Json* m : frames)
    if (truthy(get(*m, "squawk"))) sq.insert((*m)["squawk"].get<std::string>());
  if (sq.size() > 1) flags.push_back("squawk changed: " + repr_list(sq));
  std::set<std::string> em;
  for (const char* e : {"7500", "7600", "7700"})
    if (sq.count(e)) em.insert(e);
  if (!em.empty()) flags.push_back("emergency squawk " + repr_list(em));
  const Json emg = get(a, "emergency");
  if (!emg.is_null() && emg != "no emergency") flags.push_back("emergency: " + py_str(emg));
  if (get(a, "category") == Json("unmanned aerial vehicle")) {  // emitter category B6 (TC 3, CA 6): a drone on ADS-B
    a["uav"] = true;
    flags.push_back("UAV (ADS-B emitter category B6)");
  }

  // kinematics
  if (!alts.empty()) {
    const Json* mn = &alts[0].second;
    const Json* mx = &alts[0].second;
    for (const auto& x : alts) {
      if (py_less(x.second, *mn)) mn = &x.second;
      if (py_less(*mx, x.second)) mx = &x.second;
    }
    a["altitude_ft"] = Json::object({{"first", alts.front().second}, {"last", alts.back().second}, {"min", *mn}, {"max", *mx}});
    if (alts.size() > 1 && alts.back().first > alts.front().first)
      a["altitude_rate_fpm"] = py_round_int((num(alts.back().second) - num(alts.front().second)) /
                                            (alts.back().first - alts.front().first) * 60);
    const double vr = a.contains("altitude_rate_fpm") ? num(a["altitude_rate_fpm"]) : 0;
    a["vertical_trend"] = vr > 300 ? "climbing" : vr < -300 ? "descending" : "level";
  }
  for (const char* k : {"gnss_altitude_ft", "gnss_minus_baro_ft", "ground_speed_kt", "track_deg", "vrate_fpm",
                        "vrate_source", "ias_kt", "tas_kt", "heading_deg", "heading_source"}) {
    const Json v = last_of(frames, k);
    if (!v.is_null()) a[k] = v;
  }
  if (!vel.empty()) {
    const Json* mn = &vel[0].second.first;
    const Json* mx = &vel[0].second.first;
    for (const auto& x : vel) {
      if (py_less(x.second.first, *mn)) mn = &x.second.first;
      if (py_less(*mx, x.second.first)) mx = &x.second.first;
    }
    a["ground_speed_range_kt"] = Json::array({*mn, *mx});
    const double turn = py_fmod(num(vel.back().second.second) - num(vel.front().second.second) + 540, 360) - 180;
    if (std::abs(turn) >= 5) a["turn_deg"] = py_round(turn, 1);   // + right, - left, over the observation
  }
  const double alt_last = alts.empty() ? 0 : num(alts.back().second);
  const Json gmb = get(a, "gnss_minus_baro_ft");
  if (std::abs(truthy(gmb) ? num(gmb) : 0) > std::max(1500.0, 0.08 * alt_last))   // beyond what temperature explains
    flags.push_back("GNSS and baro altitude differ by " + py_str(gmb) + " ft");

  // position
  const std::vector<Json> fixes = positions(frames, ref);
  std::vector<const Json*> good;
  for (const Json& f : fixes)
    if (!f.contains("implausible_speed_kt")) good.push_back(&f);
  if (!good.empty()) {
    const Json& g = *good.back();
    a["position"] = g;
    a["n_positions"] = good.size();
    Json tr = Json::array();
    for (const Json* f : good) tr.push_back(Json::array({(*f)["t"], (*f)["lat"], (*f)["lon"], (*f)["alt_ft"]}));
    a["track"] = tr;
    if (ref.set) {
      const double plat = g["lat"].get<double>(), plon = g["lon"].get<double>();
      const double d = haversine_km(ref.lat, ref.lon, plat, plon);
      a["range_km"] = py_round(d, 2);
      a["bearing_deg"] = py_round(bearing(ref.lat, ref.lon, plat, plon), 1);
      if (!g["alt_ft"].is_null()) {
        const double dz = num(g["alt_ft"]) * 0.3048 - (ref.has_alt ? ref.alt : 0);
        a["elevation_deg"] = py_round(degrees(std::atan2(dz, d * 1000)), 1);
      }
    }
  }
  if (fixes.size() > good.size())
    flags.push_back(fmt("%zu position(s) implausible (> %.0f kt jump): possible spoofing or bad CPR",
                        fixes.size() - good.size(), kMaxKt));

  // intent (TC29, BDS40)
  for (const char* k : {"selected_altitude_ft", "selected_altitude_source", "baro_setting_hpa", "selected_heading_deg",
                        "modes", "tcas_operational"}) {
    const Json v = last_of(frames, k);
    if (!v.is_null()) {
      if (!a.contains("intent")) a["intent"] = Json::object();
      a["intent"][k] = v;
    }
  }

  // integrity and capabilities (TC29 / TC31 / DF0)
  Json integ = Json::object();
  for (const char* k : {"nac_p", "nic_a", "nic_b", "nic_baro", "sil", "sil_per", "gva", "nac_v", "sda"}) {
    const Json v = last_of(frames, k);
    if (!v.is_null()) integ[k] = v;
  }
  if (integ.contains("nac_p")) {
    const auto epu = integ["nac_p"].is_number_integer() ? nacp_epu_m(integ["nac_p"].get<int>()) : std::nullopt;
    integ["position_accuracy_m"] = "< " + (epu ? std::to_string(*epu) : std::string("unknown"));
  }
  if (!integ.empty()) a["integrity"] = integ;
  const Json lc = last_of(frames, "capabilities");
  Json cap = truthy(lc) ? lc : Json::object();
  for (const char* k : {"single_antenna", "heading_reference", "acas_reply_info", "acas_sensitivity", "crosslink_capable"}) {
    const Json v = last_of(frames, k);
    if (!v.is_null()) cap[k] = v;
  }
  if (!cap.empty()) a["capabilities"] = cap;

  // Comm-B (Mode S ground radar EHS / MRAR) and derived met
  std::vector<Json*> fr(frames.begin(), frames.end());
  const std::vector<Cb> commb = resolve_commb(fr, vel, alts);
  if (!commb.empty()) {
    Counter<std::string> regs;
    for (const auto& c : commb) regs.add(c.reg);
    Json cb = Json::object();
    cb["registers"] = counter_json(regs);
    for (const auto& c : commb)
      for (auto it = c.d.begin(); it != c.d.end(); ++it)
        if (!it.value().is_null()) cb[it.key()] = it.value();
    a["commb"] = cb;
    const Json mt = met(commb, vel, alts);
    if (!mt.empty()) a["met"] = mt;
    std::vector<std::string> bdscs;
    for (const auto& c : commb)
      if (c.reg == "BDS20" && truthy(get(c.d, "callsign"))) {
        std::string s = c.d["callsign"].get<std::string>();
        size_t i = 0, j = s.size();
        while (i < j && (s[i] == '_' || s[i] == ' ')) ++i;
        while (j > i && (s[j - 1] == '_' || s[j - 1] == ' ')) --j;
        bdscs.push_back(s.substr(i, j - i));
      }
    if (!cs.empty() && !bdscs.empty() && bdscs.back() != cs.back())
      flags.push_back("Comm-B callsign " + bdscs.back() + " != ADS-B callsign " + cs.back());
  }

  // who is interrogating this aircraft
  int64_t acas = 0, sq1718 = 0;
  std::vector<const Json*> radar;
  for (const Json* m : frames) {
    const int64_t df = (*m)["df"].get<int64_t>();
    if (df == 0 || df == 16) ++acas;
    if (df == 4 || df == 5 || df == 20 || df == 21) radar.push_back(m);
    if (df == 17 || df == 18) ++sq1718;
  }
  std::vector<std::vector<double>> dwells;   // radar beam dwells: replies within 0.2 s
  for (const Json* m : radar) {
    const double t = (*m)["t"].get<double>();
    if (dwells.empty() || t - dwells.back().back() > 0.2) dwells.push_back({t});
    else dwells.back().push_back(t);
  }
  std::set<int64_t> iis;
  for (const Json* m : radar)
    if (get(*m, "ids") == 1) iis.insert((*m)["iis"].get<int64_t>());
  Counter<std::string> d11;
  for (const Json* m : frames)
    if (truthy(get(*m, "interrogator"))) d11.add((*m)["interrogator"].get<std::string>());
  Json in = Json::object();
  in["acas_replies"] = acas;
  in["mode_s_radar_replies"] = radar.size();
  in["radar_dwells"] = dwells.size();
  in["radar_iis"] = Json::array();
  for (int64_t x : iis) in["radar_iis"].push_back(x);
  in["df11_interrogators"] = counter_json(d11);
  in["squitters_df17_18"] = sq1718;
  if (dwells.size() >= 2) {
    std::vector<double> gaps;
    for (size_t i = 0; i + 1 < dwells.size(); ++i) gaps.push_back(dwells[i + 1][0] - dwells[i][0]);
    std::sort(gaps.begin(), gaps.end());
    in["radar_revisit_s"] = py_round(gaps[gaps.size() / 2], 2);   // antenna rotation period (or a multiple)
  }
  a["interrogation"] = in;
  for (const Json* m : frames)
    if (truthy(get(*m, "acas_ra_report")) || truthy(get(*m, "tcas_ra_active"))) {
      flags.push_back("ACAS resolution advisory reported");
      break;
    }

  std::vector<double> snr;
  for (const Json* m : frames)
    if (!get(*m, "snr_db").is_null()) snr.push_back(num((*m)["snr_db"]));
  if (!snr.empty()) {
    std::sort(snr.begin(), snr.end());
    a["signal"] = Json::object({{"snr_db_min", snr.front()}, {"snr_db_median", snr[snr.size() / 2]}, {"snr_db_max", snr.back()}});
  }
  const double span = frames.back()->at("t").get<double>() - frames.front()->at("t").get<double>();
  if (span > 1) a["frame_rate_hz"] = py_round(static_cast<double>(frames.size()) / span, 2);
  if (!flags.empty()) {
    Json f = Json::array();
    for (const auto& s : flags) f.push_back(s);
    a["flags"] = f;
  }
  return a;
}

}  // namespace

Json tracks(std::vector<Json>& msgs, const Json& ref_json) {
  Ref ref;
  if (ref_json.is_array() && ref_json.size() >= 2) {
    ref.set = true;
    ref.lat = num(ref_json[0]);
    ref.lon = num(ref_json[1]);
    if (ref_json.size() > 2) {
      ref.has_alt = true;
      ref.alt = num(ref_json[2]);
    }
  }
  std::vector<std::pair<std::string, std::vector<Json*>>> by;
  for (Json& m : msgs) {
    const Json ic = get(m, "icao");
    if (!truthy(ic)) continue;
    const std::string k = ic.get<std::string>();
    auto it = std::find_if(by.begin(), by.end(), [&](const auto& e) { return e.first == k; });
    if (it == by.end()) by.emplace_back(k, std::vector<Json*>{&m});
    else it->second.push_back(&m);
  }
  std::vector<Json> out;
  for (auto& [k, v] : by) out.push_back(one(k, v, ref));
  std::stable_sort(out.begin(), out.end(), [](const Json& a, const Json& b) {
    return a["n_frames"].get<int64_t>() > b["n_frames"].get<int64_t>();
  });
  Json arr = Json::array();
  for (auto& a : out) arr.push_back(std::move(a));
  return arr;
}

}  // namespace si::rx::adsb
