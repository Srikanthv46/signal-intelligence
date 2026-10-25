// Mode S / ADS-B field encodings from ICAO Annex 10 Vol IV, DO-260B, ICAO Doc 9871 and the ISA
// (ICAO Doc 7488). See codec.hpp. Floating-point expressions keep the evaluation order of the Python
// reference output they are checked against (int * k / n as (int * k) / n, etc.).
#include "receivers/adsb/codec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace si::rx::adsb {

namespace {

constexpr double kPi = 3.141592653589793;
constexpr double kDegToRad = kPi / 180.0;   // python math.radians / np.radians
constexpr double kRadToDeg = 180.0 / kPi;   // math.degrees / np.degrees

// python float % (floored, result takes the divisor's sign)
double py_fmod(double a, double b) {
  double r = std::fmod(a, b);
  if (r != 0.0) {
    if ((b < 0) != (r < 0)) r += b;
  } else {
    r = std::copysign(0.0, b);
  }
  return r;
}
int64_t py_imod(int64_t a, int64_t b) {
  int64_t r = a % b;
  if (r != 0 && ((r < 0) != (b < 0))) r += b;
  return r;
}
// python float // (CPython float_floor_div: consistent with %, so 80.0 // (360 / 54) == 11.0)
double py_floordiv(double vx, double wx) {
  double mod = std::fmod(vx, wx);
  double div = (vx - mod) / wx;
  if (mod != 0.0) {
    if ((wx < 0) != (mod < 0)) div -= 1.0;
  }
  if (div != 0.0) {
    double fd = std::floor(div);
    if (div - fd > 0.5) fd += 1.0;
    return fd;
  }
  return std::copysign(0.0, vx / wx);
}
int64_t floor_i(double x) { return static_cast<int64_t>(std::floor(x)); }

bool bit(const std::string& d, size_t i) { return d[i] == '1'; }

// Gray code -> binary
int64_t gray2int(int64_t g) {
  g ^= g >> 8;
  g ^= g >> 4;
  g ^= g >> 2;
  g ^= g >> 1;
  return g;
}

// ------------------------------------------------------------------ ISA (ICAO Doc 7488)
constexpr double kKts = 0.514444, kFt = 0.3048;
constexpr double kR = 287.05287, kGamma = 1.40, kP0 = 101325.0, kRho0 = 1.225;

struct Atmos {
  double p, rho, T;
};
// troposphere lapse -0.0065 K/m to 11 km, isothermal 216.65 K above
Atmos atmos(double H) {
  const double T = std::max(288.15 - 0.0065 * H, 216.65);
  const double rhotrop = 1.225 * std::pow(T / 288.15, 4.256848030018761);
  const double dhstrat = std::max(0.0, H - 11000.0);
  const double rho = rhotrop * std::exp(-dhstrat / 6341.552161);
  const double p = rho * kR * T;
  return {p, rho, T};
}
double mach2tas(double mach, double H) {
  const double a = std::sqrt(kGamma * kR * atmos(H).T);
  return mach * a;
}
// subsonic compressible pitot relations (Saint-Venant) between TAS and CAS
double tas2cas(double tas, double H) {
  const Atmos A = atmos(H);
  const double qdyn = A.p * (std::pow(1.0 + A.rho * tas * tas / (7.0 * A.p), 3.5) - 1.0);
  double cas = std::sqrt(7.0 * kP0 / kRho0 * (std::pow(qdyn / kP0 + 1.0, 2.0 / 7.0) - 1.0));
  return tas < 0 ? -1 * cas : cas;
}
double cas2tas(double cas, double H) {
  const Atmos A = atmos(H);
  const double qdyn = kP0 * (std::pow(1.0 + kRho0 * cas * cas / (7.0 * kP0), 3.5) - 1.0);
  double tas = std::sqrt(7.0 * A.p / A.rho * (std::pow(1.0 + qdyn / A.p, 2.0 / 7.0) - 1.0));
  return cas < 0 ? -1 * tas : tas;
}

// ------------------------------------------------------------------ Comm-B registers (Doc 9871)
std::string mb_of(const std::string& bits) { return bits.substr(32, 56); }

// status bit sb clear but the field msb..lsb (1-based) not zero
bool wrong_status(const std::string& d, size_t sb, size_t msb, size_t lsb) {
  return d[sb - 1] == '0' && ubits(d, msb - 1, lsb) != 0;
}

int64_t signed_field(const std::string& d, size_t sign_bit, size_t i, size_t j) {
  int64_t v = ubits(d, i, j);
  if (bit(d, sign_bit)) v -= int64_t(1) << (j - i);
  return v;
}

// BDS 2,0 aircraft identification (8 characters, '#' kept)
std::string cs20(const std::string& d) {
  std::string cs;
  for (size_t i = 8; i < 56; i += 6) cs += callsign_char(static_cast<int>(ubits(d, i, i + 6)));
  return cs;
}

// field functions (null when the status bit says "not available")
Json field(const std::string& fn, const std::string& d) {
  auto st = [&](size_t i) { return bit(d, i); };
  if (fn == "cs20") return cs20(d);
  // BDS 4,0 selected vertical intention
  if (fn == "selalt40mcp") return st(0) ? Json(ubits(d, 1, 13) * 16) : Json();
  if (fn == "selalt40fms") return st(13) ? Json(ubits(d, 14, 26) * 16) : Json();
  if (fn == "p40baro") return st(26) ? Json(static_cast<double>(ubits(d, 27, 39)) * 0.1 + 800) : Json();
  // BDS 4,4 meteorological routine air report
  if (fn == "wind44") {
    if (!st(4)) return Json();
    return Json::array({ubits(d, 5, 14), static_cast<double>(ubits(d, 14, 23) * 180) / 256});
  }
  if (fn == "temp44") return static_cast<double>(signed_field(d, 23, 24, 34)) * 0.25;
  if (fn == "p44") return st(34) ? Json(ubits(d, 35, 46)) : Json();
  if (fn == "hum44") return st(49) ? Json(static_cast<double>(ubits(d, 50, 56) * 100) / 64) : Json();
  if (fn == "turb44") return st(46) ? Json(ubits(d, 47, 49)) : Json();
  // BDS 4,5 meteorological hazard report
  if (fn == "turb45") return st(0) ? Json(ubits(d, 1, 3)) : Json();
  if (fn == "ws45") return st(3) ? Json(ubits(d, 4, 6)) : Json();
  if (fn == "mb45") return st(6) ? Json(ubits(d, 7, 9)) : Json();
  if (fn == "ic45") return st(9) ? Json(ubits(d, 10, 12)) : Json();
  if (fn == "wv45") return st(12) ? Json(ubits(d, 13, 15)) : Json();
  if (fn == "temp45") return static_cast<double>(signed_field(d, 16, 17, 26)) * 0.25;   // status bit not checked
  if (fn == "p45") return st(26) ? Json(ubits(d, 27, 38)) : Json();
  if (fn == "rh45") return st(38) ? Json(ubits(d, 39, 51) * 16) : Json();
  // BDS 5,0 track and turn report
  if (fn == "roll50") return st(0) ? Json(static_cast<double>(signed_field(d, 1, 2, 11) * 45) / 256) : Json();
  if (fn == "trk50") {
    if (!st(11)) return Json();
    double t = static_cast<double>(signed_field(d, 12, 13, 23) * 90) / 512;
    return t < 0 ? t + 360 : t;
  }
  if (fn == "gs50") return st(23) ? Json(ubits(d, 24, 34) * 2) : Json();
  if (fn == "rtrk50") return st(34) ? Json(static_cast<double>(signed_field(d, 35, 36, 45) * 8) / 256) : Json();
  if (fn == "tas50") return st(45) ? Json(ubits(d, 46, 56) * 2) : Json();
  // BDS 6,0 heading and speed report
  if (fn == "hdg60") {
    if (!st(0)) return Json();
    double h = static_cast<double>(signed_field(d, 1, 2, 12) * 90) / 512;
    return h < 0 ? h + 360 : h;
  }
  if (fn == "ias60") return st(12) ? Json(ubits(d, 13, 23)) : Json();
  if (fn == "mach60") return st(23) ? Json(static_cast<double>(ubits(d, 24, 34)) * 2.048 / 512) : Json();
  if (fn == "vr60baro" || fn == "vr60ins") {
    const size_t s = fn == "vr60baro" ? 34 : 45;
    if (!st(s)) return Json();
    return signed_field(d, s + 1, s + 2, s + 11) * 32;
  }
  return Json();
}

bool is_num(const Json& v) { return v.is_number(); }
double num(const Json& v) { return v.get<double>(); }

bool all_zero(const std::string& d) { return d.find('1') == std::string::npos; }

bool is10(const std::string& d) {
  if (d.compare(0, 8, "00010000") != 0) return false;
  if (ubits(d, 9, 14) != 0) return false;
  const int64_t ov = ubits(d, 16, 23);
  if (bit(d, 14) && ov < 5) return false;
  if (!bit(d, 14) && ov > 4) return false;
  return true;
}

const std::array<const char*, 24> kCap17 = {"BDS05", "BDS06", "BDS07", "BDS08", "BDS09", "BDS0A", "BDS20", "BDS21",
                                            "BDS40", "BDS41", "BDS42", "BDS43", "BDS44", "BDS45", "BDS48", "BDS50",
                                            "BDS51", "BDS52", "BDS53", "BDS54", "BDS55", "BDS56", "BDS5F", "BDS60"};

bool is17(const std::string& d) {
  if (ubits(d, 24, 56) != 0) return false;
  return bit(d, 6);  // BDS 2,0 must be supported
}

bool is20(const std::string& d) {
  if (d.compare(0, 8, "00100000") != 0) return false;
  return cs20(d).find('#') == std::string::npos;
}

bool is30(const std::string& d) {
  if (d.compare(0, 8, "00110000") != 0) return false;
  if (d.compare(28, 2, "11") == 0) return false;   // threat type 3 not assigned
  if (ubits(d, 15, 22) >= 48) return false;        // ACAS RA bits reserved
  return true;
}

bool is40(const std::string& d) {
  if (wrong_status(d, 1, 2, 13) || wrong_status(d, 14, 15, 26) || wrong_status(d, 27, 28, 39) ||
      wrong_status(d, 48, 49, 51) || wrong_status(d, 54, 55, 56))
    return false;
  if (ubits(d, 39, 47) != 0 || ubits(d, 51, 53) != 0) return false;   // reserved
  return true;
}

bool is44(const std::string& d) {
  if (wrong_status(d, 5, 6, 23) || wrong_status(d, 35, 36, 46) || wrong_status(d, 47, 48, 49) ||
      wrong_status(d, 50, 51, 56))
    return false;
  if (ubits(d, 0, 4) > 4) return false;            // figure of merit / source
  const Json w = field("wind44", d);
  if (!w.is_array() || w[0].get<int64_t>() > 250) return false;   // wind is required
  const double t = field("temp44", d).get<double>();
  if (t > 60 || t < -80) return false;
  return true;
}

bool is45(const std::string& d) {
  if (wrong_status(d, 1, 2, 3) || wrong_status(d, 4, 5, 6) || wrong_status(d, 7, 8, 9) ||
      wrong_status(d, 10, 11, 12) || wrong_status(d, 13, 14, 15) || wrong_status(d, 16, 17, 26) ||
      wrong_status(d, 27, 28, 38) || wrong_status(d, 39, 40, 51))
    return false;
  if (ubits(d, 51, 56) != 0) return false;         // reserved
  const Json t = field("temp45", d);
  if (is_num(t) && num(t) != 0 && (num(t) > 60 || num(t) < -80)) return false;
  return true;
}

bool is50(const std::string& d) {
  if (wrong_status(d, 1, 3, 11) || wrong_status(d, 12, 13, 23) || wrong_status(d, 24, 25, 34) ||
      wrong_status(d, 35, 36, 45) || wrong_status(d, 46, 47, 56))
    return false;
  const Json roll = field("roll50", d), gs = field("gs50", d), tas = field("tas50", d);
  if (is_num(roll) && std::abs(num(roll)) > 50) return false;
  if (is_num(gs) && num(gs) > 600) return false;
  if (is_num(tas) && num(tas) > 600) return false;
  if (is_num(gs) && is_num(tas) && std::abs(num(tas) - num(gs)) > 200) return false;
  return true;
}

bool is60(const std::string& bits, const std::string& d) {
  if (wrong_status(d, 1, 2, 12) || wrong_status(d, 13, 14, 23) || wrong_status(d, 24, 25, 34) ||
      wrong_status(d, 35, 36, 45) || wrong_status(d, 46, 47, 56))
    return false;
  const Json ias = field("ias60", d), mach = field("mach60", d);
  const Json vb = field("vr60baro", d), vi = field("vr60ins", d);
  if (is_num(ias) && num(ias) > 500) return false;
  if (is_num(mach) && num(mach) > 1) return false;
  if (is_num(vb) && std::abs(num(vb)) > 6000) return false;
  if (is_num(vi) && std::abs(num(vi)) > 6000) return false;
  // with the reply's own altitude (DF20): IAS consistent with Mach
  if (is_num(mach) && is_num(ias) && ubits(bits, 0, 5) == 20) {
    const auto alt = altitude_ac13(bits.substr(19, 13));
    if (alt) {
      const double ias_ = tas2cas(mach2tas(num(mach), static_cast<double>(*alt) * kFt), static_cast<double>(*alt) * kFt) / kKts;
      if (std::abs(num(ias) - ias_) > 20) return false;
    }
  }
  return true;
}

}  // namespace

// ------------------------------------------------------------------ bits
std::string hex_bits(const std::string& hex) {
  std::string b;
  b.reserve(hex.size() * 4);
  for (char c : hex) {
    const int v = (c >= '0' && c <= '9') ? c - '0' : (c >= 'A' && c <= 'F') ? c - 'A' + 10 : c - 'a' + 10;
    for (int k = 3; k >= 0; --k) b += ((v >> k) & 1) ? '1' : '0';
  }
  return b;
}

int64_t ubits(const std::string& b, size_t i, size_t j) {
  int64_t v = 0;
  for (size_t k = i; k < j; ++k) v = (v << 1) | (b[k] == '1');
  return v;
}

// ------------------------------------------------------------------ Annex 10 codes
std::optional<int64_t> altitude_ac13(const std::string& b) {
  if (b.find('1') == std::string::npos) return std::nullopt;
  if (b[6] == '1') {   // M: metric units
    return static_cast<int64_t>(static_cast<double>(ubits(b.substr(0, 6) + b.substr(7), 0, 12)) * 3.28084);
  }
  if (b[8] == '1') {   // Q: 25 ft increments
    return ubits(b.substr(0, 6) + b[7] + b.substr(9), 0, 11) * 25 - 1000;
  }
  // Gillham (Gray) code, 100 ft increments: D2 D4 A1 A2 A4 B1 B2 B4 (500 ft) and C1 C2 C4 (100 ft)
  const std::string g500 = std::string() + b[10] + b[12] + b[1] + b[3] + b[5] + b[7] + b[9] + b[11];
  const std::string g100 = std::string() + b[0] + b[2] + b[4];
  const int64_t n500 = gray2int(ubits(g500, 0, 8));
  int64_t n100 = gray2int(ubits(g100, 0, 3));
  if (n100 == 0 || n100 == 5 || n100 == 6) return std::nullopt;
  if (n100 == 7) n100 = 5;
  if (n500 % 2) n100 = 6 - n100;
  return n500 * 500 + n100 * 100 - 1300;
}

std::string squawk_id13(const std::string& b) {
  // a code of all zeros or all ones is rejected as not a binary string (the reference raises)
  if (b.find('1') == std::string::npos || b.find('0') == std::string::npos)
    throw std::runtime_error("Input must be 13 bits binary string");
  auto oct = [&](size_t i4, size_t i2, size_t i1) { return char('0' + (bit(b, i4) << 2 | bit(b, i2) << 1 | bit(b, i1))); };
  // C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4
  return std::string() + oct(5, 3, 1) + oct(11, 9, 7) + oct(4, 2, 0) + oct(12, 10, 8);
}

char callsign_char(int c) {
  if (c >= 1 && c <= 26) return static_cast<char>('A' + c - 1);
  if (c == 32) return '_';
  if (c >= 48 && c <= 57) return static_cast<char>('0' + c - 48);
  return '#';
}

// ------------------------------------------------------------------ extended squitter (DO-260B)
int es_typecode(const std::string& bits) { return static_cast<int>(ubits(bits, 32, 37)); }

std::string es_callsign(const std::string& bits) {
  std::string cs;
  for (size_t i = 40; i < 88; i += 6) {
    const char c = callsign_char(static_cast<int>(ubits(bits, i, i + 6)));
    if (c != '#') cs += c;
  }
  return cs;
}

Json es_altitude(const std::string& bits, int tc) {
  const std::string a = bits.substr(40, 12);
  if (20 <= tc && tc <= 22) return static_cast<int64_t>(static_cast<double>(ubits(a, 0, 12)) * 3.28084);  // GNSS height, m
  const auto alt = altitude_ac13(a.substr(0, 6) + "0" + a.substr(6));
  return alt ? Json(*alt) : Json();
}

std::pair<Json, Json> es_surface_movement(const std::string& bits) {
  const std::string me = bits.substr(32, 56);
  Json trk;
  if (me[12] == '1') trk = static_cast<double>(ubits(me, 13, 20) * 360) / 128;
  const int64_t mov = ubits(me, 5, 12);
  Json spd;
  // DO-260B movement quantisation: (first code, speed at it, step); integer steps give integer speeds
  if (mov == 1) {
    spd = 0.0;
  } else if (mov == 124) {
    spd = 175.0;
  } else if (mov >= 2 && mov < 124) {
    if (mov < 9) spd = 0.125 + static_cast<double>(mov - 2) * 0.125;
    else if (mov < 13) spd = 1 + static_cast<double>(mov - 9) * 0.25;
    else if (mov < 39) spd = 2 + static_cast<double>(mov - 13) * 0.5;
    else if (mov < 94) spd = 15 + (mov - 39) * 1;
    else if (mov < 109) spd = 70 + (mov - 94) * 2;
    else spd = 100 + (mov - 109) * 5;
  }
  return {spd, trk};
}

Velocity es_velocity(const std::string& bits) {
  const std::string me = bits.substr(32, 56);
  const int64_t sub = ubits(me, 5, 8);
  Velocity v;
  if (sub == 1 || sub == 2) {
    int64_t ew = ubits(me, 14, 24), ns = ubits(me, 25, 35);
    v.kind = "GS";
    v.ang_src = "TRUE_NORTH";
    if (ew != 0 && ns != 0) {
      ew -= 1;
      ns -= 1;
      if (sub == 2) {
        ew *= 4;
        ns *= 4;
      }
      const int64_t we = (me[13] == '1' ? -1 : 1) * ew, sn = (me[24] == '1' ? -1 : 1) * ns;
      v.spd = static_cast<int64_t>(std::sqrt(static_cast<double>(sn * sn + we * we)));
      double trk = std::atan2(static_cast<double>(we), static_cast<double>(sn)) * kRadToDeg;
      v.ang = trk >= 0 ? trk : trk + 360;
    }
  } else {
    if (me[13] == '1') v.ang = static_cast<double>(ubits(me, 14, 24)) / 1024 * 360.0;
    int64_t spd = ubits(me, 25, 35);
    if (spd != 0) v.spd = (spd - 1) * (sub == 4 ? 4 : 1);
    v.kind = me[24] == '0' ? "IAS" : "TAS";
    v.ang_src = "MAGNETIC_NORTH";
  }
  v.vr_src = me[35] == '0' ? "GNSS" : "BARO";
  const int64_t vr = ubits(me, 37, 46);
  if (vr != 0) v.vr = (me[36] == '1' ? -1 : 1) * (vr - 1) * 64;
  return v;
}

Json es_gnss_baro_diff(const std::string& bits) {
  const int64_t v = ubits(bits, 81, 88);
  if (v == 0 || v == 127) return Json();
  return (bits[80] == '1' ? -1 : 1) * (v - 1) * 25;
}

std::string es_emergency_squawk(const std::string& bits) { return squawk_id13(bits.substr(43, 13)); }

// ------------------------------------------------------------------ CPR (DO-260B A.1.7)
int cpr_nl(double lat) {
  // NL(0) = 59, NL(+-87) = 2, |lat| > 87 -> 1 (with numpy isclose tolerances at 0 and 87 degrees)
  if (std::abs(lat) <= 1e-8) return 59;
  if (std::abs(std::abs(lat) - 87) <= 1e-8 + 1e-5 * 87) return 2;
  if (lat > 87 || lat < -87) return 1;
  const double nz = 15;
  const double a = 1 - std::cos(kPi / (2 * nz));
  const double c = std::cos(kPi / 180 * std::abs(lat));
  const double b = c * c;
  const double nl = 2 * kPi / std::acos(1 - a / b);
  return static_cast<int>(floor_i(nl));
}

std::optional<std::pair<double, double>> cpr_airborne_global(const std::string& even_hex, const std::string& odd_hex,
                                                             double t0, double t1) {
  std::string m0 = hex_bits(even_hex).substr(32), m1 = hex_bits(odd_hex).substr(32);
  if (m0[21] == '1' && m1[21] == '0') {
    std::swap(m0, m1);
    std::swap(t0, t1);
  } else if (!(m0[21] == '0' && m1[21] == '1')) {
    throw std::runtime_error("RuntimeError: Both even and odd CPR frames are required.");
  }
  const double lat_e = static_cast<double>(ubits(m0, 22, 39)) / 131072, lon_e = static_cast<double>(ubits(m0, 39, 56)) / 131072;
  const double lat_o = static_cast<double>(ubits(m1, 22, 39)) / 131072, lon_o = static_cast<double>(ubits(m1, 39, 56)) / 131072;
  const double dlat_e = 360.0 / 60, dlat_o = 360.0 / 59;
  const int64_t j = floor_i(59 * lat_e - 60 * lat_o + 0.5);
  double rlat_e = dlat_e * (static_cast<double>(py_imod(j, 60)) + lat_e);
  double rlat_o = dlat_o * (static_cast<double>(py_imod(j, 59)) + lat_o);
  if (rlat_e >= 270) rlat_e -= 360;
  if (rlat_o >= 270) rlat_o -= 360;
  if (cpr_nl(rlat_e) != cpr_nl(rlat_o)) return std::nullopt;
  double lat, lon;
  if (t0 > t1) {
    lat = rlat_e;
    const int nl = cpr_nl(lat);
    const int ni = std::max(nl, 1);
    const int64_t m = floor_i(lon_e * (nl - 1) - lon_o * nl + 0.5);
    lon = (360.0 / ni) * (static_cast<double>(py_imod(m, ni)) + lon_e);
  } else {
    lat = rlat_o;
    const int nl = cpr_nl(lat);
    const int ni = std::max(nl - 1, 1);
    const int64_t m = floor_i(lon_e * (nl - 1) - lon_o * nl + 0.5);
    lon = (360.0 / ni) * (static_cast<double>(py_imod(m, ni)) + lon_o);
  }
  if (lon > 180) lon -= 360;
  return std::make_pair(lat, lon);   // not rounded
}

// DO-260B A.1.7.5/7: j = floor(lat_ref / dlat) + floor(1/2 + MOD(lat_ref, dlat) / dlat - YZ / 2^17), with the
// floor of the quotient taken consistently with MOD (floored division)
std::pair<double, double> cpr_local(const std::string& hex, double lat_ref, double lon_ref, bool surface) {
  const std::string mb = hex_bits(hex).substr(32);
  const double cprlat = static_cast<double>(ubits(mb, 22, 39)) / 131072;
  const double cprlon = static_cast<double>(ubits(mb, 39, 56)) / 131072;
  const int i = mb[21] == '1';
  const double span = surface ? 90.0 : 360.0;
  const double d_lat = i ? span / 59 : span / 60;
  const int64_t j = static_cast<int64_t>(py_floordiv(lat_ref, d_lat)) + floor_i(0.5 + (py_fmod(lat_ref, d_lat) / d_lat) - cprlat);
  const double lat = d_lat * (static_cast<double>(j) + cprlat);
  const int ni = cpr_nl(lat) - i;
  const double d_lon = ni > 0 ? span / ni : span;
  const int64_t m = static_cast<int64_t>(py_floordiv(lon_ref, d_lon)) + floor_i(0.5 + (py_fmod(lon_ref, d_lon) / d_lon) - cprlon);
  const double lon = d_lon * (static_cast<double>(m) + cprlon);
  return {lat, lon};   // not rounded (tracks rounds the fix, keeps the full value as reference)
}

// ------------------------------------------------------------------ Comm-B
Json commb_infer(const std::string& bits) {
  const std::string d = mb_of(bits);
  if (all_zero(d)) return "EMPTY";
  std::string out;
  auto add = [&](bool ok, const char* name) {
    if (!ok) return;
    if (!out.empty()) out += ",";
    out += name;
  };
  add(is10(d), "BDS10");
  add(is17(d), "BDS17");
  add(is20(d), "BDS20");
  add(is30(d), "BDS30");
  add(is40(d), "BDS40");
  add(is44(d), "BDS44");
  add(is45(d), "BDS45");
  add(is50(d), "BDS50");
  add(is60(bits, d), "BDS60");
  return out.empty() ? Json() : Json(out);
}

Json commb_field(const std::string& fn, const std::string& bits) { return field(fn, mb_of(bits)); }

Json commb_cap17(const std::string& bits) {
  const std::string d = mb_of(bits);
  Json out = Json::array();
  for (size_t i = 0; i < kCap17.size(); ++i)
    if (d[i] == '1') out.push_back(kCap17[i]);
  return out;
}

Json commb_is50or60(const std::string& bits, double spd_ref, double trk_ref, double alt_ref) {
  const std::string d = mb_of(bits);
  if (!(is50(d) && is60(bits, d))) return Json();
  const Json h50 = field("trk50", d), v50 = field("gs50", d);
  if (h50.is_null() || v50.is_null()) return "BDS50,BDS60";
  const Json h60 = field("hdg60", d), m60j = field("mach60", d), i60j = field("ias60", d);
  // IAS and Mach of a BDS60 reading must agree at the aircraft's altitude (as in is60 for DF20)
  if (!m60j.is_null() && !i60j.is_null()) {
    const double ias_ = tas2cas(mach2tas(num(m60j), alt_ref * kFt), alt_ref * kFt) / kKts;
    if (std::abs(num(i60j) - ias_) > 20) return "BDS50";
  }
  if (h60.is_null() || (m60j.is_null() && i60j.is_null())) return "BDS50,BDS60";
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double m60 = m60j.is_null() ? nan : num(m60j), i60 = i60j.is_null() ? nan : num(i60j);
  auto vxy = [](double v, double ang) {
    return std::make_pair(v * std::sin(ang * kDegToRad), v * std::cos(ang * kDegToRad));
  };
  const auto xy5 = vxy(num(v50) * kKts, num(h50));
  const auto xy6m = vxy(mach2tas(m60, alt_ref * kFt), num(h60));
  const auto xy6i = vxy(cas2tas(i60 * kKts, alt_ref * kFt), num(h60));
  const auto mu = vxy(spd_ref * kKts, trk_ref);
  const std::array<std::pair<double, double>, 3> X = {xy5, xy6m, xy6i};
  const std::array<const char*, 3> names = {"BDS50", "BDS60", "BDS60"};
  int best = -1;
  double bd = 0;
  for (int k = 0; k < 3; ++k) {   // np.nanargmin of the euclidean distances
    const double dx = X[k].first - mu.first, dy = X[k].second - mu.second;
    const double dist = std::sqrt(dx * dx + dy * dy);
    if (std::isnan(dist)) continue;
    if (best < 0 || dist < bd) {
      best = k;
      bd = dist;
    }
  }
  if (best < 0) return "BDS50,BDS60";
  return names[best];
}

}  // namespace si::rx::adsb
