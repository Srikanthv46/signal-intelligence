// python receivers/rc_link/adsl.py and ogn.py (see there for the sources, the PHY and the packets).
#include "receivers/rc_link/adsl.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <tuple>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link {
namespace {
constexpr uint32_t kDelta = 0x9E3779B9u;

std::string fmt6(unsigned v) {
  char b[12];
  std::snprintf(b, sizeof b, "%06X", v);
  return b;
}

// python FL._decode(chips[odd:]): Manchester pairs -> data bits (01 -> 1)
std::vector<uint8_t> man_bits(const std::vector<uint8_t>& ch, int odd) {
  const size_t np_ = ch.size() > static_cast<size_t>(odd) ? (ch.size() - odd) / 2 : 0;
  std::vector<uint8_t> bits(np_);
  for (size_t k = 0; k < np_; ++k) bits[k] = ch[odd + 2 * k + 1] > ch[odd + 2 * k] ? 1 : 0;
  return bits;
}

int64_t s24(int64_t v) { return (v & 0x800000) ? v - (int64_t(1) << 24) : v; }
}  // namespace

namespace adsl {
namespace {
const char* kCategory[] = {"unknown", "light fixed-wing", "small-heavy fixed-wing", "light rotorcraft", "glider",
                           "lighter-than-air", "ultralight / motorised hang-glider", "paraglider", "parachutist",
                           "eVTOL / UAM", "gyrocopter", "UAS open category", "UAS specific category",
                           "UAS certified category", "model plane", "heavy rotorcraft", "hang-glider", "paramotor"};
const char* kFlight[] = {"undefined", "ground", "airborne", "reserved"};
const char* kEmergency[] = {"undefined", "none", "general", "medical", "no communications", "unlawful interference",
                            "downed", "reserved"};

uint32_t mx(uint32_t z, uint32_t y, uint32_t s) { return (((z >> 5) ^ (y << 2)) + ((y >> 3) ^ (z << 4))) ^ ((s ^ y) + z); }

std::vector<uint32_t> words(const std::vector<uint8_t>& d) {
  std::vector<uint32_t> v(d.size() / 4);
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = d[4 * i] | (d[4 * i + 1] << 8) | (d[4 * i + 2] << 16) | (static_cast<uint32_t>(d[4 * i + 3]) << 24);
  return v;
}

std::vector<uint8_t> bytes_of(const std::vector<uint32_t>& v) {
  std::vector<uint8_t> d;
  for (uint32_t w : v)
    for (int k = 0; k < 4; ++k) d.push_back(static_cast<uint8_t>(w >> (8 * k)));
  return d;
}

int64_t exp_u(int64_t raw, int n) {
  const int64_t e = raw >> n, base = raw & ((int64_t(1) << n) - 1);
  return (int64_t(1) << e) * ((int64_t(1) << n) + base) - (int64_t(1) << n);
}

int64_t exp_s(int64_t raw, int n) {
  const int64_t sign = raw >> (n + 2);
  const int64_t v = exp_u(raw & ((int64_t(1) << (n + 2)) - 1), n);
  return sign ? -v : v;
}

void traffic(const uint8_t* pl, Json& out) {
  auto f = [pl](int off, int n) {
    int64_t v = 0;
    for (int i = 0; i < n; ++i) v |= static_cast<int64_t>((pl[(off + i) / 8] >> ((off + i) % 8)) & 1) << i;
    return v;
  };
  const int64_t lat = f(16, 24), lon = f(40, 24), spd = f(64, 8), alt = f(72, 14), vr = f(86, 9), cat = f(8, 5);
  out["timestamp_qs"] = f(0, 6);
  out["flight_state"] = kFlight[f(6, 2)];
  out["aircraft_category"] = cat < 18 ? Json(kCategory[cat]) : Json("reserved " + std::to_string(cat));
  out["uas"] = cat == 11 || cat == 12 || cat == 13;
  out["emergency"] = kEmergency[f(13, 3)];
  out["ground_track_deg"] = py_round(static_cast<double>(f(95, 9)) * 360 / 512, 1);
  out["sil"] = f(104, 2);
  out["design_assurance"] = f(106, 2);
  out["nic"] = f(108, 4);
  out["hfom"] = f(112, 3);
  out["vfom"] = f(115, 2);
  out["velocity_accuracy"] = f(117, 2);
  if (lat != 0x800000) {
    out["lat"] = py_round(static_cast<double>(s24(lat)) / 93206, 5);
    out["lon"] = py_round(static_cast<double>(s24(lon)) / 46603, 5);
  }
  out["ground_speed_ms"] = spd == 0xFF ? Json(nullptr) : Json(static_cast<double>(exp_u(spd, 6)) / 4);
  out["alt_m"] = alt == 0x3FFF ? Json(nullptr) : Json(exp_u(alt, 12) - 320);
  out["climb_ms"] = vr == 0x1FF ? Json(nullptr) : Json(static_cast<double>(exp_s(vr, 6)) / 8);
}
}  // namespace

uint32_t crc24(const uint8_t* p, size_t n) {
  uint32_t c = 0;
  for (size_t k = 0; k < n; ++k) {
    c ^= static_cast<uint32_t>(p[k]) << 16;
    for (int i = 0; i < 8; ++i) c = (c & 0x800000u) ? ((c << 1) ^ 0xFFF409u) & 0xFFFFFFu : (c << 1) & 0xFFFFFFu;
  }
  return c;
}

std::vector<uint8_t> xxtea_encrypt(std::vector<uint8_t> d, int rounds) {
  auto v = words(d);
  const size_t n = v.size();
  uint32_t s = 0, z = v[n - 1];
  for (int r = 0; r < rounds; ++r) {
    s += kDelta;
    for (size_t p = 0; p < n; ++p) {
      const uint32_t y = v[(p + 1) % n];
      v[p] += mx(z, y, s);
      z = v[p];
    }
  }
  return bytes_of(v);
}

std::vector<uint8_t> xxtea_decrypt(std::vector<uint8_t> d, int rounds) {
  auto v = words(d);
  const size_t n = v.size();
  uint32_t s = static_cast<uint32_t>(rounds) * kDelta, y = v[0];
  for (int r = 0; r < rounds; ++r) {
    for (size_t p = n; p-- > 0;) {
      const uint32_t z = v[(p + n - 1) % n];
      v[p] -= mx(z, y, s);
      y = v[p];
    }
    s -= kDelta;
  }
  return bytes_of(v);
}

Json parse(const std::vector<uint8_t>& pkt) {
  if (pkt.empty()) return nullptr;
  const unsigned L = pkt[0];
  if (L < 8 || pkt.size() < L + 1u || (L - 4) % 4) return nullptr;
  const uint8_t* body = pkt.data() + 1;
  if (crc24(body, L - 3) != ((static_cast<uint32_t>(body[L - 3]) << 16) | (body[L - 2] << 8) | body[L - 1])) return nullptr;
  const unsigned hdr = body[0], key = (hdr >> 5) & 3;
  std::vector<uint8_t> data(body + 1, body + L - 3);
  if (key == 0)
    data = xxtea_decrypt(data);
  else if (key != 3)
    return nullptr;
  const unsigned typ = data[0];
  const uint32_t w = data[1] | (data[2] << 8) | (data[3] << 16) | (static_cast<uint32_t>(data[4]) << 24);
  const unsigned amt = w & 0x3F;
  Json out = Json::object();
  out["protocol"] = "ADS-L (SRD-860 M-band)";
  out["crc_ok"] = true;
  out["length"] = L;
  out["version"] = hdr & 0x0F;
  out["signed"] = (hdr & 0x10) != 0;
  out["scrambled"] = key == 0;
  out["error_control"] = (hdr & 0x80) ? "RS FEC" : "CRC-24";
  const char* types[] = {nullptr, nullptr, "traffic", "status", "traffic uplink", "FIS-B uplink", "remote identification"};
  char tb[8];
  std::snprintf(tb, sizeof tb, "0x%02X", typ);
  out["payload_type"] = typ < 7 && types[typ] ? std::string(types[typ]) : (typ == 66 ? std::string("OGN diagnostics") : std::string(tb));
  std::string at;
  switch (amt) {
    case 0: at = "random / privacy"; break;
    case 5: at = "ICAO"; break;
    case 6: at = "FLARM"; break;
    case 7: at = "OGN tracker"; break;
    case 8: at = "FANET"; break;
    default: at = (amt >= 9 && amt <= 63) ? "manufacturer page " + std::to_string(amt - 9) : "reserved " + std::to_string(amt);
  }
  out["address_table"] = at;
  out["address"] = fmt6((w >> 6) & 0xFFFFFFu);
  out["relay"] = (w >> 31) != 0;
  if (typ == 2 && data.size() >= 20) traffic(data.data() + 5, out);
  return out;
}

Json try_decode(const cfvec& y, double ofs, double chip_rate) {
  constexpr uint32_t kSync24 = 0xF5724B;
  for (const auto& [ph, chips0] : sik::demod_bits(y, ofs, chip_rate))
    for (bool inv : {false, true}) {
      std::vector<uint8_t> ch(chips0.size());
      for (size_t i = 0; i < ch.size(); ++i) ch[i] = static_cast<uint8_t>(chips0[i] ^ (inv ? 1 : 0));
      for (int odd : {0, 1}) {
        const auto bits = man_bits(ch, odd);
        if (bits.size() < 24 + 8 * 12) continue;
        for (size_t p = 0; p + 24 <= bits.size(); ++p) {
          int c = 0;
          for (int i = 0; i < 24; ++i) c += (bits[p + i] ? 1 : -1) * (((kSync24 >> (23 - i)) & 1) ? 1 : -1);
          if (c < 22) continue;
          const size_t s = p + 24;
          if (bits.size() - s < 8) continue;
          unsigned L = 0;
          for (int i = 0; i < 8; ++i) L = (L << 1) | bits[s + i];
          if (bits.size() - s < 8 * (L + 1)) continue;
          std::vector<uint8_t> pkt(L + 1);
          for (size_t k = 0; k <= L; ++k) {
            unsigned v = 0;
            for (int i = 0; i < 8; ++i) v = (v << 1) | bits[s + 8 * k + i];
            pkt[k] = static_cast<uint8_t>(v);
          }
          Json d = parse(pkt);
          if (d.is_null()) continue;
          d["sync_bit_errors"] = (24 - c) / 2;
          d["polarity_inverted"] = inv;
          return d;
        }
      }
    }
  return nullptr;
}

}  // namespace adsl

namespace ogn {
namespace {
const char* kAddrTypes[] = {"random", "ICAO", "FLARM", "OGN"};
const char* kAircraft[] = {"unknown", "glider", "tow plane", "helicopter", "skydiver", "drop plane", "hang-glider",
                           "paraglider", "powered aircraft", "jet", "UFO", "balloon", "airship", "UAV",
                           "ground vehicle", "static object"};
int64_t uvr(int64_t v, int n) {
  const int64_t r = v >> n, m = v & ((int64_t(1) << n) - 1), N = int64_t(1) << n;
  switch (r) {
    case 0: return m;
    case 1: return N + 1 + 2 * m;
    case 2: return 3 * N + 2 + 4 * m;
    default: return 7 * N + 4 + 8 * m;
  }
}
int popcount(uint32_t v) {
  int c = 0;
  for (; v; v &= v - 1) ++c;
  return c;
}
}  // namespace

std::pair<uint32_t, uint32_t> tea_key0(uint32_t v0, uint32_t v1, int loops, bool decrypt) {
  if (!decrypt) {
    uint32_t s = 0;
    for (int i = 0; i < loops; ++i) {
      s += kDelta;
      v0 += (v1 << 4) ^ (v1 + s) ^ (v1 >> 5);
      v1 += (v0 << 4) ^ (v0 + s) ^ (v0 >> 5);
    }
    return {v0, v1};
  }
  uint32_t s = kDelta * static_cast<uint32_t>(loops);
  for (int i = 0; i < loops; ++i) {
    v1 -= (v0 << 4) ^ (v0 + s) ^ (v0 >> 5);
    v0 -= (v1 << 4) ^ (v1 + s) ^ (v1 >> 5);
    s -= kDelta;
  }
  return {v0, v1};
}

Json parse(const uint8_t* d) {
  uint32_t w[5];
  for (int i = 0; i < 5; ++i)
    w[i] = d[4 * i] | (d[4 * i + 1] << 8) | (d[4 * i + 2] << 16) | (static_cast<uint32_t>(d[4 * i + 3]) << 24);
  const uint32_t h = w[0];
  if (popcount(h & 0x0FFFFFFFu) % 2) return nullptr;
  Json out = Json::object();
  out["protocol"] = "OGN tracker-like (OGNTP)";
  out["address"] = fmt6(h & 0xFFFFFFu);
  out["address_type"] = kAddrTypes[(h >> 24) & 3];
  out["position_packet"] = ((h >> 26) & 1) == 0;
  out["relay_count"] = (h >> 28) & 3;
  out["encrypted"] = ((h >> 30) & 1) != 0;
  out["emergency"] = (h >> 31) != 0;
  out["header_parity_ok"] = true;
  out["fec"] = "LDPC not checked (matrix not available)";
  if (out["encrypted"].get<bool>() || !out["position_packet"].get<bool>()) return out;
  std::tie(w[1], w[2]) = tea_key0(w[1], w[2], 8, true);
  std::tie(w[3], w[4]) = tea_key0(w[3], w[4], 8, true);
  const double lat = static_cast<double>(s24(w[1] & 0xFFFFFF) * 8 + 4) / 600000;
  const double lon = static_cast<double>(s24(w[2] & 0xFFFFFF) * 16 + 8) / 600000;
  const unsigned t = (w[1] >> 24) & 0x3F, fixq = (w[1] >> 30) & 3;
  if (std::fabs(lat) > 90 || std::fabs(lon) > 180 || t >= 60 || fixq > 2) return nullptr;
  const unsigned at = (w[4] >> 20) & 0xF;
  out["lat"] = py_round(lat, 5);
  out["lon"] = py_round(lon, 5);
  out["time_s"] = t;
  out["fix_quality"] = fixq;
  out["alt_m"] = uvr(w[3] & 0x3FFF, 12);
  out["speed_ms"] = static_cast<double>(uvr((w[3] >> 14) & 0x3FF, 8)) / 10;
  out["heading_deg"] = py_round(static_cast<double>(((w[4] & 0x3FF) * 3600 + 512) >> 10) / 10, 1);
  out["aircraft_type"] = kAircraft[at];
  out["uav"] = at == 13;
  out["stealth"] = ((w[4] >> 19) & 1) != 0;
  return out;
}

Json try_decode(const cfvec& y, double ofs, double chip_rate) {
  constexpr uint32_t kSync = 0x0AF3656Cu;
  for (const auto& [ph, chips0] : sik::demod_bits(y, ofs, chip_rate))
    for (bool inv : {false, true}) {
      std::vector<uint8_t> ch(chips0.size());
      for (size_t i = 0; i < ch.size(); ++i) ch[i] = static_cast<uint8_t>(chips0[i] ^ (inv ? 1 : 0));
      for (int odd : {0, 1}) {
        const auto bits = man_bits(ch, odd);
        if (bits.size() < 32 + 208) continue;
        for (size_t p = 0; p + 32 <= bits.size(); ++p) {
          int c = 0;
          for (int i = 0; i < 32; ++i) c += (bits[p + i] ? 1 : -1) * (((kSync >> (31 - i)) & 1) ? 1 : -1);
          if (c < 30) continue;
          const size_t a = odd + 2 * (p + 32);
          if (a + 416 > ch.size()) continue;
          bool clean = true;
          for (size_t k = 0; k < 208 && clean; ++k) clean = ch[a + 2 * k] != ch[a + 2 * k + 1];
          if (!clean) continue;
          uint8_t d20[20];
          for (size_t k = 0; k < 20; ++k) {
            unsigned v = 0;
            for (int i = 0; i < 8; ++i) v = (v << 1) | bits[p + 32 + 8 * k + i];
            d20[k] = static_cast<uint8_t>(v);
          }
          Json d = parse(d20);
          if (d.is_null()) continue;
          d["sync_bit_errors"] = (32 - c) / 2;
          d["polarity_inverted"] = inv;
          return d;
        }
      }
    }
  return nullptr;
}

}  // namespace ogn
}  // namespace si::rx::rc_link
