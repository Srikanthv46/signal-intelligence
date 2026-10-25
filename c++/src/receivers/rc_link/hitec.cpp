// python receivers/rc_link/hitec.py (see there for the sources, the PHYs and the packets).
#include "receivers/rc_link/hitec.hpp"

#include <cmath>
#include <cstdio>
#include <string>

#include "core/types.hpp"
#include "receivers/rc_link/cc2500.hpp"
#include "receivers/rc_link/radiolink.hpp"
#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::hitec {
namespace {
constexpr unsigned kBindSync = 0xD391;

std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}

std::string hx(unsigned v, int w) {
  char b[12];
  std::snprintf(b, sizeof b, "0x%0*X", w, v);
  return b;
}

const char* length_kind(unsigned L) {
  switch (L) {
    case 0x1A: return "data";
    case 0x16: return "bind";
    case 0x0C: return "telemetry";
    case 0x0A: return "bind reply";
    default: return nullptr;
  }
}

Json hitec_parse(unsigned sync, const std::vector<uint8_t>& raw) {
  if (raw.empty()) return nullptr;
  const unsigned L = raw[0];
  const char* kind = length_kind(L);
  if (!kind || raw.size() < L + 3u) return nullptr;
  const uint8_t* f = raw.data();
  if (cc2500::crc16_cc2500(f, L + 1) != ((raw[L + 1] << 8) | raw[L + 2])) return nullptr;
  if (sync != kBindSync && sync != ((static_cast<unsigned>(f[3]) << 8) | f[2])) return nullptr;
  const std::string k = kind;
  Json out = Json::object();
  out["protocol"] = "Hitec Optima / Minima";
  out["packet_name"] = k;
  out["tx_id"] = hexs(f + 1, 3);
  out["length"] = L;
  out["crc_ok"] = true;
  out["sync_matches_tx_id"] = sync != kBindSync;
  out["packet_hex"] = hexs(f, L + 1);
  if (k == "data") {
    Json ch = Json::array();
    for (int i = 0; i < 9; ++i) ch.push_back((f[4 + 2 * i] << 8) | f[5 + 2 * i]);
    out["channels"] = ch;
    out["repeat_flag"] = hx(f[23], 2);
    out["marker_ok"] = f[22] == 0xEE;
  } else if (k == "bind") {
    const unsigned ph = f[4];
    out["bind_phase"] = hx(ph, 2);
    out["variant"] = (ph & 0xF0) == 0x80 ? "Minima" : "Optima";
    if (ph == 0x72 || ph == 0x82 || ph == 0x73 || ph == 0x83) {
      Json hp = Json::array();
      const int n = (ph == 0x72 || ph == 0x82) ? 14 : 7;
      for (int i = 0; i < n; ++i) hp.push_back(2 * f[5 + i]);
      out["hop_channels_part"] = hp;
    }
  } else if (k == "telemetry") {
    out["protocol"] = "Hitec receiver telemetry";
    out["frame"] = hx(f[5], 2);
    if (f[5] == 0x00 || f[5] == 0x11) out["rx_battery_v"] = py_round(((f[6] << 8) | f[7]) / 28.0, 2);
  }
  return out;
}

// python _doubled_words(b, min_pre=16)
std::vector<int64_t> doubled_words(const std::vector<uint8_t>& b, int min_pre = 16) {
  std::vector<int64_t> out;
  const int64_t n = static_cast<int64_t>(b.size());
  if (n < min_pre + 32) return out;
  std::vector<int64_t> edges = {-1};
  for (int64_t i = 0; i + 1 < n; ++i)
    if (b[static_cast<size_t>(i)] == b[static_cast<size_t>(i + 1)]) edges.push_back(i);
  edges.push_back(n - 1);
  for (size_t e = 0; e + 1 < edges.size(); ++e) {
    const int64_t a0 = edges[e], a1 = edges[e + 1];
    if (a1 - a0 < min_pre) continue;
    for (int64_t p = std::max(a0 + 1, a1 - 14); p < std::min(a1 + 3, n - 32); ++p) {
      bool same = true, broken = false;
      for (int i = 0; i < 16; ++i) same = same && b[static_cast<size_t>(p + i)] == b[static_cast<size_t>(p + 16 + i)];
      for (int i = 0; i + 1 < 32; ++i) broken = broken || b[static_cast<size_t>(p + i)] == b[static_cast<size_t>(p + i + 1)];
      if (same && broken) out.push_back(p);
    }
  }
  return out;
}

const char* model_hint(uint8_t a, uint8_t b) {
  if (a == 0x42 && b == 0x07) return "T8J";
  if (a == 0x0F && b == 0x09) return "T10J";
  if (a == 0x9A && b == 0x06) return "TM-FH";
  return nullptr;
}

Json sfhss_parse(const uint8_t* p, unsigned crc) {
  int init = -1;
  for (unsigned i : {0xFFFFu, 0x0000u})
    if (init < 0 && crc16_cc2400(p, 13, static_cast<uint16_t>(i)) == crc) init = static_cast<int>(i);
  if (init < 0 || (p[0] != 0x80 && p[0] != 0x81 && p[0] != 0x82)) return nullptr;
  const int k = p[5] >> 3;
  if (k > 29) return nullptr;
  Json ch = Json::array({((p[5] & 7) << 9) | (p[6] << 1) | (p[7] >> 7), ((p[7] & 0x7F) << 5) | (p[8] >> 3),
                         ((p[8] & 7) << 9) | (p[9] << 1) | (p[10] >> 7), ((p[10] & 0x7F) << 5) | (p[11] >> 3)});
  const unsigned code = ((p[11] & 7u) << 2) | (p[12] >> 6), cmd = p[12] & 0x3Fu;
  const char* hint = model_hint(p[3], p[4]);
  Json out = Json::object();
  out["protocol"] = "Futaba S-FHSS";
  out["packet_name"] = "data";
  out["header"] = hx(p[0], 2);
  out["tx_id"] = hexs(p + 1, 2);
  out["model_bytes"] = hexs(p + 3, 2);
  out["model_hint"] = hint ? Json(hint) : Json(nullptr);
  out["hop_index"] = k;
  out["hop_code"] = code;
  out["channels"] = ch;
  out["command"] = hx(cmd, 2);
  out["data_half"] = (cmd & 1) ? 2 : 1;
  out["failsafe"] = (cmd & 4) != 0;
  out["channel_bank"] = (cmd & 8) ? "CH5-8" : "CH1-4";
  out["crc_ok"] = true;
  out["crc_init"] = hx(static_cast<unsigned>(init), 4);
  out["packet_hex"] = hexs(p, 13);
  return out;
}
}  // namespace

uint16_t crc16_cc2400(const uint8_t* p, size_t n, uint16_t init) {
  unsigned c = init;
  for (size_t k = 0; k < n; ++k)
    for (int i = 0; i < 8; ++i) {
      const unsigned bit = (p[k] >> (7 - i)) & 1u;
      c = (((c >> 15) & 1u) ^ bit) ? ((c << 1) ^ 0x8005u) & 0xFFFFu : (c << 1) & 0xFFFFu;
    }
  return static_cast<uint16_t>(c);
}

Json hitec_try_decode(const cfvec& y, double ofs, std::optional<double> f_hz, double rate) {
  static const std::vector<uint8_t> w = radiolink::pn9(64);
  for (const auto& [ph, bits] : sik::demod_bits_mid(y, ofs, rate))
    for (int inv = 0; inv < 2; ++inv) {
      std::vector<uint8_t> b = bits;
      if (inv)
        for (auto& v : b) v ^= 1;
      for (int64_t pos : doubled_words(b)) {
        unsigned word = 0;
        for (int i = 0; i < 16; ++i) word = (word << 1) | b[static_cast<size_t>(pos + i)];
        const size_t start = static_cast<size_t>(pos + 32);
        const size_t nbits = std::min(b.size() - start, static_cast<size_t>(8 * (0x1A + 3)));
        std::vector<uint8_t> raw(nbits / 8);
        for (size_t k = 0; k < raw.size(); ++k) {
          unsigned v = 0;
          for (int i = 0; i < 8; ++i) v = (v << 1) | b[start + 8 * k + static_cast<size_t>(i)];
          raw[k] = static_cast<uint8_t>(v ^ w[k]);
        }
        Json d = hitec_parse(word, raw);
        if (d.is_null()) continue;
        d["sync"] = hx(word, 4);
        d["polarity_inverted"] = inv != 0;
        if (f_hz) d["channr"] = static_cast<int64_t>(std::nearbyint((*f_hz - kF0) / kStep));
        return d;
      }
    }
  return nullptr;
}

Json sfhss_try_decode(const cfvec& y, double ofs, std::optional<double> f_hz, double rate) {
  const uint8_t s4[4] = {0xD3, 0x91, 0xD3, 0x91};
  for (const auto& [ph, bits] : sik::demod_bits_mid(y, ofs, rate))
    for (int ns : {4, 2})
      for (const auto& h : cc2500::find_sync(bits, s4, ns)) {
        constexpr size_t kBits = 8 * (13 + 2);
        if (bits.size() - static_cast<size_t>(h.pos) < kBits) continue;
        uint8_t raw[15];
        for (size_t k = 0; k < 15; ++k) {
          unsigned v = 0;
          for (int i = 0; i < 8; ++i) v = (v << 1) | static_cast<unsigned>(bits[static_cast<size_t>(h.pos) + 8 * k + i] ^ (h.inv ? 1 : 0));
          raw[k] = static_cast<uint8_t>(v);
        }
        Json d = sfhss_parse(raw, (static_cast<unsigned>(raw[13]) << 8) | raw[14]);
        if (d.is_null()) continue;
        d["sync_bit_errors"] = h.err;
        d["polarity_inverted"] = h.inv;
        d["sync_bits"] = 8 * ns;
        if (f_hz) {
          const double kk = (*f_hz - kSfF0) / kSfStep;
          d["channel"] = static_cast<int64_t>(std::nearbyint(kk));
          d["hop_check"] = std::fabs(kk - d["hop_index"].get<double>()) <= 0.3 ? "match" : "mismatch";
        }
        return d;
      }
  return nullptr;
}

}  // namespace si::rx::rc_link::hitec
