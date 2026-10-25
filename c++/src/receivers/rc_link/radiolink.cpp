// python receivers/rc_link/radiolink.py (see there for the sources, the PHY and the packets).
#include "receivers/rc_link/radiolink.hpp"

#include <cmath>
#include <cstdio>
#include <string>

#include "receivers/rc_link/cc2500.hpp"
#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::radiolink {
namespace {
struct Sync { unsigned word; const char* name; };
constexpr Sync kSyncs[] = {{0xABCD, "surface / air"}, {0xBADC, "DumboRC"}, {0xABA5, "RC4G"}};

uint32_t prng(uint32_t r) { return 0xA5E2A705u * r + 0x754DB79Bu; }

void shuffle(std::vector<int>& hf, uint32_t seed) {
  uint32_t r2 = seed, r1 = seed;
  for (int i = 0; i < 31; ++i) r1 = prng(r1);
  for (int i = 0; i < 16; ++i) {
    const unsigned k = (((r2 >> 16) + r2 + (r1 >> 16) + r1) & 0xFFu) & 0x0Fu;
    r2 = prng(r2);
    r1 = prng(r1);
    std::swap(hf[k], hf[static_cast<size_t>(i)]);
  }
}

std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}

std::string h2(unsigned v) {
  char b[8];
  std::snprintf(b, sizeof b, "0x%02X", v);
  return b;
}

unsigned sum(const uint8_t* p, size_t n) {
  unsigned s = 0;
  for (size_t i = 0; i < n; ++i) s += p[i];
  return s & 0xFF;
}

Json parse(unsigned sync, const char* sub, const std::vector<uint8_t>& raw) {
  if (raw.empty()) return nullptr;
  const unsigned L = raw[0];
  if ((L != 0x21 && L != 0x0F) || raw.size() < L + 3u) return nullptr;
  const uint8_t* f = raw.data();
  if (cc2500::crc16_cc2500(f, L + 1) != ((raw[L + 1] << 8) | raw[L + 2])) return nullptr;
  char sw[8];
  std::snprintf(sw, sizeof sw, "0x%04X", sync);
  Json out = Json::object();
  out["sync"] = sw;
  out["length"] = L;
  out["crc_ok"] = true;
  out["packet_hex"] = hexs(f, L + 1);
  if (sync == 0xABA5) {
    if (L != 0x0F) return nullptr;
    Json ch = Json::array();
    for (int k = 0; k < 4; ++k) ch.push_back(f[4 + k] | (((f[8 + k / 2] >> (4 * (k % 2))) & 0x0F) << 8));
    out["protocol"] = "RadioLink RC4G";
    out["sub_protocol"] = "RC4G";
    out["packet_name"] = "data";
    out["tx_id"] = hexs(f + 1, 3);
    out["channels"] = ch;
    out["next_channel"] = f[15];
    out["checksum_ok"] = nullptr;
    return out;
  }
  if (L == 0x21) {
    if (sum(f + 1, 32) != f[33]) return nullptr;
    const unsigned hdr = f[1];
    if ((hdr & ~0x02u) != 0x01 && (hdr & ~0x02u) != 0x21) return nullptr;
    Json ch = Json::array();
    for (int i = 0; i < 16; ++i) {   // 16 x 11 bits, LSB first, from byte 6
      unsigned v = 0;
      for (int b = 0; b < 11; ++b) {
        const int bit = 11 * i + b;
        v |= ((f[6 + bit / 8] >> (bit % 8)) & 1u) << b;
      }
      ch.push_back(v);
    }
    const unsigned pseudo = f[28] | (f[29] << 8);
    const std::string kind = (hdr & 0x20) ? "air" : (sync == 0xBADC ? "DumboRC" : "surface");
    out["protocol"] = "RadioLink (" + kind + ")";
    out["sub_protocol"] = kind;
    out["packet_name"] = "data";
    out["header"] = h2(hdr);
    out["telemetry_request"] = (hdr & 0x02) != 0;
    out["tx_id"] = hexs(f + 2, 4);
    out["channels"] = ch;
    out["pseudo"] = pseudo;
    out["hop_slot"] = pseudo & 0x0F;
    out["replaced_slot"] = f[32];
    out["checksum_ok"] = true;
    auto [hf, rep] = hop_table(f + 2);
    if (f[32] < 16) hf[f[32]] = rep;
    out["hop_channels"] = hf;
    out["expected_channel"] = hf[pseudo & 0x0F];
    return out;
  }
  if (sum(f + 1, 14) != f[15]) return nullptr;
  out["protocol"] = std::string("RadioLink receiver telemetry (") + sub + ")";
  out["sub_protocol"] = sub;
  out["packet_name"] = "telemetry";
  out["telemetry_type"] = f[1];
  out["tx_id"] = hexs(f + 2, 4);
  out["echoed_header"] = h2(f[6]);
  out["rssi_raw"] = f[7];
  out["checksum_ok"] = true;
  return out;
}
}  // namespace

std::vector<uint8_t> pn9(size_t n) {
  std::vector<uint8_t> out;
  unsigned st = 0x1FF;
  for (size_t i = 0; i < n; ++i) {
    out.push_back(static_cast<uint8_t>(st & 0xFF));
    for (int k = 0; k < 8; ++k) st = (st >> 1) | (((st ^ (st >> 5)) & 1u) << 8);
  }
  return out;
}

std::pair<std::vector<int>, int> hop_table(const uint8_t* id) {
  const int inc = 3 * (id[0] & 3);
  std::vector<int> hf(16);
  for (int i = 0; i < 16; ++i) hf[static_cast<size_t>(i)] = 12 * i + inc;
  for (auto [a, b] : {std::make_pair(0, 1), std::make_pair(2, 3)})
    shuffle(hf, (static_cast<uint32_t>(id[a] + 256u * id[b]) * 0xF65EF9F9u + 0x2EDDF6CAu));
  return {hf, 192 + (inc == 9 ? 6 : inc)};
}

Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz, double rate) {
  static const std::vector<uint8_t> w = pn9(64);
  for (const auto& [ph, bits] : sik::demod_bits(y, ofs, rate))
    for (const auto& sy : kSyncs) {
      const uint8_t s4[4] = {static_cast<uint8_t>(sy.word >> 8), static_cast<uint8_t>(sy.word & 0xFF),
                             static_cast<uint8_t>(sy.word >> 8), static_cast<uint8_t>(sy.word & 0xFF)};
      for (const auto& h : cc2500::find_sync(bits, s4, 4)) {
        const size_t avail = std::min(bits.size() - static_cast<size_t>(h.pos), static_cast<size_t>(8 * (kMaxLen + 3)));
        std::vector<uint8_t> raw(avail / 8);
        for (size_t k = 0; k < raw.size(); ++k) {
          unsigned v = 0;
          for (int i = 0; i < 8; ++i) v = (v << 1) | static_cast<unsigned>(bits[static_cast<size_t>(h.pos) + 8 * k + i] ^ (h.inv ? 1 : 0));
          raw[k] = static_cast<uint8_t>(v ^ w[k]);
        }
        Json d = parse(sy.word, sy.name, raw);
        if (d.is_null()) continue;
        d["sync_bit_errors"] = h.err;
        d["polarity_inverted"] = h.inv;
        if (f_hz) {
          const double n = (*f_hz - kF0) / kStep;
          d["channel"] = static_cast<int64_t>(std::nearbyint(n));
          if (d.contains("expected_channel"))
            d["hop_check"] = std::fabs(n - d["expected_channel"].get<double>()) <= 1.5 ? "match" : "mismatch";
        }
        return d;
      }
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::radiolink
