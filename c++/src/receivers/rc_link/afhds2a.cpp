// python receivers/rc_link/afhds2a.py (see there for the sources, the packet and the candidate conventions).
#include "receivers/rc_link/afhds2a.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::afhds2a {
namespace {
constexpr uint32_t kId = 0x5475C52Au;
constexpr int kIdMaxErr = 3;
const char* type_name(uint8_t t) {
  switch (t) {
    case 0x58: return "sticks";
    case 0x56: return "failsafe";
    case 0xAA: return "settings / telemetry";
    case 0xAC: return "telemetry";
    case 0xBB: case 0xBC: return "bind";
    default: return nullptr;
  }
}
constexpr const char* kLayouts[] = {"systematic", "classic"};
constexpr unsigned kCrc[][2] = {{0x1021, 0x0000}, {0x1021, 0xFFFF}, {0x8005, 0x0000}, {0x8005, 0xFFFF}};

unsigned encode7(unsigned nib, int layout) {
  const unsigned d3 = (nib >> 3) & 1, d2 = (nib >> 2) & 1, d1 = (nib >> 1) & 1, d0 = nib & 1;
  const unsigned b[7] = {layout == 0 ? d3 : d3 ^ d2 ^ d0, layout == 0 ? d2 : d3 ^ d1 ^ d0, layout == 0 ? d1 : d3,
                         layout == 0 ? d0 : d2 ^ d1 ^ d0, layout == 0 ? d3 ^ d2 ^ d0 : d2,
                         layout == 0 ? d3 ^ d1 ^ d0 : d1, layout == 0 ? d2 ^ d1 ^ d0 : d0};
  unsigned v = 0;
  for (unsigned x : b) v = (v << 1) | x;
  return v;
}

struct Tables { std::array<std::array<uint8_t, 16>, 2> enc; std::array<std::array<uint8_t, 128>, 2> dec; };
const Tables& tables() {
  static const Tables t = [] {
    Tables a{};
    for (int l = 0; l < 2; ++l) {
      for (unsigned n = 0; n < 16; ++n) a.enc[l][n] = static_cast<uint8_t>(encode7(n, l));
      for (unsigned w = 0; w < 128; ++w) {   // nearest codeword; python min() keeps the first on ties
        int best = 0, bd = 99;
        for (int n = 0; n < 16; ++n) {
          const int d = __builtin_popcount(a.enc[l][static_cast<size_t>(n)] ^ w);
          if (d < bd) bd = d, best = n;
        }
        a.dec[l][w] = static_cast<uint8_t>(best);
      }
    }
    return a;
  }();
  return t;
}

struct Hit { int64_t pos; bool inv; int err; };
std::vector<Hit> find_id(const std::vector<uint8_t>& bits) {
  std::vector<Hit> hits;
  const int64_t n = static_cast<int64_t>(bits.size());
  if (n < 40) return hits;
  std::vector<int> c(static_cast<size_t>(n - 31));
  for (int64_t p = 0; p + 32 <= n; ++p) {
    int v = 0;
    for (int i = 0; i < 32; ++i) v += (bits[static_cast<size_t>(p + i)] ? 1 : -1) * (((kId >> (31 - i)) & 1) ? 1 : -1);
    c[static_cast<size_t>(p)] = v;
  }
  for (bool inv : {false, true})
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p) {
      const int cc = inv ? -c[static_cast<size_t>(p)] : c[static_cast<size_t>(p)];
      if (cc < 32 - 2 * kIdMaxErr || p < 8) continue;
      int tr = 0;
      for (int64_t j = p - 8; j + 1 < p; ++j) tr += bits[static_cast<size_t>(j)] != bits[static_cast<size_t>(j + 1)];
      if (tr >= 6) hits.push_back({p + 32, inv, (32 - cc) / 2});
    }
  return hits;
}

std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}

constexpr int kA1Len = 21;
constexpr uint8_t kA1Table[8][4] = {{0x12, 0x34, 0x56, 0x78}, {0x18, 0x27, 0x36, 0x45}, {0x41, 0x82, 0x36, 0x57},
                                    {0x84, 0x13, 0x65, 0x72}, {0x87, 0x64, 0x15, 0x32}, {0x76, 0x84, 0x13, 0x52},
                                    {0x71, 0x62, 0x84, 0x35}, {0x71, 0x86, 0x43, 0x52}};

// python _afhds1(raw, f_hz)
Json afhds1(const std::vector<uint8_t>& raw, std::optional<double> f_hz) {
  if (raw.size() < kA1Len + 2u || (raw[0] != 0x55 && raw[0] != 0xAA)) return nullptr;
  const uint8_t* p = raw.data();
  Json ch = Json::array();
  for (int i = 0; i < 8; ++i) {
    const int v = p[5 + 2 * i] | (p[6 + 2 * i] << 8);
    if (v < 800 || v > 2200) return nullptr;
    ch.push_back(v);
  }
  const auto hops = afhds1_hops(p[1]);
  char t[8];
  std::snprintf(t, sizeof t, "0x%02X", p[0]);
  const uint8_t tid[4] = {p[4], p[3], p[2], p[1]};
  Json out = Json::object();
  out["protocol"] = "FlySky AFHDS 1 (A7105)";
  out["packet_type"] = t;
  out["packet_name"] = p[0] == 0x55 ? "data" : "bind";
  out["tx_id"] = hexs(tid, 4);
  out["channels_us"] = ch;
  out["hop_channels"] = hops;
  out["packet_hex"] = hexs(p, kA1Len);
  out["crc_ok"] = nullptr;
  for (const auto& cv : kCrc)
    if (crc16(p, kA1Len, cv[0], cv[1]) == ((raw[kA1Len] << 8) | raw[kA1Len + 1])) {
      char cp[8], ci[8];
      std::snprintf(cp, sizeof cp, "0x%04X", cv[0]);
      std::snprintf(ci, sizeof ci, "0x%04X", cv[1]);
      out["crc_ok"] = true;
      out["convention"] = {{"crc_poly", cp}, {"crc_init", ci}};
      break;
    }
  if (f_hz) {
    const int64_t n = static_cast<int64_t>(std::nearbyint((*f_hz - kA7105F0) / kA7105Step));
    out["channel"] = n;
    bool in = p[0] == 0xAA && n == 1;
    for (int h : hops) in = in || h == n;
    out["hop_check"] = in ? "match" : "mismatch";
  }
  const bool crc = out["crc_ok"].is_boolean();
  if (!crc && out.value("hop_check", std::string()) != "match") return nullptr;
  if (!crc) out["protocol"] = "FlySky AFHDS 1-like (A7105; header, channels and hop table verified, CRC convention not)";
  return out;
}
}  // namespace

std::vector<int> afhds1_hops(unsigned id3) {
  const unsigned row = id3 & 0x0F, off = id3 >> 4;
  std::vector<int> hf(16);
  for (int i = 0; i < 16; ++i) {
    const unsigned b = kA1Table[row >> 1][i >> 2];
    const unsigned n = (i & 2) ? b & 0x0F : b >> 4;
    const int t = static_cast<int>(10 * n) + ((i & 1) ? 80 : 0);
    hf[static_cast<size_t>((row & 1) ? 15 - i : i)] = t - static_cast<int>(off);
  }
  return hf;
}

uint16_t crc16(const uint8_t* p, size_t n, unsigned poly, unsigned init) {
  unsigned c = init;
  for (size_t k = 0; k < n; ++k) {
    c ^= static_cast<unsigned>(p[k]) << 8;
    for (int i = 0; i < 8; ++i) c = (c & 0x8000) ? ((c << 1) ^ poly) & 0xFFFF : (c << 1) & 0xFFFF;
  }
  return static_cast<uint16_t>(c);
}

Json try_identify(const cfvec& y, double ofs, double rate, std::optional<double> f_hz) {
  const size_t nbits = (kPayload + 2) * 2 * 7;
  const Tables& T = tables();
  Json best = nullptr;
  int best_err = 1 << 30;
  char idh[16];
  std::snprintf(idh, sizeof idh, "0x%08X", static_cast<unsigned>(kId));
  for (const auto& [ph, bits] : sik::demod_bits(y, ofs, rate))
    for (const Hit& h : find_id(bits)) {
      if (static_cast<size_t>(h.pos) + 8 * (kA1Len + 2) <= bits.size()) {   // FlySky AFHDS 1 (PLAN G5)
        std::vector<uint8_t> raw(kA1Len + 2);
        for (size_t k = 0; k < raw.size(); ++k) {
          unsigned v = 0;
          for (int i = 0; i < 8; ++i) v = (v << 1) | static_cast<unsigned>(bits[static_cast<size_t>(h.pos) + 8 * k + i] ^ (h.inv ? 1 : 0));
          raw[k] = static_cast<uint8_t>(v);
        }
        Json d1 = afhds1(raw, f_hz);
        if (!d1.is_null()) {
          d1["id_code"] = idh;
          d1["id_bit_errors"] = h.err;
          d1["polarity_inverted"] = h.inv;
          return d1;
        }
      }
      if (static_cast<size_t>(h.pos) + nbits <= bits.size()) {
        std::vector<unsigned> w(nbits / 7);
        for (size_t k = 0; k < w.size(); ++k) {
          unsigned v = 0;
          for (int i = 0; i < 7; ++i) v = (v << 1) | static_cast<unsigned>(bits[static_cast<size_t>(h.pos) + 7 * k + i] ^ (h.inv ? 1 : 0));
          w[k] = v;
        }
        for (int l = 0; l < 2; ++l)
          for (bool high_first : {true, false}) {
            std::vector<uint8_t> raw(w.size() / 2);
            int fixes = 0;
            for (size_t k = 0; k < w.size(); ++k) fixes += T.enc[l][T.dec[l][w[k]]] != w[k];
            for (size_t k = 0; k < raw.size(); ++k) {
              const unsigned a = T.dec[l][w[2 * k]], b = T.dec[l][w[2 * k + 1]];
              raw[k] = static_cast<uint8_t>(high_first ? (a << 4) | b : (b << 4) | a);
            }
            if (!type_name(raw[0])) continue;
            for (const auto& cv : kCrc) {
              if (crc16(raw.data(), kPayload, cv[0], cv[1]) != ((raw[kPayload] << 8) | raw[kPayload + 1])) continue;
              char t[8], cp[8], ci[8];
              std::snprintf(t, sizeof t, "0x%02X", raw[0]);
              std::snprintf(cp, sizeof cp, "0x%04X", cv[0]);
              std::snprintf(ci, sizeof ci, "0x%04X", cv[1]);
              Json out = Json::object();
              out["protocol"] = "FlySky AFHDS 2A (A7105)";
              out["id_code"] = idh;
              out["id_bit_errors"] = h.err;
              out["polarity_inverted"] = h.inv;
              out["crc_ok"] = true;
              out["packet_type"] = t;
              out["packet_name"] = type_name(raw[0]);
              out["tx_id"] = hexs(raw.data() + 1, 4);
              out["rx_id"] = hexs(raw.data() + 5, 4);
              out["binding"] = raw[5] == 0xFF && raw[6] == 0xFF && raw[7] == 0xFF && raw[8] == 0xFF;
              out["fec_corrections"] = fixes;
              out["packet_hex"] = hexs(raw.data(), kPayload);
              out["convention"] = {{"fec", kLayouts[l]}, {"nibble", high_first ? "high first" : "low first"},
                                   {"crc_poly", cp}, {"crc_init", ci}};
              return out;
            }
          }
      }
      if (best.is_null() || h.err < best_err) {
        best_err = h.err;
        best = Json::object();
        best["protocol"] = "FlySky AFHDS 2A-like (A7105 ID code)";
        best["id_code"] = idh;
        best["id_bit_errors"] = h.err;
        best["polarity_inverted"] = h.inv;
        best["crc_ok"] = nullptr;
        best["note"] = "ID code only: FEC / CRC convention not verified (research)";
      }
    }
  return best;
}

}  // namespace si::rx::rc_link::afhds2a
