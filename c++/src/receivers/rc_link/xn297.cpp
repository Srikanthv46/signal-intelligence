// python receivers/rc_link/xn297.py (see there for the sources, the scrambler and the Bayang check).
#include "receivers/rc_link/xn297.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "core/npcompat.hpp"
#include "receivers/rc_link/sik.hpp"
#include "receivers/rc_link/toys.hpp"

namespace si::rx::rc_link::xn297 {
namespace {
constexpr uint32_t kPrefix = 0x710F55;
constexpr int kAddr = 5, kLen = 15;
constexpr int kMinPay = 10, kMaxPay = 21;   // toy protocols (toys.cpp)
bool is_bind(uint8_t t) { return t == 0xA1 || t == 0xA2 || t == 0xA3 || t == 0xA4 || t == 0x53; }
bool is_data(uint8_t t) { return t == 0xA5 || t == 0xA6; }
uint8_t rev(uint8_t v) {
  uint8_t r = 0;
  for (int i = 0; i < 8; ++i) r = static_cast<uint8_t>((r << 1) | ((v >> i) & 1));
  return r;
}
std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}
std::string hx2(unsigned v) {
  char b[8];
  std::snprintf(b, sizeof b, "0x%02X", v);
  return b;
}
struct Hit { int64_t pos; bool inv; int err; };
}  // namespace

const std::vector<uint8_t>& scramble() {
  static const std::vector<uint8_t> t = [] {
    std::vector<uint8_t> out;
    unsigned st = 0x71;
    for (int n = 0; n < 40; ++n) {
      unsigned b = 0;
      for (int i = 0; i < 8; ++i) {
        const unsigned bit = (st >> 6) & 1;
        st = ((st << 1) | (((st >> 6) ^ (st >> 2)) & 1)) & 0x7F;
        b = (b << 1) | bit;
      }
      out.push_back(static_cast<uint8_t>(b));
    }
    return out;
  }();
  return t;
}

uint16_t crc16(const uint8_t* p, size_t n, unsigned init) {
  unsigned c = init;
  for (size_t k = 0; k < n; ++k) {
    c ^= static_cast<unsigned>(p[k]) << 8;
    for (int i = 0; i < 8; ++i) c = (c & 0x8000) ? ((c << 1) ^ 0x1021u) & 0xFFFF : (c << 1) & 0xFFFF;
  }
  return static_cast<uint16_t>(c);
}

Json try_decode(const cfvec& y, double ofs, double rate) {
  const size_t nb = kAddr + kMaxPay + 2, nb_min = kAddr + kMinPay + 2;
  const auto& S = scramble();
  Json like = nullptr;
  char pre[16];
  std::snprintf(pre, sizeof pre, "0x%06X", static_cast<unsigned>(kPrefix));
  int pat[24];
  for (int i = 0; i < 24; ++i) pat[i] = ((kPrefix >> (23 - i)) & 1) ? 1 : -1;
  // python try_decode body for one candidate stream (on-air bits after the prefix): the result to return, or null
  auto handle = [&](const std::vector<uint8_t>& bits, bool inv, int err) -> Json {
    const size_t na = bits.size() / 8;
    std::vector<uint8_t> air(na);
    for (size_t b = 0; b < na; ++b) {
      unsigned v = 0;
      for (int i = 0; i < 8; ++i) v = (v << 1) | bits[8 * b + i];
      air[b] = static_cast<uint8_t>(v);
    }
    uint8_t addr[kAddr];
    for (int j = 0; j < kAddr; ++j) addr[j] = air[kAddr - 1 - j] ^ S[kAddr - 1 - j];
    std::vector<uint8_t> unscr;
    for (size_t i = 0; i < std::min<size_t>(kMaxPay, na - kAddr); ++i)
      unscr.push_back(rev(static_cast<uint8_t>(air[kAddr + i] ^ S[kAddr + i])));
    const uint8_t* p = unscr.data();
    const bool full = na >= static_cast<size_t>(kAddr + kLen + 2);
    char resid[8] = "";
    if (full) {
      const unsigned rx = (air[kAddr + kLen] << 8) | air[kAddr + kLen + 1];
      std::snprintf(resid, sizeof resid, "0x%04X", static_cast<unsigned>(crc16(air.data(), kAddr + kLen) ^ rx));
    }
    unsigned sum = 0;
    for (int i = 0; i < 14 && i < static_cast<int>(unscr.size()); ++i) sum += p[i];
    if (unscr.size() >= static_cast<size_t>(kLen) && full && (is_bind(p[0]) || is_data(p[0])) && p[14] == (sum & 0xFF)) {
      Json out = Json::object();
      out["protocol"] = "Bayang (XN297)";
      out["prefix"] = pre;
      out["prefix_bit_errors"] = err;
      out["polarity_inverted"] = inv;
      out["address"] = hexs(addr, kAddr);
      out["checksum_ok"] = true;
      out["crc_ok"] = nullptr;
      out["crc_residual"] = resid;
      out["payload_hex"] = hexs(p, kLen);
      out["packet_type"] = hx2(p[0]);
      if (is_bind(p[0])) {
        out["packet_name"] = "bind";
        out["tx_id"] = hexs(p + 1, 5);
        out["hop_channels"] = Json::array({p[6], p[7], p[8], p[9]});
      } else {
        out["packet_name"] = "data";
        Json st = Json::object();
        const char* names[] = {"aileron", "elevator", "throttle", "rudder"};
        for (int k = 0; k < 4; ++k) st[names[k]] = ((p[4 + 2 * k] & 0x03) << 8) | p[5 + 2 * k];
        out["sticks"] = st;
        out["flags"] = Json::array({hx2(p[2]), hx2(p[3])});
        out["tx_id_byte2"] = hx2(p[12]);
      }
      return out;
    }
    {  // other toy protocols (research, PLAN G6)
      Json toy = toys::match_xn297(std::vector<uint8_t>(addr, addr + kAddr), unscr);
      if (!toy.is_null()) {
        toy["prefix"] = pre;
        toy["prefix_bit_errors"] = err;
        toy["polarity_inverted"] = inv;
        return toy;
      }
    }
    if (like.is_null() && err == 0) {
      like = Json::object();
      like["protocol"] = "XN297-like (toy drone)";
      like["prefix"] = pre;
      like["prefix_bit_errors"] = 0;
      like["polarity_inverted"] = inv;
      like["address"] = hexs(addr, kAddr);
      like["checksum_ok"] = nullptr;
      like["crc_ok"] = nullptr;
      like["note"] = "5-byte address assumed; payload protocol unknown";
    }
    return nullptr;
  };
  const size_t NB = 8 * nb, NB_MIN = 8 * nb_min;
  // 1. python _soft_streams: full-bit phase change, level from the 24 known prefix bits (PLAN L1c)
  {
    const double sps = ofs / rate;
    dvec ang(y.size());
    for (size_t i = 0; i < y.size(); ++i) ang[i] = std::atan2(static_cast<double>(y[i].imag()), static_cast<double>(y[i].real()));
    const dvec ph = np::unwrap(ang);
    const double step = std::max(sps / 8, 1.0);
    for (double p0 = 0; p0 < sps; p0 += step) {
      const int64_t n = static_cast<int64_t>((static_cast<double>(y.size()) - p0) / sps) - 1;
      if (n < 1) continue;
      dvec d(static_cast<size_t>(n));
      for (int64_t k = 0; k < n; ++k) {
        const int64_t t0 = static_cast<int64_t>(std::nearbyint(p0 + static_cast<double>(k) * sps));
        const int64_t t1 = std::min<int64_t>(static_cast<int64_t>(std::nearbyint(p0 + static_cast<double>(k + 1) * sps)),
                                             static_cast<int64_t>(ph.size()) - 1);
        d[static_cast<size_t>(k)] = ph[static_cast<size_t>(t1)] - ph[static_cast<size_t>(t0)];
      }
      dvec tmp = d;
      const double med = np::median(tmp);
      for (auto& v : d) v -= med;
      if (d.size() < 24 + NB_MIN) continue;
      for (size_t q = 0; q + 24 <= d.size(); ++q) {
        int c = 0;
        for (int i = 0; i < 24; ++i) c += (d[q + i] > 0 ? 1 : d[q + i] < 0 ? -1 : 0) * pat[i];
        if (std::abs(c) < 22) continue;
        const bool inv = c < 0;
        const double sg = inv ? -1.0 : 1.0;
        double sp = 0, sn = 0;
        int np_ = 0, nn = 0;
        for (int i = 0; i < 24; ++i) {
          const double v = sg * d[q + i];
          if (pat[i] > 0) sp += v, ++np_; else sn += v, ++nn;
        }
        const double thr = 0.5 * (sp / np_ + sn / nn);
        const size_t len = std::min(NB, d.size() - (q + 24));
        if (len < NB_MIN) continue;
        std::vector<uint8_t> bits(len);
        for (size_t i = 0; i < len; ++i) bits[i] = (sg * d[q + 24 + i] - thr) > 0 ? 1 : 0;
        Json r = handle(bits, inv, (24 - std::abs(c)) / 2);
        if (!r.is_null()) return r;
      }
    }
  }
  // 2. the discriminator slicer (python: sik.demod_bits + _find_prefix)
  for (const auto& [ph, bits] : sik::demod_bits(y, ofs, rate)) {
    const int64_t n = static_cast<int64_t>(bits.size());
    if (n < static_cast<int64_t>(24 + NB_MIN)) continue;
    std::vector<int> c(static_cast<size_t>(n - 23));
    for (int64_t p = 0; p + 24 <= n; ++p) {
      int v = 0;
      for (int i = 0; i < 24; ++i) v += (bits[static_cast<size_t>(p + i)] ? 1 : -1) * pat[i];
      c[static_cast<size_t>(p)] = v;
    }
    std::vector<Hit> hits;   // python: non-inverted hits, then inverted
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p)
      if (c[static_cast<size_t>(p)] >= 22) hits.push_back({p + 24, false, (24 - c[static_cast<size_t>(p)]) / 2});
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p)
      if (-c[static_cast<size_t>(p)] >= 22) hits.push_back({p + 24, true, (24 + c[static_cast<size_t>(p)]) / 2});
    for (const Hit& h : hits) {
      const size_t len = std::min(NB, bits.size() - static_cast<size_t>(h.pos));
      if (len < NB_MIN) continue;
      std::vector<uint8_t> b(len);
      for (size_t i = 0; i < len; ++i) b[i] = bits[static_cast<size_t>(h.pos) + i] ^ (h.inv ? 1 : 0);
      Json r = handle(b, h.inv, h.err);
      if (!r.is_null()) return r;
    }
  }
  return like;
}

}  // namespace si::rx::rc_link::xn297
