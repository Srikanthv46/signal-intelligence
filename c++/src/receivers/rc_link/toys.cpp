// python receivers/rc_link/toys.py (see there for the protocols, addresses and checksums).
#include "receivers/rc_link/toys.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <optional>
#include <string>

#include "core/npcompat.hpp"
#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::toys {
namespace {
using Bytes = std::vector<uint8_t>;

unsigned sum8(const uint8_t* p, size_t n) {
  unsigned s = 0;
  for (size_t i = 0; i < n; ++i) s += p[i];
  return s & 0xFF;
}
unsigned xor8(const uint8_t* p, size_t n) {
  unsigned v = 0;
  for (size_t i = 0; i < n; ++i) v ^= p[i];
  return v;
}
std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}
Bytes unhex(const char* h) {
  Bytes v;
  for (size_t i = 0; h[i] && h[i + 1]; i += 2) v.push_back(static_cast<uint8_t>(std::stoi(std::string(h + i, 2), nullptr, 16)));
  return v;
}

struct XnRule {
  const char* name;
  std::optional<Bytes> addr;
  size_t n;
  std::function<bool(const uint8_t*)> rule;
  int tid0, tid1;   // -1 = address
};

const std::vector<XnRule>& xn_rules() {
  static const std::vector<XnRule> r = [] {
    auto mjxq = [](const uint8_t* p) { return sum8(p, 15) == p[15]; };
    auto cx10 = [](const uint8_t* p) { return p[0] == 0x55 || p[0] == 0xAA; };
    auto jjrc = [](const uint8_t* p) {
      if (p[0] != 0) return false;
      unsigned s = 2;
      for (int i = 0; i < 13; ++i) s += i == 3 ? 0 : p[i];
      return p[13] == (s & 0xFF) && p[3] == (((((s >> 8) & 0xFF) << 2) + 2) & 0xFF);
    };
    auto h83d = [](const uint8_t* p) {
      return (p[0] == 0x13 || p[0] == 0x14) && sum8(p + 1, 4) == p[8] && sum8(p + 9, 10) == p[19];
    };
    auto e016h = [](const uint8_t* p) { return sum8(p, 9) == p[9]; };
    auto potensic = [](const uint8_t* p) {
      unsigned s = 0;
      for (int i = 1; i < 8; ++i) s += p[i];
      return p[0] == 0x64 && (s & 0x0F) == (p[8] & 0x0Fu);
    };
    return std::vector<XnRule>{
        {"MJXQ (WLH08 / E010 class)", unhex("6d6a777777"), 16, mjxq, 7, 10},
        {"MJXQ (X600 / X800)", unhex("6d6a737373"), 16, mjxq, 7, 10},
        {"JJRC345 / SkyTumbler", unhex("cccccccccc"), 16, jjrc, 14, 16},
        {"CX10 family (CX10 / DM007 / JC3015 / MK33041)", unhex("cccccccccc"), 15, cx10, 1, 5},
        {"H8_3D / CG023 / H20", unhex("c457096521"), 20, h83d, 1, 5},
        {"H8_3D (H20H)", unhex("eeddccbb11"), 20, h83d, 1, 5},
        {"E016H", std::nullopt, 10, e016h, -1, -1},
        {"Potensic A20", std::nullopt, 10, potensic, -1, -1}};
  }();
  return r;
}

struct NrfProto {
  const char* name;
  std::optional<Bytes> addr;
  size_t n;
  std::function<bool(const uint8_t*)> check;
  int tid0, tid1;
};

const std::vector<NrfProto>& nrf_protos() {
  static const std::vector<NrfProto> r = [] {
    auto sum15 = [](const uint8_t* p) { return sum8(p, 15) == p[15]; };
    auto syma = [](const uint8_t* p) { return ((xor8(p, 9) + 0x55) & 0xFF) == p[9]; };
    return std::vector<NrfProto>{{"V2x2 (WLtoys V202 / V222 / V262)", unhex("6688686868"), 16, sum15, 7, 10},
                                 {"SymaX X5C", unhex("6d6a737373"), 16, sum15, -1, -1},
                                 {"SymaX (bind)", unhex("abacadaeaf"), 10, syma, 0, 5},
                                 {"SymaX", std::nullopt, 10, syma, -1, -1}};
  }();
  return r;
}

unsigned bits_to_uint(const std::vector<uint8_t>& b, size_t a, size_t n) {
  unsigned v = 0;
  for (size_t i = 0; i < n; ++i) v = (v << 1) | b[a + i];
  return v;
}

// python _nrf_check
Json nrf_check(const std::vector<uint8_t>& b, size_t s0, const NrfProto& pr, bool inv) {
  uint8_t air[5];
  for (int i = 0; i < 5; ++i) air[i] = static_cast<uint8_t>(bits_to_uint(b, s0 + 8 * i, 8));
  const uint8_t reg[5] = {air[4], air[3], air[2], air[1], air[0]};
  for (size_t pcf : {size_t(0), size_t(9)}) {
    const size_t L = 40 + pcf + 8 * pr.n + 16;
    if (s0 + L > b.size()) continue;
    const std::vector<uint8_t> body(b.begin() + static_cast<int64_t>(s0), b.begin() + static_cast<int64_t>(s0 + L - 16));
    if (crc16_nrf_bits(body) != bits_to_uint(b, s0 + L - 16, 16)) continue;
    Bytes p(pr.n);
    for (size_t i = 0; i < pr.n; ++i) p[i] = static_cast<uint8_t>(bits_to_uint(b, s0 + 40 + pcf + 8 * i, 8));
    if (!pr.check(p.data())) continue;
    Json out = Json::object();
    out["protocol"] = std::string(pr.name) + " (nRF24)";
    out["address"] = hexs(reg, 5);
    out["payload_len"] = pr.n;
    out["payload_hex"] = hexs(p.data(), pr.n);
    out["esb_control_field"] = pcf != 0;
    out["crc_ok"] = true;
    out["checksum_ok"] = true;
    out["polarity_inverted"] = inv;
    out["tx_id"] = pr.tid0 >= 0 ? hexs(p.data() + pr.tid0, static_cast<size_t>(pr.tid1 - pr.tid0)) : hexs(reg, 5);
    return out;
  }
  return nullptr;
}

// python _fullbit_streams
std::vector<std::pair<double, std::vector<uint8_t>>> fullbit_streams(const cfvec& y, double ofs, double rate) {
  std::vector<std::pair<double, std::vector<uint8_t>>> out;
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
    const double thr = sik::two_level_threshold(d);
    std::vector<uint8_t> bits(d.size());
    for (size_t k = 0; k < d.size(); ++k) bits[k] = d[k] > thr;
    out.emplace_back(p0, std::move(bits));
  }
  return out;
}
}  // namespace

uint16_t crc16_nrf_bits(const std::vector<uint8_t>& bits, unsigned init) {
  unsigned c = init;
  for (uint8_t b : bits) c = ((((c >> 15) & 1u) ^ b) != 0) ? ((c << 1) ^ 0x1021u) & 0xFFFF : (c << 1) & 0xFFFF;
  return static_cast<uint16_t>(c);
}

Json match_xn297(const std::vector<uint8_t>& addr, const std::vector<uint8_t>& unscr) {
  for (const auto& r : xn_rules()) {
    if (r.addr && addr != *r.addr) continue;
    if (unscr.size() < r.n || !r.rule(unscr.data())) continue;
    if (!r.addr && std::string(r.name) == "E016H" && !(addr[0] == 0xA5 && addr[1] == 0x00)) continue;
    Json out = Json::object();
    out["protocol"] = std::string(r.name) + " (XN297)";
    out["address"] = hexs(addr.data(), addr.size());
    out["payload_len"] = r.n;
    out["payload_hex"] = hexs(unscr.data(), r.n);
    out["checksum_ok"] = true;
    out["tx_id"] = r.tid0 >= 0 ? hexs(unscr.data() + r.tid0, static_cast<size_t>(r.tid1 - r.tid0))
                               : hexs(addr.data(), addr.size());
    if (std::string(r.name).rfind("CX10", 0) == 0) {
      out["checksum_ok"] = nullptr;
      out["packet_name"] = unscr[0] == 0xAA ? "bind" : "data";
      out["note"] = "no payload checksum: fixed address + header byte only";
    }
    return out;
  }
  return nullptr;
}

Json nrf24_try_decode(const cfvec& y, double ofs, double rate) {
  auto streams = fullbit_streams(y, ofs, rate);
  for (auto& s : sik::demod_bits_mid(y, ofs, rate)) streams.push_back(std::move(s));
  const auto& protos = nrf_protos();
  for (const auto& [ph, bits] : streams)
    for (bool inv : {false, true}) {
      std::vector<uint8_t> b(bits.size());
      for (size_t i = 0; i < b.size(); ++i) b[i] = inv ? 1 - bits[i] : bits[i];
      if (b.size() < 40) continue;
      for (const auto& pr : protos) {
        if (!pr.addr) continue;
        int pat[40];
        for (int i = 0; i < 40; ++i) pat[i] = (((*pr.addr)[4 - i / 8] >> (7 - i % 8)) & 1) ? 1 : -1;
        for (size_t s0 = 0; s0 + 40 <= b.size(); ++s0) {
          int c = 0;
          for (int i = 0; i < 40; ++i) c += (b[s0 + i] ? 1 : -1) * pat[i];
          if (c < 38) continue;
          Json d = nrf_check(b, s0, pr, inv);
          if (!d.is_null()) return d;
        }
      }
      const NrfProto& sy = protos.back();   // SymaX data: A2 + TX ID after the preamble run
      std::vector<int64_t> edges = {-1};
      for (size_t i = 0; i + 1 < b.size(); ++i)
        if (b[i] == b[i + 1]) edges.push_back(static_cast<int64_t>(i));
      edges.push_back(static_cast<int64_t>(b.size()) - 1);
      for (size_t e = 0; e + 1 < edges.size(); ++e) {
        const int64_t a0 = edges[e], a1 = edges[e + 1];
        if (a1 - a0 < 6) continue;
        for (int64_t s0 = std::max(a0 + 1, a1 - 12); s0 < a1 + 3; ++s0) {
          if (s0 < 0 || s0 + 40 > static_cast<int64_t>(b.size())) continue;
          if (bits_to_uint(b, static_cast<size_t>(s0) + 32, 8) != 0xA2) continue;
          Json d = nrf_check(b, static_cast<size_t>(s0), sy, inv);
          if (!d.is_null()) return d;
        }
      }
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::toys
