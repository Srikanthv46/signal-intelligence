// python receivers/rc_link/ble_coded.py (see the header there for the packet and the method).
#include "receivers/rc_link/ble_coded.hpp"

#include <algorithm>
#include <array>
#include <cstdio>

#include "receivers/rc_link/gfsk_coherent.hpp"

namespace si::rx::rc_link::coded {
namespace {

constexpr int kPre[8] = {0, 0, 1, 1, 1, 1, 0, 0};

struct Trans {
  int pre[8][2], a0[8][2], a1[8][2];
  Trans() {
    for (int ns = 0; ns < 8; ++ns) {
      const int b = ns >> 2, s0 = (ns >> 1) & 1, s1 = ns & 1;
      for (int s2 = 0; s2 < 2; ++s2) {
        pre[ns][s2] = (s0 << 2) | (s1 << 1) | s2;
        a0[ns][s2] = 2 * (b ^ s0 ^ s1 ^ s2) - 1;
        a1[ns][s2] = 2 * (b ^ s1 ^ s2) - 1;
      }
    }
  }
};
const Trans kT;

std::string hex2(uint8_t v) {
  char b[3];
  std::snprintf(b, sizeof b, "%02X", v);
  return b;
}
std::string addr(const std::vector<uint8_t>& by, size_t i) {
  std::string s;
  for (int j = 5; j >= 0; --j) s += (j < 5 ? ":" : "") + hex2(by[i + static_cast<size_t>(j)]);
  return s;
}

}  // namespace

std::vector<uint8_t> conv_encode(const std::vector<uint8_t>& bits) {
  int s0 = 0, s1 = 0, s2 = 0;
  std::vector<uint8_t> out;
  for (uint8_t b : bits) {
    out.push_back(static_cast<uint8_t>(b ^ s0 ^ s1 ^ s2));
    out.push_back(static_cast<uint8_t>(b ^ s1 ^ s2));
    s2 = s1, s1 = s0, s0 = b;
  }
  return out;
}

std::vector<uint8_t> known_symbols(const std::vector<uint8_t>& aa_bits, int ci) {
  std::vector<uint8_t> out;
  for (int r = 0; r < 10; ++r)
    for (int v : kPre) out.push_back(static_cast<uint8_t>(v));
  std::vector<uint8_t> b(aa_bits);
  b.push_back(static_cast<uint8_t>(ci & 1));
  b.push_back(static_cast<uint8_t>((ci >> 1) & 1));
  for (int i = 0; i < 3; ++i) b.push_back(0);
  for (uint8_t c : conv_encode(b))
    for (int v : (c ? std::initializer_list<int>{1, 1, 0, 0} : std::initializer_list<int>{0, 0, 1, 1}))
      out.push_back(static_cast<uint8_t>(v));
  return out;
}

std::vector<uint8_t> viterbi(const dvec& llr) {
  const size_t n = llr.size() / 2;
  double pm[8];
  std::fill(pm, pm + 8, 1e18);
  pm[0] = 0.0;
  std::vector<std::array<int8_t, 8>> surv(n);
  for (size_t t = 0; t < n; ++t) {
    double npm[8];
    for (int ns = 0; ns < 8; ++ns) {
      double m[2];
      for (int k = 0; k < 2; ++k)
        m[k] = pm[kT.pre[ns][k]] - kT.a0[ns][k] * llr[2 * t] - kT.a1[ns][k] * llr[2 * t + 1];
      const int k = m[1] < m[0] ? 1 : 0;
      surv[t][ns] = static_cast<int8_t>(k);
      npm[ns] = m[k];
    }
    std::copy(npm, npm + 8, pm);
  }
  std::vector<uint8_t> out(n);
  int st = 0;
  for (size_t t = n; t-- > 0;) {
    out[t] = static_cast<uint8_t>(st >> 2);
    st = kT.pre[st][surv[t][st]];
  }
  return out;
}

std::vector<std::pair<int64_t, double>> find_coded(const cdvec& p, const std::vector<uint8_t>& aa_bits) {
  return gfsk::find_access_address(p, known_symbols(aa_bits, 0), 2, 0.4);
}

std::optional<Decoded> decode(const cfvec& y, int64_t s0, const std::vector<uint8_t>& aa_bits, int max_pdu_octets) {
  const int64_t n1 = 80 + 4 * 2 * 37;
  const int64_t nsym = n1 + 4 * 2 * (8 * static_cast<int64_t>(max_pdu_octets) + 24 + 3);
  const auto dm = gfsk::demodulate(y, s0, nsym);
  if (!dm || static_cast<int64_t>(dm->bits.size()) < n1 + 2 * 43) return std::nullopt;
  const size_t n = dm->bits.size();
  dvec soft(n);
  for (size_t k = 0; k < n; ++k)
    soft[k] = (2.0 * dm->bits[k] - 1) * std::min(dm->rel[k], dm->rel[k ? k - 1 : 0]);
  const double p8[4] = {1.0, 1.0, -1.0, -1.0};
  auto demap8 = [&](size_t from, size_t ncoded) {
    dvec l(ncoded);
    for (size_t c = 0; c < ncoded; ++c) {
      double v = 0;
      for (int j = 0; j < 4; ++j) v += soft[from + 4 * c + static_cast<size_t>(j)] * p8[j];
      l[c] = v;
    }
    return l;
  };
  const auto blk1 = viterbi(demap8(80, 74));
  Decoded out;
  for (int i = 0; i < 32; ++i) out.aa_bit_errors += blk1[i] != aa_bits[i];
  out.ci = blk1[32] | (blk1[33] << 1);
  out.h = dm->h;
  const size_t rest = n - static_cast<size_t>(n1);
  if (out.ci == 0) {
    const size_t nb = rest / 8;
    out.bits = viterbi(demap8(static_cast<size_t>(n1), 2 * nb));
    out.has_bits = true;
  } else if (out.ci == 1) {
    const size_t nb = rest / 2;
    out.bits = viterbi(dvec(soft.begin() + n1, soft.begin() + n1 + static_cast<ptrdiff_t>(2 * nb)));
    out.has_bits = true;
  }
  return out;
}

Json parse_extended(const std::vector<uint8_t>& by, std::vector<uint8_t>& data) {
  Json out = Json::object();
  data.clear();
  if (by.empty()) return out;
  const size_t hl = by[0] & 0x3F;
  out["adv_mode"] = by[0] >> 6;
  if (hl == 0) {
    data.assign(by.begin() + 1, by.end());
    return out;
  }
  const int flags = by.size() > 1 ? by[1] : 0;
  size_t i = 2;
  const size_t end = 1 + hl;
  if ((flags & 0x01) && i + 6 <= end && i + 6 <= by.size()) out["adv_addr"] = addr(by, i), i += 6;
  if ((flags & 0x02) && i + 6 <= end && i + 6 <= by.size()) out["target_addr"] = addr(by, i), i += 6;
  if ((flags & 0x04) && i + 1 <= end) i += 1;
  if ((flags & 0x08) && i + 2 <= end && i + 2 <= by.size()) {
    const int v = by[i] | (by[i + 1] << 8);
    out["adi"] = {{"did", v & 0xFFF}, {"sid", v >> 12}};
    i += 2;
  }
  if ((flags & 0x10) && i + 3 <= end && i + 3 <= by.size()) {
    const int64_t v = by[i] | (by[i + 1] << 8) | (static_cast<int64_t>(by[i + 2]) << 16);
    const int units = (v >> 7) & 1;
    const int phy = static_cast<int>((v >> 21) & 7);
    out["aux_ptr"] = {{"channel", v & 0x3F},
                      {"ca", (v >> 6) & 1},
                      {"offset_us", ((v >> 8) & 0x1FFF) * (units ? 300 : 30)},
                      {"phy", phy == 0 ? "1M" : phy == 1 ? "2M" : phy == 2 ? "Coded" : "?"}};
    i += 3;
  }
  if ((flags & 0x20) && i + 18 <= end) i += 18;
  if ((flags & 0x40) && i + 1 <= end && i < by.size()) {
    const int t = by[i];
    out["tx_power_dbm"] = t > 127 ? t - 256 : t;
    i += 1;
  }
  if (end < by.size()) data.assign(by.begin() + static_cast<ptrdiff_t>(end), by.end());
  return out;
}

}  // namespace si::rx::rc_link::coded
