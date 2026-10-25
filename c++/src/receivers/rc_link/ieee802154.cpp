// python receivers/rc_link/ieee802154.py (see there for the sources and the receiver).
#include "receivers/rc_link/ieee802154.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <string>

namespace si::rx::rc_link::ieee802154 {
namespace {
constexpr const char* kChips[16] = {
    "11011001110000110101001000101110", "11101101100111000011010100100010", "00101110110110011100001101010010",
    "00100010111011011001110000110101", "01010010001011101101100111000011", "00110101001000101110110110011100",
    "11000011010100100010111011011001", "10011100001101010010001011101101", "10001100100101100000011101111011",
    "10111000110010010110000001110111", "01111011100011001001011000000111", "01110111101110001100100101100000",
    "00000111011110111000110010010110", "01100000011101111011100011001001", "10010110000001110111101110001100",
    "11001001011000000111011110111000"};
constexpr int kSps = 4;

double sgn(double v) { return v > 0 ? 1.0 : v < 0 ? -1.0 : 0.0; }

// python modulate(chips, sps): O-QPSK half-sine, even chips on I, odd on Q one chip later
std::vector<std::complex<double>> modulate(const std::vector<int>& chips) {
  const size_t n = chips.size() * kSps + 2 * kSps;
  std::vector<double> i(n, 0.0), q(n, 0.0);
  for (size_t k = 0; k + 1 < chips.size(); k += 2) {
    const double a = 2.0 * chips[k] - 1, b = 2.0 * chips[k + 1] - 1;
    for (int t = 0; t < 2 * kSps; ++t) {
      const double p = std::sin(M_PI * t / (2.0 * kSps));
      i[k * kSps + t] += a * p;
      q[(k + 1) * kSps + t] += b * p;
    }
  }
  std::vector<std::complex<double>> y(n);
  for (size_t k = 0; k < n; ++k)   // python casts to complex64
    y[k] = std::complex<double>(static_cast<float>(i[k]), static_cast<float>(q[k]));
  return y;
}

// python _disc: angle(y[n] y*[n-1]) summed over each run of sps
template <typename C>
std::vector<double> disc(const std::vector<C>& y) {
  std::vector<double> d(y.size() > 1 ? y.size() - 1 : 0);
  for (size_t k = 1; k < y.size(); ++k) {
    const std::complex<double> a(y[k].real(), y[k].imag()), b(y[k - 1].real(), y[k - 1].imag());
    d[k - 1] = std::arg(a * std::conj(b));
  }
  std::vector<double> out(d.size() >= kSps ? d.size() - kSps + 1 : 0);
  for (size_t k = 0; k < out.size(); ++k) {
    double s = 0;
    for (int j = 0; j < kSps; ++j) s += d[k + j];
    out[k] = s;
  }
  return out;
}

const std::array<std::array<double, 31>, 16>& templates() {
  static const std::array<std::array<double, 31>, 16> T = [] {
    std::array<std::array<double, 31>, 16> t{};
    for (int s = 0; s < 16; ++s) {
      std::vector<int> ch;
      for (int r : {0, s, 0})
        for (int c = 0; c < 32; ++c) ch.push_back(kChips[r][c] - '0');
      const auto d = disc(modulate(ch));
      const int base = 32 * kSps + kSps;
      for (int k = 0; k < 31; ++k) t[s][k] = sgn(d[static_cast<size_t>(base + k * kSps - 1)]);
    }
    return t;
  }();
  return T;
}

std::string hexrev(const std::vector<uint8_t>& m, size_t i, size_t n) {
  std::string s;
  char b[3];
  for (size_t k = 0; k < n; ++k) std::snprintf(b, sizeof b, "%02x", m[i + n - 1 - k]), s += b;
  return s;
}
}  // namespace

std::optional<int> channel(double f_hz, double tol_hz) {
  const int k = static_cast<int>(py_round_int((f_hz - 2405e6) / 5e6)) + 11;
  if (k >= 11 && k <= 26 && std::fabs(f_hz - (2405e6 + 5e6 * (k - 11))) <= tol_hz) return k;
  return std::nullopt;
}

uint16_t crc16_kermit(const uint8_t* p, size_t n) {
  unsigned c = 0;
  for (size_t k = 0; k < n; ++k) {
    c ^= p[k];
    for (int i = 0; i < 8; ++i) c = (c & 1) ? (c >> 1) ^ 0x8408u : c >> 1;
  }
  return static_cast<uint16_t>(c);
}

Json parse_mac(const std::vector<uint8_t>& m) {
  Json out = Json::object();
  if (m.size() < 3) return out;
  const unsigned fcf = m[0] | (m[1] << 8);
  const unsigned ftype = fcf & 7, pan_comp = (fcf >> 6) & 1, dmode = (fcf >> 10) & 3, smode = (fcf >> 14) & 3;
  static const char* kTypes[] = {"beacon", "data", "ack", "mac command"};
  out["frame_type"] = ftype < 4 ? Json(kTypes[ftype]) : Json(ftype);
  out["seq"] = m[2];
  out["security"] = ((fcf >> 3) & 1) != 0;
  out["ack_request"] = ((fcf >> 5) & 1) != 0;
  size_t i = 3;
  auto take = [&](size_t n) -> Json {
    Json v = i + n <= m.size() ? Json(hexrev(m, i, n)) : Json(nullptr);
    i += n;
    return v;
  };
  if (dmode) {
    out["dst_pan"] = take(2);
    out["dst_addr"] = take(dmode == 2 ? 2 : 8);
  }
  if (smode) {
    if (!(pan_comp && dmode)) out["src_pan"] = take(2);
    out["src_addr"] = take(smode == 2 ? 2 : 8);
  }
  out["payload_bytes"] = m.size() > i ? static_cast<int64_t>(m.size() - i) : 0;
  for (const char* k : {"src_addr", "dst_addr"})
    if (out.contains(k) && out[k].is_string() && out[k].get<std::string>().size() == 16 &&
        out[k].get<std::string>().rfind("0013a2", 0) == 0) {
      out["vendor"] = "Digi International (XBee)";
      break;
    }
  return out;
}

Json demod(const cfvec& y, double ofs) {
  if (py_round_int(ofs / kChipRate) != kSps) return nullptr;
  const auto& T = templates();
  const std::vector<double> d = disc(y);
  const int64_t nd = static_cast<int64_t>(d.size());
  const int64_t n_sym = (nd - kSps) / (32 * kSps) - 1;
  if (n_sym < 12) return nullptr;
  const int64_t span = std::max<int64_t>(1, std::min<int64_t>(nd - 32 * kSps * 10, 32 * kSps * 4));
  int64_t off0 = 0;
  double bestv = -1e300;
  for (int64_t o = 0; o < span; ++o) {
    double v = 0;
    for (int k = 0; k < 31; ++k) v += sgn(d[static_cast<size_t>(o + k * kSps)]) * T[0][static_cast<size_t>(k)];
    if (v > bestv) bestv = v, off0 = o;   // python argmax: first maximum
  }
  if (bestv < 20) return nullptr;
  auto sym_at = [&](int64_t k) -> int {
    const int64_t b = off0 + k * 32 * kSps;
    if (b + 30 * kSps >= nd) return -1;
    int j = 0;
    double best = -1e300;
    for (int s = 0; s < 16; ++s) {
      double v = 0;
      for (int c = 0; c < 31; ++c) v += T[static_cast<size_t>(s)][static_cast<size_t>(c)] * d[static_cast<size_t>(b + c * kSps)];
      if (v > best) best = v, j = s;
    }
    return j;
  };
  int64_t k = 0;
  for (; k < 12; ++k)
    if (sym_at(k) == 7 && sym_at(k + 1) == 0xA) break;
  if (k == 12) return nullptr;
  k += 2;
  std::vector<int> syms;
  int plen = -1;
  for (;;) {
    const int s = sym_at(k);
    if (s < 0) break;
    syms.push_back(s);
    ++k;
    if (syms.size() == 2) {
      plen = (syms[0] | (syms[1] << 4)) & 0x7F;
      if (plen < 5) return nullptr;
    }
    if (plen >= 0 && static_cast<int>(syms.size()) >= 2 + 2 * plen) break;
  }
  if (plen < 0 || static_cast<int>(syms.size()) < 2 + 2 * plen) return nullptr;
  std::vector<uint8_t> by(static_cast<size_t>(plen));
  for (int i = 0; i < plen; ++i) by[static_cast<size_t>(i)] = static_cast<uint8_t>(syms[2 + 2 * i] | (syms[3 + 2 * i] << 4));
  const unsigned fcs = by[by.size() - 2] | (by[by.size() - 1] << 8);
  const bool ok = crc16_kermit(by.data(), by.size() - 2) == fcs;
  std::string hx;
  char b[3];
  for (uint8_t v : by) std::snprintf(b, sizeof b, "%02x", v), hx += b;
  Json out = Json::object();
  out["protocol"] = "IEEE 802.15.4 (O-QPSK 2.4 GHz)";
  out["length"] = plen;
  out["crc_ok"] = ok;
  out["mpdu_hex"] = hx;
  if (ok) out.update(parse_mac(std::vector<uint8_t>(by.begin(), by.end() - 2)));
  return out;
}

}  // namespace si::rx::rc_link::ieee802154
