// Bluetooth BR/EDR access codes (research, signal gaps S7; python: receivers/rc_link/bredr.py).
#include "receivers/rc_link/bredr.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::rc_link::bredr {

namespace {

constexpr uint64_t kPn = 0x83848D96BBCC54FCull;
constexpr uint64_t kG = 0x585713DA9ull;   // 0o260534236651, degree 34
constexpr int kMaxErrors = 3;
constexpr int64_t kSearchBits = 160;
constexpr double kRate = 1e6, kFilterHz = 0.6e6;
constexpr uint64_t kTails[2] = {0b101100, 0b010011};

int bitlen(uint64_t a) { return a ? 64 - __builtin_clzll(a) : 0; }

uint64_t polymod(uint64_t a, uint64_t g) {
  const int dg = bitlen(g) - 1;
  while (bitlen(a) - 1 >= dg) a ^= g << (bitlen(a) - 1 - dg);
  return a;
}

}  // namespace

uint64_t syncword(uint32_t lap) {
  const uint64_t barker = ((lap >> 23) & 1) == 0 ? 0b101100 : 0b010011;
  const uint64_t xt = (static_cast<uint64_t>(lap) | (barker << 24)) ^ (kPn >> 34);
  return ((xt << 34) | polymod(xt << 34, kG)) ^ kPn;
}

std::optional<std::tuple<int64_t, uint32_t, int>> find(const std::vector<int64_t>& b) {
  const int64_t n = static_cast<int64_t>(b.size());
  if (n < 64) return std::nullopt;
  const int64_t last = std::min(n - 64, kSearchBits);
  for (int64_t p = 0; p <= last; ++p) {
    uint64_t top = 0;
    for (int q = 0; q < 6; ++q) top |= static_cast<uint64_t>(b[static_cast<size_t>(p + 58 + q)]) << q;
    if (top != kTails[0] && top != kTails[1]) continue;
    uint64_t w = 0;
    for (int q = 0; q < 64; ++q) w |= static_cast<uint64_t>(b[static_cast<size_t>(p + q)]) << q;
    const uint32_t lap = static_cast<uint32_t>((w >> 34) & 0xFFFFFF);
    const int err = static_cast<int>(std::bitset<64>(w ^ syncword(lap)).count());
    if (err <= kMaxErrors) return std::make_tuple(p, lap, err);
  }
  return std::nullopt;
}

Json try_decode(const cfvec& y, double ofs) {
  const double sps = ofs / kRate;
  if (static_cast<double>(y.size()) < 72 * sps) return nullptr;
  const size_t n0 = y.size();
  cdvec z(n0);
  for (size_t i = 0; i < n0; ++i) z[i] = cd(y[i].real(), y[i].imag());
  const dvec f = fft::fftfreq(n0, 1 / ofs);
  cdvec Z = fft::fft(z);
  for (size_t i = 0; i < n0; ++i)
    if (!(std::fabs(f[i]) <= kFilterHz)) Z[i] = cd(0.0, 0.0);
  z = fft::ifft(Z);
  // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
  const bool elide = (n0 - 1) * sizeof(cd) >= 256 * 1024;
  dvec d(n0 - 1);
  for (size_t i = 1; i < n0; ++i) {
    const cd q = elide ? np::cmul(std::conj(z[i - 1]), z[i]) : np::cmul(z[i], std::conj(z[i - 1]));
    d[i - 1] = std::atan2(q.imag(), q.real());
  }
  const size_t k = static_cast<size_t>(std::max<int64_t>(1, py_round_int(sps)));
  dvec c(d.size() + 1, 0.0);
  for (size_t i = 0; i < d.size(); ++i) c[i + 1] = c[i] + d[i];
  if (c.size() <= k) return nullptr;
  dvec v(c.size() - k);
  for (size_t i = 0; i < v.size(); ++i) v[i] = (c[i + k] - c[i]) / static_cast<double>(k);
  dvec srt = v;
  std::sort(srt.begin(), srt.end());
  const double thr = 0.5 * (srt[static_cast<size_t>(0.1 * static_cast<double>(srt.size()))] +
                            srt[static_cast<size_t>(0.9 * static_cast<double>(srt.size()))]);
  const int64_t n = static_cast<int64_t>((static_cast<double>(v.size()) - sps) / sps);
  for (int ph = 0; ph < static_cast<int>(sps); ++ph) {
    std::vector<int64_t> b(static_cast<size_t>(std::max<int64_t>(n, 0)));
    for (int64_t i = 0; i < n; ++i)
      b[static_cast<size_t>(i)] = v[static_cast<size_t>(static_cast<int64_t>(ph + static_cast<double>(i) * sps))] > thr;
    for (int pol = 0; pol < 2; ++pol) {
      std::vector<int64_t> bits = b;
      if (pol)
        for (auto& x : bits) x = 1 - x;
      const auto r = find(bits);
      if (r) {
        const auto [p, lap, err] = *r;
        Json o = Json::object();
        o["protocol"] = "Bluetooth BR/EDR";
        o["lap"] = fmt("%06X", lap);
        o["access_code_errors"] = static_cast<int64_t>(err);
        o["inquiry"] = lap == kGiac || lap == kLiac;
        o["giac"] = lap == kGiac;
        o["access_code_ok"] = true;
        o["inverted"] = pol == 1;
        o["t_offset_us"] = py_round((ph + static_cast<double>(p) * sps) / ofs * 1e6, 2);
        return o;
      }
    }
  }
  return nullptr;
}

}  // namespace si::rx::rc_link::bredr
