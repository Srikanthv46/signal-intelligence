// python receivers/cofdm/dvbt2_p1.py (see there for the method).
#include "receivers/cofdm/dvbt2_p1.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "core/fft.hpp"
#include "core/npcompat.hpp"

namespace si::rx::cofdm::p1 {
namespace {
constexpr double kPi = 3.14159265358979323846;

std::vector<int> hexbits(const char* h) {
  std::vector<int> b;
  for (const char* c = h; *c; ++c) {
    const int v = std::stoi(std::string(1, *c), nullptr, 16);
    for (int i = 3; i >= 0; --i) b.push_back((v >> i) & 1);
  }
  return b;
}

const std::vector<int>& scr() {
  static const std::vector<int> s = [] {
    int reg[15];
    const char* init = "100111001000110";
    for (int i = 0; i < 15; ++i) reg[i] = init[i] - '0';
    std::vector<int> out(384);
    for (int i = 0; i < 384; ++i) {
      const int b = reg[13] ^ reg[14];
      out[i] = 1 - 2 * b;
      for (int k = 14; k > 0; --k) reg[k] = reg[k - 1];
      reg[0] = b;
    }
    return out;
  }();
  return s;
}

std::string bin(int v, int n) {
  std::string s;
  for (int i = n - 1; i >= 0; --i) s += ((v >> i) & 1) ? '1' : '0';
  return s;
}
}  // namespace

dvec metric(const cdvec& y) {
  const int64_t ny = static_cast<int64_t>(y.size()), n = ny - 2 * kN1;
  if (n <= 0) return {};
  cdvec c1(static_cast<size_t>(ny) + 1), c2(static_cast<size_t>(ny) + 1);
  dvec ce(static_cast<size_t>(ny) + 1);
  c1[0] = c2[0] = cd(0, 0);
  ce[0] = 0;
  for (int64_t m = 0; m < ny; ++m) {
    const cd rot = std::exp(cd(0, -2.0 * kPi * static_cast<double>(m) / kN1));
    const cd q1 = m < ny - kLC ? y[m] * std::conj(y[m + kLC]) * rot : cd(0, 0) * rot;
    const cd q2 = m >= kLB ? y[m] * std::conj(y[m - kLB]) * rot : cd(0, 0) * rot;
    c1[m + 1] = c1[m] + q1;
    c2[m + 1] = c2[m] + q2;
    ce[m + 1] = ce[m] + std::norm(y[m]);
  }
  dvec out(static_cast<size_t>(n));
  for (int64_t s = 0; s < n; ++s) {
    const double num = std::abs(c1[s + kLC] - c1[s]) + std::abs(c2[s + 2 * kN1] - c2[s + kLC + kN1]);
    const double den = 0.5 * ((ce[s + kLC] - ce[s]) + (ce[s + kLC + kN1] - ce[s + kN1]) +
                              (ce[s + 2 * kN1] - ce[s + kLC + kN1]) + (ce[s + kLC + kN1 - kLB] - ce[s + kN1 - kLB])) + 1e-30;
    out[s] = num / den;
  }
  return out;
}

std::optional<Json> decode(const cdvec& y, int64_t s0) {
  if (s0 + kLC + kN1 > static_cast<int64_t>(y.size())) return std::nullopt;
  const cdvec A = fft::fft(cdvec(y.begin() + s0 + kLC, y.begin() + s0 + kLC + kN1));
  dvec P(kN1);
  for (int i = 0; i < kN1; ++i) P[i] = std::norm(A[i]);
  auto idx = [](int k, int sh) { return ((k - 426 + sh) % kN1 + kN1) % kN1; };
  double be = -1;
  int sh = 0;
  for (int s = -100; s <= 100; ++s) {
    dvec v(384);
    for (int i = 0; i < 384; ++i) v[i] = P[idx(kKP1[i], s)];
    const double e = np::pairwise_sum(v.data(), v.size());
    if (be < 0 || e > be) be = e, sh = s;
  }
  cdvec r(384);
  for (int i = 0; i < 384; ++i) r[i] = A[idx(kKP1[i], sh)] * static_cast<double>(scr()[i]);
  dvec d(383), ad(383);
  for (int i = 1; i < 384; ++i) {
    d[i - 1] = std::real(r[i] * std::conj(r[i - 1]));
    ad[i - 1] = std::fabs(d[i - 1]);
  }
  const double sc = np::mean(ad) + 1e-30;
  dvec soft(384);
  soft[0] = 0.0;
  for (int i = 1; i < 384; ++i) soft[i] = d[i - 1] / sc;
  double s1sc[8], s2sc[16];
  for (int v = 0; v < 8; ++v) {
    const auto b = hexbits(kS1Hex[v]);
    double a = 0, c = 0;
    for (int i = 1; i < 64; ++i) a += soft[i] * (1.0 - 2.0 * b[i]);
    for (int i = 0; i < 64; ++i) c += soft[320 + i] * (1.0 - 2.0 * b[i]);
    s1sc[v] = a + c;
  }
  for (int v = 0; v < 16; ++v) {
    const auto b = hexbits(kS2Hex[v]);
    double a = 0;
    for (int i = 0; i < 256; ++i) a += soft[64 + i] * (1.0 - 2.0 * b[i]);
    s2sc[v] = a;
  }
  const int s1 = static_cast<int>(std::max_element(s1sc, s1sc + 8) - s1sc);
  const int s2 = static_cast<int>(std::max_element(s2sc, s2sc + 16) - s2sc);
  dvec t1(s1sc, s1sc + 8), t2(s2sc, s2sc + 16);
  std::sort(t1.begin(), t1.end());
  std::sort(t2.begin(), t2.end());
  static const char* kS1Names[5] = {"T2_SISO", "T2_MISO", "Non-T2 (FEF)", "T2_LITE_SISO", "T2_LITE_MISO"};
  static const char* kBase[8] = {"2K, any GI", "8K, GI 1/32 / 1/16 / 1/8 / 1/4", "4K, any GI", "1K, any GI", "16K, any GI",
                                 "32K, GI 1/32 / 1/16 / 1/8", "8K, GI 1/128 / 19/256 / 19/128", "32K, GI 1/128 / 19/256 / 19/128"};
  static const char* kLite[8] = {"2K, any GI", "8K, GI 1/32 / 1/16 / 1/8 / 1/4", "4K, any GI", "16K, GI 1/128 / 19/256 / 19/128",
                                 "16K, GI 1/32 / 1/16 / 1/8 / 1/4", nullptr, "8K, GI 1/128 / 19/256 / 19/128", nullptr};
  const int f1 = s2 >> 1;
  std::string fft_gi;
  if (s1 == 2) fft_gi = f1 == 0 ? "undefined FEF part" : "reserved";
  else if (s1 == 3 || s1 == 4) fft_gi = kLite[f1] ? kLite[f1] : "reserved";
  else fft_gi = kBase[f1];
  Json out = Json::object();
  out["s1"] = bin(s1, 3);
  out["s2"] = bin(s2, 4);
  out["preamble"] = s1 < 5 ? kS1Names[s1] : "reserved";
  out["fft_gi"] = fft_gi;
  out["mixed"] = (s2 & 1) != 0;
  out["carrier_offset"] = sh;
  out["s1_margin"] = py_round((t1[7] - t1[6]) / 127, 3);
  out["s2_margin"] = py_round((t2[15] - t2[14]) / 256, 3);
  return out;
}

std::vector<Found> find(const cdvec& y, double thr) {
  const dvec m = metric(y);
  std::vector<Found> out;
  const int64_t n = static_cast<int64_t>(m.size());
  int64_t k = 0;
  while (k < n) {
    const int64_t e = std::min(n, k + 2 * kN1);
    const int64_t j = static_cast<int64_t>(std::max_element(m.begin() + k, m.begin() + e) - m.begin());
    if (m[j] >= thr) {
      if (auto d = decode(y, j)) out.push_back({j, m[j], *d});
      k = j + 2 * kN1;
    } else {
      k += 2 * kN1;
    }
  }
  return out;
}

}  // namespace si::rx::cofdm::p1
