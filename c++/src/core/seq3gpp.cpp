#include "core/seq3gpp.hpp"

#include <cmath>
#include <functional>
#include <stdexcept>

namespace si::seq3gpp {

cdvec lte_pss_freq(int u) {
  // numpy: -1j*pi*u*n*(n+1)/63 -> phase ((-pi*u)*n*(n+1)) * (1/63) (complex division by a
  // real is a multiplication by its reciprocal); exp of a purely imaginary value
  cdvec d(62);
  const double a = -kPi * static_cast<double>(u);
  const double r63 = 1.0 / 63.0;
  for (int n = 0; n < 62; ++n) {
    const double ph = n <= 30 ? (a * n) * (n + 1) * r63 : (a * (n + 1)) * (n + 2) * r63;
    d[n] = cd(std::cos(ph), std::sin(ph));
  }
  return d;
}

namespace {
std::vector<int> mseq(const std::function<int(const std::vector<int>&, int)>& fb) {
  std::vector<int> x{0, 0, 0, 0, 1};
  for (int i = 0; i < 26; ++i) x.push_back(fb(x, i) % 2);
  for (auto& v : x) v = 1 - 2 * v;
  return x;
}
}  // namespace

const std::vector<int>& mseq_s() {
  static const std::vector<int> s = mseq([](const std::vector<int>& x, int i) { return x[i + 2] + x[i]; });
  return s;
}
const std::vector<int>& mseq_c() {
  static const std::vector<int> s = mseq([](const std::vector<int>& x, int i) { return x[i + 3] + x[i]; });
  return s;
}
const std::vector<int>& mseq_z() {
  static const std::vector<int> s =
      mseq([](const std::vector<int>& x, int i) { return x[i + 4] + x[i + 2] + x[i + 1] + x[i]; });
  return s;
}

dvec lte_sss_freq(int nid1, int nid2, int sub) {
  const int qp = nid1 / 30;
  const int q = (nid1 + qp * (qp + 1) / 2) / 30;
  const int mp = nid1 + q * (q + 1) / 2;
  const int m0 = mp % 31, m1 = (m0 + mp / 31 + 1) % 31;
  const auto &S = mseq_s(), &C = mseq_c(), &Z = mseq_z();
  dvec d(62, 0.0);
  for (int n = 0; n < 31; ++n) {
    const int s0 = S[(n + m0) % 31], s1 = S[(n + m1) % 31];
    const int c0 = C[(n + nid2) % 31], c1 = C[(n + nid2 + 3) % 31];
    const int z0 = Z[(n + m0 % 8) % 31], z1 = Z[(n + m1 % 8) % 31];
    if (sub == 0) {
      d[2 * n] = s0 * c0;
      d[2 * n + 1] = s1 * c1 * z0;
    } else {
      d[2 * n] = s1 * c0;
      d[2 * n + 1] = s0 * c1 * z1;
    }
  }
  return d;
}

const std::vector<int>& sc62() {
  static const std::vector<int> v = [] {
    std::vector<int> r;
    for (int k = -31; k < 0; ++k) r.push_back(k);
    for (int k = 1; k < 32; ++k) r.push_back(k);
    return r;
  }();
  return v;
}

const SssBank& sss_bank() {
  static const SssBank bank = [] {
    SssBank b(3);
    for (int nid2 = 0; nid2 < 3; ++nid2)
      for (int s = 0; s < 2; ++s) {
        b[nid2][s].resize(168);
        for (int n1 = 0; n1 < 168; ++n1) b[nid2][s][n1] = lte_sss_freq(n1, nid2, s == 0 ? 0 : 5);
      }
    return b;
  }();
  return bank;
}

namespace {
// x(k+31) = sum_t x(k+t) mod 2 over the taps (python _lfsr31)
std::vector<int8_t> lfsr31(const std::vector<int8_t>& init, const std::vector<int>& taps, size_t n) {
  std::vector<int8_t> x(n + 31, 0);
  for (size_t i = 0; i < 31; ++i) x[i] = init[i];
  for (size_t k = 0; k < n; ++k) {
    int8_t acc = x[k + taps[0]];
    for (size_t t = 1; t < taps.size(); ++t) acc ^= x[k + taps[t]];
    x[k + 31] = acc;
  }
  return x;
}
}  // namespace

std::vector<int8_t> gold(uint32_t c_init, size_t n, size_t nc) {
  std::vector<int8_t> i1(31, 0), i2(31, 0);
  i1[0] = 1;
  for (int i = 0; i < 31; ++i) i2[i] = static_cast<int8_t>((c_init >> i) & 1u);
  const auto x1 = lfsr31(i1, {0, 3}, n + nc);
  const auto x2 = lfsr31(i2, {0, 1, 2, 3}, n + nc);
  std::vector<int8_t> c(n);
  for (size_t k = 0; k < n; ++k) c[k] = static_cast<int8_t>(x1[nc + k] ^ x2[nc + k]);
  return c;
}

cdvec qpsk(const std::vector<int8_t>& bits) {
  if (bits.size() % 2) throw std::invalid_argument("qpsk: odd number of bits");
  const double r = 1.0 / std::sqrt(2.0);
  cdvec s(bits.size() / 2);
  for (size_t i = 0; i < s.size(); ++i)
    s[i] = cd((1 - 2 * bits[2 * i]) * r, (1 - 2 * bits[2 * i + 1]) * r);
  return s;
}

}  // namespace si::seq3gpp
