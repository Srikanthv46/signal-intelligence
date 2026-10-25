// python: receivers/adsb/modes.py (demod, _crc24, _residual).
#include "receivers/adsb/modes.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::adsb {

uint32_t crc24(const uint8_t* bits, size_t n) {
  constexpr uint32_t g = 0xFFF409;
  uint32_t reg = 0;
  for (size_t i = 0; i < n + 24; ++i) {   // the data, then 24 zero bits
    const uint32_t top = (reg >> 23) & 1;
    reg = ((reg << 1) & 0xFFFFFF) | (i < n ? bits[i] : 0);
    if (top) reg ^= g;
  }
  return reg;
}

uint32_t parity_residual(const uint8_t* bits, size_t n) {
  uint32_t p = 0;
  for (size_t i = n - 24; i < n; ++i) p = (p << 1) | bits[i];
  return crc24(bits, n - 24) ^ p;
}

namespace {

std::string hex_of(const uint8_t* bits, size_t n) {
  static const char* kHex = "0123456789ABCDEF";
  std::string s;
  for (size_t i = 0; i < n; i += 4) s += kHex[bits[i] << 3 | bits[i + 1] << 2 | bits[i + 2] << 1 | bits[i + 3]];
  return s;
}

}  // namespace

std::vector<Json> demod(const fvec& env, float thr) {
  constexpr int L = 8 * kSpu;
  double pre[L] = {};
  for (int t8 : {0, 8, 28, 36})   // pulses at 0, 1, 3.5, 4.5 us
    for (int k = 0; k < kSpu / 2; ++k) pre[t8 + k] = 1;
  const double mu = np::mean(pre, L);
  for (double& p : pre) p -= mu;
  // np.correlate(env, pre, "valid") in float64
  const size_t N = env.size() >= L ? env.size() - L + 1 : 0;
  dvec c(N);
  for (size_t s = 0; s < N; ++s) {
    double acc = 0;
    for (int k = 0; k < L; ++k) acc += static_cast<double>(env[s + k]) * pre[k];
    c[s] = acc;
  }
  double noise = static_cast<double>(np::median(env));
  if (noise == 0) noise = 1e-12;
  const float noise_f = static_cast<float>(noise);

  std::vector<Json> out;
  int64_t busy = -1;
  uint8_t bits[112];
  float a_on[112], b_on[112], mx[112];
  for (size_t s = 0; s < N; ++s) {
    if (!(c[s] > static_cast<double>(thr))) continue;
    double m = c[s];   // c == max(np.roll(c, k), k = -3..3): wraps around the ends
    for (int k = -3; k <= 3; ++k) m = std::max(m, c[(static_cast<int64_t>(s) - k + static_cast<int64_t>(N) * 4) % static_cast<int64_t>(N)]);
    if (c[s] != m) continue;
    if (static_cast<int64_t>(s) < busy) continue;
    const size_t d0 = s + 8 * kSpu;
    if (d0 + 112 * kSpu > env.size()) break;
    const float* a = env.data() + d0;
    for (int r = 0; r < 112; ++r) {
      a_on[r] = np::pairwise_sum(a + r * kSpu, kSpu / 2);
      b_on[r] = np::pairwise_sum(a + r * kSpu + kSpu / 2, kSpu / 2);
      bits[r] = a_on[r] > b_on[r];
    }
    const int df = bits[0] << 4 | bits[1] << 3 | bits[2] << 2 | bits[3] << 1 | bits[4];
    const int n = df >= 16 ? 112 : 56;
    const uint32_t res = parity_residual(bits, n);
    bool ok;
    if (df == 17 || df == 18) ok = res == 0;
    else if (df == 11) ok = res < 128;
    else ok = df == 0 || df == 4 || df == 5 || df == 16 || df == 20 || df == 21;
    if (!ok) continue;
    for (int r = 0; r < n; ++r) mx[r] = std::max(a_on[r] / 4.0f, b_on[r] / 4.0f);   // np.maximum of the means
    const float on = np::mean(mx, static_cast<size_t>(n));
    const float snr = 20.0f * std::log10(on / noise_f);
    Json f = Json::object();
    f["t_s"] = static_cast<double>(s) / kFs;
    f["df"] = df;
    f["bits"] = n;
    f["hex"] = hex_of(bits, n);
    f["parity_residual"] = fmt("%06X", res);
    f["snr_db"] = py_round(static_cast<double>(snr), 1);
    if (df == 11 || df == 17 || df == 18) f["aa"] = hex_of(bits + 8, 24);
    else f["aa"] = nullptr;
    out.push_back(std::move(f));
    busy = static_cast<int64_t>(d0) + n * kSpu;
  }
  // address/parity frames are kept only for addresses confirmed by DF11/17/18 (dump1090 rule)
  std::set<std::string> known;
  for (const Json& m : out)
    if (m["aa"].is_string()) known.insert(m["aa"].get<std::string>());
  std::vector<Json> kept;
  for (Json& m : out) {
    const int df = m["df"].get<int>();
    if (df != 11 && df != 17 && df != 18) {
      m["confirmed"] = known.count(m["parity_residual"].get<std::string>()) > 0;
      m["aa"] = m["parity_residual"];
    } else {
      m["confirmed"] = true;
    }
    if (m["confirmed"].get<bool>()) kept.push_back(std::move(m));
  }
  return kept;
}

}  // namespace si::rx::adsb
