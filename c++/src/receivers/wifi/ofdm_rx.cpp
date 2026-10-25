// Port of python receivers/wifi/ofdm_rx.py (see ofdm_rx.hpp). Arithmetic follows numpy's:
// complex products through np::cmul, divisions through np::cdiv, |z| through np::cabs, sums
// pairwise where numpy sums pairwise, np.round where Python rounds a numpy float and
// py_round where it rounds a Python float. Phase ramps exp(-2j*pi*f/fs*n) are built as
// Python builds them (a Python complex scalar times an integer array).
#include "receivers/wifi/ofdm_rx.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"

namespace si::rx::wifi::ofdm {

namespace {

constexpr double kPi = 3.141592653589793;
using np::cabs;
using np::cabs2;
using np::cdiv;
using np::cmul;

inline cd expj(double x) { return {std::cos(x), std::sin(x)}; }
inline double angle(cd z) { return std::atan2(z.imag(), z.real()); }
inline cd conjd(cd z) { return {z.real(), -z.imag()}; }
// imag part of the Python complex (-2j * pi * f / fs): the per-sample phase step
inline double phase_step(double f, double fs) { return (-2.0 * kPi * f) / fs; }

std::vector<int> range_sc(int lo, int hi, bool skip_dc, bool skip_pilots) {
  std::vector<int> v;
  for (int k = lo; k <= hi; ++k) {
    if (skip_dc && k == 0) continue;
    if (skip_pilots && (k == -21 || k == -7 || k == 7 || k == 21)) continue;
    v.push_back(k);
  }
  return v;
}

std::vector<double> pilot_polarity() {
  int s[7] = {1, 1, 1, 1, 1, 1, 1};
  std::vector<double> out;
  for (int i = 0; i < 127; ++i) {
    const int b = s[6] ^ s[3];
    out.push_back(1 - 2 * b);
    for (int j = 6; j > 0; --j) s[j] = s[j - 1];
    s[0] = b;
  }
  return out;
}

inline int parity(int v) { return __builtin_popcount(static_cast<unsigned>(v)) & 1; }
constexpr int kG0 = 0133, kG1 = 0171;

struct Trellis {
  double out0[64][2], out1[64][2];
  Trellis() {
    for (int s = 0; s < 64; ++s) {
      const int b = s >> 5;
      const int p0 = 2 * (s & 31), p1 = p0 + 1;
      const int r0 = (b << 6) | p0, r1 = (b << 6) | p1;
      out0[s][0] = 2 * parity(r0 & kG0) - 1;
      out0[s][1] = 2 * parity(r0 & kG1) - 1;
      out1[s][0] = 2 * parity(r1 & kG0) - 1;
      out1[s][1] = 2 * parity(r1 & kG1) - 1;
    }
  }
};
const Trellis& trellis() {
  static const Trellis t;
  return t;
}

// RATE bits R1..R4 -> (Mbps, N_DBPS); (0, 0) if invalid
std::pair<int, int> rate_of(const uint8_t* b) {
  const int k = b[0] << 3 | b[1] << 2 | b[2] << 1 | b[3];
  switch (k) {
    case 0b1101: return {6, 24};
    case 0b1111: return {9, 36};
    case 0b0101: return {12, 48};
    case 0b0111: return {18, 72};
    case 0b1001: return {24, 96};
    case 0b1011: return {36, 144};
    case 0b0001: return {48, 192};
    case 0b0011: return {54, 216};
    default: return {0, 0};
  }
}

int64_t bits_to_int(const std::vector<uint8_t>& b, size_t lo, size_t hi) {
  int64_t v = 0;
  for (size_t i = lo; i < hi; ++i) v |= static_cast<int64_t>(b[i]) << (i - lo);
  return v;
}

// one 3.2 us L-LTF symbol in time (python LTF_TIME)
const cdvec& ltf_time() {
  static const cdvec t = [] {
    cdvec X(64, cd(0, 0));
    for (int k = -26; k <= 26; ++k) X[bin_of(k)] = cd(kLtf[k + 26], 0);
    cdvec y = fft::ifft(X);
    const double s = std::sqrt(52.0);
    for (auto& z : y) z = cdiv(cmul(z, cd(64, 0)), cd(s, 0));
    return y;
  }();
  return t;
}

const std::vector<int>& intlv() {  // python _INTLV
  static const std::vector<int> v = [] {
    std::vector<int> o(48);
    for (int k = 0; k < 48; ++k) o[k] = 3 * (k % 16) + k / 16;
    return o;
  }();
  return v;
}

// np.convolve(a, ones(L) * w, "same") for an odd or even box (n_left = L // 2), summed in order
dvec box_same(const double* a, size_t n, size_t L, double w) {
  dvec out(n);
  const int64_t nl = static_cast<int64_t>(L / 2);
  for (int64_t k = 0; k < static_cast<int64_t>(n); ++k) {
    double s = 0;
    bool first = true;
    for (int64_t i = k - nl; i < k - nl + static_cast<int64_t>(L); ++i) {
      if (i < 0 || i >= static_cast<int64_t>(n)) continue;
      const double t = a[i] * w;
      s = first ? t : s + t;
      first = false;
    }
    out[k] = s;
  }
  return out;
}

// scipy.signal.resample_poly(x, up, down) for complex64, bit-identical to sig::resample_poly
// (same filter, same products and the same order of additions per output sample), computing
// only the kept outputs with register accumulators: output y reads input
// floor(y*down/up) - hpp + 1 .. floor(y*down/up) with filter phase (y*down) mod up.
}  // namespace

namespace detail {
cfvec resample_cf(const cfvec& x, int up_, int down_) {
  const int64_t g = std::gcd(up_, down_);
  const int64_t up = up_ / g, down = down_ / g;
  if (up == 1 && down == 1) return x;
  const int64_t n_in = static_cast<int64_t>(x.size());
  int64_t n_out = n_in * up;
  n_out = n_out / down + (n_out % down ? 1 : 0);
  const int64_t max_rate = std::max(up, down);
  const int64_t half_len = 10 * max_rate;
  const dvec hd = sig::firwin(static_cast<size_t>(2 * half_len + 1), 1.0 / static_cast<double>(max_rate), std::nullopt,
                              sig::Window::kaiser(5.0));
  auto output_len = [&](int64_t len_h) {
    const int64_t in_len_copy = n_in + (len_h + ((up - len_h % up) % up)) / up - 1;
    const int64_t nt = in_len_copy * up;
    return nt / down + (nt % down > 0 ? 1 : 0);
  };
  const int64_t n_pre_pad = down - half_len % down;
  int64_t n_post_pad = 0;
  const int64_t n_pre_remove = (half_len + n_pre_pad) / down;
  const int64_t len_h0 = static_cast<int64_t>(hd.size());
  while (output_len(len_h0 + n_pre_pad + n_post_pad) < n_out + n_pre_remove) ++n_post_pad;
  std::vector<float> h(static_cast<size_t>(n_pre_pad), 0.0f);
  for (double v : hd) h.push_back(static_cast<float>(v) * static_cast<float>(up));
  h.resize(h.size() + static_cast<size_t>(n_post_pad), 0.0f);
  const int64_t lh = static_cast<int64_t>(h.size());
  const int64_t h_padlen = lh + ((up - lh % up) % up);
  const int64_t hpp = h_padlen / up;
  std::vector<float> htf(static_cast<size_t>(h_padlen), 0.0f);
  for (int64_t ph = 0; ph < up; ++ph)
    for (int64_t j = 0; j < hpp; ++j) {
      const int64_t src = (hpp - 1 - j) * up + ph;
      htf[ph * hpp + j] = src < lh ? h[src] : 0.0f;
    }
  const int64_t len_out = output_len(lh);
  const int64_t padded_len = n_in + hpp - 1;
  cfvec out(static_cast<size_t>(n_out), cf(0, 0));
  for (int64_t k = 0; k < n_out; ++k) {
    const int64_t y = k + n_pre_remove;
    if (y >= len_out) break;
    const int64_t xi = (y * down) / up, t = (y * down) % up;
    if (xi >= padded_len) break;
    const float* hp = &htf[t * hpp];
    const int64_t x0 = xi - hpp + 1;
    const int64_t j0 = std::max<int64_t>(0, -x0), j1 = std::min<int64_t>(hpp, n_in - x0);
    float ar = 0.0f, ai = 0.0f;
    for (int64_t j = j0; j < j1; ++j) {
      const cf v = x[x0 + j];
      ar = ar + v.real() * hp[j];
      ai = ai + v.imag() * hp[j];
    }
    out[k] = cf(ar, ai);
  }
  return out;
}
}  // namespace detail

const std::vector<int> kPilotSc = {-21, -7, 7, 21};
const double kPilotBase[4] = {1, 1, 1, -1};
const std::vector<int> kDataSc = range_sc(-26, 26, true, true);
const std::vector<int> kUsedSc = range_sc(-26, 26, true, false);
const std::vector<int> kHtSc = range_sc(-28, 28, false, false);
const std::vector<int> kHtDataSc = range_sc(-28, 28, true, true);
const std::vector<int> kHtUsedSc = range_sc(-28, 28, true, false);
const double kLtf[53] = {1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0,
                         1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1};
const double kHtLtf[57] = {1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0,
                           1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1};
const std::vector<double>& pilot_pol() {
  static const std::vector<double> p = pilot_polarity();
  return p;
}

// ------------------------------------------------------------------ coding
std::vector<uint8_t> viterbi(const dvec& soft, bool terminated) {
  const Trellis& T = trellis();
  const size_t n = soft.size() / 2;
  double metric[64], next[64];
  std::fill(metric, metric + 64, -1e18);
  metric[0] = 0;
  std::vector<uint8_t> took(n * 64);
  for (size_t t = 0; t < n; ++t) {
    const double y0 = soft[2 * t], y1 = soft[2 * t + 1];
    uint8_t* tk = &took[t * 64];
    for (int s = 0; s < 64; ++s) {
      const int p0 = 2 * (s & 31);
      const double m0 = metric[p0] + (T.out0[s][0] * y0 + T.out0[s][1] * y1);
      const double m1 = metric[p0 + 1] + (T.out1[s][0] * y0 + T.out1[s][1] * y1);
      const bool win = m1 > m0;
      next[s] = win ? m1 : m0;
      tk[s] = win;
    }
    std::copy(next, next + 64, metric);
  }
  int s = terminated ? 0 : static_cast<int>(std::max_element(metric, metric + 64) - metric);
  std::vector<uint8_t> bits(n);
  for (size_t t = n; t-- > 0;) {
    bits[t] = static_cast<uint8_t>(s >> 5);
    s = 2 * (s & 31) + took[t * 64 + s];
  }
  return bits;
}

std::vector<uint8_t> conv_encode(const std::vector<uint8_t>& bits) {
  int s = 0;
  std::vector<uint8_t> out;
  out.reserve(2 * bits.size());
  for (uint8_t b : bits) {
    const int reg = (b << 6) | s;
    out.push_back(static_cast<uint8_t>(parity(reg & kG0)));
    out.push_back(static_cast<uint8_t>(parity(reg & kG1)));
    s = reg >> 1;
  }
  return out;
}

std::vector<uint8_t> crc8_sig(const uint8_t* bits, size_t n) {
  int crc = 0xFF;
  for (size_t i = 0; i < n; ++i) {
    const int fb = ((crc >> 7) & 1) ^ bits[i];
    crc = (crc << 1) & 0xFF;
    if (fb) crc ^= 0x07;
  }
  crc = ~crc & 0xFF;
  std::vector<uint8_t> out(8);
  for (int i = 0; i < 8; ++i) out[i] = static_cast<uint8_t>((crc >> (7 - i)) & 1);
  return out;
}

std::optional<HeSigA> parse_he_sig_a(const std::vector<uint8_t>& b, int lsig_length) {
  const auto c = crc8_sig(b.data(), 42);
  for (int i = 0; i < 4; ++i)
    if (c[i] != b[42 + i]) return std::nullopt;
  for (int i = 46; i < 52; ++i)
    if (b[i]) return std::nullopt;
  static const int kBw[4] = {20, 40, 80, 160};
  auto g = [&](size_t i, size_t j) { return static_cast<int>(bits_to_int(b, i, j)); };
  HeSigA h;
  const int tx = g(26, 33);  // A2 B0-B6
  if (tx != 127) h.txop_us = (tx & 1) ? 512 + 128 * (tx >> 1) : 8 * (tx >> 1);
  if (lsig_length % 3 == 2) {
    const int bwi = g(15, 18);
    h.format = "MU";
    h.uplink = b[0] != 0;
    h.bss_color = g(5, 11);
    h.bw = bwi < 4 ? kBw[bwi] : (bwi == 4 || bwi == 5 ? 80 : 160);
  } else if (b[0]) {
    h.format = "SU";
    h.uplink = b[2] != 0;
    h.mcs = g(3, 7);
    h.dcm = b[7];
    h.bss_color = g(8, 14);
    h.bw = kBw[g(19, 21)];
    h.nsts = g(23, 26) + 1;
    h.ldpc = b[26 + 7];
  } else {
    h.format = "TB";
    h.uplink = true;
    h.bss_color = g(1, 7);
    h.bw = kBw[g(24, 26)];
  }
  return h;
}

// ------------------------------------------------------------------ channels
const std::vector<StdChannel>& all_channels() {
  static const std::vector<StdChannel> c = [] {
    std::vector<StdChannel> v;
    for (int ch = 1; ch <= 14; ++ch) v.push_back({"2g" + std::to_string(ch), ch, ch == 14 ? 2484e6 : 2407e6 + 5e6 * ch});
    std::vector<int> c5;
    for (int ch = 36; ch < 65; ch += 4) c5.push_back(ch);
    for (int ch = 100; ch < 145; ch += 4) c5.push_back(ch);
    for (int ch = 149; ch < 178; ch += 4) c5.push_back(ch);
    for (int ch : c5) v.push_back({"5g" + std::to_string(ch), ch, 5000e6 + 5e6 * ch});
    for (int ch = 1; ch < 234; ch += 4) v.push_back({"6g" + std::to_string(ch), ch, 5950e6 + 5e6 * ch});
    return v;
  }();
  return c;
}

std::vector<Channel> channel_plan(double fc_hz, double fs_hz, double usable_bw_hz, bool allow_partial,
                                  double max_missing_hz) {
  const double half = std::min(usable_bw_hz, fs_hz) / 2;
  const double spill = allow_partial ? max_missing_hz : 0.0;
  std::vector<Channel> out;
  for (const auto& c : all_channels()) {
    if (!(std::fabs(c.freq_hz - fc_hz) + 8.4e6 <= half + spill)) continue;
    out.push_back({c.name, c.ch, c.freq_hz, c.freq_hz - fc_hz, half, std::fabs(c.freq_hz - fc_hz) + 8.4e6 > half});
  }
  return out;
}

std::optional<double> wide_center(const std::string& name, int bw_mhz) {
  if (name.rfind("5g", 0) != 0 || bw_mhz <= 20) return std::nullopt;
  const int ch = std::stoi(name.substr(2));
  static const std::vector<std::pair<int, int>> g40 = {{36, 40},   {44, 48},   {52, 56},   {60, 64},  {100, 104},
                                                       {108, 112}, {116, 120}, {124, 128}, {132, 136}, {140, 144},
                                                       {149, 153}, {157, 161}, {165, 169}, {173, 177}};
  static const std::vector<std::pair<int, int>> g80 = {{36, 48},   {52, 64},   {100, 112}, {116, 128},
                                                       {132, 144}, {149, 161}, {165, 177}};
  static const std::vector<std::pair<int, int>> g160 = {{36, 64}, {100, 128}, {149, 177}};
  const auto* g = bw_mhz == 40 ? &g40 : bw_mhz == 80 ? &g80 : bw_mhz == 160 ? &g160 : nullptr;
  if (!g) return std::nullopt;
  for (const auto& [lo, hi] : *g)
    if (lo <= ch && ch <= hi) return 5000e6 + 5e6 * (lo + hi) / 2;
  return std::nullopt;
}

std::pair<int64_t, int64_t> limit_denominator(double x, int64_t max_den) {
  using i128 = __int128;
  int e = 0;
  const double m = std::frexp(x, &e);
  int64_t n = static_cast<int64_t>(std::ldexp(m, 53));
  int ex = e - 53;
  while (n != 0 && n % 2 == 0 && ex < 0) {
    n /= 2;
    ++ex;
  }
  if (ex >= 0) return {n << ex, 1};
  const int64_t d = int64_t(1) << (-ex);
  if (d <= max_den) return {n, d};
  i128 p0 = 0, q0 = 1, p1 = 1, q1 = 0, nn = n, dd = d;
  while (true) {
    const i128 a = nn / dd;
    const i128 q2 = q0 + a * q1;
    if (q2 > max_den) break;
    const i128 np1 = p0 + a * p1;
    p0 = p1;
    q0 = q1;
    p1 = np1;
    q1 = q2;
    const i128 t = nn - a * dd;
    nn = dd;
    dd = t;
  }
  const i128 k = (max_den - q0) / q1;
  const i128 bp = p0 + k * p1, bq = q0 + k * q1;
  auto absd = [](i128 v) { return v < 0 ? -v : v; };
  // |p1/q1 - n/d| <= |bp/bq - n/d|  <=>  |p1 d - n q1| * bq <= |bp d - n bq| * q1
  const i128 e2 = absd(p1 * d - i128(n) * q1) * bq, e1 = absd(bp * d - i128(n) * bq) * q1;
  i128 rp = e2 <= e1 ? p1 : bp, rq = e2 <= e1 ? q1 : bq;
  const i128 g = std::gcd(static_cast<int64_t>(rp), static_cast<int64_t>(rq));
  return {static_cast<int64_t>(rp / g), static_cast<int64_t>(rq / g)};
}

// ------------------------------------------------------------------ receiver
struct WifiDetector::Det {
  int64_t start20 = 0, end20 = 0;
  int dur_us = 0;
  std::string fmt;
  int mbps = 0, length = 0, bw = 20;
  std::optional<int> mcs;
  std::optional<bool> crc_ok;
  double snr = 0, cfo = 0, stf = 0, noise_pwr = 0;
  Packet ex;  // the extra fields (sig_ok .. lsig_bit_agreement)
  std::shared_ptr<const RxState> rx;
};

WifiDetector::WifiDetector(double fs, double fc, const DetectorConfig& cfg, std::optional<std::vector<Channel>> channels)
    : fs_(fs), fc_(fc), cfg_(cfg) {
  power_margin_ = std::pow(10.0, cfg.power_margin_db / 10);
  channels_ = channels ? *channels : channel_plan(fc, fs, cfg.usable_bw);
  const auto r = limit_denominator(kFsWifi / fs, 4096);
  up_ = static_cast<int>(r.first);
  down_ = static_cast<int>(r.second);
  ratio_ = fs / kFsWifi;
}

cdvec WifiDetector::baseband(const cfvec& x, double offset_hz) const {
  const double b = phase_step(offset_hz, fs_);
  cfvec y(x.size());
  for (size_t n = 0; n < x.size(); ++n) {
    const double ph = b * static_cast<double>(n);
    const cf e(static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
    y[n] = cmul(x[n], e);
  }
  const cfvec z = detail::resample_cf(y, up_, down_);
  return cdvec(z.begin(), z.end());
}

std::vector<WifiDetector::Det> WifiDetector::detect_channel(std::shared_ptr<const cdvec> x20p, const Channel& chan,
                                                            size_t& rejects) const {
  const cdvec& r = *x20p;
  const int64_t n = static_cast<int64_t>(r.size());
  std::vector<Det> out;
  const int64_t L = 48, D = 16;
  if (n < L + D) return out;
  // stf_metric
  const int64_t nc = n - D, np_ = nc - L + 1;
  cdvec c(nc);
  dvec e(nc);
  for (int64_t i = 0; i < nc; ++i) {
    c[i] = cmul(r[i], conjd(r[i + D]));
    e[i] = cabs2(r[i + D]);
  }
  cdvec P(np_);
  dvec R(np_), M(np_), pwr(np_);
  for (int64_t k = 0; k < np_; ++k) {
    cd sp = c[k];
    double sr = e[k];
    for (int64_t j = 1; j < L; ++j) {
      sp += c[k + j];
      sr += e[k + j];
    }
    P[k] = sp;
    R[k] = sr;
    M[k] = cabs(sp) / (sr + 1e-20);
    pwr[k] = sr / 48;
  }
  const double floor = np::percentile(pwr, 10.0);
  std::vector<int64_t> starts, ends;
  bool prev = false;
  for (int64_t k = 0; k <= np_; ++k) {
    const bool ok = k < np_ && (M[k] > cfg_.stf_thresh) && (pwr[k] > floor * power_margin_);
    if (ok && !prev) starts.push_back(k);
    if (!ok && prev) ends.push_back(k);
    prev = ok;
  }
  const double thr = cfg_.min_bit_agreement;
  const cdvec& LT = ltf_time();
  const auto& DB = kDataSc;
  int64_t busy_until = -1;
  for (size_t ri = 0; ri < starts.size(); ++ri) {
    const int64_t s0 = starts[ri], e0 = ends[ri];
    if (e0 - s0 < cfg_.min_plateau || s0 < busy_until) continue;
    // ---- _decode_at
    const double cfo_c = -angle(np::mean(&P[s0], static_cast<size_t>(e0 - s0))) / (2 * kPi * 16) * kFsWifi;
    const double stf = np::mean(&M[s0], static_cast<size_t>(e0 - s0));
    const int64_t lo = s0 + 100, hi = std::min(s0 + 330, n - 700);
    if (hi <= lo) continue;
    const double bc = phase_step(cfo_c, kFsWifi);
    const int64_t ns = hi + 128 - lo;
    cdvec seg(ns);
    for (int64_t i = 0; i < ns; ++i) seg[i] = cmul(r[lo + i], expj(bc * static_cast<double>(lo + i)));
    const int64_t nx = ns - 64 + 1;
    if (nx < 65) continue;
    dvec xc(nx);
    for (int64_t k = 0; k < nx; ++k) {
      cd s(0, 0);
      for (int j = 0; j < 64; ++j) s += cmul(seg[k + j], conjd(LT[j]));
      xc[k] = cabs(s);
    }
    dvec pair(nx - 64);
    for (int64_t k = 0; k + 64 < nx; ++k) pair[k] = xc[k] + xc[k + 64];
    const int64_t t1 = lo + static_cast<int64_t>(np::argmax(pair));
    const int64_t start = t1 - 192;
    if (start < 0) continue;
    cdvec ab(64);
    for (int64_t i = 0; i < 64; ++i) {
      const cd a = cmul(r[t1 + i], expj(bc * static_cast<double>(t1 + i)));
      const cd b = cmul(r[t1 + 64 + i], expj(bc * static_cast<double>(t1 + 64 + i)));
      ab[i] = cmul(b, conjd(a));
    }
    const double cfo = cfo_c - angle(np::sum(ab)) / (2 * kPi * 64) * kFsWifi;
    if (start + 560 + 16 > n) continue;
    const int64_t nrr = std::min<int64_t>(800, n - start);
    const double bf = phase_step(cfo, kFsWifi);
    cdvec rr(nrr);
    for (int64_t i = 0; i < nrr; ++i) rr[i] = cmul(r[start + i], expj(bf * static_cast<double>(start + i)));
    const int back = 3;
    const double s52 = std::sqrt(52.0);
    auto fft_at = [&](int64_t off) {
      cdvec w(rr.begin() + (off - back), rr.begin() + (off - back + 64));
      cdvec Y = fft::fft(w);
      for (auto& z : Y) {
        z = cdiv(z, cd(64, 0));
        z = cmul(z, cd(s52, 0));
      }
      return Y;
    };
    const cdvec Y1 = fft_at(192), Y2 = fft_at(256);
    dvec Lf(64, 0.0);
    for (int k = -26; k <= 26; ++k) Lf[bin_of(k)] = kLtf[k + 26];
    std::vector<int> good_sc, used;
    for (int k : kUsedSc) {
      const double f = chan.offset_hz + k * 312.5e3;
      if (std::fabs(f) <= chan.usable_half_hz - 0.4e6) {
        good_sc.push_back(k);
        used.push_back(bin_of(k));
      }
    }
    const size_t nu = used.size();
    cdvec H(64, cd(std::numeric_limits<double>::infinity(), 0));
    for (int b : used) H[b] = cdiv(cdiv(Y1[b] + Y2[b], cd(2, 0)), cd(Lf[b], 0));
    dvec nk(64, 0.0), nku(nu);
    for (size_t i = 0; i < nu; ++i) {
      const int b = used[i];
      nk[b] = cabs2(Y1[b] - Y2[b]) / 2;
      nku[i] = nk[b];
    }
    const dvec nsm = box_same(nku.data(), nu, 3, 1.0 / 3);
    const double med_nsm = np::median(nsm);
    dvec nk_s(64, 0.0);
    for (size_t i = 0; i < nu; ++i) nk_s[used[i]] = std::max(nsm[i], med_nsm * 0.5);
    const double noise = np::median(nku);
    dvec h2(nu);
    for (size_t i = 0; i < nu; ++i) h2[i] = cabs2(H[used[i]]);
    const double sigp = np::mean(h2);
    const double ratio_sn = sigp / std::max(noise, 1e-20);
    const double snr_db = 10 * std::log10(1e-3 > ratio_sn ? 1e-3 : ratio_sn);
    dvec W(64, 0.0), wu(nu);
    for (size_t i = 0; i < nu; ++i) wu[i] = W[used[i]] = h2[i] / nk_s[used[i]];
    const double wmed4 = np::median(wu) * 4;
    for (size_t i = 0; i < nu; ++i) W[used[i]] = std::min(wu[i], wmed4);
    std::vector<uint8_t> is_used(64, 0);
    for (int b : used) is_used[b] = 1;

    auto equalize = [&](int64_t off, int n_sym) {
      const cdvec Y = fft_at(off);
      cdvec Z(64, cd(0, 0));
      for (int b : used) Z[b] = cdiv(Y[b], H[b]);
      cd pil[4];
      double w[4];
      const double pp = pilot_pol()[n_sym % 127];
      for (int j = 0; j < 4; ++j) {
        const int b = bin_of(kPilotSc[j]);
        pil[j] = cdiv(Z[b], cd(kPilotBase[j] * pp, 0));
        w[j] = W[b];
      }
      int nlive = 0;
      for (double v : w) nlive += v > 0;
      auto wsum = [&]() {
        cd t[4];
        for (int j = 0; j < 4; ++j) t[j] = cmul(cd(w[j], 0), pil[j]);
        return np::pairwise_sum(t, 4);
      };
      if (nlive >= 3) {
        const double ph0 = angle(wsum());
        double err[4];
        for (int j = 0; j < 4; ++j) err[j] = w[j] > 0 ? cabs(cmul(pil[j], expj(-ph0)) - cd(1, 0)) : -1;
        const int am = static_cast<int>(std::max_element(err, err + 4) - err);
        if (err[am] > 0.8) w[am] = 0;
      }
      const double ph = angle(wsum());
      const cd rot = expj(-ph);
      for (auto& z : Z) z = cmul(z, rot);
      return Z;
    };
    auto data_soft = [&](const cdvec& Z, bool q) {
      dvec v(48), out(48);
      for (int i = 0; i < 48; ++i) {
        const int b = bin_of(DB[i]);
        v[i] = (q ? Z[b].imag() : Z[b].real()) * W[b];
      }
      for (int k = 0; k < 48; ++k) out[k] = v[intlv()[k]];
      return out;
    };
    auto live_data = [&](const cdvec& Z) {
      cdvec d;
      for (int k : DB)
        if (is_used[bin_of(k)]) d.push_back(Z[bin_of(k)]);
      return d;
    };

    // ---- L-SIG
    const cdvec Zs = equalize(336, 0);
    const dvec soft0 = data_soft(Zs, false);
    const std::vector<uint8_t> lbits = viterbi(soft0);
    const auto [mbps, ndbps] = rate_of(lbits.data());
    const int length = static_cast<int>(bits_to_int(lbits, 5, 17));
    int par = 0;
    for (int i = 0; i < 17; ++i) par += lbits[i];
    bool tail = false;
    for (int i = 18; i < 24; ++i) tail |= lbits[i] != 0;
    if (mbps == 0 || lbits[4] != 0 || (par & 1) != lbits[17] || tail || length == 0) {
      ++rejects;
      continue;
    }
    const cdvec dl = live_data(Zs);
    dvec ev(dl.size());
    for (size_t i = 0; i < dl.size(); ++i) {
      const double re = dl[i].real(), sg = re > 0 ? 1.0 : re < 0 ? -1.0 : 0.0;
      const double a = re - sg;
      ev[i] = a * a + dl[i].imag() * dl[i].imag();
    }
    const double lsig_evm_db = 10 * std::log10(np::mean(ev) + 1e-12);
    const std::vector<uint8_t> enc = conv_encode(lbits);
    size_t nlive = 0, nagree = 0;
    for (size_t i = 0; i < 48; ++i) {
      if (soft0[i] == 0) continue;
      ++nlive;
      const double sg = soft0[i] > 0 ? 1.0 : -1.0;
      nagree += (2.0 * enc[i] - 1) * sg > 0;
    }
    const double agree = static_cast<double>(nagree) / static_cast<double>(nlive);
    if (agree < thr) {
      ++rejects;
      continue;
    }
    const int dur_us = 20 + 4 * static_cast<int>((16 + 8 * length + 6 + ndbps - 1) / ndbps);

    Det d;
    d.fmt = "NON_HT";
    if (mbps == 6 && start + 560 + 16 <= n) {
      const cdvec Z1 = equalize(416, 1), Z2 = equalize(496, 2);
      const cdvec d0 = live_data(Zs), d1 = live_data(Z1), d2 = live_data(Z2);
      auto qratio = [](const cdvec& z) {
        dvec im(z.size()), re(z.size());
        for (size_t i = 0; i < z.size(); ++i) {
          im[i] = z[i].imag() * z[i].imag();
          re[i] = z[i].real() * z[i].real();
        }
        return np::mean(im) / (np::mean(re) + 1e-12);
      };
      const double q1 = qratio(d1), q2 = qratio(d2);
      const double rl = cabs(np::vdot(d0, d1)) / (np::norm(d0) * np::norm(d1) + 1e-12);
      auto sig_bits = [&](bool q_first) {
        dvec s = data_soft(Z1, q_first);
        const dvec s2 = data_soft(Z2, true);
        s.insert(s.end(), s2.begin(), s2.end());
        return viterbi(s);
      };
      auto crc_pass = [](const std::vector<uint8_t>& bits) {
        const auto c = crc8_sig(bits.data(), 34);
        for (int i = 0; i < 8; ++i)
          if (c[i] != bits[34 + i]) return false;
        for (int i = 42; i < 48; ++i)
          if (bits[i]) return false;
        return true;
      };
      if (rl > 0.85 && (length % 3) != 0) {
        d.fmt = "HE";
        if (start + 640 + 16 <= n) {
          // python he_sig_soft: outer tones +-27 / +-28 trained on L-SIG / RL-SIG (-1 -1 -1 +1), then HE-SIG-A1 / A2
          // (FFT at 496 / 576, pilot polarity p_2 / p_3) on the 52 HT data subcarriers, N_COL 13 deinterleave
          static const int kExt[4] = {-28, -27, 27, 28};
          static const double kExtVal[4] = {-1, -1, -1, 1};
          auto cpe = [&](const cdvec& Y, int n_sym) {
            cd t[4];
            const double pp = pilot_pol()[n_sym];
            for (int j = 0; j < 4; ++j) {
              const int b = bin_of(kPilotSc[j]);
              const cd pil = is_used[b] ? cdiv(cdiv(Y[b], H[b]), cd(kPilotBase[j] * pp, 0)) : cd(0, 0);
              t[j] = cmul(cd(W[b], 0), pil);
            }
            return angle(np::pairwise_sum(t, 4));
          };
          cdvec Hf = H;
          dvec Wf = W;
          const cdvec Yl = fft_at(336), Yr = fft_at(416);
          const cd rl_ = expj(-cpe(Yl, 0)), rr_ = expj(-cpe(Yr, 1));
          dvec wu;
          for (int b : used) wu.push_back(W[b]);
          const bool edge_ok = is_used[bin_of(-26)] && is_used[bin_of(26)];
          for (int j = 0; j < 4; ++j) {
            const int b = bin_of(kExt[j]);
            Hf[b] = cdiv(cdiv(cmul(Yl[b], rl_), cd(kExtVal[j], 0)) + cdiv(cmul(Yr[b], rr_), cd(kExtVal[j], 0)), cd(2, 0));
            Wf[b] = edge_ok ? np::median(wu) : 0.0;
          }
          dvec soft;
          for (int s = 0; s < 2; ++s) {
            const cdvec Y = fft_at(496 + 80 * s);
            const cd rot = expj(-cpe(Y, 2 + s));
            dvec v(52);
            for (int i = 0; i < 52; ++i) {
              const int k = kHtDataSc[i], b = bin_of(k);
              const bool live = (is_used[b] || std::abs(k) >= 27) && Wf[b] != 0;
              v[i] = live ? cmul(cdiv(Y[b], Hf[b]), rot).real() * Wf[b] : 0.0;
            }
            for (int k = 0; k < 52; ++k) soft.push_back(v[4 * (k % 13) + k / 13]);
          }
          const auto he = parse_he_sig_a(viterbi(soft), length);
          d.crc_ok = he.has_value();
          if (he) {
            d.mcs = he->mcs;
            d.bw = he->bw;
            d.ex.he_format = he->format;
            d.ex.he_uplink = he->uplink;
            d.ex.bss_color = he->bss_color;
            d.ex.dcm = he->dcm;
            d.ex.ldpc = he->ldpc;
            d.ex.txop_us = he->txop_us;
            d.ex.he_nsts = he->nsts;
          }
        }
      } else if (q1 > 2.0 && length % 3 == 0) {
        const auto bits = sig_bits(true);
        const bool ok = crc_pass(bits);
        d.crc_ok = ok;
        d.fmt = ok ? "HT" : "HT?";
        if (ok) {
          const int mcs = static_cast<int>(bits_to_int(bits, 0, 7));
          d.mcs = mcs;
          d.bw = bits[7] ? 40 : 20;
          d.ex.sig_ok = true;
          d.ex.ht_length = static_cast<int>(bits_to_int(bits, 8, 24));
          d.ex.aggregation = bits[27];
          d.ex.smoothing = bits[24];
          d.ex.short_gi = bits[31];
          d.ex.ldpc = bits[30];
          d.ex.stbc = static_cast<int>(bits_to_int(bits, 28, 30));
          d.ex.n_ess = static_cast<int>(bits_to_int(bits, 32, 34));
          d.ex.nsts_su = mcs / 8 + 1;
        }
      } else if (q2 > 2.0 && q1 < 0.5 && length % 3 == 0) {
        const auto bits = sig_bits(false);
        const bool ok = crc_pass(bits);
        d.crc_ok = ok;
        d.fmt = ok ? "VHT" : "VHT?";
        if (ok) {
          static const int kBw[4] = {20, 40, 80, 160};
          d.bw = kBw[bits_to_int(bits, 0, 2)];
          d.mcs = static_cast<int>(bits_to_int(bits, 28, 32));
          const int gid = static_cast<int>(bits_to_int(bits, 4, 10));
          d.ex.sig_ok = true;
          d.ex.nsts_su = static_cast<int>(bits_to_int(bits, 10, 13)) + 1;
          d.ex.short_gi = bits[24];
          d.ex.ldpc = bits[26];
          d.ex.stbc = bits[3];
          d.ex.group_id = gid;
          if (gid == 0 || gid == 63) d.ex.partial_aid = static_cast<int>(bits_to_int(bits, 13, 22));
        }
      }
    }
    d.start20 = start;
    d.end20 = start + static_cast<int64_t>(dur_us) * 20;
    d.dur_us = dur_us;
    d.mbps = mbps;
    d.length = length;
    d.snr = snr_db;
    d.cfo = cfo;
    d.stf = stf;
    d.noise_pwr = floor;
    d.ex.lsig_evm_db = py_round(lsig_evm_db, 1);
    d.ex.lsig_bit_agreement = py_round(agree, 3);
    auto rx = std::make_shared<RxState>();
    rx->start = start;
    rx->cfo = cfo;
    rx->back = back;
    rx->used = used;
    rx->good_sc = good_sc;
    rx->H = H;
    rx->W_rel = W;
    rx->nk = nk_s;
    rx->noise = noise;
    rx->r = x20p;
    d.rx = rx;
    out.push_back(std::move(d));
    busy_until = out.back().end20 - 8;
  }
  return out;
}

namespace {
// python WifiDetector._energy_end
std::optional<double> energy_end(const cdvec& x20, int64_t start20, int64_t end20, double noise_pwr) {
  const int64_t k = 20;
  const int64_t stop = std::min<int64_t>(static_cast<int64_t>(x20.size()), end20 + 800);
  const int64_t a = start20 + 400;
  if (stop - a < 2 * k) return std::nullopt;
  const size_t n = static_cast<size_t>(stop - a);
  dvec e(n);
  for (size_t i = 0; i < n; ++i) e[i] = cabs2(x20[a + i]);
  const dvec p = box_same(e.data(), n, 20, 1.0 / 20);
  const size_t g = 80;
  size_t run = 0;
  for (size_t i = 0; i < n; ++i) {
    run = p[i] < noise_pwr * 4 ? run + 1 : 0;
    if (run == g) return static_cast<double>(a + static_cast<int64_t>(i + 1 - g)) / kFsWifi * 1e6;
  }
  return std::nullopt;
}
}  // namespace

bool WifiDetector::accept(Packet& p) const {
  const double ag = p.lsig_bit_agreement;
  const bool end_ok = p.truncated || (p.energy_end_us && std::fabs(*p.energy_end_us - p.end_us) <= cfg_.end_tol_us);
  if (ag >= cfg_.strong_agreement || (p.sig_crc_ok && *p.sig_crc_ok)) return true;
  return ag >= cfg_.min_bit_agreement && end_ok;
}

std::vector<Packet> WifiDetector::run(const cfvec& x, int nthreads) {
  const size_t nch = channels_.size();
  std::vector<std::vector<Det>> dets(nch);
  std::vector<size_t> rej(nch, 0);
  std::vector<std::vector<std::optional<double>>> eends(nch);
  parallel_for(
      nch,
      [&](size_t i) {
        auto x20 = std::make_shared<const cdvec>(baseband(x, channels_[i].offset_hz));
        dets[i] = detect_channel(x20, channels_[i], rej[i]);
        for (const Det& d : dets[i]) eends[i].push_back(energy_end(*x20, d.start20, d.end20, d.noise_pwr));
      },
      nthreads);
  const int64_t n_cap = static_cast<int64_t>(x.size());
  std::vector<Packet> out;
  for (size_t ci = 0; ci < nch; ++ci) {
    n_rejects_ += rej[ci];
    const Channel& chan = channels_[ci];
    for (size_t di = 0; di < dets[ci].size(); ++di) {
      const Det& d = dets[ci][di];
      const int64_t s = static_cast<int64_t>(std::nearbyint(static_cast<double>(d.start20) * ratio_));
      const int64_t e = static_cast<int64_t>(std::nearbyint(static_cast<double>(d.end20) * ratio_));
      Packet p = d.ex;
      const int bw = d.bw;
      const std::optional<double> wc = wide_center(chan.name, bw);
      const double cen = wc ? *wc : chan.freq_hz + d.cfo;
      const bool cen_np = !wc;  // numpy float64 when it carries the CFO
      const double half = bw == 20 ? 8.3e6 : bw == 40 ? 18.1e6 : bw == 80 ? 38.1e6 : 78.1e6;
      const double ch_lo = fc_ - fs_ / 2, ch_hi = fc_ + fs_ / 2;
      p.channel = chan.name;
      p.ch_freq_mhz = chan.freq_hz / 1e6;
      p.center_mhz = np::round((chan.freq_hz + d.cfo) / 1e6, 3);
      p.start_sample = s;
      p.end_sample = std::min(e, n_cap);
      p.start_us = py_round(static_cast<double>(s) / fs_ * 1e6, 3);
      p.duration_us = d.dur_us;
      p.end_us = py_round(static_cast<double>(e) / fs_ * 1e6, 3);
      p.fmt = d.fmt;
      p.lsig_rate_mbps = d.mbps;
      p.lsig_length = d.length;
      p.bandwidth_mhz = bw;
      {
        const double v = cen - half;
        if (ch_lo > v) p.f_lo_mhz = py_round(ch_lo / 1e6, 3);
        else p.f_lo_mhz = cen_np ? np::round(v / 1e6, 3) : py_round(v / 1e6, 3);
        const double w = cen + half;
        if (ch_hi < w) p.f_hi_mhz = py_round(ch_hi / 1e6, 3);
        else p.f_hi_mhz = cen_np ? np::round(w / 1e6, 3) : py_round(w / 1e6, 3);
      }
      p.snr_db = py_round(d.snr, 1);
      p.cfo_khz = py_round(d.cfo / 1e3, 2);
      p.stf_metric = py_round(d.stf, 3);
      p.mcs = d.mcs;
      p.sig_crc_ok = d.crc_ok;
      p.truncated = e > n_cap;
      const auto& ee = eends[ci][di];
      p.energy_end_us = ee ? std::optional<double>(np::round(*ee, 2)) : std::nullopt;
      p.rx = d.rx;
      p.chan_offset_hz = chan.offset_hz;
      p.chan_usable_half_hz = chan.usable_half_hz;
      out.push_back(std::move(p));
    }
  }
  std::vector<Packet> acc;
  for (auto& p : out) {
    if (accept(p)) acc.push_back(std::move(p));
    else ++n_rejects_;
  }
  // _merge_grid
  std::vector<size_t> ord(acc.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(),
                   [&](size_t a, size_t b) { return std::fabs(acc[a].cfo_khz) < std::fabs(acc[b].cfo_khz); });
  std::vector<Packet> keep;
  for (size_t i : ord) {
    const Packet& p = acc[i];
    bool dup = false;
    for (const auto& k : keep)
      if (std::fabs(k.start_us - p.start_us) < 1.0 && std::fabs(k.center_mhz - p.center_mhz) < 2.5) {
        dup = true;
        break;
      }
    if (!dup) keep.push_back(p);
  }
  for (auto& p : keep) {
    const StdChannel* best = nullptr;
    double bd = 0;
    for (const auto& c : all_channels()) {
      const double dd = std::fabs(c.freq_hz / 1e6 - p.center_mhz);
      if (!best || dd < bd) {
        best = &c;
        bd = dd;
      }
    }
    p.off_grid = std::fabs(best->freq_hz / 1e6 - p.center_mhz) > 1.0;
    if (!p.off_grid) p.channel = best->name;
  }
  // _dedupe
  auto better = [](const Packet& a, const Packet& b) {  // score(a) > score(b)
    const bool ca = a.sig_crc_ok && *a.sig_crc_ok, cb = b.sig_crc_ok && *b.sig_crc_ok;
    if (ca != cb) return ca > cb;
    if (-a.lsig_evm_db != -b.lsig_evm_db) return -a.lsig_evm_db > -b.lsig_evm_db;
    return a.snr_db > b.snr_db;
  };
  std::stable_sort(keep.begin(), keep.end(), better);
  std::vector<Packet> ded;
  for (auto& p : keep) {
    bool dup = false;
    for (const auto& k : ded)
      if (std::fabs(k.start_us - p.start_us) < 1.0 && std::fabs(k.ch_freq_mhz - p.ch_freq_mhz) >= 20) {
        dup = true;
        break;
      }
    if (!dup) ded.push_back(std::move(p));
  }
  std::stable_sort(ded.begin(), ded.end(), [](const Packet& a, const Packet& b) {
    if (a.start_sample != b.start_sample) return a.start_sample < b.start_sample;
    return a.ch_freq_mhz < b.ch_freq_mhz;
  });
  return ded;
}

// ------------------------------------------------------------------ off-grid 802.11 (python offgrid_plan)
namespace {
constexpr double kOffgridStepHz = 5e6;
constexpr double kOffgridMaxSliceHz = 20e6;  // slice width: the largest power-of-two bin count within 20 MHz
constexpr double kOffgridRefLagsS[] = {2.1e-6, 2.3e-6, 2.5e-6, 2.7e-6, 3.7e-6, 3.9e-6, 4.1e-6, 4.3e-6};
constexpr double kOffgridRatio = 4.0;
constexpr double kOffgridMinCorr = 0.02;
constexpr size_t kOffgridMinRun = 3;
constexpr double kOffgridSpanHz = 5e6;
// windowed fallback (PLAN G13): a weak, low-duty link is diluted in the whole-capture correlation
constexpr double kOffgridWinS = 250e-6;
constexpr double kOffgridWinMinCorr = 0.06;
constexpr double kOffgridWinRatio = 6.0;
constexpr int64_t kOffgridWinMin = 2;
constexpr double kOffgridWinGate = 2.0;  // cheap trigger: whole-capture ratio

// python _offgrid_windows: >= kOffgridWinMin windows whose lag-3.2 us correlation passes
bool offgrid_windows(const cdvec& y, double ofs, int64_t L, const std::vector<int64_t>& refs) {
  const int64_t wn = static_cast<int64_t>(kOffgridWinS * ofs);
  int64_t maxlag = L + 1;
  for (int64_t r : refs) maxlag = std::max(maxlag, r);
  const int64_t nw = (static_cast<int64_t>(y.size()) - maxlag) / wn;
  if (nw < kOffgridWinMin) return false;
  std::vector<int64_t> lags = {L - 1, L, L + 1};
  lags.insert(lags.end(), refs.begin(), refs.end());
  int64_t passed = 0;
  dvec rs(lags.size());
  for (int64_t w = 0; w < nw; ++w) {
    const size_t a = static_cast<size_t>(w * wn), b = static_cast<size_t>((w + 1) * wn);
    double ew = 0;
    for (size_t i = a; i < b; ++i) ew += std::norm(y[i]);
    ew += 1e-30;
    for (size_t li = 0; li < lags.size(); ++li) {
      const size_t lag = static_cast<size_t>(lags[li]);
      cd acc(0, 0);
      for (size_t i = a; i < b; ++i) acc += std::conj(y[i]) * y[i + lag];
      rs[li] = std::abs(acc) / ew;
    }
    const double pk = std::max({rs[0], rs[1], rs[2]});
    const double rm = np::median(rs.data() + 3, rs.size() - 3);
    if (pk >= kOffgridWinMinCorr && pk >= kOffgridWinRatio * rm) ++passed;
  }
  return passed >= kOffgridWinMin;
}
}  // namespace

std::vector<Channel> scan_plan(double fc_hz, double fs_hz, double usable_bw_hz, double step_hz, double max_missing_hz) {
  const auto std_ch = channel_plan(fc_hz, fs_hz, usable_bw_hz, true, max_missing_hz);
  const double half = std::min(usable_bw_hz, fs_hz) / 2;
  const double lim = half + max_missing_hz - 8.4e6;
  const double k = std::floor(lim / step_hz);
  std::vector<Channel> out;
  for (double i = -k; i < k + 1; i += 1.0) {
    const double off = i * step_hz, f = fc_hz + off;
    const Channel* near = nullptr;
    for (const auto& c : std_ch)
      if (std::fabs(c.freq_hz - f) < step_hz / 2) {
        near = &c;
        break;
      }
    Channel c;
    c.name = near ? near->name : fmt("%.1fMHz", f / 1e6);
    c.ch = near ? near->ch : std::nullopt;
    c.freq_hz = f;
    c.offset_hz = off;
    c.usable_half_hz = half;
    c.partial = std::fabs(off) + 8.4e6 > half;
    c.standard = near != nullptr;
    out.push_back(c);
  }
  return out;
}

std::vector<Channel> offgrid_plan(const cfvec& x, double fs, double fc, const std::vector<Packet>& known,
                                  double usable_bw, const cfvec* Xin) {
  std::vector<Channel> plan;
  for (auto& c : scan_plan(fc, fs, usable_bw))
    if (!c.standard) plan.push_back(c);
  const size_t n = x.size();
  const double df = fs / static_cast<double>(n);
  const int64_t M = int64_t(1) << static_cast<int>(std::floor(std::log2(kOffgridMaxSliceHz / df)));
  if (plan.empty() || static_cast<double>(M) * df < 10e6) return {};
  cfvec Xown;
  if (!Xin) Xown = fft::np_fft(x);
  const cfvec& X = Xin ? *Xin : Xown;
  const int64_t K = M / 2;
  const double ofs = static_cast<double>(M) * df;
  const int64_t L = py_round_int(3.2e-6 * ofs);
  std::vector<int64_t> refs;
  for (double t : kOffgridRefLagsS) refs.push_back(py_round_int(t * ofs));
  double lim = 0;
  for (const auto& c : plan) lim = std::max(lim, std::fabs(c.offset_hz));
  const double nk = std::floor(lim / kOffgridStepHz);
  dvec offs;
  for (double i = -nk; i < nk + 1; i += 1.0) offs.push_back(i * kOffgridStepHz);
  std::vector<char> hit(offs.size(), 0);
  const int64_t nn = static_cast<int64_t>(n);
  parallel_for(offs.size(), [&](size_t si) {  // slices are independent: results do not depend on the threads
    const double off = offs[si];
    const int64_t k0 = py_round_int(off / df);
    cfvec sl(static_cast<size_t>(2 * K));
    for (int64_t j = 0; j < 2 * K; ++j) {  // ifftshift(X[idx]): idx = (k0 + arange(-K, K)) % n, rotated by K
      const int64_t src = (((k0 + (j + K) % (2 * K) - K) % nn) + nn) % nn;
      sl[static_cast<size_t>(j)] = X[static_cast<size_t>(src)];
    }
    const cfvec yf = fft::np_ifft(sl);
    cdvec y(yf.size());
    for (size_t i = 0; i < yf.size(); ++i) y[i] = cd(yf[i].real(), yf[i].imag());
    for (const auto& p : known)
      if (p.f_hi_mhz * 1e6 > fc + off - ofs / 2 && p.f_lo_mhz * 1e6 < fc + off + ofs / 2) {
        const int64_t a = static_cast<int64_t>(static_cast<double>(p.start_sample) / fs * ofs);
        const int64_t b = std::min<int64_t>(static_cast<int64_t>(static_cast<double>(p.end_sample) / fs * ofs) + 1,
                                            static_cast<int64_t>(y.size()));
        for (int64_t i = a; i < b; ++i) y[static_cast<size_t>(i)] = cd(0, 0);
      }
    double e = 0;
    for (const auto& v : y) e += std::norm(v);
    e += 1e-30;
    auto r = [&](int64_t lag) {
      cd acc(0, 0);
      for (size_t i = 0; i + static_cast<size_t>(lag) < y.size(); ++i) acc += std::conj(y[i]) * y[i + static_cast<size_t>(lag)];
      return std::abs(acc) / e;
    };
    const double pk = std::max({r(L - 1), r(L), r(L + 1)});
    dvec ref;
    for (int64_t lag : refs) ref.push_back(r(lag));
    const double rm = np::median(ref.data(), ref.size());
    hit[si] = (pk >= kOffgridMinCorr && pk >= kOffgridRatio * rm) ||
              (pk >= kOffgridWinGate * rm && offgrid_windows(y, ofs, L, refs));
  });
  dvec centres;
  for (size_t k = 0; k < offs.size();) {
    if (!hit[k]) {
      ++k;
      continue;
    }
    size_t j = k;
    while (j + 1 < offs.size() && hit[j + 1]) ++j;
    if (j - k + 1 >= kOffgridMinRun) centres.push_back(fc + 0.5 * (offs[k] + offs[j]));
    k = j + 1;
  }
  std::vector<Channel> out;
  for (const auto& c : plan)
    for (double m : centres)
      if (std::fabs(c.freq_hz - m) <= kOffgridSpanHz + 1) {
        out.push_back(c);
        break;
      }
  return out;
}

std::vector<Packet> merge_packets(std::vector<Packet> known, std::vector<Packet> extra, double tol_us) {
  std::stable_sort(extra.begin(), extra.end(), [](const Packet& a, const Packet& b) { return a.stf_metric > b.stf_metric; });
  std::vector<Packet> added;
  auto same = [&](const Packet& p, const Packet& q) {
    return std::fabs(p.start_us - q.start_us) <= tol_us && p.f_lo_mhz < q.f_hi_mhz && q.f_lo_mhz < p.f_hi_mhz;
  };
  for (auto& p : extra) {
    bool dup = std::any_of(known.begin(), known.end(), [&](const Packet& q) { return same(p, q); }) ||
               std::any_of(added.begin(), added.end(), [&](const Packet& q) { return same(p, q); });
    if (!dup) added.push_back(std::move(p));
  }
  std::stable_sort(added.begin(), added.end(), [](const Packet& a, const Packet& b) {
    if (a.start_sample != b.start_sample) return a.start_sample < b.start_sample;
    return a.ch_freq_mhz < b.ch_freq_mhz;
  });
  for (auto& p : added) known.push_back(std::move(p));
  return known;
}

}  // namespace si::rx::wifi::ofdm
