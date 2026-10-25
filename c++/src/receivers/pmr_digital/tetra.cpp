// python: receivers/pmr_digital/tetra.py (same operations and dtypes; complex products / divisions / exp as numpy).
#include "receivers/pmr_digital/tetra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>

#include "core/fft.hpp"
#include "core/npcompat.hpp"

namespace si::rx::pmr_digital::tetra {

namespace {

constexpr int8_t kYSync[38] = {1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
                               0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1};
constexpr int8_t kNTs[22] = {1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0};
constexpr int8_t kPTs[22] = {0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0};

// python rrc(beta, sps, span)
dvec rrc(double beta, int sps, int span = 8) {
  const int n = 2 * span * sps + 1;
  dvec h(n);
  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i - span * sps) / static_cast<double>(sps);
    const double x = (4 * beta) * t;
    const double den = 1 - x * x;
    if (std::fabs(den) < 1e-8) {
      h[i] = kPi / 4 * np::sinc(1 / (4 * beta));
    } else {
      h[i] = np::sinc(t) * std::cos((kPi * beta) * t) / den;
    }
  }
  dvec h2(n);
  for (int i = 0; i < n; ++i) h2[i] = h[i] * h[i];
  const double nrm = std::sqrt(np::sum(h2));
  for (double& v : h) v = v / nrm;
  return h;
}

const cdvec& rrc_c() {
  static const cdvec k = [] {
    const dvec h = rrc(0.35, kSps);
    cdvec c(h.size());
    for (size_t i = 0; i < h.size(); ++i) c[i] = cd(h[i], 0.0);
    return c;
  }();
  return k;
}

// np.convolve(y, _RRC, "same"): complex64 -> complex128
cdvec matched(const cfvec& y) {
  cdvec yd(y.size());
  for (size_t i = 0; i < y.size(); ++i) yd[i] = cd(y[i].real(), y[i].imag());
  return np::convolve(yd, rrc_c(), np::ConvMode::Same);
}

// numpy complex ** 4 (binary exponentiation in npy_cpow)
inline cd pow4(const cd& a) {
  const cd a2 = np::cmul(a, a);
  return np::cmul(a2, a2);
}

// np.exp(1j * theta) with a zero real part (libm cexp)
inline cd expj(double theta) { return std::exp(cd(0.0, theta)); }

inline double angle(const cd& z) { return std::atan2(z.imag(), z.real()); }

std::vector<int64_t> hits(const std::vector<int8_t>& bits, const int8_t* seq, int n, int max_err) {
  std::vector<int64_t> out;
  if (bits.size() < static_cast<size_t>(n)) return out;
  for (size_t i = 0; i + n <= bits.size(); ++i) {
    int e = 0;
    for (int j = 0; j < n && e <= max_err; ++j) e += bits[i + j] != seq[j];
    if (e <= max_err) out.push_back(static_cast<int64_t>(i));
  }
  return out;
}

// ------------------------------------------------------------------------------------------- BSCH codec
constexpr int kTaps[14] = {32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1};

std::vector<int8_t> scrambling_sequence(int n) {
  std::vector<int8_t> p(32, 0);
  p[0] = p[1] = 1;  // p(-31) = p(-30) = 1, p(-29..0) = e = 0
  std::vector<int8_t> out(n);
  for (int i = 0; i < n; ++i) {
    int b = 0;
    for (int t : kTaps) b ^= p[p.size() - t];
    p.push_back(static_cast<int8_t>(b));
    out[i] = static_cast<int8_t>(b);
  }
  return out;
}

struct Codec {
  std::vector<int8_t> scr;
  std::array<int, 120> deint{};
  std::array<int, 120> depunct{};
  double out[16][2][4]{};
  int nxt[16][2]{};
  Codec() : scr(scrambling_sequence(kBschBits)) {
    for (int i = 1; i <= 120; ++i) deint[i - 1] = (11 * i) % 120;
    static const int P23[3] = {1, 2, 5};
    for (int j = 1; j <= 120; ++j) {
      const int q = (j - 1) / 3;
      depunct[j - 1] = 8 * q + P23[j - 3 * q - 1] - 1;
    }
    static const int G[4][5] = {{1, 1, 0, 0, 1}, {1, 0, 1, 1, 1}, {1, 1, 1, 0, 1}, {1, 1, 0, 1, 1}};
    for (int s = 0; s < 16; ++s)
      for (int u = 0; u < 2; ++u) {
        const int reg[5] = {u, s & 1, (s >> 1) & 1, (s >> 2) & 1, (s >> 3) & 1};
        for (int i = 0; i < 4; ++i) {
          int acc = 0;
          for (int j = 0; j < 5; ++j) acc += G[i][j] * reg[j];
          out[s][u][i] = 2.0 * (acc % 2) - 1.0;
        }
        nxt[s][u] = ((s << 1) | u) & 0xF;
      }
  }
};

const Codec& codec() {
  static const Codec c;
  return c;
}

// python viterbi(v, n): 16-state rate-1/4 mother code, soft (+ = 1, 0 = erased), terminated in state 0
std::vector<int8_t> viterbi(const dvec& v, int n) {
  const Codec& C = codec();
  double pm[16];
  std::fill(pm, pm + 16, -1e9);
  pm[0] = 0.0;
  std::vector<std::array<int16_t, 16>> prev(n);
  std::vector<std::array<int8_t, 16>> pbit(n);
  for (int k = 0; k < n; ++k) {
    double m[16][2];
    for (int s = 0; s < 16; ++s)
      for (int u = 0; u < 2; ++u) {
        double acc = 0.0;
        for (int i = 0; i < 4; ++i) acc += C.out[s][u][i] * v[4 * k + i];
        m[s][u] = pm[s] + acc;
      }
    double nw[16];
    std::fill(nw, nw + 16, -1e9);
    for (int u = 0; u < 2; ++u)
      for (int s = 0; s < 16; ++s) {
        const int t = C.nxt[s][u];
        if (m[s][u] > nw[t]) {
          nw[t] = m[s][u];
          prev[k][t] = static_cast<int16_t>(s);
          pbit[k][t] = static_cast<int8_t>(u);
        }
      }
    std::copy(nw, nw + 16, pm);
  }
  std::vector<int8_t> bits(n);
  int s = 0;
  for (int k = n - 1; k >= 0; --k) {
    bits[k] = pbit[k][s];
    s = prev[k][s];
  }
  return bits;
}

int crc16(const int8_t* bits, int n) {
  int r = 0xFFFF;
  for (int i = 0; i < n; ++i) {
    const int fb = ((r >> 15) & 1) ^ bits[i];
    r = (r << 1) & 0xFFFF;
    if (fb) r ^= 0x1021;
  }
  return r ^ 0xFFFF;
}

int field(const int8_t* b, int& p, int n) {
  int v = 0;
  for (int i = 0; i < n; ++i) v = (v << 1) | b[p + i];
  p += n;
  return v;
}

}  // namespace

double symbol_line(const cfvec& y, double fs) {
  const cdvec z = matched(y);
  dvec e(z.size());
  for (size_t i = 0; i < z.size(); ++i) {
    const double a = np::cabs(z[i]);
    e[i] = a * a;
  }
  const double mu = np::mean(e);
  for (double& v : e) v = v - mu;
  const cdvec E = fft::rfft(e);
  const double val = 1.0 / (static_cast<double>(e.size()) * (1 / fs));
  const double w_on = 3 * fs / static_cast<double>(e.size()) + 5;
  double on = -1.0;
  dvec bg;
  for (size_t k = 0; k < E.size(); ++k) {
    const double f = static_cast<double>(k) * val;
    const double a = np::cabs(E[k]);
    if (std::fabs(f - kSr) <= w_on) on = std::max(on, a);
    if (std::fabs(f - kSr) > 1e3 && std::fabs(f - kSr) < 6e3) bg.push_back(a);
  }
  if (on < 0 || bg.empty()) return 0.0;
  return on / (np::median(bg) + 1e-30);
}

Demod demod(const cfvec& y, int quarter, int block_syms) {
  cdvec z = matched(y);
  const size_t L = z.size();
  Demod out;
  if (L <= static_cast<size_t>(kSps)) {
    out.err = kPi;
    return out;
  }
  // d = z[SPS:] * conj(z[:-SPS]); numpy elides the conj temporary (>= 256 KiB): conj(prev) * cur
  const size_t nd = L - kSps;
  const bool elide = nd * sizeof(cd) >= 256 * 1024;
  cdvec d4(nd);
  for (size_t i = 0; i < nd; ++i) {
    const cd d = elide ? np::cmul(std::conj(z[i]), z[i + kSps]) : np::cmul(z[i + kSps], std::conj(z[i]));
    d4[i] = pow4(d);
  }
  const double rot = (angle(np::mean(d4)) - kPi) / 4 + quarter * kPi / 2;
  const double b = -rot / kSps;
  for (size_t n = 0; n < L; ++n) z[n] = np::cmul(z[n], expj(b * static_cast<double>(n)));
  const int64_t blk = static_cast<int64_t>(block_syms) * kSps;
  cdvec dd;
  const int64_t stop = std::max<int64_t>(1, static_cast<int64_t>(L) - blk / 4);
  for (int64_t b0 = 0; b0 < stop; b0 += blk) {
    const int64_t e0 = std::min<int64_t>(static_cast<int64_t>(L), b0 + blk + kSps);
    if (b0 >= e0) continue;
    int ph = 0;
    double best = -1.0;
    for (int p = 0; p < kSps; ++p) {
      dvec a2;
      for (int64_t i = b0 + p; i < e0; i += kSps) {
        const double a = np::cabs(z[i]);
        a2.push_back(a * a);
      }
      const double mv = a2.empty() ? std::nan("") : np::mean(a2);
      if (p == 0 || mv > best) {  // python max(range(SPS), key=...): first maximum
        best = mv;
        ph = p;
      }
    }
    cdvec v;
    for (int64_t i = b0 + ph; i < e0 && static_cast<int64_t>(v.size()) < block_syms + 1; i += kSps) v.push_back(z[i]);
    if (v.size() < 2) continue;
    cdvec blkd(v.size() - 1), b4(v.size() - 1);
    for (size_t i = 0; i + 1 < v.size(); ++i) {
      blkd[i] = np::cmul(v[i + 1], std::conj(v[i]));
      b4[i] = pow4(blkd[i]);
    }
    const double r = (angle(np::mean(b4)) - kPi) / 4;
    const cd er = std::exp(cd(0.0, -r));
    for (const cd& x : blkd) dd.push_back(np::cmul(x, er));
  }
  const size_t m = dd.size();
  dvec ph(m), ab(m);
  for (size_t i = 0; i < m; ++i) {
    ph[i] = angle(dd[i]);
    ab[i] = np::cabs(dd[i]);
  }
  out.bits.resize(2 * m);
  out.soft.resize(2 * m);
  const double med = m ? np::median(ab) : 0.0;
  for (size_t i = 0; i < m; ++i) {
    out.bits[2 * i] = ph[i] < 0;
    out.bits[2 * i + 1] = std::fabs(ph[i]) > kPi / 2;
    const cd a = np::cmul(np::cdiv(dd[i], cd(ab[i] + 1e-12, 0.0)), cd(std::sqrt(ab[i] / (med + 1e-12)), 0.0));
    out.soft[2 * i] = -a.imag();
    out.soft[2 * i + 1] = -a.real();
  }
  if (m) {
    dvec e2(m);
    for (size_t i = 0; i < m; ++i) {
      const double t = angle(expj(4 * (ph[i] - kPi / 4)));
      e2[i] = t * t;
    }
    out.err = std::sqrt(np::mean(e2)) / 4;
  } else {
    out.err = kPi;
  }
  out.foff = rot * kSr / (2 * kPi);
  return out;
}

Training find_training(const std::vector<int8_t>& bits) {
  Training t;
  t.sync = hits(bits, kYSync, 38, 4);
  t.normal = hits(bits, kNTs, 22, 2);
  const auto p = hits(bits, kPTs, 22, 2);
  t.normal.insert(t.normal.end(), p.begin(), p.end());
  return t;
}

int64_t slot_consistent(const std::vector<int64_t>& pos) {
  if (pos.empty()) return 0;
  std::map<int64_t, int64_t> cnt;
  int64_t best = 0;
  for (int64_t p : pos) best = std::max(best, ++cnt[p % kSlotBits]);
  return best;
}

std::optional<Json> decode_bsch(const double* soft120) {
  const Codec& C = codec();
  double b4[120];
  for (int i = 0; i < 120; ++i) b4[i] = C.scr[i] == 1 ? -soft120[i] : soft120[i];
  dvec v(320, 0.0);
  for (int j = 0; j < 120; ++j) v[C.depunct[j]] = b4[C.deint[j]];
  const std::vector<int8_t> b2 = viterbi(v, 80);
  int crc = 0;
  for (int i = 60; i < 76; ++i) crc = (crc << 1) | b2[i];
  if (crc16(b2.data(), 60) != crc) return std::nullopt;
  int p = 0;
  Json f = Json::object();
  f["system_code"] = field(b2.data(), p, 4);
  f["colour_code"] = field(b2.data(), p, 6);
  f["timeslot"] = field(b2.data(), p, 2) + 1;  // TN coded 0..3 for timeslots 1..4
  f["frame"] = field(b2.data(), p, 5);
  f["multiframe"] = field(b2.data(), p, 6);
  f["sharing_mode"] = field(b2.data(), p, 2);
  f["ts_reserved_frames"] = field(b2.data(), p, 3);
  f["u_plane_dtx"] = field(b2.data(), p, 1);
  f["frame18_extension"] = field(b2.data(), p, 1);
  field(b2.data(), p, 1);  // reserved
  f["mcc"] = field(b2.data(), p, 10);
  f["mnc"] = field(b2.data(), p, 14);
  f["neighbour_cell_broadcast"] = field(b2.data(), p, 2);
  f["cell_service_level"] = field(b2.data(), p, 2);
  f["late_entry"] = field(b2.data(), p, 1);
  return f;
}

}  // namespace si::rx::pmr_digital::tetra
