// Port of python receivers/wifi/dsss.py (see dsss.hpp). Arithmetic follows numpy's: complex
// products through np::cmul, |z| through np::cabs, sums pairwise where numpy sums pairwise,
// np.round (rint) where Python rounds a numpy float, py_round where it rounds a Python float.
#include "receivers/wifi/dsss.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <tuple>

#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/signal.hpp"
#include "receivers/wifi/mac80211.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::dsss {

const double kBarker[11] = {1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1};

namespace {

using np::cabs;
using np::cmul;

constexpr uint16_t kSfdLong = 0xF3A0, kSfdShort = 0x05CF;
const int kDqpskBits[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // phase change k*pi/2 -> dibit

std::optional<double> rate_of(int64_t signal) {
  switch (signal) {
    case 0x0A: return 1.0;
    case 0x14: return 2.0;
    case 0x37: return 5.5;
    case 0x6E: return 11.0;
    default: return std::nullopt;
  }
}

inline cd conjd(cd z) { return {z.real(), -z.imag()}; }
// np.exp(-1j * x) for a real x: (cos(-x), sin(-x))
inline cd expmj(double x) { return {std::cos(-x), std::sin(-x)}; }
inline cd expj(double x) { return {std::cos(x), std::sin(x)}; }
inline double angle(cd z) { return std::atan2(z.imag(), z.real()); }
inline int64_t pymod(int64_t a, int64_t m) { return ((a % m) + m) % m; }

Bits sfd_bits(uint16_t v) {
  Bits b(16);
  for (int i = 0; i < 16; ++i) b[i] = (v >> i) & 1;
  return b;
}

// W_{-1}(z) and scipy's _calc_oa_lens block size, as core/signal.cpp computes them
double lambertw_m1(double z) {
  double w;
  if (z < -0.25) {
    const double p = -std::sqrt(2.0 * (std::exp(1.0) * z + 1.0));
    w = -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p;
  } else {
    const double l1 = std::log(-z), l2 = std::log(-l1);
    w = l1 - l2 + l2 / l1;
  }
  for (int i = 0; i < 200; ++i) {
    const double ew = std::exp(w), f = w * ew - z;
    const double wn = w - f / (ew * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0));
    if (std::fabs(wn - w) <= 1e-16 * std::fabs(wn)) {
      w = wn;
      break;
    }
    w = wn;
  }
  return w;
}

}  // namespace

namespace detail {
// scipy.signal.oaconvolve(x, k, "full") for a long complex64 x and a short float64 kernel:
// the same blocks, transforms and overlap-add as sig::oaconvolve (bit-identical), with the
// kernel spectrum computed once and the blocks transformed in parallel.
cdvec oaconvolve_cf_real(const cfvec& x, const dvec& k, int nthreads) {
  const size_t s1 = x.size(), s2 = k.size();
  if (s2 < 2 || s1 <= 2 * s2) return sig::oaconvolve(x, k, np::ConvMode::Full);
  const size_t overlap = s2 - 1;
  const double e = 2.718281828459045;
  const double opt = -static_cast<double>(overlap) * lambertw_m1(-1.0 / (2 * e * static_cast<double>(overlap)));
  const size_t blk = sig::next_fast_len(static_cast<size_t>(std::ceil(opt)), false);
  if (blk >= s1) return sig::oaconvolve(x, k, np::ConvMode::Full);
  const size_t step = blk - s2 + 1, shape_final = s1 + s2 - 1;
  size_t ns = (s1 + 1 + step - 1) / step;
  if ((blk - overlap) * ns < shape_final) ns += 1;
  // kernel spectrum: scipy.fft.fft of a real input (rfft + conjugate-symmetric fill)
  const cdvec h = fft::rfft(k.data(), s2, blk);
  cdvec K(blk);
  std::copy(h.begin(), h.end(), K.begin());
  for (size_t i = blk / 2 + 1; i < blk; ++i) K[i] = std::conj(K[blk - i]);
  const size_t keep = blk - overlap;
  cdvec out(ns * keep + overlap, cd(0, 0));
  std::vector<cdvec> tails(ns);
  parallel_for(ns, [&](size_t b) {
    const size_t a0 = b * step;
    const size_t len = a0 < s1 ? std::min(step, s1 - a0) : 0;
    const cfvec S = fft::fft(x.data() + a0, len, blk);
    cdvec P(blk);
    for (size_t i = 0; i < blk; ++i) P[i] = np::cmul(cd(S[i].real(), S[i].imag()), K[i]);
    fft::c2c(P.data(), P.data(), blk, false, 1.0 / static_cast<double>(blk));
    std::copy(P.begin(), P.begin() + keep, out.begin() + b * keep);
    tails[b].assign(P.begin() + keep, P.end());
  }, nthreads);
  for (size_t b = 1; b < ns; ++b)  // overlap-add, block by block as sig::oaconvolve
    for (size_t i = 0; i < overlap; ++i) out[b * keep + i] += tails[b - 1][i];
  out.resize(shape_final);
  return out;
}
}  // namespace detail

namespace {

// oaconvolve(x, kernel, "full")[len(kernel) - 1:] for a complex64 x and a real kernel
cdvec mf(const cfvec& x, const dvec& k) {
  cdvec c = detail::oaconvolve_cf_real(x, k, 1);  // callers run per channel in parallel
  c.erase(c.begin(), c.begin() + static_cast<ptrdiff_t>(k.size() - 1));
  return c;
}

dvec repeat_barker_rev(int n) {
  dvec v;
  for (int i = 10; i >= 0; --i)
    for (int k = 0; k < n; ++k) v.push_back(kBarker[i]);
  return v;
}

}  // namespace

// ------------------------------------------------------------------------------- bits
Bits descramble(const Bits& bits) {
  Bits out = bits;
  const size_t n = bits.size();
  for (size_t i = 4; i < n; ++i) out[i] ^= bits[i - 4];
  for (size_t i = 7; i < n; ++i) out[i] ^= bits[i - 7];
  for (size_t i = 0; i < std::min<size_t>(7, n); ++i) out[i] = bits[i];
  return out;
}

bool crc16_ok(const uint8_t* bits32, const uint8_t* crc_bits) {
  unsigned crc = 0xFFFF;
  for (int i = 0; i < 32; ++i) {
    const unsigned fb = ((crc >> 15) & 1) ^ bits32[i];
    crc = (crc << 1) & 0xFFFF;
    if (fb) crc ^= 0x1021;
  }
  crc ^= 0xFFFF;
  for (int i = 0; i < 16; ++i)
    if (((crc >> (15 - i)) & 1) != crc_bits[i]) return false;
  return true;
}

int64_t bits_to_int(const uint8_t* bits, size_t n) {
  int64_t v = 0;
  for (size_t i = 0; i < n; ++i) v |= static_cast<int64_t>(bits[i]) << i;
  return v;
}

std::vector<uint8_t> bits_to_bytes(const Bits& bits) {
  std::vector<uint8_t> out(bits.size() / 8);
  for (size_t i = 0; i < out.size(); ++i) {
    uint8_t b = 0;
    for (int k = 0; k < 8; ++k) b |= static_cast<uint8_t>(bits[8 * i + k] << k);
    out[i] = b;
  }
  return out;
}

int psdu_bytes(double rate, int64_t length_us, int64_t service) {
  int n = static_cast<int>(std::floor(length_us * rate / 8 + 1e-9));
  if (rate == 11.0 && (service & 0x80)) n -= 1;
  return n;
}

// ------------------------------------------------------------------------------- CCK
namespace {
CckBook make_book(int n_bits) {
  const double pi = kPi;
  const double qpsk[2][2] = {{0.0, pi / 2}, {pi, 3 * pi / 2}};  // CCK_QPSK[(b0, b1)]
  CckBook book;
  const int nd = n_bits == 6 ? 64 : 4;
  for (int d = 0; d < nd; ++d) {
    Bits b(n_bits);
    for (int i = 0; i < n_bits; ++i) b[i] = (d >> i) & 1;
    double p2, p3, p4;
    if (n_bits == 6) {
      p2 = qpsk[b[0]][b[1]];
      p3 = qpsk[b[2]][b[3]];
      p4 = qpsk[b[4]][b[5]];
    } else {
      p2 = b[0] * pi + pi / 2;
      p3 = 0.0;
      p4 = b[1] * pi;
    }
    const cd row[8] = {expj(p2 + p3 + p4), expj(p3 + p4), expj(p2 + p4), -expj(p4),
                       expj(p2 + p3),      expj(p3),      -expj(p2),     cd(1.0, 0.0)};
    book.C.insert(book.C.end(), row, row + 8);
    book.bits.push_back(b);
  }
  return book;
}
}  // namespace

const CckBook& cck11() {
  static const CckBook b = make_book(6);
  return b;
}
const CckBook& cck55() {
  static const CckBook b = make_book(2);
  return b;
}

// ------------------------------------------------------------------------------- detection
std::vector<Candidate> candidates(const cfvec& y2) {
  const dvec k = repeat_barker_rev(kSps);
  const cdvec c = mf(y2, k);
  const int64_t nsym = static_cast<int64_t>(c.size()) / kSym;
  if (nsym < kDetMinSym) return {};
  Mat<double> Pk(nsym, kSym);
  dvec E(nsym);
  {
    double row[kSym];
    for (int64_t m = 0; m < nsym; ++m) {
      for (int p = 0; p < kSym; ++p) {
        const double a = cabs(c[m * kSym + p]);
        Pk(m, p) = a * a;
        const double ay = static_cast<double>(cabs(y2[m * kSym + p]));
        row[p] = ay * ay;
      }
      E[m] = np::pairwise_sum(row, kSym);
    }
  }
  // Pk = cumsum(P, 0); Pk[K:] -= Pk[:-K]
  for (int64_t m = 1; m < nsym; ++m)
    for (int p = 0; p < kSym; ++p) Pk(m, p) += Pk(m - 1, p);
  for (int64_t m = nsym - 1; m >= kDetK; --m)
    for (int p = 0; p < kSym; ++p) Pk(m, p) -= Pk(m - kDetK, p);
  dvec Ek = np::cumsum(E);
  for (int64_t m = nsym - 1; m >= kDetK; --m) Ek[m] -= Ek[m - kDetK];
  for (int64_t m = 0; m + 1 < nsym; ++m) Ek[m] += E[m + 1];
  dvec rho(nsym);
  std::vector<int> best(nsym);
  std::vector<int64_t> idx;
  for (int64_t m = 0; m < nsym; ++m) {
    const double* r = Pk.row(m);
    const size_t b = np::argmax(r, kSym);
    const double peak = r[b];
    const double ratio = peak / (np::mean(r, kSym) + 1e-30);
    rho[m] = peak / (kSym * Ek[m] + 1e-30);
    best[m] = static_cast<int>(b);
    if (ratio > kDetRatio) idx.push_back(m);
  }
  if (idx.empty()) return {};
  std::vector<std::pair<int64_t, int64_t>> runs;
  int64_t start = idx[0], prev = idx[0];
  for (size_t q = 1; q < idx.size(); ++q) {
    if (idx[q] - prev > 3) {
      runs.push_back({start, prev});
      start = idx[q];
    }
    prev = idx[q];
  }
  runs.push_back({start, prev});
  std::vector<Candidate> out;
  for (const auto& [a, b] : runs) {
    std::vector<int64_t> cuts = {a};
    int ph = best[a], moved = 0;
    for (int64_t m = a + 1; m <= b; ++m) {
      const int d = pymod(best[m] - ph, kSym);
      if (std::min(d, kSym - d) > 3) {
        if (++moved >= 6) {  // a new timing has taken over: new PPDU
          cuts.push_back(m - 5);
          ph = best[m];
          moved = 0;
        }
      } else {
        moved = 0;
      }
    }
    cuts.push_back(b + 1);
    for (size_t q = 0; q + 1 < cuts.size(); ++q) {
      const int64_t u = cuts[q], v = cuts[q + 1];
      if (v - u >= kDetMinSym / 2 || (u == a && v - u >= 8)) {
        double mx = rho[u];
        for (int64_t m = u + 1; m < v; ++m) mx = std::max(mx, rho[m]);
        out.push_back({std::max<int64_t>(0, (u - kDetK) * kSym), v * kSym, mx});
      }
    }
  }
  return out;
}

// ------------------------------------------------------------------------------- tracking
Chan::Chan(cfvec y4, double ofs4, double fc_) : y(std::move(y4)), fc(fc_), fs(ofs4), T(ofs4 / 1e6) {
  c = mf(y, repeat_barker_rev(4));
  ac.resize(c.size());
  for (size_t i = 0; i < c.size(); ++i) ac[i] = cabs(c[i]);
  ym = mf(y, dvec(4, 0.25));
}

Chan::Track Chan::track(double tau, int64_t nsym, std::optional<double> T_, std::optional<double> A_) const {
  double Tl = T_ ? *T_ : T;
  const dvec& a = ac;
  const double lim = static_cast<double>(c.size()) - 6;
  double A;
  if (A_) {
    A = *A_;
  } else {
    const int64_t i = static_cast<int64_t>(tau);
    A = 0.0;
    if (i >= 0 && i < static_cast<int64_t>(a.size())) {
      const int64_t lo = std::max<int64_t>(0, i - 2), hi = std::min<int64_t>(a.size(), i + 3);
      A = *std::max_element(a.begin() + lo, a.begin() + hi);
    }
  }
  Track tr;
  tr.s.reserve(nsym > 0 ? nsym : 0);
  tr.taus.reserve(nsym > 0 ? nsym : 0);
  for (int64_t n = 0; n < nsym; ++n) {
    if (tau < 4 || tau > lim) break;
    int64_t i = static_cast<int64_t>(tau);
    double f = tau - i;
    double o = a[i] * (1 - f) + a[i + 1] * f;
    if (o < 0.6 * A) {
      const int64_t j = i - 3 + static_cast<int64_t>(np::argmax(a.data() + i - 3, 7));
      if (a[j] > 1.5 * o && a[j] > 0.5 * A) {  // off the main lobe: re-acquire the peak
        tau = static_cast<double>(j);
        i = j;
        f = 0.0;
        o = a[j];
      }
    }
    tr.taus.push_back(tau);
    const cd c0 = c[i], c1 = c[i + 1];
    tr.s.push_back(cd(c0.real() * (1 - f), c0.imag() * (1 - f)) + cd(c1.real() * f, c1.imag() * f));
    double err = 0.0;
    if (o > 0.4 * A) {
      const double e = a[i - 1] * (1 - f) + a[i] * f;
      const double lt = a[i + 1] * (1 - f) + a[i + 2] * f;
      err = std::max(-1.0, std::min(1.0, 3.0 * (lt - e) / (lt + e + 1e-12)));
      A = 0.98 * A + 0.02 * o;
    }
    Tl += 0.002 * err;
    tau += Tl + 0.3 * err;
  }
  tr.T = Tl;
  return tr;
}

// ------------------------------------------------------------------------------- demod helpers
Bits dbpsk_bits(const cdvec& s) {
  if (s.size() < 2) return {};
  const size_t n = s.size() - 1;
  cdvec d(n), d2(n);
  for (size_t k = 0; k < n; ++k) {
    d[k] = cmul(s[k + 1], conjd(s[k]));
    d2[k] = np::cdiv(cmul(d[k], d[k]), cd(cabs(d[k]) + 1e-30, 0.0));
  }
  // np.convolve(d2, ones(31), "same"), summed in index order
  const int64_t W = 31, h = (W - 1) / 2;
  Bits out(n);
  for (int64_t k = 0; k < static_cast<int64_t>(n); ++k) {
    const int64_t lo = std::max<int64_t>(0, k - h), hi = std::min<int64_t>(n - 1, k + W - 1 - h);
    double re = 0, im = 0;
    for (int64_t j = lo; j <= hi; ++j) {
      re += d2[j].real();
      im += d2[j].imag();
    }
    const double rot = angle(cd(re, im)) / 2;
    out[k] = std::fabs(angle(cmul(d[k], expmj(rot)))) > kPi / 2;
  }
  return out;
}

std::vector<std::pair<int64_t, std::string>> find_sfds(const Bits& bits) {
  std::vector<std::pair<int64_t, std::string>> out;
  const int64_t n = static_cast<int64_t>(bits.size());
  if (n < 40) return out;
  const int64_t nw = n - 15;
  std::vector<int> ones(nw);
  int acc = 0;
  for (int64_t i = 0; i < 16; ++i) acc += bits[i];
  for (int64_t k = 0; k < nw; ++k) {
    ones[k] = acc;
    if (k + 16 < n) acc += bits[k + 16] - bits[k];
  }
  const std::pair<uint16_t, const char*> sfds[] = {{kSfdLong, "long"}, {kSfdShort, "short"}};
  for (const auto& [v, pre] : sfds) {
    const Bits sfd = sfd_bits(v);
    std::vector<int> err(nw);
    for (int64_t k = 0; k < nw; ++k) {
      int e = 0;
      for (int j = 0; j < 16; ++j) e += bits[k + j] != sfd[j];
      err[k] = e;
    }
    const bool is_long = std::string(pre) == "long";
    for (int64_t k = 16; k < nw; ++k) {
      if (err[k] > 1) continue;
      const int run = ones[k - 16];
      if ((is_long && run >= 13) || (!is_long && run <= 3)) {
        if (err[k] == 0 || (k + 1 < nw && err[k + 1] > 0 && err[k - 1] > 0)) out.push_back({k + 16, pre});
      }
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

std::pair<Bits, int64_t> rebuild_scrambled(const Bits& raw, const Bits& desc) {
  const size_t n = raw.size();
  Bits best_bits;
  int64_t best = -1;
  Bits out(n);
  for (int st = 0; st < 128; ++st) {
    uint8_t hist[7];  // hist[0] = s_{-7} .. hist[6] = s_{-1}
    for (int k = 0; k < 7; ++k) hist[6 - k] = (st >> k) & 1;
    int64_t agree = 0;
    for (size_t i = 0; i < n; ++i) {
      const uint8_t s = desc[i] ^ hist[3] ^ hist[0];
      out[i] = s;
      agree += s == raw[i];
      for (int k = 0; k < 6; ++k) hist[k] = hist[k + 1];
      hist[6] = s;
    }
    if (agree > best) {
      best = agree;
      best_bits = out;
    }
  }
  return {best_bits, best};
}

namespace {

struct Coherent {
  std::vector<int> q;
  double ph;
  cdvec r;
};

// python _coherent: phase-tracked decisions on a k*grid constellation
Coherent coherent(const cdvec& z, double grid, double ph0, double gain = 0.08) {
  const int step = static_cast<int>(py_round_int(grid / (kPi / 2)));
  Coherent out;
  out.q.resize(z.size());
  out.r.resize(z.size());
  double ph = ph0;
  for (size_t k = 0; k < z.size(); ++k) {
    const cd v = z[k];
    const double a = angle(cmul(v, expmj(ph)));
    const int64_t qq = static_cast<int64_t>(std::nearbyint(a / grid)) * step;
    out.q[k] = static_cast<int>(pymod(qq, 4));
    out.r[k] = cmul(v, expmj(ph + qq * kPi / 2));
    ph += gain * angle(out.r[k]);
  }
  out.ph = ph;
  return out;
}

std::vector<int64_t> dibit_flips(int64_t dq_old, int64_t dq_new, int64_t at) {
  const int* a = kDqpskBits[pymod(dq_old, 4)];
  const int* b = kDqpskBits[pymod(dq_new, 4)];
  std::vector<int64_t> out;
  for (int i = 0; i < 2; ++i)
    if (a[i] != b[i]) out.push_back(at + i);
  return out;
}

struct Alt {
  double margin;
  std::vector<int64_t> flips;
};

std::vector<Alt> qpsk_alternatives(const std::vector<int>& q, int q_prev, const cdvec& r, int rate_bits) {
  std::vector<Alt> alts;
  const size_t n = q.size();
  std::vector<int64_t> qq(n + 1);
  qq[0] = q_prev;
  for (size_t k = 0; k < n; ++k) qq[k + 1] = q[k];
  const double s2 = std::sqrt(2.0);
  for (size_t k = 0; k < n; ++k) {
    const cd v = r[k];
    if (rate_bits == 1) {
      std::vector<int64_t> fl = {static_cast<int64_t>(k)};
      if (k + 1 < n) fl.push_back(k + 1);
      alts.push_back({std::fabs(v.real()), fl});
      continue;
    }
    const std::pair<int, double> opts[2] = {{1, (v.real() - v.imag()) / s2}, {-1, (v.real() + v.imag()) / s2}};
    for (const auto& [d, m] : opts) {
      const int64_t dq0 = qq[k + 1] - qq[k];
      auto fl = dibit_flips(dq0, dq0 + d, 2 * k);
      if (k + 1 < n) {
        const int64_t dq1 = qq[k + 2] - qq[k + 1];
        auto f2 = dibit_flips(dq1, dq1 - d, 2 * k + 2);
        fl.insert(fl.end(), f2.begin(), f2.end());
      }
      alts.push_back({std::fabs(m), fl});
    }
  }
  return alts;
}

uint32_t le32(const std::vector<uint8_t>& b, size_t at) {
  return static_cast<uint32_t>(b[at]) | (static_cast<uint32_t>(b[at + 1]) << 8) |
         (static_cast<uint32_t>(b[at + 2]) << 16) | (static_cast<uint32_t>(b[at + 3]) << 24);
}

// CRC-32 of the MPDU body XOR its FCS field (0 when the FCS passes), and the bytes
std::pair<uint32_t, std::vector<uint8_t>> fcs_syndrome(const Bits& d) {
  auto b = bits_to_bytes(d);
  const size_t n = b.size() - 4;
  return {mac80211::crc32(b.data(), n) ^ le32(b, n), b};
}

// python _fcs_correct: <= kFixMaxWeight of the least reliable decisions whose flips pass the FCS
std::optional<std::pair<Bits, int>> fcs_correct(const Bits& d, const std::vector<Alt>& alts) {
  const size_t nbits = d.size();
  const uint32_t S0 = fcs_syndrome(d).first;
  if (S0 == 0) return std::make_pair(d, 0);
  std::vector<uint8_t> zeros(nbits / 8 - 4, 0);
  const uint32_t zero = mac80211::crc32(zeros.data(), zeros.size());
  std::vector<size_t> order(alts.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return alts[a].margin < alts[b].margin; });
  if (order.size() > kFixCandidates) order.resize(kFixCandidates);
  const size_t nc = order.size();
  std::vector<Bits> errs;
  std::vector<uint32_t> syns;
  for (size_t ci : order) {
    Bits e(nbits + 8, 0);
    for (int64_t p : alts[ci].flips)
      for (int64_t q : {p, p + 4, p + 7}) e[q] ^= 1;
    e.resize(nbits);
    auto b = bits_to_bytes(e);
    const size_t n = b.size() - 4;
    syns.push_back(mac80211::crc32(b.data(), n) ^ zero ^ le32(b, n));
    errs.push_back(std::move(e));
  }
  std::vector<std::vector<size_t>> groups;
  for (size_t i = 0; i < nc; ++i) groups.push_back({i});
  for (size_t i = 0; i < nc; ++i)
    for (size_t j = i + 1; j < nc; ++j) groups.push_back({i, j});
  std::map<uint32_t, size_t> table;  // syndrome -> first group
  auto syn_of = [&](const std::vector<size_t>& g) {
    uint32_t s = 0;
    for (size_t i : g) s ^= syns[i];
    return s;
  };
  for (size_t gi = 0; gi < groups.size(); ++gi) table.emplace(syn_of(groups[gi]), gi);
  std::optional<std::set<size_t>> best;
  for (const auto& g : groups) {
    const uint32_t sg = syn_of(g);
    std::vector<std::vector<size_t>> hs;
    if (sg == S0) hs.push_back({});
    auto it = table.find(S0 ^ sg);
    if (it != table.end()) hs.push_back(groups[it->second]);
    for (const auto& h : hs) {
      std::set<size_t> combo(g.begin(), g.end());
      combo.insert(h.begin(), h.end());
      if (combo.size() == g.size() + h.size() && combo.size() <= kFixMaxWeight)
        if (!best || combo.size() < best->size()) best = combo;
    }
  }
  if (!best) return std::nullopt;
  Bits out = d;
  for (size_t i : *best)
    for (size_t k = 0; k < nbits; ++k) out[k] ^= errs[i][k];
  if (fcs_syndrome(out).first != 0) return std::nullopt;
  return std::make_pair(out, static_cast<int>(best->size()));
}

Bits q_to_bits(const std::vector<int>& q, int q_prev, int rate_bits) {
  Bits out;
  int prev = q_prev;
  for (int v : q) {
    const int dq = static_cast<int>(pymod(v - prev, 4));
    prev = v;
    if (rate_bits == 1) {
      out.push_back(dq == 2);
    } else {
      out.push_back(kDqpskBits[dq][0]);
      out.push_back(kDqpskBits[dq][1]);
    }
  }
  return out;
}

}  // namespace

Bits diff_bits(const cdvec& s, double rot, int rate_bits) {
  Bits out;
  const cd er = expmj(rot);
  for (size_t k = 0; k + 1 < s.size(); ++k) {
    const double d = angle(cmul(cmul(s[k + 1], conjd(s[k])), er));
    if (rate_bits == 1) {
      out.push_back(std::fabs(d) > kPi / 2);
    } else {
      const int q = static_cast<int>(pymod(static_cast<int64_t>(std::nearbyint(d / (kPi / 2))), 4));
      out.push_back(kDqpskBits[q][0]);
      out.push_back(kDqpskBits[q][1]);
    }
  }
  return out;
}

// ------------------------------------------------------------------------------- PPDU decode
namespace {

struct PsduOut {
  std::optional<std::vector<uint8_t>> mpdu;
  std::string how;
  int nfix = 0;
};

class Ppdu {
 public:
  Ppdu(const Chan& ch, const cdvec& s, const dvec& taus, double T, int64_t e, const std::string& pre, double amp)
      : ch_(ch), s_(s), taus_(taus), T_(T), e_(e), pre_(pre), amp_(amp) {}

  std::optional<Json> run();

 private:
  PsduOut psdu(double rate, int nbytes, int64_t k_psdu);
  struct Got {
    Bits bits;
    std::vector<Alt> alts;
  };
  std::optional<Got> method(int which, double rate, size_t nbits, int64_t k_psdu);
  std::optional<cdvec> barker_syms(double rate, size_t nbits, int64_t k_psdu);
  std::optional<Got> barker_coherent(double rate, size_t nbits, int64_t k_psdu);
  std::optional<Got> barker_diff(double rate, size_t nbits, int64_t k_psdu);
  std::optional<Got> cck(double rate, size_t nbits, int64_t k_psdu, int ntaps, bool fit_clock);

  const Chan& ch_;
  const cdvec& s_;
  const dvec& taus_;
  double T_;
  int64_t e_;
  std::string pre_;
  double amp_;
  double rot_ = 0, th0_ = 0, ph_ = 0;
  int64_t p0_ = 0;
  std::vector<int> q_pre_, q_hdr_;
  Bits prefix_;
};

std::optional<Json> Ppdu::run() {
  const bool lng = pre_ == "long";
  const int64_t npre = lng ? 144 : 72;
  const int64_t e = e_;
  const int64_t p0 = std::max<int64_t>(0, e - npre + 1);
  Bits desc_full(npre - 16, lng ? 1 : 0);
  const Bits sfd = sfd_bits(lng ? kSfdLong : kSfdShort);
  desc_full.insert(desc_full.end(), sfd.begin(), sfd.end());
  const int64_t ofs = npre - (e - p0);
  Bits desc(desc_full.begin() + std::min<int64_t>(std::max<int64_t>(ofs, 0), npre), desc_full.end());
  const int64_t hi = std::min<int64_t>(e + 1, s_.size());
  const cdvec seg(s_.begin() + p0, s_.begin() + std::max(hi, p0));
  const Bits raw = dbpsk_bits(seg);
  if (static_cast<int64_t>(raw.size()) < e - p0 || e - p0 < 24) return std::nullopt;
  const auto rebuilt = rebuild_scrambled(Bits(raw.begin(), raw.begin() + (e - p0)), desc);
  const Bits& sc = rebuilt.first;
  const int64_t agree = rebuilt.second;
  // data-aided reference: symbol p0 has phase 0; symbol n+1 = symbol n + pi * scrambled bit n
  const size_t nz = seg.size();
  dvec a(nz);
  a[0] = 1.0;
  for (size_t i = 1; i < nz; ++i) a[i] = a[i - 1] * (1.0 - 2.0 * sc[i - 1]);
  cdvec z(nz);
  for (size_t i = 0; i < nz; ++i) z[i] = cmul(seg[i], cd(a[i], 0.0));
  cdvec prod(nz - 1);
  for (size_t i = 0; i + 1 < nz; ++i) prod[i] = cmul(z[i + 1], conjd(z[i]));
  const double rot = angle(np::pairwise_sum(prod.data(), prod.size()));
  cdvec zr(nz);
  for (size_t i = 0; i < nz; ++i) zr[i] = cmul(z[i], expmj(rot * static_cast<double>(i)));
  const cd mu = np::mean(zr);
  const double th0 = angle(mu);
  const double amu = cabs(mu);
  const double snr_sym = amu * amu / (np::var(zr) + 1e-30);
  rot_ = rot;
  th0_ = th0;
  p0_ = p0;
  q_pre_.resize(nz);
  for (size_t i = 0; i < nz; ++i) q_pre_[i] = static_cast<int>(pymod(static_cast<int64_t>(1 - a[i]), 4));
  Json rec = Json::object();
  rec["preamble"] = pre_;
  rec["sync_bits_agree"] = py_round(static_cast<double>(agree) / static_cast<double>(e - p0), 3);
  rec["cfo_hz"] = py_round_int(rot / (2 * kPi) * 1e6);
  rec["snr_db"] = np::round(10 * std::log10(snr_sym + 1e-12) - 10.4, 1);
  rec["start_s"] = (taus_[e] - (npre - 1) * T_) / ch_.fs;
  rec["snr_db_symbol"] = np::round(10 * std::log10(snr_sym + 1e-12), 1);
  // ---- PLCP header
  const int64_t nh = lng ? 48 : 24;
  if (static_cast<int64_t>(s_.size()) < e + 1 + nh) return std::nullopt;
  cdvec hs(s_.begin() + e + 1, s_.begin() + e + 1 + nh);
  cdvec zh(nh);
  for (int64_t k = 0; k < nh; ++k) {
    const double kk = static_cast<double>(e + 1 - p0 + k);
    zh[k] = cmul(hs[k], expmj(th0 + rot * kk));
  }
  const int hb_bits = lng ? 1 : 2;
  const int q_last = q_pre_.back();
  const Coherent ch0 = coherent(zh, hb_bits == 1 ? kPi : kPi / 2, 0.0);
  struct Try {
    const char* how;
    Bits hraw;
    bool has_q;
  };
  cdvec ds(1, s_[e]);
  ds.insert(ds.end(), hs.begin(), hs.end());
  const Try tries[2] = {{"coherent", q_to_bits(ch0.q, q_last, hb_bits), true},
                        {"differential", diff_bits(ds, rot, hb_bits), false}};
  const Try* ok = nullptr;
  Bits hb;
  auto header_bits = [&](const Bits& hraw) {
    Bits all = sc;
    all.insert(all.end(), hraw.begin(), hraw.end());
    Bits d = descramble(all);
    Bits out(48, 0);
    for (size_t i = 0; i < 48 && sc.size() + i < d.size(); ++i) out[i] = d[sc.size() + i];
    return out;
  };
  for (const auto& t : tries) {
    hb = header_bits(t.hraw);
    if (crc16_ok(hb.data(), hb.data() + 32)) {
      ok = &t;
      break;
    }
  }
  auto rate_json = [](std::optional<double> r) { return r ? Json(*r) : Json(nullptr); };
  if (!ok) {
    hb = header_bits(tries[0].hraw);
    const int64_t sigv = bits_to_int(hb.data(), 8);
    rec["plcp_crc_ok"] = false;
    rec["signal"] = sigv;
    rec["rate_mbps"] = rate_json(rate_of(sigv));
    rec["service"] = bits_to_int(hb.data() + 8, 8);
    rec["length_us"] = bits_to_int(hb.data() + 16, 16);
    return rec;
  }
  const int64_t signal = bits_to_int(hb.data(), 8), service = bits_to_int(hb.data() + 8, 8),
                length = bits_to_int(hb.data() + 16, 16);
  const auto rate = rate_of(signal);
  rec["plcp_crc_ok"] = true;
  rec["signal"] = signal;
  rec["rate_mbps"] = rate_json(rate);
  rec["service"] = service;
  rec["length_us"] = length;
  rec["locked_clocks"] = static_cast<bool>(service & 0x04);
  rec["header_demod"] = ok->how;
  if (!rate) return rec;
  if (ok->has_q) {
    q_hdr_ = ch0.q;
    ph_ = ch0.ph;
  } else {  // differential header: rebuild quadrants
    q_hdr_.resize(nh);
    int64_t acc = q_last;
    for (int64_t i = 0; i < nh; ++i) {
      int64_t dq;
      if (hb_bits == 2) {
        const int b0 = ok->hraw[2 * i], b1 = ok->hraw[2 * i + 1];
        dq = 0;
        for (int v = 0; v < 4; ++v)
          if (kDqpskBits[v][0] == b0 && kDqpskBits[v][1] == b1) {
            dq = v;
            break;
          }
      } else {
        dq = 2 * ok->hraw[i];
      }
      acc += dq;
      q_hdr_[i] = static_cast<int>(pymod(acc, 4));
    }
    ph_ = 0.0;
  }
  prefix_ = sc;
  prefix_.insert(prefix_.end(), ok->hraw.begin(), ok->hraw.end());
  const int nbytes = psdu_bytes(*rate, length, service);
  rec["psdu_bytes"] = nbytes;
  if (nbytes < 10 || nbytes > 2400) return rec;
  PsduOut p = psdu(*rate, nbytes, e + 1 + nh);
  if (!p.mpdu) {
    rec["truncated"] = true;
    return rec;
  }
  const bool fok = mac80211::fcs_ok(*p.mpdu);
  rec["fcs_ok"] = fok;
  rec["psdu_demod"] = p.how;
  rec["symbols_corrected"] = p.nfix;
  if (fok) {
    const size_t n = p.mpdu->size() - 4;
    rec["mac"] = mac80211::parse(p.mpdu->data(), n);
    rec["mpdu_hex"] = wifi::detail::hex(wifi::detail::Bytes(p.mpdu->data(), std::min<size_t>(n, 512)));
  }
  return rec;
}

PsduOut Ppdu::psdu(double rate, int nbytes, int64_t k_psdu) {
  const size_t nbits = static_cast<size_t>(nbytes) * 8;
  static const char* kNames[2][3] = {{"coherent Barker", "differential Barker", nullptr},
                                     {"CCK, LS equaliser", "CCK, 1-tap", "CCK, LS equaliser, fitted clock"}};
  const int fam = (rate == 1.0 || rate == 2.0) ? 0 : 1;
  struct Tried {
    std::string how;
    Bits d;
    std::vector<Alt> alts;
    std::vector<uint8_t> mpdu;
  };
  std::vector<Tried> tried;
  for (int w = 0; w < 3 && kNames[fam][w]; ++w) {
    auto got = method(fam * 2 + w, rate, nbits, k_psdu);
    if (!got) continue;
    Bits all = prefix_;
    all.insert(all.end(), got->bits.begin(), got->bits.end());
    Bits dd = descramble(all);
    Bits d(dd.begin() + prefix_.size(), dd.begin() + std::min(prefix_.size() + nbits, dd.size()));
    auto [S, mpdu] = fcs_syndrome(d);
    if (S == 0) return {mpdu, kNames[fam][w], 0};
    tried.push_back({kNames[fam][w], std::move(d), std::move(got->alts), std::move(mpdu)});
  }
  for (const auto& t : tried) {
    if (t.alts.empty()) continue;
    auto fixed = fcs_correct(t.d, t.alts);
    if (fixed)
      return {bits_to_bytes(fixed->first), t.how + ", FCS-guided correction of " + std::to_string(fixed->second) +
                                               " symbol(s)",
              fixed->second};
  }
  if (tried.empty()) return {std::nullopt, "", 0};
  return {tried[0].mpdu, tried[0].how, 0};
}

std::optional<Ppdu::Got> Ppdu::method(int which, double rate, size_t nbits, int64_t k_psdu) {
  switch (which) {
    case 0: return barker_coherent(rate, nbits, k_psdu);
    case 1: return barker_diff(rate, nbits, k_psdu);
    case 2: return cck(rate, nbits, k_psdu, kEqTaps, false);
    case 3: return cck(rate, nbits, k_psdu, 1, false);
    default: return cck(rate, nbits, k_psdu, kEqTaps, true);
  }
}

std::optional<cdvec> Ppdu::barker_syms(double rate, size_t nbits, int64_t k_psdu) {
  const int64_t nsym = static_cast<int64_t>(nbits / (rate == 1.0 ? 1 : 2));
  cdvec s = s_;
  const int64_t have = static_cast<int64_t>(s.size()) - k_psdu;
  if (have < nsym) {  // continue tracking past the run
    auto more = ch_.track(taus_.back() + T_, nsym - have + 1, T_, amp_);
    s.insert(s.end(), more.s.begin(), more.s.end());
  }
  const int64_t a = std::min<int64_t>(k_psdu, s.size()), b = std::min<int64_t>(k_psdu + nsym, s.size());
  if (b - a != nsym) return std::nullopt;
  return cdvec(s.begin() + a, s.begin() + b);
}

std::optional<Ppdu::Got> Ppdu::barker_coherent(double rate, size_t nbits, int64_t k_psdu) {
  auto seg = barker_syms(rate, nbits, k_psdu);
  if (!seg) return std::nullopt;
  cdvec z(seg->size());
  for (size_t k = 0; k < z.size(); ++k) {
    const double kk = static_cast<double>(k_psdu - p0_ + static_cast<int64_t>(k));
    z[k] = cmul((*seg)[k], expmj(th0_ + rot_ * kk));
  }
  const int nb = rate == 1.0 ? 1 : 2;
  const Coherent c = coherent(z, rate == 1.0 ? kPi : kPi / 2, ph_);
  return Got{q_to_bits(c.q, q_hdr_.back(), nb), qpsk_alternatives(c.q, q_hdr_.back(), c.r, nb)};
}

std::optional<Ppdu::Got> Ppdu::barker_diff(double rate, size_t nbits, int64_t k_psdu) {
  auto seg = barker_syms(rate, nbits, k_psdu);
  if (!seg) return std::nullopt;
  cdvec s(1, s_[k_psdu - 1]);
  s.insert(s.end(), seg->begin(), seg->end());
  return Got{diff_bits(s, rot_, rate == 1.0 ? 1 : 2), {}};
}

std::optional<Ppdu::Got> Ppdu::cck(double rate, size_t nbits, int64_t k_psdu, int ntaps, bool fit_clock) {
  const Chan& ch = ch_;
  const size_t bps = rate == 11.0 ? 8 : 4;
  const size_t nsym = nbits / bps;
  const int64_t p0 = p0_;
  const int64_t nk = k_psdu - p0;
  if (nk < 24) return std::nullopt;
  dvec ks(nk), taus(nk);
  for (int64_t i = 0; i < nk; ++i) {
    ks[i] = static_cast<double>(p0 + i);
    taus[i] = taus_[p0 + i];
  }
  // symbol clock: locked clocks -> the chip rate follows the carrier offset
  const double T = fit_clock ? np::polyfit(ks, taus, 1)[0] : ch.T / (1 + rot_ / (2 * kPi) * 1e6 / ch.fc);
  dvec dev(nk);
  for (int64_t i = 0; i < nk; ++i) dev[i] = taus[i] - T * static_cast<double>(i);
  const double tau0 = np::mean(dev);
  const double Tc = T / 11;
  const double w_rad = rot_ / T;
  // training chips: preamble (+-1) and header (quadrants) symbols, known after the CRC
  std::vector<int> qs(q_pre_.begin(), q_pre_.begin() + (nk - static_cast<int64_t>(q_hdr_.size())));
  qs.insert(qs.end(), q_hdr_.begin(), q_hdr_.end());
  const size_t ntr = qs.size() * 11;
  cdvec known(ntr);
  for (size_t i = 0; i < qs.size(); ++i) {
    const cd e = expj(kPi / 2 * qs[i]);
    for (int j = 0; j < 11; ++j) known[i * 11 + j] = cd(e.real() * kBarker[j], e.imag() * kBarker[j]);
  }
  const double t0_psdu = tau0 + static_cast<double>(nk) * T;
  const size_t nps = nsym * 8;
  if (nps == 0) return std::nullopt;
  if (t0_psdu + (static_cast<double>(nps - 1) + 0.375) * Tc + 8 * Tc >= static_cast<double>(ch.ym.size()))
    return std::nullopt;
  const int half = ntaps / 2;
  const double xmax = static_cast<double>(ch.ym.size()) - 1.001;
  // python chips_at(t): interpolated chip-MF samples at t + Tc*m, derotated
  auto chips_at = [&](double t, cd* row) {
    for (int m = -half; m <= half; ++m) {
      const double tt = t + Tc * m;
      double x = std::min(std::max(tt - 1.5, 0.0), xmax);
      const int64_t i = static_cast<int64_t>(x);
      const double f = x - static_cast<double>(i);
      const cd y0 = ch.ym[i], y1 = ch.ym[i + 1];
      const cd v(y0.real() * (1 - f) + y1.real() * f, y0.imag() * (1 - f) + y1.imag() * f);
      row[m + half] = cmul(v, expmj(w_rad * tt + th0_ - w_rad * tau0));
    }
  };
  Mat<cd> A(ntr, ntaps);
  for (size_t i = 0; i < ntr; ++i) chips_at(tau0 + (static_cast<double>(i) + 0.375) * Tc, A.row(i));
  double lam = 0;
  for (size_t i = 0; i < ntr; ++i) {
    const cd v = A(i, half);
    lam += v.real() * v.real() + v.imag() * v.imag();
  }
  lam *= 1e-3;
  Mat<cd> AhA(ntaps, ntaps, cd(0, 0));
  cdvec Ahb(ntaps, cd(0, 0));
  for (int r = 0; r < ntaps; ++r) {
    for (int c = 0; c < ntaps; ++c) {
      cd acc(0, 0);
      for (size_t i = 0; i < ntr; ++i) acc += cmul(conjd(A(i, r)), A(i, c));
      AhA(r, c) = acc + (r == c ? cd(lam, 0) : cd(0, 0));
    }
    cd acc(0, 0);
    for (size_t i = 0; i < ntr; ++i) acc += cmul(conjd(A(i, r)), known[i]);
    Ahb[r] = acc;
  }
  cdvec w;
  try {
    w = np::solve(AhA, Ahb);
  } catch (const std::exception&) {
    return std::nullopt;  // numpy raises LinAlgError (singular): Python would fail this PPDU too
  }
  const CckBook& book = rate == 11.0 ? cck11() : cck55();
  const size_t nc = book.ncode();
  std::vector<cd> rot4(4);
  for (int q = 0; q < 4; ++q) rot4[q] = expj(-(kPi / 2) * q);
  Bits out(nsym * bps, 0);
  std::vector<int> qsym(nsym);
  struct Second {
    double margin;
    size_t m2;
    int q2;
  };
  std::vector<Second> second(nsym);
  const int q0 = q_hdr_.back();
  int q_prev = q0;
  double ph = 0.0;
  cd row[16];
  cdvec E(8), corr(nc);
  std::vector<double> M(nc * 4);
  for (size_t i = 0; i < nsym; ++i) {
    for (int j = 0; j < 8; ++j) {
      chips_at(t0_psdu + (static_cast<double>(i * 8 + j) + 0.375) * Tc, row);
      cd acc(0, 0);
      for (int t = 0; t < ntaps; ++t) acc += cmul(row[t], w[t]);
      E[j] = acc;
    }
    const cd eph = expmj(ph);
    for (size_t m = 0; m < nc; ++m) {
      cd acc(0, 0);
      for (int j = 0; j < 8; ++j) acc += cmul(E[j], conjd(book.C[m * 8 + j]));
      corr[m] = cd(acc.real() / 8, acc.imag() / 8);
      const cd v = cmul(corr[m], eph);
      for (int q = 0; q < 4; ++q) M[m * 4 + q] = cmul(v, rot4[q]).real();
    }
    // top two of argsort(M) (numpy's quicksort order matters only for exact ties)
    size_t b1 = 0;
    for (size_t k = 1; k < M.size(); ++k)
      if (M[k] > M[b1]) b1 = k;
    size_t b2 = b1 == 0 ? 1 : 0;
    for (size_t k = 0; k < M.size(); ++k)
      if (k != b1 && M[k] > M[b2]) b2 = k;
    bool tie = false;
    for (size_t k = 0; k < M.size() && !tie; ++k)
      tie = (k != b1 && M[k] == M[b1]) || (k != b1 && k != b2 && M[k] == M[b2]);
    if (tie) {
      auto order = np::argsort_quick(M);
      b1 = static_cast<size_t>(order[order.size() - 1]);
      b2 = static_cast<size_t>(order[order.size() - 2]);
    }
    const size_t m = b1 / 4, m2 = b2 / 4;
    const int q = static_cast<int>(b1 % 4), q2 = static_cast<int>(b2 % 4);
    second[i] = {M[b1] - M[b2], m2, q2};
    const cd vm = cmul(corr[m], eph);
    const double a = angle(cmul(vm, expmj(q * kPi / 2)));
    ph += 0.05 * a;
    const int dq = static_cast<int>(pymod(q - q_prev - ((i % 2) ? 2 : 0), 4));
    out[i * bps] = kDqpskBits[dq][0];
    out[i * bps + 1] = kDqpskBits[dq][1];
    for (size_t k = 2; k < bps; ++k) out[i * bps + k] = book.bits[m][k - 2];
    qsym[i] = q;
    q_prev = q;
  }
  std::vector<Alt> alts;
  alts.reserve(nsym);
  for (size_t i = 0; i < nsym; ++i) {
    const auto& sc = second[i];
    const int qp = i ? qsym[i - 1] : q0;
    uint8_t nw[8];
    const int dq = static_cast<int>(pymod(sc.q2 - qp - ((i % 2) ? 2 : 0), 4));
    nw[0] = kDqpskBits[dq][0];
    nw[1] = kDqpskBits[dq][1];
    for (size_t k = 2; k < bps; ++k) nw[k] = book.bits[sc.m2][k - 2];
    std::vector<int64_t> fl;
    for (size_t j = 0; j < bps; ++j)
      if (nw[j] != out[i * bps + j]) fl.push_back(static_cast<int64_t>(i * bps + j));
    if (i + 1 < nsym) {
      const int nd = static_cast<int>(pymod(qsym[i + 1] - sc.q2 - (((i + 1) % 2) ? 2 : 0), 4));
      for (int j = 0; j < 2; ++j)
        if (kDqpskBits[nd][j] != out[(i + 1) * bps + j]) fl.push_back(static_cast<int64_t>((i + 1) * bps + j));
    }
    alts.push_back({sc.margin, fl});
  }
  return Got{out, alts};
}

// ------------------------------------------------------------------------------- interference
// python _excise: STFT cells above the PPDU's own (preamble) spectrum scaled down to it
cfvec excise(const cfvec& y, double fs, std::pair<int64_t, int64_t> ref, size_t nper = 128) {
  const size_t nov = nper * 3 / 4;
  auto st = sig::stft(y, fs, nper, static_cast<long>(nov));
  const Mat<cf>& X = st.Z;
  Mat<float> P(X.rows, X.cols);
  for (size_t i = 0; i < X.data.size(); ++i) {
    const float a = cabs(X.data[i]);
    P.data[i] = a * a;
  }
  const int64_t hop = static_cast<int64_t>(nper - nov);
  const int64_t k0 = ref.first / hop + 1, k1 = std::max(ref.first / hop + 5, ref.second / hop);
  std::vector<float> med(P.rows);
  for (size_t r = 0; r < P.rows; ++r) {
    const int64_t a = std::min<int64_t>(k0, P.cols), b = std::min<int64_t>(std::max(k1, a), P.cols);
    med[r] = np::median(std::vector<float>(P.row(r) + a, P.row(r) + b));
  }
  const std::vector<float> S = ndi::uniform_filter(med, 5);
  const Mat<float> U = ndi::uniform_filter(P, 3, 5);
  Mat<cf> XW(X.rows, X.cols);
  for (size_t r = 0; r < X.rows; ++r)
    for (size_t c = 0; c < X.cols; ++c) {
      const float q = S[r] / (U(r, c) + 1e-30f);
      const float w = std::isnan(q) ? q : std::min(1.0f, q);  // np.minimum propagates NaN
      const cf x = X(r, c);
      XW(r, c) = cf(x.real() * w, x.imag() * w);
    }
  auto yc = sig::istft(XW, fs, nper, static_cast<long>(nov)).x;
  cfvec out = y;
  const size_t n = std::min(y.size(), yc.size());
  std::copy(yc.begin(), yc.begin() + static_cast<ptrdiff_t>(n), out.begin());
  return out;
}

std::optional<Json> retry_excised(const Chan& ch, const Json& r) {
  const double fs = ch.fs;
  const double start = r["start_s"].get<double>(), end = r["end_s"].get<double>();
  const int64_t a = std::max<int64_t>(0, static_cast<int64_t>((start - 60e-6) * fs));
  const int64_t b = std::min<int64_t>(ch.y.size(), static_cast<int64_t>((end + 60e-6) * fs));
  const double pre_us = r["preamble"] == "long" ? 140e-6 : 70e-6;
  const double af = static_cast<double>(a) / fs;
  const std::pair<int64_t, int64_t> ref = {static_cast<int64_t>((start - af) * fs),
                                           static_cast<int64_t>((start - af + pre_us) * fs)};
  const cfvec y(ch.y.begin() + a, ch.y.begin() + std::max(a, b));
  const Chan seg(excise(y, fs, ref), fs, ch.fc);
  for (Json& g : decode_region(seg, 0, b - a, r["channel"].get<int>(), ch.fc)) {
    if (wifi::detail::truthy(wifi::detail::jget(g, "fcs_ok")) && std::fabs(g["start_s"].get<double>() + af - start) < 5e-6) {
      g["start_s"] = g["start_s"].get<double>() + af;
      g["end_s"] = g["end_s"].get<double>() + af;
      g["psdu_demod"] = g["psdu_demod"].get<std::string>() + ", after STFT interference excision";
      return g;
    }
  }
  return std::nullopt;
}

// ------------------------------------------------------------------------------- scan
struct Hit {
  int64_t n;
  std::string pre;
  double amp;
  bool operator<(const Hit& o) const {
    if (n != o.n) return n < o.n;
    if (pre != o.pre) return pre < o.pre;
    return amp < o.amp;
  }
};

// python _sfd_hits: SFD-first acquisition on a fixed 1 us grid at every sample phase
std::vector<Hit> sfd_hits(const Chan& ch, int64_t n0, int64_t n1) {
  const int64_t L = py_round_int(ch.T);
  const int64_t nsym = static_cast<int64_t>(static_cast<double>(n1 - n0) / ch.T) + 8;
  std::vector<Hit> hits;
  const int64_t limit = static_cast<int64_t>(ch.c.size()) - 1;
  for (int64_t ph = 0; ph < L; ++ph) {
    std::vector<int64_t> idx;
    for (int64_t k = 0; k < nsym; ++k) {
      const int64_t v = static_cast<int64_t>(static_cast<double>(n0 + ph) + static_cast<double>(k) * ch.T);
      if (v < limit) idx.push_back(v);
    }
    if (idx.size() < 60) break;
    cdvec cs(idx.size());
    for (size_t q = 0; q < idx.size(); ++q) cs[q] = ch.c[idx[q]];
    for (const auto& [e, pre] : find_sfds(descramble(dbpsk_bits(cs)))) {
      const int64_t lo = std::max<int64_t>(0, e - 48);
      dvec v;
      for (int64_t q = lo; q <= e; ++q) v.push_back(ch.ac[idx[q]]);
      hits.push_back({idx[e], pre, np::median(v)});
    }
  }
  std::sort(hits.begin(), hits.end());
  std::vector<Hit> out;
  for (const auto& h : hits) {  // main lobe and sidelobes of one SFD: keep the peak
    if (!out.empty() && h.n - out.back().n <= 30 && h.pre == out.back().pre) {
      if (h.amp > out.back().amp) out.back() = h;
    } else {
      out.push_back(h);
    }
  }
  return out;
}

}  // namespace

std::vector<Json> decode_region(const Chan& ch, int64_t n0, int64_t n1, int chan, double fc_ch,
                                const std::function<bool(double, double)>& skip) {
  std::vector<Json> out;
  for (const auto& h : sfd_hits(ch, n0, n1)) {
    const int64_t npre = h.pre == "long" ? 144 : 72;
    const double t_start = (static_cast<double>(h.n) - (npre - 1) * ch.T) / ch.fs;
    if (h.amp <= 0 || (skip && skip(t_start, h.amp))) continue;
    double tau = std::max(4.0, static_cast<double>(h.n) - (npre - 1) * ch.T);
    const int64_t i = static_cast<int64_t>(tau);
    const int64_t j = i >= 3 ? i - 3 + static_cast<int64_t>(np::argmax(ch.ac.data() + i - 3, 7)) : i;
    tau = static_cast<double>(j);
    const int64_t e_guess = py_round_int((static_cast<double>(h.n) - tau) / ch.T);
    const auto tr = ch.track(tau, static_cast<int64_t>((static_cast<double>(std::max(n1, h.n)) - tau) / ch.T) + 80,
                             std::nullopt, h.amp);
    if (static_cast<int64_t>(tr.s.size()) < e_guess + 30) continue;
    int64_t sf = -1;
    for (const auto& [e, p] : find_sfds(descramble(dbpsk_bits(tr.s))))
      if (p == h.pre && std::llabs(e - e_guess) <= 3) {
        sf = e;
        break;
      }
    if (sf < 0) continue;
    auto r = Ppdu(ch, tr.s, tr.taus, tr.T, sf, h.pre, h.amp).run();
    if (!r) continue;
    (*r)["channel"] = chan;
    (*r)["ch_freq_mhz"] = py_round(fc_ch / 1e6, 3);
    (*r)["barker_amplitude"] = py_round(h.amp, 1);
    const double hdr_us = h.pre == "long" ? 192 : 96;
    const bool crc = (*r)["plcp_crc_ok"].get<bool>();
    (*r)["end_s"] = (*r)["start_s"].get<double>() + hdr_us * 1e-6 +
                    (crc ? (*r)["length_us"].get<double>() * 1e-6 : 0.0);
    out.push_back(std::move(*r));
  }
  return out;
}

ScanResult scan(const Channelizer& chz, double cap_fc, double cap_fs) {
  const double lo = cap_fc - 0.47 * cap_fs, hi = cap_fc + 0.47 * cap_fs;
  struct Cand {
    double rho;
    int chan;
    double f, t0, t1;
  };
  std::vector<std::pair<int, double>> chs;
  for (const auto& [chan, f] : bands::wifi_24())
    if (!(f - 9e6 < lo || f + 9e6 > hi)) chs.push_back({chan, f});
  auto per = pmap(chs, [&](const std::pair<int, double>& cf_) {
    auto sl = chz.slice(cf_.second, 11e6, kFs);
    std::vector<Cand> out;
    for (const auto& c : candidates(sl.y))
      out.push_back({c.rho, cf_.first, cf_.second, c.a / sl.fs, c.b / sl.fs});
    return out;
  });
  std::vector<Cand> cands;
  for (auto& v : per) cands.insert(cands.end(), v.begin(), v.end());
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return -a.rho < -b.rho; });

  // every channel with a candidate, built up front in parallel (python builds them lazily: same data)
  std::vector<std::pair<int, double>> need;
  for (const auto& c : cands)
    if (std::none_of(need.begin(), need.end(), [&](const auto& p) { return p.first == c.chan; }))
      need.push_back({c.chan, c.f});
  auto built = pmap(need, [&](const std::pair<int, double>& cf_) {
    auto sl = chz.slice(cf_.second, 11e6, kFsDec);
    return std::make_shared<const Chan>(std::move(sl.y), sl.fs, cf_.second);
  });
  std::map<int, std::shared_ptr<const Chan>> chans;
  for (size_t i = 0; i < need.size(); ++i) chans[need[i].first] = built[i];

  std::vector<std::optional<Json>> pk;
  auto dup = [&](int chan, double t0) {
    std::vector<size_t> out;
    for (size_t i = 0; i < pk.size(); ++i)
      if (pk[i] && std::abs((*pk[i])["channel"].get<int>() - chan) <= 3 &&
          std::fabs((*pk[i])["start_s"].get<double>() - t0) < 4e-6)
        out.push_back(i);
    return out;
  };
  auto rank = [](const Json& r) {
    return std::make_tuple(wifi::detail::truthy(wifi::detail::jget(r, "fcs_ok")), r["plcp_crc_ok"].get<bool>(),
                           r["barker_amplitude"].get<double>());
  };
  for (const auto& c : cands) {
    bool inside = false;
    for (const auto& r : pk) {
      if (r && (*r)["plcp_crc_ok"].get<bool>() && (*r)["channel"].get<int>() == c.chan &&
          wifi::detail::truthy(wifi::detail::jget(*r, "fcs_ok")) && c.t0 >= (*r)["start_s"].get<double>() + 200e-6 &&
          c.t1 <= (*r)["end_s"].get<double>() + 10e-6) {
        inside = true;
        break;
      }
    }
    if (inside) continue;
    const Chan& ch = *chans.at(c.chan);
    auto skip = [&](double ts, double amp) {
      for (size_t i : dup(c.chan, ts))
        if ((*pk[i])["barker_amplitude"].get<double>() * 1.2 >= amp) return true;
      return false;
    };
    for (Json& r : decode_region(ch, static_cast<int64_t>(c.t0 * ch.fs), static_cast<int64_t>(c.t1 * ch.fs), c.chan,
                                 c.f, skip)) {
      r["coherence"] = py_round(c.rho, 3);
      const auto same = dup(c.chan, r["start_s"].get<double>());
      if (!same.empty()) {
        auto best = rank(*pk[same[0]]);
        for (size_t i : same) best = std::max(best, rank(*pk[i]));
        if (best >= rank(r)) continue;
      }
      for (size_t i : same) pk[i].reset();
      pk.push_back(std::move(r));
    }
  }
  std::vector<Json> recs;
  for (auto& r : pk)
    if (r) recs.push_back(std::move(*r));
  for (auto& r : recs) {  // last resort for header-ok / FCS-bad PPDUs
    const Json& f = wifi::detail::jget(r, "fcs_ok");
    if (r["plcp_crc_ok"].get<bool>() && f.is_boolean() && !f.get<bool>() && !wifi::detail::truthy(wifi::detail::jget(r, "truncated")) &&
        r["snr_db"].get<double>() > -6) {
      auto g = retry_excised(*chans.at(r["channel"].get<int>()), r);
      if (g) {
        (*g)["coherence"] = r["coherence"];
        r = std::move(*g);
      }
    }
  }
  std::stable_sort(recs.begin(), recs.end(), [](const Json& a, const Json& b) {
    return a["start_s"].get<double>() < b["start_s"].get<double>();
  });
  ScanResult res;
  for (auto& r : recs) {
    const Json& ds = wifi::detail::jget(wifi::detail::jget(r, "mac"), "ds_channel");
    r["seen_on_channels"] = Json::array({r["channel"]});
    std::optional<double> dsf;
    if (wifi::detail::truthy(ds))
      for (const auto& [chn, f] : bands::wifi_24())
        if (ds.is_number_integer() && chn == ds.get<int>()) dsf = f;
    if (dsf) {
      r["channel_source"] = "DS parameter set (beacon)";
      r["decoded_on_channel"] = r["channel"];
      r["channel"] = ds;
      r["ch_freq_mhz"] = *dsf / 1e6;
    } else {
      r["channel_source"] = "strongest Barker coherence among channels within +-15 MHz";
    }
    res.ppdus.push_back(std::move(r));
  }
  // Barker-coherent regions that no decoded PPDU explains
  struct Span {
    int c;
    double s0, e;
  };
  std::vector<Span> spans;
  for (const auto& r : res.ppdus)
    if (r["plcp_crc_ok"].get<bool>())
      spans.push_back({r.contains("decoded_on_channel") ? r["decoded_on_channel"].get<int>() : r["channel"].get<int>(),
                       r["start_s"].get<double>(), r["end_s"].get<double>()});
  auto by_t = cands;
  std::stable_sort(by_t.begin(), by_t.end(), [](const Cand& a, const Cand& b) { return a.t0 < b.t0; });
  for (const auto& c : by_t) {
    bool expl = c.rho < kBarkerReportRho;
    for (const auto& sp : spans)
      expl = expl || (std::abs(sp.c - c.chan) <= 3 && c.t0 < sp.e + 20e-6 && c.t1 > sp.s0 - 20e-6);
    if (expl) continue;
    Json* hit = nullptr;
    for (auto& r : res.unexplained)
      if (std::abs(r["channel"].get<int>() - c.chan) <= 3 && c.t0 < r["t1"].get<double>() &&
          c.t1 > r["t0"].get<double>()) {
        hit = &r;
        break;
      }
    if (hit) {
      Json& r = *hit;
      if (c.rho > r["coherence"].get<double>()) {
        r["channel"] = c.chan;
        r["coherence"] = py_round(c.rho, 3);
        r["ch_freq_mhz"] = py_round(c.f / 1e6, 3);
      }
      r["t0"] = std::min(r["t0"].get<double>(), c.t0);
      r["t1"] = std::max(r["t1"].get<double>(), c.t1);
      continue;
    }
    Json r = Json::object();
    r["t0"] = c.t0;
    r["t1"] = c.t1;
    r["channel"] = c.chan;
    r["ch_freq_mhz"] = py_round(c.f / 1e6, 3);
    r["coherence"] = py_round(c.rho, 3);
    res.unexplained.push_back(std::move(r));
  }
  res.n_candidates = cands.size();
  return res;
}

}  // namespace si::rx::wifi::dsss
