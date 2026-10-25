#include "receivers/gsm/sch.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"

namespace si::rx::gsm {

namespace {

constexpr uint8_t kTsSch[64] = {1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
                                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0,
                                0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1};

// convolutional encoder tables (python _conv_encode_tables)
struct Tables {
  int nxt[16][2];
  int out[16][2][2];
};
const Tables& tables() {
  static const Tables t = [] {
    Tables x{};
    for (int s = 0; s < 16; ++s)
      for (int b = 0; b < 2; ++b) {
        const int reg = (b << 4) | s;
        int d[5];
        for (int i = 0; i < 5; ++i) d[i] = (reg >> (4 - i)) & 1;
        x.out[s][b][0] = d[0] ^ d[3] ^ d[4];
        x.out[s][b][1] = d[0] ^ d[1] ^ d[3] ^ d[4];
        x.nxt[s][b] = reg >> 1;
      }
    return x;
  }();
  return t;
}

}  // namespace

const std::vector<std::string>& gsm_bands() {
  static const std::vector<std::string> b = {"b8_dl", "b3_dl"};
  return b;
}

std::optional<std::pair<int, double>> arfcn_freq(double f) {
  const auto a = bands::gsm_arfcn(f);
  if (!a) return std::nullopt;
  const double fr = *a < 512 ? 935e6 + 0.2e6 * *a : 1805.2e6 + 0.2e6 * (*a - 512);
  return std::make_pair(*a, fr);
}

namespace {

// python find_fcch on complex64 (slice output) or complex128 (CFO-corrected) input: |y|^2 in the input precision
template <class C>
std::optional<Fcch> find_fcch_impl(const std::vector<C>& y, int64_t win) {
  using R = typename C::value_type;
  const int64_t n = static_cast<int64_t>(y.size());
  if (n < win) return std::nullopt;
  // z = y * exp(-2j pi F n / OFS): phase ((-2 pi F) * n) * (1/OFS) (numpy complex / real)
  const double a = -(2.0 * kPi) * kFFcch, inv = 1.0 / kOfs;
  cdvec z(n);
  std::vector<R> p(n);
  for (int64_t i = 0; i < n; ++i) {
    const double ph = (a * static_cast<double>(i)) * inv;
    z[i] = np::cmul(cd(y[i].real(), y[i].imag()), cd(std::cos(ph), std::sin(ph)));
    const R m = np::cabs(y[i]);
    p[i] = m * m;
  }
  // valid convolutions with ones(win): running sums re-anchored every 64 outputs
  const int64_t nv = n - win + 1;
  dvec r(nv);
  cd cs(0, 0);
  double es = 0;
  for (int64_t i = 0; i < nv; ++i) {
    if (i % 64 == 0) {
      cs = cd(0, 0);
      es = 0;
      for (int64_t j = 0; j < win; ++j) {
        cs += z[i + j];
        es += static_cast<double>(p[i + j]);
      }
    } else {
      cs = cs + z[i + win - 1] - z[i - 1];
      es = es + static_cast<double>(p[i + win - 1]) - static_cast<double>(p[i - 1]);
    }
    const double e = es * static_cast<double>(win);
    const double c = np::cabs(cs);
    r[i] = c * c / (e + 1e-30);
  }
  const int64_t k = static_cast<int64_t>(np::argmax(r.data(), r.size()));
  if (r[k] < 0.5) return std::nullopt;
  // S = start of the >= 95% purity plateau (exact burst start)
  int64_t s0 = -1;
  for (int64_t i = std::max<int64_t>(0, k - 148 * kSps); i <= std::min(nv - 1, k + 148 * kSps); ++i)
    if (r[i] > 0.95 * r[k]) {
      s0 = i;
      break;
    }
  cdvec prod(win - 1);
  for (int64_t i = 0; i + 1 < win; ++i) prod[i] = np::cmul(z[k + i + 1], std::conj(z[k + i]));
  const cd sm = np::pairwise_sum(prod.data(), prod.size());
  const double cfo = std::atan2(sm.imag(), sm.real()) * kOfs / (2 * kPi);
  return Fcch{s0, cfo, r[k]};
}

}  // namespace

std::optional<Fcch> find_fcch(const cfvec& y, int64_t win) { return find_fcch_impl(y, win); }
std::optional<Fcch> find_fcch(const cdvec& y, int64_t win) { return find_fcch_impl(y, win); }

namespace {

// python _peak_hz: parabolic interpolation of the power-spectrum peak at bin j on log power (input precision)
template <class R>
double peak_hz(const std::vector<R>& p, int64_t j, int64_t nfft) {
  double l[3];
  for (int d = -1; d <= 1; ++d) {
    const R v = p[static_cast<size_t>((j + d + nfft) % nfft)] + static_cast<R>(1e-30);
    l[d + 1] = static_cast<double>(std::log(v));
  }
  const double den = l[0] - 2 * l[1] + l[2];
  const double d = den < 0 ? 0.5 * (l[0] - l[2]) / den : 0.0;
  const double f = j < nfft / 2 ? static_cast<double>(j) : static_cast<double>(j - nfft);
  return (f + d) * kOfs / static_cast<double>(nfft);
}

// y * exp(-2j pi cfo n / OFS) in complex128 (numpy complex / real: phase ((-2 pi cfo) * n) * (1/OFS))
cdvec rotate(const cfvec& y, double cfo) {
  const double a = -(2.0 * kPi) * cfo, inv = 1.0 / kOfs;
  cdvec out(y.size());
  for (size_t i = 0; i < y.size(); ++i) {
    const double ph = (a * static_cast<double>(i)) * inv;
    out[i] = np::cmul(cd(y[i].real(), y[i].imag()), cd(std::cos(ph), std::sin(ph)));
  }
  return out;
}

}  // namespace

std::optional<Fcch> find_fcch_cfo(const cfvec& y, double max_cfo) {
  const int64_t win = kFcchWin, hop = kFcchWin / 4, nfft = static_cast<int64_t>(kFcchNfft);
  const int64_t len = static_cast<int64_t>(y.size());
  const int64_t m = (len - win) / hop + 1;
  if (len < win || m < 1) return std::nullopt;
  const dvec f = fft::fftfreq(kFcchNfft, 1.0 / kOfs);
  std::vector<int64_t> sel;
  for (int64_t k = 0; k < nfft; ++k)
    if (std::fabs(f[k] - kFFcch) <= max_cfo) sel.push_back(k);
  if (sel.empty()) return std::nullopt;
  // per frame: |np.fft.fft(frame, 2048)|^2 (complex64 -> float32), frame energy sum(|frame|^2) * win + 1e-30 (float32)
  struct Best {
    float pur = -1;
    int64_t i = 0, j = 0;
  };
  std::vector<Best> best_of(static_cast<size_t>(m));
  parallel_for(static_cast<size_t>(m), [&](size_t i) {
    const cf* fr = y.data() + static_cast<int64_t>(i) * hop;
    fvec a2(static_cast<size_t>(win));
    for (int64_t t = 0; t < win; ++t) {
      const float v = np::cabs(fr[t]);
      a2[t] = v * v;
    }
    float e = np::pairwise_sum(a2.data(), a2.size()) * static_cast<float>(win);
    e = e + 1e-30f;
    const cfvec X = fft::np_fft(fr, static_cast<size_t>(win), kFcchNfft);
    Best b;
    for (size_t jj = 0; jj < sel.size(); ++jj) {
      const float v = np::cabs(X[sel[jj]]);
      const float pur = (v * v) / e;
      if (pur > b.pur) b = {pur, static_cast<int64_t>(i), static_cast<int64_t>(jj)};
    }
    best_of[i] = b;
  });
  Best b = best_of[0];
  for (const Best& x : best_of)
    if (x.pur > b.pur) b = x;  // first maximum in row-major order (np.argmax)
  if (b.pur < 0.3f) return std::nullopt;
  const cf* fr = y.data() + b.i * hop;
  const cfvec X = fft::np_fft(fr, static_cast<size_t>(win), kFcchNfft);
  fvec Xp(kFcchNfft);
  for (size_t k = 0; k < kFcchNfft; ++k) {
    const float v = np::cabs(X[k]);
    Xp[k] = v * v;
  }
  const double cfo0 = peak_hz(Xp, sel[b.j], nfft) - kFFcch;
  const cdvec y2 = rotate(y, cfo0);
  const int64_t a0 = std::max<int64_t>(0, b.i * hop - win);
  const int64_t a1 = std::min<int64_t>(len, b.i * hop + 2 * win + 148 * kSps);
  const auto fc = find_fcch(cdvec(y2.begin() + a0, y2.begin() + std::max(a0, a1)), win);
  if (!fc) return std::nullopt;
  const int64_t s0 = fc->start + a0;
  const int64_t e1 = std::min<int64_t>(len, s0 + win);
  const cdvec P = fft::fft(y2.data() + s0, static_cast<size_t>(std::max<int64_t>(e1 - s0, 0)), kFcchNfft);
  dvec Pp(kFcchNfft);
  for (size_t k = 0; k < kFcchNfft; ++k) {
    const double v = np::cabs(P[k]);
    Pp[k] = v * v;
  }
  int64_t kbest = -1;
  for (int64_t k = 0; k < nfft; ++k)
    if (std::fabs(f[k] - kFFcch) <= 2 * kOfs / static_cast<double>(nfft) && (kbest < 0 || Pp[k] > Pp[kbest])) kbest = k;
  return Fcch{s0, cfo0 + peak_hz(Pp, kbest, nfft) - kFFcch, fc->purity};
}

std::optional<Sch> demod_sch(const cdvec& y, int64_t start) {
  const int64_t n = static_cast<int64_t>(y.size());
  struct Best {
    int m;
    int64_t s;
    dvec dph;
    std::array<uint8_t, 149> b;
  };
  std::optional<Best> best;
  dvec dph(148);
  std::array<uint8_t, 149> b{};
  for (int64_t off = -20 * kSps; off <= 20 * kSps; ++off) {
    const int64_t s = start + off;
    if (s < 0 || s + 150 * kSps >= n) continue;
    for (int i = 0; i < 148; ++i) {
      const cd v = np::cmul(y[s + (i + 1) * kSps], std::conj(y[s + i * kSps]));
      dph[i] = std::atan2(v.imag(), v.real());
    }
    for (int pol : {1, -1}) {
      b[0] = 0;
      for (int i = 0; i < 148; ++i) {
        const int sg = dph[i] > 0 ? 1 : (dph[i] < 0 ? -1 : 0);
        const uint8_t d = pol * sg < 0;
        b[i + 1] = d ^ b[i];
      }
      int m = 0;
      for (int i = 0; i < 64; ++i) m += b[42 + i] == kTsSch[i];
      m = std::max(m, 64 - m);
      if (!best || m > best->m) best = Best{m, s, dph, b};
    }
  }
  if (!best || best->m < 56) return std::nullopt;
  auto bb = best->b;
  int same = 0;
  for (int i = 0; i < 64; ++i) same += bb[42 + i] == kTsSch[i];
  if (same < 32)
    for (auto& v : bb) v ^= 1;
  dvec rel(149);
  rel[0] = 1.0;
  for (int i = 0; i < 148; ++i) rel[i + 1] = std::clamp(std::fabs(best->dph[i]) / (kPi / 2), 0.05, 1.5);
  Sch out;
  for (int i = 3; i < 42; ++i) out.soft.push_back((1.0 - 2.0 * bb[i]) * rel[i]);
  for (int i = 106; i < 145; ++i) out.soft.push_back((1.0 - 2.0 * bb[i]) * rel[i]);
  out.ts_match = best->m;
  out.sample = best->s;
  return out;
}

std::vector<int> viterbi(const dvec& soft) {
  const Tables& T = tables();
  const size_t n = soft.size() / 2;
  constexpr double INF = 1e18;
  std::array<double, 16> pm;
  pm.fill(INF);
  pm[0] = 0;
  std::vector<std::array<int, 16>> hist(n), hbit(n);
  for (size_t i = 0; i < n; ++i) {
    hist[i].fill(0);
    hbit[i].fill(0);
    const double s0 = soft[2 * i], s1 = soft[2 * i + 1];
    std::array<double, 16> npm;
    npm.fill(INF);
    for (int s = 0; s < 16; ++s) {
      if (pm[s] >= INF) continue;
      for (int b = 0; b < 2; ++b) {
        const double cost = (T.out[s][b][0] ? s0 : -s0) + (T.out[s][b][1] ? s1 : -s1);
        const int t = T.nxt[s][b];
        if (pm[s] + cost < npm[t]) {
          npm[t] = pm[s] + cost;
          hist[i][t] = s;
          hbit[i][t] = b;
        }
      }
    }
    pm = npm;
  }
  std::vector<int> bits(n);
  int s = 0;  // 4 zero tail bits terminate in state 0
  for (size_t i = n; i-- > 0;) {
    bits[i] = hbit[i][s];
    s = hist[i][s];
  }
  return bits;
}

bool crc10_ok(const std::vector<int>& bits25, const std::vector<int>& par10) {
  const int poly = 0b10101110101;  // D10+D8+D6+D5+D4+D2+1
  std::vector<int> reg(bits25.begin(), bits25.end());
  reg.resize(35, 0);
  for (int i = 0; i < 25; ++i)
    if (reg[i])
      for (int j = 0; j < 11; ++j) reg[i + j] ^= (poly >> (10 - j)) & 1;
  for (int i = 0; i < 10; ++i)
    if ((reg[25 + i] ^ 1) != (par10[i] & 1)) return false;  // remainder is inverted
  return true;
}

Json sch_fields(const std::vector<int>& b) {
  int by[4];
  for (int k = 0; k < 4; ++k) {
    by[k] = 0;
    for (int j = 0; j < 8; ++j)
      if (8 * k + j < 25) by[k] += b[8 * k + j] << j;
  }
  const int bsic = (by[0] >> 2) & 0x3F;
  const int t1 = ((by[0] & 3) << 9) | (by[1] << 1) | (by[2] >> 7);
  const int t2 = (by[2] >> 2) & 0x1F;
  const int t3p = ((by[2] & 3) << 1) | (by[3] & 1);
  const int t3 = 10 * t3p + 1;
  const int d = ((t3 - t2) % 26 + 26) % 26;
  const int fn = 51 * d + t3 + 51 * 26 * t1;
  Json f = Json::object();
  f["bsic"] = bsic;
  f["ncc"] = bsic >> 3;
  f["bcc"] = bsic & 7;
  f["t1"] = t1;
  f["t2"] = t2;
  f["t3_prime"] = t3p;
  f["frame_number"] = fn;
  f["fields_plausible"] = t2 <= 25 && t3p <= 4;
  return f;
}

}  // namespace si::rx::gsm
