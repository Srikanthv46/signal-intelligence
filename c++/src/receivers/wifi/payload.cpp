// Port of python receivers/wifi/payload.py (see payload.hpp). numpy arithmetic as in
// ofdm_rx.cpp: np::cmul / np::cdiv / np::cabs, pairwise sums, convolutions summed in order.
#include "receivers/wifi/payload.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "receivers/wifi/ldpc.hpp"
#include "receivers/wifi/mac80211.hpp"

namespace si::rx::wifi::payload {

namespace {

using np::cabs;
using np::cabs2;
using np::cdiv;
using np::cmul;
using ofdm::bin_of;

constexpr double kPi = 3.141592653589793;
constexpr double kLogE2 = 0.693147180559945309417232121458176568;
inline cd expj(double x) { return {std::cos(x), std::sin(x)}; }
inline double angle(cd z) { return std::atan2(z.imag(), z.real()); }
inline cd conjd(cd z) { return {z.real(), -z.imag()}; }
inline double phase_step(double f, double fs) { return (-2.0 * kPi * f) / fs; }

double kmod(int nbpsc) {
  switch (nbpsc) {
    case 1: return 1.0;
    case 2: return 1 / std::sqrt(2.0);
    case 4: return 1 / std::sqrt(10.0);
    default: return 1 / std::sqrt(42.0);
  }
}

Rate nonht_mcs(int rate, int* nbpsc) {
  switch (rate) {
    case 6: *nbpsc = 1; return {1, 2};
    case 9: *nbpsc = 1; return {3, 4};
    case 12: *nbpsc = 2; return {1, 2};
    case 18: *nbpsc = 2; return {3, 4};
    case 24: *nbpsc = 4; return {1, 2};
    case 36: *nbpsc = 4; return {3, 4};
    case 48: *nbpsc = 6; return {2, 3};
    case 54: *nbpsc = 6; return {3, 4};
    default: throw std::out_of_range(std::to_string(rate));  // python KeyError
  }
}

Rate ht_mcs(int mcs, int* nbpsc) {
  static const int nb[8] = {1, 2, 2, 4, 4, 6, 6, 6};
  static const Rate r[8] = {{1, 2}, {1, 2}, {3, 4}, {1, 2}, {3, 4}, {2, 3}, {3, 4}, {5, 6}};
  *nbpsc = nb[mcs];
  return r[mcs];
}

// python pam_llr(x, bits) for one axis value
void pam_llr(double x, int per, double* out) {
  if (per == 2) {
    out[0] = x;
    out[1] = 2 - std::fabs(x);
  } else {
    out[0] = x;
    out[1] = 4 - std::fabs(x);
    out[2] = 2 - std::fabs(4 - std::fabs(x));
  }
}

// python _levels(per): PAM levels and their bit labels
struct Pam {
  std::vector<double> lv;
  std::vector<std::array<int, 3>> lab;
};
const Pam& pam(int per) {
  static const std::array<Pam, 2> P = [] {
    std::array<Pam, 2> a;
    for (int per2 = 2; per2 <= 3; ++per2) {
      Pam& p = a[per2 - 2];
      for (int v = -((1 << per2) - 1); v < (1 << per2); v += 2) {
        double l[3];
        pam_llr(v, per2, l);
        p.lv.push_back(v);
        p.lab.push_back({l[0] > 0, l[1] > 0, per2 == 3 ? l[2] > 0 : 0});
      }
    }
    return a;
  }();
  return P[per - 2];
}

inline double logaddexp(double x, double y) {
  if (x == y) return x + kLogE2;
  const double t = x - y;
  if (t > 0) return x + std::log1p(std::exp(-t));
  if (t <= 0) return y + std::log1p(std::exp(t));
  return t;
}

std::vector<uint8_t> crc8_lsb(const uint8_t* two) {  // python _crc8_delim bits
  std::vector<uint8_t> bits;
  for (int i = 0; i < 2; ++i)
    for (int b = 0; b < 8; ++b) bits.push_back((two[i] >> b) & 1);
  return ofdm::crc8_sig(bits.data(), bits.size());
}

std::vector<uint8_t> pack_lsb(const std::vector<uint8_t>& descr, size_t lo, size_t nbytes) {
  std::vector<uint8_t> out(nbytes, 0);
  for (size_t i = 0; i < nbytes; ++i) {
    int v = 0;
    for (int j = 0; j < 8; ++j) {
      const size_t k = lo + 8 * i + j;
      if (k < descr.size() && descr[k]) v |= 1 << j;
    }
    out[i] = static_cast<uint8_t>(v);
  }
  return out;
}

int64_t seed_int(const std::vector<uint8_t>& s) {
  int64_t v = 0;
  for (size_t i = 0; i < s.size(); ++i) v += static_cast<int64_t>(s[i]) << i;
  return v;
}

// fft of 64 samples of seg starting at off - back (zero padded past the end), / 64 * sqrt(52)
cdvec fft_at(const cdvec& seg, int64_t off, int back) {
  cdvec w(64, cd(0, 0));
  const int64_t a = off - back;
  for (int64_t i = 0; i < 64; ++i)
    if (a + i >= 0 && a + i < static_cast<int64_t>(seg.size())) w[i] = seg[a + i];
  cdvec Y = fft::fft(w);
  const double s52 = std::sqrt(52.0);
  for (auto& z : Y) z = cmul(cdiv(z, cd(64, 0)), cd(s52, 0));
  return Y;
}

// sum_j a[k + j] * w over a box of L, "valid", in order (np.convolve with a constant kernel)
template <class T>
std::vector<T> box_valid(const std::vector<T>& a, size_t L, double w) {
  std::vector<T> out(a.size() - L + 1);
  for (size_t k = 0; k < out.size(); ++k) {
    T s = T(0);
    for (size_t j = 0; j < L; ++j) {
      if constexpr (std::is_same_v<T, cd>) s += cmul(a[k + j], cd(w, 0));
      else s += a[k + j] * w;
    }
    out[k] = s;
  }
  return out;
}

}  // namespace

std::pair<std::vector<int>, std::vector<int>> deinterleave_indices(int ncbps, int nbpsc, int ncol) {
  const int s = std::max(nbpsc / 2, 1);
  std::vector<int> j(ncbps), inv(ncbps);
  for (int k = 0; k < ncbps; ++k) {
    const int i = (ncbps / ncol) * (k % ncol) + k / ncol;
    j[k] = s * (i / s) + (i + ncbps - (ncol * i / ncbps)) % s;
  }
  for (int k = 0; k < ncbps; ++k) inv[j[k]] = k;
  return {j, inv};
}

dvec depuncture(const dvec& soft, Rate rate) {
  std::vector<int> pat;
  if (rate == Rate{3, 4}) pat = {1, 1, 1, 0, 0, 1};
  else if (rate == Rate{2, 3}) pat = {1, 1, 1, 0};
  else if (rate == Rate{5, 6}) pat = {1, 1, 1, 0, 0, 1, 1, 0, 0, 1};
  else return soft;
  int keep_per = 0;
  for (int v : pat) keep_per += v;
  const size_t nin = soft.size();
  const size_t reps = (nin + keep_per - 1) / keep_per;
  dvec out(reps * pat.size(), 0.0);
  size_t m = 0;
  for (size_t i = 0; i < out.size() && m < nin; ++i)
    if (pat[i % pat.size()]) out[i] = soft[m++];
  out.resize(out.size() / 2 * 2);
  return out;
}

std::vector<uint8_t> scrambler_seq(const std::vector<uint8_t>& seed7, size_t n) {
  uint8_t s[7];
  for (int i = 0; i < 7; ++i) s[i] = seed7[i];
  std::vector<uint8_t> out(n);
  for (size_t k = 0; k < n; ++k) {
    const uint8_t fb = s[6] ^ s[3];
    out[k] = fb;
    for (int i = 6; i > 0; --i) s[i] = s[i - 1];
    s[0] = fb;
  }
  return out;
}

std::pair<std::vector<uint8_t>, std::vector<uint8_t>> descramble(const std::vector<uint8_t>& b) {
  std::vector<uint8_t> best;
  if (b.size() >= 7) {
    for (int seed = 1; seed < 128; ++seed) {
      std::vector<uint8_t> s7(7);
      for (int i = 0; i < 7; ++i) s7[i] = (seed >> i) & 1;
      if (scrambler_seq(s7, 7) == std::vector<uint8_t>(b.begin(), b.begin() + 7)) {
        best = s7;
        break;
      }
    }
  }
  if (best.empty()) best = {1, 0, 1, 1, 1, 0, 1};
  const auto full = scrambler_seq(best, b.size());
  std::vector<uint8_t> out(b.size());
  for (size_t i = 0; i < b.size(); ++i) out[i] = b[i] ^ full[i];
  return {out, best};
}

dvec llr_exact(const cdvec& sym, int nbpsc, const dvec& var) {
  const double k = kmod(nbpsc);
  const size_t n = sym.size();
  dvec out(n * nbpsc);
  for (size_t i = 0; i < n; ++i) {
    const double v = std::max(var[i], 1e-9);
    if (nbpsc == 1) {
      out[i] = -4 * sym[i].real() * k / v;
    } else if (nbpsc == 2) {
      out[2 * i] = -4 * sym[i].real() * k / v;
      out[2 * i + 1] = -4 * sym[i].imag() * k / v;
    } else {
      const int per = nbpsc / 2;
      const Pam& P = pam(per);
      const double s2 = v / 2 / (k * k);
      const size_t nl = P.lv.size();
      double d[8];
      for (int ax = 0; ax < 2; ++ax) {
        const double x = ax == 0 ? sym[i].real() : sym[i].imag();
        for (size_t l = 0; l < nl; ++l) {
          const double t = x / k - P.lv[l];
          d[l] = -(t * t) / (2 * s2);
        }
        for (int b = 0; b < per; ++b) {
          double l0 = -std::numeric_limits<double>::infinity(), l1 = l0;
          for (size_t l = 0; l < nl; ++l) {
            if (P.lab[l][b] == 0) l0 = logaddexp(l0, d[l]);
            else l1 = logaddexp(l1, d[l]);
          }
          out[i * nbpsc + ax * per + b] = l0 - l1;
        }
      }
    }
  }
  return out;
}

cdvec smooth_channel(const cdvec& h, const std::vector<int>& sc_in, int w) {
  if (w <= 0) return h;
  std::vector<int> sc = sc_in;
  std::sort(sc.begin(), sc.end());
  const size_t n = sc.size();
  cdvec v(n);
  for (size_t i = 0; i < n; ++i) v[i] = h[bin_of(sc[i])];
  cdvec t(n > 0 ? n - 1 : 0);
  for (size_t i = 0; i + 1 < n; ++i)
    t[i] = cmul(cmul(v[i + 1], conjd(v[i])), cd(sc[i + 1] - sc[i] == 1 ? 1.0 : 0.0, 0));
  const double tau = angle(np::sum(t));
  cdvec g(n);
  for (size_t i = 0; i < n; ++i) g[i] = cmul(v[i], expj(-tau * sc[i]));
  cdvec ext;
  for (int i = 0; i < w; ++i) ext.push_back(g[0]);
  ext.insert(ext.end(), g.begin(), g.end());
  for (int i = 0; i < w; ++i) ext.push_back(g[n - 1]);
  const cdvec c = box_valid(ext, 2 * w + 1, 1.0 / (2 * w + 1));
  cdvec out = h;
  for (size_t i = 0; i < n; ++i) out[bin_of(sc[i])] = cmul(c[i], expj(tau * sc[i]));
  return out;
}

std::vector<std::pair<std::vector<uint8_t>, bool>> deaggregate(const std::vector<uint8_t>& data) {
  std::vector<std::pair<std::vector<uint8_t>, bool>> out;
  size_t i = 0;
  const size_t n = data.size();
  while (i + 4 <= n) {
    const uint8_t* d = &data[i];
    if (d[3] != 0x4E) {
      size_t nxt = i + 4;
      while (nxt < n && data[nxt] != 0x4E) ++nxt;
      if (nxt >= n) break;
      i = nxt - 3;
      continue;
    }
    const size_t ln = ((((d[1] << 8) | d[0]) >> 4) & 0xFFF);
    const auto c = crc8_lsb(d);
    int cv = 0;
    for (int b = 0; b < 8; ++b) cv += c[b] << b;
    const bool ok = cv == d[2];
    i += 4;
    if (ln) {
      const size_t e = std::min(n, i + ln);
      out.push_back({std::vector<uint8_t>(data.begin() + std::min(i, n), data.begin() + e), ok});
      i += ln;
      i += (4 - i % 4) % 4;
    }
  }
  return out;
}

bool ht_mcs_known(std::optional<int> mcs) { return mcs && *mcs >= 0 && *mcs <= 7; }

double ht_rate_mbps(int mcs, bool short_gi) {
  int nb = 0;
  const Rate r = ht_mcs(mcs, &nb);
  return py_round(static_cast<double>(52 * nb * r.first) / r.second / (short_gi ? 3.6 : 4.0), 1);
}

std::pair<bool, std::string> ht_decodable(std::optional<int> mcs, int ldpc, int stbc, int n_ess, double bw_mhz) {
  (void)ldpc;
  if (!ht_mcs_known(mcs))
    return {false, "MCS " + (mcs ? std::to_string(*mcs) : std::string("None")) + ": 2+ spatial streams (needs as many antennas)"};
  if (bw_mhz != 20) return {false, py_str(Json(static_cast<int64_t>(bw_mhz))) + " MHz HT payload (only 20 MHz implemented)"};
  if (stbc > 1) return {false, "space-time block coding with STBC=" + std::to_string(stbc) + " (N_STS > 2 for one stream)"};
  if (n_ess) return {false, "extension spatial streams"};
  return {true, ""};
}

int vht_uplink_paid(const std::string& bssid) {
  int o[6];
  for (int i = 0; i < 6; ++i) o[i] = std::stoi(bssid.substr(3 * i, 2), nullptr, 16);
  return ((o[4] >> 7) | (o[5] << 1)) & 0x1FF;
}

// ------------------------------------------------------------------ non-HT
NonHt decode_nonht_payload(const ofdm::RxState& rx, int rate_mbps, int length) {
  const cdvec& r = *rx.r;
  const int64_t start = rx.start;
  const int back = rx.back;
  int nbpsc = 0;
  const Rate cr = nonht_mcs(rate_mbps, &nbpsc);
  const int ncbps = 48 * nbpsc;
  const int ndbps = ncbps * cr.first / cr.second;
  const int n_sym = (16 + 8 * length + 6 + ndbps - 1) / ndbps;
  std::vector<uint8_t> is_used(64, 0);
  for (int b : rx.used) is_used[b] = 1;
  const auto& DS = ofdm::kDataSc;
  const auto inv = deinterleave_indices(ncbps, nbpsc).second;
  const int64_t want = 400 + 80 * static_cast<int64_t>(n_sym) + 80;
  const int64_t nseg = std::max<int64_t>(0, std::min<int64_t>(want, static_cast<int64_t>(r.size()) - start));
  const double bs = phase_step(rx.cfo, ofdm::kFsWifi);
  cdvec seg(nseg);
  for (int64_t i = 0; i < nseg; ++i) seg[i] = cmul(r[start + i], expj(bs * static_cast<double>(start + i)));
  const double km = kmod(nbpsc);
  dvec soft;
  soft.reserve(static_cast<size_t>(n_sym) * ncbps);
  dvec coded(ncbps), llr(ncbps);
  for (int k = 0; k < n_sym; ++k) {
    const cdvec Y = fft_at(seg, 400 + 80 * k + 16, back);
    cdvec Z(64, cd(0, 0));
    for (int b : rx.used) Z[b] = cdiv(Y[b], rx.H[b]);
    cd pil[4];
    double w[4];
    const double pp = ofdm::pilot_pol()[(k + 1) % 127];
    for (int j = 0; j < 4; ++j) {
      const int b = bin_of(ofdm::kPilotSc[j]);
      pil[j] = cdiv(Z[b], cd(ofdm::kPilotBase[j] * pp, 0));
      w[j] = rx.W_rel[b];
    }
    auto wsum = [&]() {
      cd t[4];
      for (int j = 0; j < 4; ++j) t[j] = cmul(cd(w[j], 0), pil[j]);
      return np::pairwise_sum(t, 4);
    };
    int nlive = 0;
    for (double v : w) nlive += v > 0;
    if (nlive >= 3) {
      const double ph0 = angle(wsum());
      double err[4];
      for (int j = 0; j < 4; ++j) err[j] = w[j] > 0 ? cabs(cmul(pil[j], expj(-ph0)) - cd(1, 0)) : -1;
      const int am = static_cast<int>(std::max_element(err, err + 4) - err);
      if (err[am] > 0.8) w[am] = 0;
    }
    const cd rot = expj(-angle(wsum()));
    for (int i = 0; i < 48; ++i) {
      const int b = bin_of(DS[i]);
      const cd d = cmul(Z[b], rot);
      const double conf = rx.W_rel[b];
      const bool live = is_used[b];
      const double I = d.real() / km, Q = d.imag() / km;
      double v[6];
      if (nbpsc == 1) {
        v[0] = I;
      } else if (nbpsc == 2) {
        v[0] = I;
        v[1] = Q;
      } else {
        const int per = nbpsc / 2;
        pam_llr(I, per, v);
        pam_llr(Q, per, v + per);
      }
      for (int t = 0; t < nbpsc; ++t) llr[i * nbpsc + t] = live ? v[t] * conf : 0.0;
    }
    for (int t = 0; t < ncbps; ++t) coded[inv[t]] = llr[t];
    soft.insert(soft.end(), coded.begin(), coded.end());
  }
  const auto databits = ofdm::viterbi(depuncture(soft, cr));
  auto [descr, seed] = descramble(databits);
  NonHt out;
  out.n_sym = n_sym;
  out.seed = seed;
  out.mpdu = pack_lsb(descr, 16, static_cast<size_t>(length));
  if (length >= 4) {
    const size_t nb = out.mpdu.size() - 4;
    const uint32_t calc = mac80211::crc32(out.mpdu.data(), nb);
    const uint32_t fcs = out.mpdu[nb] | (out.mpdu[nb + 1] << 8) | (out.mpdu[nb + 2] << 16) |
                         (static_cast<uint32_t>(out.mpdu[nb + 3]) << 24);
    out.fcs_ok = calc == fcs;
    out.info = mac80211::parse(out.mpdu.data(), nb);
  }
  bool any = false;
  for (size_t i = 7; i < 16 && i < descr.size(); ++i) any |= descr[i] != 0;
  out.service_ok = descr.size() >= 16 && !any;
  return out;
}

// ------------------------------------------------------------------ HT
namespace {

constexpr int kHtLtfP[2][2] = {{1, -1}, {1, 1}};
constexpr int kPsi1[4] = {1, 1, 1, -1};
constexpr int kPsi2[2][4] = {{1, 1, -1, -1}, {1, -1, -1, 1}};

// python ht_pilots(n, nsts)[i][j]
double ht_pilot(int n, int nsts, int i, int j) {
  const int psi = nsts == 1 ? kPsi1[(n + j) % 4] : kPsi2[i][(n + j) % 4];
  return psi * ofdm::pilot_pol()[(n + 3) % 127];
}

struct HtCtx {
  int nsts = 1, n_sym = 0, nbpsc = 1, ncbps = 52;
  std::vector<cdvec> Ys, Yl;
  std::vector<int> used, data_bins, pb;
  std::vector<uint8_t> live_data;  // per data subcarrier
  dvec Lh, nshape;
  double base = 1.0;
};

using Hmat = std::vector<cdvec>;  // nsts x 64

struct Tracked {
  std::vector<cdvec> Zc;
  double lvl = 0;
};

Tracked tracked(const HtCtx& c, const Hmat& Hs) {
  const int n_sym = c.n_sym;
  std::vector<std::array<cd, 4>> pred(n_sym);
  dvec slope(n_sym);
  const double psc[4] = {-21, -7, 7, 21};
  for (int k = 0; k < n_sym; ++k) {
    for (int j = 0; j < 4; ++j) {
      cd s = cmul(cd(ht_pilot(k, c.nsts, 0, j), 0), Hs[0][c.pb[j]]);
      for (int i = 1; i < c.nsts; ++i) s += cmul(cd(ht_pilot(k, c.nsts, i, j), 0), Hs[i][c.pb[j]]);
      pred[k][j] = s;
    }
    // _track
    cd z[4];
    for (int j = 0; j < 4; ++j) z[j] = cmul(c.Ys[k][c.pb[j]], conjd(pred[k][j]));
    const double a0 = angle(np::pairwise_sum(z, 4));
    double num[4], den[4];
    for (int j = 0; j < 4; ++j) {
      const double t = angle(z[j]) - a0;
      const double d = angle(expj(t));
      const double w = cabs2(pred[k][j]);
      num[j] = w * d * psc[j];
      den[j] = w * (psc[j] * psc[j]);
    }
    slope[k] = np::pairwise_sum(num, 4) / (np::pairwise_sum(den, 4) + 1e-12);
  }
  if (n_sym >= 5) {
    const int w = std::min(9, n_sym | 1);
    dvec pad;
    for (int i = 0; i < w / 2; ++i) pad.push_back(slope.front());
    pad.insert(pad.end(), slope.begin(), slope.end());
    for (int i = 0; i < w / 2; ++i) pad.push_back(slope.back());
    slope = box_valid(pad, w, 1.0 / w);
  }
  Tracked out;
  dvec resid;
  resid.reserve(4 * n_sym);
  for (int k = 0; k < n_sym; ++k) {
    cdvec Yk(64);
    for (int b = 0; b < 64; ++b) {
      const int kk = b < 32 ? b : b - 64;
      Yk[b] = cmul(c.Ys[k][b], expj(-slope[k] * kk));
    }
    cd t[4];
    for (int j = 0; j < 4; ++j) t[j] = cmul(Yk[c.pb[j]], conjd(pred[k][j]));
    const cd rot = expj(-angle(np::pairwise_sum(t, 4)));
    for (auto& y : Yk) y = cmul(y, rot);
    for (int j = 0; j < 4; ++j) resid.push_back(cabs2(Yk[c.pb[j]] - pred[k][j]));
    out.Zc.push_back(std::move(Yk));
  }
  out.lvl = resid.empty() ? c.base : np::mean(resid);
  return out;
}

struct Combined {
  std::vector<cdvec> pts;
  std::vector<dvec> var;
};

Combined combine(const HtCtx& c, const Hmat& Hs, const std::vector<cdvec>& Zc, double lvl) {
  const size_t nd = c.data_bins.size();
  const double inf = std::numeric_limits<double>::infinity();
  dvec noise(nd);
  for (size_t i = 0; i < nd; ++i) noise[i] = lvl * c.nshape[c.data_bins[i]];
  Combined out;
  if (c.nsts == 1) {
    cdvec h(nd);
    dvec g(nd), var(nd);
    for (size_t i = 0; i < nd; ++i) {
      h[i] = c.live_data[i] ? Hs[0][c.data_bins[i]] : cd(1, 0);
      g[i] = c.live_data[i] ? cabs2(h[i]) : inf;
      var[i] = noise[i] / g[i];
    }
    for (int k = 0; k < c.n_sym; ++k) {
      cdvec p(nd);
      for (size_t i = 0; i < nd; ++i) p[i] = c.live_data[i] ? cdiv(Zc[k][c.data_bins[i]], h[i]) : cd(0, 0);
      out.pts.push_back(std::move(p));
      out.var.push_back(var);
    }
  } else {
    cdvec h1(nd), h2(nd);
    dvec g(nd), gs(nd), var(nd);
    for (size_t i = 0; i < nd; ++i) {
      h1[i] = Hs[0][c.data_bins[i]];
      h2[i] = Hs[1][c.data_bins[i]];
      g[i] = c.live_data[i] ? cabs2(h1[i]) + cabs2(h2[i]) : inf;
      gs[i] = c.live_data[i] ? g[i] : 1.0;
      var[i] = noise[i] / g[i];
    }
    for (int m = 0; m < c.n_sym / 2; ++m) {
      const cdvec& Za = Zc[2 * m];
      const cdvec& Zb = Zc[2 * m + 1];
      cdvec pa(nd), pbv(nd);
      for (size_t i = 0; i < nd; ++i) {
        const cd ya = Za[c.data_bins[i]], yb = Zb[c.data_bins[i]];
        const cd sa = cdiv(cmul(conjd(h1[i]), ya) + cmul(h2[i], conjd(yb)), cd(gs[i], 0));
        const cd sb = cdiv(cmul(conjd(h1[i]), yb) - cmul(h2[i], conjd(ya)), cd(gs[i], 0));
        pa[i] = c.live_data[i] ? sa : cd(0, 0);
        pbv[i] = c.live_data[i] ? sb : cd(0, 0);
      }
      out.pts.push_back(std::move(pa));
      out.var.push_back(var);
      out.pts.push_back(std::move(pbv));
      out.var.push_back(var);
    }
  }
  return out;
}

// python _qam_points (including its uint8 wrap-around: 2*0-1 == 255 for BPSK/QPSK)
cd qam_point(const uint8_t* bits, int nbpsc) {
  const double k = kmod(nbpsc);
  if (nbpsc == 1) return cd(static_cast<uint8_t>(2 * bits[0] - 1) * 1.0, 0);
  if (nbpsc == 2) {
    const double i = static_cast<uint8_t>(2 * bits[0] - 1), q = static_cast<uint8_t>(2 * bits[1] - 1);
    return cmul(cd(i, q), cd(k, 0));
  }
  const int per = nbpsc / 2;
  const Pam& P = pam(per);
  auto level = [&](const uint8_t* b) {
    for (size_t l = 0; l < P.lv.size(); ++l) {
      bool eq = true;
      for (int t = 0; t < per; ++t) eq &= P.lab[l][t] == b[t];
      if (eq) return P.lv[l];
    }
    return 0.0;
  };
  return cmul(cd(level(bits), level(bits + per)), cd(k, 0));
}

Hmat data_aided_channel(const HtCtx& c, const std::map<int64_t, std::vector<uint8_t>>& known,
                        const ldpc::Params& LP, const std::vector<int64_t>& offs, const std::vector<cdvec>& Zc,
                        const Hmat& Hs) {
  const int n_sym = c.n_sym, nsts = c.nsts, nbpsc = c.nbpsc, ncbps = c.ncbps;
  const size_t nb = static_cast<size_t>(n_sym) * ncbps;
  std::vector<uint8_t> bits(nb, 0), have(nb, 0);
  for (const auto& [i, cw] : known) {
    const auto tx = ldpc::transmitted(cw, LP, static_cast<size_t>(i));
    const int64_t lo = offs[i];
    const int64_t hi = std::min<int64_t>(lo + static_cast<int64_t>(tx.size()), static_cast<int64_t>(nb));
    for (int64_t t = lo; t < hi; ++t) {
      bits[t] = tx[t - lo];
      have[t] = 1;
    }
  }
  const int nsc = ncbps / nbpsc;
  std::vector<cd> d(static_cast<size_t>(n_sym) * nsc, cd(0, 0));
  std::vector<uint8_t> hv(static_cast<size_t>(n_sym) * nsc, 0);
  for (int k = 0; k < n_sym; ++k)
    for (int s = 0; s < nsc; ++s) {
      const size_t o = (static_cast<size_t>(k) * nsc + s) * nbpsc;
      bool all = true;
      for (int t = 0; t < nbpsc; ++t) all &= have[o + t] != 0;
      hv[k * nsc + s] = all;
      if (all) d[k * nsc + s] = qam_point(&bits[o], nbpsc);
    }
  auto Y = [&](int k, int s) { return Zc[k][c.data_bins[s]]; };
  Hmat out = Hs;
  std::vector<Mat<cd>> G(nsc, Mat<cd>(nsts, nsts, cd(0, 0)));
  std::vector<cdvec> cv(nsc, cdvec(nsts, cd(0, 0)));
  for (int n = 0; n < nsts; ++n) {
    double a[2];
    for (int i = 0; i < nsts; ++i) a[i] = nsts > 1 ? kHtLtfP[i][n] : 1;
    for (int s = 0; s < nsc; ++s) {
      for (int i = 0; i < nsts; ++i)
        for (int j = 0; j < nsts; ++j) G[s](i, j) += cd(a[i] * a[j], 0);
      const double lk = c.Lh[c.data_bins[s]];
      const cd y = cdiv(c.Yl[n][c.data_bins[s]], cd(lk != 0 ? lk : 1, 0));
      for (int i = 0; i < nsts; ++i) cv[s][i] += cmul(y, cd(a[i], 0));
    }
  }
  if (nsts == 1) {
    for (int s = 0; s < nsc; ++s) {
      double g = 0;
      cd cc(0, 0);
      for (int k = 0; k < n_sym; ++k) {
        const size_t q = static_cast<size_t>(k) * nsc + s;
        g += cabs2(d[q]) * hv[q];
        cc += cmul(cmul(Y(k, s), conjd(d[q])), cd(hv[q], 0));
      }
      G[s](0, 0) += cd(g, 0);
      cv[s][0] += cc;
    }
  } else {
    for (int m = 0; m < n_sym / 2; ++m) {
      for (int s = 0; s < nsc; ++s) {
        const size_t q0 = static_cast<size_t>(2 * m) * nsc + s, q1 = static_cast<size_t>(2 * m + 1) * nsc + s;
        const double ok = hv[q0] && hv[q1] ? 1.0 : 0.0;
        const cd d0 = d[q0], d1 = d[q1];
        const cd rows[2][3] = {{Y(2 * m, s), d0, -conjd(d1)}, {Y(2 * m + 1, s), d1, conjd(d0)}};
        for (const auto& rw : rows) {
          const cd A[2] = {cmul(rw[1], cd(ok, 0)), cmul(rw[2], cd(ok, 0))};
          const cd yo = cmul(rw[0], cd(ok, 0));
          for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) G[s](i, j) += cmul(conjd(A[i]), A[j]);
          for (int i = 0; i < 2; ++i) cv[s][i] += cmul(conjd(A[i]), yo);
        }
      }
    }
  }
  for (int s = 0; s < nsc; ++s) {
    if (!c.live_data[s]) continue;
    Mat<cd> A = G[s];
    for (int i = 0; i < nsts; ++i) A(i, i) += cd(1e-9, 0);
    const cdvec sol = np::solve(A, cv[s]);
    for (int i = 0; i < nsts; ++i) out[i][c.data_bins[s]] = sol[i];
  }
  // pilots: known on every symbol
  for (int j = 0; j < 4; ++j) {
    const int b = c.pb[j];
    dmat A(n_sym + nsts, nsts);
    cdvec vals;
    for (int k = 0; k < n_sym; ++k) {
      for (int i = 0; i < nsts; ++i) A(k, i) = ht_pilot(k, nsts, i, j);
      vals.push_back(Zc[k][b]);
    }
    for (int n = 0; n < nsts; ++n) {
      for (int i = 0; i < nsts; ++i) A(n_sym + n, i) = nsts > 1 ? kHtLtfP[i][n] : 1;
      vals.push_back(cdiv(c.Yl[n][b], cd(c.Lh[b], 0)));
    }
    const cdvec x = np::lstsq_x(A, vals);
    for (int i = 0; i < nsts; ++i) out[i][b] = x[i];
  }
  return out;
}

}  // namespace

Ht decode_ht_payload(const ofdm::RxState& rx, int mcs, int length, bool short_gi, bool aggregation, int n_ess,
                     double chan_offset_hz, double usable_half_hz, int stbc, int ldpc_on, int smoothing) {
  const auto dec = ht_decodable(mcs, ldpc_on, stbc, n_ess, 20);
  if (!dec.first) throw std::runtime_error(dec.second);
  if (length > kHtMaxLength)
    throw std::runtime_error("HT length " + std::to_string(length) + " B over the " + std::to_string(kHtMaxLength) +
                             " B decode budget");
  const cdvec& r = *rx.r;
  const int64_t start = rx.start;
  const int back = rx.back;
  HtCtx c;
  const Rate cr = ht_mcs(mcs, &c.nbpsc);
  c.ncbps = 52 * c.nbpsc;
  const int ndbps = c.ncbps * cr.first / cr.second;
  c.nsts = stbc ? 2 : 1;
  const int m_stbc = stbc ? 2 : 1;
  ldpc::Params LP;
  if (ldpc_on) {
    LP = ldpc::params(length, c.ncbps, cr, m_stbc);
    c.n_sym = static_cast<int>(LP.nsym);
  } else {
    c.n_sym = m_stbc * ((16 + 8 * length + 6 + m_stbc * ndbps - 1) / (m_stbc * ndbps));
  }
  const int n_ltf = std::max(1, c.nsts);
  const int64_t ltf1 = 640;
  const int64_t data0 = ltf1 + 80 * (n_ltf + n_ess);
  const int sym = short_gi ? 72 : 80, gi = short_gi ? 8 : 16;
  const int64_t need = data0 + static_cast<int64_t>(sym) * c.n_sym + 96;
  const int64_t nseg = std::max<int64_t>(0, std::min<int64_t>(need, static_cast<int64_t>(r.size()) - start));
  if (nseg < data0 + 64) throw std::runtime_error("PPDU truncated before HT-LTF");
  const double bs = phase_step(rx.cfo, ofdm::kFsWifi);
  cdvec seg(nseg);
  for (int64_t i = 0; i < nseg; ++i) seg[i] = cmul(r[start + i], expj(bs * static_cast<double>(start + i)));

  std::vector<int> good;
  for (int k : ofdm::kHtUsedSc)
    if (std::fabs(chan_offset_hz + k * 312.5e3) <= usable_half_hz - 0.4e6) good.push_back(k);
  for (int k : good) c.used.push_back(bin_of(k));
  std::vector<uint8_t> is_used(64, 0);
  for (int b : c.used) is_used[b] = 1;
  c.Lh.assign(64, 0.0);
  for (int k = -28; k <= 28; ++k) c.Lh[bin_of(k)] = ofdm::kHtLtf[k + 28];
  for (int i = 0; i < c.nsts; ++i) c.Yl.push_back(fft_at(seg, ltf1 + 80 * i + 16, back));
  Hmat H(c.nsts, cdvec(64, cd(0, 0)));
  for (int i = 0; i < c.nsts; ++i)
    for (int b : c.used) {
      cd acc = c.Yl[0][b];
      if (c.nsts > 1) acc = cmul(cd(kHtLtfP[i][0], 0), c.Yl[0][b]) + cmul(cd(kHtLtfP[i][1], 0), c.Yl[1][b]);
      H[i][b] = cdiv(acc, cd(c.nsts * c.Lh[b], 0));
    }
  dvec nkl;
  for (int b : rx.used) nkl.push_back(rx.nk[b]);
  c.base = np::median(nkl);
  dvec nk(64, c.base);
  for (int b : c.used) nk[b] = rx.nk[b] > 0 ? rx.nk[b] : c.base;
  for (int k : ofdm::kHtDataSc) {
    c.data_bins.push_back(bin_of(k));
    c.live_data.push_back(is_used[bin_of(k)]);
  }
  for (int k : ofdm::kPilotSc) c.pb.push_back(bin_of(k));
  for (int k = 0; k < c.n_sym; ++k) c.Ys.push_back(fft_at(seg, data0 + static_cast<int64_t>(sym) * k + gi, back));
  {
    dvec padded;
    for (int i = 60; i < 64; ++i) padded.push_back(nk[i]);
    padded.insert(padded.end(), nk.begin(), nk.end());
    for (int i = 0; i < 4; ++i) padded.push_back(nk[i]);
    dvec ns = box_valid(padded, 9, 1.0 / 9);
    dvec nu;
    for (int b : c.used) nu.push_back(ns[b]);
    const double med = np::median(nu);
    const double den = 1e-20 > med ? 1e-20 : med;
    for (auto& v : ns) v = std::min(std::max(v / den, 0.5), 4.0);
    c.nshape = ns;
  }

  struct Cand {
    double lvl;
    int w;
    Hmat Hs;
    std::vector<cdvec> Zc;
  };
  std::vector<Cand> cands;
  for (int w = 0; w <= (smoothing ? 3 : 0); ++w) {
    Hmat Hs(c.nsts);
    for (int i = 0; i < c.nsts; ++i) Hs[i] = smooth_channel(H[i], good, w);
    Tracked t = tracked(c, Hs);
    cands.push_back({t.lvl, w, std::move(Hs), std::move(t.Zc)});
  }
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.lvl < b.lvl; });

  Ht out;
  out.n_sym = c.n_sym;
  out.ldpc = ldpc_on != 0;
  const size_t nd = c.data_bins.size();
  auto llr_of = [&](const Combined& cb) {
    dvec llr;
    llr.reserve(cb.pts.size() * nd * c.nbpsc);
    for (size_t k = 0; k < cb.pts.size(); ++k) {
      const dvec l = llr_exact(cb.pts[k], c.nbpsc, cb.var[k]);
      for (size_t i = 0; i < nd; ++i)
        for (int t = 0; t < c.nbpsc; ++t) llr.push_back(c.live_data[i] ? l[i * c.nbpsc + t] : 0.0);
    }
    return llr;
  };
  std::vector<uint8_t> databits;
  if (ldpc_on) {
    const auto offs = ldpc::offsets(LP);
    std::map<int64_t, std::vector<uint8_t>> known;
    std::vector<std::vector<uint8_t>> firsts;
    struct Trial {
      double lvl;
      const Hmat* Hs;
      const std::vector<cdvec>* Zc;
    };
    std::vector<Trial> trials;
    for (const auto& cd_ : cands) trials.push_back({cd_.lvl, &cd_.Hs, &cd_.Zc});
    Hmat Hda;
    Tracked t2;
    dvec llr;
    const int64_t k0 = LP.K - ldpc::split(LP.nshrt, LP.ncw)[0];
    for (int rnd = 0; rnd < 8; ++rnd) {
      int nnew = 0;
      for (const auto& tr : trials) {
        llr = llr_of(combine(c, *tr.Hs, *tr.Zc, tr.lvl));
        auto sd = ldpc::decode_stream(llr, LP, cr);
        std::vector<uint8_t> f(static_cast<size_t>(std::min<int64_t>(k0, static_cast<int64_t>(llr.size()))));
        for (size_t i = 0; i < f.size(); ++i) f[i] = llr[i] < 0;
        firsts.push_back(std::move(f));
        for (size_t i = 0; i < sd.ok.size(); ++i)
          if (sd.ok[i] && !known.count(static_cast<int64_t>(i))) {
            known[static_cast<int64_t>(i)] = sd.hard[i];
            ++nnew;
          }
        if (static_cast<int64_t>(known.size()) == LP.ncw) break;
      }
      if (static_cast<int64_t>(known.size()) == LP.ncw || !nnew || rnd == 7) break;
      Hda = data_aided_channel(c, known, LP, offs, cands[0].Zc, cands[0].Hs);
      t2 = tracked(c, Hda);
      trials = {{t2.lvl, &Hda, &t2.Zc}};
    }
    const auto lay = ldpc::layout(LP);
    for (size_t i = 0; i < lay.size(); ++i) {
      const int64_t k_i = lay[i].k;
      auto it = known.find(static_cast<int64_t>(i));
      if (it != known.end()) {
        databits.insert(databits.end(), it->second.begin(), it->second.begin() + k_i);
      } else if (i == 0) {
        const std::vector<uint8_t>* pick = &firsts[0];
        for (const auto& f : firsts) {
          const auto d = descramble(std::vector<uint8_t>(f.begin(), f.begin() + std::min<size_t>(16, f.size()))).first;
          bool any = false;
          for (size_t q = 7; q < 16 && q < d.size(); ++q) any |= d[q] != 0;
          if (!any) {
            pick = &f;
            break;
          }
        }
        databits.insert(databits.end(), pick->begin(), pick->end());
      } else {
        const int64_t lo = offs[i];
        for (int64_t q = lo; q < lo + k_i && q < static_cast<int64_t>(llr.size()); ++q) databits.push_back(llr[q] < 0);
      }
    }
    out.n_codewords = static_cast<int>(LP.ncw);
    out.n_codewords_ok = static_cast<int>(known.size());
  } else {
    const Combined cb = combine(c, cands[0].Hs, cands[0].Zc, cands[0].lvl);
    const auto inv = deinterleave_indices(c.ncbps, c.nbpsc, 13).second;
    dvec soft;
    soft.reserve(cb.pts.size() * c.ncbps);
    dvec coded(c.ncbps);
    for (size_t k = 0; k < cb.pts.size(); ++k) {
      const dvec l = llr_exact(cb.pts[k], c.nbpsc, cb.var[k]);
      for (size_t i = 0; i < nd; ++i)
        for (int t = 0; t < c.nbpsc; ++t) {
          const size_t q = i * c.nbpsc + t;
          coded[inv[q]] = c.live_data[i] ? -l[q] : 0.0;
        }
      soft.insert(soft.end(), coded.begin(), coded.end());
    }
    databits = ofdm::viterbi(depuncture(soft, cr));
  }
  auto [descr, seed] = descramble(databits);
  out.scrambler_seed = seed_int(seed);
  bool any = false;
  for (size_t i = 7; i < 16 && i < descr.size(); ++i) any |= descr[i] != 0;
  out.service_ok = descr.size() >= 16 && !any;
  const std::vector<uint8_t> data = pack_lsb(descr, 16, static_cast<size_t>(length));
  std::vector<std::pair<std::vector<uint8_t>, bool>> parts;
  if (aggregation) parts = deaggregate(data);
  else parts.push_back({data, true});
  for (const auto& [mpdu, delim_ok] : parts) {
    (void)delim_ok;
    if (mpdu.size() < 4) continue;
    const size_t nb = mpdu.size() - 4;
    const uint32_t fcs = mpdu[nb] | (mpdu[nb + 1] << 8) | (mpdu[nb + 2] << 16) | (static_cast<uint32_t>(mpdu[nb + 3]) << 24);
    HtSubframe sf;
    sf.fcs_ok = mac80211::crc32(mpdu.data(), nb) == fcs;
    if (sf.fcs_ok) sf.mac = mac80211::parse(mpdu.data(), nb);
    if (sf.fcs_ok) ++out.n_fcs_ok;
    out.subframes.push_back(std::move(sf));
  }
  return out;
}

}  // namespace si::rx::wifi::payload
