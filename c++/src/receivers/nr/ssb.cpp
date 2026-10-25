#include "receivers/nr/ssb.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/seq3gpp.hpp"
#include "receivers/nr/bch.hpp"

namespace si::rx::nr {

namespace {

constexpr double kTwoPi = 2.0 * kPi;

inline int64_t pymod(int64_t a, int64_t b) {
  int64_t r = a % b;
  return r < 0 ? r + b : r;
}

// np.convolve(|y|^2, ones(N), "valid")[:L-N] + 1e-30 (running float64 sum re-anchored every 1024)
dvec window_energy(const cfvec& y, size_t N) {
  const size_t L = y.size();
  dvec out;
  if (L <= N) return out;
  fvec p(L);
  for (size_t i = 0; i < L; ++i) {
    const float a = np::cabs(y[i]);
    p[i] = a * a;
  }
  out.resize(L - N);
  double s = 0;
  for (size_t i = 0; i < out.size(); ++i) {
    if (i % 1024 == 0) {
      s = 0;
      for (size_t j = 0; j < N; ++j) s += static_cast<double>(p[i + j]);
    } else {
      s = s + static_cast<double>(p[i + N - 1]) - static_cast<double>(p[i - 1]);
    }
    out[i] = s + 1e-30;
  }
  return out;
}

// 7-bit m-sequence x(i+7) = x(i+a) + x(i) mod 2 from x(0..6) = init
std::array<int, 127> mseq7(const std::array<int, 7>& init, int a) {
  std::array<int, 127> x{};
  for (int i = 0; i < 7; ++i) x[i] = init[i];
  for (int i = 0; i < 120; ++i) x[i + 7] = (x[i + a] + x[i]) % 2;
  return x;
}

// z-score of the maximum of m: (m[k] - med) / (1.4826 * MAD + 1e-12)
struct ZK {
  double z;
  int64_t k;
};
ZK zscore(const dvec& m, dvec& buf, dvec& scratch) {
  const size_t k = np::argmax(m.data(), m.size());
  const double med = np::median_buf(m.data(), m.size(), buf);
  scratch.resize(m.size());
  for (size_t i = 0; i < m.size(); ++i) scratch[i] = std::fabs(m[i] - med);
  const double mad = 1.4826 * np::median_inplace(scratch.data(), scratch.size()) + 1e-12;
  return {(m[k] - med) / mad, static_cast<int64_t>(k)};
}

// SSS differential bank s_k s_{k+1} for every N_ID1 (python _sss_bank)
const std::array<std::vector<std::array<double, 126>>, 3>& sss_diff_bank() {
  static const auto bank = [] {
    std::array<std::vector<std::array<double, 126>>, 3> b;
    for (int n2 = 0; n2 < 3; ++n2) {
      b[n2].resize(336);
      for (int n1 = 0; n1 < 336; ++n1) {
        const dvec s = sss_seq(3 * n1 + n2);
        for (int i = 0; i < 126; ++i) b[n2][n1][i] = s[i] * s[i + 1];
      }
    }
    return b;
  }();
  return bank;
}

int bits_to_int(const Bits& b, size_t i, size_t n) {
  int v = 0;
  for (size_t k = i; k < i + n; ++k) v = (v << 1) | (b[k] & 1);
  return v;
}

// python _parse_mib
Json parse_mib(const BchResult& r) {
  const Bits& b = r.payload;
  const int sfn_msb = bits_to_int(b, 1, 6);
  const int lsb = bits_to_int(r.sfn_lsb, 0, 4);
  const int k_ssb_lsb = bits_to_int(b, 8, 4);
  const int k_msb = 0;
  const int pdcch = bits_to_int(b, 13, 8);
  Json m = Json::object();
  m["message"] = b[0] == 0 ? "BCCH-BCH mib" : "messageClassExtension";
  m["systemFrameNumber"] = (sfn_msb << 4) | lsb;
  m["subCarrierSpacingCommon_khz"] = b[7] == 0 ? 15 : 30;
  m["ssb_SubcarrierOffset"] = k_ssb_lsb + 16 * k_msb;
  m["dmrs_TypeA_Position"] = b[12] == 0 ? 2 : 3;
  m["pdcch_ConfigSIB1"] = Json{{"controlResourceSetZero", pdcch >> 4}, {"searchSpaceZero", pdcch & 15}};
  m["cellBarred"] = b[21] == 0 ? "barred" : "notBarred";
  m["intraFreqReselection"] = b[22] == 0 ? "allowed" : "notAllowed";
  m["half_frame"] = r.hrf;
  m["sib1_scheduled"] = k_ssb_lsb + 16 * k_msb <= 23;
  return m;
}

}  // namespace

const std::vector<NrBand>& nr_bands() {
  static const std::vector<NrBand> b = {{"n78", "n78", 30e3},   {"b28_dl", "n28", 15e3}, {"b40", "n40", 30e3},
                                        {"b41", "n41", 30e3},   {"b3_dl", "n3", 15e3},   {"b1_dl", "n1", 15e3},
                                        {"b5_dl", "n5", 15e3},  {"b8_dl", "n8", 15e3}};
  return b;
}

std::vector<std::pair<double, double>> bands_for(const std::string& key) {
  std::vector<std::pair<double, double>> out;
  for (const auto& b : bands::bands())
    if (b.key == key) out.emplace_back(b.lo, b.hi);
  return out;
}

dvec gscn(double lo, double hi) {
  dvec out;
  if (hi > 3000e6) {
    const int64_t n0 = static_cast<int64_t>(std::ceil((std::max(lo, 3000e6) - 3000e6) / 1.44e6));
    const int64_t n1 = static_cast<int64_t>((hi - 3000e6) / 1.44e6);
    for (int64_t n = n0; n <= n1; ++n) out.push_back(3000e6 + static_cast<double>(n) * 1.44e6);
  }
  if (lo < 3000e6) {
    const int64_t a = static_cast<int64_t>(lo / 1.2e6), b = static_cast<int64_t>(std::min(hi, 3000e6) / 1.2e6);
    for (int64_t n = a; n <= b; ++n)
      for (double m : {50e3, 150e3, 250e3}) {
        const double f = static_cast<double>(n) * 1.2e6 + m;
        if (lo <= f && f <= hi) out.push_back(f);
      }
  }
  return out;
}

dvec pss_seq(int nid2) {
  const auto x = mseq7({0, 1, 1, 0, 1, 1, 1}, 4);
  dvec d(127);
  for (int n = 0; n < 127; ++n) d[n] = 1 - 2 * x[(n + 43 * nid2) % 127];
  return d;
}

dvec sss_seq(int pci) {
  const auto x0 = mseq7({1, 0, 0, 0, 0, 0, 0}, 4);
  const auto x1 = mseq7({1, 0, 0, 0, 0, 0, 0}, 1);
  const int n1 = pci / 3, n2 = pci % 3;
  const int m0 = 15 * (n1 / 112) + 5 * n2, m1 = n1 % 112;
  dvec d(127);
  for (int n = 0; n < 127; ++n) d[n] = (1 - 2 * x0[(n + m0) % 127]) * (1 - 2 * x1[(n + m1) % 127]);
  return d;
}

cdvec pbch_dmrs(int pci, int ib) {
  const uint32_t c_init = (1u << 11) * (ib + 1) * (pci / 4 + 1) + (1u << 6) * (ib + 1) + (pci % 4);
  return seq3gpp::qpsk(seq3gpp::gold(c_init, 288));
}

std::shared_ptr<const TemplateBank> template_bank(double scs, double ofs, size_t nf) {
  auto b = std::make_shared<TemplateBank>();
  b->nf = nf;
  b->ofs = ofs;
  b->T = Mat<cd>(3, nf);
  // python _pss_templates(scs, nf) uses t = n / (NFFT * scs), _templates() t = n / ofs: equal when
  // the slice rate is exactly NFFT * scs
  const double tden = ofs == kNfft * scs ? kNfft * scs : ofs;
  std::vector<cd> row(127);
  for (int nid2 = 0; nid2 < 3; ++nid2) {
    const dvec p = pss_seq(nid2);
    cdvec x(kNfft);
    for (int t = 0; t < kNfft; ++t) {
      const double tt = static_cast<double>(t) / tden;
      for (int k = 0; k < 127; ++k) {
        const double ph = ((kTwoPi * (k - 64)) * scs) * tt;
        row[k] = np::cmul(cd(p[k], 0.0), cd(std::cos(ph), std::sin(ph)));
      }
      x[t] = np::pairwise_sum(row.data(), row.size());
    }
    const double nrm = np::norm(x);
    for (auto& v : x) v = np::cdiv(v, cd(nrm, 0.0));
    const cdvec F = fft::fft(x, nf);
    for (size_t i = 0; i < nf; ++i) b->T(nid2, i) = std::conj(F[i]);
  }
  return b;
}

SliceDims slice_dims(const Channelizer& ch, double out_fs) {
  const int64_t M = std::max<int64_t>(16, py_round_int(static_cast<double>(ch.L()) * out_fs / ch.fs()));
  SliceDims d;
  d.L = static_cast<size_t>(M);
  d.ofs = static_cast<double>(M) * ch.fs() / static_cast<double>(ch.L());
  d.nf = fft::next_fast_len(d.L);   // circular correlation: the lags used (< L - NFFT) never wrap (PLAN P1)
  return d;
}

SSBSlice::SSBSlice(const Channelizer& ch, double f_ss, double scs_, std::shared_ptr<const TemplateBank> bank)
    : scs(scs_) {
  cp = static_cast<int>(py_round_int(144.0 / 2048.0 * kNfft));
  auto s = ch.slice(f_ss, 0.55 * 240 * scs, kNfft * scs);
  y = std::move(s.y);
  ofs = s.fs;
  L = y.size();
  nf = fft::next_fast_len(L);
  Y_ = fft::np_fft(y, nf);
  e_ = window_energy(y, kNfft);
  bank_ = (bank && bank->nf == nf && bank->ofs == ofs) ? std::move(bank) : template_bank(scs, ofs, nf);
}

dvec SSBSlice::metric(double df, int nid2) const {
  cdvec w;
  dvec m;
  metric_into(df, nid2, w, m);
  return m;
}

void SSBSlice::metric_into(double df, int nid2, cdvec& w, dvec& m) const {
  const int64_t sh = py_round_int(df / ofs * static_cast<double>(nf));
  const int64_t NF = static_cast<int64_t>(nf);
  w.resize(nf);
  const cd* T = bank_->T.row(nid2);
  // np.roll(Y, -sh): element i is Y[(i + sh) mod nf]
  const size_t s0 = static_cast<size_t>(pymod(sh, NF));
  for (size_t i = 0, j = s0; i < nf; ++i, j = (j + 1 == nf ? 0 : j + 1)) {
    const cf v = Y_[j];
    w[i] = np::cmul(cd(v.real(), v.imag()), T[i]);
  }
  fft::c2c(w.data(), w.data(), nf, false, 1.0 / static_cast<double>(nf));
  const size_t n = L - kNfft;
  m.resize(n);
  for (size_t i = 0; i < n; ++i) {
    const double a = np::cabs(w[i]);
    m[i] = a * a / e_[i];
  }
}

PssBest SSBSlice::best(const dvec& offsets) const {
  PssBest b;
  bool have = false;
  dvec scratch, m(L - kNfft);
  cdvec yr(nf);
  dvec ms;
  Mat<cd> W(3, nf);
  const int64_t NF = static_cast<int64_t>(nf);
  for (double df : offsets) {
    const int64_t sh = py_round_int(df / ofs * static_cast<double>(nf));
    // scipy.fft.ifft(roll(Y, -sh)[None, :] * T, axis=1): the three rows in one batched transform
    const size_t s0 = static_cast<size_t>(pymod(sh, NF));
    // roll(Y, -sh) in complex128 once, then the three template products
    for (size_t i = 0, j = s0; i < nf; ++i, j = (j + 1 == nf ? 0 : j + 1)) yr[i] = cd(Y_[j].real(), Y_[j].imag());
    for (int r = 0; r < 3; ++r) np::cmul_n(yr.data(), bank_->T.row(r), W.row(r), nf);
    fft::fft_rows(W, false, 1);
    for (int nid2 = 0; nid2 < 3; ++nid2) {
      const cd* w = W.row(nid2);
      for (size_t i = 0; i < m.size(); ++i) {
        const double a = np::cabs(w[i]);
        m[i] = a * a / e_[i];
      }
      const ZK zk = zscore(m, ms, scratch);
      if (!have || zk.z > b.z) {
        have = true;
        b.z = zk.z;
        b.nid2 = nid2;
        b.pos = zk.k;
        b.df = static_cast<double>(sh) * ofs / static_cast<double>(nf);
        b.metric = m[zk.k];
      }
    }
  }
  return b;
}

cdvec SSBSlice::shifted(double df) const {
  const double a = -kTwoPi * df;
  const double inv = 1.0 / ofs;
  cdvec out(L);
  for (size_t n = 0; n < L; ++n) {
    const double ph = (a * static_cast<double>(n)) * inv;
    out[n] = np::cmul(cd(y[n].real(), y[n].imag()), cd(std::cos(ph), std::sin(ph)));
  }
  return out;
}

std::optional<Mat<cd>> SSBSlice::grid(const cdvec& sy, int64_t pos) const {
  Mat<cd> G(4, 240);
  for (int l = 0; l < 4; ++l) {
    const int64_t s = pos + l * (kNfft + cp);
    if (s + kNfft > static_cast<int64_t>(L)) return std::nullopt;
    const cdvec X = fft::fft(sy.data() + s, kNfft);
    for (int k = 0; k < 240; ++k) G(l, k) = X[pymod(k - 120, kNfft)];
  }
  return G;
}

std::vector<int64_t> peaks(const dvec& m, double frac, int64_t sep, size_t max_n) {
  if (m.empty()) return {};
  const double top = m[np::argmax(m.data(), m.size())];
  // argsort(m)[::-1]: descending; the loop stops at the first value below frac * top
  std::vector<int64_t> cand;
  for (size_t i = 0; i < m.size(); ++i)
    if (!(m[i] < frac * top)) cand.push_back(static_cast<int64_t>(i));
  std::sort(cand.begin(), cand.end(), [&](int64_t a, int64_t b) { return m[a] != m[b] ? m[a] > m[b] : a > b; });
  if (cand.size() > 20000) cand.resize(20000);
  std::vector<int64_t> out;
  for (int64_t i : cand) {
    if (out.size() >= max_n) break;
    bool ok = true;
    for (int64_t j : out)
      if (!(std::llabs(i - j) > sep)) ok = false;
    if (ok) out.push_back(i);
  }
  std::sort(out.begin(), out.end());
  return out;
}

dvec sss_scores(const Mat<cd>& G, int nid2) {
  std::array<double, 126> dr;
  for (int i = 0; i < 126; ++i) dr[i] = np::cmul(G(2, 56 + i), std::conj(G(2, 57 + i))).real();
  const auto& bank = sss_diff_bank()[nid2];
  dvec out(336);
  for (int n1 = 0; n1 < 336; ++n1) {
    double s = 0;
    for (int i = 0; i < 126; ++i) s += bank[n1][i] * dr[i];
    out[n1] = s;
  }
  return out;
}

std::pair<int, dvec> pbch_llr(const Mat<cd>& G, int pci) {
  const int v = pci % 4;
  std::vector<int> k13, k2, d13, d2;
  for (int k = v; k < 240; k += 4) k13.push_back(k);
  for (int k = v; k < 48; k += 4) k2.push_back(k);
  for (int k = 192 + v; k < 240; k += 4) k2.push_back(k);
  for (int k = 0; k < 240; ++k)
    if (k % 4 != v) d13.push_back(k);
  for (int k = 0; k < 240; ++k)
    if ((k < 48 || k >= 192) && k % 4 != v) d2.push_back(k);
  cdvec rx;
  for (int k : k13) rx.push_back(G(1, k));
  for (int k : k2) rx.push_back(G(2, k));
  for (int k : k13) rx.push_back(G(3, k));
  double best_met = 0;
  int best_ib = -1;
  cdvec h(rx.size()), hb, prod(rx.size() - 1);
  for (int ib = 0; ib < 8; ++ib) {
    const cdvec d = pbch_dmrs(pci, ib);
    for (size_t i = 0; i < rx.size(); ++i) h[i] = np::cmul(rx[i], std::conj(d[i]));
    for (size_t i = 0; i + 1 < h.size(); ++i) prod[i] = np::cmul(h[i + 1], std::conj(h[i]));
    const double met = np::cabs(np::pairwise_sum(prod.data(), prod.size()));
    if (best_ib < 0 || met > best_met) {
      best_met = met;
      best_ib = ib;
      hb = h;
    }
  }
  const dvec w(5, 1.0 / 5);
  auto smooth = [&](size_t a, size_t n) {
    dvec re(n), im(n);
    for (size_t i = 0; i < n; ++i) {
      re[i] = hb[a + i].real();
      im[i] = hb[a + i].imag();
    }
    return std::make_pair(np::convolve(re, w, np::ConvMode::Same), np::convolve(im, w, np::ConvMode::Same));
  };
  const auto h1 = smooth(0, 60), h3 = smooth(84, 60);
  dvec kk(240), xp(k13.begin(), k13.end());
  for (int i = 0; i < 240; ++i) kk[i] = i;
  const dvec H1r = np::interp(kk, xp, h1.first), H1i = np::interp(kk, xp, h1.second);
  const dvec H3r = np::interp(kk, xp, h3.first), H3i = np::interp(kk, xp, h3.second);
  auto H1 = [&](int k) { return cd(H1r[k], H1i[k]); };
  auto H3 = [&](int k) { return cd(H3r[k], H3i[k]); };
  auto H2 = [&](int k) { return np::cmul(cd(0.5, 0.0), H1(k) + H3(k)); };
  cdvec Hd, rxd;
  for (int k : d13) Hd.push_back(H1(k)), rxd.push_back(G(1, k));
  for (int k : d2) Hd.push_back(H2(k)), rxd.push_back(G(2, k));
  for (int k : d13) Hd.push_back(H3(k)), rxd.push_back(G(3, k));
  dvec p2(Hd.size());
  for (size_t i = 0; i < Hd.size(); ++i) {
    const double a = np::cabs(Hd[i]);
    p2[i] = a * a;
  }
  const double den = np::mean(p2.data(), p2.size()) + 1e-12;
  dvec llr(2 * Hd.size());
  for (size_t i = 0; i < Hd.size(); ++i) {
    const cd z = np::cmul(rxd[i], std::conj(Hd[i]));
    llr[2 * i] = z.real() / den;
    llr[2 * i + 1] = z.imag() / den;
  }
  return {best_ib, llr};
}

std::pair<std::optional<int>, SsbInfo> decode_ssb_multi(const SSBSlice& sl, const PssBest& b) {
  const dvec m = sl.metric(b.df, b.nid2);
  const std::vector<int64_t> inst = peaks(m);
  const cdvec sy = sl.shifted(b.df);
  std::vector<std::pair<int64_t, Mat<cd>>> grids;
  for (int64_t p : inst) {
    auto G = sl.grid(sy, p);
    if (G) grids.emplace_back(p, std::move(*G));
  }
  SsbInfo info;
  info.n_ssb_instances = grids.size();
  if (grids.empty()) return {std::nullopt, info};
  // (N_ID2, N_ID1) chosen jointly
  struct C {
    int nid2, nid1;
    double z, z2;
  };
  std::vector<C> cells;
  dvec scratch;
  for (int n2 = 0; n2 < 3; ++n2) {
    dvec acc(336, 0.0);
    for (const auto& [p, G] : grids) {
      const dvec s = sss_scores(G, n2);
      for (int i = 0; i < 336; ++i) acc[i] += s[i];
    }
    // argsort(acc)[::-1]: o[0], o[1] = the two largest (ties: the later index first)
    std::vector<int> o(336);
    for (int i = 0; i < 336; ++i) o[i] = i;
    std::partial_sort(o.begin(), o.begin() + 2, o.end(),
                      [&](int a, int c) { return acc[a] != acc[c] ? acc[a] > acc[c] : a > c; });
    scratch.assign(acc.begin(), acc.end());
    const double med = np::median_inplace(scratch.data(), scratch.size());
    for (int i = 0; i < 336; ++i) scratch[i] = std::fabs(acc[i] - med);
    const double mad = 1.4826 * np::median_inplace(scratch.data(), scratch.size()) + 1e-12;
    cells.push_back({n2, o[0], (acc[o[0]] - med) / mad, (acc[o[1]] - med) / mad});
  }
  std::stable_sort(cells.begin(), cells.end(), [](const C& a, const C& c) { return -a.z < -c.z; });
  const C& bc = cells[0];
  info.nid1 = bc.nid1;
  info.nid2 = bc.nid2;
  info.sss_z = bc.z;
  info.sss_z2 = bc.z2;
  for (size_t i = 1; i < cells.size(); ++i)
    if (cells[i].z >= 4.5 && cells[i].z / std::max(cells[i].z2, 1e-3) >= 1.3)
      info.other_cells.push_back({3 * cells[i].nid1 + cells[i].nid2, py_round(cells[i].z, 1)});
  const int pci = 3 * bc.nid1 + bc.nid2;
  struct Per {
    int64_t p;
    int ib;
    dvec llr;
  };
  std::vector<Per> per;
  for (const auto& [p, G] : grids) {
    auto [ib, llr] = pbch_llr(G, pci);
    per.push_back({p, ib, std::move(llr)});
  }
  for (const auto& x : per) info.ssb_indices.push_back(x.ib);
  {
    size_t bi = 0;
    for (size_t i = 1; i < per.size(); ++i)
      if (m[per[i].p] > m[per[bi].p]) bi = i;
    info.ibar_ssb = per[bi].ib;
  }
  const int64_t half = static_cast<int64_t>(5e-3 * sl.ofs);
  std::vector<std::pair<int64_t, std::vector<const Per*>>> groups;  // insertion order
  for (const auto& x : per) {
    const int64_t key = x.p / half;
    auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return g.first == key; });
    if (it == groups.end()) {
      groups.push_back({key, {}});
      it = groups.end() - 1;
    }
    it->second.push_back(&x);
  }
  std::stable_sort(groups.begin(), groups.end(),
                   [](const auto& a, const auto& c) { return -static_cast<int64_t>(a.second.size()) <
                                                             -static_cast<int64_t>(c.second.size()); });
  const std::vector<int8_t> prbs = seq3gpp::gold(static_cast<uint32_t>(pci), 8 * 864);
  int tried = 0;
  for (int lmax : {8, 4}) {
    for (const auto& [key, members] : groups) {
      dvec soft(864, 0.0);
      for (const Per* x : members) {
        const int v = x->ib % lmax;
        for (int i = 0; i < 864; ++i) soft[i] += x->llr[i] * (1.0 - 2.0 * prbs[v * 864 + i]);
      }
      const double sd = std::sqrt(np::var(soft.data(), soft.size())) + 1e-12;
      for (auto& s : soft) s = s / sd * 3;
      for (double sign : {1.0, -1.0}) {
        ++tried;
        dvec in(864);
        for (int i = 0; i < 864; ++i) in[i] = sign * soft[i];
        const BchResult r = bch_decode(in, pci);
        if (r.crc_ok) {
          info.crc_ok = true;
          info.lmax = lmax;
          info.n_combined = static_cast<int>(members.size());
          info.mib = parse_mib(r);
          return {pci, info};
        }
      }
    }
  }
  info.crc_ok = false;
  info.decode_attempts = tried;
  return {pci, info};
}

}  // namespace si::rx::nr
