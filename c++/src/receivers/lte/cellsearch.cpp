#include "receivers/lte/cellsearch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/seq3gpp.hpp"

namespace si::rx::lte {

namespace {

constexpr double kTwoPi = 2.0 * kPi;

// python: c = concatenate(([0], cumsum(abs(y.astype(complex128)) ** 2))); c[N:L] - c[:L - N] + 1e-30
dvec window_energy(const cfvec& y, size_t N) {
  const size_t L = y.size();
  dvec out;
  if (L <= N) return out;
  dvec c(L + 1);
  c[0] = 0.0;
  for (size_t i = 0; i < L; ++i) {
    const double a = np::cabs(cd(y[i].real(), y[i].imag()));
    c[i + 1] = c[i] + a * a;
  }
  out.resize(L - N);
  for (size_t i = 0; i < out.size(); ++i) out[i] = (c[N + i] - c[i]) + 1e-30;
  return out;
}

// python _smooth_half(h): convolve(h, ones(n), "same") / convolve(ones(len(h)), ones(n), "same")
void smooth_half(const cd* h, size_t len, cd* out) {
  const int n = kSssSmooth;
  const int64_t L = static_cast<int64_t>(len), off = (n - 1) / 2;
  for (int64_t i = 0; i < L; ++i) {
    // numpy correlate of the reversed (all-ones) kernel: terms in ascending input index
    double re = 0, im = 0, cnt = 0;
    for (int64_t k = 0; k < n; ++k) {
      const int64_t j = i - off + k;
      if (j < 0 || j >= L) continue;
      re += h[j].real();
      im += h[j].imag();
      cnt += 1.0;
    }
    out[i] = np::cdiv(cd(re, im), cd(cnt, 0.0));
  }
}

// np.median of a copy
double median_copy(const dvec& x, dvec& scratch) {
  scratch.assign(x.begin(), x.end());
  return np::median_inplace(scratch.data(), scratch.size());
}

double mad_of(const dvec& x, double med, dvec& scratch) {
  scratch.resize(x.size());
  for (size_t i = 0; i < x.size(); ++i) scratch[i] = std::fabs(x[i] - med);
  return 1.4826 * np::median_inplace(scratch.data(), scratch.size()) + 1e-12;
}

// Python int mod (non-negative for a positive divisor)
inline int64_t pymod(int64_t a, int64_t b) {
  int64_t r = a % b;
  return r < 0 ? r + b : r;
}

// python _templates(): the three PSS time-domain templates at 3.84 MSps, unit energy
const Mat<cd>& time_templates() {
  static const Mat<cd> tpl = [] {
    Mat<cd> out(3, kN);
    const auto& sc = seq3gpp::sc62();
    std::vector<cd> row(62);
    for (int r = 0; r < 3; ++r) {
      const cdvec d = seq3gpp::lte_pss_freq(kRoots[r]);
      cdvec x(kN);
      for (int t = 0; t < kN; ++t) {
        const double tt = static_cast<double>(t) / kOfs;
        for (int k = 0; k < 62; ++k) {
          const double ph = ((kTwoPi * sc[k]) * 15e3) * tt;
          row[k] = np::cmul(d[k], cd(std::cos(ph), std::sin(ph)));
        }
        x[t] = np::pairwise_sum(row.data(), row.size());
      }
      const double nrm = np::norm(x);
      for (int t = 0; t < kN; ++t) out(r, t) = np::cdiv(x[t], cd(nrm, 0.0));
    }
    return out;
  }();
  return tpl;
}

}  // namespace

const std::vector<std::pair<std::string, std::string>>& lte_bands() {
  static const std::vector<std::pair<std::string, std::string>> b = {
      {"b28_dl", "B28"}, {"b5_dl", "B5"}, {"b8_dl", "B8"}, {"b3_dl", "B3"},
      {"b1_dl", "B1"},   {"b40", "B40"},  {"b41", "B41"}};
  return b;
}

const std::vector<double>& bws() {
  static const std::vector<double> b = {1.4e6, 3e6, 5e6, 10e6, 15e6, 20e6};
  return b;
}

double occ(double bw) {
  static const std::vector<std::pair<double, double>> o = {{1.4e6, 1.08e6}, {3e6, 2.7e6},  {5e6, 4.5e6},
                                                           {10e6, 9e6},     {15e6, 13.5e6}, {20e6, 18e6}};
  for (const auto& [k, v] : o)
    if (k == bw) return v;
  return 0.0;
}

std::shared_ptr<const TemplateBank> template_bank(size_t nf) {
  auto b = std::make_shared<TemplateBank>();
  b->nf = nf;
  b->T = Mat<cd>(3, nf);
  const Mat<cd>& tpl = time_templates();
  for (int r = 0; r < 3; ++r) {
    const cdvec F = fft::fft(tpl.row(r), kN, nf);
    for (size_t i = 0; i < nf; ++i) b->T(r, i) = std::conj(F[i]);
  }
  return b;
}

Slice::Slice(const Channelizer& ch, double fc, std::shared_ptr<const TemplateBank> bank) : f_center(fc) {
  auto s = ch.slice(fc, 1.0e6 + 60e3, kOfs);
  y = std::move(s.y);
  ofs = s.fs;
  L = y.size();
  nf = fft::next_fast_len(L + kN);  // >= L + N - 1: circular = linear correlation
  cdvec yd(L);
  for (size_t i = 0; i < L; ++i) yd[i] = cd(y[i].real(), y[i].imag());
  Y_ = fft::fft(yd.data(), L, nf);
  bank_ = (bank && bank->nf == nf) ? std::move(bank) : template_bank(nf);
  e_ = window_energy(y, kN);
  k = L > static_cast<size_t>(kN) ? static_cast<int64_t>((L - kN) / kHalf) : 0;
}

Slice::Score Slice::score(double df, int r, cdvec& work, dvec& m) const {
  const int64_t sh = py_round_int(df / ofs * static_cast<double>(nf));
  const int64_t NF = static_cast<int64_t>(nf);
  work.resize(nf);
  const cd* T = bank_->T.row(r);
  // np.roll(Y, -sh): element i is Y[(i + sh) mod nf]
  const size_t s0 = static_cast<size_t>(pymod(sh, NF));
  for (size_t i = 0, j = s0; i < nf; ++i, j = (j + 1 == nf ? 0 : j + 1)) work[i] = np::cmul(Y_[j], T[i]);
  fft::c2c(work.data(), work.data(), nf, false, 1.0 / static_cast<double>(nf));
  const size_t n = L - kN;
  m.resize(n);
  for (size_t i = 0; i < n; ++i) {
    const double a = np::cabs(work[i]);
    m[i] = a * a / e_[i];
  }
  // fold over 5 ms (rows added in order, as numpy's sum(axis=0))
  dvec fold(m.begin(), m.begin() + kHalf);
  for (int64_t j = 1; j < k; ++j)
    for (int64_t i = 0; i < kHalf; ++i) fold[i] += m[j * kHalf + i];
  const int64_t kk = static_cast<int64_t>(np::argmax(fold.data(), fold.size()));
  dvec scratch;
  const double med = median_copy(fold, scratch);
  const double mad = mad_of(fold, med, scratch);
  double far_max = -std::numeric_limits<double>::infinity();
  for (int64_t i = 0; i < kHalf; ++i)
    if (std::llabs(pymod(i - kk + kHalf / 2, kHalf) - kHalf / 2) > kN) far_max = std::max(far_max, fold[i]);
  Score s;
  s.iso = (fold[kk] - med) / (far_max - med + 1e-12);
  s.z = (fold[kk] - med) / mad;
  s.k = kk;
  s.dfq = static_cast<double>(sh) * ofs / static_cast<double>(nf);
  return s;
}

PssBest Slice::best(const dvec& offsets, bool want_m, bool parallel) const {
  const size_t H = offsets.size() * 3;
  std::vector<Score> sc(H);
  auto run = [&](size_t h) {
    cdvec work;
    dvec m;
    sc[h] = score(offsets[h / 3], static_cast<int>(h % 3), work, m);
  };
  if (parallel) {
    parallel_for(H, run);
  } else {
    for (size_t h = 0; h < H; ++h) run(h);
  }
  size_t bi = 0;
  for (size_t h = 1; h < H; ++h)
    if (sc[h].z > sc[bi].z) bi = h;
  PssBest b;
  b.z = sc[bi].z;
  b.df = sc[bi].dfq;
  b.nid2 = static_cast<int>(bi % 3);
  b.phase = sc[bi].k;
  b.iso = sc[bi].iso;
  if (want_m) {
    cdvec work;
    score(offsets[bi / 3], b.nid2, work, b.m);
  }
  return b;
}

cdvec Slice::shifted(double df) const {
  // -2j*pi*df*arange(L)/ofs: phase ((-2 pi df) * n) * (1/ofs) (numpy complex / real)
  const double a = -kTwoPi * df;
  const double inv = 1.0 / ofs;
  cdvec out(L);
  for (size_t n = 0; n < L; ++n) {
    const double ph = (a * static_cast<double>(n)) * inv;
    out[n] = np::cmul(cd(y[n].real(), y[n].imag()), cd(std::cos(ph), std::sin(ph)));
  }
  return out;
}

Instances instances(const dvec& m, int64_t phase) {
  const int64_t n = static_cast<int64_t>(m.size());
  double pk = -std::numeric_limits<double>::infinity();
  for (int64_t p = phase; p < n; p += kHalf) pk = std::max(pk, m[p]);
  Instances out;
  int64_t j = 0;
  for (int64_t p = phase; p < n; p += kHalf, ++j) {
    double mx = -std::numeric_limits<double>::infinity();
    for (int64_t i = std::max<int64_t>(0, p - 2); i < std::min(n, p + 3); ++i) mx = std::max(mx, m[i]);
    if (mx > 0.3 * pk) out.emplace_back(j, p);
  }
  return out;
}

std::optional<SssResult> sss_decode(const cdvec& y, int nid2, const Instances& inst) {
  const auto& sc62 = seq3gpp::sc62();
  std::vector<int> sc(62);
  for (int i = 0; i < 62; ++i) sc[i] = static_cast<int>(pymod(sc62[i], kN));
  const cdvec d = seq3gpp::lte_pss_freq(kRoots[nid2]);
  const auto& bank = seq3gpp::sss_bank()[nid2];  // [sub][nid1][62]
  struct Gap {
    const char* dup;
    const char* cp;
    int64_t gap;
  };
  const Gap gaps[4] = {{"FDD", "normal", kN + kCpN},
                       {"TDD", "normal", 3 * (kN + kCpN) + (kCp0N - kCpN)},
                       {"FDD", "extended", kN + kCpE},
                       {"TDD", "extended", 3 * (kN + kCpE)}};
  const int64_t Ly = static_cast<int64_t>(y.size());
  std::vector<SssResult> res;
  cdvec H(62), Hs(62), E(62);
  for (const Gap& g : gaps) {
    dvec acc(2 * 168, 0.0);
    int used = 0;
    for (const auto& [kh, p] : inst) {
      const int64_t s = p - g.gap;
      if (s < 0 || p + kN > Ly) continue;
      const cdvec Yp = fft::fft(y.data() + p, kN);
      const cdvec Ys = fft::fft(y.data() + s, kN);
      for (int i = 0; i < 62; ++i) H[i] = np::cdiv(Yp[sc[i]], d[i]);
      // maximum-ratio combining with the PSS channel estimate smoothed over kSssSmooth subcarriers per side
      smooth_half(H.data(), 31, Hs.data());
      smooth_half(H.data() + 31, 31, Hs.data() + 31);
      dvec e2(62);
      for (int i = 0; i < 62; ++i) {
        E[i] = np::cmul(Ys[sc[i]], std::conj(Hs[i]));
        const double a = np::cabs(E[i]);
        e2[i] = a * a;
      }
      const double rms = std::sqrt(np::mean(e2.data(), e2.size())) + 1e-30;
      for (int i = 0; i < 62; ++i) E[i] = np::cdiv(E[i], cd(rms, 0.0));
      double v[2][168];
      for (int sb = 0; sb < 2; ++sb)
        for (int n1 = 0; n1 < 168; ++n1) {
          const dvec& b = bank[sb][n1];
          double s2 = 0;
          for (int i = 0; i < 62; ++i) s2 += b[i] * E[i].real();
          v[sb][n1] = s2;
        }
      for (int par = 0; par < 2; ++par) {
        const int sb = pymod(kh + par, 2) == 0 ? 0 : 1;
        for (int n1 = 0; n1 < 168; ++n1) acc[par * 168 + n1] += v[sb][n1];
      }
      ++used;
    }
    if (!used) continue;
    dvec flat = acc;
    std::sort(flat.begin(), flat.end());
    const size_t am = np::argmax(acc.data(), acc.size());
    dvec scratch;
    const double med = median_copy(flat, scratch);
    const double noise = mad_of(flat, med, scratch);
    SssResult r;
    r.duplex = g.dup;
    r.cp = g.cp;
    r.nid1 = static_cast<int>(am % 168);
    r.parity = static_cast<int>(am / 168);
    r.n = used;
    r.z = (flat[flat.size() - 1] - med) / noise;
    r.z2 = (flat[flat.size() - 2] - med) / noise;
    res.push_back(r);
  }
  if (res.empty()) return std::nullopt;
  size_t bi = 0;
  auto key = [](const SssResult& r) { return r.z / std::max(r.z2, 1e-3); };
  for (size_t i = 1; i < res.size(); ++i)
    if (key(res[i]) > key(res[bi])) bi = i;
  return res[bi];
}

bool sss_decisive(const std::optional<SssResult>& s) {
  return s && s->z >= 4.5 && s->z / std::max(s->z2, 1e-3) >= 1.35;
}

std::optional<Json> wide_ppm_search(const Channelizer& ch, const dvec& centres, double span, double step) {
  std::optional<PssBest> best;
  double best_centre = 0;
  const dvec grid = np::arange<double>(-span, span + step / 2, step);
  std::shared_ptr<const TemplateBank> bank;
  for (double c : centres) {
    Slice sl(ch, c, bank);
    bank = nullptr;
    PssBest b = sl.best(grid, false, true);
    if (!best || b.z > best->z) {
      best = b;
      best_centre = c;
    }
  }
  if (!best || best->z < 30) return std::nullopt;
  // refine: the 5 kHz grid alone leaves +-2.5 kHz (+-2.7 ppm at 900 MHz)
  Slice sl(ch, best_centre);
  const PssBest fine = sl.best(np::arange<double>(best->df - 3e3, best->df + 3.01e3, 250.0), false, true);
  if (fine.z >= best->z) {
    best->z = fine.z;
    best->df = fine.df;
  }
  const double f = best_centre + best->df;
  const double raster = static_cast<double>(py_round_int(f / 100e3)) * 100e3;
  return Json{{"ppm", (f - raster) / raster * 1e6}, {"z", best->z}, {"carrier_hz", f}, {"raster_hz", raster}};
}

const double kFlatDb = 8.0;  // python FLAT_DB

std::optional<double> bandwidth_from_energy(const DetectResult& det, double fcar) {
  const Spectrogram& C = det.coarse;
  const fvec p20 = np::percentile_axis0(C.P, 20.0);
  dvec prof(p20.size());
  for (size_t j = 0; j < p20.size(); ++j) {
    float v = std::log10(p20[j] + 1e-30f);
    v = 10.0f * v;
    v = v + 6.5f;
    prof[j] = static_cast<double>(v) - det.floor_coarse_db[j];
  }
  std::optional<double> best;
  for (double bw : bws()) {
    const double h = occ(bw) / 2;
    size_t n = 0, above = 0;
    double first = 0, last = 0;
    dvec in;
    for (size_t j = 0; j < C.f.size(); ++j) {
      if (C.f[j] > fcar - h && C.f[j] < fcar + h) {
        if (n == 0) first = C.f[j];
        last = C.f[j];
        ++n;
        above += prof[j] > 3;
        in.push_back(prof[j]);
      }
    }
    if (n < 3 || first - C.f.front() < 0 || C.f.back() - last < 0) break;
    // python: flat profile (10-90th percentile spread <= FLAT_DB): adjacent carriers do not widen the claim
    const bool flat = np::percentile(in, 90.0) - np::percentile(in, 10.0) <= kFlatDb;
    if (static_cast<double>(above) / static_cast<double>(n) >= 0.85 && flat) {
      best = bw;
    } else {
      break;
    }
  }
  return best;
}

dvec raster_grid(const Signal& sig) {
  double lo = sig.f0 + 0.55e6, hi = sig.f1 - 0.55e6;
  if (hi < lo) lo = hi = static_cast<double>(py_round_int(sig.fc() / 100e3)) * 100e3;
  return np::arange<double>(std::ceil(lo / 100e3) * 100e3, hi + 1, 100e3);
}

dvec raster_offsets(double f_hz, bool ppm_known) {
  if (ppm_known && 2.5e-6 * f_hz < 3e3) return {0.0};
  return {-5e3, 0.0, 5e3};
}

std::vector<Cell> search_cells(const Channelizer& ch, const Signal& sig, double ppm, std::set<int64_t>* done,
                               bool ppm_known, bool require_sss) {
  dvec rasters;
  for (double r : raster_grid(sig)) {
    if (done) {
      const int64_t key = py_round_int(r / 100e3);
      if (done->count(key)) continue;
      done->insert(key);
    }
    rasters.push_back(r);
  }
  if (rasters.empty()) return {};
  struct Hit {
    double r;
    std::shared_ptr<Slice> sl;
    PssBest b;
  };
  // the slice length (hence nf and the template spectra) is the same for every raster
  auto first = std::make_shared<Slice>(ch, rasters[0] * (1 + ppm * 1e-6));
  const auto bank = template_bank(first->nf);
  auto probe = [&](size_t i) -> std::optional<Hit> {
    auto sl = i == 0 ? first : std::make_shared<Slice>(ch, rasters[i] * (1 + ppm * 1e-6), bank);
    PssBest b = sl->best(raster_offsets(rasters[i], ppm_known), true);
    if (!(b.z >= kZMin)) return std::nullopt;
    return Hit{rasters[i], sl, std::move(b)};
  };
  std::vector<size_t> idx(rasters.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  auto probed = pmap(idx, [&](size_t i) { return probe(i); });
  std::vector<Hit> hits;
  for (auto& h : probed)
    if (h) hits.push_back(std::move(*h));
  // non-max suppression: one cell per (carrier, N_ID2) within +-300 kHz
  std::stable_sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return -a.b.z < -b.b.z; });
  std::vector<Hit> keep;
  for (auto& h : hits) {
    bool dup = false;
    for (const auto& k : keep)
      if (std::fabs(h.r - k.r) < 300e3 && h.b.nid2 == k.b.nid2) dup = true;
    if (!dup) keep.push_back(std::move(h));
  }
  std::vector<size_t> kidx(keep.size());
  for (size_t i = 0; i < kidx.size(); ++i) kidx[i] = i;
  auto refined = pmap(kidx, [&](size_t i) -> std::optional<Cell> {
    const Hit& h = keep[i];
    const Slice& sl = *h.sl;
    PssBest fine = sl.best(np::arange<double>(h.b.df - 3e3, h.b.df + 3.01e3, 500.0), true);
    if (fine.nid2 != h.b.nid2) fine = h.b;
    const cdvec y = sl.shifted(fine.df);
    Instances inst = instances(fine.m, fine.phase);
    std::optional<SssResult> s = inst.size() >= 2 ? sss_decode(y, fine.nid2, inst) : std::nullopt;
    const double fcar = sl.f_center + fine.df;
    const bool ok = sss_decisive(s);
    // a loaded OFDM band leaves heavy-tailed fold peaks: without a decisive SSS
    // only a very strong, isolated PSS counts
    if (!ok && (require_sss || !(fine.z >= 40 && fine.iso >= 3))) return std::nullopt;
    Cell c;
    c.raster_hz = h.r;
    c.carrier_hz = fcar;
    c.pss_z = fine.z;
    c.pss_isolation = fine.iso;
    c.nid2 = fine.nid2;
    c.inst = std::move(inst);
    c.sss = s;
    c.sss_ok = ok;
    c.err_hz = fcar - h.r;
    c.phase = fine.phase;
    return c;
  });
  std::vector<Cell> cells;
  for (auto& c : refined)
    if (c) cells.push_back(std::move(*c));
  // Zadoff-Chu time-frequency ambiguity: drop a PSS-only hit when a stronger / SSS-resolved
  // cell with the same N_ID2 and timing (+-1 symbol) lies within 20 MHz.
  std::stable_sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) {
    if (a.sss_ok != b.sss_ok) return a.sss_ok;  // (not sss_ok) False first
    return -a.pss_z < -b.pss_z;
  });
  std::vector<Cell> out;
  for (auto& c : cells) {
    bool drop = false;
    if (!c.sss_ok)
      for (const auto& o : out)
        if (o.nid2 == c.nid2 && std::fabs(o.carrier_hz - c.carrier_hz) < 20e6 &&
            std::llabs(pymod(o.phase - c.phase + kHalf / 2, kHalf) - kHalf / 2) <= kN)
          drop = true;
    if (!drop) out.push_back(std::move(c));
  }
  return out;
}

}  // namespace si::rx::lte
