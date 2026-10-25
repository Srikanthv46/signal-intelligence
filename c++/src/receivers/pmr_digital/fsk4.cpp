// python: receivers/pmr_digital/fsk4.py (float32 discriminator as numpy keeps complex64's dtype; filters in
// float64 as np.convolve with a float64 kernel).
#include "receivers/pmr_digital/fsk4.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>

#include "core/fft.hpp"
#include "core/npcompat.hpp"

namespace si::rx::pmr_digital::fsk4 {

namespace {

constexpr int kSyncErrCoarse = 6;
constexpr int kSyncErr = 4;
constexpr double kSyncCorrMin = 0.9;
constexpr double kDmrDevUnit = 648.0;
constexpr double kGainLo = 0.5, kGainHi = 1.5;
constexpr uint64_t kNxdnFsw = 0xCDF59;

struct SyncPat {
  const char* name;
  uint64_t v;
};
constexpr SyncPat kDmrSync[9] = {{"bs_voice", 0x755FD7DF75F7},      {"bs_data", 0xDFF57D75DF5D},
                                 {"ms_voice", 0x7F7D5DD57DFD},      {"ms_data", 0xD5D7F77FD757},
                                 {"ms_rc", 0x77D55F7DFD77},         {"direct1_voice", 0x5D577F7757FF},
                                 {"direct1_data", 0xF7FDD5DDFD55}, {"direct2_voice", 0x7DFFD5F55D5F},
                                 {"direct2_data", 0xD7557F5FF7F5}};
const char* kDataTypes[12] = {"PI header",     "voice LC header", "terminator with LC", "CSBK",
                              "MBC header",    "MBC continuation", "data header",       "rate 1/2 data",
                              "rate 3/4 data", "idle",             "rate 1 data",       "unified single block data"};

std::vector<int8_t> bits_of(uint64_t v, int n) {
  std::vector<int8_t> b(n);
  for (int i = 0; i < n; ++i) b[i] = static_cast<int8_t>((v >> (n - 1 - i)) & 1);
  return b;
}

// python _symbols(v, n_bits): float32 symbols of a pattern
fvec symbols(uint64_t v, int n_bits) {
  const auto b = bits_of(v, n_bits);
  fvec s(n_bits / 2);
  for (int i = 0; i < n_bits / 2; ++i) {
    const int x = b[2 * i], y = b[2 * i + 1];
    s[i] = x == 0 ? (y == 1 ? 3.0f : 1.0f) : (y == 0 ? -1.0f : -3.0f);
  }
  return s;
}

inline void dibit(int s, std::vector<int8_t>& out) {
  switch (s) {
    case 3: out.push_back(0); out.push_back(1); break;
    case 1: out.push_back(0); out.push_back(0); break;
    case -1: out.push_back(1); out.push_back(0); break;
    default: out.push_back(1); out.push_back(1); break;
  }
}

// np.clip(np.round((v + 3) / 2) * 2 - 3, -3, 3) (np.round: half to even)
inline double slice4(double v) { return std::clamp(std::nearbyint((v + 3) / 2) * 2 - 3, -3.0, 3.0); }

dvec boxcar_same(const fvec& d, int k) {
  dvec dd(d.begin(), d.end());
  return np::convolve(dd, dvec(k, 1.0 / k), np::ConvMode::Same);
}

// python _rrc(beta, sps, span)
dvec rrc(double beta, double sps, int span = 6) {
  const int64_t h = static_cast<int64_t>(span * sps);
  dvec out;
  for (int64_t i = -h; i <= h; ++i) {
    const double t = static_cast<double>(i) / sps;
    const double x = (2 * beta) * t;
    out.push_back(np::sinc(t) * std::cos((kPi * beta) * t) / (1 - x * x + 1e-9));
  }
  const double s = np::sum(out);
  for (double& v : out) v = v / s;
  return out;
}

// python _coarse_bits(d, sps) -> (bits, phase)
std::pair<std::vector<int8_t>, double> coarse_bits(const fvec& d, double sps) {
  const int w = std::max<int>(1, static_cast<int>(py_round_int(0.6 * sps)));
  const dvec z = boxcar_same(d, w);
  bool have = false;
  double best_err = 0, best_ph = 0, best_a = 0;
  dvec best_v;
  for (double ph : np::arange<double>(0, sps, sps / 10)) {
    const int64_t cnt = static_cast<int64_t>((static_cast<double>(z.size()) - ph) / sps) - 1;
    if (cnt < 24) continue;
    dvec v(cnt), av(cnt);
    for (int64_t k = 0; k < cnt; ++k) {
      v[k] = z[static_cast<size_t>(ph + sps * static_cast<double>(k))];
      av[k] = std::fabs(v[k]);
    }
    const double a = np::percentile(av, 75.0) + 1e-9;
    dvec e2(cnt);
    for (int64_t k = 0; k < cnt; ++k) {
      const double x = v[k] / a * 3;
      const double q = slice4(x);
      e2[k] = (x - q) * (x - q);
    }
    const double err = np::mean(e2);
    if (!have || err < best_err) {
      have = true;
      best_err = err;
      best_v = std::move(v);
      best_a = a;
      best_ph = ph;
    }
  }
  std::vector<int8_t> bits;
  if (!have) return {bits, 0.0};
  for (double x : best_v) dibit(static_cast<int>(slice4(x / best_a * 3)), bits);
  return {bits, best_ph};
}

std::vector<std::pair<int64_t, int>> find(const std::vector<int8_t>& bits, uint64_t pattern, int n, int max_err) {
  std::vector<std::pair<int64_t, int>> out;
  if (bits.size() < static_cast<size_t>(n)) return out;
  const auto p = bits_of(pattern, n);
  for (size_t i = 0; i + n <= bits.size(); ++i) {
    int e = 0;
    for (int j = 0; j < n; ++j) e += bits[i + j] != p[j];
    if (e <= max_err) out.emplace_back(static_cast<int64_t>(i), e);
  }
  return out;
}

// ---------------------------------------------------------------------------------- DMR FEC (B.3.1, B.3.4, B.1.1)
constexpr int8_t kGolayG[8][20] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0},
                                   {0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1},
                                   {0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1},
                                   {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0},
                                   {0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1},
                                   {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1}};
constexpr int8_t kH15P[11][4] = {{1, 0, 0, 1}, {1, 1, 0, 1}, {1, 1, 1, 1}, {1, 1, 1, 0}, {0, 1, 1, 1}, {1, 0, 1, 0},
                                 {0, 1, 0, 1}, {1, 0, 1, 1}, {1, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 1, 1}};
constexpr int8_t kH13P[9][4] = {{1, 1, 1, 1}, {1, 1, 1, 0}, {0, 1, 1, 1}, {1, 0, 1, 0}, {0, 1, 0, 1},
                                {1, 0, 1, 1}, {1, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 1, 1}};

// python golay_20_8(r) -> (data byte or none, distance of the nearest codeword)
std::pair<std::optional<int>, int> golay_20_8(const std::vector<int8_t>& r) {
  int best = 0, bd = 99;
  for (int v = 0; v < 256; ++v) {
    int dist = 0;
    for (int j = 0; j < 20; ++j) {
      int c = 0;
      for (int i = 0; i < 8; ++i) c ^= ((v >> (7 - i)) & 1) & kGolayG[i][j];
      dist += c != r[j];
    }
    if (dist < bd) {
      bd = dist;
      best = v;
    }
  }
  return {bd <= 3 ? std::optional<int>(best) : std::nullopt, bd};
}

// python _hamming_fix(word, P): single-error correction of a systematic [data | data P] word, in place
template <int K>
bool hamming_fix(int8_t* const* w, const int8_t (*P)[4]) {
  int syn[4];
  for (int c = 0; c < 4; ++c) {
    int s = 0;
    for (int i = 0; i < K; ++i) s += *w[i] * P[i][c];
    syn[c] = (s + *w[K + c]) % 2;
  }
  if (!(syn[0] | syn[1] | syn[2] | syn[3])) return true;
  int hit = -1, nhit = 0;
  for (int pos = 0; pos < K + 4; ++pos) {
    bool eq = true;
    for (int c = 0; c < 4; ++c) {
      const int col = pos < K ? P[pos][c] : (pos - K == c);
      eq = eq && col == syn[c];
    }
    if (eq) {
      hit = pos;
      ++nhit;
    }
  }
  if (nhit == 1) {
    *w[hit] ^= 1;
    return true;
  }
  return false;
}

struct Bptc {
  std::vector<int8_t> info;
  int rows = 0, cols = 0;
};

Bptc bptc_196_96(const std::vector<int8_t>& r) {
  int8_t M[13][15];
  for (int i = 1; i < 196; ++i) M[(i - 1) / 15][(i - 1) % 15] = r[(181 * i) % 196];
  Bptc b;
  for (int j = 0; j < 15; ++j) {
    int8_t* w[13];
    for (int i = 0; i < 13; ++i) w[i] = &M[i][j];
    b.cols += hamming_fix<9>(w, kH13P);
  }
  for (int i = 0; i < 9; ++i) {
    int8_t* w[15];
    for (int j = 0; j < 15; ++j) w[j] = &M[i][j];
    b.rows += hamming_fix<11>(w, kH15P);
  }
  for (int j = 3; j < 11; ++j) b.info.push_back(M[0][j]);
  for (int i = 1; i < 9; ++i)
    for (int j = 0; j < 11; ++j) b.info.push_back(M[i][j]);
  return b;
}

int crc_ccitt(const int8_t* bits, int n) {
  int r = 0;
  for (int i = 0; i < n; ++i) {
    const int fb = ((r >> 15) & 1) ^ bits[i];
    r = (r << 1) & 0xFFFF;
    if (fb) r ^= 0x1021;
  }
  return r ^ 0xFFFF;
}

std::optional<int> crc_mask(int dt) {
  switch (dt) {
    case 0: return 0x6969;
    case 3: return 0xA5A5;
    case 4: return 0xAAAA;
    case 6: return 0xCCCC;
    case 11: return 0x3333;
    default: return std::nullopt;
  }
}

bool bptc_type(int dt) { return dt <= 7 || dt == 9 || dt == 11; }

void slot_type(const std::vector<int8_t>& st, Json& rec) {
  const auto [v, nerr] = golay_20_8(st);
  if (!v) return;
  rec["colour_code"] = *v >> 4;
  const int dt = *v & 15;
  rec["data_type"] = dt < 12 ? std::string(kDataTypes[dt]) : "reserved " + std::to_string(dt);
  rec["data_type_code"] = dt;
  rec["slot_type_errors"] = nerr;
}

void data_burst(const std::vector<int8_t>& bb, Json& rec) {
  std::vector<int8_t> st(bb.begin() + 98, bb.begin() + 108);
  st.insert(st.end(), bb.begin() + 156, bb.begin() + 166);
  slot_type(st, rec);
  if (!rec.contains("data_type_code")) return;
  const int dt = rec["data_type_code"].get<int>();
  if (!bptc_type(dt)) return;
  std::vector<int8_t> r(bb.begin(), bb.begin() + 98);
  r.insert(r.end(), bb.begin() + 166, bb.begin() + 264);
  const Bptc b = bptc_196_96(r);
  rec["bptc_rows_valid"] = b.rows;
  rec["bptc_cols_valid"] = b.cols;
  const auto mask = crc_mask(dt);
  if (mask && b.rows == 9 && b.cols == 15) {
    int rx = 0;
    for (int i = 80; i < 96; ++i) rx = (rx << 1) | b.info[i];
    const bool ok = crc_ccitt(b.info.data(), 80) == (rx ^ *mask);
    rec["crc_ok"] = ok;
    if (ok && dt == 3) {
      int o[10];
      for (int k = 0; k < 10; ++k) {
        o[k] = 0;
        for (int i = 0; i < 8; ++i) o[k] = (o[k] << 1) | b.info[8 * k + i];
      }
      static const char* hx = "0123456789abcdef";
      std::string data;
      for (int k = 2; k < 10; ++k) {
        data += hx[o[k] >> 4];
        data += hx[o[k] & 15];
      }
      rec["csbk_last_block"] = o[0] >> 7;
      rec["csbk_protect"] = (o[0] >> 6) & 1;
      rec["csbko"] = o[0] & 0x3F;
      rec["feature_set_id"] = o[1];
      rec["csbk_data_hex"] = data;
    }
  }
}

}  // namespace

fvec discriminator(const cfvec& y, double fs) {
  fvec d(y.size() > 0 ? y.size() - 1 : 0);
  const float sc = static_cast<float>(fs / (2 * kPi));
  // numpy temporary elision (>= 256 KiB): conj(prev) * cur
  const bool elide = d.size() * sizeof(cf) >= 256 * 1024;
  for (size_t i = 0; i < d.size(); ++i) {
    const cf z = elide ? np::cmul(std::conj(y[i]), y[i + 1]) : np::cmul(y[i + 1], std::conj(y[i]));
    d[i] = std::atan2(z.imag(), z.real()) * sc;
  }
  if (d.empty()) return d;
  const float med = np::median(d);
  for (float& v : d) v = v - med;
  return d;
}

std::pair<std::optional<double>, double> symbol_line(const fvec& d, double fs) {
  const int k = std::max<int>(1, static_cast<int>(py_round_int(fs / 9600)));
  const dvec dl = boxcar_same(d, k);
  if (dl.size() < 65) return {std::nullopt, 0.0};
  dvec e(dl.size() - 1);
  for (size_t i = 0; i + 1 < dl.size(); ++i) {
    const double x = dl[i + 1] - dl[i];
    e[i] = x * x;
  }
  const double mu = np::mean(e);
  for (double& v : e) v = v - mu;
  const cdvec E = fft::rfft(e);
  const double val = 1.0 / (static_cast<double>(e.size()) * (1 / fs));
  std::optional<double> rate;
  double best = 0.0;
  for (double r : {2400.0, 4800.0}) {
    double on = -1.0;
    dvec bg;
    for (size_t j = 0; j < E.size(); ++j) {
      const double f = static_cast<double>(j) * val;
      const double a = np::cabs(E[j]);
      if (std::fabs(f - r) <= 15) on = std::max(on, a);
      if (std::fabs(f - r) > 40 && std::fabs(f - r) < 400) bg.push_back(a);
    }
    if (on >= 0 && !bg.empty()) {
      const double q = on / (np::median(bg) + 1e-30);
      if (q > best) {
        best = q;
        rate = r;
      }
    }
  }
  return {rate, best};
}

std::vector<Json> dmr(const fvec& d, double fs, double sr) {
  const double sps = fs / sr;
  const auto [bits, ph0] = coarse_bits(d, sps);
  (void)ph0;
  dvec dd(d.begin(), d.end());
  const dvec z = np::convolve(dd, rrc(0.2, sps), np::ConvMode::Same);
  const int64_t nz = static_cast<int64_t>(z.size());
  std::vector<Json> out;
  std::vector<int64_t> seen;
  const dvec dts = np::arange<double>(-1.5 * sps, 1.5 * sps, sps / 20);
  for (const SyncPat& sp : kDmrSync) {
    const fvec ref = symbols(sp.v, 48);
    double nref2 = 0;
    {
      float s = 0;
      for (float r : ref) s += r * r;
      nref2 = static_cast<double>(std::sqrt(s));  // np.linalg.norm(float32) -> float32
    }
    dvec refd(ref.begin(), ref.end());
    const auto pat_bits = bits_of(sp.v, 48);
    const std::string name = sp.name;
    const bool is_data = name.size() >= 4 && name.compare(name.size() - 4, 4, "data") == 0;
    for (const auto& [h, herr] : find(bits, sp.v, 48, kSyncErrCoarse)) {
      (void)herr;
      const double s0 = static_cast<double>(h / 2) * sps;
      bool have = false;
      double bc = 0, bdt = 0;
      for (double dt : dts) {
        int64_t idx[24];
        bool ok = true;
        for (int k = 0; k < 24; ++k) {
          idx[k] = static_cast<int64_t>((s0 + dt) + sps * static_cast<double>(k));
          ok = ok && idx[k] >= 0 && idx[k] < nz;
        }
        if (!ok) continue;
        dvec v(24);
        for (int k = 0; k < 24; ++k) v[k] = z[idx[k]];
        const double mu = np::mean(v);
        double dot = 0, nv = 0;
        for (int k = 0; k < 24; ++k) {
          v[k] = v[k] - mu;
          dot += v[k] * refd[k];
          nv += v[k] * v[k];
        }
        const double c = dot / (std::sqrt(nv) * nref2 + 1e-9);
        if (!have || c > bc) {
          have = true;
          bc = c;
          bdt = dt;
        }
      }
      if (!have || bc < kSyncCorrMin) continue;
      const double t_sync = s0 + bdt;
      const int64_t key = py_round_int(t_sync / sps);
      bool dup = false;
      for (int64_t k : seen) dup = dup || std::llabs(key - k) < 12;
      if (dup) continue;
      dvec zi(24);
      for (int k = 0; k < 24; ++k) zi[k] = z[static_cast<int64_t>(t_sync + sps * static_cast<double>(k))];
      const dvec pf = np::polyfit(refd, zi, 1);
      const double gain = pf[0], dc = pf[1];
      if (!(kGainLo * kDmrDevUnit <= gain && gain <= kGainHi * kDmrDevUnit)) continue;
      Json rec = Json::object();
      rec["sync"] = name;
      rec["t_s"] = py_round(t_sync / fs, 6);
      rec["sync_corr"] = py_round(bc, 3);
      rec["deviation_hz_per_level"] = py_round(gain, 1);
      auto burst = [&](double start, int nsym, std::vector<int8_t>& bb, double* mse) -> bool {
        std::vector<int64_t> idx(nsym);
        for (int k = 0; k < nsym; ++k) {
          idx[k] = static_cast<int64_t>(start + sps * static_cast<double>(k));
          if (idx[k] < 0 || idx[k] >= nz) return false;
        }
        dvec e2(nsym);
        for (int k = 0; k < nsym; ++k) {
          const double sym = (z[idx[k]] - dc) / gain;
          const double q = slice4(sym);
          dibit(static_cast<int>(q), bb);
          e2[k] = (sym - q) * (sym - q);
        }
        if (mse) *mse = np::mean(e2);
        return true;
      };
      std::vector<int8_t> bb;
      double mse = 0;
      if (burst(t_sync - 54 * sps, 132, bb, &mse)) {
        int errs = 0;
        for (int i = 0; i < 48; ++i) errs += bb[108 + i] != pat_bits[i];
        rec["sync_bit_errors"] = errs;
        if (errs > kSyncErr) continue;
        rec["symbol_mse"] = py_round(mse, 3);
        if (is_data || name == "ms_rc") data_burst(bb, rec);
      } else if (is_data) {
        bb.clear();
        if (burst(t_sync - 5 * sps, 34, bb, nullptr)) {
          std::vector<int8_t> st(bb.begin(), bb.begin() + 10);
          st.insert(st.end(), bb.begin() + 58, bb.begin() + 68);
          slot_type(st, rec);
        }
      }
      seen.push_back(key);
      out.push_back(std::move(rec));
    }
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const Json& a, const Json& b) { return a["t_s"].get<double>() < b["t_s"].get<double>(); });
  return out;
}

std::vector<std::pair<int64_t, int>> nxdn(const fvec& d, double fs, double sr) {
  const double sps = fs / sr;
  const auto [bits, ph] = coarse_bits(d, sps);
  std::vector<std::pair<int64_t, int>> out;
  for (const auto& [h, e] : find(bits, kNxdnFsw, 20, 1))
    out.emplace_back(static_cast<int64_t>(ph + static_cast<double>(h / 2) * sps), e);
  return out;
}

int64_t frame_consistent(const std::vector<int64_t>& pos, double fs, double sr) {
  if (pos.empty()) return 0;
  const double fr = 192 * fs / sr, tol = 2 * fs / sr;
  int64_t best = 0;
  for (int64_t p : pos) {
    int64_t n = 0;
    for (int64_t q : pos) {
      const double x = static_cast<double>(q - p) / fr;
      n += std::fabs(x - static_cast<double>(py_round_int(x))) * fr <= tol;
    }
    best = std::max(best, n);
  }
  return best;
}

}  // namespace si::rx::pmr_digital::fsk4
