// python receivers/cofdm/detector.py (same stages, order, keys).
#include "receivers/cofdm/detector.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"
#include "receivers/cofdm/dvbt2_p1.hpp"
#include "receivers/cofdm/tables.hpp"

namespace si::rx::cofdm {

namespace {

constexpr double kPi = 3.14159265358979323846;
const double kBandwidths[] = {8.0, 7.0, 6.0, 5.0, 4.0, 2.5, 2.0};
constexpr double kOcc = 0.951, kPeakMin = 0.2, kContrast = 3.0, kPilotDb = 1.5, kPilotMargin = 1.2;
constexpr int64_t kMinSyms = 8;

struct Mode { const char* name; int64_t N, K; const std::vector<int>* cp; const std::vector<int>* tps; };
const std::vector<Mode>& modes() {
  static const std::vector<Mode> m = {{"2K", 2048, 1705, &tables::cp_2k(), &tables::tps_2k()},
                                      {"8K", 8192, 6817, &tables::cp_8k(), &tables::tps_8k()}};
  return m;
}
const std::vector<std::pair<const char*, int64_t>> kGuard = {{"1/4", 4}, {"1/8", 8}, {"1/16", 16}, {"1/32", 32}};

struct Slice { cdvec y; double ofs; };
// python _slice (brick-wall slice of the capture FFT)
Slice slice(const cdvec& X, double fs, double f_off, double half_bw, double out_fs) {
  const int64_t L = static_cast<int64_t>(X.size());
  const int64_t M = py_round_int(static_cast<double>(L) * out_fs / fs);
  const int64_t k0 = py_round_int(f_off / fs * static_cast<double>(L));
  const int64_t kk = static_cast<int64_t>(half_bw / fs * static_cast<double>(L));
  cdvec Y(static_cast<size_t>(M), cd(0, 0));
  for (int64_t m = -kk; m <= kk; ++m) Y[static_cast<size_t>(((m % M) + M) % M)] = X[static_cast<size_t>((((k0 + m) % L) + L) % L)];
  cdvec y = fft::ifft(Y);
  const double sc = static_cast<double>(M) / static_cast<double>(L);
  for (auto& v : y) v *= sc;
  return {std::move(y), static_cast<double>(M) * fs / static_cast<double>(L)};
}

struct Cand { double off, B, width, lev; };
std::vector<Cand> candidates(const cfvec& x, double fs) {
  const size_t nper = std::min<size_t>(1 << 14, x.size());
  auto W = sig::welch(x, fs, nper, -1, false);
  const size_t n = W.f.size();
  dvec f(n);
  std::vector<float> P(n);
  for (size_t i = 0; i < n; ++i) {  // fftshift
    const size_t j = (i + (n + 1) / 2) % n;
    f[i] = W.f[j];
    P[i] = W.Pxx[j];
  }
  const float fl = np::percentile(P, 20.0);
  const float thr = static_cast<float>(1.6) * fl;
  std::vector<int> on(n, 0);
  for (size_t i = 0; i < n; ++i) {  // convolve(on, ones(5), "same") > 0
    bool any = false;
    for (int d = -2; d <= 2; ++d) {
      const int64_t k = static_cast<int64_t>(i) + d;
      if (k >= 0 && k < static_cast<int64_t>(n) && P[static_cast<size_t>(k)] > thr) any = true;
    }
    on[i] = any;
  }
  const double df = f[1] - f[0];
  std::vector<Cand> out;
  size_t i = 0;
  while (i < n) {
    if (!on[i]) { ++i; continue; }
    size_t b = i;
    while (b < n && on[b]) ++b;
    const double w = static_cast<double>(b - i) * df;
    for (double B : kBandwidths)
      if (std::fabs(w - kOcc * B * 1e6) <= 0.1 * kOcc * B * 1e6) {
        std::vector<float> seg(P.begin() + static_cast<ptrdiff_t>(i), P.begin() + static_cast<ptrdiff_t>(b));
        const float lev = 10.0f * std::log10(np::median(seg) / fl);
        out.push_back({(f[i] + f[b - 1]) / 2, B, w, static_cast<double>(lev)});
        break;
      }
    i = b;
  }
  return out;
}

struct Fold { double peak, con; int64_t ph; cd csum; int64_t nsym; };
std::optional<Fold> cp_fold(const cdvec& y, int64_t N, int64_t Ng) {
  const int64_t L = N + Ng;
  const int64_t ny = static_cast<int64_t>(y.size());
  int64_t nsym = (ny - N) / L;
  if (nsym < kMinSyms) return std::nullopt;
  const int64_t nc = ny - N;
  cdvec cs(static_cast<size_t>(nc + 1), cd(0, 0));
  dvec es(static_cast<size_t>(nc + 1), 0.0);
  for (int64_t i = 0; i < nc; ++i) {
    cs[i + 1] = cs[i] + y[i] * std::conj(y[i + N]);
    es[i + 1] = es[i] + (np::cabs2(y[i]) + np::cabs2(y[i + N])) / 2;
  }
  const int64_t nm = nc + 1 - Ng;
  if (nm < nsym * L) --nsym;
  cdvec F(static_cast<size_t>(L), cd(0, 0));
  dvec E(static_cast<size_t>(L), 0.0);
  for (int64_t s = 0; s < nsym; ++s)
    for (int64_t j = 0; j < L; ++j) {
      const int64_t q = s * L + j;
      F[j] += cs[q + Ng] - cs[q];
      E[j] += es[q + Ng] - es[q];
    }
  dvec prof(static_cast<size_t>(L));
  for (int64_t j = 0; j < L; ++j) prof[j] = std::abs(F[j]) / (E[j] + 1e-30);
  const int64_t p = static_cast<int64_t>(np::argmax(prof));
  return Fold{prof[p], prof[p] / (np::median(prof) + 1e-12), p, F[p], nsym};
}

}  // namespace

Json tps_decode(const std::vector<int>& bits) {
  const auto& sync = tables::sync();
  Json res = nullptr;
  const int64_t nb = static_cast<int64_t>(bits.size());
  for (int64_t p = 1; p < nb - 66; ++p) {
    bool eq = true, inv = true;
    for (int i = 0; i < 16; ++i) {
      eq = eq && bits[p + i] == sync[i];
      inv = inv && bits[p + i] == 1 - sync[i];
    }
    if (!eq && !inv) continue;
    std::vector<int> s(68);
    s[0] = 0;
    for (int i = 0; i < 67; ++i) s[i + 1] = bits[p + i];
    const std::vector<int> info(s.begin() + 1, s.begin() + 54), par(s.begin() + 54, s.end());
    const bool ok = tables::bch_parity(info) == par;
    auto v = [&](int a, int n) {
      int r = 0;
      for (int i = 0; i < n; ++i) r = (r << 1) | s[a + i];
      return r;
    };
    auto bin = [](int val, int n) {
      std::string o;
      for (int i = n - 1; i >= 0; --i) o += ((val >> i) & 1) ? '1' : '0';
      return o;
    };
    const int len = v(17, 6), hier = v(27, 3), hp = v(30, 3), lp = v(33, 3);
    static const char* kCon[] = {"QPSK", "16-QAM", "64-QAM", "reserved"};
    static const char* kHier[] = {"non-hierarchical", "alpha 1", "alpha 2", "alpha 4"};
    static const char* kCr[] = {"1/2", "2/3", "3/4", "5/6", "7/8"};
    static const char* kGi[] = {"1/32", "1/16", "1/8", "1/4"};
    static const char* kMode[] = {"2K", "8K", "4K", "reserved"};
    std::string length = len == 0b010111 ? "23 bits (no cell ID)" : len == 0b011111 ? "31 bits (cell ID)"
                         : len == 0b100001 ? "33 bits (DVB-H)" : bin(len, 6);
    Json r = {{"bch_ok", ok},
              {"frame_number", v(23, 2) + 1},
              {"length_indicator", length},
              {"constellation", kCon[v(25, 2)]},
              {"hierarchy", hier < 4 ? std::string(kHier[hier]) : bin(hier, 3)},
              {"code_rate_hp", hp < 5 ? std::string(kCr[hp]) : std::string("reserved")},
              {"code_rate_lp", lp < 5 ? std::string(kCr[lp]) : std::string("reserved")},
              {"guard_interval", kGi[v(36, 2)]},
              {"mode", kMode[v(38, 2)]},
              {"frame_start_symbol", p - 1}};
    if (len == 0b011111) {
      const int fn = r["frame_number"].get<int>();
      r["cell_id_byte"] = fmt("0x%02X", v(40, 8)) + ((fn == 1 || fn == 3) ? " (high)" : " (low)");
    }
    if (ok) return r;
    if (res.is_null()) res = r;
  }
  return res;
}

std::vector<Emission> detect(const cfvec& x, double fs, double fc) {
  std::vector<Emission> out;
  const cdvec X = fft::fft(cdvec(x.begin(), x.end()));
  for (const auto& c : candidates(x, fs)) {
    const double clk = 8e6 * c.B / 7;
    const Slice sl = slice(X, fs, c.off, 0.5 * c.B * 1e6, clk);
    const cdvec& y = sl.y;
    const double ofs = sl.ofs;
    {  // DVB-T2 P1 preamble (python analyse_t2; research, PLAN D5)
      // python P1_MARGIN: a real P1 decodes S1 and S2 with margins >= 0.74 even at -6 dB; correlation peaks from ADS-B
      // pulses, analog video and clipped captures decoded with <= 0.13 on S1 or S2 (regression 2026-10-22)
      auto found = p1::find(y, 0.35);
      found.erase(std::remove_if(found.begin(), found.end(), [](const auto& f) {
                    return !(f.decoded["s1_margin"].template get<double>() >= 0.4 &&
                             f.decoded["s2_margin"].template get<double>() >= 0.4 && f.decoded["preamble"] != "reserved");
                  }), found.end());
      if (!found.empty()) {
        const auto& f0 = found[0];
        const double scs = ofs / p1::kN1;
        Json F = Json::object();
        F["standard"] = "DVB-T2 (EN 302 755) P1 preamble";
        F["channel_bw_mhz"] = c.B;
        F["occupied_mhz"] = np::round(c.width / 1e6, 3);
        F["elementary_rate_mhz"] = py_round(ofs / 1e6, 6);
        F["p1_count"] = found.size();
        F["p1_metric"] = py_round(f0.metric, 3);
        F["p1_start_us"] = py_round(static_cast<double>(f0.start) / ofs * 1e6, 1);
        F["p1"] = f0.decoded;
        F["carrier_mhz"] = np::round((fc + c.off + f0.decoded["carrier_offset"].get<double>() * scs) / 1e6, 6);
        F["level_db"] = py_round(c.lev, 1);
        const bool ok = f0.decoded["s1_margin"].get<double>() > 0.2 && f0.decoded["s2_margin"].get<double>() > 0.2;
        const double lo = (fc + c.off - c.width / 2) / 1e6, hi = (fc + c.off + c.width / 2) / 1e6;
        Emission e;
        e.family = "cofdm";
        e.label = ok ? "dvbt2" : "dvbt2_like";
        e.start_us = 0.0;
        e.end_us = py_round(static_cast<double>(x.size()) / fs * 1e6, 1);
        e.f_lo_mhz = np::round(lo, 3);
        e.f_hi_mhz = np::round(hi, 3);
        e.center_mhz = np::round((lo + hi) / 2, 3);
        e.snr_db = F["level_db"].get<double>();
        e.confidence = ok ? 0.95 : 0.7;
        e.features = std::move(F);
        out.push_back(std::move(e));
        continue;
      }
    }
    // numerology: 2K / 8K x GI
    std::optional<Fold> best;
    const Mode* bm = nullptr;
    const char* bgi = nullptr;
    int64_t bNg = 0;
    for (const auto& m : modes())
      for (const auto& [gi, d] : kGuard) {
        auto r = cp_fold(y, m.N, m.N / d);
        if (!r) continue;
        if (r->peak >= kPeakMin && r->con >= kContrast && (!best || r->peak > best->peak)) {
          best = r;
          bm = &m;
          bgi = gi;
          bNg = m.N / d;
        }
      }
    if (!best) continue;
    const int64_t N = bm->N, K = bm->K, Ng = bNg, L = N + Ng;
    const double eps = -std::arg(best->csum) / (2 * kPi * static_cast<double>(N));
    const int64_t ny = static_cast<int64_t>(y.size());
    cdvec yc(y.size());
    for (int64_t i = 0; i < ny; ++i) yc[i] = y[i] * std::exp(cd(0, (-2.0 * kPi * eps) * static_cast<double>(i)));
    std::vector<int64_t> starts;
    for (int64_t s = 0; s < best->nsym; ++s) {
      const int64_t st = best->ph + s * L;
      if (st + Ng + N <= ny) starts.push_back(st);
    }
    const int64_t offw = Ng / 2;
    std::vector<cdvec> Z;
    cdvec w(static_cast<size_t>(N));
    for (int64_t st : starts) {
      std::copy(yc.begin() + st + offw, yc.begin() + st + offw + N, w.begin());
      Z.push_back(fft::fft(w));
    }
    const int64_t kc = (K - 1) / 2;
    dvec Pm(static_cast<size_t>(N), 0.0);
    for (const auto& z : Z)
      for (int64_t k = 0; k < N; ++k) Pm[k] += np::cabs2(z[k]);
    for (auto& v : Pm) v /= static_cast<double>(Z.size());
    const double pmean = np::mean(Pm);
    const int64_t S = N / 8;
    dvec sc;
    for (int64_t s = -S; s <= S; ++s) {
      dvec vals;
      for (int k : *bm->cp) vals.push_back(Pm[static_cast<size_t>((((k - kc + s) % N) + N) % N)]);
      sc.push_back(np::mean(vals) / (pmean + 1e-30));
    }
    const size_t j = np::argmax(sc);
    const int64_t si = static_cast<int64_t>(j) - S;
    const double med = np::median(sc);
    const double pilot_db = 10 * std::log10(sc[j] / (med + 1e-30));
    dvec dev(sc.size());
    for (size_t i = 0; i < sc.size(); ++i) dev[i] = std::fabs(sc[i] - med);
    const double mad = 1.4826 * np::median(dev) + 1e-12;
    const double pilot_z = (sc[j] - med) / mad;
    double other = -1e300;
    for (size_t i = 0; i < sc.size(); ++i)
      if (i != j) other = std::max(other, sc[i]);
    const double margin = sc[j] / (other + 1e-30);
    const double scs = ofs / static_cast<double>(N);
    Json F = Json::object();
    F["standard"] = "DVB-T (EN 300 744) framing";
    F["channel_bw_mhz"] = c.B;
    F["occupied_mhz"] = np::round(c.width / 1e6, 3);
    F["mode"] = bm->name;
    F["guard_interval"] = bgi;
    F["fft_size"] = N;
    F["elementary_rate_mhz"] = py_round(ofs / 1e6, 6);
    F["scs_hz"] = py_round(scs, 2);
    F["tu_us"] = py_round(static_cast<double>(N) / ofs * 1e6, 3);
    F["symbol_us"] = py_round(static_cast<double>(L) / ofs * 1e6, 3);
    F["cp_peak"] = py_round(best->peak, 3);
    F["cp_contrast"] = py_round(best->con, 1);
    F["symbols"] = starts.size();
    F["carrier_mhz"] = np::round((fc + c.off + static_cast<double>(si) * scs + eps * ofs) / 1e6, 6);
    F["integer_offset_carriers"] = si;
    F["continual_pilot_db"] = py_round(pilot_db, 1);
    F["continual_pilot_z"] = py_round(pilot_z, 1);
    F["continual_pilot_margin"] = py_round(margin, 2);
    F["level_db"] = py_round(c.lev, 1);
    if (pilot_db < kPilotDb || margin < kPilotMargin) {
      continue;  // DVB-T framing always carries its continual pilots (python analyse returns None)
    }
    {
      std::vector<int> bits = {0};
      for (size_t l = 1; l < Z.size(); ++l) {
        double d = 0;
        for (int k : *bm->tps) {
          const size_t b = static_cast<size_t>((((k - kc + si) % N) + N) % N);
          d += std::real(Z[l][b] * std::conj(Z[l - 1][b]));
        }
        bits.push_back(d < 0 ? 1 : 0);
      }
      const Json r = tps_decode(bits);
      F["tps"] = r;
      if (r.is_null())
        F["tps_note"] = fmt("no TPS sync word in %d symbols (a TPS frame is 68 symbols, %.1f ms)", static_cast<int>(bits.size()),
                            68.0 * static_cast<double>(L) / ofs * 1e3);
    }
    const double lo = (fc + c.off - c.width / 2) / 1e6, hi = (fc + c.off + c.width / 2) / 1e6;
    const bool bch = F["tps"].is_object() && F["tps"]["bch_ok"].get<bool>();
    Emission e;
    e.family = "cofdm";
    e.label = bch ? "dvbt" : "dvbt_like";
    e.start_us = 0.0;
    e.end_us = py_round(static_cast<double>(x.size()) / fs * 1e6, 1);
    e.f_lo_mhz = np::round(lo, 3);
    e.f_hi_mhz = np::round(hi, 3);
    e.center_mhz = np::round((lo + hi) / 2, 3);
    e.snr_db = F["level_db"].get<double>();
    e.confidence = bch ? 0.95 : F["tps"].is_object() ? 0.75 : 0.6;
    e.features = std::move(F);
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace si::rx::cofdm
