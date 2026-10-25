// 315 / 433 MHz remotes and TPMS (research, signal gaps S6; python: receivers/ism_devices/detector.py): bursty
// carriers from a per-bin-normalised STFT, envelope pulses -> OOK decoders, FSK symbols -> TPMS.
#include "receivers/ism_devices/detector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "receivers/ism_devices/ook.hpp"
#include "receivers/ism_devices/tpms.hpp"

namespace si::rx::ism_devices {

const std::vector<std::pair<double, double>> kBands = {{314.0e6, 316.0e6}, {433.05e6, 434.79e6}};

namespace {

constexpr double kScanFs = 2e6;
constexpr size_t kNstft = 256;
constexpr double kOnDb = 12.0;
constexpr double kClusterHz = 100e3;
constexpr int64_t kMinFrames = 3;
constexpr size_t kMaxCand = 4;
constexpr double kSliceHalf = 100e3;
constexpr double kSliceFs = 1e6;
constexpr size_t kSmooth = 8;
constexpr double kMinRunUs = 30.0;
constexpr double kEnvDb = 12.0;
constexpr double kFskMinUs = 2000.0;
constexpr double kFskDevHz = 10e3;
constexpr double kFskFrac = 0.6;
constexpr double kSymbolUs = 52.0;
constexpr double kPi = 3.14159265358979323846;

// python _mavg(a, n): moving average over the "valid" range via a float64 cumsum
dvec mavg(const dvec& a, size_t n) {
  if (a.size() < n) return {};
  dvec c(a.size() + 1, 0.0);
  for (size_t i = 0; i < a.size(); ++i) c[i + 1] = c[i] + a[i];
  dvec out(a.size() + 1 - n);
  for (size_t i = 0; i < out.size(); ++i) out[i] = (c[i + n] - c[i]) / static_cast<double>(n);
  return out;
}

struct Run {
  size_t start, len;
  bool v;
};

std::vector<Run> runs(const std::vector<char>& on) {
  std::vector<Run> out;
  size_t s = 0;
  for (size_t i = 1; i <= on.size(); ++i)
    if (i == on.size() || on[i] != on[i - 1]) {
      out.push_back({s, i - s, on[s] != 0});
      s = i;
    }
  return out;
}

struct Pulses {
  std::vector<ook::Pulse> p;
  std::vector<std::pair<size_t, size_t>> marks;
  double env_db;
};

// python pulses(y, fs)
std::optional<Pulses> pulses(const cfvec& y, double fs) {
  dvec a(y.size());
  for (size_t i = 0; i < y.size(); ++i) a[i] = static_cast<double>(np::cabs(y[i]));
  const dvec env = mavg(a, kSmooth);
  if (env.empty()) return std::nullopt;
  const double noise = np::percentile(env, 10.0), peak = np::percentile(env, 99.5);
  if (noise <= 0 || 20 * std::log10(peak / noise) < kEnvDb) return std::nullopt;
  const double thr = std::sqrt(noise * peak);
  std::vector<char> on(env.size());
  for (size_t i = 0; i < env.size(); ++i) on[i] = env[i] > thr;
  const size_t shortn = static_cast<size_t>(py_round_int(kMinRunUs * 1e-6 * fs));
  for (const Run& r : runs(on))   // debounce: noise spikes / dropouts join their neighbours
    if (r.len < shortn) std::fill(on.begin() + static_cast<ptrdiff_t>(r.start),
                                  on.begin() + static_cast<ptrdiff_t>(r.start + r.len), static_cast<char>(!r.v));
  const double us = 1e6 / fs;
  const auto rr = runs(on);
  Pulses out;
  for (size_t i = 0; i < rr.size(); ++i) {
    if (!rr[i].v) continue;
    const double sp = i + 1 < rr.size() ? static_cast<double>(rr[i + 1].len) * us
                                        : std::numeric_limits<double>::infinity();
    out.p.push_back({static_cast<double>(rr[i].start) * us, static_cast<double>(rr[i].len) * us, sp});
    out.marks.emplace_back(rr[i].start + kSmooth / 2, rr[i].start + rr[i].len + kSmooth / 2);
  }
  out.env_db = py_round(20 * std::log10(peak / noise), 1);
  return out;
}

// python symbols(seg, fs)
std::optional<std::vector<int>> symbols(const cf* seg, size_t n, double fs) {
  if (n < 2 * kSmooth) return std::nullopt;
  // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
  const bool elide = (n - 1) * sizeof(cf) >= 256 * 1024;
  const double sc = fs / (2 * kPi);
  dvec d(n - 1);
  for (size_t i = 1; i < n; ++i) {
    const cf z = elide ? np::cmul(std::conj(seg[i - 1]), seg[i]) : np::cmul(seg[i], std::conj(seg[i - 1]));
    d[i - 1] = static_cast<double>(std::atan2(z.imag(), z.real())) * sc;
  }
  d = mavg(d, kSmooth);
  const double mid = 0.5 * (np::percentile(d, 10.0) + np::percentile(d, 90.0));   // mid-way between the tones
  dvec a(d.size());
  for (size_t i = 0; i < d.size(); ++i) {
    d[i] -= mid;
    a[i] = std::fabs(d[i]);
  }
  const double dev = np::median(a);   // two-level: |d| sits near one deviation (a chirp spreads it)
  size_t in = 0;
  for (double v : a) in += v >= 0.6 * dev && v <= 1.4 * dev;
  if (dev < kFskDevHz || static_cast<double>(in) / static_cast<double>(a.size()) < kFskFrac) return std::nullopt;
  const double spb = kSymbolUs * 1e-6 * fs;
  std::vector<char> pos(d.size());
  for (size_t i = 0; i < d.size(); ++i) pos[i] = d[i] > 0;
  std::vector<int> out;
  for (const Run& r : runs(pos)) {
    const int64_t k = py_round_int(static_cast<double>(r.len) / spb);
    for (int64_t q = 0; q < k; ++q) out.push_back(r.v ? 1 : 0);
  }
  return out;
}

Emission make(const std::string& label, double t0, double t1, double f, double env_db, double conf, Json decoded,
              Json features) {
  Emission e;
  e.family = "ism_devices";
  e.label = label;
  e.start_us = t0;
  e.end_us = t1;
  e.f_lo_mhz = (f - kSliceHalf) / 1e6;
  e.f_hi_mhz = (f + kSliceHalf) / 1e6;
  e.center_mhz = f / 1e6;
  e.snr_db = env_db;
  e.confidence = conf;
  e.decoded = std::move(decoded);
  e.features = std::move(features);
  return e;
}

}  // namespace

std::vector<Candidate> candidates(const Channelizer& ch) {
  std::vector<Candidate> out;
  for (const auto& [lo, hi] : kBands) {
    const double a = std::max(lo, ch.fc() - 0.45 * ch.fs()), b = std::min(hi, ch.fc() + 0.45 * ch.fs());
    if (b - a < 50e3) continue;
    const double mid = (a + b) / 2, half = (b - a) / 2;
    const auto sl = ch.slice(mid, half, std::max(kScanFs, b - a));
    const size_t nf = sl.y.size() / kNstft;
    if (nf < 8) continue;
    const dvec win = np::hanning(kNstft);
    const dvec fr = fft::fftfreq(kNstft, 1 / sl.fs);
    std::vector<dvec> P(nf, dvec(kNstft));
    cdvec in(kNstft), F(kNstft);
    for (size_t r = 0; r < nf; ++r) {
      for (size_t i = 0; i < kNstft; ++i) {
        const cf v = sl.y[r * kNstft + i];
        in[i] = cd(static_cast<double>(v.real()) * win[i], static_cast<double>(v.imag()) * win[i]);
      }
      fft::c2c(in.data(), F.data(), kNstft, true);
      for (size_t i = 0; i < kNstft; ++i) {
        const double m = np::cabs(F[i]);
        P[r][i] = m * m;
      }
    }
    dvec med(kNstft), col(nf);
    for (size_t i = 0; i < kNstft; ++i) {
      for (size_t r = 0; r < nf; ++r) col[r] = P[r][i];
      med[i] = np::median(col);
    }
    dvec f_on, d_on;
    for (size_t r = 0; r < nf; ++r) {
      size_t k = 0;
      double best = -1;
      for (size_t i = 0; i < kNstft; ++i) {
        const double v = std::fabs(fr[i]) > half ? 0.0 : P[r][i] / (med[i] + 1e-30);
        if (v > best) {
          best = v;
          k = i;
        }
      }
      if (10 * std::log10(best + 1e-30) >= kOnDb) {
        f_on.push_back(fr[k] + mid);
        d_on.push_back(10 * std::log10(best));
      }
    }
    if (f_on.empty()) continue;
    const auto order = np::argsort_stable(f_on);
    dvec fs_(f_on.size()), ds_(f_on.size());
    for (size_t i = 0; i < order.size(); ++i) {
      fs_[i] = f_on[static_cast<size_t>(order[i])];
      ds_[i] = d_on[static_cast<size_t>(order[i])];
    }
    std::vector<size_t> cuts = {0};
    for (size_t i = 0; i + 1 < fs_.size(); ++i)
      if (fs_[i + 1] - fs_[i] > kClusterHz) cuts.push_back(i + 1);
    cuts.push_back(fs_.size());
    for (size_t c = 0; c + 1 < cuts.size(); ++c) {
      const size_t c0 = cuts[c], c1 = cuts[c + 1];
      if (static_cast<int64_t>(c1 - c0) < kMinFrames) continue;
      const double f = np::median(dvec(fs_.begin() + static_cast<ptrdiff_t>(c0), fs_.begin() + static_cast<ptrdiff_t>(c1)));
      const double mx = *std::max_element(ds_.begin() + static_cast<ptrdiff_t>(c0), ds_.begin() + static_cast<ptrdiff_t>(c1));
      out.push_back({f, static_cast<int64_t>(c1 - c0), py_round(mx, 1)});
    }
  }
  std::stable_sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) { return a.frames > b.frames; });
  if (out.size() > kMaxCand) out.resize(kMaxCand);
  return out;
}

std::vector<Emission> detect(const Channelizer& ch) {
  std::vector<Emission> out;
  for (const Candidate& c : candidates(ch)) {
    const double f = c.f;
    const auto sl = ch.slice(f, kSliceHalf, kSliceFs);
    const auto pr = pulses(sl.y, sl.fs);
    if (!pr) continue;
    const auto& p = pr->p;
    Json base = Json::object();
    base["carrier_mhz"] = py_round(f / 1e6, 4);
    base["frames"] = c.frames;
    base["scan_db"] = c.peak_db;
    base["envelope_db"] = pr->env_db;
    base["n_pulses"] = static_cast<int64_t>(p.size());
    bool found = false;
    int64_t fsk_marks = 0;
    for (size_t i = 0; i < p.size(); ++i) {
      if (p[i].mark_us < kFskMinUs) continue;
      const auto [s0, s1] = pr->marks[i];
      const size_t e1 = std::min(s1, sl.y.size());
      if (s0 >= e1) continue;
      const auto sy = symbols(sl.y.data() + s0, e1 - s0, sl.fs);
      if (!sy) continue;
      ++fsk_marks;
      if (auto r = tpms::decode(*sy)) {
        found = true;
        Json feats = base;
        feats["symbols"] = static_cast<int64_t>(sy->size());
        feats["mark_us"] = py_round(p[i].mark_us, 1);
        out.push_back(make("tpms", p[i].t_us, p[i].t_us + p[i].mark_us, f, pr->env_db, 0.9, *r, feats));
      }
    }
    if (auto r = ook::fixed_code(p)) {
      found = true;
      out.push_back(make("ook_fixed_code", (*r)["t_first_us"].get<double>(), (*r)["t_end_us"].get<double>(), f,
                         pr->env_db, (*r)["repeats"].get<int64_t>() >= 2 ? 0.75 : 0.45, *r, base));
    }
    if (auto r = ook::keeloq(p)) {
      found = true;
      out.push_back(make("keeloq", (*r)["t_first_us"].get<double>(), p.back().t_us + p.back().mark_us, f, pr->env_db,
                         (*r)["matching_packets"].get<int64_t>() >= 2 ? 0.75 : 0.6, *r, base));
    }
    if (found) continue;
    dvec mk;
    for (const auto& q : p) mk.push_back(q.mark_us);
    const bool ook_like = p.size() >= 16 && std::all_of(mk.begin(), mk.end(), [](double m) { return 50.0 <= m && m <= 2000.0; });
    if (fsk_marks) {
      Json feats = base;
      feats["fsk_bursts"] = fsk_marks;
      out.push_back(make("fsk_device", p.front().t_us, p.back().t_us + p.back().mark_us, f, pr->env_db, 0.55, nullptr,
                         feats));
    } else if (ook_like) {
      Json feats = base;
      feats["mark_short_us"] = py_round(np::percentile(mk, 10.0), 1);
      feats["mark_long_us"] = py_round(np::percentile(mk, 90.0), 1);
      out.push_back(make("ook_device", p.front().t_us, p.back().t_us + p.back().mark_us, f, pr->env_db, 0.55, nullptr,
                         feats));
    }
  }
  return out;
}

}  // namespace si::rx::ism_devices
