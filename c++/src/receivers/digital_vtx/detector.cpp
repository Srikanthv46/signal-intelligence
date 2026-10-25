// digital_vtx detector entry point and shared DSP (python: receivers/digital_vtx/detector.py
// _slice, _pss_time, _notch, detect).
#include "receivers/digital_vtx/detector.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/seq3gpp.hpp"
#include "core/spectral.hpp"
#include "receivers/digital_vtx/common.hpp"
#include "receivers/digital_vtx/hdzero.hpp"

namespace si::rx::digital_vtx {
namespace detail {

cdvec slice(const cdvec& X, double fs, double f_off, double half_bw, double out_fs, double* ofs_out) {
  const int64_t L = static_cast<int64_t>(X.size());
  const int64_t M = py_round_int(static_cast<double>(L) * out_fs / fs);
  const int64_t k0 = py_round_int(f_off / fs * static_cast<double>(L));
  const int64_t kk = trunc_i(half_bw / fs * static_cast<double>(L));
  cdvec Y(static_cast<size_t>(M), cd(0.0, 0.0));
  for (int64_t m = -kk; m <= kk; ++m) Y[static_cast<size_t>(pymod(m, M))] = X[static_cast<size_t>(pymod(k0 + m, L))];
  cdvec y = fft::ifft(Y);
  const cd s(static_cast<double>(M) / static_cast<double>(L), 0.0);
  for (auto& v : y) v = np::cmul(v, s);
  if (ofs_out) *ofs_out = static_cast<double>(M) * fs / static_cast<double>(L);
  return y;
}

cdvec pss_time(int u, double ofs) {
  const int64_t n = py_round_int(ofs / 15e3);
  const cdvec d = seq3gpp::lte_pss_freq(u);
  const auto& sc = seq3gpp::sc62();
  // theta = ((2pi*sc)*15e3)*t, as numpy's complex multiplies leave it (real part 0)
  double w[62];
  for (size_t k = 0; k < 62; ++k) w[k] = (2 * kPi * static_cast<double>(sc[k])) * 15e3;
  cdvec out(static_cast<size_t>(n));
  cdvec row(62);
  for (int64_t j = 0; j < n; ++j) {
    const double t = static_cast<double>(j) / ofs;
    for (size_t k = 0; k < 62; ++k) row[k] = np::cmul(d[k], expj(w[k] * t));
    out[static_cast<size_t>(j)] = np::pairwise_sum(row.data(), row.size());
  }
  return out;
}

// python hdzero_detect(x, fs, fc, X, ctx, info): one Emission per HDZero link over the whole capture
std::vector<Emission> hdzero_detect(const cfvec& x, double fs, double fc, const cdvec& X) {
  const hdzero::Slicer sl = [&](double f, double h, double r, double* o) { return slice(X, fs, f, h, r, o); };
  const auto links = hdzero::hdzero_links(sl, fs, fc, std::min(fs, 56e6) / 2, static_cast<int64_t>(x.size()));
  std::vector<Emission> out;
  if (links.empty()) return out;
  const double noise = noise_floor_psd(x).floor_power_per_sample;
  for (const auto& lk : links) {
    const double p_noise = noise * 2 * lk.half_bw / fs;
    const double snr = 10 * std::log10(std::max(lk.p_band / std::max(p_noise, 1e-30) - 1, 1e-3));
    const double lo = (fc + lk.lo) / 1e6, hi = (fc + lk.hi) / 1e6;
    const double cm = lk.channel_mhz + lk.cfo_hz / 1e6;
    Json f = Json::object();
    f["channel"] = lk.channel;
    f["channel_mhz"] = lk.channel_mhz;
    f["carrier_mhz"] = py_round(cm, 4);
    f["cfo_hz"] = py_round(lk.cfo_hz, 1);
    f["cfo_ppm"] = py_round(lk.ppm, 2);
    f["bw_mhz"] = py_round(hi - lo, 2);
    f["bandwidth_mode"] = "27 MHz (normal)";
    f["scs_khz"] = py_round(hdzero::SCS / 1e3, 3);
    f["tu_us"] = py_round(static_cast<double>(hdzero::N) / hdzero::FS * 1e6, 3);
    f["cp_us"] = py_round(static_cast<double>(hdzero::CP) / hdzero::FS * 1e6, 3);
    f["symbol_us"] = py_round(static_cast<double>(hdzero::TS) / hdzero::FS * 1e6, 3);
    f["fft_size"] = hdzero::N;
    f["sample_rate_msps"] = hdzero::FS / 1e6;
    f["used_subcarriers"] = 2 * hdzero::HALF + 1;
    f["continual_pilots"] = hdzero::N_CONTINUAL;
    f["pilot_coherence"] = py_round(lk.pilot_coherence, 3);
    f["data_coherence"] = py_round(lk.data_coherence, 3);
    f["cp_fold_contrast"] = py_round(lk.fold_contrast, 1);
    f["n_symbols"] = lk.n_symbols;
    f["preamble_ms"] = lk.preamble_ms;
    f["frame_symbols"] = lk.frame_symbols ? Json(*lk.frame_symbols) : Json(nullptr);
    f["frame_period_ms"] = lk.frame_ms ? Json(py_round(*lk.frame_ms, 3)) : Json(nullptr);
    f["frame_rate_fps"] = lk.frame_ms ? Json(py_round(1e3 / *lk.frame_ms, 2)) : Json(nullptr);
    f["link_snr_db"] = py_round(snr, 1);
    Emission e;
    e.family = "hdzero";
    e.label = "hdzero_video";
    e.start_us = py_round(lk.t0 * 1e6, 1);
    e.end_us = py_round(lk.t1 * 1e6, 1);
    e.f_lo_mhz = py_round(lo, 3);
    e.f_hi_mhz = py_round(hi, 3);
    e.center_mhz = py_round((lo + hi) / 2, 3);
    e.snr_db = py_round(snr, 1);
    e.confidence = 0.9;
    e.features = std::move(f);
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace detail

std::vector<Emission> detect(const cfvec& x, double fs, double fc, bool /*report_wifi: Tu < 10 us is never reported*/) {
  using namespace detail;
  const auto spurs = find_spurs(x, fs, fc);
  cdvec xd(x.size());
  for (size_t i = 0; i < x.size(); ++i) xd[i] = cd(x[i].real(), x[i].imag());
  cdvec X = fft::fft(xd);
  xd.clear();
  xd.shrink_to_fit();
  if (!spurs.empty()) {  // python _notch(X, fs, spurs)
    const dvec f = fft::fftfreq(X.size(), 1 / fs);
    for (const auto& s : spurs)
      for (size_t i = 0; i < X.size(); ++i)
        if (std::fabs(f[i] - s.offset_hz) <= 30e3) X[i] = cd(0.0, 0.0);
  }
  std::vector<Emission> out = herelink_detect(x.size(), fs, fc, X);
  auto hz = hdzero_detect(x, fs, fc, X);
  out.insert(out.end(), hz.begin(), hz.end());
  std::vector<Emission> excl;
  for (const auto& e : out)
    if (e.family == "herelink" || e.family == "hdzero") excl.push_back(e);
  auto gen = generic_detect(x, fs, fc, X, excl);
  out.insert(out.end(), gen.begin(), gen.end());
  std::stable_sort(out.begin(), out.end(), [](const Emission& a, const Emission& b) {
    return a.start_us != b.start_us ? a.start_us < b.start_us : a.center_mhz < b.center_mhz;
  });
  return out;
}

}  // namespace si::rx::digital_vtx
