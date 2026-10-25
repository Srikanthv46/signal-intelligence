// Shared DSP of the research receivers (python: core/spectral.py): FFT channelizer with
// spur notches (complex128), CW spur finder, CP metrics, PSD noise floor, spectrogram,
// robust threshold, in-band SNR. add_awgn (evaluation only) is not ported.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si {

// python core.spectral.Channelizer: np.fft.fft of the whole capture (complex128), spur
// bins |f - f0| <= half zeroed; slice(f_off, half_bw, out_fs) is the complex baseband of
// [f_off +- half_bw] (offset from the capture centre) at M = round(L*out_fs/fs) samples.
// Unlike si::Channelizer (core/dsp.hpp) M has no minimum and offsets are relative.
class SpectralChannelizer {
 public:
  SpectralChannelizer(const cdvec& x, double fs, const std::vector<std::pair<double, double>>& notch_hz = {});
  cdvec slice(double f_off, double half_bw, double out_fs) const;

  const cdvec& X() const { return X_; }
  const dvec& f() const { return f_; }
  size_t L() const { return L_; }
  double fs() const { return fs_; }

 private:
  double fs_;
  size_t L_;
  cdvec X_;
  dvec f_;
};

// Indices i (ascending, as np.where) with |f[i] - f_off| <= half_bw for f = fftfreq(L, 1/fs),
// without scanning all L bins; the boundary bins are tested with the same expression.
std::vector<size_t> select_bins(const dvec& f, double f_off, double half_bw);

struct Spur {
  double offset_hz = 0;
  double freq_mhz = 0;       // round(.., 4)
  double prominence_db = 0;  // round(.., 1)
  std::string note;          // "LO leakage (DC)" | "CW tone"
};
Json to_json(const Spur& s);

// python find_spurs(x, fs, fc, nfft=8192, min_prom_db=15.0)
std::vector<Spur> find_spurs(const cfvec& x, double fs, double fc, size_t nfft = 8192, double min_prom_db = 15.0);
std::vector<Spur> find_spurs(const cdvec& x, double fs, double fc, size_t nfft = 8192, double min_prom_db = 15.0);

// python _cp_metric(y, N, cp) -> (m, c): normalised CP correlation over a cp window
struct CpMetric {
  dvec m;
  cdvec c;
};
CpMetric cp_metric(const cdvec& y, size_t N, size_t cp);
CpMetric cp_metric(const cd* y, size_t len, size_t N, size_t cp);

// python _folded_scores(m, sym, W, hop)
dvec folded_scores(const dvec& m, size_t sym, size_t W, size_t hop);

// python noise_floor_psd(x, nfft=1024)
struct NoiseFloor {
  dvec floor_bin;
  double floor_power_per_sample = 0;
  size_t nfft = 0;
};
NoiseFloor noise_floor_psd(const cfvec& x, size_t nfft = 1024);
NoiseFloor noise_floor_psd(const cdvec& x, size_t nfft = 1024);

// python spectrogram(x, fs, fc, nfft=1024, hop=None): P (frames x nfft, fftshifted,
// linear power), frame-centre times (us), bin frequencies (MHz)
struct SpectralSpectrogram {
  dmat P;
  dvec t_us;
  dvec f_mhz;
};
SpectralSpectrogram spectrogram(const cfvec& x, double fs, double fc, size_t nfft = 1024, size_t hop = 0);
SpectralSpectrogram spectrogram(const cdvec& x, double fs, double fc, size_t nfft = 1024, size_t hop = 0);

// python robust_threshold(v, k_sigma) -> (median + k*1.4826*MAD, median)
std::pair<double, double> robust_threshold(const dvec& v, double k_sigma);

// python measure_snr(x, fs, fc, t0_us, t1_us, f_lo_mhz, f_hi_mhz, nfft=1024)
double measure_snr(const cfvec& x, double fs, double fc, double t0_us, double t1_us, double f_lo_mhz,
                   double f_hi_mhz, size_t nfft = 1024);
double measure_snr(const cdvec& x, double fs, double fc, double t0_us, double t1_us, double f_lo_mhz,
                   double f_hi_mhz, size_t nfft = 1024);

}  // namespace si
