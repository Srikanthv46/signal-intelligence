// scipy.signal equivalents (scipy 1.13 semantics) for the calls the receivers make:
// find_peaks, firwin, butter (sos), sosfilt / sosfiltfilt / sosfilt_zi, lfilter /
// lfilter_zi, fftconvolve / oaconvolve, resample_poly, welch, stft / istft, get_window.
//
// Precision follows scipy's dtype rules: float32 / complex64 inputs are processed in
// single precision where scipy does (FFT convolution, resample_poly, welch), stft works in
// double and rounds to complex64 like scipy, IIR filtering is always double.
#pragma once

#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "core/npcompat.hpp"
#include "core/types.hpp"

namespace si::sig {

// ------------------------------------------------------------------ windows
struct Window {
  enum Kind { Hann, Hamming, Kaiser, Boxcar } kind = Hann;
  double beta = 0.0;  // Kaiser
  static Window hann() { return {Hann, 0}; }
  static Window hamming() { return {Hamming, 0}; }
  static Window kaiser(double b) { return {Kaiser, b}; }
  static Window boxcar() { return {Boxcar, 0}; }
};
// scipy.signal.get_window(window, N, fftbins): fftbins=true -> periodic (DFT-even)
dvec get_window(const Window& w, size_t N, bool fftbins = true);
// scipy.special.i0 (Cephes)
double i0(double x);

// ------------------------------------------------------------------ find_peaks
struct PeakOpts {
  std::optional<double> height_min, height_max;        // height=h -> min only; (None, h) -> max only
  std::optional<double> threshold_min, threshold_max;
  std::optional<double> distance;                      // >= 1
  bool prominence = false;                             // prominence given (even (None, None))
  std::optional<double> prominence_min, prominence_max;
  bool width = false;                                  // width given
  std::optional<double> width_min, width_max;
  double wlen = -1;                                    // None
  double rel_height = 0.5;
};
struct Peaks {
  std::vector<int64_t> peaks;
  // properties, filled when the matching condition was given (as scipy's dict)
  dvec peak_heights, left_thresholds, right_thresholds;
  dvec prominences;
  std::vector<int64_t> left_bases, right_bases;
  dvec widths, width_heights, left_ips, right_ips;
};
// scipy.signal.find_peaks(x, ...) (x converted to float64 as scipy does)
Peaks find_peaks(const dvec& x, const PeakOpts& opt = {});
Peaks find_peaks(const fvec& x, const PeakOpts& opt = {});

// ------------------------------------------------------------------ FIR / IIR design
// scipy.signal.firwin(numtaps, cutoff, window=..., pass_zero, scale, fs); fs = None ->
// cutoff relative to Nyquist = 1.
dvec firwin(size_t numtaps, const dvec& cutoff, std::optional<double> fs = std::nullopt,
            const Window& window = Window::hamming(), bool pass_zero = true, bool scale = true);
inline dvec firwin(size_t numtaps, double cutoff, std::optional<double> fs = std::nullopt,
                   const Window& window = Window::hamming(), bool pass_zero = true, bool scale = true) {
  return firwin(numtaps, dvec{cutoff}, fs, window, pass_zero, scale);
}

enum class BType { Lowpass, Highpass, Bandpass };
// scipy.signal.butter(N, Wn, btype, fs=fs, output="sos"): (n_sections x 6) rows
// [b0 b1 b2 1 a1 a2], scipy's zpk2sos pairing ('nearest') and section order.
dmat butter_sos(int N, const dvec& Wn, BType btype = BType::Lowpass, std::optional<double> fs = std::nullopt);
inline dmat butter_sos(int N, double Wn, BType btype = BType::Lowpass, std::optional<double> fs = std::nullopt) {
  return butter_sos(N, dvec{Wn}, btype, fs);
}

// ------------------------------------------------------------------ filtering (float64)
// scipy.signal.sosfilt(sos, x, zi): zi (n_sections x 2) used as initial state and
// replaced by the final state when given.
dvec sosfilt(const dmat& sos, const dvec& x, dmat* zi = nullptr);
// scipy.signal.sosfilt_zi(sos)
dmat sosfilt_zi(const dmat& sos);
// scipy.signal.sosfiltfilt(sos, x, padtype='odd', padlen=None). The odd extension is
// built in the input dtype (float32 stays float32), filtering is double.
enum class PadType { Odd, Even, Constant, None };
template <class T>
dvec sosfiltfilt(const dmat& sos, const std::vector<T>& x, PadType padtype = PadType::Odd, long padlen = -1);
// scipy.signal.lfilter_zi(b, a)
dvec lfilter_zi(const dvec& b, const dvec& a);
// scipy.signal.lfilter(b, a, x, zi): FIR (len(a) == 1) via np.convolve as scipy, IIR
// via the direct form II transposed recursion. zi in/out (length max(len a, len b) - 1).
dvec lfilter(const dvec& b, const dvec& a, const dvec& x, dvec* zi = nullptr);

// ------------------------------------------------------------------ FFT convolution
// scipy.signal.fftconvolve / oaconvolve (1-D), mode full / same / valid. Each input is
// transformed in its own precision, the product and inverse in the promoted type.
// Supported pairs: (float, float), (double, double), (cf, cf), (cd, cd), (cd, double),
// (double, cd), (cf, double) -> cd, (cf, float) -> cf, (float, cf) -> cf.
template <class A, class B>
struct conv_out {
  using R = std::conditional_t<std::is_same_v<np::real_t<A>, double> || std::is_same_v<np::real_t<B>, double>, double, float>;
  using type = std::conditional_t<np::is_cplx<A>::value || np::is_cplx<B>::value, std::complex<R>, R>;
};
template <class A, class B>
using conv_out_t = typename conv_out<A, B>::type;
template <class A, class B>
std::vector<conv_out_t<A, B>> fftconvolve(const std::vector<A>& in1, const std::vector<B>& in2,
                                          np::ConvMode mode = np::ConvMode::Full);
template <class A, class B>
std::vector<conv_out_t<A, B>> oaconvolve(const std::vector<A>& in1, const std::vector<B>& in2,
                                         np::ConvMode mode = np::ConvMode::Full);
size_t next_fast_len(size_t n, bool real);  // scipy.fft.next_fast_len(n, real)

// ------------------------------------------------------------------ resampling
// scipy.signal.resample_poly(x, up, down) with the default window ('kaiser', 5.0) and
// padtype 'constant'. T = float, double, cf, cd (output same dtype).
template <class T>
std::vector<T> resample_poly(const std::vector<T>& x, int up, int down, const Window& window = Window::kaiser(5.0));

// ------------------------------------------------------------------ spectra
template <class R>
struct Welch {
  dvec f;
  std::vector<R> Pxx;
};
// scipy.signal.welch(x, fs, nperseg=..., noverlap=None, return_onesided=...) with the
// defaults hann / detrend='constant' / scaling='density' / average='mean'. Complex input
// is always two-sided. Output float32 for float / cf input, float64 otherwise.
template <class T>
Welch<np::real_t<T>> welch(const std::vector<T>& x, double fs, size_t nperseg, long noverlap = -1,
                           bool return_onesided = true);

template <class C>
struct Stft {
  dvec f, t;
  Mat<C> Z;  // (nfft frequencies) x (segments), as scipy's Zxx
};
// scipy.signal.stft(x, fs, nperseg, noverlap, return_onesided=False) for complex x
// (hann window, boundary='zeros', padded=True, scaling='spectrum').
template <class C>
Stft<C> stft(const std::vector<C>& x, double fs, size_t nperseg, long noverlap = -1);
template <class C>
struct Istft {
  dvec t;
  std::vector<C> x;
};
// scipy.signal.istft(Z, fs, nperseg, noverlap, input_onesided=False) (hann, boundary=True)
template <class C>
Istft<C> istft(const Mat<C>& Z, double fs, size_t nperseg, long noverlap = -1);

}  // namespace si::sig
