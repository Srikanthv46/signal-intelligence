// Analog FPV video (FM PAL / NTSC composite) detector (python: receivers/analog_vtx/detector.py).
//
//   detect(x, fs, fc) -> Emission family "analog_vtx", label "pal" | "ntsc" | "analog_fm_video"
//   video_frame(x, fs, fc, em) -> uint8 picture [lines, 720] (or nullptr)
//
// A line-period fold of a ~1 us x ~1 MHz time-frequency picture finds the sync tips far
// below the FM threshold; the emitter is then characterised (line rate, sync-tip frequency,
// occupied band, field timing, channel table, FM-discriminator features). Depends on core/ only.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/spectral.hpp"
#include "core/types.hpp"

namespace si::rx::analog_vtx {

constexpr double kSyncS = 4.7e-6;
constexpr double kThreshold = 10.0;
constexpr double kPpmMax = 500.0;
constexpr int kNB = 128;
constexpr float kClip = 8.0f;
constexpr int kMaxEmitters = 4;
constexpr double kDevCorrHz = 1.0e6;
constexpr size_t kMaxSamples = 4000000;

struct Nominal {
  const char* name;
  double fl;
  int lines;
  double fsc;
};
// python NOMINAL, in its (dict) order: pal, ntsc
const std::vector<Nominal>& nominal();
const Nominal& nominal(const std::string& std_name);

// python match_channel(f_mhz, tol_mhz): (name, nominal MHz) nearest first
std::vector<std::pair<std::string, int>> match_channel(double f_mhz, double tol_mhz = 3.0);

// ---------------------------------------------------------------- time-frequency fold
struct TfPower {
  fmat q, P;      // median-normalised clipped power, raw power (frames x nfft), float32
  dvec t, f;      // frame times (s), bin offsets (Hz)
};
// python tf_power(x[:n], fs)
TfPower tf_power(const cf* x, size_t n, double fs);

struct Fold {
  dmat F;         // nb x bins
  dvec cnt;
};
// python fold(q, t, rate, nb): float32 accumulation in time order (scipy sparse product), / count
Fold fold(const fmat& q, const dvec& t, double rate, int nb = kNB);

// python zscore(F, q_sd, n_per_bin, nb_box) -> (zb, zs)
std::pair<dmat, dmat> zscore(const dmat& F, const fvec& q_sd, double n_per_bin, int nb_box);

struct ScanBest {
  double score = 0;
  std::string std_name;
  double ppm = 0, rate = 0;
  int b0 = 0, k0 = 0;
  dmat zb;
};

// python scan(): every (standard, ppm) hypothesis is folded once (in parallel) and its
// z-score map kept, so the repeated scans of detect() (with growing column masks) only
// redo the argmax.
class Scanner {
 public:
  Scanner(const TfPower& tf, double ppm_max = kPpmMax);
  ScanBest best(const std::vector<uint8_t>* col_mask) const;
  size_t n_hypotheses() const { return hyp_.size(); }
  // python best["ppm_profile"]: the unmasked maximum per hypothesis (pal ppms, then ntsc)
  dvec profile() const;

 private:
  struct Hyp {
    int s;          // index into nominal()
    double ppm, rate;
    dmat zs;
  };
  const TfPower& tf_;
  fvec q_sd_;
  std::vector<Hyp> hyp_;
};

// python emitter_columns(zb, b0, k0, f, nbins_win, thr)
std::pair<int, int> emitter_columns(const dmat& zb, int b0, int k0, int nbins_win, double thr = 3.5);

// ---------------------------------------------------------------- refinement
// python slice_band(ch, f_off, half_bw, out_fs): bins mapped by signed frequency
struct Band {
  cdvec y;
  double fs = 0;
};
Band slice_band(const SpectralChannelizer& ch, double f_off, double half_bw, double out_fs);

double grid_line_rate(const dvec& e, double efs, double rate, double span_ppm = 160.0);
std::pair<double, double> drift_line_rate(const dvec& e, double efs, double rate, int nblk, int nbf = 1024,
                                          int iters = 2, int smooth = 8);
// python sync_tip_freq -> (f, prominence, f_is_numpy)
struct SyncTip {
  double f = 0, prom = 0;
};
SyncTip sync_tip_freq(const SpectralChannelizer& ch, double fs, double f_guess, double rate, double t0,
                      double half = 2.0e6);

struct Timing {
  double rate = 0, sig = 0, t0 = 0, f_sync = 0, sync_prom = 0;
  bool f_sync_np = false;   // f_sync is a numpy float64 (rounding semantics in the report)
  dvec e;
  double efs = 0;
};
Timing refine_timing(const SpectralChannelizer& ch, double fs, double f_sync, double rate, double score,
                     double eclip = 24.0);

struct Vbi {
  int64_t k;
  double z, t;
};
std::vector<Vbi> field_timing(const dvec& e, double efs, double rate, double t0, const std::string& std_name);

struct Spec {
  double centre = 0, e_lo = 0, e_hi = 0, rms = 0, snr_db = 0, f_white = 0;
  int lo = 0, hi = 0;
  bool edge_cut = false;
  std::string source;
};
Spec spectrum_features(const fmat& P, const dvec& f, int k_sync, double max_half_hz = 25e6);
Spec line_locked_spectrum(const fmat& P, const dvec& t, const dvec& f, double rate, double t0, int k_sync,
                          double max_half_hz = 25e6, double clip = 50.0);

// python fm_discriminate(ch, f_off, half_bw, fs) -> (d, ofs)
struct Disc {
  dvec d;
  double fs = 0;
};
Disc fm_discriminate(const SpectralChannelizer& ch, double f_off, double half_bw, double fs);

// python chroma_standard(): colour standard from the line-locked burst phasors (research)
Json chroma_standard(const dvec& d, double ofs, double rate, double t0, const std::string& std_name);
// python audio_features(d, ofs, fa): audio FM subcarrier -> tone / deviation / flatness (research)
Json audio_features(const dvec& d, double ofs, double fa);
// python secam_test(): SECAM line-alternating back-porch subcarrier (research)
Json secam_test(const dvec& d, double ofs, double rate, double t0);
// python demod_features(): ordered as Python inserts them; values unrounded
Json demod_features(const dvec& d, double ofs, double rate, double t0, const std::string& std_name, int polarity);

// ---------------------------------------------------------------- entry points
// python detect(x, fs, fc, ctx={}) (builds its own spur-notched channelizer)
std::vector<Emission> detect(const cfvec& x, double fs, double fc);

// python video_frame(x, fs, fc, em, ncols=720, max_lines=1250)
std::shared_ptr<Image> video_frame(const cfvec& x, double fs, double fc, const Emission& em, int ncols = 720,
                                   int max_lines = 1250, const SpectralChannelizer* ch = nullptr);

}  // namespace si::rx::analog_vtx
