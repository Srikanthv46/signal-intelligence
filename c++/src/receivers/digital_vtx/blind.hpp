// Vendor-independent CP-OFDM link detector (python: receivers/digital_vtx/blind.py): a coherent lag-Tu
// autocorrelation scan of the band (time-AGC'd), symbol-rate line / per-window fold confirmation,
// CP cross-spectrum band, bursts on local symbol grids, TDD frame completion, band and burst extent from
// in-band power; hopping narrowband CP-OFDM links (RC uplinks). See blind.py.
#pragma once

#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "core/bursts.hpp"
#include "core/types.hpp"

namespace si::rx::digital_vtx::blind {

// sl(f_off, half_bw, out_fs, &actual_rate) -> slice (python: the caller's _slice)
using Slicer = std::function<cdvec(double, double, double, double*)>;

constexpr double TU_MIN = 10e-6;

// python rate_for(tu, B, cp) -> (rate, k)
std::pair<double, int64_t> rate_for(double tu, double B, double cp = 0.0);

struct Hit {
  double W = 0, off = 0;
  int64_t N = 0;
  double ofs = 0, tu = 0, z = 0;
};
// skip(off): true leaves a candidate / scan hit at that offset unconfirmed (cheap band gate, PLAN G8); empty = none
using SkipFn = std::function<bool(double)>;
// python scan(sl, usable, skip): a slice whose whole band (sampled every 1 MHz) is skipped is not scanned
std::vector<Hit> scan(const Slicer& sl, double usable, const SkipFn& skip = {});

struct Confirmed {
  double tu = 0, ts = 0, line_z = 0;
  int64_t nwin = 0;
  bool ts_np = true;  // ts is a numpy float in python (decides how round() treats it downstream)
};
std::optional<Confirmed> confirm(const Slicer& sl, double off, double W, int64_t N0, double ofs0);

struct Run {
  double t0 = 0, t1 = 0;
  int64_t n = 0;
};
struct Segmented {
  std::vector<Run> runs;
  dvec tk, zr;  // symbol times (s), per-symbol z
  double ts = 0;
  int64_t ngrid = 0;
  double frac = 0;
};
Segmented segment(const Slicer& sl, double off, double B, double tu, double ts);

// -> (T_s or nullopt, contrast)
std::pair<std::optional<double>, double> frame_period(const dvec& tk, const dvec& zr, double ts, double lo = 0.9e-3,
                                                      double hi = 12e-3);
std::vector<Run> frame_complete(const std::vector<Run>& runs, const dvec& tk, const dvec& zr, double ts, double T,
                                double min_con = 6.0);
// -> (frequency offsets, CP cross-spectrum per symbol)
std::pair<dvec, dvec> band_cp(const Slicer& sl, double off, double W, double tu, double ts, const dvec& tsym);
std::optional<std::pair<double, double>> measure_band(const Slicer& sl, double off, double W, double tu, double ts,
                                                      const dvec& tsym, double usable);

struct Link {
  double tu = 0, ts = 0, lo = 0, hi = 0, cen = 0, bw = 0, line_z = 0, acf_z = 0;
  int64_t win_z = 0;
  std::vector<Run> runs;
  int64_t n_grids = 0;
  std::optional<double> frame_s;
  double frame_contrast = 0, snr_db = 0, p_on = 0, duty = 0, cfo_frac_hz = 0;
  bool ts_np = true;
  cdvec slice;              // kept when snr_db >= keep_slice_snr: centred, band B, rate slice_rate
  double slice_rate = 0;
};
// tf: the capture's core::bursts::TFPower (band and burst extent from in-band power), or nullptr;
// keep_slice_snr: links at or above this SNR keep their slice (python lk['slice'])
std::vector<Link> blind_links(const Slicer& sl, double fs, double fc, double usable, double noise_pps,
                              const bursts::TFPower* tf = nullptr, double keep_slice_snr = 1e300,
                              const SkipFn& skip = {});

// ---- hopping narrowband CP-OFDM links (RC uplinks), python blind_hops
struct Hop {
  double t0 = 0, t1 = 0, lo = 0, hi = 0, z = 0, contrast = 0, channel = 0;
  bool tracked = false;
  size_t num = 0;
};
struct HopLink {
  double tu = 0, ts = 0, cpf = 0;
  std::vector<Hop> hops;
  int64_t n_channels = 0;
  std::optional<double> frame_s;
  dvec paired;  // centre offsets of the blind links of the same numerology
  double dur_s = 0, bw = 0;
};
double hop_width(double tu);
// python refine_hop(sl, f_off, t0_s, t1_s, tu, cpf, W); nullopt = None
std::optional<Hop> refine_hop(const Slicer& sl, double f_off, double t0_s, double t1_s, double tu, double cpf, double W);
std::optional<double> hop_period(const dvec& starts, double lo = 1e-3, double hi = 20e-3, double tol = 30e-6);
std::vector<HopLink> blind_hops(const Slicer& sl, double fs, double fc, double usable, const std::vector<Link>& links,
                                const SkipFn& skip = {});

}  // namespace si::rx::digital_vtx::blind
