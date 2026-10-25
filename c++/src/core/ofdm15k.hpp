// 15 kHz-SCS OFDM burst scan and hop tracker (python: core/ofdm15k.py): LTE-like links
// such as the DJI RC uplink. Blind CP scan over a (centre, time-window) grid, burst
// refinement (edges, -10 dB bandwidth, CP-phase CFO, CP-peak train check), duplicate
// merge, and the hop-train confirm pass. Used by the dji and emitters receivers.
// ofdm_numerology / classify_numerology / hop_summary have no C++ user and are not ported.
#pragma once

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::ofdm15k {

// slice(f_off, half_bw, out_fs): complex baseband of [f_off +- half_bw] (offset from the
// capture centre) at out_fs, as python Channelizer.slice / dji _FastChannelizer.slice.
using Slicer = std::function<cdvec(double f_off, double half_bw, double out_fs)>;

struct Profile {
  std::string name;
  double ofs;    // output rate
  size_t N;      // FFT size at ofs
  size_t cp;     // normal CP at ofs
  double occ;    // occupied bandwidth
  double step;   // scan step
};
const std::vector<Profile>& profiles();  // PROFILES: 1.4MHz, 10MHz, 20MHz

struct Hit {
  std::string profile;
  double ofs = 0;
  size_t N = 0, cp = 0;
  double occ = 0, f_off = 0, t0_s = 0, t1_s = 0, score = 0, thr = 0;
};
struct Stats {
  double median = 0, sigma = 0;
};
struct Scan {
  std::vector<Hit> hits;
  std::map<std::string, Stats> stats;
};

struct Burst {
  double start_s = 0, end_s = 0, f_off = 0, cfo_hz = 0, bw_hz = 0;
  int64_t nsym = 0;
  double score = 0;
  std::string profile;
  double occ = 0;
  int64_t cp_peaks = 0;
  double period_us = 0, cp_scatter_us = 0;
  bool tracked = false;  // python: b.get("tracked", False)
};

// python ofdm15k_scan(ch, fc, profiles, win_s=0.5e-3, k_sigma=7.0, usable_half=28e6)
Scan ofdm15k_scan(const Slicer& ch, const std::vector<Profile>& profs, double win_s = 0.5e-3, double k_sigma = 7.0,
                  double usable_half = 28e6);
// python refine_ofdm15k(ch, hit, pad_s=1e-3); nullopt = None
std::optional<Burst> refine_ofdm15k(const Slicer& ch, const Hit& hit, double pad_s = 1.0e-3);
// python merge_ofdm15k(bursts): one record per burst, narrowest profile first
std::vector<Burst> merge_ofdm15k(const std::vector<Burst>& bursts);
// python track_hops(ch, fc, bursts, stats, dur_s, k_sigma=4.0, usable_half=28e6)
std::vector<Burst> track_hops(const Slicer& ch, const std::vector<Burst>& bursts,
                              const std::map<std::string, Stats>& stats, double dur_s, double k_sigma = 4.0,
                              double usable_half = 28e6);
// python label_ofdm15k(b) -> (label, why)
std::pair<std::string, std::string> label_ofdm15k(const Burst& b);
// python _cycle(starts_s): repeat period of a burst train, or nullopt
std::optional<double> cycle(const dvec& starts_s);

}  // namespace si::ofdm15k
