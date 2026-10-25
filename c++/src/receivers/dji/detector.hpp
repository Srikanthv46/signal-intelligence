// DJI OcuSync / Lightbridge links (python: receivers/dji/detector.py detect()).
//
//   detect(x, fs, fc) -> one Emission per burst:
//     video   15/30 kHz CP-OFDM, CP metric folded on the 0.5 ms slot over the whole capture,
//             burst segmentation on the known grid, DJI test (15 kHz: length-601/1201 Zadoff-Chu
//             slot symbol in DJI's mapping; 30 kHz O4: reference symbol repeated 4 symbols later);
//     uplink  1.4 MHz 15 kHz-OFDM hops (core/ofdm15k scan + hop tracker), two length-73 ZC symbols
//             of one root 6 symbols apart;
//     droneid 9-symbol 15.36 MSps burst: full CP schedule + ZC roots 600 / 147.
//   DroneID payloads are not decoded in C++ (see README.md): decode_status is what Python
//   reports without its external decoder.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::dji {

// python VIDEO_PROFILES: (name, SCS kHz, output rate, FFT N, normal CP, occupied BW, scan step, ZC length)
struct VideoProfile {
  std::string name;
  int scs;
  double ofs;
  size_t N, cp;
  double occ, step;
  size_t L;  // 0 = None (30 kHz O4: no ZC symbol)
};
const std::vector<VideoProfile>& video_profiles();

// x: complex64 capture at fs, centre fc. decode_droneid as Python ctx["decode_droneid"].
std::vector<Emission> detect(const cfvec& x, double fs, double fc, bool decode_droneid = true);

// python zc_freq(root, L): ZC as DJI maps it (middle value dropped, empty DC)
cdvec zc_freq(int64_t root, size_t L);
// python zc_root_search(X, L, roots=None) -> (corr, root); root -1 = None
std::pair<double, int64_t> zc_root_search(const cdvec& X, size_t L, const std::vector<int64_t>& roots = {});

namespace detail {
// python _FastChannelizer.slice (stateless, thread-safe): X = fft(x) complex128 (notched)
class FastSlicer {
 public:
  FastSlicer(const cdvec& X, double fs) : X_(X), fs_(fs) {}
  cdvec slice(double f_off, double half_bw, double out_fs) const;
  int64_t slice_len(double out_fs) const;  // M = round(L * out_fs / fs)

 private:
  const cdvec& X_;
  double fs_;
};

// python _cp_metric(y, N, cp) -> (|movsum(p)| / movsum(e), movsum(p))
struct CpOut {
  dvec m;
  cdvec cs;
};
// a2: |y|**2 precomputed (optional; shared by profiles that use the same slice)
CpOut cp_metric(const cd* y, size_t len, size_t N, size_t cp, const double* a2 = nullptr);
inline CpOut cp_metric(const cdvec& y, size_t N, size_t cp) { return cp_metric(y.data(), y.size(), N, cp); }
struct FoldIdx {
  std::vector<int64_t> idx, cnt;
};
FoldIdx fold_index(size_t n, double P);       // python _FOLD_IDX entry
dvec fold(const dvec& v, const FoldIdx& f);
dvec fold(const dvec& v, double P);           // python _fold
dvec robust_z(const dvec& prof);       // python _robust_z

struct Segments {
  std::vector<std::pair<size_t, size_t>> runs;
  dvec t, on, off;
};
// python _segment(m, phase, P, n_avg=4, k_sigma=4.0)
Segments segment(const dvec& m, int64_t phase, double P, size_t n_avg = 4, double k_sigma = 4.0);

std::pair<cd, double> zc_dd(const cdvec& X, size_t L);  // python _zc_dd
int64_t root_from_dd(cd acc, size_t L);                 // python _root_from_dd
dvec droneid_cp_score(const cdvec& y);                  // python _droneid_cp_score

// python _frame_period(starts_ms); from_cycle: the value is a Python float (else np.float64)
struct FramePeriod {
  double T = 0;
  bool from_cycle = false;
};
std::optional<FramePeriod> frame_period(const dvec& starts_ms);
}  // namespace detail

}  // namespace si::rx::dji
