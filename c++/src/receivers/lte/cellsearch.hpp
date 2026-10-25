// LTE downlink cell search (python: receivers/lte/cellsearch.py): PSS fold over 5 ms, ZC
// ambiguity NMS, SSS (FDD/TDD, normal/extended CP) -> PCI, receiver ppm search, occupied
// bandwidth from energy. Depends on core/ only.
#pragma once

#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/spectrogram.hpp"
#include "core/types.hpp"

namespace si::rx::lte {

constexpr double kOfs = 3.84e6;           // slice rate
constexpr int kN = 256;                   // 15 kHz useful symbol at 3.84 MSps
constexpr int kCpN = 18, kCp0N = 20, kCpE = 64;
constexpr int kRoots[3] = {25, 29, 34};
constexpr int64_t kHalf = 19200;          // int(5e-3 * OFS)
constexpr double kZMin = 15.0;
constexpr int kSssSmooth = 15;            // channel-estimate smoothing (subcarriers) per side of DC

// python LTE_BANDS (band key -> name), in dict order
const std::vector<std::pair<std::string, std::string>>& lte_bands();
// python BWS / OCC
const std::vector<double>& bws();
double occ(double bw);

// One PSS hypothesis result (python dict z, df, nid2, phase, m, iso).
struct PssBest {
  double z = 0, df = 0, iso = 0;
  int nid2 = 0;
  int64_t phase = 0;
  dvec m;                                 // normalised correlation (filled when asked for)
};

// conj(fft(templates, nf)) of the three PSS roots (python _Slice.T), shared by slices of one capture
struct TemplateBank {
  size_t nf = 0;
  Mat<cd> T;                              // 3 x nf
};
std::shared_ptr<const TemplateBank> template_bank(size_t nf);

// One 3.84 MSps slice; frequency offsets are integer-bin rotations of its FFT (python _Slice).
class Slice {
 public:
  Slice(const Channelizer& ch, double f_center, std::shared_ptr<const TemplateBank> bank = nullptr);
  // python _Slice.best(offsets): first maximum of z over (offset, root). want_m fills .m.
  PssBest best(const dvec& offsets, bool want_m, bool parallel = false) const;
  cdvec shifted(double df) const;         // y * exp(-2j pi df n / ofs), complex128

  double f_center = 0, ofs = 0;
  cfvec y;
  size_t L = 0, nf = 0;
  int64_t k = 0;

 private:
  struct Score {
    double z, iso, dfq;
    int64_t k;
  };
  Score score(double df, int r, cdvec& work, dvec& m) const;
  cdvec Y_;
  dvec e_;
  std::shared_ptr<const TemplateBank> bank_;
};

struct SssResult {
  std::string duplex, cp;
  int nid1 = 0, parity = 0, n = 0;
  double z = 0, z2 = 0;
};
using Instances = std::vector<std::pair<int64_t, int64_t>>;  // (half-frame index, sample)

// python _instances(m, phase)
Instances instances(const dvec& m, int64_t phase);
// python sss_decode(y, nid2, inst)
std::optional<SssResult> sss_decode(const cdvec& y, int nid2, const Instances& inst);
// python sss_decisive(s)
bool sss_decisive(const std::optional<SssResult>& s);

// python wide_ppm_search(ch, centres): {"ppm", "z", "carrier_hz", "raster_hz"} or nullopt
std::optional<Json> wide_ppm_search(const Channelizer& ch, const dvec& centres, double span = 500e3,
                                    double step = 5e3);
// python bandwidth_from_energy(det, fcar): channel bandwidth (Hz) or nullopt
std::optional<double> bandwidth_from_energy(const DetectResult& det, double fcar);

struct Cell {
  double raster_hz = 0, carrier_hz = 0, pss_z = 0, pss_isolation = 0, err_hz = 0;
  int nid2 = 0;
  Instances inst;
  std::optional<SssResult> sss;
  bool sss_ok = false;
  int64_t phase = 0;
};
// python raster_grid(sig): EARFCN raster carriers whose PSS (+-0.55 MHz) lies inside the energy signal
dvec raster_grid(const Signal& sig);
// python raster_offsets(f_hz, ppm_known)
dvec raster_offsets(double f_hz, bool ppm_known);
// python search_cells(ctx, sig, ppm, done, ppm_known, require_sss). done: raster keys round(r / 100 kHz) already
// searched for another energy signal (skipped here and updated); nullptr = search all.
std::vector<Cell> search_cells(const Channelizer& ch, const Signal& sig, double ppm, std::set<int64_t>* done = nullptr,
                               bool ppm_known = false, bool require_sss = false);

}  // namespace si::rx::lte
