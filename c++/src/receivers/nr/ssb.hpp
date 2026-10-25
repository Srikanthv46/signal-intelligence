// 5G NR SSB receiver (python: receivers/nr/ssb.py): GSCN search, PSS, differential SSS,
// DMRS i_SSB, PBCH LLRs, polar decode + CRC-24C (bch.hpp), MIB fields. PSS/SSS m-sequences
// from TS 38.211 7.4.2.2 / 7.4.2.3, PBCH DMRS from 7.4.1.4.1. Depends on core/ only.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/types.hpp"

namespace si::rx::nr {

constexpr int kNfft = 256;

struct NrBand {
  std::string key, name;
  double scs;
};
const std::vector<NrBand>& nr_bands();                       // python NR_BANDS, dict order
std::vector<std::pair<double, double>> bands_for(const std::string& key);  // python _bands_for
dvec gscn(double lo, double hi);                              // python _gscn

dvec pss_seq(int nid2);   // TS 38.211 7.4.2.2.1, +-1 (127)
dvec sss_seq(int pci);    // TS 38.211 7.4.2.3.1, +-1 (127)
cdvec pbch_dmrs(int pci, int ib);  // TS 38.211 7.4.1.4.1 (144 QPSK symbols)

// conj(fft(PSS templates, nf)) for one slice rate (python _pss_templates / _templates)
struct TemplateBank {
  size_t nf = 0;
  double ofs = 0;
  Mat<cd> T;  // 3 x nf
};
std::shared_ptr<const TemplateBank> template_bank(double scs, double ofs, size_t nf);
// the slice length / rate / FFT size Channelizer::slice gives for out_fs (no slicing done)
struct SliceDims {
  size_t L, nf;
  double ofs;
};
SliceDims slice_dims(const Channelizer& ch, double out_fs);

struct PssBest {
  double z = 0, df = 0, metric = 0;
  int nid2 = 0;
  int64_t pos = 0;
};

class SSBSlice {
 public:
  SSBSlice(const Channelizer& ch, double f_ss, double scs, std::shared_ptr<const TemplateBank> bank = nullptr);
  PssBest best(const dvec& offsets) const;                  // python _SSBSlice.best (without m)
  dvec metric(double df, int nid2) const;                   // |ifft(roll(Y) * T[nid2])|^2 / e
  cdvec shifted(double df) const;                           // y * exp(-2j pi df n / ofs)
  std::optional<Mat<cd>> grid(const cdvec& shifted, int64_t pos) const;  // 4 x 240

  double scs = 0, ofs = 0;
  int cp = 0;
  cfvec y;
  size_t L = 0, nf = 0;

 private:
  void metric_into(double df, int nid2, cdvec& w, dvec& m) const;
  cfvec Y_;
  dvec e_;
  std::shared_ptr<const TemplateBank> bank_;
};

std::vector<int64_t> peaks(const dvec& m, double frac = 0.2, int64_t sep = kNfft, size_t max_n = 24);
dvec sss_scores(const Mat<cd>& G, int nid2);                       // python _sss_scores (336)
std::pair<int, dvec> pbch_llr(const Mat<cd>& G, int pci);           // python _pbch_llr: (i_SSB, 864 LLRs)

struct OtherCell {
  int pci;
  double sss_z;
};
struct SsbInfo {
  size_t n_ssb_instances = 0;
  std::optional<int> nid1, ibar_ssb;
  int nid2 = 0;
  double sss_z = 0, sss_z2 = 0;
  std::vector<OtherCell> other_cells;
  std::vector<int> ssb_indices;
  bool crc_ok = false;
  int lmax = 0, n_combined = 0, decode_attempts = 0;
  Json mib;
};
// python decode_ssb_multi(sl, b): (pci or nullopt, info)
std::pair<std::optional<int>, SsbInfo> decode_ssb_multi(const SSBSlice& sl, const PssBest& b);

}  // namespace si::rx::nr
