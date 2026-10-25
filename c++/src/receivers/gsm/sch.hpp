// GSM BCCH carrier receiver (python: receivers/gsm/sch.py): FCCH tone search, SCH
// differential-GMSK demodulation, Viterbi, CRC-10, BSIC / frame number fields.
// Signals are at kOfs (4 samples per symbol). Depends on core/ only.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::gsm {

constexpr double kRsym = 1625e3 / 6;  // 270.833 ksym/s
constexpr int kSps = 4;
constexpr double kOfs = kRsym * kSps;
constexpr double kFFcch = kRsym / 4;  // +67.708 kHz

const std::vector<std::string>& gsm_bands();  // python GSM_BANDS
// python _arfcn_freq(f): (ARFCN, downlink carrier Hz) or nullopt
std::optional<std::pair<int, double>> arfcn_freq(double f);

struct Fcch {
  int64_t start;
  double cfo, purity;
};
constexpr int64_t kFcchWin = static_cast<int64_t>(0.4e-3 * kOfs);  // python FCCH_WIN
constexpr size_t kFcchNfft = 2048;                                   // python FCCH_NFFT
// python find_fcch(y): complex64 (slice) or complex128 (CFO-corrected) input
std::optional<Fcch> find_fcch(const cfvec& y, int64_t win = kFcchWin);
std::optional<Fcch> find_fcch(const cdvec& y, int64_t win = kFcchWin);
// python find_fcch_cfo(y, max_cfo): FCCH with an unknown carrier offset (|cfo| <= max_cfo)
std::optional<Fcch> find_fcch_cfo(const cfvec& y, double max_cfo);

struct Sch {
  dvec soft;   // 78 coded-bit soft values, positive = 0
  int ts_match = 0;
  int64_t sample = 0;
};
// python demod_sch(y, start) on the CFO-corrected (complex128) signal
std::optional<Sch> demod_sch(const cdvec& y, int64_t start);
std::vector<int> viterbi(const dvec& soft);       // 39 bits
bool crc10_ok(const std::vector<int>& bits25, const std::vector<int>& par10);
Json sch_fields(const std::vector<int>& b);       // python sch_fields (dict order)

}  // namespace si::rx::gsm
