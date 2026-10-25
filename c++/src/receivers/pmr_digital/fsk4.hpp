// 4-level FSK land-mobile radio (python: receivers/pmr_digital/fsk4.py): DMR (ETSI TS 102 361-1) and NXDN from the
// FM discriminator output. Dibits: +3 -> 01, +1 -> 00, -1 -> 10, -3 -> 11. Depends on core/ only.
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::pmr_digital::fsk4 {

// python discriminator(y, fs): float32 instantaneous frequency (Hz), median removed
fvec discriminator(const cfvec& y, double fs);
// python symbol_line(d, fs) -> (rate or None, ratio)
std::pair<std::optional<double>, double> symbol_line(const fvec& d, double fs);
// python dmr(d, fs, sr): one record per burst whose SYNC is in the block (python dict, keys in python order)
std::vector<Json> dmr(const fvec& d, double fs, double sr = 4800.0);
// python nxdn(d, fs, sr): (sample position, bit errors) of NXDN frame-sync words
std::vector<std::pair<int64_t, int>> nxdn(const fvec& d, double fs, double sr);
// python frame_consistent(pos, fs, sr)
int64_t frame_consistent(const std::vector<int64_t>& pos, double fs, double sr);

}  // namespace si::rx::pmr_digital::fsk4
