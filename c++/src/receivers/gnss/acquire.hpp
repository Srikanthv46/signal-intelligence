// GPS L1 C/A acquisition, C/N0, GNSS jamming / spoofing indicators (research, signal gaps S3;
// python: receivers/gnss/codes.py + acquire.py).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/dsp.hpp"
#include "core/types.hpp"

namespace si::rx::gnss {

constexpr double kFL1 = 1575.42e6;

// python codes.ca_code(prn): 1023 chips as 0/1
std::vector<uint8_t> ca_code(int prn);
// python acquire.chi2_sf_even / threshold
double chi2_sf_even(double x, int dof);
double threshold(int m, double p_cell);

// python acquire.swept_ppd(ch, f_l1): a fast swept carrier (GNSS personal privacy device), or nullopt
std::optional<Json> swept_ppd(const Channelizer& ch, double f_l1);
// python acquire.analyse(ch, ppm, exclude_hz): nullopt when L1 is not in the block
std::optional<Json> analyse(const Channelizer& ch, std::optional<double> ppm, const dvec& exclude_hz = {});

}  // namespace si::rx::gnss
