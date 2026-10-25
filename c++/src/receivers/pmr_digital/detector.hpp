// Digital land-mobile radio detector (python: receivers/pmr_digital/detector.py): TETRA, DMR, NXDN.
//
//   detect(ch, band, exclude_hz) -> Emission families "tetra", "dmr", "nxdn"
//
// Narrow carriers (1 kHz bins, >= 8 dB, <= 30 kHz) in the land-mobile ranges; TETRA: 18 ksym/s line, pi/4-DQPSK,
// training sequences, BSCH (CRC); 4FSK: DMR SYNC / slot type / BPTC / CSBK CRC, NXDN frame sync.
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::pmr_digital {

const std::vector<std::pair<double, double>>& tetra_ranges();
const std::vector<std::pair<double, double>>& fsk_ranges();

// python detect(ch, {"band": band, "exclude_hz": exclude_hz}); band nullopt = fc +- 0.45 fs
std::vector<Emission> detect(const Channelizer& ch, std::optional<std::pair<double, double>> band = std::nullopt,
                             const dvec& exclude_hz = {});

}  // namespace si::rx::pmr_digital
