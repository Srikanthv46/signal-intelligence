// Research catalog of drone video / data-link systems (python receivers/digital_vtx/catalog.py): ranked
// system candidates for a blind CP-OFDM link, each with what matched and what conflicted.
#pragma once

#include <optional>
#include <string>

#include "core/types.hpp"

namespace si::rx::digital_vtx::catalog {

// tx_antennas: 1 / 2 / nullopt; frame_ms: nullopt when not measured. -> Json array (python list of dicts)
Json candidates(double tu_us, std::optional<double> cp_frac, std::optional<double> bw_mhz, std::optional<double> carrier_mhz,
                std::optional<int> tx_antennas, std::optional<double> frame_ms, size_t top = 3);

// python model_candidates(system_candidates, carrier_mhz, bw_mhz, tx_antennas): the model block or null
// (research, PLAN M3 / M4)
Json model_candidates(const Json& system_candidates, std::optional<double> carrier_mhz, std::optional<double> bw_mhz,
                      std::optional<int> tx_antennas);

// python hdzero_model(bandwidth_mode): model block for a link the hdzero path identified (family certain, model not)
Json hdzero_model(const std::string& bandwidth_mode);

// python hop_model(hop_bw_mhz, hop_us, period_ms): model block from DJI uplink hop timing (PLAN M6) or null
Json hop_model(std::optional<double> hop_bw_mhz, std::optional<double> hop_us, std::optional<double> period_ms);

}  // namespace si::rx::digital_vtx::catalog
