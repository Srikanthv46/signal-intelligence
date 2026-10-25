// Transmit streams per subcarrier of a blind CP-OFDM link from one receive channel (python:
// receivers/digital_vtx/streams.py): per-subcarrier kurtosis and the Alamouti pair power correlation on runs of
// contiguous full-band symbols -> 2tx_stbc_like / 2tx_sfbc_like / 1_stream_like / unknown. See streams.py.
#pragma once

#include <vector>

#include "core/types.hpp"
#include "receivers/digital_vtx/blind.hpp"

namespace si::rx::digital_vtx::streams {

constexpr double SNR_MIN = 20.0;

// y: the link's slice (centred, `ofs` samples/s, Tu ~ an integer number of samples); runs: the link's bursts;
// returns the decoded-field object (python tx_streams). y may be empty when snr_db < SNR_MIN.
Json tx_streams(const cdvec& y, double ofs, double tu, double ts, double bw, const std::vector<blind::Run>& runs,
                double snr_db, double dur_s);

}  // namespace si::rx::digital_vtx::streams
