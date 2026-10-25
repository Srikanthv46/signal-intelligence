// 802.11 OFDM data-plane decoder (python receivers/wifi/payload.py): non-HT PPDUs (all eight
// legacy rates) and single-stream 20 MHz HT PPDUs (MCS 0-7, BCC or LDPC, with or without
// Alamouti STBC), from the receiver state of ofdm::WifiDetector, down to FCS-checked MPDUs
// dissected by mac80211. The CLI driver / CSV writer of the Python module is not ported.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"
#include "receivers/wifi/ofdm_rx.hpp"

namespace si::rx::wifi::payload {

using Rate = std::pair<int, int>;

// python deinterleave_indices(ncbps, nbpsc, ncol) -> (j, inv)
std::pair<std::vector<int>, std::vector<int>> deinterleave_indices(int ncbps, int nbpsc, int ncol = 16);
dvec depuncture(const dvec& soft, Rate rate);
std::vector<uint8_t> scrambler_seq(const std::vector<uint8_t>& seed7, size_t n);
// python descramble(bits) -> (descrambled, seed bits)
std::pair<std::vector<uint8_t>, std::vector<uint8_t>> descramble(const std::vector<uint8_t>& bits);
// python llr_exact(sym, nbpsc, var): (n, nbpsc) row-major, log P(0)/P(1)
dvec llr_exact(const cdvec& sym, int nbpsc, const dvec& var);
// python smooth_channel(h, sc, w)
cdvec smooth_channel(const cdvec& h, const std::vector<int>& sc, int w);
// python deaggregate(data) -> [(mpdu, delimiter_crc_ok)]
std::vector<std::pair<std::vector<uint8_t>, bool>> deaggregate(const std::vector<uint8_t>& data);

// python decode_nonht_payload
struct NonHt {
  std::vector<uint8_t> mpdu;
  bool fcs_ok = false;
  std::vector<uint8_t> seed;
  Json info = Json::object();  // mac80211.parse of the body (without FCS), when length >= 4
  bool service_ok = false;
  int n_sym = 0;
};
NonHt decode_nonht_payload(const ofdm::RxState& rx, int rate_mbps, int length);

// python HT_MCS, ht_rate_mbps, ht_decodable, HT_MAX_LENGTH
bool ht_mcs_known(std::optional<int> mcs);
double ht_rate_mbps(int mcs, bool short_gi);
std::pair<bool, std::string> ht_decodable(std::optional<int> mcs, int ldpc, int stbc, int n_ess, double bw_mhz);
constexpr int kHtMaxLength = 8192;

// python decode_ht_payload (raises std::runtime_error with Python's ValueError messages)
struct HtSubframe {
  bool fcs_ok = false;
  std::optional<Json> mac;
};
struct Ht {
  int n_sym = 0;
  std::vector<HtSubframe> subframes;
  int n_fcs_ok = 0;
  int64_t scrambler_seed = 0;
  bool service_ok = false;
  bool ldpc = false;
  int n_codewords = 0, n_codewords_ok = 0;
};
Ht decode_ht_payload(const ofdm::RxState& rx, int mcs, int length, bool short_gi, bool aggregation, int n_ess,
                     double chan_offset_hz, double usable_half_hz, int stbc, int ldpc, int smoothing);

// python vht_uplink_paid(bssid)
int vht_uplink_paid(const std::string& bssid);

}  // namespace si::rx::wifi::payload
