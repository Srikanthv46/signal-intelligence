// 802.11a/g/n/ac/ax legacy-preamble receiver (python receivers/wifi/ofdm_rx.py): every 20 MHz
// channel of the capture is mixed to baseband and resampled to 20 MSps, the L-STF found by its
// 16-sample autocorrelation plateau, the L-LTF gives timing, fine CFO, the channel and a noise
// estimate per subcarrier, L-SIG is Viterbi-decoded with parity, tail and bit-agreement checks,
// and HT-SIG / VHT-SIG-A are decoded and CRC-8 checked. The undecoded-burst / plotting / CLI
// helpers of the Python module are not part of the pipeline and are not ported.
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::ofdm {

constexpr double kFsWifi = 20e6;

// subcarrier sets (python SC, PILOT_SC, DATA_SC, USED_SC, HT_*)
extern const std::vector<int> kPilotSc;     // -21 -7 7 21
extern const double kPilotBase[4];          // 1 1 1 -1
extern const std::vector<int> kDataSc;      // 48
extern const std::vector<int> kUsedSc;      // 52
extern const std::vector<int> kHtSc;        // -28..28
extern const std::vector<int> kHtDataSc;    // 52
extern const std::vector<int> kHtUsedSc;    // 56
extern const double kLtf[53];               // L-LTF, subcarriers -26..26
extern const double kHtLtf[57];             // HT-LTF, -28..28
const std::vector<double>& pilot_pol();     // p_n, 127 long
inline int bin_of(int sc) { return ((sc % 64) + 64) % 64; }

// python viterbi(soft, terminated): soft > 0 means bit 1, pairs (A, B)
std::vector<uint8_t> viterbi(const dvec& soft, bool terminated = false);
std::vector<uint8_t> conv_encode(const std::vector<uint8_t>& bits);
std::vector<uint8_t> crc8_sig(const uint8_t* bits, size_t n);  // 8 bits, MSB first

// python parse_he_sig_a(bits, lsig_length): the 52 decoded HE-SIG-A bits -> the fields the report uses;
// nullopt when the CRC-4 or the tail fails. L-SIG LENGTH mod 3 = 2 -> MU, else SU (A1 B0 = 1) / TB.
struct HeSigA {
  std::string format;
  bool uplink = false;
  int bss_color = 0, bw = 20, ldpc = 0, dcm = 0;
  std::optional<int> mcs, nsts, txop_us;
};
std::optional<HeSigA> parse_he_sig_a(const std::vector<uint8_t>& bits, int lsig_length);

struct Channel {
  std::string name;
  std::optional<int> ch;
  double freq_hz = 0, offset_hz = 0, usable_half_hz = 1e12;
  bool partial = false;
  bool standard = true;  // python scan_plan "standard": on a standard 20 MHz channel centre
};
struct StdChannel {
  std::string name;
  int ch;
  double freq_hz;
};
const std::vector<StdChannel>& all_channels();
std::vector<Channel> channel_plan(double fc_hz, double fs_hz, double usable_bw_hz, bool allow_partial = true,
                                  double max_missing_hz = 6e6);
std::optional<double> wide_center(const std::string& ch_name, int bw_mhz);
// python scan_plan: candidate 20 MHz centres every step_hz (standard ones carry their channel name)
std::vector<Channel> scan_plan(double fc_hz, double fs_hz, double usable_bw_hz, double step_hz = 1e6,
                               double max_missing_hz = 6e6);

// receiver state a payload decoder needs (python _decode_at's rx dict)
struct RxState {
  int64_t start = 0;
  double cfo = 0;
  int back = 3;
  std::vector<int> used;        // FFT bins of the usable legacy subcarriers
  std::vector<int> good_sc;     // the same as subcarrier indices
  cdvec H;                      // 64, inf where erased
  dvec W_rel, nk;               // 64
  double noise = 0;
  std::shared_ptr<const cdvec> r;  // the 20 MSps baseband of the channel
};

// python Packet (+ the extra fields the plugin reads)
struct Packet {
  std::string channel;
  double ch_freq_mhz = 0, center_mhz = 0;
  int64_t start_sample = 0, end_sample = 0;
  double start_us = 0, duration_us = 0, end_us = 0;
  std::string fmt;
  int lsig_rate_mbps = 0, lsig_length = 0;
  double bandwidth_mhz = 20;
  double f_lo_mhz = 0, f_hi_mhz = 0, snr_db = 0, cfo_khz = 0, stf_metric = 0;
  std::optional<int> mcs;
  std::optional<bool> sig_crc_ok;
  bool truncated = false;
  std::optional<double> energy_end_us;
  // extra
  bool sig_ok = false;  // HT/VHT fields below are set
  int ht_length = 0, aggregation = 0, smoothing = 0, short_gi = 0, ldpc = 0, stbc = 0, n_ess = 0, nsts_su = 0;
  std::optional<int> group_id, partial_aid;
  // HE-SIG-A (set when sig_crc_ok on an HE PPDU): format SU / MU / TB, UL/DL, BSS colour, TXOP; SU: streams, DCM
  std::string he_format;
  bool he_uplink = false;
  int bss_color = 0, dcm = 0;
  std::optional<int> txop_us, he_nsts;
  double lsig_evm_db = 0, lsig_bit_agreement = 0;
  bool off_grid = false;
  // payload state: (rx, channel offset, usable half band)
  std::shared_ptr<const RxState> rx;
  double chan_offset_hz = 0, chan_usable_half_hz = 1e12;
};

struct DetectorConfig {
  double usable_bw = 60e6, stf_thresh = 0.75;
  int min_plateau = 50;
  double power_margin_db = 4.0, min_bit_agreement = 0.85, strong_agreement = 0.90, end_tol_us = 5.0;
};

// python WifiDetector
class WifiDetector {
 public:
  WifiDetector(double fs, double fc, const DetectorConfig& cfg = {},
               std::optional<std::vector<Channel>> channels = std::nullopt);
  // python run(x); channels run on up to nthreads threads (0 = si::threads())
  std::vector<Packet> run(const cfvec& x, int nthreads = 0);
  // python baseband(x, offset_hz): complex64 mix, resample_poly to 20 MSps, as complex128
  cdvec baseband(const cfvec& x, double offset_hz) const;
  size_t n_rejects() const { return n_rejects_; }
  const std::vector<Channel>& channels() const { return channels_; }

  struct Det;  // one decoded preamble (python _decode_at's dict)

 private:
  std::vector<Det> detect_channel(std::shared_ptr<const cdvec> x20, const Channel& chan, size_t& rejects) const;
  bool accept(Packet& p) const;
  double fs_, fc_;
  DetectorConfig cfg_;
  double power_margin_;
  std::vector<Channel> channels_;
  int up_ = 1, down_ = 1;
  double ratio_;
  size_t n_rejects_ = 0;
};

// python offgrid_plan: off-grid blind-scan channels around 20 MHz slices (every 5 MHz) whose 3.2 us CP correlation
// stands out once the known PPDUs are blanked (>= 3 adjacent slices; run midpoint +- 5 MHz)
// Xin: the block's FFT if the caller has it (Channelizer::X()); slices are a power-of-two number of bins <= 20 MHz
std::vector<Channel> offgrid_plan(const cfvec& x, double fs, double fc, const std::vector<Packet>& known,
                                  double usable_bw = 60e6, const cfvec* Xin = nullptr);
// python merge_packets: known (order kept) + the extra PPDUs that are not the same PPDU again
std::vector<Packet> merge_packets(std::vector<Packet> known, std::vector<Packet> extra, double tol_us = 2.0);

namespace detail {
// sig::resample_poly(x, up, down) for complex64, bit-identical, faster (only the kept outputs)
cfvec resample_cf(const cfvec& x, int up, int down);
}  // namespace detail

// Fraction(x).limit_denominator(max_den) -> (numerator, denominator)
std::pair<int64_t, int64_t> limit_denominator(double x, int64_t max_den);

}  // namespace si::rx::wifi::ofdm
