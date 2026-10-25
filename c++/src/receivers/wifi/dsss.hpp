// 802.11b DSSS/CCK receiver (python receivers/wifi/dsss.py; IEEE 802.11-2020 Clauses 15/16).
//
// Per 22 MHz channel: Barker-energy detection at 2 samples/chip, SFD-first acquisition and an
// early/late tracked Barker grid at 4 samples/chip, data-aided CFO/phase/SNR from the rebuilt
// scrambled preamble, PLCP header CRC-16 (coherent, differential fallback), PSDU as coherent
// Barker DBPSK/DQPSK (1/2 Mb/s) or equalised CCK (5.5/11 Mb/s), descrambling, FCS with a
// Chase-style correction, STFT interference excision as a last resort, MAC dissection.
//
// scan() returns the PPDU records exactly as the Python dicts (Json, Python key order), so
// they can go to the report verbatim. The building blocks are exposed for the unit tests.
#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/types.hpp"

namespace si::rx::wifi::dsss {

constexpr double kFs = 22e6;      // detection rate: 2 samples per chip
constexpr double kFsDec = 44e6;   // decoding rate: 4 samples per chip
constexpr int kSps = 2;
constexpr int kSym = 11 * kSps;   // detection samples per 1 us symbol
constexpr int kDetK = 16;
constexpr double kDetRatio = 2.8;
constexpr int kDetMinSym = 32;
constexpr double kBarkerReportRho = 0.2;
constexpr int kEqTaps = 9;
constexpr size_t kFixCandidates = 20;
constexpr size_t kFixMaxWeight = 4;
extern const double kBarker[11];

using Bits = std::vector<uint8_t>;

// ---- bit-level helpers (python _descramble, _crc16_ok, _bits_to_int, _bytes, psdu_bytes)
Bits descramble(const Bits& bits);
bool crc16_ok(const uint8_t* bits32, const uint8_t* crc_bits16);
int64_t bits_to_int(const uint8_t* bits, size_t n);
std::vector<uint8_t> bits_to_bytes(const Bits& bits);  // LSB first, len // 8 bytes
int psdu_bytes(double rate, int64_t length_us, int64_t service);

// ---- CCK codebooks (phi1 = 0): rows of 8 chips, and the data bits d2.. of each row
struct CckBook {
  std::vector<cd> C;   // ncode x 8
  std::vector<Bits> bits;
  size_t ncode() const { return bits.size(); }
};
const CckBook& cck11();
const CckBook& cck55();

// ---- detection (python _candidates): [(sample0, sample1, peak coherence)] at 2 sps
struct Candidate {
  int64_t a, b;
  double rho;
};
std::vector<Candidate> candidates(const cfvec& y2);

// ---- demodulation helpers
Bits dbpsk_bits(const cdvec& s);                                       // _dbpsk_bits
std::vector<std::pair<int64_t, std::string>> find_sfds(const Bits& bits);  // _find_sfds
std::pair<Bits, int64_t> rebuild_scrambled(const Bits& raw, const Bits& desc);
Bits diff_bits(const cdvec& s, double rot, int rate_bits);             // _diff_bits

// One 22 MHz channel at 4 samples/chip with its Barker matched filter (python _Chan).
struct Chan {
  cfvec y;
  double fc = 0, fs = 0, T = 0;
  cdvec c;     // Barker MF: c[n] = Barker starting at n
  dvec ac;     // |c|
  cdvec ym;    // chip MF (4 samples), centred at n + 1.5
  Chan(cfvec y4, double ofs4, double fc);
  struct Track {
    cdvec s;
    dvec taus;
    double T;
  };
  Track track(double tau, int64_t nsym, std::optional<double> T = std::nullopt,
              std::optional<double> A = std::nullopt) const;
};

// Every PPDU whose SFD lies in [n0, n1) (python _decode_region); skip(t_start, amp) drops one.
std::vector<Json> decode_region(const Chan& ch, int64_t n0, int64_t n1, int chan, double fc_ch,
                                const std::function<bool(double, double)>& skip = nullptr);

// python scan(ctx): the capture's PPDU records, the candidate count and the Barker-coherent
// regions no PPDU explains (ctx.cache "wifi_dsss", "wifi_dsss_candidates", "wifi_dsss_unexplained").
struct ScanResult {
  std::vector<Json> ppdus;
  size_t n_candidates = 0;
  std::vector<Json> unexplained;
};
ScanResult scan(const Channelizer& ch, double cap_fc, double cap_fs);

namespace detail {
// scipy.signal.oaconvolve(x, k, "full") for a complex64 x and a short float64 kernel,
// bit-identical to sig::oaconvolve (same blocks and transforms), faster (one kernel FFT,
// blocks on nthreads threads, 0 = si::threads())
cdvec oaconvolve_cf_real(const cfvec& x, const dvec& k, int nthreads = 0);
}  // namespace detail

}  // namespace si::rx::wifi::dsss
