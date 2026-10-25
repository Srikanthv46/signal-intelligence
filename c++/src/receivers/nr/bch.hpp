// NR BCH decoding for L_max = 8 (python: receivers/nr/ssb.py _bch_decode), implemented from
// the standards (no py3gpp code):
//   TS 38.212 5.4.1   rate recovery (repetition combining, sub-block de-interleaving), E = 864, N = 512
//   TS 38.212 5.3.1   polar code: reliability sequence (Table 5.3.1.2-1), frozen set for K = 56,
//                     n_PC = 0, input bit interleaver (Table 5.3.1.1-1, I_IL = 1)
//   successive-cancellation decoding with the min-sum f function (what the Python reference's
//   decoder does; it accepts a list size but decodes a single path), then
//   TS 38.212 5.1     CRC-24C check
//   TS 38.212 7.1.2 / 7.1.1   PBCH payload descrambling and de-interleaving.
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::nr {

using Bits = std::vector<int>;

const std::vector<int>& polar_reliability();      // Q_0^{1023}, Table 5.3.1.2-1 (most reliable last)
std::vector<int> polar_info_positions(int K, int N);  // ascending
std::vector<int> polar_input_interleaver(int K);  // Table 5.3.1.1-1 pattern for K

// TS 38.212 5.4.1 inverse for E >= N: y[k mod N] += llr[k], then sub-block de-interleave
dvec rate_recover_polar(const dvec& llr, int N);
// SC decode (min-sum) of N = 2^nmax LLRs (positive = 0), K info bits incl. CRC, input-deinterleaved
Bits polar_decode(const dvec& rec, int K, int nmax = 9);
// CRC-24C remainder over data + 24 parity bits == 0
bool crc24c_ok(const Bits& blk);
const std::vector<int>& pbch_g();                 // TS 38.212 Table 7.1.1-1 G(j), j < 32

struct BchResult {
  bool crc_ok = false;
  Bits scrblk;    // 32 decoded (still scrambled) bits
  Bits payload;   // 24 MIB bits
  Bits sfn_lsb;   // 4 SFN LSBs
  int hrf = 0;    // half-frame bit
};
// python _bch_decode(softbits, lssb, ncellid) for 864 soft bits
BchResult bch_decode(const dvec& soft, int pci);

}  // namespace si::rx::nr
