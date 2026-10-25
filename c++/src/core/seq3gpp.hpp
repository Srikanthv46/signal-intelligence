// 3GPP synchronisation / scrambling sequences (python: core/seq3gpp.py): LTE PSS
// (Zadoff-Chu), SSS (m-sequences), the 62 PSS/SSS subcarriers, the SSS bank, the Gold
// sequence c(n) of 36.211 7.2 / 38.211 5.2.1 and QPSK mapping.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::seq3gpp {

// 36.211 6.11.1.1 PSS: 62-length ZC with root u in {25, 29, 34} (N_ID2 0/1/2)
cdvec lte_pss_freq(int u);
// 36.211 6.11.2.1 SSS for subframe 0 (sub = 0) or 5 (sub = 5): +-1 values
dvec lte_sss_freq(int nid1, int nid2, int sub);
// _SC62: subcarrier indices -31..-1, 1..31 (DC skipped)
const std::vector<int>& sc62();
// _SSS_BANK[nid2][s][nid1][k], s = 0 (subframe 0) / 1 (subframe 5), k < 62
using SssBank = std::vector<std::array<std::vector<dvec>, 2>>;  // [nid2][s][nid1] -> 62 values
const SssBank& sss_bank();
// m-sequences _S, _C, _Z (+-1, length 31)
const std::vector<int>& mseq_s();
const std::vector<int>& mseq_c();
const std::vector<int>& mseq_z();

// gold(c_init, n, nc=1600): c(0..n-1) as 0/1 int8 (bit-identical to py3gpp nrPRBS)
std::vector<int8_t> gold(uint32_t c_init, size_t n, size_t nc = 1600);
// qpsk(bits): ((1-2b(2i)) + j(1-2b(2i+1))) / sqrt(2) (numpy: times the rounded 1/sqrt(2))
cdvec qpsk(const std::vector<int8_t>& bits);

}  // namespace si::seq3gpp
