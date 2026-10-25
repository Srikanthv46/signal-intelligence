// IEEE 802.11n/ac LDPC code (python receivers/wifi/ldpc.py, 802.11-2020 19.3.11.7): the twelve
// lifted parity-check matrices, the PPDU codeword parameters (shortening, puncturing,
// repetition) and the normalised min-sum decoder. LLRs are log P(0)/P(1). The prototype
// matrices are those published in Saurabh Tavildar's LDPC library (MIT licence), as in Python.
// The systematic encoder (only the Python test transmitter uses it) is not ported.
#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::ldpc {

using Rate = std::pair<int, int>;

// python params(): codeword count and size, shortening, puncturing, repetition, N_SYM
struct Params {
  int64_t npld = 0, navbits = 0, nsym = 0, ncw = 0, L = 0, K = 0, nshrt = 0, npunc = 0, nrep = 0;
};
Params params(int64_t length_bytes, int64_t ncbps, Rate rate, int64_t m_stbc = 1);

std::vector<int64_t> split(int64_t n, int64_t ncw);  // python _split
struct Seg {
  int64_t k, par, rep;  // information bits sent, parity bits sent, repeated bits
};
std::vector<Seg> layout(const Params& P);          // python _layout
std::vector<int64_t> offsets(const Params& P);     // python offsets
// python transmitted(cw, P, i): the bits codeword i (length L) puts on the air
std::vector<uint8_t> transmitted(const std::vector<uint8_t>& cw, const Params& P, size_t i);

// python decode(llr, L, rate) for one codeword: hard bits, syndrome satisfied, iterations used
struct Decoded {
  std::vector<uint8_t> hard;
  bool ok = false;
  int used = 0;
};
Decoded decode(const dvec& llr, int64_t L, Rate rate, int iters = 60, double alpha = 0.8);

// python decode_stream(llr, P, rate, full=True)
struct StreamDecoded {
  std::vector<uint8_t> info;
  std::vector<bool> ok;
  std::vector<int> used;
  std::vector<std::vector<uint8_t>> hard;  // (ncw, L)
};
StreamDecoded decode_stream(const dvec& llr, const Params& P, Rate rate, int iters = 60);

}  // namespace si::rx::wifi::ldpc
