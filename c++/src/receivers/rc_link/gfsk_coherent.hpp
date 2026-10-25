// Coherent detection of h ~ 0.5 GFSK (BLE 1M) for low SNR + CRC-24 Chase correction
// (python receivers/rc_link/gfsk_coherent.py; research, synthetic validation only, PLAN L1).
#pragma once

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::gfsk {

constexpr int kSps = 8;  // samples per bit (the BLE slice is 8 MSps)

// python diff_products(y): LPF + one-bit differential products p[i] = z[i + T] z*[i]
cdvec diff_products(const cfvec& y);
// python find_access_address(p, aa_bits, top=3, min_ncorr=0.55): [(sample, normalised correlation)]
std::vector<std::pair<int64_t, double>> find_access_address(const cdvec& p, const std::vector<uint8_t>& aa_bits,
                                                            int top = 3, double min_ncorr = 0.55);

struct Demod {
  std::vector<uint8_t> bits;
  dvec rel;          // reliability of each b_k
  double h = 0.5;    // modulation index
  double cfo_hz = 0;
  int dt = 0;
};
// python demodulate(y, s0, n, iters=2)
std::optional<Demod> demodulate(const cfvec& y, int64_t s0, int64_t n, int iters = 2);

// python crc_syndrome_parts(dw, n_pdu) -> (syndrome, contributions)
std::pair<int64_t, std::vector<int64_t>> crc_syndrome_parts(const std::vector<uint8_t>& dw, size_t n_pdu);
// python chase(dw, n_pdu, rel, off, k_least=16, max_units=3, protect=range(plo, phi)) -> (bits, units) or nullopt
std::optional<std::pair<std::vector<uint8_t>, int>> chase(const std::vector<uint8_t>& dw, size_t n_pdu,
                                                          const dvec& rel, int64_t off, size_t plo = 8,
                                                          size_t phi = 16, int k_least = 16, int max_units = 3);

}  // namespace si::rx::rc_link::gfsk
