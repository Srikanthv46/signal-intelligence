// 2 Mb/s GFSK packet links of the Enhanced ShockBurst family (python: receivers/rc_link/esb.py): per-burst bit
// slicing behind a 0101 preamble, and the joint CRC-16/CCITT check of the packets of one link (same sync word),
// which finds the CRC span and init from the packets alone.
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::rc_link::esb {

constexpr double kRate = 2e6;            // bit/s
constexpr double kHalfBw = 1.25e6;       // slice half-width
constexpr double kOutFs = 16e6;          // 8 samples per bit
constexpr int kMinPreamble = 14;         // alternating bits (a 16-bit preamble measures 14-15)
constexpr int kPreambleStart = 32;       // ... starting within the first 32 bits of the burst
constexpr int kSyncBits = 16;            // packets of one link share their first 16 bits after the preamble
constexpr uint32_t kCrcPoly = 0x1021;
constexpr int kMinSpan = 24;             // CRC span >= address
constexpr int kMinVerified = 3;          // distinct packets on one (E, K)
constexpr uint32_t kShortByOne = 0x8000u | (kCrcPoly >> 1);  // K(E_true - 1) classes differ by this

struct Packet {
  std::vector<uint8_t> bits;   // from the sync word (the preamble's last bit) to the end of the burst
  int64_t preamble_bits = 0;
  double dev_hz = 0, f_off_hz = 0;
  double f_carrier = 0;        // set by the caller (fc + slice offset + f_off_hz)
};

// python slice_packet(y, ofs) on y[0..n)
std::optional<Packet> slice_packet(const cf* y, size_t n, double ofs);

// python verify_links(frames, starts_us): one decoded dict per frame
std::vector<Json> verify_links(const std::vector<std::vector<uint8_t>>& frames, const std::vector<double>& starts_us);

// python frame_period_us(starts_us), nullopt for < 3 packets
std::optional<double> frame_period_us(std::vector<double> starts_us);

}  // namespace si::rx::rc_link::esb
