// Toy-drone protocols on nRF24L01 / XN297 radios (python receivers/rc_link/toys.py; research, synthetic only, PLAN G6).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::toys {

// python match_xn297(addr, unscrambled): XN297 address (register order) + de-scrambled payload -> dict or null
Json match_xn297(const std::vector<uint8_t>& addr, const std::vector<uint8_t>& unscrambled);
// python crc16_nrf(b"", bits=bits): nRF24 CRC-16 over a bit list
uint16_t crc16_nrf_bits(const std::vector<uint8_t>& bits, unsigned init = 0xFFFF);
// python nrf24_try_decode(y, ofs, rate=RATE): CRC-verified plain nRF24 toy protocol or null
Json nrf24_try_decode(const cfvec& y, double ofs, double rate = 1e6);

}  // namespace si::rx::rc_link::toys
