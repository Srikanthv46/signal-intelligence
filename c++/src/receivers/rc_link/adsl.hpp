// ADS-L 4 SRD-860 M-band and the OGN tracker protocol on the FLARM PHY (python receivers/rc_link/adsl.py and ogn.py;
// research, synthetic only, PLAN G10).
#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::adsl {

// python crc24(data): Mode-S CRC-24 (0xFFF409, init 0)
uint32_t crc24(const uint8_t* p, size_t n);
// python xxtea_encrypt / xxtea_decrypt (zero key, 6 rounds) over little-endian words
std::vector<uint8_t> xxtea_encrypt(std::vector<uint8_t> d, int rounds = 6);
std::vector<uint8_t> xxtea_decrypt(std::vector<uint8_t> d, int rounds = 6);
// python parse(pkt): length + header + data + CRC-24 -> dict or null
Json parse(const std::vector<uint8_t>& pkt);
// python try_decode(y, ofs): the dict or null (CRC-24 verified)
Json try_decode(const cfvec& y, double ofs, double chip_rate = 100e3);

}  // namespace si::rx::rc_link::adsl

namespace si::rx::rc_link::ogn {

// python tea_key0(v0, v1, loops=8, decrypt)
std::pair<uint32_t, uint32_t> tea_key0(uint32_t v0, uint32_t v1, int loops = 8, bool decrypt = false);
// python parse(data20): dict or null
Json parse(const uint8_t* d20);
// python try_decode(y, ofs): the dict or null
Json try_decode(const cfvec& y, double ofs, double chip_rate = 100e3);

}  // namespace si::rx::rc_link::ogn
