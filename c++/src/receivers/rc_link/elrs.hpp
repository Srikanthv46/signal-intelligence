// ExpressLRS over-the-air packets: OTA4 / OTA8, self-verifying SYNC (python receivers/rc_link/elrs.py; research,
// synthetic validation only, PLAN R2c).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::elrs {

int crc(const std::vector<uint8_t>& data, int init, int bits);
// python parse(payload): the dict or null
Json parse(const std::vector<uint8_t>& payload);

// LR1121 GFSK modes (research, PLAN G3): python gfsk_try_decode(y, ofs) -> verified SYNC packet dict or null
constexpr double kGfskRate = 300e3, kGfskDev = 100e3;
std::vector<uint8_t> gfsk_whiten(const std::vector<uint8_t>& b, const char* variant);
Json gfsk_try_decode(const cfvec& y, double ofs, double rate = kGfskRate);

// FHSS hop sequence (research, PLAN R2d; python fhss_seed / fhss_sequence / uid23_candidates / recover_hops)
uint32_t fhss_seed(int uid2, int uid3, int uid4, int uid5, int ver = 4);
std::vector<int> fhss_sequence(uint32_t seed, int count);
std::vector<std::pair<int, int>> uid23_candidates(int uid4, int uid5, int count,
                                                  const std::vector<std::pair<int64_t, int>>& observations, int ver = 4,
                                                  int uid3 = -1);
// packets: (t_us, channel index); the dict or null
Json recover_hops(const Json& sync, double t_sync_us, int sync_ch, double interval_us, int hop_every, int count,
                  const std::vector<std::pair<double, int>>& packets);

}  // namespace si::rx::rc_link::elrs
