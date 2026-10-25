// RDS codec (research, signal gaps S20; python: receivers/broadcast_fm/rds.py).
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "core/types.hpp"

namespace si::rx::broadcast_fm::rds {

int syndrome(uint32_t word26);
uint32_t encode_block(uint16_t data, const std::string& offset);   // offset A / B / C / Cp / D
std::vector<int> encode_group(uint16_t a, uint16_t b, uint16_t c, uint16_t d, bool version_b = false);

struct Block {
  int64_t pos;
  std::string offset;
  int data;
};
std::vector<Block> blocks(const std::vector<int>& bits);

struct Group {
  int64_t pos;
  int a;
  std::optional<int> b, c, d;
  bool version_b = false;
};
std::vector<Group> groups(const std::vector<Block>& blks);
Json decode(const std::vector<Group>& grps);

}  // namespace si::rx::broadcast_fm::rds
