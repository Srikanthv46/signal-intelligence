// python receivers/rc_link/mlrs.py (see there for the sources and the conventions tried).
#include "receivers/rc_link/mlrs.hpp"

#include <cstdio>
#include <string>

namespace si::rx::rc_link::mlrs {

const std::vector<Mode>& modes() {
  static const std::vector<Mode> t = {
      {"mlrs_50hz_2g4", "2G4", 5, 812.5e3, 5, 12, kFrameLen, 20000, false},
      {"mlrs_31hz_2g4", "2G4", 6, 812.5e3, 5, 12, kFrameLen, 32000, false},
      {"mlrs_19hz_2g4", "2G4", 7, 812.5e3, 5, 12, kFrameLen, 53000, false},
      {"mlrs_31hz_900", "900", 5, 500e3, 5, 12, kFrameLen, 32000, false},
      {"mlrs_19hz_900", "900", 6, 500e3, 5, 12, kFrameLen, 53000, false},
      {"mlrs_19hz_7x", "900", 6, 500e3, 5, 12, kFrameLen, 53000, true}};
  return t;
}

uint16_t crc16_mcrf4xx(const uint8_t* p, size_t n, unsigned init) {
  unsigned c = init;
  for (size_t k = 0; k < n; ++k) {
    c ^= p[k];
    for (int i = 0; i < 8; ++i) c = (c & 1) ? (c >> 1) ^ 0x8408u : c >> 1;
  }
  return static_cast<uint16_t>(c);
}

Json parse(const std::vector<uint8_t>& f) {
  if (static_cast<int>(f.size()) != kFrameLen) return nullptr;
  std::vector<uint8_t> all(f.begin(), f.begin() + 89), skip;
  for (size_t i = 0; i < 89; ++i)
    if (i != 13 && i != 14) skip.push_back(f[i]);   // python _CRC1_POS = (13, 14)
  const std::pair<const char*, const std::vector<uint8_t>*> vars[] = {{"all bytes 0-88", &all}, {"crc1 bytes skipped", &skip}};
  for (const auto& [cov, data] : vars) {
    const unsigned c = crc16_mcrf4xx(data->data(), data->size());
    for (const char* order : {"little", "big"}) {
      const bool le = std::string(order) == "little";
      const unsigned rx = le ? (f[89] | (f[90] << 8)) : ((f[89] << 8) | f[90]);
      if (c != rx) continue;
      const unsigned sw = le ? (f[0] | (f[1] << 8)) : ((f[0] << 8) | f[1]);
      char sws[8], st[12];
      std::snprintf(sws, sizeof sws, "0x%04X", sw);
      std::snprintf(st, sizeof st, "%02x%02x%02x%02x%02x", f[2], f[3], f[4], f[5], f[6]);
      Json out = Json::object();
      out["protocol"] = "mLRS";
      out["crc_ok"] = true;
      out["frame_sync_word"] = sws;
      out["frame_bytes"] = kFrameLen;
      out["convention"] = {{"crc_coverage", cov}, {"byte_order", order}};
      out["status_hex"] = st;
      return out;
    }
  }
  return nullptr;
}

}  // namespace si::rx::rc_link::mlrs
