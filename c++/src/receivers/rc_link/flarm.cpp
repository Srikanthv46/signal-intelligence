// python receivers/rc_link/flarm.py (see there for the sources; the obfuscated payload is not decoded).
#include "receivers/rc_link/flarm.hpp"

#include <cstdio>
#include <tuple>
#include <vector>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::flarm {
namespace {
constexpr uint32_t kSync = 0xF531FAB6u;
}

bool in_band(double f) { return (865.5e6 <= f && f <= 866.5e6) || (867.9e6 <= f && f <= 868.7e6); }

Json try_identify(const cfvec& y, double ofs, double chip_rate, int max_err) {
  int pat[32];
  for (int i = 0; i < 32; ++i) pat[i] = ((kSync >> (31 - i)) & 1) ? 1 : -1;
  Json best = nullptr;   // python: the hit with the fewest (header, total) Manchester violations, then sync errors
  std::tuple<int, int, int> best_key{1 << 30, 0, 0};
  for (const auto& [ph, chips0] : sik::demod_bits(y, ofs, chip_rate))
    for (bool inv : {false, true}) {
      std::vector<uint8_t> ch(chips0.size());
      for (size_t i = 0; i < ch.size(); ++i) ch[i] = static_cast<uint8_t>(chips0[i] ^ (inv ? 1 : 0));
      for (int odd : {0, 1}) {
        if (ch.size() < static_cast<size_t>(odd)) continue;
        const size_t np_ = (ch.size() - odd) / 2;
        std::vector<uint8_t> bits(np_);
        for (size_t k = 0; k < np_; ++k) bits[k] = ch[odd + 2 * k + 1] > ch[odd + 2 * k] ? 1 : 0;
        if (bits.size() < 64) continue;
        for (size_t p = 0; p + 32 <= bits.size(); ++p) {
          int c = 0;
          for (int i = 0; i < 32; ++i) c += (bits[p + i] ? 1 : -1) * pat[i];
          if (c < 32 - 2 * max_err) continue;
          if (p + 64 > bits.size()) continue;
          int by[4];
          for (int b = 0; b < 4; ++b) {
            int v = 0;
            for (int i = 0; i < 8; ++i) v = (v << 1) | bits[p + 32 + 8 * b + i];
            by[b] = v;
          }
          int viol = 0;
          for (size_t k = 0; k < 64; ++k) {
            const size_t a = odd + 2 * (p + k);
            if (a + 1 < ch.size() && ch[a] == ch[a + 1]) ++viol;
          }
          if (viol > 4) continue;
          const int rid = by[0] | (by[1] << 8) | (by[2] << 16);
          int hv = 0;   // python: Manchester violations in the 32 header bits -> the ID is not reliable
          for (size_t k = 32; k < 64; ++k) {
            const size_t a = odd + 2 * (p + k);
            if (a + 1 < ch.size() && ch[a] == ch[a + 1]) ++hv;
          }
          char rs[8], ss[12], at[8];
          std::snprintf(rs, sizeof rs, "%06X", rid);
          std::snprintf(ss, sizeof ss, "0x%08X", static_cast<unsigned>(kSync));
          std::snprintf(at, sizeof at, "0x%02X", by[3]);
          Json out = Json::object();
          out["protocol"] = "FLARM-like (868 MHz beacon)";
          out["sync"] = ss;
          out["sync_bit_errors"] = (32 - c) / 2;
          out["manchester_violations"] = viol;
          out["radio_id"] = hv == 0 ? Json(rs) : Json(nullptr);
          out["address_type"] = hv != 0 ? Json(nullptr)
                                : by[3] == 0x20 ? Json("FLARM ID") : by[3] == 0x10 ? Json("ICAO 24-bit address") : Json(at);
          out["polarity_inverted"] = inv;
          out["payload"] = "obfuscated (not decoded)";
          if (hv) out["note"] = std::to_string(hv) + " Manchester violation(s) in the header: ID not reliable";
          const std::tuple<int, int, int> key{hv, viol, (32 - c) / 2};
          if (best.is_null() || key < best_key) best_key = key, best = std::move(out);
        }
      }
    }
  return best;
}

}  // namespace si::rx::rc_link::flarm
