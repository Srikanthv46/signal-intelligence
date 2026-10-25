// python receivers/rc_link/fsk_frame.py (see there).
#include "receivers/rc_link/fsk_frame.hpp"

#include <algorithm>
#include <cstdio>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::fsk_frame {
namespace {
struct Known { const char* name; uint32_t word; int n; };
// python KNOWN_SYNCS (longest first)
constexpr Known kKnown[] = {{"FrSky D16-like (CC2500)", 0xD391D391u, 32},
                            {"FlySky AFHDS 2A-like (A7105 ID code)", 0x5475C52Au, 32},
                            {"TBS Crossfire-like", 0x3E6B3Eu, 24},
                            {"FrSky D8-like (CC2500)", 0xD391u, 16},
                            {"SiK-like (Si1000 / Si4432)", 0x2DD4u, 16}};

Json known_sync(const std::vector<uint8_t>& bits, int64_t e) {
  const int64_t nb = static_cast<int64_t>(bits.size());
  for (const auto& k : kKnown)
    for (int d : {0, -1, 1, -2, 2, -3, 3, -4, 4}) {
      const int64_t p = e + d;
      if (p < 0 || p + k.n > nb) continue;
      for (bool inv : {false, true}) {
        int err = 0;
        for (int i = 0; i < k.n; ++i)
          err += ((bits[static_cast<size_t>(p + i)] ^ (inv ? 1 : 0)) & 1) != ((k.word >> (k.n - 1 - i)) & 1);
        if (err > (k.n >= 24 ? 1 : 0)) continue;
        char hx[16];
        std::snprintf(hx, sizeof hx, "%0*X", k.n / 4, static_cast<unsigned>(k.word));
        return {{"name", k.name}, {"sync_hex", hx}, {"offset_bits", d}, {"inverted", inv}, {"bit_errors", err}};
      }
    }
  return nullptr;
}
}  // namespace

Json frame(const cfvec& y, double ofs, double rate) {
  int64_t best_ln = -1, best_s = 0;
  std::vector<uint8_t> best_bits;
  for (auto& [ph, bits] : sik::demod_bits(y, ofs, rate)) {
    const int64_t n = static_cast<int64_t>(bits.size());
    if (n < kMinPre + kSyncBits) continue;
    // python _alt_runs: maximal runs where consecutive bits differ; the first long one is the preamble
    int64_t i = 0;
    const int64_t nd = n - 1;
    while (i < nd) {
      if (bits[static_cast<size_t>(i)] != bits[static_cast<size_t>(i + 1)]) {
        int64_t j = i;
        while (j < nd && bits[static_cast<size_t>(j)] != bits[static_cast<size_t>(j + 1)]) ++j;
        const int64_t s = i, ln = j - i + 1;
        i = j;
        if (ln < kMinPre || s > std::max<int64_t>(16, n / 3)) continue;
        if (best_ln < 0 || ln > best_ln) best_ln = ln, best_s = s, best_bits = bits;
        break;
      }
      ++i;
    }
  }
  if (best_ln < 0) return nullptr;
  const int64_t e = best_s + best_ln;
  if (e + kSyncBits > static_cast<int64_t>(best_bits.size())) return nullptr;
  uint32_t w = 0;
  for (int k = 0; k < kSyncBits; ++k) w = (w << 1) | best_bits[static_cast<size_t>(e + k)];
  char hx[16];
  std::snprintf(hx, sizeof hx, "%08X", static_cast<unsigned>(w));
  Json out = Json::object();
  out["preamble_bits"] = best_ln;
  out["preamble_last_bit"] = static_cast<int>(best_bits[static_cast<size_t>(e - 1)]);
  out["sync32_hex"] = hx;
  out["payload_bits"] = static_cast<int64_t>(best_bits.size()) - e - kSyncBits;
  out["polarity"] = "bit 1 = higher frequency";
  Json ks = known_sync(best_bits, e);
  if (!ks.is_null()) out["known_sync"] = std::move(ks);
  return out;
}

}  // namespace si::rx::rc_link::fsk_frame
