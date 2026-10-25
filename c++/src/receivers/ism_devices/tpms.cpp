// TPMS frames (research, signal gaps S6; python: receivers/ism_devices/tpms.py): Continental / VDO (Manchester,
// CRC-8 0x07 / 0xAA) and Toyota PMV-C210 (differential Manchester, CRC-8 0x07 / 0x80), bit conventions of rtl_433.
#include "receivers/ism_devices/tpms.hpp"

#include <algorithm>

#include "core/pyfmt.hpp"

namespace si::rx::ism_devices::tpms {

namespace {

const std::vector<int> kVdoPre = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
                                  1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1};
const std::vector<int> kToyotaPre = {1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0};

long find(const std::vector<int>& s, const std::vector<int>& pat, size_t start) {
  if (s.size() < pat.size()) return -1;
  for (size_t i = start; i + pat.size() <= s.size(); ++i)
    if (std::equal(pat.begin(), pat.end(), s.begin() + static_cast<ptrdiff_t>(i))) return static_cast<long>(i);
  return -1;
}

std::vector<uint8_t> to_bytes(const std::vector<int>& bits, size_t n) {
  std::vector<uint8_t> b;
  for (size_t k = 0; k + 7 < n; k += 8) {
    int v = 0;
    for (int q = 0; q < 8; ++q) v |= bits[k + q] << (7 - q);
    b.push_back(static_cast<uint8_t>(v));
  }
  return b;
}

std::string hex(const std::vector<uint8_t>& b) {
  std::string s;
  for (uint8_t v : b) s += fmt("%02x", v);
  return s;
}

// python diff_manchester(s, start, max_bits) -> (bits, symbols consumed)
std::pair<std::vector<int>, size_t> diff_manchester(const std::vector<int>& s, size_t start, size_t max_bits) {
  std::vector<int> out;
  size_t i = start;
  const size_t end = std::min(s.size(), start + 2 * max_bits);
  int b2 = 0;
  while (i + 2 < end) {
    const int b1 = s[i], c2 = s[i + 1], b3 = s[i + 2];
    b2 = c2;
    i += 2;
    if (b1 != b2) {
      if (b2 != b3) {
        out.push_back(0);
      } else {
        b2 = b1;
        i -= 1;
        break;
      }
    } else {
      b2 = 1 - b1;
      i -= 2;
      break;
    }
  }
  while (i + 1 < end) {
    const int b1 = s[i];
    i += 1;
    if (b1 == b2) break;
    b2 = s[i];
    i += 1;
    out.push_back(b1 == b2 ? 1 : 0);
  }
  return {out, i - start};
}

}  // namespace

uint8_t crc8(const std::vector<uint8_t>& data, uint8_t poly, uint8_t init) {
  unsigned c = init;
  for (uint8_t b : data) {
    c ^= b;
    for (int k = 0; k < 8; ++k) c = (c & 0x80) ? ((c << 1) ^ poly) & 0xFF : (c << 1) & 0xFF;
  }
  return static_cast<uint8_t>(c);
}

std::vector<int> manchester(const std::vector<int>& s, size_t start, size_t max_bits) {
  std::vector<int> out;
  size_t i = start;
  const size_t end = std::min(s.size(), start + 2 * max_bits);
  while (i + 1 < end) {
    if (s[i] == s[i + 1]) break;
    out.push_back(s[i + 1]);
    i += 2;
  }
  return out;
}

std::optional<Json> vdo(const std::vector<int>& s) {
  std::vector<int> inv(s.size());
  for (size_t i = 0; i < s.size(); ++i) inv[i] = 1 - s[i];
  size_t pos = 0;
  while (true) {
    const long p = find(inv, kVdoPre, pos);
    if (p < 0 || static_cast<size_t>(p) + 32 + 160 > inv.size()) return std::nullopt;
    pos = static_cast<size_t>(p);
    const auto bits = manchester(inv, pos + 32, 80);
    if (bits.size() >= 80) {
      const auto b = to_bytes(bits, 80);
      if (crc8(std::vector<uint8_t>(b.begin(), b.begin() + 9), 0x07, 0xAA) == b[9]) {
        Json r = Json::object();
        r["protocol"] = "Continental / VDO";
        r["id"] = fmt("%08x", (unsigned(b[1]) << 24) | (unsigned(b[2]) << 16) | (unsigned(b[3]) << 8) | b[4]);
        r["state"] = static_cast<int64_t>(b[0]);
        r["flags"] = static_cast<int64_t>(b[5] >> 4);
        r["repeat"] = static_cast<int64_t>(b[5] & 0x0F);
        r["pressure_kpa"] = py_round(b[6] * 1.375, 3);
        r["temperature_c"] = static_cast<int64_t>(b[7]) - 50;
        r["battery_raw"] = static_cast<int64_t>(b[8]);
        r["crc_ok"] = true;
        r["frame_hex"] = hex(b);
        return r;
      }
    }
    pos += 2;
  }
}

std::optional<Json> toyota(const std::vector<int>& s) {
  size_t pos = 0;
  while (true) {
    const long p = find(s, kToyotaPre, pos);
    if (p < 0 || static_cast<size_t>(p) + 156 > s.size()) return std::nullopt;
    pos = static_cast<size_t>(p);
    const auto [bits, used] = diff_manchester(s, pos + 11, 80);
    if (used >= 144 && bits.size() >= 72) {
      const auto b = to_bytes(bits, 72);
      const int p1 = ((b[4] & 0x7F) << 1) | (b[5] >> 7);
      if (crc8(std::vector<uint8_t>(b.begin(), b.begin() + 8), 0x07, 0x80) == b[8] && p1 == (b[7] ^ 0xFF)) {
        const int temp = ((b[5] & 0x7F) << 1) | (b[6] >> 7);
        Json r = Json::object();
        r["protocol"] = "Toyota / Pacific PMV-C210";
        r["id"] = fmt("%08x", (unsigned(b[0]) << 24) | (unsigned(b[1]) << 16) | (unsigned(b[2]) << 8) | b[3]);
        r["status"] = static_cast<int64_t>((b[4] & 0x80) | (b[6] & 0x7F));
        r["pressure_psi"] = p1 * 0.25 - 7.0;
        r["pressure_kpa"] = py_round((p1 * 0.25 - 7.0) * 6.894757, 3);
        r["temperature_c"] = static_cast<int64_t>(temp) - 40;
        r["crc_ok"] = true;
        r["frame_hex"] = hex(b);
        return r;
      }
    }
    pos += 1;
  }
}

std::optional<Json> decode(const std::vector<int>& s) {
  for (int pol = 0; pol < 2; ++pol) {
    std::vector<int> t = s;
    if (pol)
      for (int& v : t) v = 1 - v;
    if (auto r = vdo(t)) return r;
    if (auto r = toyota(t)) return r;
  }
  return std::nullopt;
}

}  // namespace si::rx::ism_devices::tpms
