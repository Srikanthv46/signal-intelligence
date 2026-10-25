#include "receivers/analog_vtx/png.hpp"

#include <algorithm>
#include <fstream>

namespace si::rx::analog_vtx {

namespace {

uint32_t crc32(const uint8_t* p, size_t n, uint32_t c = 0) {
  static uint32_t table[256];
  static bool init = [] {
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t v = i;
      for (int k = 0; k < 8; ++k) v = (v & 1) ? 0xEDB88320u ^ (v >> 1) : v >> 1;
      table[i] = v;
    }
    return true;
  }();
  (void)init;
  c = ~c;
  for (size_t i = 0; i < n; ++i) c = table[(c ^ p[i]) & 0xFF] ^ (c >> 8);
  return ~c;
}

void be32(std::vector<uint8_t>& o, uint32_t v) {
  for (int s = 24; s >= 0; s -= 8) o.push_back(static_cast<uint8_t>(v >> s));
}

void chunk(std::vector<uint8_t>& o, const char* type, const std::vector<uint8_t>& data) {
  be32(o, static_cast<uint32_t>(data.size()));
  const size_t start = o.size();
  o.insert(o.end(), type, type + 4);
  o.insert(o.end(), data.begin(), data.end());
  be32(o, crc32(o.data() + start, o.size() - start));
}

}  // namespace

std::vector<uint8_t> png_gray8(const uint8_t* px, size_t rows, size_t cols) {
  std::vector<uint8_t> raw;   // filter byte 0 + row
  raw.reserve(rows * (cols + 1));
  for (size_t r = 0; r < rows; ++r) {
    raw.push_back(0);
    raw.insert(raw.end(), px + r * cols, px + (r + 1) * cols);
  }
  std::vector<uint8_t> z = {0x78, 0x01};   // zlib header, no compression
  uint32_t a = 1, b = 0;
  for (uint8_t v : raw) {
    a = (a + v) % 65521;
    b = (b + a) % 65521;
  }
  size_t pos = 0;
  do {
    const size_t len = std::min<size_t>(65535, raw.size() - pos);
    const bool last = pos + len == raw.size();
    z.push_back(last ? 1 : 0);
    z.push_back(static_cast<uint8_t>(len & 0xFF));
    z.push_back(static_cast<uint8_t>(len >> 8));
    z.push_back(static_cast<uint8_t>(~len & 0xFF));
    z.push_back(static_cast<uint8_t>((~len >> 8) & 0xFF));
    z.insert(z.end(), raw.begin() + pos, raw.begin() + pos + len);
    pos += len;
  } while (pos < raw.size());
  be32(z, (b << 16) | a);
  std::vector<uint8_t> o = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
  std::vector<uint8_t> ihdr;
  be32(ihdr, static_cast<uint32_t>(cols));
  be32(ihdr, static_cast<uint32_t>(rows));
  ihdr.insert(ihdr.end(), {8, 0, 0, 0, 0});   // 8-bit grayscale, deflate, no filter set, no interlace
  chunk(o, "IHDR", ihdr);
  chunk(o, "IDAT", z);
  chunk(o, "IEND", {});
  return o;
}

std::vector<uint8_t> png_gray8(const Image& img) {
  std::vector<uint8_t> px(img.pixels.size());
  for (size_t i = 0; i < px.size(); ++i) px[i] = static_cast<uint8_t>(std::min(std::max(img.pixels[i], 0.f), 255.f));
  return png_gray8(px.data(), img.rows, img.cols);
}

bool write_png(const std::string& path, const Image& img) {
  const auto bytes = png_gray8(img);
  std::ofstream f(path, std::ios::binary);
  f.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  return f.good();
}

}  // namespace si::rx::analog_vtx
