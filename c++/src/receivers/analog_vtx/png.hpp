// Minimal grayscale PNG writer for the demodulated picture (python: cli._write_pictures,
// PIL Image.fromarray(uint8).save). zlib-free: the IDAT stream uses stored (uncompressed)
// deflate blocks, so the file is larger than PIL's and its bytes differ, but it decodes to
// the identical 8-bit pixels (PIL, libpng and browsers read stored deflate).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/model.hpp"

namespace si::rx::analog_vtx {

// 8-bit grayscale PNG of a rows x cols picture (pixels 0..255, row-major)
std::vector<uint8_t> png_gray8(const uint8_t* pixels, size_t rows, size_t cols);
std::vector<uint8_t> png_gray8(const Image& img);   // Image pixels hold 0..255
bool write_png(const std::string& path, const Image& img);

}  // namespace si::rx::analog_vtx
