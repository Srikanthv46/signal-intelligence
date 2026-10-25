// Python semantics the Wi-Fi MAC dissector relies on: a bytes view (slices clamp, indexing
// raises IndexError), bytes.decode("utf-8"/"ascii", "replace"), str.strip/split/lower and
// truthiness of dict values. Shared by the wifi receiver (its detail namespace re-exports these) and
// core/remote_id (Open Drone ID over Wi-Fi and Bluetooth).
#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::pybytes {

// python IndexError (caught where the Python code catches it)
struct IndexError : std::out_of_range {
  IndexError() : std::out_of_range("index out of range") {}
};

// python bytes: b[i] raises past the end, b[a:b] clamps.
class Bytes {
 public:
  Bytes() = default;
  Bytes(const uint8_t* p, size_t n) : p_(p), n_(n) {}
  explicit Bytes(const std::vector<uint8_t>& v) : p_(v.data()), n_(v.size()) {}
  size_t size() const { return n_; }
  bool empty() const { return n_ == 0; }
  const uint8_t* data() const { return p_; }
  int operator[](size_t i) const {
    if (i >= n_) throw IndexError();
    return p_[i];
  }
  Bytes sub(size_t a, size_t b) const {  // b[a:b]
    if (a > n_) a = n_;
    if (b > n_) b = n_;
    if (b < a) b = a;
    return Bytes(p_ + a, b - a);
  }
  Bytes from(size_t a) const { return sub(a, n_); }
  std::vector<uint8_t> vec() const { return std::vector<uint8_t>(p_, p_ + n_); }

 private:
  const uint8_t* p_ = nullptr;
  size_t n_ = 0;
};

int u16le(const Bytes& b, size_t i);            // b[i] | b[i+1] << 8
int be16(const Bytes& b, size_t i);             // b[i] << 8 | b[i+1]
uint64_t le_uint(const Bytes& b);               // int.from_bytes(b, "little") (<= 8 bytes)
uint64_t be_uint(const Bytes& b);               // int.from_bytes(b, "big") (<= 8 bytes)
std::string hex(const Bytes& b);                // b.hex()
std::string mac_str(const Bytes& b);            // ":".join(f"{x:02x}")

std::string decode_utf8(const Bytes& b);        // b.decode("utf-8", "replace")
std::string decode_ascii(const Bytes& b);       // b.decode("ascii", "replace")

std::string py_strip(const std::string& s);                          // str.strip()
std::string py_strip_chars(const std::string& s, const std::string& ascii_chars);  // str.strip(chars)
std::vector<std::string> py_split_ws(const std::string& s);          // str.split()
std::string ascii_lower(std::string s);
std::string ascii_upper(std::string s);
std::string head_codepoints(const std::string& s, size_t n);         // s[:n]
size_t len_codepoints(const std::string& s);

bool truthy(const Json& v);                                          // bool(v)
// d.get(k) (null if d is not a dict or k is absent)
const Json& jget(const Json& d, const char* k);
// str(v) as an f-string prints it (strings unquoted, None, True, floats as repr)
std::string pstr(const Json& v);
// a Python dict key as json.dumps writes it (str: itself, float: repr, int: str, None: "null")
std::string json_key(const Json& v);

}  // namespace si::pybytes
