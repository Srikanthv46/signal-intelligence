// Python-dict semantics over Json values, for the device tracker (python: plain dict / Counter
// code in receivers/devices). Header-only; used by tracker.cpp, threats.cpp and devices.cpp.
#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "core/pyfmt.hpp"
#include "core/types.hpp"

namespace si::rx::devices::detail {

// bool(x)
inline bool truthy(const Json& v) {
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number_integer() || v.is_number_unsigned()) return v.get<int64_t>() != 0;
  if (v.is_number_float()) return v.get<double>() != 0.0;
  if (v.is_string()) return !v.get_ref<const std::string&>().empty();
  return !v.empty();
}

// d.get(k) (None if d is not a dict or k is absent)
inline Json get(const Json& d, const char* k) {
  if (!d.is_object()) return nullptr;
  auto it = d.find(k);
  return it == d.end() ? Json(nullptr) : *it;
}

// a or b
inline Json py_or(const Json& a, const Json& b) { return truthy(a) ? a : b; }

inline bool is_num(const Json& v) { return v.is_number() || v.is_boolean(); }
inline double num(const Json& v) { return v.is_boolean() ? (v.get<bool>() ? 1.0 : 0.0) : v.get<double>(); }

// a == b as Python compares dict keys / set members (1 == 1.0 == True)
inline bool py_eq(const Json& a, const Json& b) {
  if (is_num(a) && is_num(b)) return num(a) == num(b);
  return a == b;
}

// the string json.dumps writes for a dict key
inline std::string key_str(const Json& k) {
  if (k.is_string()) return k.get<std::string>();
  if (k.is_boolean()) return k.get<bool>() ? "true" : "false";
  if (k.is_null()) return "null";
  if (k.is_number_float()) return py_float_repr(k.get<double>());
  return k.dump();
}

inline bool str_startswith(const Json& v, const std::string& p) {
  return v.is_string() && v.get_ref<const std::string&>().compare(0, p.size(), p) == 0;
}

// collections.Counter with insertion order (most_common() is a stable sort on the counts)
class Counter {
 public:
  void add(const Json& k, int64_t n = 1) {
    for (auto& e : v_)
      if (py_eq(e.first, k)) {
        e.second += n;
        return;
      }
    v_.emplace_back(k, n);
  }
  bool empty() const { return v_.empty(); }
  const std::vector<std::pair<Json, int64_t>>& items() const { return v_; }
  std::vector<std::pair<Json, int64_t>> most_common(size_t n = SIZE_MAX) const {
    auto s = v_;
    std::stable_sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    if (s.size() > n) s.resize(n);
    return s;
  }
  Json to_json() const {  // dict(counter) after a json round trip
    Json o = Json::object();
    for (const auto& e : v_) o[key_str(e.first)] = e.second;
    return o;
  }

 private:
  std::vector<std::pair<Json, int64_t>> v_;
};

// sorted(set(strings))
inline std::vector<std::string> sorted_unique(std::vector<std::string> v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

inline Json str_array(const std::vector<std::string>& v) {
  Json a = Json::array();
  for (const auto& s : v) a.push_back(s);
  return a;
}

}  // namespace si::rx::devices::detail
