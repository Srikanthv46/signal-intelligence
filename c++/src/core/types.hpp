// Basic types shared by every module: complex samples, a row-major 2-D array, JSON,
// and the Python-compatible rounding the report numbers go through.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace si {

using cf = std::complex<float>;     // numpy complex64
using cd = std::complex<double>;    // numpy complex128
using Json = nlohmann::ordered_json;  // keeps insertion order, like a Python dict

using fvec = std::vector<float>;
using dvec = std::vector<double>;
using cfvec = std::vector<cf>;
using cdvec = std::vector<cd>;
using ivec = std::vector<int64_t>;

constexpr double kPi = 3.14159265358979323846;

// Row-major 2-D array (numpy (rows, cols) C order).
template <class T>
struct Mat {
  size_t rows = 0, cols = 0;
  std::vector<T> data;

  Mat() = default;
  Mat(size_t r, size_t c, T v = T()) : rows(r), cols(c), data(r * c, v) {}
  T& operator()(size_t r, size_t c) { return data[r * cols + c]; }
  const T& operator()(size_t r, size_t c) const { return data[r * cols + c]; }
  T* row(size_t r) { return data.data() + r * cols; }
  const T* row(size_t r) const { return data.data() + r * cols; }
  size_t size() const { return data.size(); }
  bool empty() const { return data.empty(); }
};

using fmat = Mat<float>;
using dmat = Mat<double>;
using cfmat = Mat<cf>;
using bmat = Mat<uint8_t>;  // boolean mask
using imat = Mat<int32_t>;  // labels

// Python round(x, ndigits): the decimal closest to the exact binary value, ties to even
// (printf with the C library's correctly rounded conversion, then back to double).
double py_round(double x, int ndigits);
// Python round(x) with no ndigits: nearest integer, ties to even.
int64_t py_round_int(double x);

// A float that the report shows as Python would after round(x, n); NaN/Inf become null.
Json jnum(double x, int ndigits);
inline Json jnum(double x) { return std::isfinite(x) ? Json(x) : Json(nullptr); }

}  // namespace si
