// scipy.ndimage equivalents (scipy 1.13 semantics): boundary modes, origin rule for even
// sizes, double accumulation with the output in the input dtype, label numbering.
//
// Window of a 1-D filter of size k at sample i (origin o): x[i - k/2 - o ... i - k/2 - o + k-1]
// (k/2 integer division), with x extended past the ends by the mode:
//   reflect  (d c b a | a b c d | d c b a)   default
//   mirror   (d c b | a b c d | c b a)
//   nearest  (a a a | a b c d | d d d)
//   wrap     (a b c d | a b c d | a b c d)
//   constant (k k k | a b c d | k k k), k = cval
// grey_dilation / binary_dilation / convolve1d flip the footprint and origin as scipy does
// (even sizes then shift by one sample the other way).
#pragma once

#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::ndi {

enum class Mode { Reflect, Nearest, Mirror, Wrap, Constant };
Mode mode_from(const std::string& s);  // "reflect" | "nearest" | "mirror" | "wrap" | "constant"

// index of the extended sample i in [0, n), or -1 for a constant-mode pad sample
int64_t ext_index(int64_t i, int64_t n, Mode m);

// ------------------------------------------------------------------ uniform (box mean)
// scipy.ndimage.uniform_filter1d: running sum in double (tmp += (in - out)/size), output
// cast to T. T = float or double.
template <class T>
std::vector<T> uniform_filter1d(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect,
                                double cval = 0.0, int origin = 0);
template <class T>
Mat<T> uniform_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode = Mode::Reflect,
                        double cval = 0.0, int origin = 0);
// scipy.ndimage.uniform_filter(x, size=(s0, s1)): axis 0 then axis 1 (sizes of 1 skipped),
// the intermediate stored in T as scipy does.
template <class T>
Mat<T> uniform_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode = Mode::Reflect, double cval = 0.0);
template <class T>
std::vector<T> uniform_filter(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect, double cval = 0.0) {
  return uniform_filter1d(x, size, mode, cval, 0);
}

// ------------------------------------------------------------------ rank filters
// scipy.ndimage.median_filter: the element of rank size//2 of the window (for an even
// size that is the upper median, no averaging), exact for any T.
template <class T>
std::vector<T> median_filter(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect,
                             double cval = 0.0, int origin = 0);
// 2-D, size (s0, s1); rank (s0*s1)//2. (1, k) / (k, 1) run as 1-D lines.
template <class T>
Mat<T> median_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode = Mode::Reflect, double cval = 0.0);

// scipy.ndimage.minimum_filter1d / maximum_filter1d
template <class T>
std::vector<T> minimum_filter1d(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect,
                                double cval = 0.0, int origin = 0);
template <class T>
std::vector<T> maximum_filter1d(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect,
                                double cval = 0.0, int origin = 0);
template <class T>
Mat<T> minimum_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode = Mode::Reflect,
                        double cval = 0.0, int origin = 0);
template <class T>
Mat<T> maximum_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode = Mode::Reflect,
                        double cval = 0.0, int origin = 0);
// scipy.ndimage.minimum_filter / maximum_filter with size=(s0, s1) (separable)
template <class T>
Mat<T> minimum_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode = Mode::Reflect, double cval = 0.0);
template <class T>
Mat<T> maximum_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode = Mode::Reflect, double cval = 0.0);

// scipy.ndimage.grey_erosion / grey_dilation / grey_opening / grey_closing (size=k, 1-D)
template <class T>
std::vector<T> grey_erosion(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect, double cval = 0.0);
template <class T>
std::vector<T> grey_dilation(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect, double cval = 0.0);
template <class T>
std::vector<T> grey_opening(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect, double cval = 0.0);
template <class T>
std::vector<T> grey_closing(const std::vector<T>& x, size_t size, Mode mode = Mode::Reflect, double cval = 0.0);

// ------------------------------------------------------------------ correlate / convolve
// scipy.ndimage.correlate1d / convolve1d (double weights; scipy's symmetric /
// antisymmetric fast paths and their summation order). T = float, double, cf, cd
// (complex: real and imaginary parts filtered separately, as scipy).
template <class T>
std::vector<T> correlate1d(const std::vector<T>& x, const dvec& w, Mode mode = Mode::Reflect,
                           double cval = 0.0, int origin = 0);
template <class T>
std::vector<T> convolve1d(const std::vector<T>& x, const dvec& w, Mode mode = Mode::Reflect,
                          double cval = 0.0, int origin = 0);
template <class T>
Mat<T> correlate1d(const Mat<T>& x, const dvec& w, int axis, Mode mode = Mode::Reflect,
                   double cval = 0.0, int origin = 0);
template <class T>
Mat<T> convolve1d(const Mat<T>& x, const dvec& w, int axis, Mode mode = Mode::Reflect,
                  double cval = 0.0, int origin = 0);

// ------------------------------------------------------------------ labelling
// generate_binary_structure(2, 1) (cross) and np.ones((3, 3))
bmat cross3();
bmat ones(size_t r, size_t c);

struct Labels {
  imat lab;   // 0 = background, 1..n
  int n = 0;
};
// scipy.ndimage.label(mask, structure) for a 3x3 structure (default: cross). Labels are
// numbered in raster order of each component's first pixel, as scipy.
Labels label(const bmat& mask, const bmat& structure);
inline Labels label(const bmat& mask) { return label(mask, cross3()); }

// 1-D label (structure [1,1,1]): runs of nonzero samples
struct Labels1d {
  std::vector<int32_t> lab;
  int n = 0;
};
Labels1d label(const std::vector<uint8_t>& mask);

// scipy.ndimage.find_objects: bounding slices [r0, r1) x [c0, c1) of labels 1..n
struct Box2 {
  size_t r0 = 0, r1 = 0, c0 = 0, c1 = 0;
  bool empty() const { return r1 == r0; }
};
std::vector<Box2> find_objects(const imat& lab, int n);
std::vector<std::pair<size_t, size_t>> find_objects(const std::vector<int32_t>& lab, int n);

// ------------------------------------------------------------------ binary morphology
// scipy.ndimage.binary_erosion / dilation / opening / closing with a given structure,
// iterations=1, border_value=0, origin=0. Outside the array counts as 0 (so erosion, and
// therefore closing, clears pixels within the structure's reach of the border).
bmat binary_erosion(const bmat& x, const bmat& structure);
bmat binary_dilation(const bmat& x, const bmat& structure);
bmat binary_opening(const bmat& x, const bmat& structure);
bmat binary_closing(const bmat& x, const bmat& structure);
// 1-D with a structure of `size` ones (np.ones(size)), or an explicit 1-D structure
std::vector<uint8_t> binary_erosion(const std::vector<uint8_t>& x, const std::vector<uint8_t>& structure);
std::vector<uint8_t> binary_dilation(const std::vector<uint8_t>& x, const std::vector<uint8_t>& structure);
std::vector<uint8_t> binary_opening(const std::vector<uint8_t>& x, const std::vector<uint8_t>& structure);
std::vector<uint8_t> binary_closing(const std::vector<uint8_t>& x, const std::vector<uint8_t>& structure);
inline std::vector<uint8_t> ones1(size_t k) { return std::vector<uint8_t>(k, 1); }

}  // namespace si::ndi
