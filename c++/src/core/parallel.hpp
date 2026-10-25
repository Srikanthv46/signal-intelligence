// Thread-level parallelism for independent work items (receivers, carrier hypotheses,
// raster channels). Results come back in input order, so a parallel run gives exactly
// what the sequential run gives. Mirrors python/spectrum_intel/core/parallel.py:
// SPECTRUM_INTEL_THREADS sets the count (1 = sequential), default min(4, cores).
#pragma once

#include <atomic>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace si {

int threads();                  // SPECTRUM_INTEL_THREADS or min(4, hardware threads)
void set_threads(int n);        // override for this process (0 = back to the default)

// for i in [0, n): fn(i), on up to nthreads threads (0 = threads()); exceptions are
// rethrown (the first one) after every started item has finished.
void parallel_for(size_t n, const std::function<void(size_t)>& fn, int nthreads = 0);

// [fn(items[0]), fn(items[1]), ...] in input order.
template <class In, class Fn>
auto pmap(const std::vector<In>& items, Fn fn, int nthreads = 0) -> std::vector<decltype(fn(items[0]))> {
  using Out = decltype(fn(items[0]));
  std::vector<Out> out(items.size());
  parallel_for(items.size(), [&](size_t i) { out[i] = fn(items[i]); }, nthreads);
  return out;
}

}  // namespace si
