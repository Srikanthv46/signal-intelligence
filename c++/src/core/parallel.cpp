#include "core/parallel.hpp"

#include <algorithm>
#include <cstdlib>
#include <string>

namespace si {
namespace {
std::atomic<int> g_override{0};
}

int threads() {
  if (int o = g_override.load()) return o;
  if (const char* v = std::getenv("SPECTRUM_INTEL_THREADS")) {
    int n = std::atoi(v);
    if (n >= 1) return n;
  }
  unsigned hw = std::thread::hardware_concurrency();
  return static_cast<int>(std::max(1u, std::min(4u, hw ? hw : 1u)));
}

void set_threads(int n) { g_override.store(std::max(0, n)); }

void parallel_for(size_t n, const std::function<void(size_t)>& fn, int nthreads) {
  if (n == 0) return;
  size_t nt = static_cast<size_t>(nthreads > 0 ? nthreads : threads());
  nt = std::min(nt, n);
  if (nt <= 1) {
    for (size_t i = 0; i < n; ++i) fn(i);
    return;
  }
  std::atomic<size_t> next{0};
  std::exception_ptr err;
  std::mutex err_mu;
  auto work = [&] {
    for (;;) {
      size_t i = next.fetch_add(1);
      if (i >= n) return;
      try {
        fn(i);
      } catch (...) {
        std::lock_guard<std::mutex> lk(err_mu);
        if (!err) err = std::current_exception();
      }
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nt - 1);
  for (size_t t = 0; t + 1 < nt; ++t) pool.emplace_back(work);
  work();
  for (auto& th : pool) th.join();
  if (err) std::rethrow_exception(err);
}

}  // namespace si
