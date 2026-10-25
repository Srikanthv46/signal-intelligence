// python core/rs.py (see there).
#include "core/rs.hpp"

#include <algorithm>

namespace si::rs {

RS::RS(unsigned gfpoly, int fcr, int prim, int nroots) : nroots_(nroots), fcr_(fcr), prim_(prim) {
  unsigned x = 1;
  for (int i = 0; i < 255; ++i) {
    exp_[i] = static_cast<int>(x);
    log_[x] = i;
    x <<= 1;
    if (x & 0x100) x ^= gfpoly;
  }
  log_[0] = 0;
  for (int i = 255; i < 512; ++i) exp_[i] = exp_[i - 255];
  gen_ = {1};
  for (int i = 0; i < nroots; ++i) {
    const int r = exp_[(prim * (fcr + i)) % 255];
    std::vector<int> ng(gen_.size() + 1, 0);
    for (size_t j = 0; j < gen_.size(); ++j) {
      ng[j] ^= gen_[j];
      ng[j + 1] ^= mul(gen_[j], r);
    }
    gen_ = ng;
  }
}

std::vector<uint8_t> RS::encode(const std::vector<uint8_t>& data) const {
  std::vector<int> par(static_cast<size_t>(nroots_), 0);
  for (uint8_t d : data) {
    const int fb = d ^ par[0];
    par.erase(par.begin());
    par.push_back(0);
    if (fb)
      for (int j = 0; j < nroots_; ++j) par[static_cast<size_t>(j)] ^= mul(fb, gen_[static_cast<size_t>(j) + 1]);
  }
  return std::vector<uint8_t>(par.begin(), par.end());
}

std::optional<std::pair<std::vector<uint8_t>, int>> RS::decode(const std::vector<uint8_t>& block) const {
  const int n = static_cast<int>(block.size()), nr = nroots_;
  std::vector<uint8_t> r = block;
  std::vector<int> S(static_cast<size_t>(nr));
  bool any = false;
  for (int i = 0; i < nr; ++i) {
    const int a = exp_[(prim_ * (fcr_ + i)) % 255];
    int s = 0;
    for (uint8_t v : r) s = mul(s, a) ^ v;
    S[static_cast<size_t>(i)] = s;
    any = any || s;
  }
  if (!any) return std::make_pair(r, 0);
  std::vector<int> lam(static_cast<size_t>(nr) + 1, 0), b(static_cast<size_t>(nr) + 1, 0);
  lam[0] = b[0] = 1;
  int L = 0, m = 1, bb = 1;
  for (int k = 0; k < nr; ++k) {
    int d = S[static_cast<size_t>(k)];
    for (int i = 1; i <= L; ++i) d ^= mul(lam[static_cast<size_t>(i)], S[static_cast<size_t>(k - i)]);
    if (d == 0) {
      ++m;
      continue;
    }
    const int coef = mul(d, exp_[255 - log_[bb]]);
    const std::vector<int> t = lam;
    for (int i = m; i <= nr; ++i) lam[static_cast<size_t>(i)] ^= mul(coef, b[static_cast<size_t>(i - m)]);
    if (2 * L <= k) {
      L = k + 1 - L;
      b = t;
      bb = d;
      m = 1;
    } else {
      ++m;
    }
  }
  int deg = 0;
  for (int i = 0; i <= nr; ++i)
    if (lam[static_cast<size_t>(i)]) deg = i;
  if (deg != L || L > nr / 2) return std::nullopt;
  std::vector<int> pos, Xs;
  for (int p = 0; p < n; ++p) {
    const int e = (prim_ * (n - 1 - p)) % 255;
    const int xinv = exp_[(255 - e) % 255];
    int v = 0, xp = 1;
    for (int i = 0; i <= L; ++i) {
      v ^= mul(lam[static_cast<size_t>(i)], xp);
      xp = mul(xp, xinv);
    }
    if (v == 0) {
      pos.push_back(p);
      Xs.push_back(exp_[e]);
    }
  }
  if (static_cast<int>(pos.size()) != L) return std::nullopt;
  std::vector<int> om(static_cast<size_t>(nr), 0);
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j <= std::min(i, L); ++j) om[static_cast<size_t>(i)] ^= mul(S[static_cast<size_t>(i - j)], lam[static_cast<size_t>(j)]);
  for (size_t q = 0; q < pos.size(); ++q) {
    const int X = Xs[q];
    const int xi = exp_[255 - log_[X]];
    int num = 0, xp = 1;
    for (int c : om) {
      num ^= mul(c, xp);
      xp = mul(xp, xi);
    }
    int den = 0;
    for (int i = 1; i <= L; i += 2) den ^= mul(lam[static_cast<size_t>(i)], i > 1 ? exp_[(log_[xi] * (i - 1)) % 255] : 1);
    if (den == 0) return std::nullopt;
    const int e1 = ((log_[X] * (1 - fcr_)) % 255 + 255) % 255;
    const int sc = exp_[e1];
    r[static_cast<size_t>(pos[q])] ^= static_cast<uint8_t>(mul(mul(num, sc), exp_[255 - log_[den]]));
  }
  return std::make_pair(r, static_cast<int>(pos.size()));
}

}  // namespace si::rs
