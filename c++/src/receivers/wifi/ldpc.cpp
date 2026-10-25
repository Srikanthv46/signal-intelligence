// Port of python receivers/wifi/ldpc.py (see ldpc.hpp). The decoder reproduces the Python
// arithmetic exactly: the variable-node sums add the check messages in edge order (scipy's
// csc_matvecs over the check-major edge list), messages are alpha * sign * min.
#include "receivers/wifi/ldpc.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace si::rx::wifi::ldpc {

namespace {

// (648, 1/2)
const int8_t kB_648_12[12][24] = {
    {0,-1,-1,-1,0,0,-1,-1,0,-1,-1,0,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {22,0,-1,-1,17,-1,0,0,12,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {6,-1,0,-1,10,-1,-1,-1,24,-1,0,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1},
    {2,-1,-1,0,20,-1,-1,-1,25,0,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1},
    {23,-1,-1,-1,3,-1,-1,-1,0,-1,9,11,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1},
    {24,-1,23,1,17,-1,3,-1,10,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1},
    {25,-1,-1,-1,8,-1,-1,-1,7,18,-1,-1,0,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1},
    {13,24,-1,-1,0,-1,8,-1,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1},
    {7,20,-1,16,22,10,-1,-1,23,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1},
    {11,-1,-1,-1,19,-1,-1,-1,13,-1,3,17,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1},
    {25,-1,8,-1,23,18,-1,14,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0},
    {3,-1,-1,-1,16,-1,-1,2,25,5,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0}};
// (648, 2/3)
const int8_t kB_648_23[8][24] = {
    {25,26,14,-1,20,-1,2,-1,4,-1,-1,8,-1,16,-1,18,1,0,-1,-1,-1,-1,-1,-1},
    {10,9,15,11,-1,0,-1,1,-1,-1,18,-1,8,-1,10,-1,-1,0,0,-1,-1,-1,-1,-1},
    {16,2,20,26,21,-1,6,-1,1,26,-1,7,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1},
    {10,13,5,0,-1,3,-1,7,-1,-1,26,-1,-1,13,-1,16,-1,-1,-1,0,0,-1,-1,-1},
    {23,14,24,-1,12,-1,19,-1,17,-1,-1,-1,20,-1,21,-1,0,-1,-1,-1,0,0,-1,-1},
    {6,22,9,20,-1,25,-1,17,-1,8,-1,14,-1,18,-1,-1,-1,-1,-1,-1,-1,0,0,-1},
    {14,23,21,11,20,-1,24,-1,18,-1,19,-1,-1,-1,-1,22,-1,-1,-1,-1,-1,-1,0,0},
    {17,11,11,20,-1,21,-1,26,-1,3,-1,-1,18,-1,26,-1,1,-1,-1,-1,-1,-1,-1,0}};
// (648, 3/4)
const int8_t kB_648_34[6][24] = {
    {16,17,22,24,9,3,14,-1,4,2,7,-1,26,-1,2,-1,21,-1,1,0,-1,-1,-1,-1},
    {25,12,12,3,3,26,6,21,-1,15,22,-1,15,-1,4,-1,-1,16,-1,0,0,-1,-1,-1},
    {25,18,26,16,22,23,9,-1,0,-1,4,-1,4,-1,8,23,11,-1,-1,-1,0,0,-1,-1},
    {9,7,0,1,17,-1,-1,7,3,-1,3,23,-1,16,-1,-1,21,-1,0,-1,-1,0,0,-1},
    {24,5,26,7,1,-1,-1,15,24,15,-1,8,-1,13,-1,13,-1,11,-1,-1,-1,-1,0,0},
    {2,2,19,14,24,1,15,19,-1,21,-1,2,-1,24,-1,3,-1,2,1,-1,-1,-1,-1,0}};
// (648, 5/6)
const int8_t kB_648_56[4][24] = {
    {17,13,8,21,9,3,18,12,10,0,4,15,19,2,5,10,26,19,13,13,1,0,-1,-1},
    {3,12,11,14,11,25,5,18,0,9,2,26,26,10,24,7,14,20,4,2,-1,0,0,-1},
    {22,16,4,3,10,21,12,5,21,14,19,5,-1,8,5,18,11,5,5,15,0,-1,0,0},
    {7,7,14,14,4,16,16,24,24,10,1,7,15,6,10,26,8,18,21,14,1,-1,-1,0}};
// (1296, 1/2)
const int8_t kB_1296_12[12][24] = {
    {40,-1,-1,-1,22,-1,49,23,43,-1,-1,-1,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {50,1,-1,-1,48,35,-1,-1,13,-1,30,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {39,50,-1,-1,4,-1,2,-1,-1,-1,-1,49,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1},
    {33,-1,-1,38,37,-1,-1,4,1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1},
    {45,-1,-1,-1,0,22,-1,-1,20,42,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1},
    {51,-1,-1,48,35,-1,-1,-1,44,-1,18,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1},
    {47,11,-1,-1,-1,17,-1,-1,51,-1,-1,-1,0,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1},
    {5,-1,25,-1,6,-1,45,-1,13,40,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1},
    {33,-1,-1,34,24,-1,-1,-1,23,-1,-1,46,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1},
    {1,-1,27,-1,1,-1,-1,-1,38,-1,44,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1},
    {-1,18,-1,-1,23,-1,-1,8,0,35,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0},
    {49,-1,17,-1,30,-1,-1,-1,34,-1,-1,19,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0}};
// (1296, 2/3)
const int8_t kB_1296_23[8][24] = {
    {39,31,22,43,-1,40,4,-1,11,-1,-1,50,-1,-1,-1,6,1,0,-1,-1,-1,-1,-1,-1},
    {25,52,41,2,6,-1,14,-1,34,-1,-1,-1,24,-1,37,-1,-1,0,0,-1,-1,-1,-1,-1},
    {43,31,29,0,21,-1,28,-1,-1,2,-1,-1,7,-1,17,-1,-1,-1,0,0,-1,-1,-1,-1},
    {20,33,48,-1,4,13,-1,26,-1,-1,22,-1,-1,46,42,-1,-1,-1,-1,0,0,-1,-1,-1},
    {45,7,18,51,12,25,-1,-1,-1,50,-1,-1,5,-1,-1,-1,0,-1,-1,-1,0,0,-1,-1},
    {35,40,32,16,5,-1,-1,18,-1,-1,43,51,-1,32,-1,-1,-1,-1,-1,-1,-1,0,0,-1},
    {9,24,13,22,28,-1,-1,37,-1,-1,25,-1,-1,52,-1,13,-1,-1,-1,-1,-1,-1,0,0},
    {32,22,4,21,16,-1,-1,-1,27,28,-1,38,-1,-1,-1,8,1,-1,-1,-1,-1,-1,-1,0}};
// (1296, 3/4)
const int8_t kB_1296_34[6][24] = {
    {39,40,51,41,3,29,8,36,-1,14,-1,6,-1,33,-1,11,-1,4,1,0,-1,-1,-1,-1},
    {48,21,47,9,48,35,51,-1,38,-1,28,-1,34,-1,50,-1,50,-1,-1,0,0,-1,-1,-1},
    {30,39,28,42,50,39,5,17,-1,6,-1,18,-1,20,-1,15,-1,40,-1,-1,0,0,-1,-1},
    {29,0,1,43,36,30,47,-1,49,-1,47,-1,3,-1,35,-1,34,-1,0,-1,-1,0,0,-1},
    {1,32,11,23,10,44,12,7,-1,48,-1,4,-1,9,-1,17,-1,16,-1,-1,-1,-1,0,0},
    {13,7,15,47,23,16,47,-1,43,-1,29,-1,52,-1,2,-1,53,-1,1,-1,-1,-1,-1,0}};
// (1296, 5/6)
const int8_t kB_1296_56[4][24] = {
    {48,29,37,52,2,16,6,14,53,31,34,5,18,42,53,31,45,-1,46,52,1,0,-1,-1},
    {17,4,30,7,43,11,24,6,14,21,6,39,17,40,47,7,15,41,19,-1,-1,0,0,-1},
    {7,2,51,31,46,23,16,11,53,40,10,7,46,53,33,35,-1,25,35,38,0,-1,0,0},
    {19,48,41,1,10,7,36,47,5,29,52,52,31,10,26,6,3,2,-1,51,1,-1,-1,0}};
// (1944, 1/2)
const int8_t kB_1944_12[12][24] = {
    {57,-1,-1,-1,50,-1,11,-1,50,-1,79,-1,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {3,-1,28,-1,0,-1,-1,-1,55,7,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {30,-1,-1,-1,24,37,-1,-1,56,14,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1},
    {62,53,-1,-1,53,-1,-1,3,35,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1},
    {40,-1,-1,20,66,-1,-1,22,28,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,-1},
    {0,-1,-1,-1,8,-1,42,-1,50,-1,-1,8,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1},
    {69,79,79,-1,-1,-1,56,-1,52,-1,-1,-1,0,-1,-1,-1,-1,-1,0,0,-1,-1,-1,-1},
    {65,-1,-1,-1,38,57,-1,-1,72,-1,27,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1,-1},
    {64,-1,-1,-1,14,52,-1,-1,30,-1,-1,32,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1,-1},
    {-1,45,-1,70,0,-1,-1,-1,77,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,-1},
    {2,56,-1,57,35,-1,-1,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0},
    {24,-1,61,-1,60,-1,-1,27,51,-1,-1,16,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0}};
// (1944, 2/3)
const int8_t kB_1944_23[8][24] = {
    {61,75,4,63,56,-1,-1,-1,-1,-1,-1,8,-1,2,17,25,1,0,-1,-1,-1,-1,-1,-1},
    {56,74,77,20,-1,-1,-1,64,24,4,67,-1,7,-1,-1,-1,-1,0,0,-1,-1,-1,-1,-1},
    {28,21,68,10,7,14,65,-1,-1,-1,23,-1,-1,-1,75,-1,-1,-1,0,0,-1,-1,-1,-1},
    {48,38,43,78,76,-1,-1,-1,-1,5,36,-1,15,72,-1,-1,-1,-1,-1,0,0,-1,-1,-1},
    {40,2,53,25,-1,52,62,-1,20,-1,-1,44,-1,-1,-1,-1,0,-1,-1,-1,0,0,-1,-1},
    {69,23,64,10,22,-1,21,-1,-1,-1,-1,-1,68,23,29,-1,-1,-1,-1,-1,-1,0,0,-1},
    {12,0,68,20,55,61,-1,40,-1,-1,-1,52,-1,-1,-1,44,-1,-1,-1,-1,-1,-1,0,0},
    {58,8,34,64,78,-1,-1,11,78,24,-1,-1,-1,-1,-1,58,1,-1,-1,-1,-1,-1,-1,0}};
// (1944, 3/4)
const int8_t kB_1944_34[6][24] = {
    {48,29,28,39,9,61,-1,-1,-1,63,45,80,-1,-1,-1,37,32,22,1,0,-1,-1,-1,-1},
    {4,49,42,48,11,30,-1,-1,-1,49,17,41,37,15,-1,54,-1,-1,-1,0,0,-1,-1,-1},
    {35,76,78,51,37,35,21,-1,17,64,-1,-1,-1,59,7,-1,-1,32,-1,-1,0,0,-1,-1},
    {9,65,44,9,54,56,73,34,42,-1,-1,-1,35,-1,-1,-1,46,39,0,-1,-1,0,0,-1},
    {3,62,7,80,68,26,-1,80,55,-1,36,-1,26,-1,9,-1,72,-1,-1,-1,-1,-1,0,0},
    {26,75,33,21,69,59,3,38,-1,-1,-1,35,-1,62,36,26,-1,-1,1,-1,-1,-1,-1,0}};
// (1944, 5/6)
const int8_t kB_1944_56[4][24] = {
    {13,48,80,66,4,74,7,30,76,52,37,60,-1,49,73,31,74,73,23,-1,1,0,-1,-1},
    {69,63,74,56,64,77,57,65,6,16,51,-1,64,-1,68,9,48,62,54,27,-1,0,0,-1},
    {51,15,0,80,24,25,42,54,44,71,71,9,67,35,-1,58,-1,29,-1,53,0,-1,0,0},
    {16,29,36,41,44,56,59,37,50,24,-1,65,4,65,52,-1,4,-1,73,52,1,-1,-1,0}};

struct Code {
  int64_t L = 0, Z = 0, M = 0, K = 0, dmax = 0;
  std::vector<int64_t> idx;  // (M, dmax), check-major, padded with L
  std::vector<int64_t> deg;  // real edges per check (they come first in each row)
};

Code build(int64_t L, Rate rate) {
  const int8_t(*B)[24] = nullptr;
  int rows = 0;
  const int p = rate.first, q = rate.second;
#define SI_LDPC_CASE(LL, PP, QQ, NR)                \
  if (L == LL && p == PP && q == QQ) {             \
    B = kB_##LL##_##PP##QQ;                        \
    rows = NR;                                     \
  }
  SI_LDPC_CASE(648, 1, 2, 12) SI_LDPC_CASE(648, 2, 3, 8) SI_LDPC_CASE(648, 3, 4, 6) SI_LDPC_CASE(648, 5, 6, 4)
  SI_LDPC_CASE(1296, 1, 2, 12) SI_LDPC_CASE(1296, 2, 3, 8) SI_LDPC_CASE(1296, 3, 4, 6) SI_LDPC_CASE(1296, 5, 6, 4)
  SI_LDPC_CASE(1944, 1, 2, 12) SI_LDPC_CASE(1944, 2, 3, 8) SI_LDPC_CASE(1944, 3, 4, 6) SI_LDPC_CASE(1944, 5, 6, 4)
#undef SI_LDPC_CASE
  if (!B) throw std::invalid_argument("no LDPC code for this length and rate");
  Code c;
  c.L = L;
  c.Z = L / 24;
  c.M = rows * c.Z;
  c.K = L - c.M;
  std::vector<std::vector<int64_t>> per(c.M);  // columns of each check, sorted
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < 24; ++j) {
      if (B[i][j] < 0) continue;
      for (int64_t t = 0; t < c.Z; ++t) per[i * c.Z + t].push_back(j * c.Z + (t + B[i][j]) % c.Z);
    }
  for (auto& v : per) {
    std::sort(v.begin(), v.end());
    c.dmax = std::max<int64_t>(c.dmax, static_cast<int64_t>(v.size()));
  }
  c.idx.assign(c.M * c.dmax, L);
  for (const auto& v : per) c.deg.push_back(static_cast<int64_t>(v.size()));
  for (int64_t m = 0; m < c.M; ++m)
    for (size_t k = 0; k < per[m].size(); ++k) c.idx[m * c.dmax + k] = per[m][k];
  return c;
}

const Code& code(int64_t L, Rate rate) {
  static std::mutex mu;
  static std::map<std::array<int64_t, 3>, std::unique_ptr<Code>> cache;
  std::lock_guard<std::mutex> lk(mu);
  auto& p = cache[{L, rate.first, rate.second}];
  if (!p) p = std::make_unique<Code>(build(L, rate));
  return *p;
}

}  // namespace

Params params(int64_t length_bytes, int64_t ncbps, Rate rate, int64_t m_stbc) {
  const int64_t p = rate.first, q = rate.second;
  const double R = static_cast<double>(p) / static_cast<double>(q);
  Params P;
  const int64_t npld = 8 * length_bytes + 16;
  int64_t navbits = ncbps * m_stbc *
                    static_cast<int64_t>(std::ceil(static_cast<double>(npld) / (static_cast<double>(ncbps) * R * static_cast<double>(m_stbc))));
  int64_t ncw, L;
  const double nav = static_cast<double>(navbits), npl = static_cast<double>(npld);
  if (navbits <= 648) {
    ncw = 1;
    L = nav >= npl + 912 * (1 - R) ? 1296 : 648;
  } else if (navbits <= 1296) {
    ncw = 1;
    L = nav >= npl + 1464 * (1 - R) ? 1944 : 1296;
  } else if (navbits <= 1944) {
    ncw = 1;
    L = 1944;
  } else if (navbits <= 2592) {
    ncw = 2;
    L = nav >= npl + 2916 * (1 - R) ? 1944 : 1296;
  } else {
    ncw = static_cast<int64_t>(std::ceil(npl / (1944 * R)));
    L = 1944;
  }
  const int64_t K = L * p / q;
  const int64_t nshrt = std::max<int64_t>(0, ncw * K - npld);
  int64_t npunc = std::max<int64_t>(0, ncw * L - navbits - nshrt);
  const double ncwd = static_cast<double>(ncw), Ld = static_cast<double>(L);
  if ((static_cast<double>(npunc) > 0.1 * ncwd * Ld * (1 - R) &&
       static_cast<double>(nshrt) < 1.2 * static_cast<double>(npunc) * R / (1 - R)) ||
      static_cast<double>(npunc) > 0.3 * ncwd * Ld * (1 - R)) {
    navbits += ncbps * m_stbc;
    npunc = std::max<int64_t>(0, ncw * L - navbits - nshrt);
  }
  const double rep = static_cast<double>(navbits) - static_cast<double>(ncw * L) * (1 - R) - npl;
  P.npld = npld;
  P.navbits = navbits;
  P.nsym = navbits / ncbps;
  P.ncw = ncw;
  P.L = L;
  P.K = K;
  P.nshrt = nshrt;
  P.npunc = npunc;
  P.nrep = std::max<int64_t>(0, static_cast<int64_t>(std::nearbyint(rep)));
  return P;
}

std::vector<int64_t> split(int64_t n, int64_t ncw) {
  std::vector<int64_t> v(ncw);
  const int64_t b = n / ncw;
  for (int64_t i = 0; i < ncw; ++i) v[i] = b + (i < n % ncw ? 1 : 0);
  return v;
}

std::vector<Seg> layout(const Params& P) {
  const auto s = split(P.nshrt, P.ncw), u = split(P.npunc, P.ncw), r = split(P.nrep, P.ncw);
  std::vector<Seg> out;
  for (int64_t i = 0; i < P.ncw; ++i) out.push_back({P.K - s[i], P.L - P.K - u[i], r[i]});
  return out;
}

std::vector<int64_t> offsets(const Params& P) {
  std::vector<int64_t> o;
  int64_t a = 0;
  for (const auto& s : layout(P)) {
    o.push_back(a);
    a += s.k + s.par + s.rep;
  }
  return o;
}

std::vector<uint8_t> transmitted(const std::vector<uint8_t>& cw, const Params& P, size_t i) {
  const Seg s = layout(P)[i];
  std::vector<uint8_t> tx(cw.begin(), cw.begin() + s.k);
  tx.insert(tx.end(), cw.begin() + P.K, cw.begin() + P.K + s.par);
  if (s.rep && !tx.empty()) {
    const size_t n = tx.size();
    for (int64_t j = 0; j < s.rep; ++j) tx.push_back(tx[static_cast<size_t>(j) % n]);
  }
  return tx;
}

Decoded decode(const dvec& llr, int64_t L, Rate rate, int iters, double alpha) {
  const Code& cd = code(L, rate);
  const int64_t E = cd.M * cd.dmax;
  dvec lam(L + 1), tot(L + 1), acc(L + 1), R(E, 0.0), Q(cd.dmax), a(cd.dmax), s(cd.dmax);
  for (int64_t v = 0; v < L; ++v) lam[v] = std::min(std::max(llr[v], -1e4), 1e4);
  lam[L] = 1e6;
  Decoded out;
  out.hard.resize(L);
  for (int64_t v = 0; v < L; ++v) out.hard[v] = lam[v] < 0;
  out.used = iters;
  std::vector<uint8_t> hard(L + 1, 0);
  const double inf = std::numeric_limits<double>::infinity();
  // Padding edges (index L, the dummy "sure 0" variable) add 0 to the dummy's sum, carry
  // Q = 1e6 and get R = 0, so only their 1e6 enters the check update (as Python's dense rows).
  const double kPad = 1e6;
  for (int it = 0; it <= iters; ++it) {
    std::fill(acc.begin(), acc.end(), 0.0);
    for (int64_t m = 0; m < cd.M; ++m) {
      const int64_t* id = &cd.idx[m * cd.dmax];
      const double* r = &R[m * cd.dmax];
      for (int64_t k = 0; k < cd.deg[m]; ++k) acc[id[k]] += r[k];
    }
    for (int64_t v = 0; v < L; ++v) {
      tot[v] = lam[v] + acc[v];
      hard[v] = tot[v] < 0;
    }
    bool ok = true;
    for (int64_t m = 0; m < cd.M && ok; ++m) {
      const int64_t* id = &cd.idx[m * cd.dmax];
      uint8_t x = 0;
      for (int64_t k = 0; k < cd.deg[m]; ++k) x ^= hard[id[k]];
      ok = x == 0;
    }
    if (ok) {
      std::copy(hard.begin(), hard.begin() + L, out.hard.begin());
      out.used = it;
      out.ok = true;
      return out;
    }
    if (it == iters) break;
    for (int64_t m = 0; m < cd.M; ++m) {
      const int64_t* id = &cd.idx[m * cd.dmax];
      double* r = &R[m * cd.dmax];
      const int64_t dg = cd.deg[m], npad = cd.dmax - dg;
      double sp = 1;
      for (int64_t k = 0; k < dg; ++k) {
        Q[k] = tot[id[k]] - r[k];
        a[k] = std::fabs(Q[k]);
        s[k] = Q[k] < 0 ? -1.0 : 1.0;
        sp *= s[k];
      }
      int64_t am = 0;
      for (int64_t k = 1; k < dg; ++k)
        if (a[k] < a[am]) am = k;
      double m1 = a[am], m2 = inf;
      if (npad > 0 && kPad < m1) {  // the first padding edge is the minimum
        am = dg;
        m1 = kPad;
        for (int64_t k = 0; k < dg; ++k)
          if (a[k] < m2) m2 = a[k];
        if (npad > 1 && kPad < m2) m2 = kPad;
      } else {
        for (int64_t k = 0; k < dg; ++k)
          if (k != am && a[k] < m2) m2 = a[k];
        if (npad > 0 && kPad < m2) m2 = kPad;
      }
      for (int64_t k = 0; k < dg; ++k) r[k] = alpha * sp * s[k] * (k == am ? m2 : m1);
    }
  }
  for (int64_t v = 0; v < L; ++v) out.hard[v] = tot[v] < 0;
  return out;
}

StreamDecoded decode_stream(const dvec& llr, const Params& P, Rate rate, int iters) {
  const int64_t L = P.L, K = P.K;
  const auto lay = layout(P);
  StreamDecoded out;
  int64_t pos = 0;
  for (const auto& sg : lay) {
    const int64_t n = sg.k + sg.par + sg.rep;
    dvec seg(n, 0.0);
    for (int64_t i = 0; i < n && pos + i < static_cast<int64_t>(llr.size()); ++i) seg[i] = llr[pos + i];
    pos += n;
    dvec cw(L, 0.0);
    for (int64_t i = 0; i < sg.k; ++i) cw[i] = seg[i];
    for (int64_t i = sg.k; i < K; ++i) cw[i] = 1e3;
    for (int64_t i = 0; i < sg.par; ++i) cw[K + i] = seg[sg.k + i];
    if (sg.rep) {
      const int64_t ns = sg.k + sg.par;
      for (int64_t j = 0; j < sg.rep; ++j) {
        const int64_t t = j % ns;
        cw[t < sg.k ? t : K + (t - sg.k)] += seg[sg.k + sg.par + j];
      }
    }
    Decoded d = decode(cw, L, rate, iters);
    out.info.insert(out.info.end(), d.hard.begin(), d.hard.begin() + sg.k);
    out.ok.push_back(d.ok);
    out.used.push_back(d.used);
    out.hard.push_back(std::move(d.hard));
  }
  return out;
}

}  // namespace si::rx::wifi::ldpc
