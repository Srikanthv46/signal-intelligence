// Reed-Solomon codes over GF(2^8), systematic and shortened (python core/rs.py; clean-room: field from the primitive
// polynomial, generator prod (x - alpha^(prim (fcr + i))), syndromes, Berlekamp-Massey, Chien search, Forney).
#pragma once

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace si::rs {

class RS {
 public:
  RS(unsigned gfpoly, int fcr, int prim, int nroots);
  // parity bytes of data (codeword = data + parity)
  std::vector<uint8_t> encode(const std::vector<uint8_t>& data) const;
  // codeword -> (corrected block, corrected byte count); nullopt when the errors exceed the capacity
  std::optional<std::pair<std::vector<uint8_t>, int>> decode(const std::vector<uint8_t>& block) const;

 private:
  int nroots_, fcr_, prim_;
  int exp_[512], log_[256];
  std::vector<int> gen_;
  int mul(int a, int b) const { return (a == 0 || b == 0) ? 0 : exp_[log_[a] + log_[b]]; }
};

}  // namespace si::rs
