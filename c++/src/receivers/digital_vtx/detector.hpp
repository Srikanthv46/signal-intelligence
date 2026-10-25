// Digital video links that are not DJI (python: receivers/digital_vtx/detector.py):
// Herelink (LTE Rel-8 DL PHY in 2.4 GHz ISM, PSS fold search -> carrier -> SSS -> PCI ->
// subframe activity), HDZero (27 MHz mode, continual pilots: hdzero.hpp) and every other CP-OFDM link with Tu >= 10 us, found blind (blind.hpp).
//
//   detect(x, fs, fc, report_wifi) -> Emissions, family herelink | hdzero | digital_vtx | wifi_ofdm
#pragma once

#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::digital_vtx {

// python detect(x, fs, fc, {"report_wifi": report_wifi})
std::vector<Emission> detect(const cfvec& x, double fs, double fc, bool report_wifi = true);

// python classify_numerology(tu_us, None, bw_mhz) -> (name, hits); numerology_role(name): the table's role
std::pair<std::string, std::vector<std::string>> classify_numerology(double tu_us,
                                                                     std::optional<double> bw_mhz = std::nullopt);
std::string numerology_role(const std::string& name);

namespace detail {
// python _slice(X, fs, f_off, half_bw, out_fs) -> (y, M*fs/L)
cdvec slice(const cdvec& X, double fs, double f_off, double half_bw, double out_fs, double* ofs_out = nullptr);
// python _pss_time(u, ofs)
cdvec pss_time(int u, double ofs);
// python pss_fold_search(X, fs, usable_half)
struct FoldSearch {
  dvec offs;
  std::map<int, dvec> Z;
  std::map<int, std::vector<int64_t>> T;
  double ofs = 0, P = 0, global_median = 0;
  int64_t M = 0;
  size_t nt = 0;
};
FoldSearch pss_fold_search(const cdvec& X, double fs, double usable_half, double step_hz = 10e3,
                           double period_s = 5e-3, double clip = 24.0,
                           std::optional<std::pair<double, double>> f_range = std::nullopt);
// python _carrier_search(X, fs, f0, t_inst_s) -> (f_c, root, score)
struct CarrierHit {
  double f_c = 0;
  int root = 25;
  double score = 0;
};
CarrierHit carrier_search(const cdvec& X, double fs, double f0, const dvec& t_inst_s);
// python hdzero_detect(x, fs, fc, X, ctx, info)
std::vector<Emission> hdzero_detect(const cfvec& x, double fs, double fc, const cdvec& X);
}  // namespace detail

}  // namespace si::rx::digital_vtx
