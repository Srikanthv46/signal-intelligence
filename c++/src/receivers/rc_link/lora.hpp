// LoRa / ExpressLRS blind detection (python: receivers/rc_link/detector.py, "LoRa: blind
// scan", "LoRa: refine one packet", lora_emission, detect_lora, _dedup).
#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"
#include "receivers/rc_link/dsp_util.hpp"

namespace si::rx::rc_link {

// python detect_lora(X, fs, fc, f_range=None): X = fft(x) (complex64); f_range: (lo, hi) Hz offsets from fc
// that the scanned sub-band centres must lie in
std::vector<Emission> detect_lora(const cfvec& X, double fs, double fc,
                                  const std::pair<double, double>* f_range = nullptr);

namespace detail {

struct LoraBw {
  double bw;
  std::vector<int> sfs;
};
std::vector<LoraBw> lora_bws_for(double fc);  // python _lora_bws_for

struct Cand {
  double f_sb = 0;
  int sf = 0;
  double bw = 0;
  int pol = 1;
  int64_t win = 0;
  double ratio = 0;
  int64_t bin = 0;
  int K = 8;
};
// python lora_scan(X, fs, fc, bw, sfs, K=8, eta=24.0, chunk=48, edge=0.45, f_range=None) (candidates only)
std::vector<Cand> lora_scan(const cfvec& X, double fs, double bw, const std::vector<int>& sfs, int K = 8,
                            double eta = 24.0, size_t chunk = 48, double edge = 0.45,
                            const std::pair<double, double>* f_range = nullptr);
std::vector<Cand> nms_candidates(std::vector<Cand> cands, double bw);  // python _nms_candidates

struct Partial {
  int sf = 0;
  double bw = 0;
  int pol = 1;
  int64_t N = 0;
  double ofs = 0, f_sb = 0;
  float k_u = 0;
  int64_t pre_start = 0, pre_end = 0, n_pre_est = 0;
  double ratio1 = 0;
  int K = 8;
  bool sfd_window_in_capture = true;  // python: int(pre_end + 3.75 N) + 2 N <= n
};

struct Walk {
  double start_samp = 0, end_samp = 0;
  int64_t n_pre = 0, n_pay = 0;
  std::vector<int64_t> sync_vals, pay_vals, all_vals;
  double snr_db = 0, slope = 0, slope_err = 0;
  double pre_ratio = 0, pay_ratio = 0, z = 0, z_sfd = 0;
  fvec sfd_r;
  double f_carrier = 0, ofs = 0, bw = 0;
  int sf = 0, pol = 1;
  int64_t N = 0;
  double sfd_ratio = 0, ratio1 = 0;
};

struct Refined {
  bool partial = true;  // python w["partial"]
  Partial prem;         // the preamble-only view (w itself when partial, else w["prem"])
  Walk walk;            // valid when !partial
};
// slices of one sub-band (cand.f_sb, bw) by t_adv, shared by its candidates (python lora_refine cache)
using SliceCache = std::map<double, SliceOut>;
// python lora_refine(X, fs, fc, cand, cache=None): nullopt = None
std::optional<Refined> lora_refine(const cfvec& X, double fs, double fc, const Cand& cand, int n_iter = 3,
                                   double min_sfd_ratio = 5.0, SliceCache* cache = nullptr);

struct Streams {
  cfvec z_up, z_dn;
  double ofs = 0, f_act = 0;
};
// python _streams(X, fs, f_off, bw, N, pol, t_adv, cache)
Streams streams(const cfvec& X, double fs, FOff f_off, double bw, int64_t N, int pol, double t_adv = 0.0,
                SliceCache* cache = nullptr);
// python _circ_mean(v, period, w)
double circ_mean(const dvec& v, int64_t period, const dvec& w);
// python _walk_packet(z_up, z_dn, b_sfd, N, cfo, pol) (start/end samples without the frac shift)
std::optional<Walk> walk_packet(const cfvec& z_up, const cfvec& z_dn, double b_sfd, int64_t N, double cfo, int pol,
                                int max_pre = 40, int max_pay = 400);

Emission lora_emission(const Walk& w, double fc, double fs);
Emission lora_preamble_emission(const Partial& p, double fc);
std::vector<Emission> dedup(std::vector<Emission> ems, double f_tol_frac = 0.5);

}  // namespace detail
}  // namespace si::rx::rc_link
