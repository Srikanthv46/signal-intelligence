#include "receivers/rc_link/lora.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <tuple>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "receivers/rc_link/dsp_util.hpp"
#include "receivers/rc_link/elrs.hpp"
#include "receivers/rc_link/mlrs.hpp"
#include "receivers/rc_link/lora_phy.hpp"
#include "receivers/rc_link/tables.hpp"

namespace si::rx::rc_link {
namespace detail {

// python _streams(X, fs, f_off, bw, N, pol, t_adv, cache)
Streams streams(const cfvec& X, double fs, FOff f_off, double bw, int64_t N, int pol, double t_adv,
                SliceCache* cache) {
  SliceOut own;
  if (!cache) {
    own = slice(X, fs, f_off, bw / 2, bw, t_adv);
  } else if (!cache->count(t_adv)) {
    (*cache)[t_adv] = slice(X, fs, f_off, bw / 2, bw, t_adv);
  }
  const SliceOut& s = cache ? cache->at(t_adv) : own;
  const cfvec u = upchirp(N);
  Streams o;
  o.z_up.resize(s.y.size());
  o.z_dn.resize(s.y.size());
  // y * np.conj(uu): the conj temporary may be elided (operands swapped); y * uu never is
  const bool el = elides(s.y.size(), sizeof(cf));
  for (size_t i = 0; i < s.y.size(); ++i) {
    const cf uu = u[i % static_cast<size_t>(N)];
    const cf zc = mul_conj(s.y[i], uu, el), zp = np::cmul(s.y[i], uu);
    o.z_up[i] = pol > 0 ? zc : zp;
    o.z_dn[i] = pol > 0 ? zp : zc;
  }
  o.ofs = s.ofs;
  o.f_act = s.f_act;
  return o;
}

namespace {

// max of S over (kb + [-tol, tol]) mod n
float max_around(const fvec& S, int64_t kb, int tol) {
  const int64_t n = static_cast<int64_t>(S.size());
  float m = -std::numeric_limits<float>::infinity();
  for (int d = -tol; d <= tol; ++d) m = std::max(m, S[static_cast<size_t>(pymod(kb + d, n))]);
  return m;
}

// python r_at(z, st, kexp, tol) inside lora_refine (st already rounded)
double r_at(const cfvec& z, int64_t st, float kexp, int64_t N, int tol = 2) {
  auto p = sym_spec(z, st, N);
  if (!p) return 0.0;
  const int64_t kb = pymod(rint_i(kexp * 8.0f), 8 * N);
  return static_cast<double>(p->per_noise(max_around(p->S, kb, tol), N));
}

fvec power_spec(const cf* seg, size_t len, size_t nfft) {
  cfvec F = fft::fft(seg, len, nfft);
  fvec S(nfft);
  for (size_t i = 0; i < nfft; ++i) S[i] = np::cabs2(F[i]);
  return S;
}

}  // namespace

// python _circ_mean(v, period, w)
double circ_mean(const dvec& v, int64_t period, const dvec& w) {
  cdvec t(v.size());
  const double inv = 1.0 / static_cast<double>(period);
  for (size_t i = 0; i < v.size(); ++i) {
    const double th = (2 * kPi * v[i]) * inv;
    t[i] = np::cmul(cd(w[i], 0.0), cd(std::cos(th), std::sin(th)));
  }
  const cd s = np::pairwise_sum(t.data(), t.size());
  const double a = std::atan2(s.imag(), s.real());
  return np::mod<double>(a / (2 * kPi) * static_cast<double>(period), static_cast<double>(period));
}

// python _walk_packet(z_up, z_dn, b_sfd, N, cfo, pol, max_pre=40, max_pay=400)
std::optional<Walk> walk_packet(const cfvec& z_up, const cfvec& z_dn, double b, int64_t N, double cfo, int pol,
                                int max_pre, int max_pay) {
  const double Nd = static_cast<double>(N);
  const double base = cfo - pol * np::mod<double>(b, Nd);
  struct Pre {
    int64_t m;
    double dv;
    float r;
    bool good;
  };
  struct Sym {
    double g, d;
    float w;
  };
  std::vector<Pre> pre;
  std::vector<Sym> syms;
  int miss = 0;
  const int64_t kb_base = pymod(rint_i(base * 8), 8 * N);
  for (int64_t m = 3; m < max_pre + 3; ++m) {
    auto p = sym_spec(z_up, rint_i(b - static_cast<double>(m * N)), N);
    if (!p) break;
    const float r_exp = p->per_noise(max_around(p->S, kb_base, 1), N);
    if (r_exp < 4.0f) {
      if (++miss > 1) break;
      continue;
    }
    miss = 0;
    const double dv = wrap<double>(static_cast<double>(p->pk) - base, Nd);
    pre.push_back({m, std::fabs(dv) < 0.5 ? dv : 0.0, r_exp, std::fabs(dv) < 0.5 && p->r > 6});
  }
  if (pre.size() < 2) return std::nullopt;
  int64_t mmax = 0;
  for (const auto& q : pre) mmax = std::max(mmax, q.m);
  Walk w;
  w.n_pre = mmax - 2;
  for (const auto& q : pre)
    if (q.good) syms.push_back({static_cast<double>(-q.m), q.dv, q.r});
  const int64_t kb_sfd = pymod(rint_i((cfo + pol * np::mod<double>(b, Nd)) * 8), 8 * N);
  for (int m = 0; m < 2; ++m) {
    auto p = sym_spec(z_dn, rint_i(b + static_cast<double>(m * N)), N);
    if (p) w.sfd_r.push_back(p->per_noise(max_around(p->S, kb_sfd, 1), N));
  }
  for (int m : {2, 1}) {
    auto p = sym_spec(z_up, rint_i(b - static_cast<double>(m * N)), N);
    if (p) w.sync_vals.push_back(pymod(rint_i(pol * wrap<double>(static_cast<double>(p->pk) - base, Nd)), N));
  }
  const double p0 = b + 2.25 * Nd;
  const double base_pay = base - static_cast<double>(pol * N) / 4.0;
  struct Pay {
    int64_t j;
    float pk, r;
  };
  std::vector<Pay> pay;
  for (int64_t j = 0; j < max_pay; ++j) {
    auto p = sym_spec(z_up, rint_i(p0 + static_cast<double>(j * N)), N);
    if (!p) break;
    pay.push_back({j, p->pk, p->r});
  }
  if (pay.empty()) return std::nullopt;
  fvec r(pay.size());
  for (size_t i = 0; i < pay.size(); ++i) r[i] = pay[i].r;
  const double c_noise = std::log(Nd) + 2.0;
  fvec r_pre(pre.size());
  for (size_t i = 0; i < pre.size(); ++i) r_pre[i] = pre[i].r;
  const float sp = 0.85f * mean_f(r_pre);
  const double sig = (c_noise + 1.0 > static_cast<double>(sp)) ? c_noise + 1.0 : static_cast<double>(sp);
  const double c = c_noise + 0.5 * (sig - c_noise);
  dvec a(r.begin(), r.end());
  a.insert(a.end(), 3, c_noise);
  const dvec ma = np::convolve(a, dvec(4, 0.25), np::ConvMode::Valid);
  size_t run_end = r.size();
  for (size_t i = 0; i < ma.size(); ++i)
    if (ma[i] < c_noise + 0.25 * (sig - c_noise)) {
      run_end = i;
      break;
    }
  run_end = std::min(run_end + 3, r.size());
  std::vector<int64_t> vals(pay.size());
  for (size_t i = 0; i < pay.size(); ++i)
    vals[i] = pymod(rint_i(pol * wrap<double>(static_cast<double>(pay[i].pk) - base_pay, Nd)), N);
  // a following packet's preamble repeats one value; random payload does so `run` times with
  // probability (3/N)^(run-1) per symbol, so small N needs a longer run (python _walk_packet)
  const size_t run = N >= 128 ? 3 : N >= 64 ? 4 : 5;
  for (size_t j = 0; j + run <= vals.size(); ++j) {
    bool same = true;
    double rmin = static_cast<double>(r[j]);
    for (size_t i = 1; i < run; ++i) {
      same = same && std::fabs(wrap<double>(static_cast<double>(vals[j + i] - vals[j]), Nd)) <= 1;
      rmin = std::min(rmin, static_cast<double>(r[j + i]));
    }
    if (same && rmin > c) {
      run_end = std::min(run_end, j);
      break;
    }
  }
  const size_t ncs = std::max<size_t>(run_end, 1);
  dvec cs(ncs);
  double acc = 0;
  for (size_t i = 0; i < ncs; ++i) {
    acc = i == 0 ? static_cast<double>(r[0]) - c : acc + (static_cast<double>(r[i]) - c);
    cs[i] = acc;
  }
  const int64_t n_pay = static_cast<int64_t>(np::argmax(cs)) + 1;
  for (int64_t j = 0; j < n_pay; ++j) {
    const double v = pol * wrap<double>(static_cast<double>(pay[j].pk) - base_pay, Nd);
    const double vi = std::nearbyint(v);
    w.pay_vals.push_back(pymod(static_cast<int64_t>(vi), N));
    syms.push_back({2.25 + static_cast<double>(j), pol * (v - vi), pay[j].r});
  }
  // in-band SNR from the preamble power at the known bin: r = 1 + N*snr
  const float t = (mean_f(r_pre) - 1.0f) / static_cast<float>(N);
  w.snr_db = (1e-4 > static_cast<double>(t)) ? 10 * std::log10(1e-4)
                                             : static_cast<double>(10.0f * std::log10(t));
  const size_t n_ev = r_pre.size() + w.sfd_r.size();
  float ssum = 0.0f;
  for (size_t i = 0; i < w.sfd_r.size(); ++i) ssum = i == 0 ? w.sfd_r[0] : ssum + w.sfd_r[i];
  const float zp = (np::pairwise_sum(r_pre.data(), r_pre.size()) + ssum) - static_cast<float>(n_ev);
  w.z = static_cast<double>(zp) / std::sqrt(static_cast<double>(n_ev));
  w.z_sfd = w.sfd_r.empty() ? 0.0
                            : static_cast<double>(ssum - static_cast<float>(w.sfd_r.size())) /
                                  std::sqrt(static_cast<double>(w.sfd_r.size()));
  // chip-clock drift: fractional peak deviation vs symbol index (weighted LS)
  w.slope = w.slope_err = std::numeric_limits<double>::quiet_NaN();
  dvec g, d;
  fvec ww;
  for (const auto& s : syms)
    if (s.w > 8.0f) {
      g.push_back(s.g);
      d.push_back(s.d);
      ww.push_back(s.w);
    }
  if (g.size() >= 6 && *std::max_element(g.begin(), g.end()) - *std::min_element(g.begin(), g.end()) > 8) {
    const float wmax = *std::max_element(ww.begin(), ww.end());
    const size_t n = g.size();
    dvec W(n);
    for (size_t i = 0; i < n; ++i) W[i] = static_cast<double>(ww[i] / wmax);
    dmat A(2, 2, 0.0);
    dvec bb(2, 0.0);
    for (size_t i = 0; i < n; ++i) {
      const double gw = g[i] * W[i], wd = W[i] * d[i];
      A(0, 0) += W[i];
      A(0, 1) += gw;
      A(1, 0) += g[i] * W[i];
      A(1, 1) += g[i] * gw;
      bb[0] += wd;
      bb[1] += g[i] * wd;
    }
    try {
      const dvec coef = np::solve(A, bb);
      dvec q(n);
      for (size_t i = 0; i < n; ++i) {
        const double res = d[i] - (coef[0] + g[i] * coef[1]);
        q[i] = W[i] * (res * res);
      }
      const double s2 = np::pairwise_sum(q.data(), n) / static_cast<double>(std::max<int64_t>(n - 2, 1));
      w.slope = coef[1];
      w.slope_err = std::sqrt(np::inv(A)(1, 1) * s2);
    } catch (const std::exception&) {
      // python: LinAlgError would propagate; a singular 2x2 needs all g equal (excluded by ptp > 8)
    }
  }
  w.start_samp = b - static_cast<double>((w.n_pre + 2) * N);
  w.end_samp = p0 + static_cast<double>(n_pay * N);
  w.n_pay = n_pay;
  w.pre_ratio = static_cast<double>(np::median(r_pre));
  w.pay_ratio = static_cast<double>(np::median(fvec(r.begin(), r.begin() + n_pay)));
  w.all_vals = vals;
  return w;
}

namespace {

enum class WalkStatus { Lost, NoPacket, Ok };

// python _refine_walk(...): cfo and b_sfd arrive as numpy float32 and become float64 after
// the first refinement step (float32 + float64)
WalkStatus refine_walk(const cfvec& X, double fs, double fc, double f_act, double ofs, float cfo32, float b32,
                       int64_t N, double bw, int pol, int n_iter, Walk& out) {
  const int sf = static_cast<int>(std::log2(static_cast<double>(N)));
  const double Nd = static_cast<double>(N);
  double f_ref = f_act, cfo = cfo32, b = b32;
  bool f32 = true;
  for (int it = 0; it < n_iter; ++it) {
    Streams st;
    if (f32) {
      const float f_new = static_cast<float>(f_ref) + (cfo32 * static_cast<float>(ofs)) / static_cast<float>(N);
      st = streams(X, fs, FOff{f_new, true}, bw, N, pol);
      cfo32 = cfo32 - static_cast<float>((st.f_act - f_ref) / (st.ofs / Nd));
    } else {
      const double f_new = f_ref + (cfo * ofs) / Nd;
      st = streams(X, fs, FOff{f_new, false}, bw, N, pol);
      cfo -= (st.f_act - f_ref) / (st.ofs / Nd);
    }
    ofs = st.ofs;
    f_ref = st.f_act;
    auto start = [&](int64_t m) {  // b_sfd + m*N, rounded
      return f32 ? rint_i(b32 + static_cast<float>(m * N)) : rint_i(b + static_cast<double>(m * N));
    };
    dvec ku_v, ku_w, kd_v, kd_w;
    for (int64_t m = 3; m < 9; ++m)
      if (auto p = sym_spec(st.z_up, start(-m), N)) {
        ku_v.push_back(p->pk);
        ku_w.push_back(p->r);
      }
    for (int64_t m = 0; m < 2; ++m)
      if (auto p = sym_spec(st.z_dn, start(m), N)) {
        kd_v.push_back(p->pk);
        kd_w.push_back(p->r);
      }
    if (ku_v.size() < 3 || kd_v.empty()) return WalkStatus::Lost;
    const double ku = circ_mean(ku_v, N, ku_w), kd = circ_mean(kd_v, N, kd_w);
    if (f32) {
      const float tb = np::mod<float>(b32, static_cast<float>(N));
      const double du = wrap<double>(ku - static_cast<double>(cfo32 - static_cast<float>(pol) * tb), Nd);
      const double dd = wrap<double>(kd - static_cast<double>(cfo32 + static_cast<float>(pol) * tb), Nd);
      cfo = static_cast<double>(cfo32) + (du + dd) / 2;
      b = static_cast<double>(b32) + pol * (dd - du) / 2;
      f32 = false;
    } else {
      const double tb = np::mod<double>(b, Nd);
      const double du = wrap<double>(ku - (cfo - pol * tb), Nd);
      const double dd = wrap<double>(kd - (cfo + pol * tb), Nd);
      cfo += (du + dd) / 2;
      b += pol * (dd - du) / 2;
    }
  }
  if (f32) {
    cfo = cfo32;
    b = b32;
  }
  // final pass on the chip-aligned stream (fractional timing removed)
  const double frac = b - std::floor(b);
  Streams st = streams(X, fs, FOff{f_ref, false}, bw, N, pol, frac);
  auto w = walk_packet(st.z_up, st.z_dn, std::floor(b), N, cfo, pol);
  if (!w) return WalkStatus::NoPacket;
  w->start_samp += frac;
  w->end_samp += frac;
  w->f_carrier = fc + st.f_act + cfo * st.ofs / Nd;
  w->ofs = st.ofs;
  w->sf = sf;
  w->bw = bw;
  w->pol = pol;
  w->N = N;
  out = std::move(*w);
  return WalkStatus::Ok;
}

int64_t lora_payload_symbols(int sf, int pl_bytes, int cr_den) {  // crc=False, implicit=True, ldro=False
  const int num = 8 * pl_bytes - 4 * sf + 28 - 20;
  return 8 + std::max<int64_t>(static_cast<int64_t>(std::ceil(static_cast<double>(num) / (4 * sf))) * cr_den, 0);
}

struct ElrsMatch {
  std::string mode;
  bool pre_ok, pay_ok;
};

std::vector<ElrsMatch> elrs_match(int sf, double bw, int64_t n_pre, int64_t n_pay) {
  std::vector<ElrsMatch> out;
  for (const auto& m : elrs_lora()) {
    if (m.sf != sf || std::fabs(m.bw - bw) > 1e3) continue;
    const int64_t exp_pay = lora_payload_symbols(sf, m.payload, m.cr_den);
    out.push_back({m.name, std::llabs(n_pre - m.preamble) <= 3, std::llabs(n_pay - exp_pay) <= 4});
  }
  return out;
}

std::string kfmt(double bw) { return fmt("%g", bw / 1e3); }

// python _mlrs_match(sf, bw, n_pre, n_pay, fc_hz): mLRS air modes consistent with the measured LoRa parameters
struct MlrsMatch { std::string mode; bool decodable; };
std::vector<MlrsMatch> mlrs_match(int sf, double bw, int64_t n_pre, int64_t n_pay, double fc_hz) {
  std::vector<MlrsMatch> out;
  const std::string band = fc_hz > 1.5e9 ? "2G4" : "900";
  for (const auto& m : mlrs::modes()) {
    if (band != m.band || m.sf != sf || std::fabs(m.bw - bw) > 1e3) continue;
    const int64_t exp_pay = lora_payload_symbols(sf, m.frame, m.cr_den);
    if (std::llabs(n_pre - m.preamble) <= 3 &&
        static_cast<double>(std::llabs(n_pay - exp_pay)) <= 0.1 * static_cast<double>(exp_pay) + 4)
      out.push_back({m.name, m.decodable});
  }
  return out;
}

// python _lora_payload(w, sf, bw); SF5: implicit header, CRC-verified
Json lora_payload(const Walk& w, int sf, double bw) {
  auto d = sf == 5 ? lora_phy::decode_sf5_implicit(w.all_vals, w.n_pay) : lora_phy::decode(w.all_vals, sf, bw);
  if (!d) return nullptr;
  const bool crc_ok = (*d).contains("crc_ok") && (*d)["crc_ok"].is_boolean() && (*d)["crc_ok"].get<bool>();
  if (!crc_ok && std::llabs((*d)["n_symbols"].get<int64_t>() - w.n_pay) > 4) return nullptr;
  const auto sw = lora_phy::sync_word(w.sync_vals);
  (*d)["sync_word"] = sw ? Json(fmt("0x%02X", *sw)) : Json(nullptr);
  (*d)["network"] = lora_phy::sync_word_network(sw);
  if ((*d).value("complete", false)) {
    const std::string h = (*d)["payload_hex"].get<std::string>();
    std::vector<uint8_t> pay;
    for (size_t i = 0; i + 1 < h.size(); i += 2) pay.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
    (*d)["payload_ascii"] = lora_phy::ascii_view(pay);
    if (crc_ok) {
      if (auto hdr = lora_phy::identify(sw, pay)) {
        (*d)["protocol"] = (*hdr)["protocol"];
        hdr->erase("protocol");
        (*d)["fields"] = *hdr;
      }
    }
    if (!(*d).contains("protocol")) {  // MAVLink in a transparent LoRa serial link (research, PLAN G2)
      Json mv = lora_phy::mavlink_in(pay);
      if (!mv.is_null()) {
        (*d)["protocol"] = "MAVLink over LoRa (transparent serial link)";
        (*d)["fields"] = mv;
      }
    }
  }
  return *d;
}

bool truthy_crc(const Json& dec) {
  return dec.is_object() && dec.contains("crc_ok") && dec["crc_ok"].is_boolean() && dec["crc_ok"].get<bool>();
}

}  // namespace

std::vector<LoraBw> lora_bws_for(double fc) {
  if (fc > 1.5e9) return {{812.5e3, {5, 6, 7, 8}}, {1625e3, {6, 7, 8}}};
  return {{500e3, {5, 6, 7, 8, 9}}, {250e3, {7, 8, 9, 10}}, {125e3, {7, 8, 9, 10, 11, 12}}};
}

// python ETA_SF5: SF5 stage-1 threshold (N = 32 opens many more windows; ~-8 dB, SF5 decodes from ~0 dB)
constexpr double kEtaSf5 = 40.0;

std::vector<Cand> lora_scan(const cfvec& X, double fs, double bw, const std::vector<int>& sfs, int K, double eta,
                            size_t chunk, double edge, const std::pair<double, double>* f_range) {
  const int64_t L = static_cast<int64_t>(X.size());
  const int64_t M = py_round_int(static_cast<double>(L) * bw / fs);
  const int64_t hb = static_cast<int64_t>(bw / 2 / fs * static_cast<double>(L));
  dvec centers = np::arange<double>(-edge * fs + bw / 2, edge * fs - bw / 2 + 1, bw / 4);
  if (f_range) {
    dvec keep;
    for (double c : centers)
      if (c >= f_range->first && c <= f_range->second) keep.push_back(c);
    centers = std::move(keep);
  }
  std::vector<cfvec> refs;
  std::vector<cfvec> refs_conj;
  for (int sf : sfs) {
    refs.push_back(upchirp(int64_t(1) << sf));
    refs_conj.push_back(refs.back());
    for (auto& v : refs_conj.back()) v = std::conj(v);
  }
  struct Hit {
    int64_t w;
    float r;
    int64_t bin;
  };
  // hits[center][sf index][pol index]
  using PerCenter = std::vector<std::vector<std::vector<Hit>>>;
  std::vector<PerCenter> hits(centers.size());
  parallel_for(centers.size(), [&](size_t ci) {
    PerCenter& H = hits[ci];
    H.assign(sfs.size(), std::vector<std::vector<Hit>>(2));
    cfvec Y(static_cast<size_t>(M), cf(0.f, 0.f));
    const int64_t k0 = py_round_int(centers[ci] / fs * static_cast<double>(L));
    for (int64_t k = -hb; k <= hb; ++k) Y[static_cast<size_t>(pymod(k, M))] = X[static_cast<size_t>(pymod(k0 + k, L))];
    const cfvec y = fft::ifft(Y);
    cfvec z, F;
    fvec S, Sa, S2;
    for (size_t si = 0; si < sfs.size(); ++si) {
      const int64_t N = int64_t(1) << sfs[si];
      const int64_t nsym = M / N;
      if (nsym < K) continue;
      const int64_t KN = K * N, hop = N * (K / 2);
      const int64_t nwin = (nsym * N - KN) / hop + 1;
      z.resize(static_cast<size_t>(nsym * N));
      F.resize(static_cast<size_t>(KN));
      S.resize(KN);
      Sa.resize(KN);
      S2.resize(KN);
      for (int pi = 0; pi < 2; ++pi) {
        const int pol = pi == 0 ? 1 : -1;
        // y[:, :nsym*N] * (conj(u) or u)[None, :] (a broadcast view: never elided)
        const cfvec& u = pol > 0 ? refs_conj[si] : refs[si];
        for (int64_t s0 = 0; s0 < nsym * N; s0 += N) np::cmul_n(y.data() + s0, u.data(), z.data() + s0, N);
        for (int64_t w = 0; w < nwin; ++w) {
          fft::c2c(z.data() + w * hop, F.data(), static_cast<size_t>(KN), true);
          for (int64_t k = 0; k < KN; ++k) S[k] = F[k].real() * F[k].real() + F[k].imag() * F[k].imag();
          // Sa = S + roll(S, -1); S2 = Sa + roll(Sa, -K)
          for (int64_t k = 0; k + 1 < KN; ++k) Sa[k] = S[k] + S[k + 1];
          Sa[KN - 1] = S[KN - 1] + S[0];
          for (int64_t k = 0; k < KN - K; ++k) S2[k] = Sa[k] + Sa[k + K];
          for (int64_t k = KN - K; k < KN; ++k) S2[k] = Sa[k] + Sa[k + K - KN];
          const size_t kb = np::argmax(S2);
          const float pk = S2[kb];
          const float r = pk * static_cast<float>(KN - 4) /
                          (np::pairwise_sum(S.data(), S.size()) - pk + static_cast<float>(1e-30));
          if (r > static_cast<float>(sfs[si] == 5 ? std::max(eta, kEtaSf5) : eta))
            H[si][pi].push_back({w * hop, r, static_cast<int64_t>(kb)});
        }
      }
    }
  });
  std::vector<Cand> cands;
  for (size_t c0 = 0; c0 < centers.size(); c0 += chunk) {
    const size_t c1 = std::min(centers.size(), c0 + chunk);
    for (size_t si = 0; si < sfs.size(); ++si)
      for (int pi = 0; pi < 2; ++pi)
        for (size_t ci = c0; ci < c1; ++ci)
          for (const auto& h : hits[ci][si][pi]) {
            Cand c;
            c.f_sb = centers[ci];
            c.sf = sfs[si];
            c.bw = bw;
            c.pol = pi == 0 ? 1 : -1;
            c.win = h.w;
            c.ratio = h.r;
            c.bin = h.bin;
            c.K = K;
            cands.push_back(c);
          }
  }
  return cands;
}

std::vector<Cand> nms_candidates(std::vector<Cand> cands, double bw) {
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return -a.ratio < -b.ratio; });
  std::vector<Cand> keep;
  for (const auto& c : cands) {
    const double T = std::pow(2.0, c.sf) / c.bw, t0 = static_cast<double>(c.win) / c.bw;
    bool dup = false;
    for (const auto& k : keep) {
      const double Tk = std::pow(2.0, k.sf) / k.bw;
      if (std::fabs(k.f_sb - c.f_sb) < 0.8 * std::max(bw, k.bw) &&
          std::fabs(static_cast<double>(k.win) / k.bw - t0) < 16 * std::max(T, Tk)) {
        dup = true;
        break;
      }
    }
    if (!dup) keep.push_back(c);
  }
  return keep;
}

std::optional<Refined> lora_refine(const cfvec& X, double fs, double fc, const Cand& cand, int n_iter,
                                   double min_sfd_ratio, SliceCache* cache) {
  const int sf = cand.sf, pol = cand.pol, K = cand.K;
  const double bw = cand.bw;
  const int64_t N = int64_t(1) << sf;
  const double Nd = static_cast<double>(N);
  const Streams s = streams(X, fs, FOff{cand.f_sb, false}, bw, N, pol, 0.0, cache);
  const double ofs = s.ofs, f_act = s.f_act;
  const int64_t n = static_cast<int64_t>(s.z_up.size());
  const int64_t s0 = cand.win;
  float k_u;
  {
    const int64_t len = std::max<int64_t>(0, std::min<int64_t>(K * N, n - s0));
    const fvec S = power_spec(s.z_up.data() + s0, static_cast<size_t>(len), static_cast<size_t>(K * N * 8));
    k_u = np::mod<float>(peak_interp(S.data(), S.size(), np::argmax(S)) / static_cast<float>(8 * K),
                         static_cast<float>(N));
  }
  // preamble extent: the de-chirped preamble tone (bin k_u) is present in any N-window inside it
  constexpr int64_t kLo = -20, kHi = 33;
  dvec prof(kHi - kLo + 1);
  for (int64_t j = kLo; j <= kHi; ++j) prof[j - kLo] = r_at(s.z_up, s0 + j * N, k_u, N);
  auto P = [&](int64_t j) { return prof[j - kLo]; };
  int64_t j0 = 0;
  {
    int best_cnt = -1;
    int64_t best_abs = 0;
    for (int64_t j = -(K / 2); j <= K / 2; ++j) {
      int cnt = 0;
      for (int i = 0; i < K; ++i) cnt += P(j + i) >= 3.0;
      if (best_cnt < 0 || cnt > best_cnt || (cnt == best_cnt && -std::llabs(j) > -best_abs)) {
        best_cnt = cnt;
        best_abs = std::llabs(j);
        j0 = j;
      }
    }
  }
  {
    dvec rin(K);
    int cnt = 0;
    for (int i = 0; i < K; ++i) {
      rin[i] = P(j0 + i);
      cnt += rin[i] >= 3.0;
    }
    const double mx = *std::max_element(rin.begin(), rin.end());
    if (cnt < K - 2 || mx > 6 * std::max(np::median(rin), 1e-9)) return std::nullopt;
  }
  int64_t j_end = j0 + K - 1;
  while (j_end + 2 <= kHi && (P(j_end + 1) + P(j_end + 2)) / 2 >= 3.0) ++j_end;
  const int64_t pre_end = s0 + (j_end + 1) * N;
  double r_d = 0.0, p_d = 0;
  bool have = false;
  fvec Sd;
  for (double adv : {0.0, 0.5}) {
    Streams sa;
    if (adv != 0.0) sa = streams(X, fs, FOff{cand.f_sb, false}, bw, N, pol, adv, cache);
    const cfvec& zd = adv == 0.0 ? s.z_dn : sa.z_dn;
    const int64_t p_lo = static_cast<int64_t>(static_cast<double>(pre_end) + 0.5 * Nd);
    const int64_t p_hi = static_cast<int64_t>(static_cast<double>(pre_end) + 3.75 * Nd);
    const int64_t step = std::max<int64_t>(N / 4, 1);
    for (int64_t p = p_lo; p < p_hi; p += step) {
      if (p + 2 * N > n) break;
      fvec S = power_spec(zd.data() + p, static_cast<size_t>(2 * N), static_cast<size_t>(2 * N * 4));
      const float mx = *std::max_element(S.begin(), S.end());
      const float med = np::median(S);
      const double r = static_cast<double>(mx) / (static_cast<double>(med) / std::log(2.0));
      if (r > r_d) {
        r_d = r;
        p_d = static_cast<double>(p) + adv;
        Sd = std::move(S);
        have = true;
      }
    }
  }
  int64_t j_start = j0;
  while (j_start - 2 >= kLo && (P(j_start - 1) + P(j_start - 2)) / 2 >= 3.0) --j_start;
  Refined out;
  Partial& pr = out.prem;
  pr.sf = sf;
  pr.bw = bw;
  pr.pol = pol;
  pr.N = N;
  pr.ofs = ofs;
  pr.f_sb = fc + cand.f_sb;
  pr.k_u = k_u;
  pr.pre_start = s0 + j_start * N;
  pr.pre_end = pre_end;
  pr.n_pre_est = j_end - j_start + 1;
  pr.ratio1 = cand.ratio;
  pr.K = K;
  pr.sfd_window_in_capture =
      static_cast<int64_t>(static_cast<double>(pre_end) + 3.75 * Nd) + 2 * N <= static_cast<int64_t>(n);
  if (!have || r_d < min_sfd_ratio) return out;
  const float k_d = np::mod<float>(peak_interp(Sd.data(), Sd.size(), np::argmax(Sd)) / 8.0f, static_cast<float>(N));
  const float cfo0 = wrap<float>((k_u + k_d) / 2.0f, Nd / 2);
  // cfo is ambiguous by N/2 (tau with it); the SFD grid point by +-1 symbol
  struct Hyp {
    double sc;
    float cfo, b;
  };
  std::vector<Hyp> hyps;
  const float half = static_cast<float>(Nd / 2);
  for (float cfo : {cfo0, cfo0 < 0 ? cfo0 + half : cfo0 - half}) {
    const float Nf = static_cast<float>(N), pf = static_cast<float>(pol);
    const float tau = np::mod<float>(pf * (k_d - cfo), Nf);
    const float f_off = static_cast<float>(f_act) + (cfo * static_cast<float>(ofs)) / Nf;
    const Streams s2 = streams(X, fs, FOff{f_off, true}, bw, N, pol);
    const float c2 = cfo - static_cast<float>((s2.f_act - f_act) / (ofs / Nd));
    const float kdx = c2 + pf * tau, kux = c2 - pf * tau;
    const float jj0 = std::floor((static_cast<float>(p_d) - tau) / Nf);
    bool any = false;
    Hyp best{0, cfo, 0};
    for (float j : {jj0 - 1.0f, jj0, jj0 + 1.0f}) {
      const float b = tau + j * Nf;
      const double r3 = r_at(s2.z_up, rint_i(b - 3.0f * Nf), kux, N, 1);
      const double r4 = r_at(s2.z_up, rint_i(b - 4.0f * Nf), kux, N, 1);
      const double r5 = r_at(s2.z_up, rint_i(b - 5.0f * Nf), kux, N, 1);
      const double sc = (((r_at(s2.z_dn, rint_i(b), kdx, N, 1) + r_at(s2.z_dn, rint_i(b + Nf), kdx, N, 1)) -
                          r_at(s2.z_dn, rint_i(b - Nf), kdx, N, 1)) -
                         r_at(s2.z_dn, rint_i(b + 2.0f * Nf), kdx, N, 1)) +
                        ((r3 + r4) + r5);
      if (!any || sc > best.sc) {
        best = {sc, cfo, b};
        any = true;
      }
    }
    hyps.push_back(best);
  }
  // walk both hypotheses; keep the one with the stronger preamble + SFD
  std::stable_sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) { return -a.sc < -b.sc; });
  std::vector<Walk> walks;
  bool failed_walk = false;
  for (const auto& h : hyps) {
    Walk wk;
    const WalkStatus st = refine_walk(X, fs, fc, f_act, ofs, h.cfo, h.b, N, bw, pol, n_iter, wk);
    if (st == WalkStatus::Ok) walks.push_back(std::move(wk));
    else if (st == WalkStatus::NoPacket) failed_walk = true;
  }
  if (walks.empty()) {
    if (failed_walk) return out;
    return std::nullopt;
  }
  size_t bi = 0;
  for (size_t i = 1; i < walks.size(); ++i)
    if (walks[i].z > walks[bi].z) bi = i;
  out.partial = false;
  out.walk = std::move(walks[bi]);
  out.walk.sfd_ratio = r_d;
  out.walk.ratio1 = cand.ratio;
  return out;
}

Emission lora_emission(const Walk& w0, double fc, double fs) {
  Walk w = w0;
  const int sf = w.sf;
  const double bw = w.bw, ofs = w.ofs, f = w.f_carrier;
  const int64_t N = w.N;
  Json clock_ppm = nullptr, clock_err = nullptr;
  if (std::isfinite(w.slope)) {
    const double delta = -w.slope / static_cast<double>(N) * (w.pol > 0 ? 1 : -1);
    const double eps = ofs / (bw * (1 + delta)) - 1;
    clock_ppm = eps * 1e6;
    clock_err = w.slope_err / static_cast<double>(N) * 1e6;
  }
  Json dec = lora_payload(w, sf, bw);
  if (dec.is_object() && (dec.value("header_ok", false) || dec.value("implicit_header", false))) {
    const int64_t ns = dec["n_symbols"].get<int64_t>();
    w.end_samp = w.end_samp - static_cast<double>(w.n_pay * N) + static_cast<double>(ns * N);
    w.n_pay = ns;
  }
  const auto nom = nominal_channel(f);
  std::vector<std::string> matches_all, full;
  if (dec.is_null())
    for (const auto& m : elrs_match(sf, bw, w.n_pre, w.n_pay)) {
      matches_all.push_back(m.mode);
      if (m.pre_ok && m.pay_ok) full.push_back(m.mode);
    }
  if (dec.is_null() && !matches_all.empty() && f < 1.5e9) {
    // ExpressLRS sub-GHz SYNC packets verify themselves (python _elrs_ota; research, PLAN R2c)
    std::vector<std::pair<int, int>> seen;
    for (const auto& name : matches_all)
      for (const auto& m : elrs_lora()) {
        if (name != m.name) continue;
        const std::pair<int, int> k{m.cr_den, m.payload};
        if (std::find(seen.begin(), seen.end(), k) != seen.end()) continue;
        seen.push_back(k);
        const auto d = lora_phy::decode_implicit(w.all_vals, sf, m.cr_den - 4, m.payload, lora_phy::default_ldro(sf, bw));
        if (!d || !dec.is_null()) continue;
        Json e = elrs::parse(d->payload);
        if (e.is_object() && e["crc_ok"].is_boolean() && e["crc_ok"].get<bool>()) {
          e["cr"] = d->cr;
          e["ldro"] = false;
          e["sync_word"] = nullptr;
          e["implicit_header"] = true;
          e["fec_corrected"] = d->fec_corrected;
          e["fec_detected_errors"] = d->fec_detected_errors;
          dec = e;
        }
      }
  }
  // mLRS (research, PLAN R3 / T6): never IQ-inverted, preamble 12, 91-byte frames; the SX127x-framed SF6 mode decodes
  std::vector<MlrsMatch> ml;
  if (dec.is_null() && full.empty() && !(w.pol < 0)) ml = mlrs_match(sf, bw, w.n_pre, w.n_pay, f);
  if (std::any_of(ml.begin(), ml.end(), [](const MlrsMatch& m) { return m.decodable; })) {
    const auto d = lora_phy::decode_implicit(w.all_vals, sf, 1, mlrs::kFrameLen, false);
    if (d) {
      Json e = mlrs::parse(d->payload);
      if (e.is_object()) {
        e["cr"] = d->cr;
        e["ldro"] = false;
        e["sync_word"] = nullptr;
        e["implicit_header"] = true;
        e["fec_corrected"] = d->fec_corrected;
        dec = e;
      }
    }
  }
  std::string family, label;
  Json cfo = nullptr;
  double conf;
  const bool near = nom && std::fabs(f - nom->f) < 150e3;
  if (!ml.empty()) {
    family = "rc_fhss";
    label = std::string(truthy_crc(dec) ? "mlrs" : "mlrs_like") + "_lora_sf" + std::to_string(sf) + "_bw" + kfmt(bw) + "k";
    conf = truthy_crc(dec) ? 0.95 : 0.7;
  } else if (near && !full.empty()) {
    family = "rc_fhss";
    label = "elrs_lora_sf" + std::to_string(sf) + "_bw" + kfmt(bw) + "k";
    cfo = f - nom->f;
    conf = truthy_crc(dec) ? 0.97 : 0.9;
  } else if (near && !matches_all.empty()) {
    family = "rc_fhss";
    label = "elrs_like_lora_sf" + std::to_string(sf) + "_bw" + kfmt(bw) + "k";
    cfo = f - nom->f;
    conf = 0.7;
  } else {
    family = "lora";
    label = "lora_sf" + std::to_string(sf) + "_bw" + kfmt(bw) + "k";
    conf = 0.6;
    if (truthy_crc(dec)) {
      conf = 0.97;
      if (dec.contains("protocol") && dec["protocol"].is_string() && !dec["protocol"].get<std::string>().empty()) {
        std::string p = dec["protocol"].get<std::string>();
        std::transform(p.begin(), p.end(), p.begin(), [](unsigned char ch) { return std::tolower(ch); });
        label = p + "_" + label;
      }
    } else if (!dec.is_null()) {
      conf = 0.8;  // header checksum passed, payload unverified
    }
  }
  const double t0 = w.start_samp / ofs * 1e6, t1 = w.end_samp / ofs * 1e6;
  Json feats = Json::object();
  feats["modulation"] = "lora";
  feats["sf"] = sf;
  feats["bw_hz"] = bw;
  feats["bw_mhz"] = bw / 1e6;
  feats["iq_inverted"] = w.pol < 0;
  feats["preamble_syms"] = w.n_pre;
  feats["payload_syms"] = w.n_pay;
  feats["sync_syms"] = w.sync_vals;
  feats["symbol_us"] = static_cast<double>(N) / bw * 1e6;
  feats["burst_us"] = t1 - t0;
  feats["f_carrier_hz"] = f;
  feats["cfo_hz"] = cfo;
  feats["cfo_ppm"] = cfo.is_null() ? Json(nullptr) : Json(cfo.get<double>() / f * 1e6);
  feats["clock_ppm"] = clock_ppm;
  feats["clock_ppm_err"] = clock_err;
  feats["elrs_channel"] = (nom && !cfo.is_null()) ? Json::array({nom->domain, nom->index}) : Json(nullptr);
  feats["elrs_modes"] = full.empty() ? matches_all : full;
  {
    Json mm = Json::array();
    for (const auto& m : ml) mm.push_back(m.mode);
    feats["mlrs_modes"] = mm;
  }
  feats["preamble_ratio"] = w.pre_ratio;
  feats["sfd_ratio"] = w.sfd_ratio;
  feats["z_pre"] = w.z;
  feats["z_sfd"] = w.z_sfd;
  if (!dec.is_null()) {
    feats["sync_word"] = dec.contains("sync_word") ? dec["sync_word"] : Json(nullptr);
    feats["coding_rate"] = dec["cr"];
    feats["ldro"] = dec["ldro"];
  }
  Emission e;
  e.family = family;
  e.label = label;
  e.start_us = t0;
  e.end_us = t1;
  e.f_lo_mhz = (f - bw / 2) / 1e6;
  e.f_hi_mhz = (f + bw / 2) / 1e6;
  e.center_mhz = f / 1e6;
  e.snr_db = w.snr_db;
  e.confidence = conf;
  e.features = std::move(feats);
  e.decoded = dec;
  return e;
}

Emission lora_preamble_emission(const Partial& p, double fc) {
  const int sf = p.sf;
  const double bw = p.bw, ofs = p.ofs, f = p.f_sb;
  const double t0 = static_cast<double>(p.pre_start) / ofs * 1e6, t1 = static_cast<double>(p.pre_end) / ofs * 1e6;
  const auto nom = nominal_channel(f);
  bool elrs = false;
  if (nom && std::fabs(f - nom->f) < bw / 4)
    for (const auto& m : elrs_lora())
      if (m.sf == sf && std::fabs(m.bw - bw) < 1e3) elrs = true;
  Json feats = Json::object();
  feats["modulation"] = "lora";
  feats["sf"] = sf;
  feats["bw_hz"] = bw;
  feats["bw_mhz"] = bw / 1e6;
  feats["iq_inverted"] = p.pol < 0;
  feats["preamble_syms"] = p.n_pre_est;
  feats["symbol_us"] = static_cast<double>(p.N) / bw * 1e6;
  feats["burst_us"] = t1 - t0;
  feats["sfd_found"] = false;
  feats["cfo_hz"] = nullptr;
  feats["cfo_ppm"] = nullptr;
  feats["clock_ppm"] = nullptr;
  feats["f_uncertainty_hz"] = bw / 4;
  feats["stage1_ratio"] = p.ratio1;
  feats["elrs_channel"] = elrs ? Json::array({nom->domain, nom->index}) : Json(nullptr);
  Emission e;
  e.family = elrs ? "rc_fhss" : "lora";
  e.label = std::string(elrs ? "elrs_like_" : "") + "lora_preamble_sf" + std::to_string(sf) + "_bw" + kfmt(bw) + "k";
  e.start_us = t0;
  e.end_us = t1;
  e.f_lo_mhz = (f - bw / 2) / 1e6;
  e.f_hi_mhz = (f + bw / 2) / 1e6;
  e.center_mhz = f / 1e6;
  e.snr_db = std::nullopt;
  e.confidence = 0.5;
  e.features = std::move(feats);
  return e;
}

std::vector<Emission> dedup(std::vector<Emission> ems, double f_tol_frac) {
  // a decoded packet beats a copy of itself walked one symbol off
  auto key = [](const Emission& e) {
    const double pr = e.features.contains("preamble_ratio") && e.features["preamble_ratio"].is_number()
                          ? e.features["preamble_ratio"].get<double>()
                          : 0.0;
    const bool dec = e.decoded.is_object() && !e.decoded.empty();
    return std::make_tuple(-static_cast<int>(truthy_crc(e.decoded)), -static_cast<int>(dec), -(pr != 0 ? pr : 0.0));
  };
  std::stable_sort(ems.begin(), ems.end(), [&](const Emission& a, const Emission& b) { return key(a) < key(b); });
  std::vector<Emission> keep;
  for (auto& e : ems) {
    const double bw = e.f_hi_mhz - e.f_lo_mhz;
    const bool no_sfd = e.features.contains("sfd_found") && e.features["sfd_found"].is_boolean() &&
                        !e.features["sfd_found"].get<bool>();
    const double tol = no_sfd ? 1.0 : f_tol_frac;
    bool dup = false;
    for (const auto& k : keep)
      if (std::fabs(k.center_mhz - e.center_mhz) < tol * bw &&
          std::min(k.end_us, e.end_us) - std::max(k.start_us, e.start_us) > 0.3 * (e.end_us - e.start_us)) {
        dup = true;
        break;
      }
    if (!dup) keep.push_back(std::move(e));
  }
  std::stable_sort(keep.begin(), keep.end(), [](const Emission& a, const Emission& b) { return a.start_us < b.start_us; });
  return keep;
}

}  // namespace detail

std::vector<Emission> detect_lora(const cfvec& X, double fs, double fc, const std::pair<double, double>* f_range) {
  using namespace detail;
  std::vector<Emission> out;
  for (const auto& [bw, sfs] : lora_bws_for(fc)) {
    std::vector<Cand> cands = nms_candidates(lora_scan(X, fs, bw, sfs, 8, 24.0, 48, 0.45, f_range), bw);
    // refine sub-band by sub-band (parallel over sub-bands), sharing each one's slices
    std::map<double, std::vector<size_t>> by_sb;
    for (size_t i = 0; i < cands.size(); ++i) by_sb[cands[i].f_sb].push_back(i);
    std::vector<std::vector<size_t>> groups;
    for (auto& [f, idx] : by_sb) groups.push_back(std::move(idx));
    const auto parts = pmap(groups, [&](const std::vector<size_t>& idx) {
      SliceCache cache;
      std::vector<std::optional<Refined>> r;
      for (size_t i : idx) r.push_back(lora_refine(X, fs, fc, cands[i], 3, 5.0, &cache));
      return r;
    });
    std::vector<std::optional<Refined>> refined(cands.size());
    for (size_t g = 0; g < groups.size(); ++g)
      for (size_t k = 0; k < groups[g].size(); ++k) refined[groups[g][k]] = std::move(parts[g][k]);
    for (size_t i = 0; i < cands.size(); ++i) {
      const auto& w = refined[i];
      if (!w) continue;
      const Cand& c = cands[i];
      // a sync word's two symbols are its nibbles x 8: anything else is not LoRa
      bool sync_ok = true;
      if (!w->partial)
        for (int64_t v : w->walk.sync_vals) sync_ok = sync_ok && std::min(v % 8, 8 - v % 8) <= 1;
      const bool ok = !w->partial && w->walk.z >= 12.0 && w->walk.z_sfd >= 4.0 && w->walk.n_pre >= 4 &&
                      w->walk.n_pay >= 4 && sync_ok;
      // preamble-only: a real preamble is followed by its SFD, so either the SFD lies past the end of the
      // capture, or it was seen (z_sfd) and only the payload is cut short
      const bool sfd_ok = !w->prem.sfd_window_in_capture || (!w->partial && w->walk.z_sfd >= 4.0);
      const bool ok_pre = !ok && c.ratio >= 60.0 && 5 <= w->prem.n_pre_est && w->prem.n_pre_est <= 40 && sfd_ok &&
                          sync_ok;
      if (ok) out.push_back(lora_emission(w->walk, fc, fs));
      else if (ok_pre) out.push_back(lora_preamble_emission(w->prem, fc));
    }
  }
  return dedup(std::move(out));
}

}  // namespace si::rx::rc_link
