#include "core/ofdm15k.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/signal.hpp"
#include "core/spectral.hpp"

namespace si::ofdm15k {

namespace {

constexpr double kPi = 3.141592653589793;

// int(x) of a Python float: truncation toward zero
int64_t py_int(double x) { return static_cast<int64_t>(x); }

dvec abs2(const cd* y, size_t n) {
  dvec out(n);
  for (size_t i = 0; i < n; ++i) {
    const double a = np::cabs(y[i]);
    out[i] = a * a;
  }
  return out;
}

double median_all(const dmat& G) { return np::median(G.data); }

}  // namespace

const std::vector<Profile>& profiles() {
  static const std::vector<Profile> p = {
      {"1.4MHz", 1.92e6, 128, 9, 1.08e6, 0.25e6},
      {"10MHz", 15.36e6, 1024, 72, 9.0e6, 1.0e6},
      {"20MHz", 30.72e6, 2048, 144, 18.0e6, 2.0e6},
  };
  return p;
}

Scan ofdm15k_scan(const Slicer& ch, const std::vector<Profile>& profs, double win_s, double k_sigma,
                  double usable_half) {
  Scan out;
  for (const Profile& pr : profs) {
    const size_t sym = pr.N + pr.cp;
    const size_t W = static_cast<size_t>(py_int(win_s * pr.ofs)), hop = W / 2;
    const dvec offs = np::arange(-usable_half + pr.occ / 2, usable_half - pr.occ / 2 + 1, pr.step);
    // one CP-metric row per scan centre (independent: parallel, ordered)
    const std::vector<dvec> rows = pmap(offs, [&](double off) {
      const cdvec y = ch(off, pr.occ / 2 * 1.05, pr.ofs);
      return folded_scores(cp_metric(y, pr.N, pr.cp).m, sym, W, hop);
    });
    const size_t nc = rows.empty() ? 0 : rows[0].size();
    dmat G(rows.size(), nc);
    for (size_t r = 0; r < rows.size(); ++r) {
      if (rows[r].size() != nc) throw std::runtime_error("ofdm15k_scan: ragged score grid");
      std::copy(rows[r].begin(), rows[r].end(), G.row(r));
    }
    const double med = median_all(G);
    dvec dev(G.size());
    for (size_t i = 0; i < G.size(); ++i) dev[i] = std::fabs(G.data[i] - med);
    const double mad = np::median(dev) * 1.4826;
    const double thr = med + k_sigma * mad;
    bmat mask(G.rows, G.cols, 0);
    for (size_t i = 0; i < G.size(); ++i) mask.data[i] = G.data[i] > thr;
    const auto lab = ndi::label(mask, ndi::ones(3, 3));
    const auto objs = ndi::find_objects(lab.lab, lab.n);
    for (int i = 0; i < lab.n; ++i) {
      const auto& o = objs[static_cast<size_t>(i)];
      // argmax over the box with other labels at -inf (first maximum, row-major)
      double best = -std::numeric_limits<double>::infinity();
      size_t ci = o.r0, wi = o.c0;
      bool first = true;
      for (size_t r = o.r0; r < o.r1; ++r)
        for (size_t c = o.c0; c < o.c1; ++c) {
          const double v = lab.lab(r, c) == i + 1 ? G(r, c) : -std::numeric_limits<double>::infinity();
          if (first || v > best) {
            best = v;
            ci = r;
            wi = c;
            first = false;
          }
        }
      Hit h;
      h.profile = pr.name;
      h.ofs = pr.ofs;
      h.N = pr.N;
      h.cp = pr.cp;
      h.occ = pr.occ;
      h.f_off = offs[ci];
      h.t0_s = static_cast<double>(wi * hop) / pr.ofs;
      h.t1_s = static_cast<double>(wi * hop + W) / pr.ofs;
      h.score = G(ci, wi);
      h.thr = thr;
      out.hits.push_back(h);
    }
    out.stats[pr.name] = Stats{med, mad};
  }
  return out;
}

std::optional<Burst> refine_ofdm15k(const Slicer& ch, const Hit& hit, double pad_s) {
  const double ofs = hit.ofs, occ = hit.occ;
  const size_t N = hit.N, cp = hit.cp;
  cdvec y = ch(hit.f_off, occ / 2 * 1.05, ofs);
  const int64_t ny = static_cast<int64_t>(y.size());
  const int64_t a = std::max<int64_t>(0, py_int((hit.t0_s - pad_s) * ofs));
  const int64_t b = std::min<int64_t>(ny, py_int((hit.t1_s + pad_s) * ofs));
  const size_t nseg = b > a ? static_cast<size_t>(b - a) : 0;
  const int64_t k = std::max<int64_t>(py_int(20e-6 * ofs), 4);
  const dvec ones_k(static_cast<size_t>(k), 1.0 / static_cast<double>(k));
  const dvec env = np::convolve(abs2(y.data() + a, nseg), ones_k, np::ConvMode::Same);
  if (env.empty()) throw std::runtime_error("refine_ofdm15k: empty burst segment");
  double pkf = ((hit.t0_s + hit.t1_s) / 2) * ofs - static_cast<double>(a);
  pkf = std::min(std::max(pkf, 0.0), static_cast<double>(env.size() - 1));
  const int64_t pk = py_int(pkf);
  const double hi = np::percentile(env.data(), env.size(), 90.0), lo = np::percentile(env.data(), env.size(), 10.0);
  const double gate = std::sqrt(hi * lo);
  int64_t s = pk, e = pk;
  while (s > 0 && env[static_cast<size_t>(s - 1)] > gate) --s;
  while (e < static_cast<int64_t>(env.size()) - 1 && env[static_cast<size_t>(e + 1)] > gate) ++e;
  const int64_t s0 = a + s, e0 = a + e;
  if (e0 - s0 < static_cast<int64_t>(2 * (N + cp))) return std::nullopt;  // shorter than two symbols

  // occupied band (-10 dB from the in-band level) at 3x the profile width: bandwidth and
  // coarse centre
  const cdvec wide = ch(hit.f_off, occ * 1.5, ofs * 2);
  const int64_t ws = s0 * 2, we = e0 * 2;
  const size_t nf = 4 * N;
  const dvec win = np::hanning(nf);
  std::vector<int64_t> starts;
  for (int64_t j = ws; j < std::max(ws + 1, we - static_cast<int64_t>(nf)); j += static_cast<int64_t>(nf / 2))
    starts.push_back(j);
  dmat frames(starts.size(), nf);
  cdvec seg(nf);
  for (size_t r = 0; r < starts.size(); ++r) {
    const int64_t j = starts[r];
    if (j + static_cast<int64_t>(nf) > static_cast<int64_t>(wide.size()))
      throw std::runtime_error("refine_ofdm15k: frame beyond the slice (numpy broadcast error)");
    for (size_t i = 0; i < nf; ++i) {
      const cd v = wide[static_cast<size_t>(j) + i];
      seg[i] = cd(v.real() * win[i], v.imag() * win[i]);
    }
    fft::c2c(seg.data(), seg.data(), nf, true);
    double* row = frames.row(r);
    for (size_t i = 0; i < nf; ++i) {
      const double m = np::cabs(seg[i]);
      row[(i + nf / 2) % nf] = m * m;  // fftshift
    }
  }
  double bw = occ, shift = 0.0;
  if (!starts.empty()) {
    dvec ps = ndi::uniform_filter1d(np::median_axis0(frames), 9);
    for (auto& v : ps) v = 10 * std::log10(std::max(v, 1e-20));
    dvec ff(nf);
    std::vector<uint8_t> inb(nf);
    dvec ps_in;
    for (size_t i = 0; i < nf; ++i) {
      ff[i] = static_cast<double>(static_cast<int64_t>(i) - static_cast<int64_t>(nf / 2)) / static_cast<double>(nf) *
              ofs * 2;
      inb[i] = std::fabs(ff[i]) < occ * 0.3;
      if (inb[i]) ps_in.push_back(ps[i]);
    }
    const double level = np::median(ps_in) - 10;
    dvec masked(nf);
    for (size_t i = 0; i < nf; ++i) masked[i] = inb[i] ? ps[i] : -std::numeric_limits<double>::infinity();
    size_t L = np::argmax(masked), R = L;  // strongest in-band bin
    while (L > 0 && ps[L - 1] > level) --L;
    while (R < nf - 1 && ps[R + 1] > level) ++R;
    bw = ff[R] - ff[L];
    shift = (ff[R] + ff[L]) / 2;
  }
  double f_off = hit.f_off + shift;
  // fine CFO from the CP correlation phase over the burst, at the corrected centre
  y = ch(f_off, occ / 2 * 1.05, ofs);
  const CpMetric cm = cp_metric(y.data() + s0, static_cast<size_t>(e0 - s0), N, cp);
  const double p90 = np::percentile(cm.m.data(), cm.m.size(), 90.0);
  const cdvec cc = np::convolve(cm.c, cdvec(cp, cd(1.0, 0.0)), np::ConvMode::Valid);
  cdvec sel;
  for (size_t i = 0; i < cm.m.size(); ++i)
    if (cm.m[i] > p90) sel.push_back(cc[i]);
  const cd tot = np::sum(sel.data(), sel.size());
  const double ph = std::atan2(tot.imag(), tot.real());
  const double cfo = -ph / (2 * kPi * static_cast<double>(N) / ofs);
  f_off += cfo;
  if (bw > 1.5 * occ) return std::nullopt;  // band edges not found
  // verification: a train of CP peaks one symbol period apart, fitted to a uniform grid
  const size_t sym = N + cp;
  sig::PeakOpts po;
  po.height_min = 0.5;
  po.distance = static_cast<double>(py_int(0.8 * static_cast<double>(sym)));
  const auto pks = sig::find_peaks(cm.m, po).peaks;
  if (pks.size() < 3) return std::nullopt;
  dvec kk(pks.size()), pp(pks.size());
  std::set<double> uniq;
  for (size_t i = 0; i < pks.size(); ++i) {
    kk[i] = std::nearbyint(static_cast<double>(pks[i] - pks[0]) / static_cast<double>(sym));
    pp[i] = static_cast<double>(pks[i]);
    uniq.insert(kk[i]);
  }
  if (uniq.size() < 3) return std::nullopt;
  const dvec fit = np::polyfit(kk, pp, 1);
  const double period = fit[0], a0 = fit[1];
  dvec res(pks.size());
  for (size_t i = 0; i < pks.size(); ++i) {
    const double d = pp[i] - (a0 + period * kk[i]);
    res[i] = d * d;
  }
  const double scatter = std::sqrt(np::mean(res.data(), res.size()));
  const double symd = static_cast<double>(sym);
  if (std::fabs(period - symd) > 0.03 * symd || scatter > 0.12 * symd) return std::nullopt;
  Burst out;
  out.start_s = static_cast<double>(s0) / ofs;
  out.end_s = static_cast<double>(e0) / ofs;
  out.f_off = f_off;
  out.cfo_hz = cfo;
  out.bw_hz = bw;
  out.nsym = py_round_int(static_cast<double>(e0 - s0) / static_cast<double>(N + cp));
  out.score = hit.score;
  out.profile = hit.profile;
  out.occ = occ;
  out.cp_peaks = static_cast<int64_t>(pks.size());
  out.period_us = py_round(period / ofs * 1e6, 2);
  out.cp_scatter_us = py_round(scatter / ofs * 1e6, 2);
  return out;
}

std::vector<Burst> merge_ofdm15k(const std::vector<Burst>& bursts) {
  std::map<std::string, size_t> order;
  for (size_t i = 0; i < profiles().size(); ++i) order[profiles()[i].name] = i;
  std::vector<Burst> sorted = bursts;
  std::stable_sort(sorted.begin(), sorted.end(), [&](const Burst& x, const Burst& y) {
    const size_t ox = order.at(x.profile), oy = order.at(y.profile);
    if (ox != oy) return ox < oy;
    return -x.score < -y.score;
  });
  std::vector<Burst> keep;
  for (const Burst& b : sorted) {
    bool dup = false;
    for (const Burst& k : keep) {
      const double ov = std::min(b.end_s, k.end_s) - std::max(b.start_s, k.start_s);
      if (ov <= 0.5 * std::min(b.end_s - b.start_s, k.end_s - k.start_s)) continue;
      if (std::fabs(b.f_off - k.f_off) <= (b.occ + k.occ) / 2) {
        dup = true;
        break;
      }
    }
    if (!dup) keep.push_back(b);
  }
  return keep;
}

std::optional<double> cycle(const dvec& starts_s) {
  dvec s = starts_s;
  std::sort(s.begin(), s.end());
  const dvec iv = np::diff(s);
  if (iv.size() >= 3) {
    dvec pair(iv.size() - 1);
    for (size_t i = 0; i + 1 < iv.size(); ++i) pair[i] = iv[i] + iv[i + 1];
    if (np::stddev(pair) < 0.15 * np::mean(pair) && np::stddev(iv) > 0.1 * np::mean(iv)) return np::mean(pair);
  }
  if (iv.size() >= 2 && np::stddev(iv) < 0.1 * np::mean(iv)) return np::mean(iv);
  return std::nullopt;
}

std::vector<Burst> track_hops(const Slicer& ch, const std::vector<Burst>& bursts,
                              const std::map<std::string, Stats>& stats, double dur_s, double k_sigma,
                              double usable_half) {
  const Profile& pr = profiles()[0];
  std::vector<Burst> narrow;
  for (const Burst& b : bursts)
    if (b.profile == pr.name) narrow.push_back(b);
  if (narrow.size() < 4) return {};
  dvec starts(narrow.size());
  for (size_t i = 0; i < narrow.size(); ++i) starts[i] = narrow[i].start_s;
  const auto T = cycle(starts);
  if (!T) return {};
  std::set<double> cs;  // round(np.float64, 5): numpy rounding
  for (double t : starts)
    for (int k = -4; k <= 4; ++k)
      if (k) cs.insert(np::round(t + k * *T, 5));
  std::vector<double> merged;
  for (double t : cs) {
    if (!(0 <= t && t <= dur_s - 0.5e-3)) continue;
    double mn = std::numeric_limits<double>::infinity();
    for (double s : starts) mn = std::min(mn, std::fabs(s - t));
    if (!(mn > 1e-3)) continue;
    if (merged.empty() || t - merged.back() > 0.3e-3) merged.push_back(t);  // coinciding predictions
  }
  if (merged.empty()) return {};
  const Stats& st = stats.at(pr.name);
  const double thr = st.median + k_sigma * st.sigma;
  const size_t sym = pr.N + pr.cp, W = static_cast<size_t>(py_int(0.5e-3 * pr.ofs));
  const dvec offs = np::arange(-usable_half + pr.occ / 2, usable_half - pr.occ / 2 + 1, pr.step / 2);
  dvec dur(narrow.size());
  for (size_t i = 0; i < narrow.size(); ++i) dur[i] = narrow[i].end_s - narrow[i].start_s;
  const double dur0 = np::median(dur);
  const dvec dts = np::arange(-0.25e-3, 0.11e-3, 0.05e-3);  // hop timing jitter

  // best (score, t0) per (prediction, centre). Python slices every centre up front; here
  // each centre is sliced once and scored for all predictions (same values, bounded memory).
  struct Best {
    double score = -1.0, t0 = 0.0;
  };
  const std::vector<std::vector<Best>> grid = pmap(offs, [&](double off) {
    const cdvec y = ch(off, pr.occ / 2 * 1.05, pr.ofs);
    std::vector<Best> row(merged.size());
    for (size_t ti = 0; ti < merged.size(); ++ti) {
      Best best;
      for (double dt : dts) {
        const int64_t a = py_int((merged[ti] + dt) * pr.ofs);
        if (a < 0 || a + static_cast<int64_t>(W + pr.N + pr.cp) >= static_cast<int64_t>(y.size())) continue;
        const CpMetric cm = cp_metric(y.data() + a, W + pr.N + pr.cp, pr.N, pr.cp);
        const dvec sc = folded_scores(cm.m, sym, W, W);
        if (!sc.empty() && sc[0] > best.score) best = Best{sc[0], merged[ti] + dt};
      }
      row[ti] = best;
    }
    return row;
  });

  std::vector<Burst> found;
  for (size_t ti = 0; ti < merged.size(); ++ti) {
    struct Cand {
      double score, off, t0;
    };
    std::vector<Cand> cands;
    for (size_t oi = 0; oi < offs.size(); ++oi)
      if (grid[oi][ti].score > thr) cands.push_back({grid[oi][ti].score, offs[oi], grid[oi][ti].t0});
    // sorted(cands, reverse=True): tuples compared element-wise, descending
    std::stable_sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) {
      if (x.score != y.score) return x.score > y.score;
      if (x.off != y.off) return x.off > y.off;
      return x.t0 > y.t0;
    });
    if (cands.size() > 8) cands.resize(8);
    for (const Cand& c : cands) {
      Hit h;
      h.profile = pr.name;
      h.ofs = pr.ofs;
      h.N = pr.N;
      h.cp = pr.cp;
      h.occ = pr.occ;
      h.f_off = c.off;
      h.t0_s = c.t0;
      h.t1_s = c.t0 + 0.5e-3;
      h.score = c.score;
      h.thr = thr;
      auto b = refine_ofdm15k(ch, h);
      if (!b || b->bw_hz > 2.5e6 || std::fabs((b->end_s - b->start_s) - dur0) > 0.25 * dur0) continue;
      b->tracked = true;
      found.push_back(*b);
      break;
    }
  }
  return found;
}

std::pair<std::string, std::string> label_ofdm15k(const Burst& b) {
  const double dur_us = (b.end_s - b.start_s) * 1e6;
  const double bw = b.bw_hz;
  if (bw < 2.5e6 && 300 <= dur_us && dur_us <= 1200)
    return {"ocusync_uplink_like", "LTE-like 15 kHz OFDM, ~1.4 MHz, slot-length burst (RC->drone hop)"};
  if (7e6 <= bw && bw <= 11e6 && 500 <= dur_us && dur_us <= 800)
    return {"droneid_candidate", "LTE-like 15 kHz OFDM, ~10 MHz, ~9-symbol burst (DJI DroneID shape)"};
  if (bw >= 7e6) return {"ocusync_video_like", "LTE-like 15 kHz OFDM, wideband (drone->RC video downlink shape)"};
  return {"ofdm15k_other", "15 kHz-SCS OFDM, unrecognised shape"};
}

}  // namespace si::ofdm15k
