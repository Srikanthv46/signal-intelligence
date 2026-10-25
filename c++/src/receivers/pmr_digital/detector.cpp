// python: receivers/pmr_digital/detector.py
#include "receivers/pmr_digital/detector.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"
#include "receivers/pmr_digital/fsk4.hpp"
#include "receivers/pmr_digital/tetra.hpp"

namespace si::rx::pmr_digital {

namespace {

constexpr double kRes = 1000.0;
constexpr double kPeakDb = 8.0;
constexpr double kMaxWidth = 30e3;
constexpr double kTetraLine = 6.0;      // worth demodulating
constexpr double kTetraLikeLine = 9.0;  // enough on its own (no training sequence) for a "-like" carrier
constexpr double kMaxEnvCv = 0.35;      // 4FSK is constant-envelope
constexpr double kFskFs = 48000.0;
constexpr size_t kMaxCand = 24;

bool in(double f, const std::vector<std::pair<double, double>>& r) {
  for (const auto& [a, b] : r)
    if (a <= f && f <= b) return true;
  return false;
}

struct Cand {
  double f_hz, peak_db, width_hz, floor;
};

// python candidates(ch, band)
std::vector<Cand> candidates(const Channelizer& ch, std::pair<double, double> band) {
  const double fs = ch.fs(), fc = ch.fc();
  const size_t L = ch.L();
  const double lo = band.first, hi = band.second;
  std::vector<std::pair<double, double>> rng;
  for (const auto* R : {&tetra_ranges(), &fsk_ranges()})
    for (const auto& r : *R)
      if (r.first < hi && r.second > lo) rng.push_back(r);
  if (rng.empty()) return {};
  const size_t g = static_cast<size_t>(std::max<int64_t>(1, py_round_int(kRes * static_cast<double>(L) / fs)));
  const size_t start = (L + 1) / 2;  // np.argsort(fftfreq): negative frequencies first
  fvec P(L);
  dvec f(L);
  for (size_t i = 0; i < L; ++i) {
    const size_t j = (start + i) % L;
    P[i] = np::cabs2(ch.X()[j]);
    f[i] = ch.f()[j];
  }
  const size_t n = L / g;
  fvec Pg(n);
  dvec fa(n);
  std::vector<uint8_t> inb(n);
  fvec inner;
  for (size_t j = 0; j < n; ++j) {
    Pg[j] = np::mean(P.data() + j * g, g);
    fa[j] = fc + np::mean(f.data() + j * g, g);
    inb[j] = fa[j] >= lo && fa[j] <= hi;
    if (inb[j]) inner.push_back(Pg[j]);
  }
  if (inner.size() < 20) return {};
  const double floor = static_cast<double>(np::median(inner));
  const float floorf = static_cast<float>(floor);
  fvec Ldb(n), Lw(n);
  for (size_t j = 0; j < n; ++j) {
    Ldb[j] = 10.0f * std::log10(Pg[j] / floorf + 1e-30f);
    Lw[j] = inb[j] ? Ldb[j] : 0.0f;
  }
  sig::PeakOpts o;
  o.height_min = kPeakDb;
  o.prominence = true;
  o.prominence_min = kPeakDb;
  o.width = true;
  o.width_max = kMaxWidth / kRes;
  o.rel_height = 0.5;
  o.distance = std::max(1.0, std::floor(6e3 / kRes));
  const sig::Peaks pk = sig::find_peaks(Lw, o);
  std::vector<Cand> out;
  for (size_t q = 0; q < pk.peaks.size(); ++q) {
    const int64_t i = pk.peaks[q];
    const double f_abs = fa[i];
    if (!in(f_abs, rng)) continue;
    const int64_t a = std::max<int64_t>(0, i - 12), b = std::min<int64_t>(static_cast<int64_t>(n), i + 13);
    fvec p;
    dvec fp;
    for (int64_t k = a; k < b; ++k) {
      p.push_back(std::max(Pg[k] - floorf, 0.0f));
      fp.push_back(fa[k]);
    }
    const float ps = np::sum(p);
    double fcen = f_abs;
    if (ps > 0) {
      dvec w(p.size());
      for (size_t k = 0; k < p.size(); ++k) w[k] = fp[k] * static_cast<double>(p[k]);
      fcen = np::sum(w) / static_cast<double>(ps);
    }
    out.push_back({fcen, static_cast<double>(Ldb[i]), pk.widths[q] * kRes, floor / (fs * static_cast<double>(L))});
  }
  std::stable_sort(out.begin(), out.end(), [](const Cand& x, const Cand& y) { return -x.peak_db < -y.peak_db; });
  std::vector<Cand> keep;
  for (const Cand& c : out) {
    bool ok = true;
    for (const Cand& k : keep) ok = ok && std::fabs(c.f_hz - k.f_hz) >= 12.5e3;
    if (ok) keep.push_back(c);
  }
  if (keep.size() > kMaxCand) keep.resize(kMaxCand);
  return keep;
}

// python _snr(y, floor_density, bw)
Json snr_of(const cfvec& y, double floor_density, double bw) {
  const double noise = floor_density * bw;
  fvec a2(y.size());
  for (size_t i = 0; i < y.size(); ++i) a2[i] = np::cabs2(y[i]);
  const double p = static_cast<double>(np::mean(a2));
  if (!(noise > 0 && p > noise)) return nullptr;
  return np::round(10 * std::log10(std::max(p - noise, 1e-30) / noise), 1);
}

std::optional<Emission> tetra_carrier(const Channelizer& ch, const Cand& c, double dur) {
  const auto sl = ch.slice(c.f_hz, 12.5e3, tetra::kSr * tetra::kSps);
  const double line = tetra::symbol_line(sl.y, sl.fs);
  if (line < kTetraLine) return std::nullopt;
  struct Best {
    int64_t n_bsch = -1, n_sync_slot = 0, norm_ok = 0;
    double err = 0, foff = 0;
    std::vector<int64_t> sync, normal;
    std::vector<Json> bsch;
  } best;
  bool have = false;
  for (int q : {0, -1, 1}) {  // quarter-turn ambiguity of the carrier estimate
    const tetra::Demod dm = tetra::demod(sl.y, q);
    const tetra::Training tr = tetra::find_training(dm.bits);
    std::vector<Json> bsch;
    for (int64_t h : tr.sync)
      if (h >= tetra::kBschBits)
        if (auto pdu = tetra::decode_bsch(dm.soft.data() + (h - tetra::kBschBits))) {
          Json b = *pdu;
          b["t_s"] = py_round(static_cast<double>(h) / 2 / tetra::kSr, 6);
          bsch.push_back(std::move(b));
        }
    const int64_t n_norm = tetra::slot_consistent(tr.normal);
    const int64_t norm_ok =
        static_cast<double>(n_norm) >= std::max(3.0, 0.5 * static_cast<double>(tr.normal.size())) ? n_norm : 0;
    const int64_t n_sync_slot = tetra::slot_consistent(tr.sync);
    const auto score = std::make_tuple(static_cast<int64_t>(bsch.size()), n_sync_slot, norm_ok);
    if (!have || score > std::make_tuple(best.n_bsch, best.n_sync_slot, best.norm_ok)) {
      have = true;
      best = {static_cast<int64_t>(bsch.size()), n_sync_slot, norm_ok, dm.err, dm.foff, tr.sync, tr.normal, bsch};
    }
  }
  const int64_t n_slot = std::max(best.n_sync_slot, best.norm_ok);
  const bool framed = !best.bsch.empty() || best.n_sync_slot >= 2 || best.norm_ok >= 3;
  const double f_car = c.f_hz + best.foff;
  Json F = Json::object();
  F["carrier_hz"] = py_round_int(f_car);
  F["symbol_rate_line"] = py_round(line, 1);
  F["phase_error_rad"] = py_round(best.err, 3);
  F["n_sync_bursts"] = best.sync.size();
  F["n_normal_training"] = best.normal.size();
  F["slot_consistent_hits"] = n_slot;
  F["n_bsch_crc_ok"] = best.bsch.size();
  F["modulation"] = "pi/4-DQPSK 18 ksym/s";
  F["snr_db"] = snr_of(sl.y, c.floor, 25e3);
  std::string label;
  double conf;
  if (!best.bsch.empty()) {
    std::vector<std::pair<std::string, int>> ids;  // identity -> count, in order of first appearance
    for (const Json& b : best.bsch) {
      const std::string k = b["mcc"].dump() + "/" + b["mnc"].dump() + "/" + b["colour_code"].dump();
      auto it = std::find_if(ids.begin(), ids.end(), [&](const auto& x) { return x.first == k; });
      if (it == ids.end()) ids.emplace_back(k, 1);
      else ++it->second;
    }
    size_t ti = 0;
    for (size_t i = 1; i < ids.size(); ++i)
      if (ids[i].second > ids[ti].second) ti = i;
    const Json* first = nullptr;
    for (const Json& b : best.bsch)
      if (b["mcc"].dump() + "/" + b["mnc"].dump() + "/" + b["colour_code"].dump() == ids[ti].first) {
        first = &b;
        break;
      }
    for (const char* k : {"mcc", "mnc", "colour_code", "system_code", "sharing_mode", "ts_reserved_frames",
                          "u_plane_dtx", "frame18_extension", "cell_service_level", "late_entry",
                          "neighbour_cell_broadcast"})
      F[k] = (*first)[k];
    Json lst = Json::array();
    for (size_t i = 0; i < best.bsch.size() && i < 20; ++i) {
      const Json& b = best.bsch[i];
      lst.push_back({{"t_s", b["t_s"]}, {"timeslot", b["timeslot"]}, {"frame", b["frame"]},
                     {"multiframe", b["multiframe"]}});
    }
    F["bsch"] = lst;
    F["bsch_identities"] = ids.size();
    label = "tetra_downlink_bsch";
    conf = 0.97;
  } else if (framed) {
    label = "tetra_downlink";
    conf = 0.85;
  } else if (line >= kTetraLikeLine) {
    label = "tetra_like_carrier";
    conf = 0.55;
  } else {
    return std::nullopt;
  }
  Emission e;
  e.family = "tetra";
  e.label = label;
  e.start_us = 0.0;
  e.end_us = dur * 1e6;
  e.f_lo_mhz = (f_car - 12.5e3) / 1e6;
  e.f_hi_mhz = (f_car + 12.5e3) / 1e6;
  e.center_mhz = f_car / 1e6;
  if (!F["snr_db"].is_null()) e.snr_db = F["snr_db"].get<double>();
  e.confidence = conf;
  if (!best.bsch.empty()) e.decoded = {{"crc_ok", true}};
  e.features = F;
  return e;
}

Emission dmr_emission(const Cand& c, double dur, const std::vector<Json>& bursts, const Json& base) {
  Json kinds = Json::object();
  for (const Json& b : bursts) {
    const std::string k = b["sync"].get<std::string>();
    kinds[k] = kinds.contains(k) ? kinds[k].get<int>() + 1 : 1;
  }
  std::vector<int> cc;
  Json dts = Json::object();
  int64_t crc_checked = 0, crc_ok = 0;
  std::vector<const Json*> csbk;
  for (const Json& b : bursts) {
    if (b.contains("colour_code")) cc.push_back(b["colour_code"].get<int>());
    if (b.contains("data_type")) {
      const std::string k = b["data_type"].get<std::string>();
      dts[k] = dts.contains(k) ? dts[k].get<int>() + 1 : 1;
    }
    if (b.contains("crc_ok")) {
      ++crc_checked;
      crc_ok += b["crc_ok"].get<bool>();
      if (b["crc_ok"].get<bool>() && b.contains("csbko")) csbk.push_back(&b);
    }
  }
  int bs = 0, ms = 0, direct = 0, voice = 0;
  for (const auto& [k, v] : kinds.items()) {
    const int n = v.get<int>();
    if (k.rfind("bs_", 0) == 0) bs += n;
    if (k.rfind("ms_", 0) == 0) ms += n;
    if (k.rfind("direct", 0) == 0) direct += n;
    if (k.size() >= 5 && k.compare(k.size() - 5, 5, "voice") == 0) voice += n;
  }
  const std::string role = bs >= std::max(ms, direct) ? "bs" : (ms >= direct ? "ms" : "direct");
  Json F = base;
  F["modulation"] = "4FSK 4800 sym/s";
  F["channel_bw_khz"] = 12.5;
  F["n_bursts"] = bursts.size();
  F["sync_types"] = kinds;
  F["source"] = role == "bs" ? "base station / repeater" : (role == "ms" ? "mobile station" : "direct mode (TDMA)");
  F["voice_bursts"] = voice;
  F["data_types"] = dts;
  Json lst = Json::array();
  for (size_t i = 0; i < bursts.size() && i < 12; ++i) lst.push_back(bursts[i]);
  F["bursts"] = lst;
  if (!cc.empty()) {
    int top = -1, topn = 0;
    for (int v = 0; v < 16; ++v) {  // most frequent; ties: the lowest
      const int n = static_cast<int>(std::count(cc.begin(), cc.end(), v));
      if (n > topn) {
        topn = n;
        top = v;
      }
    }
    F["colour_code"] = top;
    F["colour_code_votes"] = topn;
  }
  if (crc_checked) {
    F["crc_checked"] = crc_checked;
    F["crc_ok"] = crc_ok;
  }
  if (!csbk.empty()) {
    Json cl = Json::array();
    for (size_t i = 0; i < csbk.size() && i < 8; ++i)
      cl.push_back({{"csbko", (*csbk[i])["csbko"]},
                    {"feature_set_id", (*csbk[i])["feature_set_id"]},
                    {"csbk_last_block", (*csbk[i])["csbk_last_block"]},
                    {"csbk_data_hex", (*csbk[i])["csbk_data_hex"]}});
    F["csbk"] = cl;
  }
  double t0 = bursts[0]["t_s"].get<double>();
  for (const Json& b : bursts) t0 = std::min(t0, b["t_s"].get<double>());
  Emission e;
  e.family = "dmr";
  e.label = "dmr_" + role;
  e.start_us = t0 * 1e6;
  e.end_us = dur * 1e6;
  e.f_lo_mhz = (c.f_hz - 6.25e3) / 1e6;
  e.f_hi_mhz = (c.f_hz + 6.25e3) / 1e6;
  e.center_mhz = c.f_hz / 1e6;
  if (!base["snr_db"].is_null()) e.snr_db = base["snr_db"].get<double>();
  e.confidence = cc.empty() ? 0.8 : 0.9;
  if (!csbk.empty()) e.decoded = {{"crc_ok", true}};
  e.features = F;
  return e;
}

std::vector<Emission> fsk_carrier(const Channelizer& ch, const Cand& c, double dur) {
  const auto sl = ch.slice(c.f_hz, 6.25e3, kFskFs);
  fvec env(sl.y.size());
  for (size_t i = 0; i < env.size(); ++i) env[i] = np::cabs(sl.y[i]);
  if (static_cast<double>(np::stddev(env) / (np::mean(env) + 1e-30f)) > kMaxEnvCv) return {};
  const fvec d = fsk4::discriminator(sl.y, sl.fs);
  const auto [rate, line] = fsk4::symbol_line(d, sl.fs);
  (void)rate;
  const Json snr = snr_of(sl.y, c.floor, 12.5e3);
  fvec ad(d.size());
  for (size_t i = 0; i < d.size(); ++i) ad[i] = std::fabs(d[i]);
  Json base = Json::object();
  base["symbol_rate_line"] = py_round(line, 1);
  base["snr_db"] = snr;
  base["deviation_hz_p99"] = py_round_int(static_cast<double>(np::percentile(ad, 99.0)));
  const std::vector<Json> bursts = fsk4::dmr(d, sl.fs, 4800.0);
  if (!bursts.empty()) {
    Json b2 = base;
    b2["symbol_rate"] = 4800.0;
    return {dmr_emission(c, dur, bursts, b2)};
  }
  for (double sr : {4800.0, 2400.0}) {
    const auto hits = fsk4::nxdn(d, sl.fs, sr);
    std::vector<int64_t> pos;
    bool exact = false;
    for (const auto& [p, e] : hits) {
      pos.push_back(p);
      exact = exact || e == 0;
    }
    const int64_t n_frame = fsk4::frame_consistent(pos, sl.fs, sr);
    if (n_frame >= 2 || exact) {
      const double half = sr == 2400 ? 3.125e3 : 6.25e3;
      Emission e;
      e.family = "nxdn";
      e.label = "nxdn";
      e.start_us = static_cast<double>(pos[0]) / sl.fs * 1e6;
      e.end_us = dur * 1e6;
      e.f_lo_mhz = (c.f_hz - half) / 1e6;
      e.f_hi_mhz = (c.f_hz + half) / 1e6;
      e.center_mhz = c.f_hz / 1e6;
      if (!snr.is_null()) e.snr_db = snr.get<double>();
      e.confidence = 0.85;
      Json F = base;
      F["symbol_rate"] = sr;
      F["n_frame_sync"] = pos.size();
      F["frame_consistent_hits"] = n_frame;
      F["modulation"] = fmt("4FSK %g sym/s", sr);
      F["channel_bw_khz"] = sr == 2400 ? 6.25 : 12.5;
      e.features = F;
      return {e};
    }
  }
  return {};
}

}  // namespace

const std::vector<std::pair<double, double>>& tetra_ranges() {
  static const std::vector<std::pair<double, double>> r = {
      {380e6, 400e6}, {410e6, 430e6}, {450e6, 470e6}, {870e6, 876e6}, {915e6, 921e6}};
  return r;
}

const std::vector<std::pair<double, double>>& fsk_ranges() {
  static const std::vector<std::pair<double, double>> r = {{136e6, 174e6}, {300e6, 520e6}, {806e6, 870e6}};
  return r;
}

std::vector<Emission> detect(const Channelizer& ch, std::optional<std::pair<double, double>> band,
                             const dvec& exclude_hz) {
  const auto bnd = band ? *band : std::make_pair(ch.fc() - 0.45 * ch.fs(), ch.fc() + 0.45 * ch.fs());
  const double dur = static_cast<double>(ch.L()) / ch.fs();
  std::vector<Emission> out;
  for (const Cand& c : candidates(ch, bnd)) {
    bool skip = false;
    for (double f : exclude_hz) skip = skip || std::fabs(c.f_hz - f) < 3e3;
    if (skip) continue;
    if (in(c.f_hz, tetra_ranges()))
      if (auto e = tetra_carrier(ch, c, dur)) {
        out.push_back(std::move(*e));
        continue;
      }
    if (in(c.f_hz, fsk_ranges()))
      for (auto& e : fsk_carrier(ch, c, dur)) out.push_back(std::move(e));
  }
  return out;
}

}  // namespace si::rx::pmr_digital
