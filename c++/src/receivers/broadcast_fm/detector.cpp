// python receivers/broadcast_fm/detector.py (see there for the method and the sources).
#include "receivers/broadcast_fm/detector.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/signal.hpp"
#include "receivers/broadcast_fm/rds.hpp"

namespace si::rx::broadcast_fm {

namespace {

constexpr double kBandLo = 87.5e6, kBandHi = 108.0e6;
constexpr double kMpxFs = 240e3;
constexpr double kChanDb = 10.0, kPilotDb = 10.0;
constexpr double kRdsRate = 1187.5;
constexpr double kMinRdsS = 0.045;

// python Channelizer._bins(f_off, half_bw): ascending unique FFT indices within half_bw of f_off
std::vector<size_t> bins(const Channelizer& ch, double f_off, double half) {
  const size_t L = ch.L();
  const double df = ch.fs() / static_cast<double>(L);
  int64_t lo = static_cast<int64_t>(std::floor((f_off - half) / df)) - 2;
  int64_t hi = static_cast<int64_t>(std::ceil((f_off + half) / df)) + 2;
  lo = std::max<int64_t>(lo, -static_cast<int64_t>(L / 2) - 1);
  hi = std::min<int64_t>(hi, static_cast<int64_t>((L - 1) / 2) + 1);
  std::vector<size_t> idx;
  for (int64_t k = lo; k <= hi; ++k) idx.push_back(static_cast<size_t>(((k % static_cast<int64_t>(L)) + L) % L));
  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  std::vector<size_t> out;
  for (size_t i : idx)
    if (std::fabs(ch.f()[i] - f_off) <= half) out.push_back(i);
  return out;
}

double power(const Channelizer& ch, double f0, double half) {
  const auto k = bins(ch, f0 - ch.fc(), half);
  if (k.empty()) return 0.0;
  dvec p(k.size());
  for (size_t i = 0; i < k.size(); ++i) {
    const double a = static_cast<double>(np::cabs(ch.X()[k[i]]));
    p[i] = a * a;
  }
  return np::mean(p);
}

double line_db(const dvec& spec, const dvec& fr, double f0) {
  double on = -1;
  dvec bg;
  for (size_t i = 0; i < fr.size(); ++i) {
    const double d = std::fabs(fr[i] - f0);
    if (d <= 15) on = std::max(on, spec[i]);
    if (d > 500 && d < 2000) bg.push_back(spec[i]);
  }
  if (on < 0 || bg.empty()) return 0.0;
  return 10 * std::log10(on / (np::median(bg) + 1e-30) + 1e-30);
}

std::optional<Json> measure(const Channelizer& ch, double f_mhz) {
  const auto sl = ch.slice(f_mhz * 1e6, 120e3, kMpxFs);
  const cfvec& y = sl.y;
  const double fs = sl.fs;
  if (y.size() < 1024) return std::nullopt;
  const float sc = static_cast<float>(fs / (2 * kPi));
  fvec d(y.size() - 1);
  // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
  const bool elide = (y.size() - 1) * sizeof(cf) >= 256 * 1024;
  for (size_t i = 1; i < y.size(); ++i) {
    const cf z = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
    d[i - 1] = std::atan2(z.imag(), z.real()) * sc;
  }
  const float offf = np::mean(d);
  const double off = static_cast<double>(offf);
  fvec mpx(d.size());
  dvec ab(d.size()), win(d.size());
  const dvec hann = np::hanning(d.size());
  for (size_t i = 0; i < d.size(); ++i) {
    mpx[i] = d[i] - offf;
    ab[i] = static_cast<double>(std::fabs(mpx[i]));
    win[i] = static_cast<double>(mpx[i]) * hann[i];
  }
  const double dev = np::percentile(ab, 99.5);
  const cdvec S = fft::rfft(win);
  dvec spec(S.size());
  for (size_t i = 0; i < S.size(); ++i) {
    const double a = np::cabs(S[i]);
    spec[i] = a * a;
  }
  const dvec fr = fft::rfftfreq(mpx.size(), 1 / fs);
  const double pilot_db = line_db(spec, fr, 19e3);
  const bool stereo = pilot_db >= kPilotDb;
  double f_pilot = 19e3;
  if (stereo) {
    double best = -1;
    for (size_t i = 0; i < fr.size(); ++i)
      if (std::fabs(fr[i] - 19e3) <= 15 && spec[i] > best) best = spec[i], f_pilot = fr[i];
  }
  Json out = Json::object();
  out["carrier_mhz"] = py_round((f_mhz * 1e6 + off) / 1e6, 6);
  out["deviation_khz"] = py_round(dev / 1e3, 1);
  out["stereo"] = stereo;
  out["pilot_db"] = py_round(pilot_db, 1);
  out["carrier_offset_hz"] = py_round_int(off);
  if (static_cast<double>(mpx.size()) / fs >= kMinRdsS) {
    const auto bits = rds_bits(mpx, fs, 3 * f_pilot);
    const auto blks = bits.empty() ? std::vector<rds::Block>{} : rds::blocks(bits);
    if (!blks.empty()) {
      Json info = rds::decode(rds::groups(blks));
      info["blocks_ok"] = static_cast<int64_t>(blks.size());
      info["bits"] = static_cast<int64_t>(bits.size());
      out["rds"] = info;
    }
  }
  return out;
}

}  // namespace

std::vector<int> rds_bits(const fvec& mpx, double fs, double f57) {
  const size_t n = mpx.size();
  dvec re(n), im(n);
  const double w = -2 * kPi * f57;
  for (size_t i = 0; i < n; ++i) {
    const double th = w * (static_cast<double>(i) / fs);
    const double m = mpx[i];
    re[i] = m * std::cos(th);
    im[i] = m * std::sin(th);
  }
  const dmat sos = sig::butter_sos(4, 2400.0, sig::BType::Lowpass, fs);
  const dvec zr = sig::sosfiltfilt(sos, re), zi = sig::sosfiltfilt(sos, im);
  cdvec z2(n);
  for (size_t i = 0; i < n; ++i) {
    const cd z(zr[i], zi[i]);
    z2[i] = np::cmul(z, z);
  }
  const cd s = np::sum(z2);
  const double th = 0.5 * std::atan2(s.imag(), s.real());
  const cd e(std::cos(th), -std::sin(th));
  dvec r(n);
  for (size_t i = 0; i < n; ++i) r[i] = np::cmul(cd(zr[i], zi[i]), e).real();
  const double spb = fs / kRdsRate;
  const int64_t nbits = static_cast<int64_t>(static_cast<double>(n) / spb) - 1;
  if (nbits < 52) return {};
  const double half = spb / 2;
  auto seg = [&](int64_t a, int64_t b) {
    a = std::max<int64_t>(0, std::min<int64_t>(a, static_cast<int64_t>(n)));
    b = std::max<int64_t>(a, std::min<int64_t>(b, static_cast<int64_t>(n)));
    return np::sum(r.data() + a, static_cast<size_t>(b - a));
  };
  double best = 0.0, best_e = -1.0;
  const double step = spb / 32;
  for (int q = 0; q < 32; ++q) {
    const double off = static_cast<double>(q) * step;
    double e = 0.0;
    for (int64_t k = 0; k < nbits; ++k) {
      const int64_t a = static_cast<int64_t>(off + static_cast<double>(k) * spb);
      const int64_t m = static_cast<int64_t>(off + static_cast<double>(k) * spb + half);
      const int64_t b = static_cast<int64_t>(off + static_cast<double>(k + 1) * spb);
      e += std::fabs(seg(a, m) - seg(m, b));
    }
    if (e > best_e) best = off, best_e = e;
  }
  std::vector<int> enc;
  for (int64_t k = 0; k < nbits; ++k) {
    const int64_t a = static_cast<int64_t>(best + static_cast<double>(k) * spb);
    const int64_t m = static_cast<int64_t>(best + static_cast<double>(k) * spb + half);
    const int64_t b = static_cast<int64_t>(best + static_cast<double>(k + 1) * spb);
    enc.push_back(seg(a, m) - seg(m, b) > 0 ? 1 : 0);
  }
  std::vector<int> out;
  for (size_t k = 1; k < enc.size(); ++k) out.push_back(enc[k] ^ enc[k - 1]);
  return out;
}

std::vector<std::pair<double, double>> candidates(const Channelizer& ch) {
  const double lo = ch.fc() - 0.45 * ch.fs(), hi = ch.fc() + 0.45 * ch.fs();
  std::vector<double> chans;
  for (int k = 0; k < 204; ++k) {
    const double f = py_round(87.6 + 0.1 * k, 1);
    if (lo + 120e3 < f * 1e6 && f * 1e6 < hi - 120e3) chans.push_back(f);
  }
  if (chans.empty()) return {};
  dvec sel;
  for (size_t i = 0; i < ch.L(); ++i) {
    const double fa = ch.f()[i] + ch.fc();
    if (std::fabs(ch.f()[i]) < 0.45 * ch.fs() && fa >= kBandLo && fa <= kBandHi) {
      const double a = static_cast<double>(np::cabs(ch.X()[i]));
      sel.push_back(a * a);
    }
  }
  if (sel.size() < 64) return {};
  const double floor = np::median(sel) / std::log(2.0);
  std::vector<std::pair<double, double>> out, keep;
  for (double f : chans) {
    const double db = 10 * std::log10(power(ch, f * 1e6, 100e3) / floor + 1e-30);
    if (db >= kChanDb) out.emplace_back(f, db);
  }
  for (const auto& [f, db] : out) {
    bool top = true;
    for (const auto& [g, dg] : out)
      if (std::fabs(g - f) < 0.15 && g != f && db < dg) top = false;
    if (top) keep.emplace_back(f, db);
  }
  return keep;
}

std::vector<Emission> detect(const Channelizer& ch) {
  if (!(ch.fc() - ch.fs() / 2 < kBandHi && ch.fc() + ch.fs() / 2 > kBandLo)) return {};
  const double dur_us = static_cast<double>(ch.L()) / ch.fs() * 1e6;
  std::vector<Emission> out;
  for (const auto& [f, db] : candidates(ch)) {
    const auto m = measure(ch, f);
    if (!m) continue;
    const double dk = (*m)["deviation_khz"].get<double>();
    if (!(15.0 <= dk && dk <= 110.0)) continue;
    const bool rds_ok = m->contains("rds") && (*m)["rds"].contains("pi");
    const double conf = 0.7 + ((*m)["stereo"].get<bool>() ? 0.1 : 0.0) + (rds_ok ? 0.15 : 0.0);
    Json F = Json::object();
    F["channel_mhz"] = f;
    F["level_db"] = py_round(db, 1);
    for (auto it = m->begin(); it != m->end(); ++it) F[it.key()] = it.value();
    Emission e;
    e.family = "broadcast_fm";
    e.label = "fm_broadcast";
    e.start_us = 0.0;
    e.end_us = dur_us;
    e.f_lo_mhz = f - 0.1;
    e.f_hi_mhz = f + 0.1;
    e.center_mhz = (*m)["carrier_mhz"].get<double>();
    e.snr_db = py_round(db, 1);
    e.confidence = py_round(std::min(conf, 0.97), 2);
    e.decoded = m->contains("rds") ? (*m)["rds"] : Json(nullptr);
    e.features = std::move(F);
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace si::rx::broadcast_fm
