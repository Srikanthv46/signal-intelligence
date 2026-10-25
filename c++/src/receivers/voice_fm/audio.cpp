// python: receivers/voice_fm/audio.py (and the WAV format of cli._write_audio)
#include "receivers/voice_fm/audio.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>

#include "core/npcompat.hpp"
#include "core/signal.hpp"

namespace si::rx::voice_fm {

VoiceDemodulator::VoiceDemodulator(double center_hz, double deviation_hz, double channel_bw_hz, int audio_fs,
                                   bool squelch, bool am)
    : center_hz_(center_hz), channel_bw_hz_(channel_bw_hz), audio_fs_(audio_fs), squelch_(squelch), am_(am) {
  ref_dev_ = std::max(deviation_hz, channel_bw_hz <= 12.5e3 ? 2500.0 : 5000.0);
  mid_ = 3.0 * audio_fs_;
  bp_ = sig::butter_sos(4, dvec{300, 3400}, sig::BType::Bandpass, mid_);
  zbp_ = dmat(bp_.rows, 2, 0.0);
  const double a = std::exp(-1 / (mid_ * kDeemphS));
  de_b_ = {1 - a};
  de_a_ = {1, -a};
  zde_ = dvec(1, 0.0);
  dec_ = sig::firwin(63, 0.9 * audio_fs_ / 2, mid_);
  zdec_ = dvec(dec_.size() - 1, 0.0);
  const double w = 2 * kPi * 1000 * kDeemphS;
  const double g1k = 1 / std::sqrt(1 + w * w);
  gain_ = 0.5 / (ref_dev_ * g1k);
}

VoiceDemodulator VoiceDemodulator::from_fields(const Json& f, int audio_fs) {
  auto num = [&](const char* k, double dflt) {
    if (!f.contains(k) || !f[k].is_number()) return dflt;
    const double v = f[k].get<double>();
    return v != 0 ? v : dflt;   // python: f.get(k) or default
  };
  const bool am = f.contains("audio") && f["audio"] == "am";
  return VoiceDemodulator(f.at("center_mhz").get<double>() * 1e6, num("deviation_hz", 2500),
                          am ? 8e3 : num("channel_bw_khz", 12.5) * 1e3, audio_fs, true, am);
}

bool VoiceDemodulator::covers(double fs, double fc) const {
  return std::fabs(center_hz_ - fc) + channel_bw_hz_ / 2 < 0.48 * fs;
}

dvec VoiceDemodulator::gate(const dvec& snr_db, size_t n_frame, size_t n) {
  dvec g(snr_db.size());
  bool st = open_;
  for (size_t i = 0; i < snr_db.size(); ++i) {   // 3 dB hysteresis
    st = snr_db[i] >= (st ? kSquelchDb - 3 : kSquelchDb);
    g[i] = st;
  }
  open_ = st;
  dvec r;
  r.reserve(n);
  for (double v : g) r.insert(r.end(), n_frame, v);
  const double last = g.empty() ? 0.0 : g.back();
  while (r.size() < n) r.push_back(last);
  const size_t k = std::max<size_t>(1, static_cast<size_t>(0.005 * mid_));
  return np::convolve(r, dvec(k, 1.0 / static_cast<double>(k)), np::ConvMode::Same);   // 5 ms ramps
}

std::optional<AudioChunk> VoiceDemodulator::process(const cfvec& iq, double fs, double fc, const Channelizer* chp) {
  if (!covers(fs, fc) || static_cast<double>(iq.size()) < fs * 0.005) return std::nullopt;
  std::unique_ptr<Channelizer> own;
  if (!chp) {
    own = std::make_unique<Channelizer>(iq, fs, fc);
    chp = own.get();
  }
  const Channelizer::Slice sl = chp->slice(center_hz_, channel_bw_hz_ / 2, mid_);
  const cfvec& y = sl.y;
  const double fs_y = sl.fs;
  // noise density beside the channel (median: robust to neighbouring carriers)
  fvec side;
  const double c_off = center_hz_ - fc;
  for (size_t i = 0; i < chp->L(); ++i) {
    const double fi = chp->f()[i];
    const double off = std::fabs(fi - c_off);
    if (off > 15e3 && off < 100e3 && std::fabs(fi) < 0.48 * fs) side.push_back(np::cabs2(chp->X()[i]));
  }
  const bool have_side = !side.empty();
  const double nd = have_side ? static_cast<double>(np::median(side)) / std::log(2.0) / (fs * static_cast<double>(iq.size())) : 0.0;
  own.reset();
  const size_t n_frame = std::max<size_t>(1, static_cast<size_t>(0.01 * fs_y));
  const size_t nf = y.size() / n_frame;
  dvec snr(nf);
  {
    fvec a(n_frame);
    for (size_t i = 0; i < nf; ++i) {
      for (size_t j = 0; j < n_frame; ++j) a[j] = np::cabs2(y[i * n_frame + j]);
      const float p = np::mean(a);
      if (have_side) {   // nd is a numpy float64 (strong): float64 arithmetic
        snr[i] = 10 * std::log10(static_cast<double>(p) / (nd * channel_bw_hz_ + 1e-30) + 1e-30);
      } else {           // nd = 0.0 (python float, weak): float32 arithmetic
        const float r = p / static_cast<float>(nd * channel_bw_hz_ + 1e-30) + 1e-30f;
        snr[i] = static_cast<double>(10.0f * std::log10(r));
      }
    }
  }
  const dvec g = squelch_ ? gate(snr, n_frame, y.size()) : dvec(y.size(), 1.0);
  dvec v;
  if (am_) {   // envelope; the band-pass removes the carrier (DC); 100 % modulation -> 0.5 FS
    dvec e(y.size());
    for (size_t i = 0; i < y.size(); ++i) e[i] = static_cast<double>(np::cabs(y[i]));
    v = sig::sosfilt(bp_, e, &zbp_);
    const double k = 0.5 / (np::mean(e) + 1e-30);
    for (size_t i = 0; i < v.size(); ++i) v[i] = v[i] * k * g[i];
  } else {
  const double sc = static_cast<double>(static_cast<float>(fs_y / (2 * kPi)));
  dvec d(y.size(), 0.0);
  // numpy temporary elision: in `y[1:] * np.conj(y[:-1])` the conj() temporary, when it holds
  // >= 256 KiB, receives the product in place, i.e. numpy computes conj(prev) * cur; the
  // fma-based complex multiply is not symmetric in its last bit, so the order is mirrored.
  const bool elide = y.size() > 1 && (y.size() - 1) * sizeof(cf) >= 256 * 1024;
  for (size_t i = 1; i < y.size(); ++i) {
    const cf z = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
    d[i] = static_cast<double>(std::atan2(z.imag(), z.real()) * static_cast<float>(sc));
  }
  if (y.size() > 1) d[0] = d[1];
  v = sig::sosfilt(bp_, d, &zbp_);
  v = sig::lfilter(de_b_, de_a_, v, &zde_);
  for (size_t i = 0; i < v.size(); ++i) v[i] = v[i] * gain_ * g[i];
  }
  v = sig::lfilter(dec_, dvec{1.0}, v, &zdec_);
  AudioChunk out;
  out.samples.reserve(v.size() / 3 + 1);
  for (size_t i = 0; i < v.size(); i += 3) out.samples.push_back(static_cast<float>(std::min(std::max(v[i], -1.0), 1.0)));
  out.fs = fs_y / 3;
  size_t open = 0;
  for (double x : g) open += x > 0.5;
  out.meta["squelch_open"] = g.empty() ? std::nan("") : static_cast<double>(open) / static_cast<double>(g.size());
  out.meta["snr_db"] = snr.empty() ? Json(nullptr) : Json(py_round(np::median(snr), 1));
  out.meta["center_mhz"] = center_hz_ / 1e6;
  return out;
}

std::optional<AudioChunk> demodulate(const cfvec& iq, double fs, double fc, const Json& fields, int audio_fs,
                                     const Channelizer* ch) {
  return VoiceDemodulator::from_fields(fields, audio_fs).process(iq, fs, fc, ch);
}

std::vector<uint8_t> wav_bytes(const fvec& x, double fs) {
  std::vector<uint8_t> o;
  auto le = [&](uint32_t v, int n) {
    for (int i = 0; i < n; ++i) o.push_back(static_cast<uint8_t>(v >> (8 * i)));
  };
  const uint32_t data = static_cast<uint32_t>(x.size() * 2);
  const uint32_t rate = static_cast<uint32_t>(py_round_int(fs));
  o.insert(o.end(), {'R', 'I', 'F', 'F'});
  le(36 + data, 4);
  o.insert(o.end(), {'W', 'A', 'V', 'E', 'f', 'm', 't', ' '});
  le(16, 4);
  le(1, 2);          // PCM
  le(1, 2);          // mono
  le(rate, 4);
  le(rate * 2, 4);   // byte rate
  le(2, 2);          // block align
  le(16, 2);         // bits
  o.insert(o.end(), {'d', 'a', 't', 'a'});
  le(data, 4);
  for (float v : x) {
    const float c = std::min(std::max(v, -1.0f), 1.0f) * 32767.0f;
    le(static_cast<uint16_t>(static_cast<int16_t>(c)), 2);
  }
  return o;
}

bool write_wav(const std::string& path, const fvec& x, double fs) {
  const auto b = wav_bytes(x, fs);
  std::ofstream f(path, std::ios::binary);
  f.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
  return f.good();
}

}  // namespace si::rx::voice_fm
