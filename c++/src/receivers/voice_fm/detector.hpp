// Narrowband FM voice (walkie-talkies, PMR, land-mobile) detector
// (python: receivers/voice_fm/detector.py).
//
//   detect(ch, ppm) -> Emission family "voice_fm", label nbfm_voice[_<service>_ch<n>]
//
// Candidates are narrow peaks of the capture spectrum (1 kHz bins) in the voice ranges;
// each is measured on its +-10 kHz channel at 48 kSps (keying, envelope, discriminator,
// voice band, CTCSS/DCS); clipping products / IQ images of stronger carriers are dropped.
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::voice_fm {

const std::vector<std::pair<double, double>>& voice_ranges();   // VOICE_RANGES
// python air_channel(f_hz, ppm): aeronautical service, nearest 25 kHz channel, emergency frequency ({} outside)
Json air_channel(double f_hz, double ppm = 0.0);

struct Candidate {
  double f_hz, peak_db, width_hz, noise_density;
};
std::vector<Candidate> candidates(const Channelizer& ch, size_t max_n = 12);

struct Measurement {
  double f_hz = 0;
  std::vector<std::pair<int64_t, int64_t>> segs;   // keyed segments, in 5 ms frames
  double keyed_s = 0, duration_s = 0, envelope_cv = 0, deviation_hz = 0;
  double voice_rms_hz = 0, voice_kurtosis = 0, voice_activity = 0;
  double fsk_line_ratio = 0;                        // python fsk_line["ratio"]
  double am_voice_rms = 0, am_activity = 0;         // AM: envelope speech band / carrier (research, S13)
  std::vector<Json> packets;                        // AX.25 / APRS frames (research, S19)
  std::optional<Json> sub;
  std::optional<double> snr_db;
  double peak_db = 0;
  Json artifacts = Json::array();
};
std::optional<Measurement> measure(const Channelizer& ch, const Candidate& c);

// python detect(x, fs, fc, {"ch": ch, "ppm": ppm, "exclude_hz": exclude_hz})
std::vector<Emission> detect(const Channelizer& ch, double ppm = 0.0, const dvec& exclude_hz = {});
std::vector<Emission> detect(const cfvec& x, double fs, double fc, double ppm = 0.0);

}  // namespace si::rx::voice_fm
