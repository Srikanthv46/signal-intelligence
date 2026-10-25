// Narrowband FM voice -> audio, one block at a time (python: receivers/voice_fm/audio.py,
// VoiceDemodulator / AudioChunk / demodulate).
//
//   auto demod = VoiceDemodulator::from_fields(signal["decoded"]["fields"]);
//   if (auto chunk = demod.process(iq, fs, fc)) play(chunk->samples, chunk->fs);
//
// Per block: channel slice at 3 x audio_fs; carrier squelch per 10 ms frame (>= 8 dB over
// the noise density 15-100 kHz beside the channel, 3 dB hysteresis, 5 ms ramps); FM
// discriminator; 300-3400 Hz band-pass; 750 us de-emphasis; fixed gain from the rated
// deviation; decimate by 3. Filter states carry across blocks.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/dsp.hpp"
#include "core/types.hpp"

namespace si::rx::voice_fm {

constexpr double kDeemphS = 750e-6;
constexpr double kSquelchDb = 8.0;

struct AudioChunk {
  fvec samples;     // float32 mono, [-1, 1]
  double fs = 0;
  Json meta = Json::object();   // squelch_open, snr_db, center_mhz
};

class VoiceDemodulator {
 public:
  VoiceDemodulator(double center_hz, double deviation_hz = 2500.0, double channel_bw_hz = 12.5e3, int audio_fs = 16000,
                   bool squelch = true, bool am = false);
  static VoiceDemodulator from_fields(const Json& fields, int audio_fs = 16000);

  bool covers(double fs, double fc) const;
  // one block -> audio (nullopt if the block does not cover the channel or is < 5 ms).
  // ch: an existing Channelizer of this block (optional; the pipeline passes its own).
  std::optional<AudioChunk> process(const cfvec& iq, double fs, double fc, const Channelizer* ch = nullptr);

  double center_hz() const { return center_hz_; }
  double gain() const { return gain_; }
  int audio_fs() const { return audio_fs_; }

 private:
  double center_hz_, channel_bw_hz_, ref_dev_;
  int audio_fs_;
  bool squelch_;
  bool am_;                       // envelope (AM voice, research) instead of the FM discriminator
  double mid_;
  dmat bp_, zbp_;
  dvec de_b_, de_a_, zde_, dec_, zdec_;
  double gain_;
  bool open_ = false;
  dvec gate(const dvec& snr_db, size_t n_frame, size_t n);
};

// python demodulate(iq, fs, fc, fields): whole capture -> audio (one shot)
std::optional<AudioChunk> demodulate(const cfvec& iq, double fs, double fc, const Json& fields, int audio_fs = 16000,
                                     const Channelizer* ch = nullptr);

// 16-bit mono PCM WAV, byte-identical to python cli._write_audio (wave module, int16 of
// clip(x, -1, 1) * 32767 truncated toward zero, frame rate round(fs))
std::vector<uint8_t> wav_bytes(const fvec& samples, double fs);
bool write_wav(const std::string& path, const fvec& samples, double fs);

}  // namespace si::rx::voice_fm
