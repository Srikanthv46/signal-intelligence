# voice_fm (C++) — narrowband FM voice (walkie-talkies, PMR)

Port of `python/spectrum_intel/receivers/voice_fm` (detector.py, subaudio.py, channels.py, audio.py, plugin.py).

## What it does
- **Detects** analog NBFM voice carriers in 136–174, 300–520 and 806–870 MHz: narrow peaks of the 1 kHz-bin
  spectrum, each measured on its ±10 kHz channel at 48 kSps (keyed time, envelope constancy, discriminator
  carrier and deviation, voice-band activity / kurtosis); FSK data and unmodulated carriers are rejected;
  clipping products and IQ images of stronger carriers are dropped and listed (`suppressed_artifacts`).
- **Decodes** DCS (Golay(23,12)-verified, with aliases) or CTCSS, the channel (PMR446 / LPD433 / FRS-GMRS, with
  the receiver ppm from `ctx.cache` "ppm" / "sweep_ppm"), and the **audio** (`Signal::audio`, float32 at 16 kHz).

## C++ API (`namespace si::rx::voice_fm`, depends on core/ only)
```cpp
#include "receivers/voice_fm/detector.hpp"
#include "receivers/voice_fm/audio.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc, ppm);   // or detect(channelizer, ppm)

// streaming audio, filter states carried across blocks (python VoiceDemodulator)
VoiceDemodulator demod = VoiceDemodulator::from_fields(report_signal["decoded"]["fields"]);
if (auto chunk = demod.process(block, fs, fc)) play(chunk->samples, chunk->fs);   // chunk->meta as Python
auto whole = demodulate(iq, fs, fc, fields);                 // one shot
write_wav("voice.wav", whole->samples, whole->fs);           // 16-bit PCM, byte-identical to Python's wave output
```
`subaudio.hpp`: `golay_syndrome`, `dcs_encode`, `dcs_aliases`, `to_subaudio`, `dcs`, `ctcss`, `subaudible`;
`channels.hpp`: plans, CTCSS tones, DCS codes, `match_channel`.

## Implementation notes
- Same dtypes as numpy (complex64 channel from `core/dsp` Channelizer, float32 discriminator, float64 filters);
  Python-float thresholds compared in float32 (NEP 50).
- numpy temporary elision in `y[1:] * np.conj(y[:-1])` (≥ 256 KiB → numpy computes `conj(prev) * cur`) is
  mirrored; without it ~75 % of discriminator samples differed in the last bit.
- The plugin passes the pipeline's Channelizer to the audio demodulator (Python builds an identical one).

## Parity / tests
- `tools/parity.py --receivers analog_vtx,voice_fm`: 16/16 inputs identical; `real_pmr446_ch8_dcs`
  (PMR446 ch 8, D053N, 5 suppressed clipping products) matches in every field, incl. audio_s / squelch.
- `tests/unit/test_rx_voice_fm.cpp` (vectors: `tests/vectors/voice_fm/gen_voice_fm.py`): DCS tables, aliases
  and syndromes exact; candidates, per-candidate measurements (to 1e-9), Emission features, one-shot and
  streaming audio (max |diff| < 1e-5; np.convolve/lfilter use BLAS dot order in Python) on the real capture and a
  synthetic speech + CTCSS 88.5 Hz transmitter. WAV bytes verified identical to Python's `wave` writer.
- Timing (receiver_s, PMR446 capture): C++ 0.04 s vs Python 0.08 s.

## 2026-10-10 (ported with Python)
- Rejects a carrier whose discriminator shows a 2400 / 4800 sym/s FSK line (>= 6x its neighbourhood): digital PMR.
- Without a sub-audible code, rejects an envelope cv > max(0.12, 1.6 / sqrt(2 SNR)) (two-tone beats, clock lines).
- Skips the pipeline's receiver lines (`ctx.cache["receiver_lines"]`).

## AM voice and APRS / AX.25 (research, signal gaps S13 / S19)
- AM voice in the aeronautical bands (`air_channel`, envelope speech band vs carrier, phase quiet) → "Analog AM voice
  (...)"; `audio.cpp` AM mode (envelope, 300–3400 Hz).
- `afsk.{hpp,cpp}` (python `afsk.py`): Bell 202 AFSK on any FM channel, HDLC / CRC-16/X.25 / AX.25 / APRS through
  `core/hdlc` → "APRS (...)" / "AX.25 packet radio (...)". Vectors: `tests/vectors/signal_gaps/am*`, `aprs*`.
