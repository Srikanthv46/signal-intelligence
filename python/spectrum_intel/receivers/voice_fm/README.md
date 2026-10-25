# voice_fm — narrowband FM voice (walkie-talkies, PMR)

**Detects** analog narrowband FM voice carriers in 136–174, 300–520 and 806–870 MHz: walkie-talkies (PMR446,
LPD433, FRS/GMRS) and land-mobile radios.
**Decodes:**
- the channel;
- the sub-audible squelch code: a **DCS** code verified by its Golay(23,12) check, or a **CTCSS** tone;
- deviation, keyed time and speech activity;
- the **audio**.

## How
1. **Candidates.** The capture's spectrum (one FFT, shared with the pipeline) is averaged into 1 kHz bins. A
   candidate is a narrow peak at least 10 dB above the floor and above its surroundings, no more than 40 kHz wide
   at half prominence. Prominence matters: a strong or clipped carrier's skirt stays above the floor far wider
   than its channel.
2. **Per candidate**, the ±10 kHz channel at 48 kSps gives:
   - keyed time: 5 ms frames at least 10 dB over the floor, with gaps under 30 ms bridged;
   - envelope constancy: coefficient of variation ≤ 0.35 (FM);
   - from the FM discriminator: carrier frequency, peak deviation (≤ 6 kHz, i.e. narrowband) and the 300–3400 Hz
     voice band, measured as activity per 20 ms frame and as kurtosis.
3. **Rejected:**
   - an unmodulated carrier (left to the CW rule);
   - two- or four-level FSK data such as POCSAG, DMR or dPMR (a platykurtic voice band with no sub-audible
     code), left to the band rules;
   - digital PMR: a 2400 or 4800 sym/s line ≥ 6× its neighbourhood in the squared derivative of the low-passed
     discriminator output (DMR / NXDN repeaters gave ×11–20, analog FM ×1; their 4FSK had passed as "FM voice
     with CTCSS"); `pmr_digital` identifies them;
   - without a sub-audible code, an envelope that varies more than the noise explains: cv > max(0.12,
     1.6 / sqrt(2 SNR)). Real FM carriers measured 0.12–0.17 against 0.12–0.16 expected; clock / EMI lines with a
     second tone in the channel (two-tone beat) 0.33–0.43 against 0.09–0.20 — they had been "Analog FM voice";
   - the receiver's own clock lines (`ctx.cache["receiver_lines"]`: k × fs, reference-locked lines).
   Speech sensitivity is unchanged (synthetic NBFM speech without a sub-audible code: 12/12 at ≥ 12 dB, before
   and after).
4. **Sub-audible signalling** on the longest keyed segment (0–300 Hz):
   - **DCS:** 134.4 bit/s NRZ, bit clock chosen by eye opening. Every 23-bit window is checked against the
     Golay syndrome (generator 0xC75) and the `100` marker. The result is the code plus its on-air aliases
     (rotations of the same word), and the number of verified words.
   - **CTCSS:** the 60–260 Hz peak, at least 15 dB over the band's median, matched within 0.6 Hz to the EIA tone
     list. A peak that matches nothing is reported as non-standard.
   - DCS is tried first, because its line comb (k × 5.84 Hz) would otherwise read as a tone.
5. **Receiver products.** A stronger carrier at offset f from the LO, when ADC-clipped, produces copies at k·f
   (k = −3, 5, −7, 9, …) with |k| times its deviation. Its IQ image sits at −f. These are dropped and listed on
   the carrier as `suppressed_artifacts`.
6. **Channel plan:** PMR446 (16 × 12.5 kHz), LPD433 (69 × 25 kHz) and FRS/GMRS (22). The match uses the
   pipeline's receiver ppm when it has one. Tolerance is 0.4 channel spacings or 12 ppm, whichever is smaller.
7. **Audio** (`audio.py`, `VoiceDemodulator`):
   - the channel is cut out of one FFT of the block;
   - carrier squelch: channel power at least 8 dB over the noise density 15–100 kHz either side, with 3 dB
     hysteresis and 5 ms ramps;
   - FM discriminator, then 300–3400 Hz band-pass (removes CTCSS/DCS) and 750 µs de-emphasis;
   - fixed gain from the rated deviation, so silence stays silent;
   - decimate to 16 kHz.

## Labels
- **"Analog FM voice (PMR446 ch 8, D053N)":** a sub-audible code or speech was seen. Status is `success` when the
  code is verified (DCS Golay or a standard CTCSS tone).
- **"Analog FM voice-like (…)":** a modulated narrowband FM carrier with neither.

DCS polarity: a 1 bit is taken as a positive frequency shift (N). The equivalent codes with swapped polarity are
always listed in `dcs_aliases`.

## Not decoded
- **Scrambled voice** (inversion or rolling-code scramblers) and **digital voice** (DMR, dPMR, NXDN, P25, TETRA)
  are not turned into audio.
- Digital voice is not FM-voice-shaped and falls to the band rules (`PMR land-mobile-like`).
- Inverted audio is still detected as an FM voice carrier, with its channel and squelch code, but plays back
  unintelligible.

**Short snapshots:** 32.5 ms sweep snapshots are too short for DCS (at least 0.35 s) or CTCSS (at least 0.3 s).
The carrier, channel and deviation are still reported.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.voice_fm import detect, VoiceDemodulator, demodulate
for e in detect(x, fs, fc):          # Emission(family="voice_fm", label="nbfm_voice_pmr446_ch8", ...)
    print(e.center_mhz, e.features.get("channel"), e.features.get("dcs_code") or e.features.get("ctcss_hz"))
demod = VoiceDemodulator.from_fields(report_signal["decoded"]["fields"])   # streaming, one block at a time
chunk = demod.process(iq_block, fs, fc)   # AudioChunk(samples float32 [-1, 1], fs 16 kHz, meta) or None
```
- signal-vault's Signal Intelligence tab uses `VoiceDemodulator` for its **🔊 Listen** button, live and on files.
- `spectrum-intel run` saves `<name>_voice<MHz>_audio.wav` for every voice signal and records the file name in
  `decoded.fields.audio_wav`.
- `synth.py` generates synthetic NBFM transmitters (speech-like audio, tones, CTCSS, DCS, keying, noise) for
  the tests.

## Pipeline plugin (`plugin.py`)
- Runs when the capture overlaps one of the voice ranges.
- Produces one signal per carrier, with `decoded.fields.audio = "nbfm"`: the viewer's cue to offer Listen.
- Owns its band.

## Validation
`tests/receivers/test_det_voice_fm.py`:
- **DCS convention:** reproduces the published equivalences 023N = 340N = 766N = 047I, which the reciprocal
  polynomial does not, and all 104 standard codes are valid codewords.
- **Synthetic CTCSS:** tones 67.0, 88.5, 151.4 and 254.1 Hz on PMR446 channels 1, 8, 12 and 16 at 20 dB SNR, with
  exact tone and channel.
- **Synthetic DCS:** D023N, D053N, D411I and D754N, each with at least 10 verified words.
- **Synthetic keying and audio:** the keyed burst is timed within 30 ms and the deviation measured. Audio of a
  1 kHz tone peaks at 1 kHz, the CTCSS is suppressed by at least 30 dB, and the squelch closes when the carrier
  drops.
- **Real PMR446 walkie-talkie** (`testdata/voice`; Pluto, 2026-10-06, keyed next to the antenna):
  - the capture is hard-clipped;
  - one carrier, PMR446 channel 8 (measured −2.76 kHz, −6.2 ppm, which is the Pluto's +7.4 ppm);
  - DCS D053N (alias D452I) with 14–16 Golay-verified words in every one of 25 blocks;
  - four clipping products suppressed.
  - The code is not yet confirmed from the radio's display.
- **No false alarms:**
  - noise;
  - an unmodulated carrier;
  - 2400 Bd 2-FSK data;
  - the real 434 MHz LoRa capture;
  - 25 blocks of 433.92 MHz with no walkie-talkie.

**Known gaps:**
- no real capture with speech yet;
- DCS polarity not confirmed against a radio's display;
- DCS turn-off code and reverse-burst tail not reported.

## AM voice (research, signal gaps S13)

**Validation: research** (synthetic only; `research/signal_gaps/05_aviation_radar_jammers.md` section 6).
- **Where:** only the aeronautical bands, VHF civil 118–137 MHz and UHF military 225–400 MHz (`air_channel`: service,
  nearest 25 kHz channel, offset; 121.5 / 243.0 MHz emergency). Elsewhere OOK remotes and pulsed data also move the
  envelope.
- **Test:** on the ±10 kHz channel at 48 kSps, the envelope's 300–3400 Hz rms relative to the carrier is ≥ 4 % and
  3 × what noise gives, in ≥ 10 % of 20 ms frames, while the discriminator's speech-band rms stays at the noise level
  (phase quiet: not FM); SNR ≥ 10 dB; no FSK symbol-rate line.
- **Label** "Analog AM voice (<service>, <channel> MHz)" or the emergency name, 0.7 (0.8 with speech in ≥ 20 % of the
  frames), status `partial`; `fields.audio = "am"`: the viewer's Listen and `demodulate` give envelope audio
  (300–3400 Hz, carrier squelch).

## APRS / AX.25 packet radio (`afsk.py`; research, signal gaps S19)

**Validation: research** (`research/signal_gaps/06_vhf_below_300.md` section 7; published APRS examples and
synthetic AFSK).
- On every FM voice channel: Bell 202 AFSK (mark 1200 / space 2200 Hz) tone energy over one bit at 8 bit phases,
  NRZI, HDLC frames whose CRC-16/X.25 is valid (`core/hdlc`), AX.25 UI address / control / PID checks.
- **APRS** fields: uncompressed, compressed and Mic-E positions (lat / lon, symbol, course / speed, comment),
  timestamps, status, messages, objects.
- **Label** "APRS (<callsigns>)" or "AX.25 packet radio (<callsigns>)", 0.97, status `success` (every packet
  CRC-verified); no audio is offered for packet channels.
