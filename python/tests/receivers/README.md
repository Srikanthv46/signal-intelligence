# receiver validation scripts

Moved from `wifi-decoding/test_*.py`. They check the receivers in `spectrum_intel/receivers/`
against the real captures in `testdata/` (copied from wifi-decoding) and synthetic signals.
Run them from the package root:

    cd tests/receivers && ../../.venv/bin/python test_det_digital_vtx.py

| script | checks | runtime |
|---|---|---|
| `test_ofdm15k.py` | 15 kHz OFDM hop finder, Wi-Fi-numerology rejection, noise -> 0 | seconds |
| `test_det_digital_vtx.py` | Herelink carrier / PCI 411, no false alarms on DJI + Wi-Fi + noise | ~1 min |
| `test_det_dji.py` | OcuSync / DroneID on the DJI captures (+ DroneID test vectors if the decoder repo is present) | minutes |
| `test_drone_assoc.py` | burst -> emitter association (also collected by pytest) | ~1 min |
| `test_det_analog_vtx.py` | synthetic PAL/NTSC Pd vs SNR, no false alarms on real captures | minutes |
| `test_analog_color.py` | colour video: synthetic NTSC/PAL colour bars (clean, 12 dB, ±ppm, audio subcarrier), colour killer, the real NTSC VTX in colour with luma unchanged (also collected by pytest) | ~10 s |
| `test_det_rc_link.py` | LoRa / GFSK Pd vs SNR sweeps (`--quick`, `--procs 2`) | long |
| `test_det_lora.py` | LoRa payloads byte-exact vs gr-lora_sdr's transmitter (SF7–11, all CRs, CRC on/off, LDRO) and a real 434 MHz packet decoded by gr-lora_sdr; CRC/whitening check values; LoRaWAN/Meshtastic header parsing (also collected by pytest). `gen_lora_vectors.py` regenerates the vectors (radioconda Python) | ~10 s |
| `test_det_voice_fm.py` | FM voice: DCS convention vs published code equivalences, CTCSS tones / DCS codes (both polarities) / PMR446 channels on synthetic NBFM, keyed-burst timing and deviation, demodulated audio (tone, CTCSS removal, squelch), the real clipped PMR446 radio (one carrier, ch 8, D053N), no detection on noise / CW / 2-FSK data / real LoRa (also collected by pytest) | ~1 s |
| `test_det_wifi_dsss.py` | 802.11b all rates/preambles byte-exact, sensitivity, CFO/clock/multipath, FCS correction, excision, MAC dissector, Phantom 4 regression (also collected by pytest; run directly for the Pd-vs-SNR table) | ~20 s (table: ~3 min) |
| `test_det_wifi_ht.py` | 802.11n single stream: MCS 0–7 × both guard intervals byte-exact, A-MPDU de-aggregation, sensitivity, refusal of 2+ streams / 40 MHz / extension streams, **real MCS1 frames at 5240 MHz (8/11 FCS; rotating HT pilots)**, Viterbi rewrite equivalence, non-HT through the same reference transmitter (run directly for the MCS-vs-SNR table) | ~3 s (table: ~2 min) |
| `test_det_wifi_scaled.py` | 802.11 at non-standard clocks + LDPC/STBC + transmitter attribution: all LDPC codes and PPDU layouts; HT BCC/LDPC × STBC byte-exact, 40 ppm clock drift; synthetic 5 / 10 (mirrored) / 11 MHz (mirrored, MCS7 LDPC STBC) links with two stations fully attributed; no link from standard Wi-Fi, LTE-like OFDM or noise; attribution rules; real Microhard 11 MHz mirrored link and DJI 5 MHz link (also collected by pytest) | ~30 s |
| `test_det_droneid.py` | Remote ID: Open Drone ID (ASTM F3411) field by field, DJI DroneID cross-checked against the external reference decoder, refusal to report a guess, platform signatures, end to end through `mac80211.parse` | ~1 s |
| `test_det_devices.py` | device tracking and findings: AP/station roles, randomised-MAC linking, deauthentication floods, evil twin vs one radio's virtual BSSIDs, proportionate security findings, drone findings, plugin against the generic report contract | ~1 s |
| `test_det_jammer.py` | swept-carrier jammer: triangle / sawtooth / sine / noise-FM sweeps (1.2–50 kHz, 20–50 MHz, jittered + clipped beyond the band edge, over Wi-Fi traffic, down to +4 dB/MHz) with rate and extent; no detection on noise, CW, pulsed and pulsed-LFM radar, TDD OFDM, periodic / random hoppers, Wi-Fi traffic, FM video, a 4 MHz sweep, a +2 dB/MHz sweep; the two bench jammer captures, and the whole 10-block signal-vault session when present (also collected by pytest) | ~10 s |

`_data.py` holds the paths. The DroneID payload decoder is a separate repo
(`droneid-decoder-reference`, found automatically or via `$SPECTRUM_INTEL_DRONEID_SRC`).
