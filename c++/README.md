# spectrum_intel — C++ library

The production build of `spectrum_intel` for edge devices (ARM64 / x86-64, Linux or macOS). It does what the
Python library in `../python` does, with the same results: an IQ block in, the JSON response of
`../schema/response.schema.json` out. Python is the reference; a capability enters this library only after it is
verified in Python, and only with a passing parity check (`tools/parity.py`). How to port: `PORTING.md`.

## Build

Needs CMake ≥ 3.16 and a C++17 compiler (GCC ≥ 9, Clang ≥ 10). There are no external dependencies: pocketfft
(BSD-3), nlohmann/json (MIT) and doctest (MIT) are vendored in `third_party/`.

```
cmake -S c++ -B c++/build -G Ninja -DCMAKE_BUILD_TYPE=Release      # add -DSI_NATIVE=ON on the target itself
ninja -C c++/build -j2
c++/build/si_tests                                                  # unit tests
c++/build/si_analyze testdata/DJI_Mini_2_20251119-165139            # response JSON on stdout
```

| Option | Default | |
|---|---|---|
| `SI_RECEIVERS` | all | receivers to build, e.g. `"lte;nr;gsm"` (`none`: energy pipeline only); a receiver left out is not in the binary |
| `SI_NATIVE` | OFF | `-mcpu=native` / `-march=native` for the build machine |
| `SI_SHARED` | OFF | `libspectrum_intel` as a shared library (C API for other languages) |

Cross-compiling for a Jetson / Raspberry Pi: a standard CMake toolchain file; nothing in the build runs target
code.

## Use

```cpp
#include <spectrum_intel/spectrum_intel.hpp>

si::Engine engine;                      // receivers are created once; analyse many blocks with it
si::IQBlock b;
b.data = samples;                       // borrowed for the call
b.n_samples = n;                        // complex samples
b.format = si::IqFormat::CI16;          // CF32, CF64, CI16, CI8, CU8; raw ADC counts are fine
b.sample_rate_hz = 61.44e6;
b.center_frequency_hz = 2.437e9;
b.bandwidth_hz = 56e6;                  // analog bandwidth (optional)
b.receiver_ppm = 7.4;                   // SDR frequency error, if known (optional, recommended)
b.start_time_utc = "2026-10-09T10:00:00.000000Z";
si::Json r = engine.analyze(b);         // never throws for bad input: see r["status"], r["errors"]
```

- **Small SDR buffers** (2048, 4096, ... samples): `si::StreamAnalyzer stream(engine);` then
  `for (auto& resp : stream.push(block)) ...`. It analyses every 2^21 samples of one tuning.
- **C API** (`include/spectrum_intel/c_api.h`): `si_engine_create`, `si_analyze` (returns the JSON text),
  `si_free`.
- **Command line:** `si_analyze PATH|--raw FILE ... [--enable/--disable R] [--ppm X] [--compact] [--report]`
  (header of `apps/si_analyze.cpp`).
- **Threads:** one block uses `SPECTRUM_INTEL_THREADS` threads (default min(4, cores)). Receivers and independent
  hypotheses run in parallel, and results do not depend on the thread count. Several engines may run in parallel
  threads.

## Layout

See `PORTING.md`. In short: `src/core` holds the numpy/scipy equivalents and shared records; `src/detect`,
`features`, `classify`, `report` and `analyze.cpp` form the pipeline; `src/receivers/<name>` has one folder per
receiver, mirroring the Python package; `src/api` holds the Engine and the response; `apps/` the CLI.

## Status (parity with Python on testdata)

| Part | Status |
|---|---|
| core numerics (numpy/scipy equivalents) | ported; FFTs bit-identical to numpy/scipy (`test_core_fft`), `ofdm15k` verified on D_20260916-171343 (`test_core_ofdm15k`), `bursts` (band / burst extent from in-band power, `test_core_bursts` + parity); the rest is verified per receiver |
| detection pipeline, features, band rules, report | **parity**: all 25 testdata inputs identical (`parity.py --receivers none`), including the receiver-spur and LO-harmonic-image flags (`core/spurs`); pulse-train features (PRI, staggered PRI, LFM) and the pulsed-radar / DME gates (`features/measure`, `classify/rules`, 2026-10-16); pulse-train merge (`detect/merge` `merge_pulse_trains`) and the short-pulse OFDM gate (`test_core_pulse_merge`) |
| all receivers together | **parity**: 47/47 testdata inputs identical with every receiver enabled (2026-10-22, after the coherent BLE path, frame-period R², ADS-B UAV flag, SECAM, analog audio; `si_tests` 91/91); earlier (`parity.py --refresh`, 2026-10-21, after the drone-transmitter research work: validation tag, colour standard, Remote ID BLE, SiK/MAVLink, catalogs, cofdm); `si_tests` 87/87 |
| API: Engine, response, StreamAnalyzer, C API, CLI | responses validate against the schema and equal Python's (except engine, timing, ids) |

Receivers: see the table below, filled in as each port passes `tools/parity.py`.

| Receiver | Status | Parity on testdata |
|---|---|---|
| wifi | ported: wifi_ofdm (802.11a/g/n/ac/ax preambles, non-HT + HT MCS 0-7 BCC/LDPC/STBC payloads), wifi_ofdm_scaled (5/10/11 MHz clocks, mirrored spectrum), wifi_dsss; MAC, OUI, DroneID, links, wfb-ng / ESP-NOW / ELRS Backpack / MAVLink-over-Wi-Fi in the MAC parser (2026-10-22, `vectors/drone_links`); windowed off-grid CP trigger (PLAN G13) and China GB 46750 Remote ID in `core/remote_id` (PLAN G14, `vectors/g10_g14`) 2026-10-23; batch W 2026-10-25 (`vectors/wifi_field`): HE-SIG-A (W2, `ofdm_rx.parse_he_sig_a`, scene + records + HE summary), drone hotspot SSIDs (W3), Tello SDK / video datagrams (W4, `mac80211`), BSS link profiles + drone-link verdicts (W1, `field.cpp`, links + end-to-end summary, label, evidence) == python | 28/28 identical with all receivers; 16/16 alone (`--receivers wifi`): wifi_ofdm 12/12, wifi_ofdm_scaled 3/3, wifi_dsss 14/14 |
| rc_link | ported (LoRa + PHY with the SFD / sync-word checks, GFSK/FSK, BLE incl. the advertising-channel label rule, GSM-burst gate, label check, 2 Mb/s ESB-like links with the joint CRC-16 check) | 30/30 (2026-10-16); 28/28 identical with all receivers; 25/25 alone (`--receivers rc_link`); 2026-10-17 ESB: `LIVE_ESB2M_2440MHz` 15/15 identical, full re-run pending; 2026-10-20 LoRa SF5 + LoRa-only 860–930 MHz: `LIVE_LORA_SF5_890MHz`, `LIVE_PMR_860.4MHz`, `LIVE_ELRS_2G4_2440MHz` identical; the FSK/BLE unit test's vectors (`ref.json`, Oct 9) predate the Oct 17 BLE change and need regenerating; 2026-10-21 coherent low-SNR BLE path (`gfsk_coherent`, research) == python on `vectors/ble_coherent` (identical bits, h, CFO, decode); 2026-10-22 BT5 Long Range (`ble_coded`, research) == python on `vectors/ble_coded` incl. the whole rc_link pipeline; 2026-10-22 FrSky D8 / D16 (`cc2500`, research) == python on `vectors/frsky` (labels, decoded packets, CRCs); blind FSK frame (`fsk_frame`, research) == python on `vectors/fsk_frame` (labels, frames); FlySky AFHDS 2A (`afhds2a`, research) == python on `vectors/afhds2a`; XN297 / Bayang (`xn297`, research) == python on `vectors/xn297`; IEEE 802.15.4 (`ieee802154`, research) == python on `vectors/ieee802154` incl. the pipeline; FLARM (`flarm`, research) == python on `vectors/flarm`; ELRS FHSS (`elrs` fhss / uid search / recover_hops) == python on `vectors/elrs_fhss`; model level (`core/drone_models`, `rc_link/models`, `digital_vtx` catalog, analog / dji / wifi plugins) == python on `vectors/drone_models`, `vectors/rc_models`, `vectors/digital_vtx`; 2026-10-22 RadioLink / AFHDS 1 / Hubsan / toy drones / ELRS GFSK / MAVLink over LoRa (research, PLAN G2-G6) == python on `vectors/rc_more`, `vectors/drone_links`; 2026-10-23 HoTT (CC2500 FEC Viterbi) / Hitec / S-FHSS (`hott`, `hitec`), ADS-L / OGN on the FLARM PHY (`adsl`), FANET in the LoRa identify() (research, PLAN G10 / G12) == python on `vectors/g10_g14` |
| dji | ported (video 15/30 kHz incl. further symbol grids, band / bursts from in-band power (`core/bursts`), RC uplink, DroneID detection; DroneID payload decoder not ported) | 28/28 identical (dji + digital_vtx), 2026-10-13 (bursts / hops); tmpdata 8/8 identical |
| analog_vtx | ported (detector, field demodulator incl. colour, low-SNR line clock / weighted discriminator / median / noise killer, PNG; >= 20 MHz band, capture-locked line-clock rejection) | 28/28 identical with all receivers; 16/16 alone (`--receivers analog_vtx,voice_fm`); 2026-10-23 3.3 GHz range + 3G3A / 3G3B channel groups (PLAN G9) |
| voice_fm | ported (channels, detector, audio, CTCSS/DCS, envelope and 4FSK-line rejection) | 28/28 identical with all receivers; 25/25 alone (`--receivers voice_fm`) |
| digital_vtx | ported (Herelink PSS fold / SSS -> PCI / subframe bursts, HDZero 27 MHz-mode pilot path `hdzero.cpp` 2026-10-22, blind CP-OFDM link detector and blind hop (RC uplink) detector `blind.cpp`, transmit-antenna verdict `streams.cpp` 2026-10-21) | 28/28 identical (dji + digital_vtx; digital_vtx 3/3), 2026-10-13 (bursts / hops); tmpdata 8/8 identical (D_20260609 ch0–3: O4 video + 2–3 O4 hop channels each); tx_antennas: LIVE_OFDM_VTX_5805MHz identical 2026-10-21 (full set not yet re-run); 2026-10-23 wider search bands + `band_context` + cellular band gate (`bands::digital_link_ranges`, `blind::SkipFn`; PLAN G8) == python on `vectors/digital_vtx` |
| cofdm | ported 2026-10-21 (DVB-T COFDM: candidates, CP fold, continual pilots, TPS + BCH; research capability) | synthetic vectors identical (tests/vectors/cofdm, 136 assertions); ~30 ms per 2^18-sample block |
| lte | ported (cell search, PSS/SSS with MRC, EARFCN raster, bursty candidates, ppm, flat-energy bandwidth) | 28/28 identical with all receivers; 25/25 alone (`--receivers lte`); real cells: B8 PCI 278/350 (SWEEP_317), PCI 351/169/415 on the live blocks; B3 PCI 101 (20 MHz) / 110 and B8 PCI 96 / 101 on LIVE_1840/940MHz |
| pmr_digital | ported (TETRA π/4-DQPSK + BSCH Viterbi/CRC, DMR SYNC/Golay/BPTC/CSBK CRC, NXDN frame sync) | 28/28 identical with all receivers; 25/25 alone (`--receivers pmr_digital`): 7 signals (TETRA x2 MCC 404 MNC 4040, DMR x2, NXDN) |
| nr | ported (SSB, PBCH, polar SC decoder from TS 38.212, MIB) | 28/28 identical with all receivers; 16/16 alone (n28 PCI 75 MIB, n78 PCI 82) |
| gsm | ported (FCCH searched over the CFO, block-wide receiver-ppm consensus, SCH search) | 28/28 identical with all receivers (24 GSM signals on live + cellular, incl. 14 on LIVE_940MHz); 16/16 alone (8 carriers at 928 MHz, before the CFO search) |
| adsb | ported (preamble/PPM/CRC-24 demod, every DF/typecode, Comm-B BDS 1,0–6,0 inference, CPR, per-aircraft tracks, met); encodings from Annex 10 / DO-260B / Doc 9871, pyModeS only as test oracle | 28/28 identical with all receivers; 16/16 alone (`--receivers adsb`): adsb 1/1 (real_1090_commb); 2026-10-23 UAT 978 (`uat.cpp`, plugin `uat`) + `core/rs` Reed-Solomon (research, PLAN G14) == python on `vectors/g10_g14` (RS cases, demod, dump978 sample messages) |
| jammer | ported (occupancy STFT, cyclic trigger, period by fold contrast, drift-corrected rate, sweep confirmation, band extent); 2026-10-20 | `LIVE_JAMMER_SWEEP_2460MHz` + `_jitter` identical with all receivers (jammer 2/2, rc_link 1/1); the 10 blocks of the signal-vault jamming session give the same sweep as Python; unit tests on synthetic sweeps / noise / CW; 40–50 ms per 2M-sample block; full re-run on testdata pending |
| emitters | ported (fingerprints, pairwise LLR, clustering, links, burst source) | 28/28 identical (all receivers) |
| devices | ported (frames from decoded.fields.ppdus, tracker, randomised-MAC linking, all threat rules, finalize_report signal ids); batch W 2026-10-25: sequence-counter continuity in the linking (W7), beacon flood and behaviour drone-link findings (W6 / W1) == python on `vectors/devices` case `w_batch` | 28/28 identical with all receivers; 16/16 alone (`--receivers devices`, and `--receivers devices,wifi`: devices/findings on 8 inputs) |
| navaids | ported 2026-10-25 (ILS glide slope / localizer DDM / SDM, VOR radial; research, signal gaps S9 / S21) | synthetic vectors identical (`vectors/signal_gaps/navaids`); no testdata yet; regression + parity pending (end of batch) |
| gnss | ported 2026-10-25 (GPS L1 C/A acquisition, C/N0, jamming / spoofing indicators, swept PPD; research, S3 / S12) | synthetic vectors identical (`vectors/signal_gaps/gnss*`); parity pending |
| broadcast_fm | ported 2026-10-25 (FM stations, stereo pilot, RDS; research, S20) | synthetic vectors identical (`vectors/signal_gaps/broadcast_fm`); parity pending |
| marine | ported 2026-10-25 (AIS GMSK / HDLC / CRC, message types 1/2/3, 4, 5, 18, 21, 24; research, S15) | synthetic + golden vectors identical (`vectors/signal_gaps/ais`); parity pending |
| radalt | ported 2026-10-25 (FMCW radio altimeters: ridge crossings, rate, shape, bandwidth, period; research, S8) | synthetic vectors identical (`vectors/signal_gaps/radalt*`, bit-exact); ~45 ms per 2^21 block on 2 threads; parity pending |
| ism_devices | ported 2026-10-25 (315 / 433 MHz: TPMS VDO / Toyota, EV1527 / PT2262, KeeLoq clear fields, unrecognised OOK / FSK; arbitrate; research, S6) | synthetic vectors identical (`vectors/signal_gaps/ism_*`, incl. ism_devices + rc_link on one scene); parity pending |

## Not in the C++ library (by design)

- Sweep tooling (sweep summaries, noise templates from neighbouring tunings, occupancy maps), PNG rendering and
  the research CLIs: development and offline tools, kept in Python.
- Synthetic transmitters (`synth.py`, `*_tx.py`): test-vector generators, kept in Python.

**Signal gaps batch (2026-10-25, research; `research/signal_gaps/PLAN.md`):** besides the receivers above, these
additions to existing ones are ported with Python-generated vectors in `tests/vectors/signal_gaps/`
(`tests/unit/test_signal_gaps.cpp`): adsb `modeac` (Mode A/C), voice_fm AM voice and `afsk` (APRS / AX.25 via
`core/hdlc`), rc_link `core/ble_trackers` and `bredr` (BR/EDR LAP), jammer `barrage` (plugin `jammer_barrage`),
`classify/rules` radar classes, `core/deinterleave` (`pulse_emitters` in `features/measure` pulse_train, and the
pulsed-radar evidence), and the 70–300 MHz band plan. The one expected failure until the end-of-batch run is
`test_rx_rc_link` (`vectors/rc_link/ref.json` predates the `ad.device` key of real adverts).

**Signal descriptions (schema 1.4.0, 2026-10-25):** `contract::description(label, decoder, flags)` reads the embedded
`schema/descriptions.json`; the report writer and the response carry `description`. `tests/unit/test_descriptions.cpp`
== python on every label of the latest regression plus one per variant (`vectors/descriptions`, 184 cases).
