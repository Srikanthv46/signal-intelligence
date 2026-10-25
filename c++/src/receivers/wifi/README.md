# wifi (C++)

Port of `python/spectrum_intel/receivers/wifi/`. Plugins, in Python's order: `wifi_ofdm`,
`wifi_ofdm_scaled` (802.11 OFDM) and `wifi_dsss` (802.11b DSSS/CCK). Shared pieces: `mac80211`
(MAC dissector, OUI table), `droneid` (Remote ID / DJI DroneID), `links` (transmitter attribution).

## 802.11 OFDM (`wifi_ofdm`, `wifi_ofdm_scaled`)

| File | Python | What |
|---|---|---|
| `ofdm_rx.*` | `ofdm_rx.py` | `WifiDetector`: every standard 20 MHz channel mixed down (complex64) and resampled to 20 MSps; L-STF plateau (16-sample autocorrelation, 48-sample window), L-LTF timing / fine CFO / channel / per-subcarrier noise, L-SIG Viterbi + parity / tail / bit agreement, HT-SIG / VHT-SIG-A CRC-8; accept, merge over scan centres, dedupe 40/80/160 MHz duplicates. `viterbi`, `conv_encode`, `crc8_sig`, `channel_plan`, `limit_denominator` (Python `Fraction.limit_denominator`) |
| `payload.*` | `payload.py` | non-HT payloads (all 8 rates) and single-stream 20 MHz HT (MCS 0-7, BCC or LDPC, with or without Alamouti STBC): HT-LTF channel, smoothing candidates, pilot phase/slope tracking, log-MAP LLRs, LDPC rounds with data-aided channel re-estimation, descrambling, A-MPDU de-aggregation, FCS, MAC dissection |
| `ldpc.*` | `ldpc.py` | the 12 lifted 802.11 LDPC codes, PPDU codeword parameters (shortening/puncturing/repetition), normalised min-sum decoder (the encoder is only used by the Python test transmitter: not ported) |
| `ofdm_scaled.*` | `ofdm_scaled.py` | `find_bands` (60 kHz spectrogram, repeatable edges), `numerology` (CP correlation -> T_u, its periodicity -> T_s, guard, clock, occupied width), `scan` (time-scaled decode at the round clocks and the measured clock, plain and mirrored) |
| `ofdm_report.*` | `plugin_ofdm.py` | PPDU records with payloads, VHT Partial-AID attribution, the decode results of both plugins |

The plugins cache the records in `ctx.cache` (`wifi_ofdm`, `wifi_ofdm_scaled`); `links::attribute`
mutates them as the Python dicts are mutated. The undecoded-burst finder, SigMF annotation writer,
plot and CLI of the Python modules are not part of the pipeline and are not ported.

Arithmetic follows numpy (np::cmul / np::cdiv / np::cabs, pairwise sums, np.round vs Python round
per call site, phase ramps built as Python builds them). Two pieces are re-implemented locally for
speed, bit-identical to what they replace: `ofdm::detail::resample_cf` (= `sig::resample_poly` for
complex64, only the kept outputs, register accumulators; unit-tested bit for bit) and the LDPC check
update (padding edges handled analytically). Payloads of all PPDUs are decoded in parallel; the
channels of a scan and the (clock, orientation) trials of `wifi_ofdm_scaled` run in parallel with
results in input order.

Faithfully reproduced Python behaviour: `_qam_points` maps a 0 bit of BPSK/QPSK to 255 (uint8
wrap-around of `2*b - 1`); it only affects the data-aided channel of LDPC MCS 0-2, which no testdata
PPDU uses.

### Verification

- `tests/unit/test_rx_wifi_ofdm.cpp` (6 cases, 2032 assertions) against
  `tests/vectors/wifi/gen_wifi_ofdm.py`: Viterbi / encoder / CRC-8, limit_denominator, channel plans,
  resampler (bit-exact vs core), deinterleaver, depuncture, descrambler, log-MAP LLRs, channel
  smoothing, A-MPDU de-aggregation, LDPC parameters and decoding (a converging codeword, one that
  runs all 60 iterations, a two-codeword stream: same bits, same iteration counts); a synthetic
  30.72 MSps scene (non-HT 6/24/54, HT MCS7 BCC A-MPDU, HT MCS7 LDPC+STBC, HT MCS3 short GI, CFOs):
  every packet field and every record identical; a 10 MHz-clock link off the grid: bands,
  numerology (1e-9), link and records identical.
- Parity (`tools/parity.py --receivers wifi`): 16/16 inputs identical; wifi_ofdm 12/12,
  wifi_ofdm_scaled 3/3, wifi_dsss 14/14 signals identical, energy 127/127.

### Timing (receiver_s of the report; C++ `--threads 2`, Python from the parity cache)

| Input | Python wifi_ofdm / scaled (s) | C++ wifi_ofdm / scaled (s) |
|---|---|---|
| D_20251119-182450 (Phantom 4, 2.4 GHz) | 1.25 / 2.36 | 0.36 / 1.00 |
| D_20251119-183758 (2.4 GHz + 5 MHz link) | 0.95 / 3.37 | 0.28 / 1.19 |
| D_20260211-161016 (Microhard, 11 MHz) | – / 4.78 | – / 1.72 |
| D_20260211-161341 (Microhard, LDPC STBC) | – / 12.15 | – / 2.95 |
| D_20260916-171343 (5 GHz) | 0.59 / 2.39 | 0.15 / 0.93 |
| real_ht_mcs1_5240 | 0.26 / 0.29 | 0.05 / 0.27 |

wifi_ofdm is ~3-4x faster; wifi_ofdm_scaled 2-4x (its band search and numerology are FFT-bound:
4M-point transforms per candidate band, as in Python).

## 802.11b DSSS (`wifi_dsss`)

`dsss.*` (Barker / CCK receiver), `dsss_report.*` (`plugin_dsss.py` results): 14/14 signals
identical on testdata (`tests/unit/test_rx_wifi.cpp`).

## Drone links inside frames (2026-10-22, PLAN G1 / G2 / G7)
`mac80211.cpp`: wfb-ng / OpenIPC (`wfb_ng`), ESP-NOW and the ExpressLRS Backpack (`esp_now`, `elrs_backpack`, MSP v2
CRC-8), MAVLink over open-Wi-Fi UDP / TCP (`ipv4.mavlink`) and the matching `drone` flags, as in Python.
`tests/unit/test_drone_links.cpp` (vectors `tests/vectors/drone_links/gen_drone_links.py`): identical to Python.
