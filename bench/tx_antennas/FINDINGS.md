# Transmit antennas from one receive channel (bench, 2026-10-21)

In the library since 2026-10-21 as an evidence-only verdict, at the user's request: `receivers/digital_vtx/streams.py`, C++ `streams.cpp`, `decoded.fields.tx_antennas`. Not checked against Walksnail units of known antenna count.

Question: how many antennas does an unknown CP-OFDM transmitter use, seen by one Pluto channel? Asked to narrow
the Walksnail-like link in `testdata/live/LIVE_OFDM_VTX_5805MHz` to a model. Per the reported specs, the Mini 1S,
Mini 1S Lite and Nano V3 air units have one antenna and the Moonlight has two (the Pro and GT not confirmed).

## Tests (`streams.py`)
- **kappa**: per-subcarrier kurtosis over symbols, E|X|^4 / (E|X|^2)^2, after per-symbol power normalisation.
  A static one-tap channel cancels. One stream per subcarrier: QPSK 1.00, 16-QAM 1.32, 64-QAM 1.38,
  256-QAM 1.40 (CDD and antenna switching included). Two streams (STBC or spatial multiplexing): 1.45-1.65.
  Gaussian-like (DFT-spread): ~1.92. Overlap: 2-stream QPSK (~1.40) vs 1-stream 256-QAM.
- **delta** (Alamouti fingerprint): in an Alamouti pair |y1|^2 + |y2|^2 = (|h1|^2 + |h2|^2)(|a|^2 + |b|^2),
  so the pair's powers are anti-correlated. The statistic is the time-pair power correlation on the paired parity
  minus the other parity (the difference cancels the estimator's -1/(n-1) bias on short PPDUs).
  Synthetic STBC: -0.44 (equal branches) to -0.16 (2nd branch -15 dB). 1 TX, CDD, SM: 0. SFBC shows on the
  frequency pairs instead.

## Validation
Synthetic (`test_streams.py`, `test_pairs.py`, `test_imbalance.py`): 50 symbols x 1728 subcarriers, 15-31 dB SNR,
random multipath per antenna.

Real 802.11 PPDUs with SIG-stated modes (`wifi_reference.py`, raw FFT of the data symbols, pilots excluded):

| Capture | Mode (from SIG) | n | kappa | delta |
|---|---|---|---|---|
| D_20260211-161341 (Microhard) | HT MCS6/7 **STBC, 2 TX** | 12 | 1.52-1.61 (median 1.60) | **-0.33** |
| D_20260211-161341 | HT MCS15, **2 spatial streams** | 1 | 1.58 | ~0 |
| wifi/real_ht_mcs1_5240 | HT MCS1, 1 stream | 11 | 1.27 | ~0 (both parities -0.06/-0.14: short PPDUs) |
| real_ht_mcs1_5240, D_20260916, LIVE_OFDM_VTX_5805MHz | NON_HT, **FCS ok** | 8 | 1.07-1.20 | ~0 |
| same captures, LIVE_WIFI_5755MHz | "NON_HT" with **FCS failing** | 38 | 1.65-1.80 | ~0 |

The FCS-failing "NON_HT" rows are not legacy frames: they are PPDUs whose HT/VHT/HE SIG the detector did
not recognise. Only FCS-verified frames count as single-stream ground truth.

## Result on LIVE_OFDM_VTX_5805MHz (`link_streams.py --offset-hz 5.07e6`)
Frame (bench decomposition): ~50 full-band data symbols (1729 subcarriers) | ~0.7 ms gap | 7-8 symbols at 24 %
occupancy on the same grid (the 4.8 MHz return bursts) | STF (128-sample period) | LTF: ~4 x 1024-sample
repetitions, BPSK, odd comb only (Pluto integer CFO ~ +5 subcarriers), identical repetitions, nothing on the other
comb (-21 dB, incoherent: noise) | data.

| Run | kappa | time pairs even / odd | delta |
|---|---|---|---|
| 0.05 ms (truncated by the capture start) | 1.61 | -0.04 / **-0.20** | parity flipped, as expected |
| 4.23 ms | 1.60 | **-0.21** / -0.03 | -0.18 |
| 11.39 ms | 1.56 | **-0.23** / -0.04 | -0.19 |
| 18.55 ms | 1.55 | **-0.22** / -0.03 | -0.19 |
| 25.72 ms | 1.68 | **-0.13** / +0.01 | -0.14 |

**Two transmit antennas, Alamouti space-time code over OFDM symbol pairs that start at each frame's first data
symbol.** The training field is single-port (one comb, no orthogonal cover), so STBC applies to the data only.
That fits 802.16e/WiMAX-style matrix-A STC more than 802.11n-style per-stream training.

Open points:
- delta (-0.18) is weaker than pure STBC at this kappa (synthetic -0.34..-0.43). Possible causes: part of each
  frame not STBC-coded, or the second branch ~10 dB weaker. Not resolved.
- A blind Alamouti decode (h2/h1 by kurtosis minimisation per subcarrier, h1 from the LTF) did not converge:
  the per-pair phase drift and the LTF's antenna mapping are unknown. A decode would prove the scheme.
- No Walksnail reference recordings yet. Needed: a known single-antenna unit (Mini 1S / Nano V3) and a known
  dual-antenna unit (Moonlight / Pro / V2) on this bench. The test enters `python/` only after it separates them.
