# cofdm — DVB-T-framed COFDM links (drone / UAV video downlinks)

**Research capability** (2026-10-21): built from ETSI EN 300 744 and vendor information
(`research/drone_transmitters/06_dvbt_cofdm.md`). Validated on synthetic signals only, and tagged `research` in
reports (`schema/capabilities.json`, technology `cofdm_video_link`).

**Detects:**
- COFDM video transmitters using DVB-T framing, on standard 8 / 7 / 6 / 5 MHz channels and on the scaled
  4 / 2.5 / 2 MHz modes of COFDM modules. Scaling: elementary period T = 7 / (8·B) µs.
- Users include professional UAV downlinks (Domo/DTC "DVB-T mode", IMT/Vislink 2K), consumer COFDM drone links
  (HawkEye), cheap "COFDM 2K" modules (200 MHz–4.4 GHz) and DVB-T broadcast.

**Decodes:** the TPS (transmission parameter signalling), BCH(67,53)-verified:
- mode (2K / 8K), guard interval, constellation (QPSK / 16 / 64-QAM), hierarchy, HP / LP code rate, frame number;
- the cell-ID byte (high byte in frames 1 / 3, low byte in frames 2 / 4).

Encryption (AES / ABS) scrambles only the transport stream, so TPS reads even on encrypted links. The video itself
is not decoded.

## How (`detector.py`)
1. **Candidates:** Welch spectrum runs above floor + 2 dB whose width is 0.951·B (±10 %).
2. **Numerology:** the band sliced at 8B/7 MHz; the CP metric folded on N(1+g) for 2K / 8K × GI 1/4, 1/8, 1/16, 1/32.
   Best peak ≥ 0.2 and ≥ 3× its median.
3. **Carrier:** fractional offset from the CP phase. Integer offset from the continual pilots (45 / 177, +2.5 dB), which
   must be ≥ 1.5 dB over the median shift and ≥ 1.2× the runner-up (3k shifts partly align the scattered-pilot raster).
4. **TPS:** Re{z_l z*_(l−1)} summed over the 17 / 68 TPS carriers (the channel cancels) gives s_l. Then the sync word
   or its inverse, then the BCH check.

## Synthetic results (`python/tests/test_cofdm.py`, generator `synth.py`)
- **2K at 8 MHz (GI 1/8, QPSK 1/2) and 2K at 7 MHz (GI 1/32, 16-QAM 7/8):** detected, carrier to ~4 Hz, TPS decoded
  with BCH ok, down to 0 dB in-band SNR.
- **2K at 2 / 4 MHz and 8K at 8 MHz:** numerology, bandwidth and carrier correct. A TPS frame (38–68 ms) is longer than
  a 34 ms block, so the label is "-like" with a note (stream several blocks to read TPS).
- **False alarms:** none on noise, and none on an LTE-like 4.5 MHz OFDM signal (width matches the 5 MHz window).
- **Timing:** 0.07–0.2 s per 2^21-sample block in Python (dominated by the capture FFT and Welch).

## DVB-T2 P1 (research, 2026-10-22, `dvbt2_p1.py`; EN 302 755 clause 9.8, `research/drone_transmitters/src_etsi/dvbt2.txt`)
- **Detection:** per candidate band, the slice at the elementary rate is searched for the P1 preamble. The C and B
  parts are frequency-shifted (f_SH = 1/1024T) copies of A, so y[m] y*[m+542] and y[m] y*[m−482], de-rotated, peak at
  the P1 start whatever the carrier offset. Threshold 0.35; noise, DVB-T and 1K CP-OFDM give ≤ 0.11, a P1 at −3 dB ~0.55.
- **Decode gate (regression 2026-10-22):** a peak counts as a P1 only if S1 and S2 both decode with a margin ≥ 0.4 and
  S1 is defined. Synthetic P1s give ≥ 0.74 even at −6 dB. On real testdata, ADS-B pulses, a weak analog VTX and
  clipped 940 / 899 MHz blocks passed the 0.35 correlation but decoded with ≤ 0.13 on S1 or S2 (six false
  "DVB-T2" links, all removed).
- **Decoding:**
  1. FFT of A; the integer carrier offset from the 384-carrier mask (±100 carriers).
  2. Descrambling (PRBS 1+X^14+X^15, init 100111001000110).
  3. DBPSK soft bits between consecutive active carriers.
  4. Correlation with the 8 S1 / 16 S2 patterns → `p1`: S1 (T2_SISO / T2_MISO / T2-Lite / FEF), S2 field 1 (FFT
     size + guard-interval class), mixed bit, pattern margins.
- **Labels:** "COFDM video link (DVB-T2 …)"; `dvbt2_like` when a pattern margin is ≤ 0.2.
- **Synthetic results:** matches the spec's worked example (MSS_SCR, scrambler); S1/S2 correct to −3 dB with a
  13-carrier offset; no false alarm on noise / DVB-T.
- **Limit:** a T2 frame lasts up to 250 ms, so a 34 ms block may hold no P1. T2 body symbols alone are not reported
  (their continual-pilot tables are not implemented, and a CP fold alone gave false links on real captures).

## Not yet
- **DVB-T2 L1 (P2) and body continual pilots.** Without a P1 in the block, DVB-T2 is `dvb_t2`-like by band rule only.
- **Proprietary narrowband COFDM** (Domo NB 625 kHz–2.5 MHz, UMVL; iWave LDPC): not reported. Links without the DVB-T
  continual pilots are dropped (2026-10-21: before this rule, the CP fold at DVB-T periods gave 5 false links on real
  captures, on DJI OcuSync and 2.4 GHz ISM traffic). The blind CP-OFDM scan in `digital_vtx` covers Tu ≤ 120 µs only.
- **TPS across blocks** (StreamAnalyzer) for the 2 / 2.5 / 4 MHz and 8K modes.
- **The scaling rule for cheap modules** is an inference: confirm on a real module.

## API
```python
from spectrum_intel.receivers.cofdm import detect
ems = detect(x, fs, fc)          # Emission family "cofdm", label dvbt | dvbt_like
```
