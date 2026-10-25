# digital_vtx — Herelink, HDZero and any CP-OFDM video/data link (drones of any make)

**Detects**
- **HDZero (Divimath DM6300) digital FPV video, 27 MHz mode** (`hdzero.py`, measured 2026-10-22 on a user-confirmed
  HDZero Freestyle V2, `testdata/live/LIVE_HDZERO_5658MHz`): CP-OFDM, 1024 FFT at 27 MSps (Tu 37.926 µs, SCS
  26.367 kHz), CP 32 (1/32), 889 used subcarriers (−444…444) = 23.44 MHz, 26 continual pilots (carriers −444, −400,
  −360, −312, −280, −240, −224, −211, −200, −168, −136, −96, −56, 0, 56, 104, 128, 160, 184, 216, 264, 280, 312, 364,
  404, 444), scattered pilots every 32 carriers shifting 8 per symbol (4-symbol cycle), a ~17-symbol block of
  repeated preamble symbols once per video frame (426 symbols = 16.66 ms at 60 fps), one TX, continuous.
  - **How:** on each HDZero channel inside the capture (R1–R8, E1, F1, F2, F4, L1–L8) the band is cut from the block
    FFT at nw/1024 × 27 MSps, where `slice_plan` picks a fast (13-smooth) slice length whose per-symbol FFT size nw is
    an integer (2^n blocks: exactly 27 MSps; signal-vault's 2,000,000-sample blocks: 894,348 → nw 1042, 15 ms instead
    of a 0.21 s prime-length IFFT). The CP metric folded on the fractional symbol grid gives the timing and the
    fractional CFO; each symbol is transformed with its half-sample start offset corrected (otherwise the edge carriers
    lose 1.3 rad). Over the integer carrier offsets within ±25 ppm the 26 continual pilots must be coherent symbol to
    symbol: mean coherence ≥ 0.2 and ≥ 4× that of the data carriers ±2 beside them (real link 0.95 vs 0.09 at 14 dB;
    Wi-Fi 0.30 vs 0.25; Walksnail-like / DJI / noise ≤ 0.14). Preamble = symbols that repeat one of the next four
    (|ρ| > 0.3; data ~0.035); two preambles give the frame period.
  - **Result:** label `HDZero digital FPV video link` (0.9, partial): channel, carrier, CFO / ppm, frame period and
    fps, SNR, `model` (every HDZero VTX, unresolved: identical on air apart from power). Sensitivity on the real
    capture with added noise: found to about −6 dB link SNR (32 ms), frame rate to −3 dB; 65 k-sample blocks (1 ms)
    are found in 3 ms. Cost ~0.1 s per in-band channel on a 2 M-sample block. Blind links of the same Tu over the
    HDZero band are left to it. The 17 MHz narrow mode (Nano 90) and Shark Byte are not measured and not searched.
- **CubePilot Herelink:** an LTE Rel-8 downlink PHY used in the 2.4 GHz ISM band (and searched at 5 GHz), with a 5 ms TDD rhythm.
- **Any hopping narrowband CP-OFDM link (RC uplinks), found blind** (`blind.py` `blind_hops`): DJI O4 (30 kHz,
  4.4 MHz × 0.41 ms hops), OcuSync 2 / Autel (15 kHz, 1.1 MHz), or a hop train at any numerology found in the capture.
- **Any CP-OFDM link with Tu ≥ 10 µs (SCS ≤ 100 kHz), found blind** (`blind.py`): no vendor, numerology or band
  knowledge. Every drone video link in `testdata/` and `tmpdata/` is found this way, including links no protocol
  receiver knows. 802.11-family numerologies (Tu 3.2 / 6.4 / 12.8 µs, 5.82 µs) are left to the Wi-Fi receivers.

**Decodes** Herelink's physical-layer identity: carrier, PCI from PSS + SSS, CFO and clock ppm, and the
control / video / uplink subframes. The payload is proprietary or encrypted, so the status is `partial`. A blind
link is measured, not decoded (`partial`, label `-like` or `unknown`).

## How
- **Herelink:** an LTE PSS (ZC roots 25/29/34) fold on the 5 ms half-frame across the band, then a carrier search,
  SSS decode → PCI, and CFO from the CP phase. Subframe activity separates control / video / uplink.
  - **PSS purity check:** the PSS correlation is compared with what a pure LTE PSS gives at the measured SNR,
    and a link is rejected below 0.5. DJI OcuSync's sync symbol also repeats every 5 ms and partly matches the PSS.
  - Measured purity: real Herelink stays at 0.62 or more down to its sensitivity limit; DJI scores 0.10–0.35.
    At the old limit of 0.35, a far 10 MHz Mini 2 (D_20260916-171218, purity 0.351) was labelled Herelink.
- **Blind CP-OFDM links** (`blind.py`). Every OFDM symbol ends with a copy of its first Tu seconds; that is found
  without knowing the vendor, and below the noise floor:
  1. **scan:** each 10 / 20 MHz slice (half-overlapping; one slice for a narrow capture) is time-AGC'd (20 µs
     envelope, so strong Wi-Fi bursts weigh no more than the quiet time between them) and its autocorrelation over
     lags 8–120 µs is summed coherently over the whole capture (one FFT); peaks over the detrended lag profile
     (z > 7) are candidate Tu.
  2. **confirm:** a symbol-rate line in the cyclic spectrum of the CP metric (z ≥ 10; noise / Wi-Fi ≤ 5), or a
     sharp per-1 ms-window CP fold in ≥ 3 windows after the period refinement (strong LTE-like links have many
     slot-rate lines that hide the symbol line: Herelink scores line z 6.6 but 19 windows). Tu is refined to a
     fraction of a sample around the ACF maximum, and every later slice is taken at a rate where Tu is an integer
     number of samples: a one-sample lag error decorrelates the band edges of an 18 MHz link.
  3. **band:** the CP cross-spectrum (per symbol, FFT of the CP and of the same window Tu later, summed
     coherently) holds only what repeats at lag Tu: overlapping Wi-Fi, RC hops at another timing and a mirrored IQ
     image 17 dB down do not widen it. Within 6 dB of the in-band level, gaps < 1 MHz closed; the slice is widened
     and re-centred while the band touches its edge.
  4. **bursts:** the symbol period is refined from the window-to-window drift of the grid phase (a 0.17 % period
     error smears a 2 ms fold), then the grid phase is found per 2 ms window (1 ms hop), so a link that re-times
     mid-capture keeps its bursts (D_20260609-151042: the grid jumps by 10.3 µs at 16.8 ms). Per-symbol CP
     evidence → runs → the TDD frame (activity ACF, then the frame-fold contrast) completes weak bursts frame by
     frame (as the `dji` receiver). The rate is ≥ 16 samples per CP (a short CP on a narrow link is too noisy).
  4b. **band and bursts from in-band power** (`core/bursts.py`, 2026-10-13): the band is refitted as the best
     contiguous rectangle of the per-bin on-vs-off power ratio (medians over the link's burst / idle frames), and the
     bursts become the in-band power runs (per 16.7 µs frame, the 25th-percentile rank of the sub-band powers over each
     sub-band's floor; a quarter of the sub-bands next to a lit guard band left out; a much stronger signal over the
     band is "unknown" and bridged, as are dips of ≤ 2 symbols) that the CP runs touch; CP runs where the power says
     clearly off are dropped or trimmed. D_20260609-151042 ch0: band 5761.5–5779.6 MHz (the CP cross-spectrum gave
     5761.3–5782.2), 13 whole bursts.
  5. **gates:** Tu ≥ 10 µs, ≥ 40 symbols in bursts and one burst ≥ 10 symbols (an FM video's line rate folds onto
     a short "Tu" with 1–8-symbol blips); duplicates (same Tu, bands overlapping > 50 %) and IQ images (same Tu,
     mirrored about fc, ≥ 10 dB weaker) are dropped.
  - **numerology table** (`KNOWN_NUMEROLOGIES`): what the link resembles and its role — `drone` (label
    `OFDM drone video link-like (… numerology)`, confidence 0.6), `infrastructure` (`OFDM fixed-wireless
    link-like`, 0.6) or unknown (`unknown OFDM (TDD) link`, 0.45).

- **Hopping links** (`blind_hops`). An RC uplink hop (O4: a one-symbol preamble, an 8 µs gap, 10 symbols) is too
  short and too spread in frequency for the capture-level scan:
  1. numerologies: every blind link's, then the drone families' (15 kHz and 30 kHz, normal CP) — a drone's uplink
     shares its video numerology;
  2. scan (`core.ofdm15k_scan`): slices of ~75 subcarriers (≥ 1 MHz), CP metric folded in 0.25 ms windows, z > 5;
  3. each hit: edges from the in-band envelope (10 µs blocks, ±1 ms), band from the on-vs-off spectrum (rectangle
     fit over ±1.5 W, re-centred once; rejected if still at the window edge — a slice of a wider signal), then the CP
     fold z at that numerology over the burst must exceed 6; 4 symbols to 1 ms, 0.5–8 MHz, ≥ 3 dB contrast; not
     inside a burst of the link with that numerology;
  4. the hop cadence (the longest period most hop-start differences are multiples of, ≥ 3 pairs) predicts the
     missed hops, re-tested at a lower scan threshold with the same verification;
  5. one numerology per hop (stronger fold), the hops of a link share dwell (±20 %) and width (±30 %), ≥ 2 hops;
     hops whose bands overlap form a channel (one report signal each).
  Labels: `OFDM hopping RC uplink-like (… numerology)` (0.6; `RC uplink-like` without "hopping" on one channel),
  `OFDM fixed-wireless uplink-like` for an infrastructure numerology, else `unknown OFDM hopping / narrowband burst
  link` (0.45). A hop the `dji` uplink tracker covers, or that sits inside a `dji` video link's bursts and band, is
  that receiver's (`explained_receiver_detections`).

## Drone link waveforms (what the blind detector sees)
Measured here on real captures (blind, no protocol receiver), plus public information:

| Link | SCS / Tu | CP | Band | Frame / bursts | Seen on |
|---|---|---|---|---|---|
| DJI OcuSync 2 (Mini 2) video | 15 kHz / 66.67 µs (LTE-like) | 4.7 µs (1/14, LTE normal) | 9 / 18 MHz | TDD, bursts of 1–5 ms, ZC sync | DJI_Mini_2, D_20260515-*, test_signal 2484 MHz |
| DJI O4 (Mini 4 Pro) video | 30 kHz / 33.33 µs | 2.34 µs (1/14) | 9 / 18 MHz (40 MHz mode exists) | 2.5 ms frame, ~1.7 ms bursts; comb-4 repeated reference; can re-time mid-capture | DJI_Mini_Pro_4, D_20260609-151042, test_signal ch0 |
| HDZero (Freestyle V2) video, 27 MHz mode | 26.367 kHz / 37.926 µs (1024 at 27 MSps) | 1.185 µs (1/32): too short for the blind scan, found by the `hdzero` pilot path | 23.44 MHz (889 sc) | continuous; per-frame preamble every 16.66 ms (60 fps) | LIVE_HDZERO_5658MHz |
| DJI Lightbridge 2 (Phantom 4) video — inferred | 10.94 kHz / 91.43 µs (802.16e: 1024 FFT at 11.2 MSps) | 11.4 µs (1/8) | 10.3 MHz | long bursts (4–10 ms), 61 % duty | DJI_Phantom_4 (2426.1 MHz true) |
| DJI OcuSync 2 RC uplink | 15 kHz, 72 subcarriers | LTE | 1.1 MHz, hopping | 0.51 ms hops | `dji` receiver; blind: `blind_hops` (D_20260916, 4 of 5) |
| DJI O4 RC uplink | 30 kHz / 33.33 µs | 2.34 µs (72/1024) | 4.4 MHz, hopping (5747.5 / 5762.4 / 5792.5 MHz seen) | 0.41 ms hops (1 preamble symbol, 8 µs gap, 10 symbols), 2.5 ms cadence, TDD with the video | `blind_hops` (D_20260609-151042 ch0–3, test_signal 2484 MHz at ~2510 MHz) |
| DJI Wi-Fi link (Spark / Mavic Air class) | 802.11a at a 5 MHz clock: 78.125 kHz / 12.8 µs | 1/4 | 4.4 MHz | CSMA | `wifi_ofdm_scaled` (D_20251119-183758) |
| Herelink | 15 kHz (LTE Rel-8 DL) | LTE normal | 18 MHz | 5 ms TDD | Herelink_20251125 |
| Microhard pDDL / data links | 802.11n at an 11 MHz clock: 171.9 kHz / 5.82 µs | 1/4 | 9.8 MHz | CSMA | `wifi_ofdm_scaled` (D_20260211-*) |
| Walksnail Avatar (public; matched live, device unconfirmed; data 2-TX STBC-like) | 10.94 kHz (802.16e 20 MHz, 2048 FFT) | 1/8 | 18.9 MHz (1729 subcarriers, DC used) | TDD 6.67 ms frame, 5.8–6.4 ms bursts; 4.8 MHz BPSK-like uplink bursts (0.94 ms) in the gaps | `testdata/live/LIVE_OFDM_VTX_5805MHz` (5805.07 MHz) |
| Autel SkyLink (public) | 15 kHz | LTE | 9 MHz video, 1.1 MHz uplink | – | no sample |
| Cambium PMP 450 (fixed wireless, not a drone) | 43.75 kHz / 22.86 µs | 1/16 | 18.5 MHz | 2.5 ms TDD, Tu/3 preamble | D_20260515-192148 / -192835 (5765 MHz) |

Common to every drone video link here: CP-OFDM, 1–40 MHz, a TDD frame (2.5–10 ms) with long bursts, in the 2.4 /
5.1–5.8 GHz ISM bands. That is what the blind detector keys on; the numerology only names the family.

- **Transmit antennas** (`streams.py`, `tx_streams`; bench: `bench/tx_antennas/FINDINGS.md`). On links with
  SNR ≥ 20 dB and ≥ 200 subcarriers, it reuses the link's slice from `blind_links` (`lk['slice']`, rate N / Tu).
  Per burst it folds the symbol grid from the CP metric and takes the fractional CFO from the CP phase.
  Symbols with CP correlation < 0.5 or < 75 % occupancy (training, narrowband bursts) are dropped. Then, per run of
  ≥ 20 contiguous symbols on the central 90 % of the band:
  - **kurtosis**: per-subcarrier E|X|⁴/(E|X|²)², after per-symbol normalisation. One stream (QPSK … 256-QAM, CDD,
    switching) ≤ 1.40; two streams ~1.45–1.65; Gaussian-like ~1.9.
  - **pair delta** (Alamouti fingerprint): the power correlation of symbol pairs (2i, 2i+1) minus (2i+1, 2i+2), counted
    from the run start. The same over subcarrier pairs, per block between null subcarriers, detects SFBC.
  Verdicts in `decoded.fields.tx_antennas`:
  - `2tx_stbc_like` / `2tx_sfbc_like`: |delta| ≥ 0.10 with the same sign on every run the capture does not cut off
    (`tx_antennas: 2`, label suffix `2-TX STBC-like`);
  - `1_stream_like`: kurtosis ≤ 1.43 and no pairing;
  - otherwise `unknown`, with a note.
  About 50 ms per link in Python.

- **System candidates** (`catalog.py`; **research**: `research/drone_transmitters/02_digital_video.md`, synthetic /
  catalog tests only). Each blind link gets a ranked `system_candidates` list (top 3), and the evidence names the first.
  - **Inputs:** Tu (a known Tu must match, ±1.5 %), CP fraction, nearest bandwidth mode (±12 %), band, channel raster,
    antennas from `tx_antennas` (2 TX vs a 1-antenna system is a conflict), and frame period.
  - **Rasters:** DJI x.5 MHz; the DJI FPV / Walksnail channel set; Raceband 37 MHz steps; AR8030 extra centres.
  - **Frame period** only ever adds a match: the blind frame estimate can be off (Walksnail: 6.673 ms reported, the
    data runs are 7.163 ms apart).
  - **Systems:** DJI O4 (2T) / O4 Air Unit Lite (1T), O3, OcuSync 2, OcuSync 1, Herelink, Autel SkyLink, Walksnail
    dual- / single-antenna air units, the AR8030 family, Lightbridge 2, HDZero / Shark Byte (numerology not public:
    needs BW + raster), DJI O1, Cambium PMP 450.
  - **Example:** the 5805 MHz bench link ranks "Walksnail Avatar HD (dual-antenna air unit)" first (6 matches, 0
    conflicts).

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.digital_vtx import detect
ems = detect(x, fs, fc, {"report_wifi": False})   # Emission families: herelink | digital_vtx
```

## Pipeline plugin
- Runs on captures in 2.33–2.53 or 4.9–6.0 GHz.
- Herelink is reported only at 2400 MHz and above. Below that, the same PHY is LTE band 40 and is left to `lte`.
- A blind link inside a licensed cellular band is left to the `lte` / `nr` receivers.
- **Arbitration:** a blind link covered over more than half its band by a DJI video link or Herelink of the same
  subcarrier spacing is that receiver's signal: it is dropped and recorded on the host
  (`detection.explained_receiver_detections`), corroborating it.

## Validation
`tests/receivers/test_det_digital_vtx.py`:
- the Herelink carrier and PCI 411 are found; no Herelink false alarm on the DJI Mini 2, Mini 4 Pro, Phantom 4,
  D_20260916 Wi-Fi captures or pure noise;
- **blind links** (no protocol receiver): Herelink (lte_like, Tu 66.65 µs, 17.9 MHz at 2400.70 MHz), DJI Mini 2
  (lte_like, 18.1 MHz at 5770.61 MHz), Mini 4 Pro (dji_o4_like, 9.1 MHz at 5770.58 MHz, frame 2.486 ms), Phantom 4
  (lightbridge_like, CP 1/8, 10.3 MHz); none in noise or in the D_20260916 Wi-Fi + RC-uplink capture.
- On `tmpdata/`: D_20260609-151042 ch0–3 (O4, 18 MHz, ch0 at about −5 dB in-band: 14 bursts, 2 symbol grids),
  D_20260515-192148 / -192835 (OcuSync 2, 18.3 MHz; PMP 450-like link on -192835), test_signal 2484 MHz ch3
  (OcuSync 2 10 MHz at 2466.46 MHz) and ch0 (O4 10 MHz at 2467.5 MHz).
- **hops** (`blind_hops_check`): D_20260916 OcuSync 2 uplink at 5.589 / 9.591 / 15.589 / 19.589 ms (the 25.59 ms hop
  next to the LO tone is not found blind); D_20260609-151042 ch0 O4 uplink at 1.37 / 3.87 / 8.88 / 17.08 / 24.57 ms
  (plus 6.38 ms); no hops on noise, Mini 2, Mini 4 Pro, Phantom 4.
- Sensitivity (added noise): Mini 4 Pro and Mini 2 found to +5 dB, lost at +10 dB — the same limit as the `dji`
  receiver.
- **tx_antennas** (`python/tests/test_tx_streams.py`): synthetic time-domain CP-OFDM links with 2×6-tap channels
  and a fractional CFO give one stream for 1 TX at QPSK / 64-QAM and for CDD; `2tx_stbc_like` for STBC at 16 / 64-QAM
  and with the second branch at −10.6 dB; `2tx_sfbc_like` for SFBC (DC null); `unknown` for spatial multiplexing
  and below 20 dB SNR. `LIVE_OFDM_VTX_5805MHz`: `2tx_stbc_like` (kurtosis 1.55, delta −0.18, 4 runs), as on the bench.
- `bench/regress.py` on all testdata: the only change is the new Phantom 4 link; no false link in the live /
  sweep captures (LTE / NR carriers are found blind but are outside the plugin's range or in cellular bands).

The Herelink capture's frequency label is 10 MHz low; the pipeline corrects it from BLE adverts, giving a true
carrier of 2410.6 MHz.

## Model level (`fields.model`, research, PLAN M3 / M4)
`catalog.model_candidates` turns the ranked `system_candidates` into the common model block
(`core/drone_models.model_block`). It takes the best-scoring systems with no conflict and their published model
lists (`catalog.MODELS`, from `research/drone_transmitters/11_model_level_identification.md`). The DJI O3 / O4
channel centres narrow the generation, and the HDZero 17 MHz mode implies the RunCam Nano 90 camera. A 2-TX
measurement already conflicts the single-antenna systems. The block names one model only when exactly one fits.
Otherwise it lists the candidates with `unresolved_reason`. Synthetic and table tests only:
`tests/test_drone_models.py`, and C++ == Python on the `vectors/digital_vtx` catalog cases.

## Wider search bands (2026-10-23, PLAN G8; research)
- `plugin.in_range` uses `core.bands_india.DIGITAL_LINK_RANGES`: 1.35–1.444, 1.6–2.53 and 4.4–6.0 GHz. A blind link
  outside the ISM / FPV bands carries `band_context` (1.4 GHz UAV band 1430–1444 MHz, L-band 1.35–1.44 GHz, 1.6–2.33 GHz
  and 4.4–4.9 GHz government / professional bands) and the label suffix "in the … band".
- Cheap band gate (`blind.scan` / `blind_links` / `blind_hops` `skip`): a slice, candidate or hop hit whose band and IQ
  mirror are both inside an LTE / NR band is not processed (the plugin dropped such links anyway); the hop scan is
  skipped only when the whole usable band is cellular, because its threshold is a whole-grid statistic. On
  LIVE_1840MHz digital_vtx fell from 3.8 s to 0.2 s (Python) and ~2.6 s to ~0.65 s (C++).
