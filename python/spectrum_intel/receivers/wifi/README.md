# wifi — IEEE 802.11 OFDM (a/g/n/ac/ax) and 802.11b DSSS/CCK

**Detects / decodes**
- **OFDM** (`ofdm_rx.py`, `payload.py`):
  - every 802.11 legacy preamble on the 20 MHz channel plan;
  - L-SIG (rate and length, with parity and bit-agreement check), HT-SIG / VHT-SIG-A (CRC), HE-SIG-A (CRC-4,
    802.11ax; see below);
  - **non-HT MPDUs:** demodulated, deinterleaved, Viterbi-decoded and descrambled, then checked with the CRC-32
    FCS; the MAC frame is dissected (`mac80211.py`).
  - **HT (802.11n) MPDUs, MCS 0–7, 20 MHz, BCC or LDPC, with or without STBC:** the single-stream set one
    antenna can receive. LDPC (`ldpc.py`: the 12 codes, the PPDU codeword layout, min-sum decoding) and
    Alamouti STBC from two transmit antennas are decoded. An A-MPDU is split on its delimiters and every MPDU
    checked on its own FCS. 2+ spatial streams and 40 MHz payloads are declined **with the reason**, not guessed.
- **802.11 at a non-standard clock, frequency or spectrum orientation** (`ofdm_scaled.py`, any band): 5 and
  10 MHz channels, vendor links that run an 802.11 chipset on another clock and convert it out of the Wi-Fi
  bands (Microhard at 11 MHz in 2025–2110 MHz, DJI Wi-Fi links at 5 MHz), off-grid 20 MHz channels (MikroTik),
  and mirrored spectra. The numerology is measured blind, then the capture is time-scaled into a standard
  channel and everything above applies unchanged.
- **802.11 at the standard clock but off the channel grid or outside the Wi-Fi bands** (`ofdm_rx.offgrid_plan`,
  2026-10-16): fixed-wireless backhaul that tunes in 1–5 MHz steps (Mimosa at 4880 MHz, and at 5755 / 5790 MHz
  inside the 5.8 GHz band). After the standard-channel scan, 20 MHz slices every 5 MHz are taken from one FFT of
  the block with the decoded PPDUs blanked; a slice "triggers" when its correlation at the 3.2 us symbol lag is
  >= 0.02 and >= 4x the median at 8 reference lags (off the 0.8 us STF period). A run of >= 3 adjacent triggered
  slices (a 20 MHz signal lights every slice within +-10 MHz) puts blind 1 MHz centres at its midpoint +- 5 MHz,
  skipping the standard centres; their PPDUs are merged with the standard ones (same start within 2 us and
  overlapping band = the same PPDU). Cost: ~1 s per >= 20 MSps block for the trigger (Python), plus ~0.12 s per
  blind centre where it fires. CW lines, LTE/NR (66.7 / 33.3 us symbols), DJI OFDM, analog video and radar do not
  trigger (bench: ratio <= 2.5). Live: `testdata/live/LIVE_MIMOSA_4880MHz` (was 2–3 `unknown` pieces; now 83 PPDUs,
  78 FCS-valid, BSSID 20:b5:c6:46:e8:3e).
- **Transmitter of every PPDU** (`links.py`): from the MAC header, from control responses (ACK / Block Ack /
  CTS name the station they answer), from per-station scrambler-seed counters, and as a weak hint from
  received power — so even PPDUs whose payload is out of reach (2-stream MIMO, failed LDPC) are attributed.
  - **VHT (802.11ac) headers, when the payload is out of reach:** VHT-SIG-A is on the primary 20 MHz,
    so even an 80 MHz PPDU gives its **direction** (Group ID 0 = station to AP, 63 = AP to station,
    1–62 = multi-user MIMO), spatial streams, MCS, coding and airtime. Each signal reports
    `vht_airtime_ms_by_direction` and the exact reason each payload is not decoded
    (`payload_not_decodable`). An uplink is tied to its AP through the Partial AID only when the
    value is non-zero and matches exactly one decoded BSSID (`bssid_from_partial_aid`) — see Limitations.
- **DSSS/CCK** (`dsss.py`): 802.11b at 1/2/5.5/11 Mb/s, long and short preamble; PLCP header CRC-16; PSDU; FCS;
  full MAC dissection.
- **Remote ID** (`droneid.py`), read out of any vendor element:
  - **ASTM F3411-22a / ASD-STAN prEN 4709-002 ("Open Drone ID")**, OUI `fa:0b:bc` type `0x0d`: UAS serial
    or registration, aircraft type, flight status, position, altitude, height, speed, track, accuracy,
    the operator ID and **the operator's own position**;
  - **DJI's proprietary DroneID**: the 91-byte payload, *located by its CRC-16/Kermit* rather than by an
    assumed offset, so a parse means the checksum passed. Gives serial, model, position, height, home
    point and operator position;
  - platform signatures from the IEEE registry and the SSID (DJI, Parrot, Autel Robotics, Skydio, Yuneec,
    Hubsan, PowerVision, OpenHD/wfb-ng), marked "-like" because they identify equipment, not a flight.
- **MAC dissector** (`mac80211.py`), shared by both:
  - frame control, address roles (RA/TA/SA/DA/BSSID, direction), sequence, QoS;
  - beacons and probe responses: TSF (gives the AP's uptime), interval and capabilities;
  - elements: SSID (hidden or wildcard), rates, DS channel, TIM, country, ERP, BSS load, RSN/WPA (gives a
    security summary such as "WPA2/WPA3-Personal (transition) CCMP-128"), HT/VHT/HE/EHT capabilities and
    operation (gives the PHY generation, 802.11b/g/n/ac/ax/be), 802.11k/r/v/u, WMM;
  - device identity: WPS (device name, manufacturer, model, serial) and Wi-Fi Direct device name;
  - auth/assoc status, deauth/disassoc reason, action category, Block Ack;
  - data frames: CCMP/TKIP/WEP header (PN, key ID); when unencrypted, LLC/SNAP → EAPOL (4-way handshake
    message number, EAP identity), ARP, IPv4/IPv6 + ports;
  - vendor of every address from the IEEE registry (`data/oui.tsv.gz`, 58k MA-L/MA-M/MA-S blocks, via
    Wireshark's `manuf`, 2026-10-05); randomised (locally administered) addresses are flagged;
  - `remote_id` and `drone` when a frame carries them (see `droneid.py` above);
  - a one-line `summary` per frame.

Status is `success` when an FCS passes. Payloads out of reach (VHT/HE, 2+ spatial streams, 40/80 MHz) leave
the signal `partial`.
For 802.11b, Barker-11 energy without a decodable PLCP header is reported as `partial`:
"Wi-Fi 802.11b DSSS-like (preamble before capture start)" or "Barker-11 DSSS (no 802.11b PLCP header)-like",
with the coherence as evidence.

## How — 802.11b (`dsss.py`)
1. **Channelise** every 2.4 GHz channel to 22 MSps (2 samples per chip, detection) and to 44 MSps (4 samples per
   chip, decoding).
2. **Detect:**
   - compute the Barker-11 matched filter, then |·|² summed over 16 symbols at each of the 22 sample phases;
   - DSSS is where the best phase stands out: peak/mean > 2.8 (noise stays below 2.6 over 180 ms), at any level;
   - a run splits where its best phase moves, so a second PPDU gets its own acquisition.
3. **Acquire SFD-first:**
   - Barker symbols on a fixed 1 µs grid at each of the 44 sample phases (drift over a preamble < 0.3 sample);
   - DBPSK with a sliding CFO estimate, descramble (x⁷+x⁴+1), then search for the SFD: long 0xF3A0 or short
     0x05CF, at most 1 bit error, after a SYNC run;
   - merge the hits of one SFD at neighbouring phases (main lobe and sidelobes).
4. **Track** each PPDU from its own SYNC start, phase and amplitude. The early/late loop is second order and
   interpolated, and updates only on symbols with at least 40% of the PPDU's amplitude, so noise, CCK chips and
   other PPDUs cannot pull it.
5. **Train:**
   - rebuild the scrambled SYNC+SFD exactly (best of the 128 scrambler states);
   - this gives a data-aided CFO, phase and SNR (in the 22 MHz channel: `snr_db`).
6. **PLCP header:** coherent and phase-tracked, with a differential fallback; checked with its CRC-16.
7. **PSDU:**
   - 1/2 Mb/s: coherent Barker symbols (differential fallback);
   - 5.5/11 Mb/s: chips on a grid whose rate follows the CFO (locked clocks; a fitted clock is the fallback),
     through a 9-tap LS equaliser trained on the preamble and header chips, then ML codeword and φ1 decisions
     with phase tracking.
8. **FCS:**
   - check the CRC-32;
   - on failure, FCS-guided correction: CRC-32 and the scrambler are linear, so ≤ 4 of the 20 least reliable
     symbol decisions are flipped, with singles and pairs matched by syndrome (about 2·10⁴ patterns, so a
     false pass has probability about 5·10⁻⁶);
   - still failing: STFT Wiener excision of partial-band interference (BLE, an adjacent OFDM edge), then decode
     again.
   - `symbols_corrected` and `psdu_demod` record what was used.
9. **Dedup:**
   - a PPDU also correlates on neighbouring channels, so the decode with the highest Barker amplitude wins;
   - a beacon's DS parameter set names the channel (`channel_source`).

## How — 802.11n (`payload.decode_ht_payload`)

The legacy preamble is detected exactly as before; HT-SIG gives MCS, length, guard interval,
aggregation, STBC, coding and the extension-stream count. The data field then differs from non-HT in
five structural ways, all of them from Clause 19:

1. **56 subcarriers, not 52** (−28..28), so there are 52 data carriers and `N_CBPS = 52·N_BPSC`.
2. **The channel is re-estimated from the HT-LTF**, the only field that trains the four extra
   subcarriers the legacy preamble never touches. Per-subcarrier noise still comes from the L-LTF pair
   (the only place two identical symbols exist); the four HT-only carriers take the legacy median.
3. **The interleaver uses `N_COL = 13`**, not 16. (The third HT permutation is a frequency rotation
   between spatial streams, so at one stream it is the identity.)
4. **Pilots start at `p_3`** — L-SIG, HT-SIG1 and HT-SIG2 consumed `p_0..p_2`.
5. **A short guard interval** makes the symbol 72 samples instead of 80, and **MCS7 needs the rate-5/6
   puncturing pattern**, which non-HT never uses.

The payload is then Viterbi-decoded and descrambled as before. Because almost all 802.11n data is
aggregated, it is split on its A-MPDU delimiters — `reserved(4) | length(12) | CRC-8 | 0x4E` — and each
MPDU is checked on **its own FCS**; a corrupted delimiter resynchronises on the next signature byte
instead of abandoning the rest of the frame.

The Viterbi forward recursion was rewritten for this: the trellis has exactly two predecessors per
state and both arrive on the same input bit, so each step is two 64-wide adds and a compare instead of
a 128-element sort. Identical output (tested against the textbook form), 3.3× faster — which matters
because an A-MPDU is thousands of bits, not the 24 of an L-SIG.

## How — LDPC and STBC (`ldpc.py`, `payload.decode_ht_payload`)
1. **Channel per space-time stream.** STBC sends one stream as two space-time streams from two antennas, the
   second with a −400 ns cyclic shift. The HT-LTFs carry them with P = [[1, −1], [1, 1]], which separates
   h₁ and h₂ per subcarrier. A single HT-LTF per stream is a noisy estimate, so it is smoothed over 1–3
   neighbouring subcarriers after removing its linear phase (the delay and cyclic shift), unless the HT-SIG
   smoothing bit says not to; the width is the one whose channel best predicts the pilots.
2. **Pilots rotate.** Pilot j of data symbol n is Ψ[(n + j) mod 4]·p_{n+3} (802.11-2020 Eq. 19-58), with
   Ψ = (1, 1, 1, −1) for one stream and (1, 1, −1, −1) / (1, −1, −1, 1) for two. The common phase **and the
   phase slope** (sampling-clock drift: 10 ppm on the Microhard link, 0.2 samples over a 2 ms PPDU) are
   measured on the pilots of every symbol and the slope is smoothed over 9 symbols.
3. **Alamouti combining:** symbols 2m, 2m+1 carry (d₀, d₁) on stream 1 and (−d₁*, d₀*) on stream 2, so
   d₀ = (h₁* y₂ₘ + h₂ y₂ₘ₊₁*)/(|h₁|² + |h₂|²) and d₁ = (h₁* y₂ₘ₊₁ − h₂ y₂ₘ*)/(|h₁|² + |h₂|²).
4. **Soft bits:** exact log-MAP per axis from the equalised point and its noise variance (pilot-residual
   level × the L-LTF's per-subcarrier noise shape, smoothed over 9 and clipped to [0.5, 4]).
5. **LDPC** (no interleaver; bits map to subcarriers in order): codeword count and length, shortening,
   puncturing and repetition from the PSDU length (19.3.11.7.5); shortened bits are known zeros, punctured
   ones erasures, repeated ones combined; normalised min-sum (α = 0.8, up to 60 iterations), batched over
   all codewords of the PPDU. Each smoothing width is a candidate receiver and a codeword is taken from the
   first that satisfies all its parity checks.
6. **Data-aided channel re-estimation:** the codewords that passed are known exactly, so their symbols (and
   every pilot) re-train h₁, h₂ by least squares over the whole PPDU; the remaining codewords are decoded
   again, up to 8 rounds. On the real Microhard frames this takes MCS7 codewords from 2/30 to 20/30.
7. Descramble (seed from the 7 first SERVICE bits; the other 9 must descramble to zero, `service_ok`), split
   the A-MPDU, check every FCS. BCC + STBC goes through the same combining, then deinterleave + Viterbi.

The LDPC base matrices are those of the standard as published in Saurabh Tavildar's MIT-licensed LDPC library
(github.com/tavildar/LDPC); their lifting convention, the codeword layout and the repetition rule were confirmed
on real frames: repeated bits match the start of the codeword (83 % agreement at 9 % raw BER, against 49 % for
the parity), and Microhard MCS7 STBC PPDUs decode to CRC-valid MPDUs.

## How — 802.11 at a non-standard clock (`ofdm_scaled.py`)
1. **Bands:** a ~60 kHz spectrogram; per frame, runs of bins ≥ 10 dB over a per-bin floor (20th percentile
   over time, morphologically opened). Runs of 1.2–45 MHz whose edges recur in ≥ 150 µs of frames are
   candidates; adjacent fragments active in the same frames (a burst split by a fade) are also offered as one.
2. **Numerology:** the band is shifted, filtered, decimated and gated to its active frames (gaps of up to 3
   frames closed: a flickering gate imprints the frame rate). The autocorrelation peak over lags of 2.4–26 µs
   gives T_u (CP correlation, ≥ 0.06 and ≥ 5× the median). The CP-correlation metric over a quarter symbol
   is periodic at the symbol rate: its strongest line (frame-rate harmonics excluded) gives T_s. 802.11 needs
   T_s/T_u − 1 = 0.25 (0.2–0.3 accepted) and 64 subcarriers: clock = 80/T_s. The occupied width at −10 dB from
   the in-band level must be 0.7–1.0 × the clock.
3. **Decode:** for each clock candidate (the 0.5 MHz grid values within 1.5 % of the measurement, then the
   measurement) and each orientation (the capture and its complex conjugate), `WifiDetector` runs with the
   capture's sample rate declared as fs·20 MHz/clock, three scan centres 1 MHz (scaled) apart, slightly more
   sensitive thresholds (STF 0.6, 2.5 dB, 40-sample plateau). Every PPDU still needs a parity-checked,
   bit-agreeing L-SIG or a CRC-valid HT-SIG; the clock and orientation with the most such PPDUs win.
4. Packets are converted back to real units (centre, CFO, times, rates × clock/20, span) and the payload
   decoders run on their time-scaled state unchanged.
5. A standard 20 MHz channel on the grid is left to the channel-plan scan; two band candidates of one link
   are merged.

## How — transmitter attribution (`links.py`)
Applied to the PPDUs of every OFDM signal (standard and scaled):
1. **mac header:** an FCS-valid MPDU's TA.
2. **response:** the next PPDU starts within 84 µs × (SIFS/16 µs) of the end of this one (vendors stretch the
   response delay: the Microhard link answers after 46–61 µs at an 11 MHz clock) and is an ACK / Block Ack /
   CTS addressed to X → this PPDU is X's. An ACK / CTS after a frame addressed to Y was sent by Y.
3. **scrambler counter:** a station is "sequential" when ≥ 80 % of the steps between its known seeds fit the
   PPDUs in between. Then a gap that equals exactly the number of unattributed PPDUs in it gives them all to
   the station, and a single PPDU whose own (`service_ok`) seed fits the counter on both sides (allowing two
   transmissions the capture missed; probability of a random fit < (5/127)²) or the next 1–2 steps at its end
   is the station's. Iterated until nothing changes. A station with random seeds is never used this way.
4. Per transmitter: vendor, counts by evidence, formats, frame types, median received power (L-STF/L-LTF) and
   CFO, counter behaviour and seeds. A PPDU still unattributed gets `ta_likely` when exactly one transmitter is
   within 3 dB in power and every other is more than 6 dB away (power says where, not who).

## Drone links inside 802.11 frames (`mac80211.py`; research, PLAN G1 / G2 / G7, synthetic frames only)
- **wfb-ng / OpenIPC / OpenHD:** broadcast data frames (FC 08 01) whose TA = BSSID = 57:42 + channel_id (big-endian,
  link_id << 8 | radio port) → `wfb_ng` {link_id, default link domain (0x7505D6, "default"), radio port, stream (video /
  MAVLink / tunnel / bind, reverse direction), packet type data / session, FEC block / fragment}; flagged as a drone
  link. The payload is ChaCha20-Poly1305 encrypted and is not decrypted.
- **ESP-NOW / ExpressLRS Backpack:** vendor action frames (category 127, OUI 18:fe:34) → `esp_now` {version, body};
  an MSP v2 body with a valid CRC-8 (DVB-S2) sent between identical UID-derived MACs (or a bind to broadcast) →
  `elrs_backpack` {bind UID, command, VTX band / channel from the default Betaflight VTX table}. ESP32 pairs may use
  Espressif Long Range, which is not 802.11 and is not decoded.
- **MAVLink over open Wi-Fi:** any CRC_EXTRA-verified MAVLink frame in a UDP / TCP payload (ports 14550 / 14555 /
  14540 / 5760 named) → `ipv4.mavlink` (the SiK summary); a vehicle sysid flags the frame as a drone.

## HE-SIG-A (802.11ax; PLAN W2, `ofdm_rx.parse_he_sig_a`)
An HE PPDU repeats L-SIG (RL-SIG) and sets L-SIG LENGTH mod 3 (2 = MU, 1 = SU / ER SU / TB). HE-SIG-A follows: two
BPSK r1/2 symbols on the 52 data subcarriers -28..28 (the HT set, interleaved with N_COL 13), pilot polarity p_2 / p_3.
The outer tones +-27 / +-28, which the L-LTF never trains, carry -1 -1 -1 +1 in L-SIG and RL-SIG, so their channel is
estimated there. 52 bits after Viterbi: A1 B0-B25, A2 B0-B25, checked with the 4 MSBs of the HT-SIG CRC-8 over the
first 42 bits and a zero tail (802.11ax-2021 27.3.11.7). Fields: format (SU when A1 B0 = 1, else TB; MU from the
LENGTH), UL/DL, BSS colour, MCS / DCM / streams / LDPC / STBC (SU), SIG-B MCS (MU), bandwidth (punctured MU modes
included), GI + LTF, TXOP (7 bits: 127 none, B0 the 8 / 128 us granularity). The reference preamble is
`ofdm_tx.he_sig_a_bits` + `he_preamble`.
- Real: all 44 HE PPDUs of `testdata/live/LIVE_WIFI_5755MHz` pass CRC-4 (downlink: colour 14, MCS 6, 2 streams,
  80 MHz, TXOP 56 us; uplinks colour 0, MCS 5-11).
- Synthetic (`tests/receivers/test_det_wifi_he.py`): SU / MU / TB round trips; 100 % from 8 dB in-band SNR (6 dB 15 %).
- Fields: `he_airtime_ms_by_direction`, `he_formats`, `he_bss_colors`, `he_bss_from_color` (a colour announced by exactly
  one BSSID's HE Operation element in the capture). The payload is not decoded.

## Drone Wi-Fi links by behaviour (`field.py`; research, PLAN W1, research file 18)
Per signal, from every PPDU (decoded or not, with its transmitter from `links.attribute`): one profile per BSS
(`wifi_links`): AP, clients (the transmitter of to-DS frames, the receiver of from-DS frames, either end of
association / authentication; never SA / DA, which can be wired hosts), downlink / uplink PPDUs, airtime, data
airtime, bytes (non-HT L-SIG length, HT-SIG length, 802.11b PSDU), frames >= 600 B, median size, downlink gaps,
video payloads. A `drone_link` verdict comes from identity (any frame of the BSS flagged `drone`: Remote ID, Tello
SDK, MAVLink, wfb-ng, ESP-NOW, SSID / OUI signature) or else from behaviour: one client, downlink data airtime >= 3x
the uplink with >= 5 data PPDUs, >= 5 % airtime or >= 1 Mb/s, uplink median < 200 B -> "drone / FPV video link-like"
(0.55; 0.65 with a video payload), with the caveat that an IP camera or a phone hotspot streaming video looks the same
in one block. The signal gets `drone_links` and the label suffix ", drone link-like (...)"; `receivers/devices` turns
behaviour-only links into a `drone_link_behaviour` finding. Tests: `tests/receivers/test_wifi_field.py` (encrypted
drone-style link vs a 3-client home AP, Tello identity, the OFDM plugin end to end on a synthetic capture).

**Open links (PLAN W3 / W4, `mac80211`):** Tello SDK command / state datagrams (UDP 8889 / 8890; the `wifi` / `ap`
commands' credentials are not reported), video datagrams (H.264 / H.265 Annex B, RTP dynamic payload with H.264 NAL
units), and vendor-documented default hotspot names as platform signatures (model level for Wi-Fi drones: Spark / Spark-RC, MavicAir-, Mavic-Mini-, PHANTOM3_, ANAFI-, CGO3_ / CGO3P_ (Yuneec), RUKO-F11; research file 18): RMTT- (RoboMaster TT), HolyStoneFPV-,
WiFi-720P- (Eachine E58 class), ExpressLRS TX / RX (Wi-Fi setup mode), ArduPilot / PixRacer (MAVESP8266), CUAVWLINK,
DroneBridge.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.wifi import WifiDetector, decode_nonht_payload
from spectrum_intel.receivers.wifi import dsss, mac80211, dsss_tx, droneid, ofdm_tx, payload
from spectrum_intel.core.dsp import IQContext
packets = WifiDetector(fs, fc).run(x)
ppdus = dsss.scan(IQContext(x, fs, fc))        # dicts: rate, preamble, channel, snr_db, cfo_hz, fcs_ok, mac, ...
info = mac80211.parse(mpdu_without_fcs)        # info["summary"]; info["remote_id"], info["drone"] if present
x = dsss_tx.scene([dict(mpdu=dsss_tx.beacon(ssid="lab"), f_hz=2437e6, t_s=1e-3, snr_db=5)], n=200_000)

# 802.11n payload, from a detected HT packet p
r = payload.decode_ht_payload(p.extra["_rx"][0], p.mcs, p.extra["ht_length"],
                              short_gi=bool(p.extra["short_gi"]), aggregation=bool(p.extra["aggregation"]))
r["subframes"][0]["mac"]["summary"]            # each A-MPDU subframe, with its own fcs_ok

# LDPC / STBC (pass the HT-SIG flags through)
r = payload.decode_ht_payload(rx, p.mcs, e["ht_length"], stbc=e["stbc"], ldpc=e["ldpc"], smoothing=e["smoothing"])
r["scrambler_seed"], r["service_ok"], r["n_codewords_ok"]

# 802.11 at another clock / frequency / orientation, anywhere in the capture
from spectrum_intel.receivers.wifi import ofdm_scaled, links
for ln in ofdm_scaled.scan(x, fs, fc):          # Link(center_hz, clock_hz, inverted, numerology, packets)
    print(ln.center_hz, ln.clock_hz, ln.inverted, len(ln.packets))
rows = links.attribute(records, sifs_us=16.0)   # records as built by plugin_ofdm._records; sets r["ta"]

# reference OFDM transmitter (validation)
s = ofdm_tx.ppdu(ofdm_tx.ampdu([mpdu_a, mpdu_b]), mcs=5, aggregation=True)
s2 = ofdm_tx.ppdu(mpdu_a, mcs=7, ldpc=True, stbc=True)        # two antennas, standard cyclic shifts
y = ofdm_tx.timescale(s2, 11e6, 61.44e6)                      # the same PPDU at an 11 MHz clock
x = ofdm_tx.scene([dict(samples=s, f_hz=2437e6, t_s=1e-3, snr_db=20)], n=200_000)
```

## Pipeline fields (802.11b, one signal per channel)
- **Counts:** `n_ppdus` (CRC-valid headers), `n_fcs_ok`, `n_fcs_corrected`, `rates_mbps`, `preambles`,
  `snr_db_inband_median`, `cfo_hz_median`.
- **Legacy keys:** `bssids`, `ssids`, `frame_types`.
- **`bss`:** one row per BSSID: SSID, channel, vendor, security, PHY, uptime, beacon interval, country,
  BSS load, HT operation, WPS identity, the stations seen talking to it, and the frame types.
- **`stations`:** non-AP transmitters: vendor, randomised or not, associated BSSID, probed SSIDs, device name.
- **`frames`:** one line per PPDU: time, rate, preamble, SNR, CFO, FCS, corrections, `summary`.
- **`ppdus`:** the full records, including the dissected `mac` and `mpdu_hex`.

The signal's own `snr_db` comes from the energy detections it absorbs, on the same scale as every other signal.
DSSS decodes below 0 dB in-band SNR because of its 10.4 dB spreading gain.

## Pipeline fields (OFDM: `wifi_ofdm` per channel, `wifi_ofdm_scaled` per link)
- **Counts:** `n_ppdus`, `formats`, `channels`, `bandwidths_mhz`, `n_fcs_ok`, `n_mpdus`, `bssids`, `ssids`,
  `n_ht_ppdus` / `n_ht_decoded`, `payload_not_decodable` (reason → count).
- **`transmitters`:** one row per TA — `vendor`, `n_ppdus`, `evidence` (count per kind), `formats`,
  `frame_types`, `power_db_median` / `power_db_spread`, `cfo_khz_median`, `scrambler_counter`
  (sequential / not sequential), `seeds`, first/last time, `n_likely_by_power`; and `n_attributed`.
- **`ppdus`:** per PPDU, besides the PHY fields: `ta`, `ta_evidence`, `ta_likely`, `scrambler_seed`, `seed_ok`,
  `power_db`, and for HT `psdu_bytes`, `coding`, `stbc`, `ldpc_codewords_ok`, `n_subframes` / `n_fcs_ok`; for HE
  `he_format`, `bss_color`, `txop_us`, `direction`, `nss` / `ldpc` / `dcm` (SU).
- **HE:** `he_airtime_ms_by_direction`, `he_formats`, `he_bss_colors`, `he_bss_from_color`.
- **Drone links (both PHYs, also 802.11b):** `wifi_links` (top 10 BSS profiles) and `drone_links` (bssid, ssid,
  clients + the verdict: basis identity / behaviour, platform, confidence, behaviour checks, caveat).
- **Scaled links also:** `center_mhz`, `clock_mhz`, `channel_width_mhz`, `subcarrier_spacing_khz`,
  `spectrum_inverted`, `numerology` (T_u, T_s, guard, CP correlation, measured clock, occupied width) and
  `standard_equivalent`. Rates are real (L-SIG rate × clock/20); the label names the vendor of the FCS-verified
  transmitters, e.g. "Wi-Fi 802.11n-type link, 11 MHz channel at 2069.02 MHz, spectrum-inverted (Microhard
  Systems Inc.)".

## Validation
- **`tests/receivers/test_det_wifi_dsss.py`** (pytest, about 20 s): uses the synthetic transmitter (`dsss_tx.py`:
  raised-cosine chips, CFO, chip-clock ppm, multipath, AWGN). It checks:
  - all 4 rates × both preambles, byte-exact;
  - sensitivity;
  - CFO ±120 kHz, −25 ppm clock on 1500-byte CCK, and multipath to 350 ns;
  - every FCS correction is byte-exact;
  - no headers on noise;
  - two BSSs at once;
  - a BLE-like collision rescued by excision;
  - the MAC dissector;
  - Phantom 4 capture regression: ≥ 13 FCS-valid frames, with the known SSIDs.
- **Sensitivity:** FCS-valid, byte-exact fraction, by in-band SNR (22 MHz), for 230-byte MPDUs at 50 kHz CFO:

  | rate | −4 dB | −2 dB | 0 dB | 2 dB | 4 dB | 6 dB |
  |---|---|---|---|---|---|---|
  | 1 Mb/s long | 0.8 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
  | 2 Mb/s short | 0.0 | 0.9 | 1.0 | 1.0 | 1.0 | 1.0 |
  | 5.5 Mb/s short | 0.0 | 0.0 | 1.0 | 1.0 | 1.0 | 1.0 |
  | 11 Mb/s short | 0.0 | 0.0 | 0.0 | 0.7 | 1.0 | 1.0 |

- **Real captures**, against the previous receiver on the same IQ:
  - fresh 2.4 GHz snapshots of the 2026-10-04 sweep (141128_910, 141130_208, 141131_506, plus 141130_534, which
    holds 2437.2 MHz IQ):

    | | PLCP headers | FCS-valid frames |
    |---|---|---|
    | previous receiver | 41 | 23 |
    | this receiver | 74 | 32 |

    No frame the old receiver decoded is lost.
  - Phantom 4 capture (really 2430 MHz): 10 → 16 FCS-valid frames.
- **Remaining FCS failures** on real captures are mostly collisions. Mid-frame, a co-channel OFDM frame lands
  7–15 dB above the beacon, which is beyond the DSSS spreading gain. They are recorded with `fcs_ok: false`.
- **802.11n (HT)** — `tests/receivers/test_det_wifi_ht.py`, against the reference transmitter
  `ofdm_tx.py` (which takes its constellation and code straight from the receiver's own tables, so a
  failure is a real failure and not the two sides disagreeing about a convention):
  - all 8 single-stream rates × both guard intervals, byte-exact;
  - A-MPDU de-aggregation byte-exact over 4 subframes, every delimiter CRC valid;
  - MCS 0–7 all MPDUs of a 1952-byte, 8-subframe A-MPDU byte-exact by in-band SNR (20 MHz):

    | MCS | modulation | Mb/s | needs |
    |---|---|---|---|
    | 0 | BPSK 1/2 | 6.5 | 10 dB |
    | 3 | 16-QAM 1/2 | 26.0 | 15 dB |
    | 5 | 64-QAM 2/3 | 52.0 | 25 dB |
    | 7 | 64-QAM 5/6 | 65.0 | 25 dB |

    A single 61-byte MCS1 frame (the shape of a real 802.11n beacon) decodes from 8 dB.
  - MIMO / LDPC / STBC / 40 MHz are refused with a reason; an oversized payload is refused on budget;
  - the rewritten Viterbi is tested to decode identically to the textbook formulation.
  - **Not yet confirmed on a real 802.11n frame.** 116 Wi-Fi-band snapshots of the 2026-10-04 sweep hold
    only 13 HT PPDUs with a valid HT-SIG CRC: 2 are MCS14 (two streams, correctly declined) and 11 are
    repeats of one 61-byte MCS1 beacon at 5235 MHz — inside a documented ADC-clipping band, with the ADC
    railing on 0.04 % of samples and L-LTF coherence 0.84–0.92. None passed FCS, nor under 2880 alternative
    timing/pilot/interleaver conventions — none of which was the right one. **HT pilots rotate by one
    position every symbol** (802.11-2020 Eq. 19-58); the decoder and the reference transmitter both used a
    fixed [1, 1, 1, −1], so the synthetic tests agreed with themselves. Fixed 2026-10-07: **8 of the 11 now
    pass FCS** (0 with the fixed pattern, same decoder otherwise), vendor action frames broadcast by
    bc:a9:93:de:bf:1b, whose scrambler seeds step 49 → 62. They are kept as `testdata/wifi/real_ht_mcs1_5240`
    and `test_det_wifi_ht.py::test_real_802_11n_frames` pins them.
- **802.11 at non-standard clocks, LDPC, STBC, attribution** — `tests/receivers/test_det_wifi_scaled.py`
  (~30 s):
  - all 12 LDPC codes; the PPDU codeword layout round-trips for every branch (1–2+ codewords, 648/1296/1944,
    shortening, puncturing, repetition, with/without STBC);
  - HT byte-exact for BCC/LDPC × with/without STBC, MCS 0/3/5/7, and an MCS7 LDPC STBC A-MPDU with every
    codeword passing; a 7 kB MCS7 LDPC STBC A-MPDU with a 40 ppm sampling-clock offset;
  - synthetic links at 5 MHz (12 MHz off centre), 10 MHz mirrored and 11 MHz mirrored with MCS7 LDPC STBC data:
    clock, orientation and centre found, every MPDU FCS-valid, every PPDU attributed to the right one of two
    stations (data station / ACK station 10 dB stronger, both with sequential counters);
  - no link from a standard 20 MHz channel (left to the channel-plan scan), LTE-like OFDM bursts or noise;
  - attribution rules on hand-made records, and none from a station with random seeds;
  - **real Microhard link** (`testdata/D_20260211-161341`, 2069.02 MHz, 11 MHz clock, mirrored): 50 PPDUs,
    29 with CRC-valid MPDUs (26 control frames, 3 MCS7 LDPC STBC A-MPDUs); 49 of 50 PPDUs attributed —
    00:0f:92:fd:58:12 sends all 24 HT PPDUs and the RTS (19 by response, 5 by MAC header, 1 2-stream MCS15 by its
    counter gap), 00:0f:92:fd:3e:bf the 24 ACK/BA/CTS (counter 69 → 92 unbroken), 16 dB apart in power; the
    MCS14 PPDU is left with `ta_likely` 58:12;
  - **real DJI link** (`testdata/D_20251119-183758`, 2432.03 MHz, 5 MHz clock): 13 PPDUs, 8 FCS-valid (WEP data
    34:d2:62:53:20:37 → 48:1c:b9:15:c4:15, ad hoc, and one ACK), all 13 attributed to 34:d2:62:53:20:37 (counter
    91 → 104, one transmission hidden by a clipping event).
  - Sweep sample (80 snapshots, 300–6000 MHz, not a test): one off-grid link, a MikroTik station at
    5101.4 MHz (63 FCS-valid frames); nothing else; ~0.8 s per snapshot.
- **OFDM:**
  - **16-QAM / 64-QAM bit-mapping bug, fixed 2026-10-05.** The soft demapper interleaved the I and Q
    bits of each subcarrier (`b0I b0Q b1I b1Q`) instead of splitting them (`b0I b1I | b0Q b1Q`, 802.11-2020
    Figures 17-10/17-11). BPSK and QPSK carry one bit per axis, so 6/9/12/18 Mb/s were unaffected and the
    L-SIG always decoded perfectly — but **every 24/36/48/54 Mb/s frame failed FCS**, as would HT MCS 3–7.
    It survived because the reference transmitter takes its constellation from the receiver; it was found on
    live 5 GHz captures, where 24 Mb/s decoded 0 of 31 frames at every gain. A test now pins the mapping to
    the standard's figures written out by hand, independently of `ofdm_tx` (it fails on the old code).
    Effect on the same IQ:

    | capture | FCS-valid non-HT frames, before → after |
    |---|---|
    | live 5 GHz, 2026-10-05, UNII-1 + UNII-3 (8 blocks) | 39 → 70 |
    | live 5200 MHz at gain 40, 24 Mb/s only | 0/13 → 11/13 |
    | D_20260916 (5770 MHz, the regression reference) | 8 → 13 |
    | sweep 5235.36 MHz | 3 → 17 |
    | 2.4 GHz snapshots (beacons, 1–6 Mb/s) | unchanged |

    The 24 Mb/s frames still failing on D_20260916 are weak (median 9.6 dB against 11.3 dB for those that
    pass), at the 16-QAM threshold.
  - **Checked against gr-ieee802-11** (radioconda's base Python; `.venv` has NumPy 2, which it cannot load;
    wiring copied from its `examples/ieee802_11/wifi_phy_hier.grc`), on every 20 MHz channel fully inside
    the capture, comparing unique FCS-valid MPDUs byte for byte:

    | | gr-ieee802-11 | ours | both | only gr-ieee802-11 | only ours |
    |---|---|---|---|---|---|
    | live 5 GHz, 2026-10-05 | 52 | 59 | 52 | **0** | 7 |
    | D_20260916 + sweep 5235 MHz | 14 | 23 | 14 | **0** | 9 |

    Every frame the reference decodes, this receiver decodes too; the extra 16 are beacons, RTS/CTS, null and
    Block Ack frames. (Before the fix this could not have held: the reference decodes 16-QAM.)
  - `usable_bw` was raised from 56 to 60 MHz so that weak edge subcarriers are down-weighted by their L-LTF SNR
    instead of erased (+1 FCS-valid long frame on a channel at the capture edge, no loss elsewhere);
  - **VHT-SIG-A fields** are checked on live data: the two reserved bits (B2 and B23, fixed at 1 by the
    standard) land where expected in all 44 valid words, and Group ID reads only 0 or 63.
  - D_20260916 (5770 MHz) is the OFDM regression reference.
- **`tests/test_core.py`:** the DSSS descrambler inverts the scrambler.

## Limitations
- 802.11b PPDUs whose preamble lies before the capture start cannot be decoded (they are reported as
  DSSS-like). Neither can collided PPDUs whose interferer covers the whole band.
- A few weak Barker-11 1 Msym/s bursts per capture have no 802.11b SFD even with 3 bit errors (non-802.11
  framing). They are reported as "Barker-11 DSSS (no 802.11b PLCP header)-like".
- OFDM: VHT and HE payloads are not decoded, nor HT with 2+ spatial streams or 40 MHz — they need more than
  one antenna or more bandwidth (VHT LDPC at 20 MHz would only need the VHT framing around `ldpc.py`). Each is
  reported with the reason (`payload_not_decodable` in the signal's fields) rather than silently skipped.
- LDPC at the edge: MCS7 needs ~18 dB after equalisation; the Microhard HT PPDUs arrive at ~16.5 dB SINR, so
  most of them decode only partly (codewords that pass are kept, the MPDU needs all of its own).
- Non-standard clocks: the search needs ≥ 150 µs of bursts with sharp edges, the 1/4 (long) guard interval
  present (every 802.11 PPDU's legacy preamble has it) and a 2.5–26 MHz clock. 40 MHz-type links and
  continuous OFDM are not searched; links on the standard 20 MHz grid are the channel-plan scan's.
- Scrambler-counter attribution needs a chipset that steps its seed (Atheros/Qualcomm-style); with random seeds
  it simply does not apply (reported as "not sequential").
- **At 5 GHz most data cannot be captured, whatever the software.** Every VHT PPDU in the live 2026-10-05
  captures is 80 MHz (21 single-stream LDPC, 5 two-stream). 80 MHz does not fit in a 61.44 MSps capture,
  so those payloads are out of reach of an ADALM-Pluto. What is decodable at 5 GHz is the non-HT traffic —
  beacons, probe responses, RTS/CTS/ACK and Block Ack — which carries the MAC identities, plus the VHT
  headers above. Reaching the data needs a receiver with ≥ 80 MHz of bandwidth and LDPC decoding.
- **Partial AID is usually 0 in practice.** On the live captures every uplink VHT PPDU carried Partial
  AID 0 or 1 rather than the BSSID bits the standard specifies, so BSS attribution through it rarely
  succeeds; it is attempted only on a non-zero value with a unique match.
- **Gain at 5 GHz:** at 50 dB a nearby AP rails the Pluto's ADC (0.3–0.7 % of samples) and 16-QAM suffers
  (24 Mb/s: 4/15 decoded at 50 dB against 11/13 at 40 dB, live). 40 dB gave no clipping and the best decode
  rate; this is an acquisition setting in signal-vault, not something the library corrects.
- **Remote ID is self-declared.** It is broadcast in clear and unauthenticated, so the values are what
  the aircraft chose to say. An aircraft that does not broadcast it is invisible to this; the platform
  signatures only identify equipment. Neither is a flight detection.

## Off-grid sensitivity and China Remote ID (2026-10-23, PLAN G13 / G14)
- `offgrid_plan`: a slice that fails the whole-capture CP test but has a whole-capture ratio ≥ 2 is retested per 250 µs
  window (≥ 2 windows with correlation ≥ 0.06 and ≥ 6 × the reference lags). The weak, low-duty Mimosa link of
  `LIVE_MIMOSA_4880MHz_weak` (~6.6 dB, ~20 % duty; the whole-capture test lit 2 of the 3 required slices) is now planned
  and decoded (21 PPDUs, 14 FCS-ok). Cost: ~25 ms per gated slice in Python, nothing on noise.
- The F3411 beacon element (FA:0B:BC / 0x0D) may carry a China GB 46750-2025 packet (first byte 0xFF);
  `core.remote_id.parse_gb46750` decodes it (research; checked against a published DJI Mini 5 Pro frame).
