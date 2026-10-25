# ism_devices — 315 / 433 MHz remotes and TPMS (research)

**Validation: research.** Built from `research/signal_gaps/04_short_range.md` section 1, with the bit conventions
of the rtl_433 decoders read as a fact source only (GPL code, nothing copied). Checked on synthetic devices only.

**Bands:** 314–316 MHz (imported cars, US-market fobs) and 433.05–434.79 MHz (India SRD, GSR 1048(E)).

**Detects and decodes:**

| Device | Label | Fields | Verified by |
|---|---|---|---|
| TPMS, Continental / VDO (Hyundai-Kia, BMW, FCA, PSA, Mazda families) | `TPMS (Continental / VDO, ID …)` | `id`, `state`, `flags`, `repeat`, `pressure_kpa` (raw × 1.375), `temperature_c` (raw − 50), `battery_raw`, `frame_hex` | CRC-8 (poly 0x07, init 0xAA) over bytes 0–8 → `success`, 0.9 |
| TPMS, Toyota / Pacific PMV-C210 | `TPMS (Toyota / Pacific PMV-C210, ID …)` | `id`, `status`, `pressure_psi` (raw × 0.25 − 7), `pressure_kpa`, `temperature_c` (raw − 40) | CRC-8 (0x07, init 0x80) over bytes 0–7 and the inverted pressure copy → `success`, 0.9 |
| Fixed-code remotes: EV1527 / HS1527, PT2262 / SC2262 (gates, shutters, doorbells, alarm PIRs, RF sockets) | `OOK fixed-code remote (EV1527 / PT2262-like, ID …)` | `code` (24 bits), `id` (20 bits), `data` (4 bits), `tristate` (PT2262, when every bit pair is 00 / 11 / 01), `repeats`, `words`, `te_us` | no checksum: ≥ 2 identical words → 0.75; one word → 0.45 (stays `unknown`) |
| KeeLoq HCS200 / 300 / 301 (code-hopping car, garage and gate fobs) | `KeeLoq code-hopping remote-like (serial …)` | `serial` (28 bits), `button`, `learn`, `battery_low`, `repeat`, `packets`, `matching_packets`, `te_us` | structure (12 warm-up pulses, 10 Te header, exactly 66 bits); 2 matching packets → 0.75, one → 0.6 |
| anything else keyed | `OOK device-like` / `FSK device-like (unrecognised protocol)` | pulse counts and mark widths | 0.55 |

The KeeLoq **encrypted 32-bit hopping code is neither decoded nor reported**: only the clear serial, buttons and
status bits are.

**How:**
1. Candidates: the band at ≥ 2 MSps in 256-point frames; per frame the bin furthest over its own median across the
   frames (≥ 12 dB); those frequencies grouped at 100 kHz gaps (an FSK pair stays one device), ≥ 3 frames per
   group, its median frequency; the 4 groups with the most frames.
2. Per carrier: ±100 kHz at 1 MSps; envelope (8-sample average) on above the geometric mean of its 10th and 99.5th
   percentiles (≥ 12 dB apart); runs < 30 µs debounced → pulses (start, mark, space).
3. OOK decoders on the pulse list (`ook.py`): EV1527 / PT2262 (sync 1 Te + 31 Te, 24 PWM bits of 4 Te: 0 = 1 + 3 Te,
   1 = 3 + 1 Te); KeeLoq (12 × (1 Te, 1 Te), 10 Te header, 66 PWM bits of 3 Te, short mark = 1, LSB first).
4. FSK: a mark ≥ 2 ms whose discriminator is two-level (deviation ≥ 10 kHz, ≥ 60 % of samples within 0.6–1.4 of
   it: a LoRa chirp or FM voice is not); sign runs rounded to 52 µs symbols; `tpms.py` on both polarities: VDO
   (inverted stream, preamble aa aa aa a9, Manchester = the second symbol of a pair, 80 bits) and Toyota (preamble
   a9 e, differential Manchester from its last symbol, ≥ 72 bits).
5. `arbitrate`: a device decoded here (TPMS CRC, repeated code, KeeLoq) explains another receiver's unverified guess
   on the same burst (e.g. rc_link's FSK-rate "SiK-like"); it is recorded on the device's signal
   (`explained_receiver_detections`).

**Accuracy (synthetic, 2 MSps, SNR per sample over 2 MHz):** OOK remotes and KeeLoq 6/6 at 6 dB, 0/6 at 3 dB
(the 12 dB envelope gate); TPMS 12/12 down to 0 dB; a VDO frame in a 2^21 block at 61.44 MSps decodes (~10 ms).
No detection on noise, a LoRa SF9 packet or NBFM speech.

**Limitations:**
- No real capture tested yet; the KeeLoq PWM sense and the TPMS layouts come from rtl_433 and must be confirmed on a
  real fob / sensor.
- A 34 ms block (2^21 at 61.44 MSps) holds a TPMS frame but not a whole EV1527 frame (≈ 45 ms at Te 350 µs) or a
  KeeLoq packet (≈ 95 ms): those need longer blocks at a lower rate.
- Only two TPMS families (Schrader, Jansite and others are listed in the research file); the symbol time is fixed at
  52 µs.
- Weather sensors and other rtl_433 protocols are not decoded: they appear as `OOK / FSK device-like`.
