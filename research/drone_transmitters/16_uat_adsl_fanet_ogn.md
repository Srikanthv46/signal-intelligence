# 16. UAT 978, ADS-L SRD-860, FANET, OGN tracker, China GB 46750-2025

Research 2026-10-22/23 for PLAN G10 / G14 (file 14 §4). The facts below are precise enough to implement.

**Tags:**
- **[D]** = from the standard or from a permissively licensed source;
- **[S]** = from a GPL or unlicensed source (read as a spec only, never copied);
- **[calc]** = derived here;
- **[U]** = uncertain or unverified.

**Bit numbering:** each section states its own convention, because the systems differ.

**What our existing blocks give us:**
- `rc_link/flarm.py`: the FLARM PHY demodulator (100 kchip/s Manchester, ±50 kHz, sync search, both polarities). It
  serves ADS-L M-band and OGNTP once the sync word changes.
- `rc_link/lora_phy.py`: LoRa payload bytes plus the sync-word byte (`sync_word()`, `identify()`). It serves FANET.
- `core/remote_id.py` (F3411 over the Wi-Fi beacon vendor IE and BLE 0xFFFA): it serves GB 46750 after one branch
  on the first byte 0xFF.
- `receivers/adsb`: the 1090 decoder. UAT is a separate PHY.

---

## A. UAT 978 MHz (ADS-B on the Universal Access Transceiver; US, and large UAS in the US)

**Primary sources:**
- FlightAware dump978 (BSD-2) [D]: `uat_protocol.h`, `fec.cc`, `demodulator.cc`, `uat_message.cc`.
- ICAO UAT SARPs (ACP-WGW01 Appendix B) [D] for the PHY and RS constants.
- The UAT Implementation Manual v1.5 [D] for frame timing.

**Licence trap:** dump978's Reed–Solomon code `libs/fec/` is Phil Karn's fec-3.0.1, which is **LGPL**. Do not vendor it
into `c++/`. Write our own RS(255,k) decoder (Berlekamp–Massey + Chien + Forney, errors and erasures) from the
constants below; it is about 150 lines. The rest of dump978 is BSD-2 and may be followed directly.

### PHY [D: SARPs §12.4]
| Item | Value |
|---|---|
| Carrier | 978 MHz. Tolerance ±20 ppm aircraft, ±2 ppm ground |
| Modulation | binary CPFSK; 1 = shift **up**, 0 = shift **down** |
| Rate | 1.041667 Mb/s (bit period 0.96 µs) |
| Deviation | nominal ±312.5 kHz at the optimum sampling point (eye ≥ 560 kHz) |
| Modulation index | h ≥ 0.6 (2·312.5 / 1041.667 = 0.6) |
| Filtering | not specified. Tx filtering exists for spectral containment, so treat the signal as filtered CPFSK |
| Occupied BW | about 1.3 MHz [calc]. The Pluto at 61.44 MSps covers it trivially |

### Frame timing [D: UAT Implementation Manual §2.3]
- **UAT frame:** 1 UTC second, holding **3952 Message Start Opportunities (MSOs)** spaced 250 µs.
  - MSO 0 is at 6 ms.
  - The **Ground segment** (uplinks, 32 slots of 5.5 ms) spans 176 ms from MSO 0 to MSO 704.
  - After a 12 ms guard comes the **ADS-B segment**, 800 ms = MSO 752…3951.
- **Aircraft:** each aircraft sends exactly **one** ADS-B message per second, at a pseudorandom MSO.
- **Durations** [calc]:

| Message | Bits | Duration |
|---|---|---|
| Basic | 36 + 240 = 276 | 265 µs |
| Long | 36 + 384 = 420 | 403 µs |
| Uplink | 36 + 4416 = 4452 | 4.27 ms |

### Sync words (36 bits, left-most bit first) [D: SARPs + dump978]
| Link | Binary | Hex (dump978) |
|---|---|---|
| ADS-B (downlink) | `111010101100110111011010010011100010` | `0xEACDDA4E2` |
| Ground uplink | `000101010011001000100101101100011101` | `0x153225B1D` (bitwise complement of the downlink) |

### Framing and FEC [D: SARPs §12.4.4; dump978 `uat_protocol.h` / `fec.cc`]

**Frame sizes:**

| Frame | Data | RS code | Air bytes after sync |
|---|---|---|---|
| Basic ADS-B | 144 bits (18 B) | RS(30,18): 12 parity, corrects 6 | 30 |
| Long ADS-B | 272 bits (34 B) | RS(48,34): 14 parity, corrects 7 | 48 |
| Ground uplink | 3456 bits = 6 × 72 B | 6 × RS(92,72): 20 parity each, corrects 10 | 552, byte-interleaved |

**RS field and code:**
- **Field:** GF(256), primitive polynomial **p(x) = x⁸+x⁷+x²+x+1 (0x187)**.
- **Generator:** g(x) = ∏(x − αⁱ) for i = 120 … 120+nroots−1, so the **first consecutive root is 120** with step 1.
  Karn parameters: `symsize 8, gfpoly 0x187, fcr 120, prim 1, nroots 12/14/20`.
- **Shortening:** the codes are shortened. Pad = 255 − n, so 225 (basic), 207 (long) and 163 (uplink) leading zeros.
- **Order:** parity bytes follow the data, most significant coefficient first, and bits within a byte go MSB first.

**Basic or long?** The receiver does not know the length in advance. dump978 demodulates 48 bytes, then:
1. Try RS(48,34). It is accepted only if `payload_type = byte0 >> 3` is **≠ 0**.
2. Otherwise try RS(30,18) on the first 30 bytes. It is accepted only if `payload_type == 0`.

**Uplink de-interleave:** `block[b][i] = raw[i*6 + b]`, for b = 0..5 and i = 0..91. The transmission order is
1, 73, 145, 217, 289, 361, 2, 74, … (column by column of the 6×92 matrix). After correction, concatenate the six 72-byte
data parts (432 B). The payload is FIS-B (DO-358B) and is not needed for drone ID.

**Erasures:** dump978 marks a byte as an erasure when any bit's |Δφ| falls between the 0 and 1 slicer levels. This is
used only with AUTO_CENTER; by default both thresholds are 0, so there are no erasures.

**No whitening or scrambling:** the RS-corrected bytes are the payload.

### dump978 demodulator (BSD-2; may be followed directly) [D: `demodulator.cc`, `convert.cc`]
1. **Input:** **2.083333 MSps = exactly 2 samples per bit**. Each sample is converted to a 16-bit phase:
   `atan2(Q, I)` scaled so that 65536 = 2π. The cu8/cs8 inputs use a lookup table.
2. **Phase difference:** `PhaseDifference(a, b)` = wrapped int16 `(b − a)`. A positive value means a positive
   frequency, i.e. bit 1.
3. **Two sync shift registers:** dump978 runs two 36-bit registers, one per sample phase:
   - `sync0` takes the bit `Δφ(s[2k], s[2k+1]) > 0`;
   - `sync1` takes the bit `Δφ(s[2k+1], s[2k+2]) > 0`.
4. **Sync match:** a match is **Hamming distance ≤ 4** to the downlink or uplink word (`SyncWordMatch`). On a match,
   the decoder demodulates from `start` and from `start+1` (`DemodBest`) and keeps the candidate with fewer RS
   corrections.
5. **Bit slicing:** bit = `Δφ(s[2i], s[2i+1]) > 0`, MSB first into bytes, starting right after the 72 sync samples.
6. **Optional `AUTO_CENTER`:** the slicer threshold is the mean of the sync's 1-bit and 0-bit Δφ, which absorbs the
   CFO. The sync is then rechecked with ≤ 4 errors. We should use this, because the Pluto TCXO is off by about 10 ppm
   (≈ 10 kHz at 978 MHz, ≈ 3 % of the deviation, so harmless anyway).
7. **Tail kept between buffers:** (36 + 4416) × 2 samples.
8. **Porting:** take the 978 MHz channel from the block FFT / Channelizer at a convenient rate, e.g. 4.1667 MSps
   (4 sps). The same Δφ-sign slicer works with the best-phase pick over 4 phases. From 61.44 MSps the rational
   resample to 2.0833 MSps is 625/18432 [calc]; a decimating channelizer to about 4.17 MSps plus interpolated timing
   is cheaper.

### ADS-B payload (DO-282B §2.2.4.5) [D: dump978 `uat_message.cc`]
**Bit addressing** is DO-282 style: `(byte, bit)`, where byte 1 is the first and bit 1 is the **MSB**.
`Bits(b1, i1, b2, i2)` is MSB-first inclusive.

**Header (HDR, bytes 1–4):**

| Field | Bits | Meaning |
|---|---|---|
| Payload type | (1,1)–(1,5), 5 bits | which elements follow; see the table below |
| Address qualifier | (1,6)–(1,8), 3 bits | 0 ADS-B ICAO · 1 ADS-B other (self-assigned) · 2 TIS-B ICAO · 3 TIS-B track file · 4 surface vehicle · 5 fixed beacon · 6 ADS-R other · 7 reserved |
| Address | (2,1)–(4,8), 24 bits | ICAO 24-bit or self-assigned |

**Payload type → elements** (DO-282B Table 2-10, as dump978 decodes it):

| Type | Contents |
|---|---|
| 0 | HDR + SV (basic, 18 B) |
| 1 | SV + MS + AUXSV |
| 2 | SV + AUXSV |
| 3 | SV + MS + TS(@30) |
| 4 | SV + TS(@30) |
| 5 | SV + AUXSV |
| 6 | SV + TS(@25) + AUXSV |
| 7–10 | SV |
| 11–31 | HDR only |

Long messages (types 1–10) are 34 B. **The emitter category and the callsign appear only in MS, i.e. in types 1 and
3.** An aircraft rotates through payload types, so the category arrives every few seconds [U: exact rotation is in
DO-282B §3].

**State vector (SV), bytes 5–17:**

| Field | Bits | Decode |
|---|---|---|
| Latitude | (5,1)–(7,7), **23 bits** | `lat = raw·360/2²⁴`; if lat > 90 then lat −= 180 |
| Longitude | (7,8)–(10,7), 24 bits | `lon = raw·360/2²⁴`; if lon > 180 then lon −= 360 |
| Altitude type | (10,8) | 1 = geometric, 0 = pressure |
| Altitude | (11,1)–(12,4), 12 bits | 0 = unavailable, else `(raw − 41)·25` ft (−1000 … +101 350 ft) |
| NIC | (12,5)–(12,8), 4 bits | navigation integrity category |
| A/G state | (13,1)–(13,2) | 0 airborne subsonic · 1 airborne supersonic (×4 on velocities) · 2 on ground · 3 reserved |

**SV velocity, airborne:**
- **North:** sign (13,4), 1 = south; magnitude (13,5)–(14,6), 10 bits; value = raw − 1 kt (0 = unavailable).
- **East:** sign (14,7), 1 = west; magnitude (14,8)–(16,1), 10 bits; value = raw − 1 kt.
- **Vertical:** source (16,2) (0 = geometric, 1 = baro); sign (16,3), 1 = down; magnitude (16,4)–(17,4), 9 bits;
  value = (raw − 1)·64 ft/min.

**SV velocity, on ground:**
- **Ground speed:** (13,5)–(14,6), value = raw − 1 kt.
- **Track/heading type:** (14,7)–(14,8): 0 none · 1 true track · 2 magnetic heading · 3 true heading.
- **Angle:** (15,1)–(16,1), 9 bits × 360/512.
- **Aircraft size:** (16,2)–(16,5), length/width code per DO-282B Table 2-35.
- **GPS antenna offset:** fields from (16,7).

**Trailer, by address qualifier:**
- ADS-B and vehicle qualifiers: UTC-coupled (17,5) and uplink feedback (17,6)–(17,8).
- TIS-B and ADS-R qualifiers: site ID (17,5)–(17,8).

**Mode status (MS), bytes 18–29:**
- **First word:** bytes 18–19 form a 16-bit word `w1 = emitter_category·1600 + c1·40 + c2`. So
  **emitter category = (w1 / 1600) % 40**.
- **Callsign words:** bytes 20–21 and 22–23 hold `c3..c5` and `c6..c8`, as `w = a·1600 + b·40 + c`.
- **Base-40 alphabet:** `0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ *??` (index 36 = space, 37 = `*`). Trim trailing
  space or `*`.
- **CSID** (27,7): 1 = callsign, 0 = the 8 characters are a flight-plan ID (4-digit octal squawk).
- **Status fields:**
  - emergency (24,1–3); MOPS version (24,4–6); SIL (24,7–8);
  - transmit MSO (25,1–6); SDA (25,7–8);
  - NACp (26,1–4); NACv (26,5–7); NICbaro (26,8);
  - capability: UAT-in, 1090ES-in, TCAS (27,1–3);
  - operational modes: TCAS RA, IDENT, ATC services (27,4–6);
  - SIL supplement (27,8); GVA (28,1–2); single antenna (28,3); NIC supplement (28,4).

**Other elements:**
- **AUXSV:** bytes 30–31 hold the *other* altitude, `(30,1)–(31,4)`, decoded like the SV altitude.
- **TS** (target state), from byte 25 or 30: selected altitude, baro setting, selected heading and mode bits, as in
  `DecodeTS`.

**Emitter category values** [D: DO-282B Table 2-40 as also published in the public Garmin GDL90 ICD §3.5.1.10;
dump978 renders them as `'A'+(c>>3)`, `'0'+(c&7)`, the 1090 set/letter form]:

| Code | Category |
|---|---|
| 0 | no information |
| 1 | light |
| 2 | small |
| 3 | large |
| 4 | high-vortex large |
| 5 | heavy |
| 6 | highly manoeuvrable |
| 7 | rotorcraft |
| 9 | glider/sailplane |
| 10 | lighter-than-air |
| 11 | parachutist/skydiver |
| 12 | ultralight/hang/paraglider |
| **14** | **unmanned aerial vehicle (= "B6")** |
| 15 | space/trans-atmospheric |
| 17 | surface emergency vehicle |
| 18 | surface service vehicle |
| 19 | point obstacle |
| 20 | cluster obstacle |
| 21 | line obstacle |
| 8, 13, 16, 22–39 | reserved |

[U: the code numbers 0–21 are confirmed through the GDL90/dump978 mapping, not from the DO-282B text itself, which
is paywalled.]

**Drone relevance:**
- Most small UAS in the US are **prohibited** from ADS-B Out (14 CFR 89 and 91.225). Large or BVLOS UAS sometimes
  carry UAT transceivers (uAvionix ping-class); see the UK CAA 2025 note: UAT on 978 MHz for UAS under DO-282B.
- Category 14 is the UAV flag. A self-assigned address (qualifier 1) is common on such units.
- **India:** UAT is unused.

### Verification for us
- **Validation data:** dump978 ships `sample-data.txt.gz` (decoded hex lines). RS decoding of those hex frames is a
  unit test.
- **Synthetic IQ:** CPFSK h = 0.6, 1.041667 Mb/s, plus RS encode, gives end-to-end tests.
- **A real capture needs a US location** [U: no 978 traffic in India].

---

## B. ADS-L 4 SRD-860 (EASA; Issue 2, 2025-12-01)

**Sources:**
- **Primary [D]:** the EASA PDF https://www.easa.europa.eu/en/downloads/142821/en (ED Decision 2022/024/R,
  Issue 2, 52 pp).
- **Cross-checks [S]:**
  - pjalocha `esp32-ogn-tracker/main/adsl.h`, `rf.cpp`. The repo has no licence file (all rights reserved).
  - SoftRF `ADSL.cpp/.h` (GPL-3).
- Pawel Jalocha (OGN) is a co-author of the spec, and his `adsl.h` agrees with it.

**Conventions (spec §B.3):**
- Multi-byte integers are **little-endian** (LS byte first on air).
- **Each byte is sent MSB first.**
- Bit fields are numbered from the LSB of the little-endian stream: field offset 0 sits in bit 0 of byte 0, and a
  14-bit field A then a 10-bit field B gives byte0 = A[7:0], byte1 = B[1:0]·A[13:8], byte2 = B[9:2].

### PHY (Subpart C) [D]
| Item | M-band | O-band LDR | O-band HDR |
|---|---|---|---|
| Frequency | **868.2 and 868.4 MHz**. Traffic and Status alternate channel on every transmission | 869.525 MHz | 869.525 MHz |
| Channel BW | 200 kHz | 250 kHz | 250 kHz |
| Modulation | 2-GFSK BT 0.5, **±50 kHz** (0 = −, 1 = +) | 2-GFSK BT 1.0, **±12.5 kHz** | GMSK BT 0.5 (±50 kHz implicit) |
| Chip rate | **100 kchip/s, Manchester → 50 kb/s** (1 → `01`, 0 → `10`) | 38.4 kb/s, no Manchester | 200 kb/s |
| Preamble | `01` × 5…12 (default 8), then `1001 1001`, **not** Manchester-coded | `10` × 40 | `10` × 12 |
| Sync | **0x72 0x4B**, Manchester-coded | **0xB4 0x2B** | **0x2D 0xD4** |
| ERP | 14 dBm | 27 dBm | 27 dBm |
| Backoff | 15–250 ms | 10–50 ms | 10–50 ms |

**Notes:**
- **M-band sync as decoded bytes:** the raw preamble tail `1001 1001` equals Manchester-decoded bits `0101`.
  Together with the `01` run, the decoded stream before the sync ends in `…F5`. OGN tracker programs the 8-byte chip
  pattern `55 99 95 A6 9A 65 A9 6A`, which is **Manchester(F5 72 4B 18)** [calc, matches `rf.cpp`: "0xF5724B18 encoded
  in Manchester (fixed packet length 0x18 is included)"].
  - **FLARM:** F5 31 FA B6. **OGNTP:** 0A F3 65 6C. So `flarm.py` needs only a different 32-bit sync constant:
    **0xF5724B18** for 24-byte traffic packets, or 0xF5724B plus a length byte in general.
- **M-band timing:** there are no time slots. CSMA/LBT applies (ETSI EN 300 220-1), with a forced transmission after
  3 s and then 2 s of silence.
- **Rates:** Traffic is ≥ 1 Hz airborne and 0.1 Hz on ground. Status is ≥ once per 20 s.
- **HDR time slots** (relative to the UTC second): Direct 450–1000 ms (Traffic, Status, Remote ID); Uplink 200–450 ms
  (Traffic Uplink, FIS-B).
- **M-band signature:** the optional signature (68 B) is sent **without Manchester**.
- **O-band LDR/HDR:** Issue 2 corrected the O-band LDR modulation. **Do not use the Issue 1 O-band parameters.**
- **LDR sync bit order:** 0xB4 0x2B is the bit-reverse of 0x2D 0xD4 (the PilotAware/P3I heritage). OGN tracker's
  PilotAware sync is `B4 2B 00 00 00 00 18 71` [S] [U: confirm the LDR bit order on a capture].

### Data link (Subpart D) [D]
**M-band and O-band HDR** (Manchester on M-band):

| Field | Size | Content |
|---|---|---|
| Packet Length | 1 B | bytes in Network header + ADS-L Data + parity, **excluding** the length byte itself and any signature. Traffic: **0x18 = 24** |
| Data | Length bytes | the network packet below |
| Signature | 68 B | optional, not counted in the length |

**O-band LDR** (legacy and PilotAware compatible; no FEC, no signature, max payload 15 B):

| Field | Size | Content |
|---|---|---|
| Header | 6 B | fixed `00 00 00 00 18 71` |
| Data | 24 B | network header + ADS-L data + CRC24 |
| CRC8 | 1 B | poly 0x07, **init 0x71**, MSB-first shift |

[U: whether the CRC8 covers only the 24 data bytes; the spec gives only the function. The 0x71 init equals the last
header byte, which hints that the header is treated as a seed.]

### Network layer (Subpart E) [D]
**Byte 0, the network header** (LSB-first fields):

| Bits | Field | Values |
|---|---|---|
| 0–3 | Protocol version | 0 = Issue 1, 1 = Issue 2 (senders use the lowest version that encodes the payload identically) |
| 4 | Secure-signature flag | |
| 5–6 | Key index | 0 = key 0x00000000 (scrambling), 1–2 reserved, **3 = scrambling disabled** |
| 7 | Error-control mode | 0 = CRC24, 1 = RS FEC |

It is followed by the **ADS-L Data** (scrambled), then the parity.

**Scrambling (E.2):** XXTEA with a **zero key** and **6 rounds** (not the standard 6 + 52/n), DELTA = 0x9E3779B9, over
the ADS-L Data block taken as little-endian uint32 words.
- With key 0, `MX = (((z>>5) ^ (y<<2)) + ((y>>3) ^ (z<<4))) ^ ((sum ^ y) + z)`.
- **Receiver:** run standard XXTEA *decrypt* with rounds = 6, starting at `sum = 6·DELTA`.
- **Word count:** for Traffic, n = 5 words (20 B = 5 header + 15 payload), as in OGN `Descramble(): XXTEA_Decrypt_Key0
  (Word, 5, 6)`.
- [U: the spec's parameter table says `num_data_words = 6`, which contradicts the 20-byte Traffic block. Use
  `len(ADS-L Data)/4`. The OGN implementation uses 5.]
- Data sizes are always 4·k bytes (payload = 3, 7, 11, 15 … bytes).

**CRC24 (E.3.1):**
- It is the **Mode-S / 1090 ES CRC**: G(x) = x²⁴+x²³+…+x¹²+x¹⁰+x³+1, i.e. **0xFFF409**.
- Init 0. It is computed over the **network header + scrambled ADS-L Data** (not the length byte).
- It is appended **big-endian** (MS byte first).
- Reference implementation: register `poly = 0xFFFA0480`, feed bytes into the low 8 bits, then three zero bytes,
  then `>>8`.
- **Receiver check:** run over header + data + 3 CRC bytes; the result must be 0.
- **We can reuse our 1090 CRC.** OGN's `adsl.h` also lists 1-bit syndromes for single-bit correction and soft
  (erasure) flips [S: the algorithm only].

**FEC (E.3.2), for long payloads (> 64 B, uplinks):**
- RS over GF(256) with field polynomial **1+x²+x³+x⁴+x⁸ (0x11D)**.
- **32 parity bytes**, primitive element **11**, first consecutive root **121**.
- The data is zero-prefixed to 223 B; parity covers header + scrambled data.
- [U: map "primitive element 11" to Karn `prim = 11`, `fcr = 121`. Verify with a test vector if one appears.]

**Signature:** a 4-byte UNIX time followed by a 64-byte Ed25519 (R, S) signature. It covers header + scrambled data +
timestamp, not the parity. It is sent at most every 10 s.

### Presentation layer (Subpart F) [D]
**ADS-L Data** = 5-byte header + payload.

| Bit offset | Field |
|---|---|
| 0–7 | Payload type: 2 Traffic · 3 Status · 4 Traffic Uplink · 5 FIS-B Uplink · **6 Remote Identification** · 66 OGN diagnostics · ≥128 unicast (reserved) |
| 8–13 | Address mapping table (AMT): 0 random/privacy · 5 ICAO · 6 FLARM · 7 OGN tracker · 8 FANET · 9–63 manufacturer pages 0–54 (the address is then an 8-bit manufacturer prefix + a 16-bit base) |
| 14–37 | Address, 24 bits |
| 38 | reserved |
| 39 | relay/forward |

**In bytes:** byte 0 = type; bytes 1–4 = LE uint32 with AMT = bits 0–5, address = bits 6–29, relay = bit 31.

### Application layer (Subpart G) [D]
**Traffic payload (type 2, 15 B, payload bit offsets):**

| Offset | Bits | Field | Encoding |
|---|---|---|---|
| 0 | 6 | Timestamp | quarter-seconds since the full hour, mod 60 (60–63 unused); time of the PVT fix |
| 6 | 2 | Flight state | 0 undefined · 1 ground · 2 airborne |
| 8 | 5 | **Aircraft category** | see below |
| 13 | 3 | Emergency | 0 undefined · 1 none · 2 general · 3 medical · 4 no comms · 5 unlawful · 6 downed |
| 16 | 24 | Latitude | signed, **1/93206 °**; 0x800000 = no fix |
| 40 | 24 | Longitude | signed, **1/46603 °** |
| 64 | 8 | Ground speed | exponential, unsigned N = 6, × 0.25 m/s; 0xFF invalid; max 236 m/s |
| 72 | 14 | Altitude (WGS-84 ellipsoid) | exponential, unsigned N = 12, × 1 m, **−320 m offset**; 0x3FFF invalid |
| 86 | 9 | Vertical rate | exponential, signed (sign + e2 + N = 6), × 0.125 m/s; 0x1FF invalid |
| 95 | 9 | Ground track | × 360/512 ° |
| 104 | 2 | SIL | 0–3 |
| 106 | 2 | Design assurance | 0 none · 1 D · 2 C · 3 B |
| 108 | 4 | NIC | Rc table, 0–12 |
| 112 | 3 | HFOM (NACp-like) | 0 ≥ 0.5 NM … 5 < 30 m, 6 < 10 m, 7 < 3 m |
| 115 | 2 | VFOM | 0 ≥ 150 m, 1 < 150, 2 < 45, 3 < 10 m |
| 117 | 2 | Velocity accuracy | 0 ≥ 10, 1 < 10, 2 < 3, 3 < 1 m/s |
| 119 | 1 | reserved | |

**Exponential decode:**
- Unsigned: `value = 2^e·(2^N + base) − 2^N`, where `e` = the top 2 bits and `base` = the low N bits.
- Signed: the top bit is the sign.
- Checks from the spec: altitude 0 m → 0x0140, 1000 m → 0x0528; speed 120 m/s → 0xC4; climb 10 m/s → 0x048,
  −10 m/s → 0x148.

**Aircraft category (5 bits):**

| Code | Category |
|---|---|
| 0 | unknown |
| 1 | light fixed-wing |
| 2 | small–heavy fixed-wing |
| 3 | light rotorcraft |
| 4 | glider |
| 5 | lighter-than-air |
| 6 | ultralight / motorised hang-glider |
| 7 | paraglider |
| 8 | parachutist |
| 9 | eVTOL/UAM |
| 10 | gyrocopter |
| **11** | **UAS Open category** |
| **12** | **UAS Specific category** |
| **13** | **UAS Certified category** |
| 14 | model plane |
| 15 | heavy rotorcraft |
| 16 | hang-glider |
| 17 | paramotor |
| 18–31 | reserved |

(Issue 1 OGN code maps OGN UAV → 11.)

**Status payload (type 3, 7 B):**
- Max protocol version (4 bits).
- Receive capability (2 bits each) for M-band, O-LDR and O-HDR.
- ADS-L 4 Mobile (2 bits).
- E-conspicuity Rx bitmap and Tx bitmap (8 bits each). Bits: FLARM, PilotAware, FANET, OGN, 1090, UAT; Rx also has
  Mode-S/PCAS.
- Transponder (2 bits): none, Mode-C or Mode-S.
- Traffic-uplink client flag and FIS-B client flag.
- 24 bits reserved.

**Traffic Uplink (type 4):** up to 10 × (30-bit address + 2 reserved bits + 120-bit Traffic payload). It is sent on
HDR in the uplink slot.

**FIS-B Uplink (type 5):** a DO-358B UAT ground-uplink APDU, up to 217 B.

**Remote ID (type 6):** **reserved for future use.** It is meant to encapsulate an EN 4709-002 Direct-RID message, but
**no format is defined in Issue 2** [D]. A UAS today shows up only through category 11–13 in Traffic.

### Receiver chain (M-band)
1. Run the `flarm.py` PHY with sync 0xF5724B (+ 0x18).
2. Read the length L (24 for Traffic).
3. Manchester-decode L + 1 bytes and count the violations.
4. Check CRC24 over bytes [1 .. L]; it must be 0. Single-bit syndrome correction is optional.
5. If the key index is 0, XXTEA-decrypt (key 0, 6 rounds) the ADS-L Data; if it is 3, skip this step.
6. Parse the fields.

This is **fully implementable now**, with no secrets. A clean public spec is a good candidate for Python then C++.
Testdata: synthesise from the spec. Real captures need an ADS-L device (SkyEcho 2, recent FLARM firmware, OGN trackers
and SoftRF all transmit ADS-L) [U: availability in India nil; 868 is outside the Indian 865–867 MHz band].

---

## C. FANET / FANET+ (Skytraxx; paragliders; LoRa)

**Sources:**
- **Primary [S]:** `3s1d/fanet-stm32/Src/fanet/radio/protocol.txt`, `frame.cpp`, `fmac.cpp`.
- **Licence:** the GitHub API reports **no licence** (null) and the files carry no licence headers, so the code is
  all rights reserved. **Use protocol.txt as a spec only; implement our own parser.**
- **Cross-check [S]:** SoftRF `FANET.cpp` (GPL-3).

### PHY
| Zone | Frequency | LoRa |
|---|---|---|
| EU868 (default) | **868.2 MHz** | **SF7, BW 250 kHz** |
| US920 / AU920 | 920.8 MHz | SF7, BW 500 kHz |
| **IN866** | **866.2 MHz** | **SF7, BW 250 kHz, 14 dBm** (fmac.h zone table, bounding box lat 5–40°, lon 69–89°) |
| KR923 / AS920 | 923.2 MHz | SF7, BW 125 kHz |
| IL918 | 918.5 MHz | SF7, BW 125 kHz |

**LoRa parameters, common to all zones:**
- Explicit header, payload CRC **on**.
- **CR 4/8 when there are few neighbours, 4/5 otherwise.** `fmac.cpp`: `neighbors < MAC_CODING48_THRESHOLD ? CR_8 :
  CR_5`.
- **Sync word 0xF1.** On SX127x that is register value 0xF1, i.e. sync symbols nibble F and nibble 1. The SX126x
  equivalent is 0xF414.
- Preamble: SX127x default 8 symbols [U].
- No LDRO at SF7.
- **Duty cycle:** < 1 %. A tracking packet is about 15–17 B, ≈ 25–30 ms on air [calc].
- SoftRF also has a `FANET_DEPRECATED` build with sync **0x12**, the Semtech default. Accept 0x12 and 0xF1 when the
  frame parses [S].

`lora_phy.sync_word()` already returns the byte, so add **0xF1: "FANET"** to `SYNC_WORDS` and a FANET branch to
`identify()`.

### MAC header (bytes 0–3, plus optional fields)
| Byte | Bits | Field |
|---|---|---|
| 0 | 7 | extended header present |
| 0 | 6 | forward (relay) |
| 0 | 5–0 | **message type** |
| 1 | | **manufacturer ID** |
| 2–3 | | **unique device ID, little-endian** |
| 4 (if ext) | 7–6 ACK (0 none, 1 requested, 2 requested via forward); 5 unicast; 4 signature; 3 geo-forwarded | extended header |
| 5–7 (if unicast) | | destination: manufacturer + ID (LE) |
| next 4 (if signature) | | SHA1-based 32-bit signature over a pseudo-header + payload + PSK |

The payload follows. There is **no extra CRC** beyond the LoRa CRC, and **no encryption**.

**Manufacturer IDs:**

| ID | Manufacturer |
|---|---|
| 0x01 | Skytraxx |
| 0x03 | BitBroker |
| 0x04 | AirWhere |
| 0x05 | Windline |
| 0x06 | Burnair |
| 0x07 | SoftRF |
| 0x08 | GXAirCom |
| 0x09 | Airtribune |
| 0x0A | FLARM |
| 0x0B | FlyBeeper |
| 0x0C | Leaf Vario |
| 0x0D | LYNKX |
| 0x0E | Syride |
| 0x0F | FlyBrief |
| 0x10 | alfapilot |
| **0x11** | **FANET+** (FLARM + FANET devices: Skytraxx, Naviter, Skybean) |
| 0x20 | XC Tracer |
| 0xBA | Burnair |
| 0xCB | Cloudbuddy |
| 0xE0 | OGN Tracker |
| 0xE4 | 4aviation |
| 0xFA | various |
| 0xFB | Espressif base stations |
| 0xFC, 0xFD | unregistered |
| 0xFE | multicast |
| 0x00, 0xFF | reserved |

The IDs 0x0000 and 0xFFFF are reserved.

### Payload types
**Types:**

| Type | Name | Interval |
|---|---|---|
| 0 | ACK | |
| 1 | **Tracking** | about 5 s |
| 2 | Name (UTF-8 string, no terminator) | 4 min |
| 3 | Message (subheader byte + string) | |
| 4 | Service: weather or gateway | 40 s |
| 5 | Landmarks | |
| 6 | Remote config | in development |
| 7 | **Ground tracking** | |
| 8 | HW info (deprecated) | |
| 9 | Thermal | |
| 0xA | HW info v2 | |

**Absolute position:** 6 bytes.
- **Latitude:** bytes 0–2, 24-bit two's-complement LE, `lat = v / 93206`.
- **Longitude:** bytes 3–5, `lon = v / 46603`.
- The encoder rounds `lat·93206` and `lon·46603`. These are the same units as ADS-L.

**Type 1, Tracking:**

| Bytes | Field | Decode |
|---|---|---|
| 0–5 | Position | absolute |
| 6–7 (LE u16) | bit 15 online-tracking; **bits 14–12 aircraft type**; bit 11 altitude ×4 scale; bits 10–0 altitude (m) | type: 0 other · 1 paraglider · 2 hang-glider · 3 balloon · 4 glider · 5 powered aircraft · 6 helicopter · **7 UAV** |
| 8 | Speed | bit 7: ×5 scale; bits 6–0 in 0.5 km/h (max 317.5 km/h) |
| 9 | Climb | bit 7: ×5 scale; bits 6–0 **7-bit two's complement** in 0.1 m/s (±31.5 m/s) |
| 10 | Heading | × 360/256 ° |
| 11 (opt) | Turn rate | bit 7: ×4; bits 6–0 signed × 0.25 °/s |
| 12 (opt) | QNE offset | bit 7: ×4; bits 6–0 signed m |

**Encoder checks** (`app.cpp`):
- Altitude: if > 2047, write `(alt+2)/4 | 0x800`; altitude is clamped to 0–8190 m.
- Speed: in 0.5 km/h units; if > 127, write `(s+2)/5 | 0x80`.

**Other types:**
- **Type 7, Ground tracking:** position (6 bytes), then byte 6: bits 7–4 = type (0 other, 1 walking, 2 vehicle,
  3 bike, 4 boat, 8 need a ride, 9 landed well, 12 need tech support, 13 need medical help, 14 distress, 15 distress
  automatic); bit 0 = online tracking.
- **Type 4, Service:**
  - **Header byte flags:** 7 internet gateway; 6 temperature (+1 B, 0.5 °C signed); 5 wind (+3 B: heading
    360/256, speed and gust in 0.2 km/h, each with bit 7 = ×5); 4 humidity (+1 B, 0.4 %); 3 pressure (+2 B,
    (hPa − 430)·10 LE); 2 remote-config support; 1 state of charge (+1 B, low nibble × 6.67 %); 0 extended header.
  - **Layout:** the position follows the header, then the extra fields in bit order 6 → 1.
- **Type 2, Name:** pilot or device name, useful for identity.

### FANET+
FANET+ is not a separate PHY. Devices with manufacturer 0x11 (Skytraxx 3/4/5, Naviter Oudie/Blade, Skybean Strato)
carry a FLARM radio *and* a FANET radio. Our FLARM detector and the FANET parser report them separately; correlate them
by time and position.

**Implementability:** **now.** The LoRa demodulator already outputs bytes, so this is a parser of about 100 lines.
- **India:** 866.2 MHz is in the IN866 zone, so local paragliding (Bir Billing) may carry FANET. **Priority up for
  India.**
- A bench capture is possible with any SX127x board running GXAirCom or SoftRF in FANET mode (not done).

---

## D. OGN tracker radio protocol (OGNTP, "OGN1")

**Sources [S]:**
- pjalocha `esp32-ogn-tracker`: `main/ogn1.h`, `ognconv.h/.cpp`, `ldpc.cpp`, `rf.cpp`, `freqplan.h`. Its **licence is
  null**, so the code is all rights reserved.
- The same files are copied in SoftRF `libraries/OGN` under **GPL-3**.
- `glidernet/ogn-rf` (the receiver) and `glidernet/diy-tracker` are **GPL-3**.

**Treat all of these as spec or oracle only.**

### PHY
**FLARM-identical:**
- 2-FSK ±50 kHz, **100 kchip/s Manchester (50 kb/s)**.
- `rfm.h`: bit-rate register 0x0140 = 100 kb/s, Fdev 0x0333 = 50 kHz; the RFM chips do the Manchester coding in
  software.

**Frequencies and slots (`freqplan.h`):**
- **Europe/Africa:** 868.2 and 868.4 MHz (2 channels).
- **Slots:** slot 0 at 400–800 ms and slot 1 at 800–1200 ms after the UTC second, as FLARM uses
  (SoftRF `OGNTP.cpp` slot1 = {800, 1200}).
- **Other plans:** US 902.2 + n·0.4 MHz, 65 channels (FLARM-style hopping hash); AU/SA 917.0 MHz + 24 channels;
  NZ 869.25; IL 916.2; EU-433 433.2 MHz × 8.
- **India:** the SoftRF band plan "IN" uses 866.0 MHz (file 10).

**Sync:** 8 chip bytes **`AA 66 55 A5 96 99 96 5A`** = Manchester(**0A F3 65 6C**) [calc; `rf.cpp` comment
"OGNv1 SYNC: 0x0AF3656C"].
- OGNv2 (experimental) uses `55 99 96 59 A5 95 A5 6A` = Manchester(**F5 6D 37 38**).

**Packet on air:** 26 bytes (20 data + 6 LDPC parity), Manchester → 52 bytes, ~4.2 ms; SoftRF says 6 ms air time.
- Bytes go in **memory order of the little-endian uint32 words**, MSB first per byte (as the radio chip shifts them).
  SoftRF `recvBytes()` copies them straight in [U: confirm on a capture].

### Packet layout (20 bytes = 5 LE uint32 words; bit fields LSB-first per word)
**Header word 0:**

| Bits | Field |
|---|---|
| 0–23 | **Address** |
| 24–25 | **Address type**: 0 random · 1 ICAO · 2 FLARM · 3 OGN |
| 26 | NonPos: 0 = position, 1 = status/info/other |
| 27 | **Parity**: even parity over bits 0–27 (`Count1s(word & 0x0FFFFFFF)` even) |
| 28–29 | Relay count |
| 30 | **Encrypted** |
| 31 | Emergency |

**Position packet** (NonPos = 0), words 1–4:

| Word | Bits | Field | Decode |
|---|---|---|---|
| 1 | 0–23 | Latitude (signed) | `lat_deg = ((v<<3)+4) / 600000` (units of 0.0001/60 °) |
| 1 | 24–29 | Time | UTC second mod 60 |
| 1 | 30–31 | Fix quality | 0 none, 1 GPS, 2 DGPS |
| 2 | 0–23 | Longitude (signed) | `lon_deg = ((v<<4)+8) / 600000` |
| 2 | 24–29 | DOP | unsigned VR, 4-bit mantissa |
| 2 | 30 | BaroMSB | negated bit 8 of the baro − GPS altitude difference |
| 2 | 31 | Fix mode | 0 = 2-D, 1 = 3-D |
| 3 | 0–13 | Altitude | `UnsVRdecode<12>` metres |
| 3 | 14–23 | Speed | `UnsVRdecode<8>`, 0.1 m/s |
| 3 | 24–31 | Turn rate | SR2V5, 0.1 °/s; 0x80 = absent |
| 4 | 0–9 | Heading | `(v·3600+512)>>10`, 0.1 ° |
| 4 | 10–18 | Climb rate | SR2V6, 0.1 m/s; 0x100 = absent |
| 4 | 19 | Stealth | |
| 4 | 20–23 | **Aircraft type** | see below |
| 4 | 24–31 | Baro altitude difference | low 8 bits, sign from BaroMSB |

**Variable-range (VR) decoders** (`ognconv.h`):
- **`UnsVRdecode<N>(v)`:** let range r = v >> N and m = v & (2^N − 1). Then:
  - r = 0 → m;
  - r = 1 → 2^N + 1 + 2m;
  - r = 2 → 3·2^N + 2 + 4m;
  - r = 3 → 7·2^N + 4 + 8m.
- **SR2V5 (8-bit):** bit 7 = sign; 7-bit magnitude with range = bits 6–5 and m = bits 4–0:
  - r = 0 → m;
  - r = 1 → 0x21 + 2m;
  - r = 2 → 0x62 + 4m;
  - r = 3 → 0xE4 + 8m.
- **SR2V6 (9-bit):** bit 8 = sign; 8-bit magnitude with range = bits 7–6 and m = bits 5–0:
  - r = 0 → m;
  - r = 1 → 0x41 + 2m;
  - r = 2 → 0xC2 + 4m;
  - r = 3 → 0x1C4 + 8m.

**Aircraft type** (same as FLARM):

| Code | Type |
|---|---|
| 0 | unknown |
| 1 | glider |
| 2 | tow plane |
| 3 | helicopter |
| 4 | skydiver |
| 5 | drop plane |
| 6 | hang-glider |
| 7 | paraglider |
| 8 | powered aircraft |
| 9 | jet |
| 10 | UFO |
| 11 | balloon |
| 12 | airship |
| **13** | **UAV** |
| 14 | ground vehicle |
| 15 | static object |

**Non-position packets** (NonPos = 1; `ReportType` in the top bits of word 4):
- 0 = device status: time, satellites, SNR, pressure, temperature, humidity, voltage, Tx power, noise, firmware and
  hardware versions.
- 1 = info: up to 16 seven-bit characters packed. Fields are "Pilot, Manuf, Model, Type, SN, Reg, ID, Class, Task,
  Base, ICE, PilotID, Hard, Soft, Crew", each terminated by a type code < 0x20. Byte 15 is an XOR check of bytes 0–14.
- 2 = wind/thermal.
- 15 = manufacturer-specific.

### Whitening and encryption
- **The payload is not encrypted by default, but it is whitened:** if `Encrypted = 0`, the transmitter runs
  `Whiten()` on words 1–4.
- **Whiten():** **TEA encrypt with a zero key, 8 loops**, applied separately to the pairs (word1, word2) and
  (word3, word4). The key-0 TEA round is `sum += 0x9E3779B9; v0 += (v1<<4) ^ (v1+sum) ^ (v1>>5); v1 += (v0<<4) ^ (v0+sum) ^ (v0>>5)`.
- **Receiver:** run TEA-decrypt with key 0 and 8 loops on both pairs, starting at sum = 8·DELTA (`Dewhiten()`).
- The header word is never whitened, so the **address, address type and flags are readable before dewhitening**.
- If `Encrypted = 1`, words 1–4 are XXTEA with an owner key (8 loops). Only the header is readable; report
  "encrypted".
- Whitening happens **before** the FEC: the LDPC covers the whitened words.

### FEC: LDPC (208,160)
- **Code:** 160 data bits + **48 parity bits** = 208 (26 bytes).
- **Parity-check matrix:** `LDPC_ParityCheck_n208k160[48][7]` (48 rows × 208 bits). Each row ANDed with the codeword
  must have even weight.
- **Generator:** `LDPC_ParityGen_n208k160[48][5]`.
- **Where defined:** in **pjalocha's `ldpc.cpp`** (esp32-ogn-tracker, no licence) and the GPL-3 copies (SoftRF
  `libraries/OGN/ldpc.cpp`, diy-tracker, ogn-rf). **It is not published in any permissive source or document.**
  [The matrix is numeric data. Whether we may re-key it into a permissive codebase is a licence or legal question for
  the user. Default: don't.]
- **A receiver can verify a packet without the matrix**, more weakly:
  1. Exact **sync match** (32 Manchester-decoded bits, 64 chips).
  2. **Zero Manchester violations** over the 208 data chip pairs. On noise each pair is invalid with probability ½,
     so this is a strong check.
  3. The **header address-parity bit** (1 bit).
  4. Plausibility after dewhitening: time < 60, latitude and longitude in range and near the receiver, fix quality
     ≤ 2, VR fields consistent.
  5. **Consistency across repeats:** the same address at 1 Hz with a moving track.
- **Optional clean-room route [U]:** the code is linear and systematic, parity = G·data over GF(2). About 200 clean,
  linearly independent received packets let us **solve G (48 × 160) from our own captures**. That gives full LDPC
  checking (and parity-based correction) without copying the matrix.

**Implementability:**
- **Now:** detection, address, address type, aircraft type (UAV = 13), position and dewhitening. This reuses the FLARM
  PHY with sync 0x0AF3656C.
- **Pending:** the LDPC check needs the licence decision or the learned G.

---

## E. China GB 46750-2025 drone Remote ID ("运行识别", operational identification)

**Sources:**
- **The standard [D]:** GB 46750—2025 《民用无人驾驶航空器系统运行识别规范》, published 2025-10-31 and in force
  **2026-05-01**. It is a mandatory national standard, and CAAC hosts the full PDF:
  http://www.caac.gov.cn/PHONE/XXGK_17/XXGK/BZGF/GJBZ/202601/P020260120357575326873.pdf
- **Implementation [D]:** `opendroneid/opendroneid-core-c` now contains **`libopendroneidcn/`**, Apache-2.0 (the
  repo LICENSE). It has an encoder, a decoder and a beacon locator; see its README.
- **Oracle and test vectors [S]:** `isaacbentley/orecchino-esp32` `firmware/common/gb46750_decode.h` (GPL-3), whose
  layout and DJI captures come from "Light RID Scanner" (GPL-3).

### Broadcast requirements (§6.1)
- **Transport:** Bluetooth **5.0 or later** broadcast and/or Wi-Fi (2.4 or 5.8 GHz) broadcast. **One transmission must
  carry the whole packet** ("单次完整传输"). The packet is 78 B with all mandatory fields, so BLE needs **BT5
  extended advertising**; a 31-byte legacy advert cannot carry it.
- **Interval:** ≤ 1 s.
- **Network RID:** also mandatory (cellular).
- **ADS-B:** forbidden as a RID means (§5.1.10).
- **Minimum EIRP:**

| Link | Micro | Light and above |
|---|---|---|
| BT | ≥ 1 dBm (or ≥ 3 dBm average) | ≥ 4 dBm (or ≥ 6 dBm average) |
| Wi-Fi 2.4 GHz | ≥ 1 dBm | ≥ 11 dBm |
| Wi-Fi 5.8 GHz | ≥ 9 dBm | ≥ 18 dBm |

**The standard does NOT specify the BLE AD type / UUID or the Wi-Fi IE / OUI.** §5.2.2 says extensions follow
"protocols published by the civil aviation authority".

**Observed in the field [S: Orecchino and Light RID Scanner, DJI 2026 firmware, e.g. Mini 5 Pro]:** the GB packet
rides the **same Wi-Fi beacon vendor IE as ASTM F3411**:
- element 221, OUI **FA:0B:BC**, type **0x0D**, one counter byte, then the GB packet instead of an ODID message pack;
- the first byte **0xFF** cannot start an ODID pack, so the two formats do not collide.

[U: the BLE carriage. Presumably service data UUID 0xFFFA / app code 0x0D in BT5 extended advertising; not confirmed.]

### Packet (§5.2.1, Tables 1–3) [D]
| Byte | Field | Value |
|---|---|---|
| 0 | Data type | **0xFF** |
| 1 | Version | bits 7–5 = **001**; bits 4–0 = sub-version (V1.X). [U: the standard says 0–63 for a 5-bit field; libopendroneidcn uses 0–31] |
| 2 | Data length | bytes of content items only (1–200) |
| 3… | Flag bytes (≥ 3) | 7 presence bits per byte, MSB first. **Bit 0 (0x01) = another flag byte follows** |
| … | Content items | the present items, in ascending ID order, at fixed lengths |

**All multi-byte values are little-endian** (the serial and registration strings are ASCII).

| ID | Flag (byte/mask) | M/O | Bytes | Field | Encoding |
|---|---|---|---|---|---|
| 001 | 1/0x80 | M | 20 | Unique product ID (serial, GB/T 41300) | ASCII, NUL-padded at the end. Pre-2024 aircraft send their serial |
| 002 | 1/0x40 | M | 8 | Real-name registration mark | last 8 characters, ASCII, NUL-padded |
| 003 | 1/0x20 | O | 1 | Operation category | 0 undefined, 1 open, 2 specific, 3 certified |
| 004 | 1/0x10 | M | 1 | UA class | 0 micro, 1 light, 2 small, 3 medium, 4 large |
| 005 | 1/0x08 | M | 1 | Remote-station position type | 0 take-off point, 1 remote-station position |
| 006 | 1/0x04 | M | 8 | Remote-station position | **int32 longitude, then int32 latitude**, × 1e-7 °; 0xFFFFFFFF = unknown |
| 007 | 1/0x02 | M | 2 | Remote-station altitude (ellipsoid) | (h + 1000)·2, 0.5 m; 0 = unknown |
| 008 | 2/0x80 | M | 8 | **UA position** | **longitude, then latitude**, int32 × 1e-7 ° |
| 009 | 2/0x40 | M | 2 | Track | ×10, 0–3599; 0xFFFF = unknown |
| 010 | 2/0x20 | M | 2 | Ground speed | ×10 m/s; 0xFFFF = unknown |
| 011 | 2/0x10 | O | 2 | Height above take-off | (h + 9000)·2; 0 = unknown |
| 012 | 2/0x08 | O | 1 | Vertical speed | bit 7 = sign (1 = descending); magnitude × 0.5 m/s (0–127); 0xFF = unknown |
| 013 | 2/0x04 | M | 2 | Geodetic (ellipsoid) altitude | (h + 1000)·2 |
| 014 | 2/0x02 | O | 2 | Pressure altitude | (h + 1000)·2 |
| 015 | 3/0x80 | M | 1 | Operating status | 0 not reported, 1 ground, 2 airborne, 3 emergency, 4 RID failure (no emergency), 5 RID failure + emergency |
| 016 | 3/0x40 | M | 1 | Coordinate system | 0 WGS-84, 1 CGCS2000 |
| 017 | 3/0x20 | M | 1 | Horizontal accuracy | NACp-like, 0–12 (12 < 1 m) |
| 018 | 3/0x10 | M | 1 | Vertical accuracy | 0–6 (6 < 1 m) |
| 019 | 3/0x08 | M | 1 | Speed accuracy | 0–4 (4 < 0.3 m/s) |
| 020 | 3/0x04 | M | 6 | Timestamp | Unix **milliseconds**, uint48 |
| 021 | 3/0x02 | M | 1 | Timestamp accuracy | 0 > 0.5 s … 8 ≤ 10 ms |

With every item present: 72 content bytes, 78 B in total. The flags are `FF FF FE`.

**Test vector** [S: Light RID Scanner capture of a DJI Mini 5 Pro, quoted in orecchino `tests/odid_test.c`]:
- Hex: `ff2048fffffe 3135383146414e4c433235385530323952544e36 3030303030303030 00 01 01 d1823b48 3bf2eb11 …`.
- Decoded: serial `1581FANLC258U029RTN6`, class light, take-off point; remote station at lon 121.1859665, lat
  30.0675643.
- It confirms **longitude-first order**. Use it to check our parser; regenerate our own capture for `testdata/`
  (licence of the capture data [U]).

### Differences from ASTM F3411 / EN 4709-002
| | F3411 | GB 46750 |
|---|---|---|
| Message | six 25-byte message types and a message pack | **one variable-length packet** with a flag-bitmask chain |
| ID | serial **or** session ID | serial **and** registration mark, both mandatory |
| Operator position | optional | mandatory, plus its type and altitude |
| Coordinates | lat then lon, 1e-7 °, WGS-84 | **lon then lat**, 1e-7 °, WGS-84 or CGCS2000 |
| Altitude | (h + 1000)·2 | same for geodetic, pressure and remote station; height above take-off uses +9000 |
| Timestamp | 1/10 s past the hour | Unix ms (6 bytes) |
| Authentication | Auth message | none |
| UA type | 16-value type list | UA weight class (micro … large) instead |
| BLE | BT4 legacy and BT5 coded/extended | BT5 or later only |
| Wi-Fi | NAN and beacon | "Wi-Fi broadcast"; observed in the beacon vendor IE FA:0B:BC/0x0D |

**Implementability:** **now** (Apache-2.0 reference plus the standard text).
- **Python:** in `core/remote_id.parse_odid`, branch when `body[0] == 0xFF and body[1] >> 5 == 1`. The Wi-Fi beacon IE
  path is shared.
- **BLE:** needs BT5 extended advertising (AUX_ADV_IND on secondary channels) in our BLE receiver [U: check whether
  `rc_link`'s BLE path handles extended adverts].
- **Global relevance:** high. DJI's 2026 firmware emits it, and Chinese-market drones may appear anywhere.

---

## Summary: what is implementable now
| System | PHY in hand | Verification | Decode | Status |
|---|---|---|---|---|
| **A. UAT 978** | new: CPFSK Δφ slicer (dump978 BSD-2) | RS(30,18)/(48,34)/(92,72), 0x187, fcr 120: **write our own RS** (Karn's is LGPL) | full: ICAO/self address, position, altitude, velocity, **emitter category 14 = UAV**, callsign | **implementable now**; test on dump978 sample hex plus synthetic IQ. No live 978 traffic in India |
| **B. ADS-L M-band** | `flarm.py` with sync 0xF5724B(18) | CRC24 = Mode-S 0xFFF409, XXTEA key-0 descramble | full Traffic and Status; **UAS category 11/12/13**; Remote ID payload type 6 **not yet defined** | **implementable now** (public spec). O-band LDR/HDR: PHY clear, minor [U] on LDR CRC8 span and bit order; RS-FEC "primitive element 11" mapping [U] |
| **C. FANET** | `lora_phy.py` (SF7/BW250, sync 0xF1; India 866.2 MHz) | LoRa CRC | full: manufacturer + ID, **type 7 = UAV**, position, altitude, speed, name | **implementable now**; spec repo has no licence, so write our own parser |
| **D. OGNTP** | `flarm.py` with sync 0x0AF3656C | Manchester violations + address parity + plausibility; LDPC(208,160) matrix only in no-licence/GPL code | address, type (**13 = UAV**), position after TEA key-0 dewhitening | **implementable now without LDPC**; full LDPC needs a licence decision or a matrix learned from captures |
| **E. GB 46750** | existing Wi-Fi beacon vendor IE (FA:0B:BC/0x0D) and BLE (needs BT5 extended advertising) | none (plain) | full field set | **implementable now** (standard plus Apache libopendroneidcn); BLE carriage [U] |

## Sources
- dump978 (BSD-2; `libs/fec` is LGPL): https://github.com/flightaware/dump978 (`uat_protocol.h`, `fec.cc`,
  `demodulator.cc`, `uat_message.cc`, `faup978_reporter.cc`)
- ICAO UAT SARPs (ACP-WGW01 Appendix B):
  https://antena.fe.uni-lj.si/literatura/Razno/Avionika/uat/ACP-WGW01-Report%20on%20A.I.1-Appendix%20B%20-%20UAT%20SARPs.pdf
- UAT Implementation Manual v1.5:
  https://antena.fe.uni-lj.si/literatura/Razno/Avionika/uat/WGC09-WP04%20-%20UAT%20Implementation%20Manual-V1-5.pdf
- UK CAA, UAT for UAS under DO-282B: https://www.caa.co.uk/publication/download/26426
- EASA ADS-L 4 SRD-860 Issue 2: https://www.easa.europa.eu/en/downloads/142821/en
- OGN tracker (no licence; spec only): https://github.com/pjalocha/esp32-ogn-tracker (`main/adsl.h`, `ogn1.h`,
  `ognconv.*`, `ldpc.cpp`, `rf.cpp`, `freqplan.h`, `rfm.h`)
- SoftRF (GPL-3; oracle): https://github.com/lyusupov/SoftRF (`src/protocol/radio/{ADSL,OGNTP,FANET,UAT978}.cpp`,
  `libraries/OGN`)
- glidernet ogn-rf (GPL-3): https://github.com/glidernet/ogn-rf · diy-tracker (GPL-3):
  https://github.com/glidernet/diy-tracker
- FANET (no licence; spec only): https://github.com/3s1d/fanet-stm32/blob/master/Src/fanet/radio/protocol.txt
- GB 46750-2025 full text (CAAC): http://www.caac.gov.cn/PHONE/XXGK_17/XXGK/BZGF/GJBZ/202601/P020260120357575326873.pdf
- opendroneid-core-c `libopendroneidcn` (Apache-2.0):
  https://github.com/opendroneid/opendroneid-core-c/tree/master/libopendroneidcn
- Orecchino (GPL-3; oracle and test vectors): https://github.com/isaacbentley/orecchino-esp32 · Light RID Scanner
  (GPL-3): https://github.com/luyii-code-1/Light_RID_Scanner
