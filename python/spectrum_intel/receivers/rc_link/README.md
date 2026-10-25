# rc_link — LoRa / ExpressLRS, GFSK/FSK RC links, BLE advertising, 2 Mb/s ESB-like links

**Detects**
- LoRa / ExpressLRS CSS packets;
- GFSK/FSK frequency-hopping RC and telemetry links (FrSky-like, ELRS FLRC/GFSK-like);
- 1 Mb/s GFSK bursts (BLE);
- 2 Mb/s GFSK packets of the Enhanced ShockBurst family (nRF24 / nRF52 ESB and clones: wireless mice,
  keyboards, game controllers, toy RC), "2.4 GHz GFSK link (ESB-like, 2 Mb/s)".

**Decodes**
- **LoRa:** SF, BW, preamble and payload length; for explicit-header packets also the payload (`lora_phy.py`):
  length, coding rate, CRC-16, the sync word (0x12 private, 0x34 LoRaWAN, 0x2B Meshtastic), and the plain header
  fields of LoRaWAN (message type, DevAddr, FCnt, FPort; JoinEUI/DevEUI of a join request) and Meshtastic
  (sender, destination, packet id, hop limit). Their encrypted payloads are left as they are. `success` when the
  CRC passes.
- **LoRa SF5 (SX126x), implicit header** (`lora_phy.decode_sf5_implicit`): the SFD is followed by 2 extra symbols;
  the first block is CR 4/8 at the full SF (5 nibbles of payload, no reduced rate); then blocks of 4+CR symbols. The
  length and CR are not on air: every (CR, length) whose symbol count is within ±8 of the walked packet is tried and
  only a CRC-16 pass without an uncorrectable FEC error is kept (~36 hypotheses, false pass ~5e-4); two passing
  payloads -> nothing; a last block of zero nibbles (valid at every CR) -> `cr_ambiguous`. Measured on the bench
  (2026-10-19, 890.03 MHz telemetry link: 12/12 packets CRC-valid, and `synth.lora_sf5_implicit_symbols`
  re-encodes each payload to the exact on-air symbols). Explicit-header SF5 and SF6 are not verified and not decoded.
- **Bands:** the full receiver runs when the capture touches an RC band (ISM 433 / 865 / 2.4 GHz, 915–925 MHz);
  a capture that touches none of them but 860–930 MHz runs the LoRa scan only, on sub-bands centred in 860–930 MHz
  (`plugin.LORA_HZ`, `detect(..., {"lora_only": (lo, hi)})`): EU868, FCC915 drone links and off-band telemetry
  radios (the live SF5 link at 890 MHz) sit in GSM/LTE bands, where the FSK detector would mistake uplink bursts
  for RC links.
- **ExpressLRS 2.4 GHz payloads are not decoded:** ELRS sets the SX1280 to the long-interleaving coding rates
  (`SX1280_LORA_CR_LI_4_8` / `4_6`, ExpressLRS `common.cpp`), whose interleaver Semtech does not document; the
  standard chain fails even on the reduced-rate first block (bench 2026-10-20, `LIVE_ELRS_2G4_2440MHz`).
- **BLE advertising:** access address 0x8E89BED6, de-whitening and CRC-24, giving the PDU type and advertiser
  address (`success` when the CRC passes). Without a CRC a 1 Mb/s GFSK burst is "Bluetooth LE advertising-like"
  only on channels 37/38/39 (2402/2426/2480 MHz ± 0.35 MHz); elsewhere it is "Bluetooth (BLE/BR)-like", a BLE
  data channel or a BR/EDR hop (live 2026-10-16, `testdata/live/LIVE_ISM_2452MHz`: 4 hops at 2463–2472 MHz were
  "advertising-like"). The modulation index measured on CRC-valid adverts spans 0.27–0.65, so it cannot tell
  BLE (h 0.5) from BR (h 0.32).
- **Low-SNR BLE decoding (research, 2026-10-21, PLAN L1; `gfsk_coherent.py`, bench `bench/ble_lowsnr`):** at
  8 samples/bit the coherent detector runs first and the hard-decision slicer is the fallback (for the other 39
  whitening channels, i.e. shifted frequency labels).
  1. The access address is found by soft correlation of the one-bit differential phase; the bit timing comes from
     the peak of |Im z[n+T] z*[n]|.
  2. First pass: differential decisions, then h = the value whose modulation-removed products are most coherent
     (per-bit ISI table).
  3. Laurent C0 matched filter. The h walk and the matched-filter ISI are removed; s_k = r_k^2 (-1)^(k+1) is then a
     tone at twice the carrier offset. Its frequency track comes from segment FFT peaks (8, 4, 1 segments, quadratic
     fit): carrier offset plus drift. Phase from a ±16-bit moving sum; decisions on each symbol's axis. Two
     iterations.
  4. CRC-24 Chase correction: up to 3 b-symbol flips (each flips 2 adjacent bits) among the 16 least reliable, plus
     the last symbol; the length field is never flipped.
  - Bit errors at 8 dB in-band SNR: 2e-4, against 2.3e-2 (one-bit differential) and 5e-2 (old slicer).
  - It holds for h 0.45–0.55 and drift up to 100 Hz/µs: ≤ 2.3e-3 at 8 dB.
  - Remote ID adverts decode at **8 dB** (1000-bit message pack) and **6 dB** (376-bit BLE 4 advert); the old
    slicer needed 12 dB.
  - Cost: 5–11 ms per advert in Python, against ~100 ms for the slicer's 16 passes.
  - Decoded keys: `demod` "coherent", `bit_corrections`, `aa_bit_errors`, `mod_index` (measured h, a transmitter
    property), `aa_ncorr`. Plugin field `ble_coherent` (adverts, corrections, median h); tagged research
    (`ble_coherent_demod`).
  - A 2.4 GHz burst 0.5–2.5 MHz wide and ≥ 120 µs long with no measurable symbol rate (below ~6 dB slice SNR) gets
    one coherent attempt too. Nothing is reported without the CRC.
  - Not yet: combining repeated adverts (the counter byte changes between them).

- **ESB-like links (`esb.py`):** CRC-16/CCITT verified frames. Products differ in what the CRC covers and how it
  starts (address bytes, a fixed scrambling pattern), so the packets of one link (same 16-bit sync word) are
  checked together: K = CRC0(m[:E]) ^ m[E:E+16] is the same for every valid packet (CRC0 = init 0), so the CRC end
  E is the span where the most packets of distinct content share one K, and the init at the sync word follows
  (I = E inverse register steps of K). `success` for the packets on that CRC (≥ 3 distinct packets, ≥ half the
  link's packets); fields: `sync_hex` (first 32 bits: one per node), `frame_hex`, `crc_hex`, `crc_span_bits`,
  `crc_init_hex`, `node_frame_period_us`, `link_packets_crc_ok`. Fewer than 3 packets: `partial`.

- **Remote ID over BLE (research, 2026-10-21; `research/drone_transmitters/04_telemetry_remoteid.md`):** every
  CRC-valid advert's AD structures are parsed (`core/remote_id.parse_ble_ad`: flags, local name, manufacturer company
  ID, 16-bit UUIDs, TX power, appearance, service data).
  - Service Data with UUID 0xFFFA and app code 0x0D carries ASTM F3411 (`[counter][message | message pack]`). It is
    parsed by `core/remote_id.parse_odid`, shared with the Wi-Fi Beacon transport.
  - The signal is labelled "Drone Remote ID (Bluetooth LE, <PDU>)". Fields: `remote_id` (merged over the adverts),
    `remote_id_adverts`, `uas_id`, `ble_ad`.
  - **Synthetic results** (`python/tests/test_remote_id_ble.py`, `test_gfsk_coherent.py`; synth `ble_adv_burst` /
    `remote_id_adv_data`): full decode at 20, 12 and 8 dB SNR per advert (coherent path). Realistic BLE 4 adverts
    (one 25-byte message + counter) decode at 6 dB.
  - **BT5 Long Range (research, 2026-10-22; `ble_coded.py`, bench `bench/ble_coded`,
    `research/drone_transmitters/09_ble_coded_droneid_legacy.md`):** LE Coded PHY.
    - Packet: 80-symbol 00111100 preamble; AA + CI + TERM1 at S=8; PDU + CRC at S=8 or S=2.
    - Code: K=4, rate 1/2 (17 / 13 octal); S=8 mapping 0→0011, 1→1100.
    - Receiver: the coded preamble + coded AA (376 known symbols) by soft correlation, coherent soft symbols
      (`gfsk_coherent`), de-mapping, soft Viterbi, de-whitening, CRC-24.
    - Extended advertising PDUs are parsed: `extended` = AdvMode, AdvA, ADI (DID / SID), AuxPtr (channel, offset,
      PHY), TX power; the AdvData carries the F3411 message pack (transport "Bluetooth LE extended advertising (BT5)").
    - Label "Drone Remote ID (Bluetooth LE Long Range, AUX_ADV_IND)".
    - Pipeline gates:
      - S=8 looks like 500 kBd at ~1 MHz width, so a 2.4 GHz burst with an apparent rate of 0.45–0.55 MHz gets
        a coded attempt.
      - A long box that splits into ≥ 2 pieces (low SNR, or an S=2 block after the S=8 one) gets one decode over
        the whole box.
    - Synthetic results (`tests/test_ble_coded.py`): ADV_EXT_IND + an 8 ms AUX_ADV_IND message pack at 12 and 1 dB
      (S=8), 8 dB (S=2).

- **Protocol candidates (research, 2026-10-21; `protocols.py`, `research/drone_transmitters/03_rc_links.md` /
  `04_telemetry_remoteid.md`):** for FSK-family signals, the median symbol rate, deviation and carrier plus the
  shortest burst spacing (≥ 3 bursts) are matched against a table of RC / telemetry protocols. The output is a ranked
  `protocol_candidates` field, each candidate with matched and conflicting criteria (rate must agree within 5 %).
  - **Protocols:** FrSky D8 / D16 / v2 / LBT, RadioLink, Crossfire 150 Hz, mLRS FSK, FLARM, S-FHSS, HoTT, Hitec,
    ELRS GFSK / FLRC, AFHDS 2A, nRF24 / XN297, OpenLRS-ng, XBee 900HP, and SiK at all 13 air rates (deviation ≈ rate).
  - **Labels unchanged:** the existing `*_like` labels stay as they were.
  - **Known tie:** mLRS FSK vs FLARM (same rate / deviation / band); packet length and period would separate them.

- **SiK telemetry + MAVLink (research, 2026-10-21; `sik.py`, `core/mavlink.py`; `research/drone_transmitters/05_sik_mavlink_details.md`).**
  - **Gate:** a sub-GHz GFSK burst whose rate is within 5 % of a SiK air rate (2…250 kb/s) and whose deviation is
    0.6–1.4 × the rate.
  - **Slice:** re-sliced at 8 samples/bit from the burst start over the longest SiK packet, so a packet split by the
    energy detector still decodes.
  - **Decode:** sync 0x2DD4 after a 1010 preamble, in either polarity; PN9 de-whitening, LSB- or MSB-first per byte;
    NETID, length, data.
  - **Verification:** the MAVLink frames inside must pass X.25 + CRC_EXTRA. Five IBM CRC-16 variants are tried on the
    radio CRC and the matching one is reported, since the seed and bit order are not in any reachable source.
  - **Output:** label "SiK telemetry radio (MAVLink)"; fields `sik` (netid, air rate, whitening, hw_crc_variant,
    packets) and `mavlink` (vehicle type, autopilot, armed, flight mode, position, altitude, home, speeds, battery,
    attitude, status text, radio RSSI).
  - **Synthetic results** (`python/tests/test_sik_telemetry.py`): full decode at 20 / 12 / 9 dB; ~0.5 ms per attempt.
    C++ = Python (`c++/tests/vectors/sik`).
  - **Not decoded:** Golay mode (ECC=1) and AES-encrypted links.

**Frequency-label check:** a CRC-valid BLE advert only passes when de-whitened with its true channel index, so
adverts reveal a wrong capture centre-frequency label (`label_check.check_frequency_label`). The pipeline then
re-analyses the snapshot at the corrected frequency.

## How
- **LoRa:** a blind dechirp scan per sub-band (critically sampled, K-symbol coherent FFT), 2.25-downchirp SFD
  verification, a symbol walk, and a match to the ELRS channel grid. The payload chain follows gr-lora_sdr:
  symbol − 1, reduced rate (÷4) in the header block and with LDRO, Gray, diagonal deinterleaving, Hamming
  syndrome correction (4/7, 4/8), the 5-bit header checksum, de-whitening, CRC-16/CCITT. A header without a CRC
  pass is kept only when its symbol count matches the measured burst (the checksum alone passes 1 in 40 random
  implicit-header packets). The two sync symbols must be a sync word's nibbles × 8 (± 1): every real packet in
  testdata reads [8, 16] (0x12). A **preamble-only** report (no SFD confirmed) needs the SFD search window to run
  past the end of the capture, or the SFD itself to be seen with only the payload cut short: a real preamble is
  followed by its SFD. (2026-10-10: these removed 6 "LoRa / ExpressLRS preamble SF6" made of the Pluto's capture
  self-noise at 424–429 MHz, a "LoRa SF6" with sync symbols [0, 5] in the LTE B5 uplink and a "LoRa preamble SF7"
  inside an LTE downlink carrier; all 13 real packets in testdata are unchanged.)
- **GFSK/FSK:**
  - matched-area energy at 3 widths;
  - narrowband flank and constant-envelope tests;
  - instantaneous-frequency clusters for the deviation;
  - symbol rate from transition timing.
- **BLE:** preamble + access address, de-whitening for channels 37/38/39, CRC-24.
- **ESB-like:** a measured 2.4 GHz burst 0.5–3 MHz wide (40 µs – 1.1 ms, ≥ 6 dB) is re-sliced ±1.25 MHz at
  16 MSps from an FFT of the burst ± 50 µs only, cropped to the burst −2 / +6 µs (neighbouring packets would bias
  the threshold); bits at the widest-eye phase; a run of ≥ 14 alternating bits starting in the first 32 bits is
  the preamble (a 16-bit preamble measures 14–15; the bench link sends ~20 lead-in bits). The carrier is the
  discriminator's mid-level (the power centroid wanders ± 0.2 MHz); signals are grouped on the 1 MHz raster. A
  span one bit short of the true one has a signature (K in two classes 0x8810 apart) and is skipped; ties go to
  the shortest span (a fixed bit after the CRC makes E + 1 constant too).
- **Short packets in long boxes:** `fsk_scan`'s 800 µs integration smears a 1 kHz train of ~70–100 µs packets into
  ~870 µs boxes, where `fsk_refine` took the packet's level (75th percentile of the box) for noise and dropped it.
  The level then comes from the 98th percentile, if that is ≥ 4 × the noise (a noise-only box stays below 1.3 ×).

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.receivers.rc_link import detect, check_frequency_label
ems = detect(x, fs, fc)
print(check_frequency_label(ems))    # {"offset_mhz": ..., "suspect": bool, ...} or None
```
`synth.py` generates synthetic LoRa / GFSK links.

## Validation
- `tests/receivers/test_det_rc_link.py`: LoRa and GFSK Pd vs SNR on noise and on real backgrounds. It is long;
  use `--quick --procs 2`.
- `tests/receivers/test_det_lora.py`: payloads byte-exact against gr-lora_sdr's **transmitter** (SF7–SF11, CR
  4/5–4/8, CRC on/off, LDRO, two carrier offsets; `gen_lora_vectors.py` regenerates them under radioconda) and a
  real Pluto packet at 434.0 MHz whose payload gr-lora_sdr's **receiver** decoded (`testdata/lora/`).
- Live 433.92 MHz capture (2026-10-05, 25 × 0.8 s): a SF9/125 kHz, sync 0x12, CR 4/7 device at 434.000 MHz
  sending 18-byte packets; both decoders agree byte for byte. Only SF9 + 0x12 decodes; SF8/SF10 and sync words
  0x34/0x2B give nothing (controls).
- `tests/receivers/test_det_lora.py`: the live SF5 packet (`testdata/live/LIVE_LORA_SF5_890MHz`, exact payload, CR,
  length) and synthetic implicit SF5 packets at every CR (1–101 bytes, carrier offsets, noise). Synthesise SF5 at
  >= 4 MSps: at 2 MSps the chirp wraps alias in-band and the payload peaks drop ~60 % (an artefact of the synth).
- Walk end rule (2026-10-20): a following packet is recognised by `run` repeated symbol values, run = 3 for
  N >= 128, 4 for N = 64, 5 for N = 32 (a 3-run cut about half of all 85-symbol SF5 packets short).
- `synth.lora_packet` sends sync symbols as nibble x 8 at every SF (it used N/16, wrong except at SF7).
- Latency (2026-10-20): SF5 stage-1 threshold `ETA_SF5` = 40 (with 24, DMR/NXDN at 860 MHz gave ~170 false SF5
  candidates per 1.75 s dwell, 12 s of refinement; 40 is still ~-8 dB in-band SNR, SF5 decodes from ~0 dB where
  the stage-1 ratio is >= 108). `lora_refine(cache=)`: candidates of one sub-band share its full-capture slices
  (`_slice` is most of a refinement; 200 -> 87 slices on `LIVE_PMR_860.4MHz`, results bit-identical).
- Fixed on the way (2026-10-06): payload symbol values were off by N/4 (the payload grid sits 2.25 symbols after
  the SFD, a quarter symbol off the preamble grid); the CFO/timing N/2 ambiguity is now settled by walking both
  hypotheses; a stage-1 window straddling the preamble start no longer loses the packet; preamble-only ghosts one
  BW from a confirmed packet are dropped; FSK refinement and the image test work on a local segment (a 2.5 MSps
  block went from 389 s to under 40 s); the BLE decode is tried only on 2.37–2.52 GHz bursts of >= 120 µs (noise
  fragments at 433 MHz triggered 20,000 8-MSps re-slices per block; the P4/Herelink label corrections are
  unchanged); LoRa preamble-only detections were shown as "Bluetooth LE advertising-like" ("preamble" contains
  "ble"); FSK protocol labels need at least 64 symbols (90 µs blips were "SiK-like": 423 per 0.8 s block at 433 MHz,
  now 0). Trade-off, measured with `test_det_rc_link.py --quick`: FrSky-like is labelled from 12.5 dB instead of
  7.5 dB, because below that a packet is measured in 9–22-symbol fragments that rate agreement cannot tell from
  the artifact blips; it is still detected from 5 dB (as `nb_burst`). The other 15 checks are unchanged, and
  "LoRa 2g4_150hz embedded in p4" now passes. The three remaining FAILs (2g4_50hz Pd at −10 dB, ELRS in the
  Herelink capture, cfo_ppm p90 1.76 ppm) fail identically with the previous code.
- GSM bursts are not RC links (2026-10-16): GSM-900 traffic carriers (no FCCH, so the GSM receiver does not claim
  them) gave 554 µs constant-envelope GMSK bursts that the FSK classes labelled "sik_like" (957.0 / 958.0 MHz,
  `testdata/live/LIVE_GSM_958MHz_fsk`). A burst of 1–8 GSM time slots (k × 576.9 µs ± 60 µs) centred on the 200 kHz
  raster (± 15 kHz) inside a GSM-900 / 1800 uplink or downlink band is dropped (`_gsm_burst`). 865–867 MHz and the
  2.4 GHz band are not GSM bands, so RC and telemetry links there are unaffected.
- Label check: corrected the Phantom 4 capture by +20.0 MHz (8 adverts) and the Herelink capture by −10.0 MHz.
- ESB-like links (2026-10-17, `testdata/live/LIVE_ESB2M_2440MHz`, live LO 2440): a 1 kHz link on 2426 / 2450 MHz,
  two nodes 145 µs apart, found and CRC-checked with no per-device setting: CRC end at bit 124, init 0xDD19, frame
  period 999.97 µs; 44 of 47 packets CRC-valid (3 are inside Wi-Fi ch 10). On the bench the same CRC closed 86 of
  88 packets over three blocks. Before: 1–4 MHz `unknown` entries. Synthetic (`tests/test_core.py`): standard nRF
  ESB frames (5-byte address, PCF, CRC init 0xFFFF) give init 0xFFFF at the address and the exact span; random CRC
  fields never pass; a 1 kHz link on one channel for a whole block passes 16/16 at 15 dB in-band SNR, 7/16 at
  12 dB, none at 9 dB (before the box-level fix, a link that stays on one channel for the whole block was lost).

## Research additions 2026-10-22 (synthetic validation; `schema/capabilities.json` research)
- **ExpressLRS identity:**
  - **2.4 GHz:** the SX1280 inverts IQ when UID[5] bit 0 = 1 (and always while binding) → `elrs_identity`. The
    sub-GHz SX127x driver never inverts.
  - **Sub-GHz SYNC** (`elrs.py`, `lora_phy.decode_implicit`):
    - Implicit-header LoRa payload: the first block is CR 4/8 at reduced rate, then 4+CR blocks. Validated against
      gr-lora_sdr implicit-header packets: SF7–9, CR 4/7 and 4/8, 8 and 13 B (`tests/receivers/gen_lora_implicit.py`,
      `testdata/lora/tx_implicit_*`).
    - OTA4 / OTA8: a SYNC packet's CRC-14 (0x2E57) / CRC-16 (0x3D65) is seeded by its own UID4 / UID5, for OTA
      version 4 or 3 and in binding mode. A pass verifies the packet and gives the hop index, nonce, rate, TLM ratio,
      UID4 / UID5 (+UID3 on 3.x) → `elrs_sync`.
    - Model-match SYNC packets are not searched (they would raise the false-pass rate) and stay unverified.
    - Synthetic: SF7 / 8 / 9 SYNC packets in a capture decode at 10 / 3 / 0 dB (`tests/test_elrs_sync.py`).
    - **Hop recovery** (`elrs.fhss_sequence`, `uid23_candidates`, `recover_hops`; detector post-pass `_elrs_hops`;
      research file 07 §4, PLAN R2d):
      - **Sequence:** seed = UID[2..5] ^ OTA version. The LCG (214013, 2531011, mod 2³¹, output >> 16) shuffles each
        block of `count` entries, whose first entry is the sync channel. It matches a transliteration of the
        reference on 350 seeds (`bench/elrs_fhss/oracle_check.py`).
      - **Search:** each packet after the SYNC gets its nonce from the packet interval and its sequence index from
        the hop interval. UID2 / UID3 are then searched (vectorised; the LCG jumps to each observed block; UID2 bit 7
        never matters).
      - **Domains:** overlapping grids (IN866 inside EU868, AU915 inside FCC915) are all tried with the on-grid
        packets, and a domain is kept when candidates remain.
      - **Output:** `elrs_sync.hop_recovery` = observations, candidate count, the first candidates, the domain, the
        next channels when the candidates agree.
      - **Limit:** one 34 ms block holds few hops, so the UID stays ambiguous (tens of candidates on FCC915);
        several blocks would narrow it.
    - **SF6 (900 200 Hz, PLAN R2e):** the SX127x SF6 first block is reduced-rate as at SF ≥ 7 (SX1276 datasheet, research
      file 12 §6), so `decode_implicit` accepts SF6. gr-lora_sdr SF6 packets (`testdata/lora/tx_implicit_sf6_cr3_len8`)
      decode byte-exact, and a 200 Hz SYNC packet decodes in a capture at 10 dB. The SX126x SF5 / 6 framing (2 extra
      symbols, full-rate first block) is not decoded.
- **FLARM beacons** (`flarm.py`; research file 10 §5, PLAN T5): detection and the clear header only.
  - **Gate:** a 100 kchip/s (±6 %) FSK burst with a deviation of 30–70 kHz in 865.5–866.5 or 867.9–868.7 MHz, at least
    128 chips long → ±150 kHz slice at 8 samples per chip.
  - **Decode:** Manchester pairs (01 → 1, 10 → 0), sync F5 31 FA B6 (≤ 2 bit errors; it matches the published
    Manchester bytes 55 99 A5 A9 55 66 65 96), then the 24-bit radio ID (little-endian) and the address type.
  - **Reliability:** no CRC verifies the header, so the ID is reported only when its 32 bits have no Manchester
    violation (a single chip error always shows as 00 / 11). Every phase, polarity and pair alignment is tried and
    the hit with the fewest violations wins. The obfuscated payload is not decoded.
  - **Results:** synthetic (`tests/test_flarm.py`): ID 40/40 at 10 dB, 39/40 at 6 dB, 3/40 at 3 dB (detected without
    an ID otherwise); never a wrong ID; 0.5 ms per burst.
- **mLRS** (`mlrs.py`; research file 12 §1, PLAN R3 / T6):
  - **Identification:** a LoRa burst with preamble 12 (±3), normal IQ and a payload length near 91 bytes for its SF / BW
    (`_mlrs_match`; the SX127x formula, ±10 % + 4 symbols for the SX126x / SX1280 framings) → `mlrs_like_lora_…`.
    ExpressLRS sub-GHz uses preamble 8 / 10 and 8 / 13-byte payloads, and inverts IQ for half of the UIDs.
  - **Decode:** the 19 Hz "7x" mode (SF6, SX127x framing) decodes with `lora_phy.decode_implicit` (91 bytes, CR 4/5).
    CRC-16/MCRF4XX over bytes 0–88 is tried with both byte orders and both uplink crc1 coverages, and the passing
    convention is reported → `mlrs_lora_…`, FrameSyncWord (16-bit per-binding ID), status bytes. SX126x SF5 / 6 and the
    SX1280 long interleaving are not decoded.
  - **Results:** synthetic (`tests/test_mlrs.py`): the SF6 frame decodes at 10 dB, burst 22.3 ms (as the firmware's air
    time); the SF5 2.4 GHz frame is identified as mLRS-like.
- **TBS Crossfire 150 Hz** (`crossfire.py`):
  - **Detection:** an 85.1 kBd, ±42.48 kHz FSK burst in 850–935 MHz → 24-bit sync 0x3E6B3E (≤ 2 bit errors) after
    the alternating preamble, either polarity.
  - **Packet read:** 23-byte uplink or 13-byte downlink, raw and PN9 de-whitened (MSB first), plus the 260 kHz grid
    channel and link. The read starts at the piece holding the sync and extends over the longest packet.
  - **Fragmentation:** a fragmented box gets one whole-box attempt.
  - **Limit:** the per-hop CRC init is unknown, so packets are reported but not verified.
  - **Results:** label "TBS Crossfire-like (150 Hz FSK, sync 0x3E6B3E)", field `crossfire`; synthetic 3/3 packets
    at 20 and 12 dB (≤ 3 bit errors at 12 dB).
- **FrSky D8 / D16, ACCST** (`cc2500.py`; research file 10, PLAN R6a):
  - **Gate:** a 2.37–2.52 GHz FSK burst whose rate is within 8 % of 31.04 k (D8), 70 k (D16 v1), 76.97 k (D16 v2) or
    99.98 k (D16 LBT) and whose deviation is 15–95 kHz → re-slice at 8 samples / bit, read over the longest packet.
  - **Identification:** CC2500 sync (D16 0xD391D391 tried first, then D8 0xD391; ≤ 2 bit errors, either polarity,
    after the preamble), length byte (D8 0x11; D16 0x1D or LBT 0x20), CRC. D8: hardware CRC-16 0x8005 / 0xFFFF over
    length + payload. D16: software CRC (reflected-0x1021 table applied MSB first, seed 0) over bytes 3 … L−3; v2
    adds the hardware CRC after it. Every timing phase and hit is tried; the first CRC pass wins.
  - **Output:** variant (D8, D16 v1 / v2 / LBT v1), TX ID (bytes 1–2), packet bytes → label `frsky_d8` / `frsky_d16`,
    pretty "FrSky D16 v2 (ACCST)" ("FrSky D16-like" without a CRC pass), field `frsky` (packets, CRC passes,
    variants, TX IDs). An exact 32-bit D16 sync with another length is reported as a sync-only match.
  - **Results:** synthetic, `tests/test_frsky.py`: 20/20 down to 8 dB (D16 18–19/20), 14/20 at 6 dB; 0.2–0.4 ms per
    burst; noise gives nothing. Not yet seen on a real FrSky transmitter.
- **FlySky AFHDS 2A, A7105** (`afhds2a.py`; research file 10 §3, PLAN R7):
  - **Gate:** a 2.37–2.52 GHz FSK burst at 450–550 kBd with 120–260 kHz deviation, ≥ 64 symbols, not claimed by
    BLE → re-slice at 4 MHz, read over the longest packet (32 + 32 + 560 bits).
  - **Identification:** the 32-bit ID code 0x5475C52A after the preamble, ≤ 3 bit errors (the chip's own tolerance),
    either polarity → label `afhds2a_id`, "FlySky AFHDS 2A-like (A7105 ID code)".
  - **Decode:** 38 payload bytes + CRC-16, (7,4) Hamming per nibble. The bit layout and the CRC polynomial / init are
    not verified in the datasheet, so 2 layouts × 2 nibble orders × 4 CRC variants are tried; a decode is reported
    only on a CRC pass with a known packet type (0x58 / 0x56 / 0xAA / 0xAC / 0xBB / 0xBC) and names the convention
    that passed → label `afhds2a`, "FlySky AFHDS 2A (A7105)", packet type, TX ID (bytes 1–4), RX ID (bytes 5–8,
    FF×4 = binding), FEC corrections; field `afhds2a`.
  - **Results:** synthetic (`tests/test_afhds2a.py`): decode 20/20 at 12 dB, 16–20/20 at 8 dB; 0.3–0.8 ms per burst;
    noise gives nothing. A real capture will confirm (or rule out) the conventions.
- **XN297 toy drones and Bayang** (`xn297.py`; research file 10 §4, PLAN R8):
  - **Gates:** the 1 Mb/s 2.4 GHz branch when BLE does not decode. Also any undecoded constant-envelope 2.4 GHz burst
    of 100–600 µs and 0.3–1.6 MHz, because the symbol-rate estimate fails at h ≈ 0.3 (one attempt, ~0.2 ms).
  - **Decode:** prefix 0x710F55 (≤ 1 bit error, either polarity), then the 5-byte address (reverse byte order,
    scrambled) and the 15-byte payload (bit-reversed, scrambled). The scramble table comes from a 7-bit Fibonacci LFSR
    (taps 6 and 2, seed 0x71). It was found independently: it reproduces the published start E3 B1 4B EA 85 BC E5 66.
  - **Verification:** Bayang verifies itself (header A1–A6 / 0x53, byte 14 = 8-bit sum) → label `bayang`, "Bayang
    (XN297)". Bind packets give the TX ID and 4 hop channels; data packets give the sticks and flags. The CRC length
    constant is unverified, so `crc_residual` is reported; it should be the same for every packet of one length.
    Another payload protocol → `xn297_like` with the address.
  - **Demodulation (PLAN L1c):** first the phase change over each whole bit (CPFSK at h ≈ 0.3), at every timing
    phase, with the decision level set from the 24 known prefix bits; then the discriminator slicer as a fallback.
  - **Results:** synthetic (`tests/test_xn297.py`): clean 20/20; 40/38/30/10 of 40 at 9/7/5/3 dB in 8 MHz (≈ +7 dB
    in-band). The half-bit slicer alone gave 38/33/10/1. 0.3 ms per packet, 0.7 ms per noise burst.
- **IEEE 802.15.4, XBee / Zigbee** (`ieee802154.py`; research file 10 §6, PLAN T4):
  - **Gate:** a 2.4 GHz burst on an 802.15.4 channel (2405 + 5k MHz, ±300 kHz), 1.2–4 MHz wide, ≥ 250 µs, ≥ 3 dB. It
    is tried before the 2 Mb/s ESB slicer; a frame found there skips ESB.
  - **Receiver:** a ±1.6 MHz slice at 8 MHz (4 samples per chip). The FM discriminator is summed over each chip. The
    16 symbol templates are made by modulating each 32-chip sequence (half-sine O-QPSK) through the same
    discriminator, so no analytic MSK sign rule is assumed. Timing comes from the preamble template, then the SFD
    (symbols 7, A), the PHR length and the symbols by template correlation.
  - **Output:** FCS CRC-16 KERMIT pass → label `ieee802154`, MAC header (frame type, sequence, PAN IDs, addresses);
    a 64-bit address with the Digi OUI 00 13 A2 → vendor "Digi International (XBee)". No CRC pass →
    `ieee802154_like`. Field `ieee802154` (frames, CRC passes, channels, PANs, addresses, vendors); model block
    (XBee 802.15.4 / Zigbee 3.0 / DigiMesh share the PHY → unresolved).
  - **Results:** synthetic (`tests/test_ieee802154.py`): 10/10 frames down to ~7 dB in-band, 8/10 at ~4 dB; 0.7 ms per
    frame; 3/3 XBee frames through `detect` at 12 dB.
- **Blind frame structure of any other FSK packet** (`fsk_frame.py`; PLAN F1):
  - **Gate:** an FSK burst with a symbol rate of 20 kBd–1.2 MBd and a deviation, ≥ 10 dB, long enough for 12 + 32 + 16
    bits, that no protocol decoder (BLE, SiK, Crossfire, FrSky) claimed.
  - **Method:** hard decisions at 8 samples / bit at every timing phase; the longest alternating run near the start is
    the preamble. The run can absorb the sync word's first bits when they continue the alternation, so the 32-bit
    candidate may start a few bits into the sync word.
  - **Output:** `features.frame` = preamble bits, last preamble bit, 32-bit sync candidate, payload bits; `known_sync`
    when a published sync word (FrSky D16 / D8, A7105 ID code, Crossfire, SiK) is found within ±4 bits of the
    preamble end, either polarity (≤ 1 bit error for 24 / 32-bit words, exact for 16-bit). Per signal, field
    `fsk_frame`: packets, the top four 24-bit sync candidates with counts, preamble lengths, known sync names.
  - **Results:** synthetic (`tests/test_fsk_frame.py`): the same unknown sync gives the same 24-bit candidate on 3/3
    hops; the A7105 ID code is named at 18 dB. Nothing is decoded or verified ("-like" names only).


## Market-gap batch 2026-10-22 (PLAN G2-G6; synthetic validation, `schema/capabilities.json` research)
Specs from research file 14 (DIY Multiprotocol / ExpressLRS sources read as GPL references only; implemented
independently).
- **RadioLink** (`radiolink.py`): CC2500 GFSK 37.99 kBd ±38.09 kHz on f = 2401.99994 + n × 0.3999 MHz; 4 × 0xAA, the
  sync word twice (0xABCD surface / air, 0xBADC DumboRC, 0xABA5 RC4G), PN9-whitened length + payload + CRC-16 (0x8005 /
  0xFFFF). Data packet (33 bytes): header (0x01 surface / 0x21 air, | 0x02 telemetry request), TX ID, 16 × 11-bit
  channels, hop counter, replaced slot, sum checksum; the TX ID's 16-channel hop table (12 i + inc, two LCG shuffles)
  predicts the packet's channel (`hop_check`). Receiver telemetry and RC4G parsed. Detector: a 38 kBd GFSK burst at
  2.4 GHz; the 9.3 ms packet's later pieces join the decoded emission. Model: air → T8FB / T8S, surface → RC4GS V3 /
  RC6GS V3 / RC8X (the AT9S / AT10 protocol is not published).
- **FlySky AFHDS 1** (in `afhds2a.py`): same A7105 PHY and ID code as AFHDS 2A; 21 bytes, no FEC: 0x55 / 0xAA, TX ID
  (LE), 8 channels (µs); reported when a CRC convention passes, or when the header, all channel values and the TX ID's
  hop table (8 × 4 nibble table) containing the measured channel agree ("-like"). TX IDs in the field come from
  CRC-verified packets when there are any.
- **Hubsan** (`hubsan.py`): A7105 2-FSK 100 kb/s ±186 kHz, 12 channels 2410–2465 MHz; the session ID code is read where
  the preamble ends; 16 bytes summing to 0 mod 256 plus a known header and, for data, a model's constant bytes (H107,
  H107 Plus, H301, H501 class). Gated on deviation and burst length (the symbol-rate estimate reads 35–100 kHz at h ≈ 3.7).
- **Toy drones** (`toys.py`): XN297 protocols (MJXQ, CX10, JJRC345, H8_3D, E016H, Potensic) matched on the XN297
  decoder's address + payload (10–21 bytes now read); plain nRF24 (V2x2, SymaX, SymaX X5C) by correlation with the
  fixed address, then nRF CRC-16 over address [+ 9-bit ESB control field] + payload and the protocol checksum. The
  nRF24 path uses a full-bit phase demodulator (an isolated bit is missed by the half-bit discriminator at h ≈ 0.3).
- **ExpressLRS LR1121 GFSK** (`elrs.gfsk_try_decode`): 300 kb/s ±100 kHz, 16-bit 0x55 preamble, sync UID5 UID4, fixed
  OTA4 / OTA8; a SYNC packet verifies itself (CRC with its own UID bytes) and its UID bytes must equal the sync word.
  The whitening bit order is not published: the CC2500-style PN9 bytes and their bit-reversed form are tried.
- **MAVLink over transparent LoRa** (`lora_phy.mavlink_in`): any explicit-header LoRa payload with CRC_EXTRA-verified
  MAVLink frames → "MAVLink over LoRa (transparent serial link)" (Ebyte E22 / E32 / E220-class links).
- Shared: `sik.demod_bits_mid` / `two_level_threshold` (a 2-means decision level for unbalanced payloads, used by
  Hubsan and the nRF24 toys); decodes of the new protocols found on one piece of a split burst are carried to the
  merged burst (`_CARRY_KEYS`).
- Not covered yet: 250 kb/s XN297 toys, HS6200 / RF2500 toy chips, ELRS Gemini Xrossband pairing (PLAN G3b / G6b).

## Second market-gap batch 2026-10-23 (PLAN G10 / G12; synthetic validation, `schema/capabilities.json` research)
Specs from research files 16 and 17 (DIY Multiprotocol, SoftRF and the OGN tracker read as GPL / unlicensed references
only; the EASA ADS-L spec, TI DN504 and the FANET protocol text as specifications; implemented independently).
- **Graupner HoTT** (`hott.py`): CC2500 MSK 249.94 kBd, sync D391 twice; the CC2500 FEC (K = 4 rate-½ convolutional
  code, 4 × 4 symbol interleaver, terminator 0x0B) decoded by a hard-decision Viterbi that reproduces the TI DN504 test
  vector; PN9 de-whitening; CRC-16. TX packet (50 bytes): hop-table tag (one of 16 tables), hop index, 12 channels in
  0.5 µs, TX ID, RX ID (bind: FA × 5); receiver telemetry (22 bytes, sync 2C6E in "Sync" mode). Both MSK polarities.
- **Hitec Optima / Minima** (`hitec.py`): CC2500 MSK 249.94 kBd, 24-byte preamble, sync word = the TX ID (id3 id2 sent
  twice; bind D391): every 16-bit word sent twice after the preamble is a candidate; PN9, variable length, CRC-16, and
  the sync must equal the TX ID's bytes 3, 2. Data (9 channels, repeat flag), bind (hop table / 2), telemetry.
- **Futaba S-FHSS** (`hitec.py`): CC2500 2-FSK 128.14 kBd ±38 kHz, sync D391 (32 or 16 bits), 13 bytes, no whitening,
  CRC-16 in CC2400 mode (init undocumented: 0xFFFF then 0x0000; the passing one is reported). TX ID, model-hint bytes
  (T8J / T10J / TM-FH), hop index checked against the channel, hop code, 4 × 12-bit channels, command bits.
- **ADS-L SRD-860 M-band** (`adsl.py`): the FLARM PHY with sync F5 72 4B + length; CRC-24 (Mode-S polynomial) over
  header + data; XXTEA with the zero key, 6 rounds (checked against a reference implementation); Traffic fields with the
  spec's exponential encodings (checked against its examples); aircraft category 11–13 = UAS.
- **OGN tracker** (`ogn.py`): the FLARM PHY with sync 0A F3 65 6C; the LDPC (208,160) matrix is not used (unlicensed /
  GPL only): a packet needs zero Manchester violations over all 208 coded bits, the header's even parity and a
  plausible de-whitened position (TEA zero key, 8 loops; checked against a reference implementation): "-like".
- **FANET** (`lora_phy.fanet_parse`): CRC-valid LoRa payloads with sync 0xF1 (the deprecated 0x12 only for a strictly
  valid tracking packet): header, manufacturer + device ID, tracking (aircraft type 7 = UAV, position, altitude, speed,
  climb, heading), name, ground tracking, service position.
- Detector: the FLARM-PHY decoders try a piece of ≥ 40 chips (a packet split into pieces is read from the piece holding
  the preamble + sync; later pieces join the decoded emission); HoTT / Hitec run on 200–310 kBd, 30–110 kHz bursts of
  0.6–4.4 ms at 2.4 GHz, S-FHSS on 100–160 kBd, 20–60 kHz bursts of 0.4–1.9 ms below 2450 MHz.
- Not covered: FrSky R9 (PLAN R4, dropped), DSM / DEVO / M-LINK decode (PLAN R5: timing only), ADS-L O-band.

## Model level (`fields.model`, research, PLAN M5)
`models.py` applies the rules of `research/drone_transmitters/11_model_level_identification.md` §4–5 to what the
detector already decoded or measured. ELRS: the radio chip class from the air mode, plus the domain from the channel
grid. FrSky: the variant (v1 FCC / v2 / LBT = firmware generation or region build). AFHDS 2A, Crossfire: transmitters
differ only in power, so the candidates are listed and the model stays unresolved. SiK: the air-rate set separates
RFD 3.x radios from the 100 mW boards. Tests: `tests/test_drone_models.py`; C++ == Python on `vectors/rc_models`.

## Signal gaps S7 (2026-10-25; research, synthetic validation)

**Bluetooth LE trackers and notable devices** (`core/ble_trackers.py`, called by `core/remote_id.parse_ble_ad`):
byte-level recognizers on a CRC-valid advert's AD structures with strict length / constant checks: Apple Find My
(AirTag, AirPods, Apple devices; battery, maintained), iBeacon, AirPods proximity pairing (model), Apple Continuity,
Samsung SmartTag (state, battery, UWB), Google Find My Device Network (FMDN) and Eddystone, Google Fast Pair, Tile,
Chipolo, Microsoft Swift Pair, IETF DULT. `ad.device` = {vendor, kind, label, tracker, …}; the signal carries
`ble_devices` and an item tracker is labelled "Bluetooth LE tracker (<label>)". No stable identity is recovered: the
trackers broadcast rotating tokens.

**Bluetooth BR/EDR access codes** (`bredr.py`): on a 1 Mb/s GFSK hop that is not a BLE advert (nor XN297 / nRF24),
±0.6 MHz filter, discriminator integrated over one symbol, threshold between its 10th / 90th percentiles, 8 phases,
both polarities; a 64-bit window within 3 bits of the (64,30) BCH codeword of its own LAP (PN 0x83848D96BBCC54FC,
generator 0o260534236651; the LAP 0 and GIAC codewords reproduce the published ones) is an access code (chance
match ~2.5e-5 per burst). Label "Bluetooth BR/EDR (LAP <hex>[, inquiry])" (GIAC / LIAC = inquiry), 0.9, status
`partial`; fields `bredr` {laps, packets, inquiry, max_access_code_errors}. The UAP needs the master clock and
several packets, and the NAP is never sent: no vendor is claimed. Synthetic: 6/6 at 20 dB, 4/6 at 12 dB in-band
(below that the 1 Mb/s rate gate drops the hop, as for BLE); 0 LAPs on random GFSK bursts; ~1 ms per hop. Real: hops in the Phantom 4 testdata capture
(`testdata/DJI_Phantom_4_20251119-182450`) carry BR/EDR access codes that match with 0 bit errors (LAP DA7476,
B4BDD9): Bluetooth Classic devices near the bench.
