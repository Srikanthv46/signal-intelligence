# 09 — BLE Coded PHY (Remote ID BT5 Long Range) and DJI DroneID 8 vs 9 symbols

Research agent, 2026-10-21; saved from its report (the agent's own write failed on a full disk). [U] = uncertain.

## Part 1: BLE LE Coded PHY (Core 5.x Vol 6 Part B 2.2 / 3.3) for ASTM F3411 BT5 Long Range
- **Packet:**
  - Preamble: 80 uncoded symbols (`00111100` × 10).
  - FEC block 1 = AA (32 b) + CI (2 b) + TERM1 (3 b), always S=8.
  - FEC block 2 = PDU + CRC24 + TERM2: S=8 if CI=00, S=2 if CI=01.
  - Modulation: 1 Msym/s GFSK, BT 0.5, h~0.5.
  - Air time: 592 + 64N µs at S=8, 430 + 16N µs at S=2 (N = PDU octets incl. the 2-octet header).
- **FEC:**
  - Rate 1/2, K=4: G0 = 1+x+x²+x³ (octal 17), G1 = 1+x²+x³ (octal 13). Initial state 0; TERM = 3 zero bits.
  - Output a0 (G0) first [U].
  - Pattern mapper: S=8: 0→`0011`, 1→`1100`; S=2: identity.
  - TX order: CRC → whitening (7-bit LFSR seeded with the channel index) → FEC → mapper. AA, CI and TERM1 are not
    whitened. The SIG 2019 slides draw whitening after the mapper; the spec text puts it before FEC [U].
- **Remote ID primary:** ADV_EXT_IND (type 0111) on channels 37/38/39 at S=8, AA 0x8E89BED6. F3411-19 example:
  - header 0x700D;
  - extended-header flags 0x19 (AdvA + ADI + AuxPtr), no AdvData;
  - 1552 µs.
- **AuxPtr (LSB first):** channel 6 b | CA 1 b | units 1 b (0 = 30 µs, 1 = 300 µs) | offset 13 b | PHY 3 b (010 = Coded).
  - ADI = DID 12 b + SID 4 b.
  - Extended-header flag bits 0–7: AdvA, TargetA, CTEInfo, ADI, AuxPtr, SyncInfo, TxPower.
  - F3411: aux channel = (ch+9) % 36; example offsets 166 / 114 / 62 × 30 µs.
  - F3411-19 Table 15 prints the ext-header-length octet as 0x30 (MSB-first drawing): 0x0C on air [U].
- **AD data:**
  - Layout: `len, 0x16, FA FF, 0x0D, counter`, then the message pack `0xF2, 0x19, N, N × 25 B`.
  - Up to 10 messages (opendroneid caps at 9: ≤ 255 octets).
  - BT4 legacy adverts are mandatory too.
  - An AUX_ADV_IND with 4 messages is ~121 B: 8.3 ms at S=8.
- **Not verified:** F3411-22a / ASD-STAN BT5 differences (download 404); the F3411-19 Table 17 AUX_ADV_IND example
  (PDF cache full).

## Part 2: DJI DroneID
- **OFDM:** 15.36 MSps, FFT 1024, 15 kHz, 601 used carriers (600 + DC); CP long 80, short 72 samples.
  - ZC: length 601 with the middle element removed; roots 600 then 147.
- **9-symbol frame (Mini 2, Mavic Air 2; S1 SoC):**
  - CP [80, 72, 72, 72, 72, 72, 72, 72, 80]: 9880 samples = 643 µs.
  - ZC in symbols 4 and 6. Symbol 1 descrambles to all zeros (a known reference).
  - Payload in symbols 2, 3, 5, 7, 8, 9.
- **8-symbol "legacy" (Mavic 2, older OcuSync):** the first symbol is dropped.
  - ~572 µs (the paper: 576 µs).
  - ZC in symbols 3 and 5; payload in 1, 2, 4, 6, 7, 8.
  - First CP 80 (RUB code) or 72 (proto17) [U].
- **Coding (both variants):**
  - LTE Gold scrambler: Nc 1600, x2 seed 0x12345678, 7200 bits.
  - LTE turbo rate 1/3: K 1408 (176 B), CRC24A, rate matching E=7200 rv0.
  - Record 91 B: byte 0 = length 0x58, byte 1 = type 0x10, byte 2 = version (2 on Mini 2 / Air 2).
  - Inner CRC16: poly 0x1021 reflected, init 0x3692, over bytes 0..88.
  - Version 1: pitch / roll / yaw, no pilot position or time. Version 2: drops pitch / roll, adds GPS time and pilot
    lat / lon.
  - Message 0x11 = flight purpose.
  - Wi-Fi DroneID: vendor IE OUI 26:37:12, then 58 62 13, then the subcommand.
- **Frequencies and timing:**
  - 2399.5–2459.5 MHz in 15 MHz steps (2474.5 [U]); 5741.5–5831.5 MHz in 15 MHz steps.
  - One burst every ~600–640 ms, 12–20 per frequency before hopping, only while the motors run.
  - O4: encrypted.
- **Telling the variants apart:** burst length 643 vs 572 µs, or symbols before the ZC600: 3 vs 2.
- **Conflicts:**
  - Field order (altitude / height) differs between RUB and proto17.
  - UUID: 20 B (RUB) vs 19 B + pad (proto17).
  - Bender's figures show 2-byte coordinates (Kismet, RUB and proto17 use 4).
  - The 8-symbol model list beyond Mavic 2 is inferred.

## Sources
- Bluetooth Core 5.4 Vol 6 Part B: https://www.bluetooth.com/wp-content/uploads/Files/Specification/HTML/Core-54/out/en/low-energy-controller/link-layer-specification.html
- SIG BT Asia 2019 slides 30–34: https://www.bluetooth.com/wp-content/uploads/2019/05/BTAsia/Bluetooth-Asia-2019-Longer-Range-Higher-speed-and-Increased-Capacity.pdf
- MathWorks packet structure: https://www.mathworks.com/help/bluetooth/ug/bluetooth-packet-structure.html
- USPTO 11057855 (generators): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11057855
- Nordic DevZone 75320: https://devzone.nordicsemi.com/f/nordic-q-a/75320/transmitting-data-via-le-coded-phy
- ASTM F3411-19 §5.4.7: https://thedroneprofessor.com/wp-content/uploads/2022/11/F3411.40165-UAS-Remote-ID.pdf
- opendroneid-core-c, sxjack/remote_id_bt5, opendroneid receiver-android (GitHub)
- NDSS 2023 (§III-B, Fig. 13): https://www.ndss-symposium.org/wp-content/uploads/2023-217-paper.pdf
- RUB-SysSec/DroneSecurity, proto17/dji_droneid, alphafox02/antsdr_dji_droneid (GitHub)
- Bender arXiv 2207.10795; Kismet DJI DroneID KSY
