# marine — AIS ship and station reports (research)

**Validation: research.** Built from ITU-R M.1371-5 via `research/signal_gaps/06_vhf_below_300.md` section 5.
Checked on the research file's golden vector (the gpsd sentence encoded to on-air NRZI bits, and decoded back) and on
synthetic GMSK bursts only. It needs a receiver that tunes to 162 MHz (Pluto with the AD9364 range).

**Detects and decodes** AIS on both channels, AIS 1 (161.975 MHz) and AIS 2 (162.025 MHz), with the receiver ppm
applied. Per channel:
1. ±10 kHz at 48 kSps (5 samples per bit), FM discriminator, the carrier offset removed (32-bit moving average),
   integrated over one bit.
2. For each sample phase: NRZ levels → NRZI decoding → training + flag 0x7E → de-stuffing up to the closing flag →
   octets (LSB first).
3. **CRC-16/X.25**: only frames whose residue is 0xF0B8 are kept. Every reported message is CRC-verified.

**Fields**, per message: `msg_type`, `mmsi`, `mid` / `flag` (419 = India), `kind`, plus:

| Type | Fields |
|---|---|
| 1 / 2 / 3 class A position | `nav_status`, `sog_kn`, `lat`, `lon`, `cog_deg`, `heading_deg`, `utc_second`, `accuracy` |
| 4 / 11 base station | `utc`, `lat`, `lon` |
| 5 static and voyage | `imo`, `callsign`, `name`, `ship_type`, `length_m`, `beam_m`, `draught_m`, `destination` |
| 18 class B position | `sog_kn`, `lat`, `lon`, `cog_deg`, `heading_deg` |
| 21 aid to navigation | `aid_type`, `name`, `lat`, `lon`, `virtual` |
| 24 class B static | part A `name`; part B `ship_type`, `callsign` |

The signal's `vessels` list merges every message per MMSI. Status is `success`. The label is `AIS (n stations)`.

**Sensitivity (synthetic):** about 10–12 dB SNR in a 25 kHz channel, with carrier offsets of ±1.5 kHz.

**Limitations:**
- Hard-decision discriminator demodulation: a coherent / Viterbi GMSK receiver would gain 2–3 dB.
- No single-bit error correction.
- Message 5 spans two slots (51 ms) and needs a block longer than 34 ms.
- The bit layouts of types 18 / 21 / 24 come from memory in the research file and still need checking against the
  standard. Types 1 / 2 / 3 are checked against the gpsd golden vector.
