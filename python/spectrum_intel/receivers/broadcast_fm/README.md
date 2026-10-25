# broadcast_fm — FM broadcast stations and RDS (research)

**Validation: research.** Built from IEC 62106 / EN 50067 via `research/signal_gaps/06_vhf_below_300.md` section 8.
Checked on the research file's worked RDS group and on synthetic stations only. It needs a receiver that tunes to
88–108 MHz (Pluto with the AD9364 range).

**Detects** FM broadcast stations on the 100 kHz raster (87.6–107.9 MHz). A channel counts when its power over
±100 kHz is ≥ 10 dB over the band floor and it is a local maximum on the raster.

**Measures:** the station's channel and carrier, and its peak deviation (FM broadcast: 15–110 kHz). The station is
`stereo` when the 19 kHz pilot line of the demodulated MPX stands ≥ 10 dB over its surroundings.

**Decodes RDS** when the block holds ≥ 45 ms (two RDS blocks):
- 57 kHz (3 × the measured pilot) down to baseband, BPSK phase from the squared signal;
- biphase matched filter at 1187.5 bit/s, with the bit timing searched;
- differential decoding, then block sync on the offset words A / B / C / C′ / D;
- **only error-free (26,16) blocks** are used (no burst-error correction).

Fields (`rds`): `pi` (and its country nibble: India is 5), `pty` / `pty_name`, `tp`, `ta`, `music`, `ps` (the
8-character programme service name, `_` for characters not yet seen), `radiotext`, `clock_time` (UTC and the local
offset), `group_types`, `blocks_ok`.

**Status:** `success` when RDS blocks are checkword-verified, else `partial`. Many Indian stations send no RDS, so a
station without it is still reported.

**Label:** `FM broadcast (<MHz>[, stereo][, RDS '<name>' | RDS PI <pi>])`.

**Limitations:**
- One 34 ms block at 61.44 MSps is too short for RDS. A lower sample rate or a longer dwell (`StreamAnalyzer`)
  gives PI and the programme name; RadioText needs seconds.
- No burst-error correction: weak stations give fewer blocks.
- The ECC (extended country code, group 1A) is not decoded yet.
