# pmr_digital — digital land-mobile radio (TETRA, DMR, NXDN)

**Detects** digital land-mobile carriers and **decodes** what the air interface carries in the clear:

| System | Detected from | Decoded | Status |
|---|---|---|---|
| TETRA downlink | 18 ksym/s line of the matched-filtered envelope; training sequences on the 510-bit slot grid | **BSCH of every synchronisation burst (CRC-CCITT): MCC, MNC, colour code**, system code, TN/FN/MN, sharing mode, cell service level | `success` with a CRC-valid BSCH; `partial` with training sequences only; "TETRA-like carrier" from the symbol-rate line alone |
| DMR (Tier II/III) | 48-bit SYNC words (BS/MS/direct, voice/data) | slot type Golay(20,8): **colour code, data type**; BPTC(196,96) Hamming check; CSBK opcode / FID when the CRC passes; base station vs mobile vs direct mode; voice activity | `success` with a CRC-valid CSBK, otherwise `partial` |
| NXDN (4800 or 2400 sym/s) | 20-bit frame sync word (0xCDF59) | — (LICH/SACCH not decoded) | `partial` |

## How
1. **Candidates.** The capture spectrum (the pipeline's Channelizer FFT) in 1 kHz bins; peaks ≥ 8 dB over the
   floor, ≤ 30 kHz wide, in 380–400 / 410–430 / 450–470 / 870–876 / 915–921 MHz (TETRA) and 136–174 / 300–520 /
   806–870 MHz (4FSK). The receiver's own clock lines (pipeline `receiver_lines`) are skipped.
2. **TETRA** (ETSI EN 300 392-2). The 25 kHz channel at 144 kSps (8 samples/symbol). The 18 ksym/s line must be
   ≥ 6× its neighbourhood to demodulate (≥ 9× for an unframed "-like" report). π/4-DQPSK: RRC 0.35, the carrier
   rotation from the 4th power of the 1-symbol differential, re-estimated per 2000 symbols with the timing phase
   (tracks a sample-clock error over long captures). The 4th power is ambiguous by a quarter turn per symbol
   (4.5 kHz): all three hypotheses are demodulated and ranked by CRC-valid BSCH, then slot-consistent sync
   training sequences (38-bit, ≤ 4 errors), then normal training sequences (22-bit, ≤ 2 errors, only when most of
   them share one slot phase). BSCH: the 120 bits before each sync training sequence → descramble (e = 0) →
   de-interleave (K = 120, a = 11) → de-puncture RCPC 2/3 → soft Viterbi of the 16-state rate-1/4 mother code →
   CRC-CCITT → SYNC PDU + D-MLE-SYNC.
3. **4FSK** (DMR: ETSI TS 102 361-1; NXDN). The 12.5 kHz channel at 48 kSps, constant envelope (cv ≤ 0.35), FM
   discriminator. DMR: the 9 SYNC patterns on a coarse slice (≤ 6 errors); per burst, RRC 0.2 matched filter,
   timing from the sync itself (1/20-symbol grid), a sync correlation ≥ 0.9 and a level fit within 0.5–1.5× the
   DMR deviation (648 Hz per level) — the fit would otherwise scale anything into the pattern — then the burst is
   sliced: SYNC ≤ 4 errors, slot type Golay(20,8) (≤ 3 errors), BPTC(196,96) (×181 de-interleave, Hamming(13,9)
   columns and (15,11) rows with single-error correction), CRC-CCITT with the data-type mask (B.3.12). NXDN: the
   frame sync word at 4800 then 2400 sym/s; one exact word, or several on the 192-symbol frame grid.
4. **Arbitration.** A carrier verified here is digital: another receiver's analog reading of it (voice_fm) is
   dropped.

The symbol-rate line is reported for 4FSK but not required: in a 34 ms block a DMR burst at 11 dB gives only
×2.6, while its sync word is decisive.

## Validation
- **TETRA**, live Pluto captures (2026-10-10, `testdata/live`): `LIVE_TETRA_390.78MHz` (1.75 s, 2.4 MSps): two
  carriers, 390.0288 and 391.5289 MHz, BSCH CRC-valid 19/21 and 52/57: MCC 404 (India), MNC 4040, colour code 7,
  system code 2, frame/multiframe counters consistent. `LIVE_400MHz` (one 34 ms block at 61.44 MSps): both
  carriers decoded (2 BSCH each). No external oracle (osmo-tetra is AGPL; implemented from the standard; the CRC
  and the consistent counters are the check).
- **DMR**, `LIVE_PMR_860.4MHz` (1.75 s): repeaters at 859.8209 MHz (voice call: voice LC header, voice sync,
  terminators) and 860.8210 MHz (control: CSBK in every burst), colour code 0, sync correlation ≥ 0.988, level
  622–652 Hz. Slot types are exact Golay codewords; BPTC rows/columns pass Hamming. **Their CSBKs fail the
  standard CRC with the CSBK mask** (the CRC routine reproduces the reference vector): a vendor / non-standard
  control channel, so no opcode is claimed. `LIVE_868MHz` (34 ms): 860.82 MHz found from one burst.
- **NXDN**, `LIVE_PMR_860.4MHz`: 859.3455 MHz, 43 frame syncs, all on the 384-bit frame grid.
- **Null tests**: the PMR446 capture's clipping products (analog FM) and the drone / Wi-Fi / cellular testdata give
  no DMR/NXDN/TETRA (a clipping product passed the sync search before the correlation and deviation checks).

## Limitations
- A 34 ms block holds ~2.4 TETRA timeslots: a BSCH is decoded when a synchronisation burst falls in the block
  (the main carrier sends one at least every multiframe, 1.02 s; this network sends them in many slots).
- NXDN: a 40 ms frame (80 ms at 2400 sym/s) does not always fit a 34 ms block; without its sync word the carrier
  stays with voice_fm / the band rules.
- TETRA MAC/MLE broadcast (SYSINFO: main carrier, location area), DMR voice (AMBE+2) and NXDN LICH/SACCH are not
  decoded; TETRA traffic is TEA-encrypted in practice.

## API (depends on numpy, scipy and `spectrum_intel.core` only)
```python
from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.pmr_digital import detect, tetra, fsk4
emissions = detect(Channelizer(x, fs, fc))                 # core.emission.Emission per carrier
bits, soft, err, f_off = tetra.demod(y)                    # y: carrier at 144 kSps
ok, sync_pdu = tetra.decode_bsch(soft[h - 120:h])          # h: a sync training sequence position
bursts = fsk4.dmr(fsk4.discriminator(y48k, 48e3), 48e3)
```

## Pipeline plugin
`plugin.py` (`PmrDigitalDecoder`, an `EmissionReceiver`): one capture-wide detection, one signal per carrier,
`owns_band` (its energy detections are absorbed), and `arbitrate` (drops voice_fm's reading of a verified digital
carrier). Technologies: `tetra`, `dmr`, `nxdn` (schema/technologies.json); identifiers `mcc`, `mnc`,
`colour_code`.
