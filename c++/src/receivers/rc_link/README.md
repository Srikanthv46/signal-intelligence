# rc_link (C++) — LoRa / ExpressLRS, GFSK/FSK RC links, BLE advertising, 2 Mb/s ESB-like links

Port of `python/spectrum_intel/receivers/rc_link` (detector.py, esb.py, lora_phy.py, label_check.py, plugin.py). Same
method; see the Python README for the details.

## What it does
- **LoRa / ExpressLRS:** blind dechirp scan per sub-band, SFD verification, symbol walk, ELRS channel grid; payload
  chain as gr-lora_sdr (Gray, deinterleave, Hamming, header checksum, de-whitening, CRC-16) with the LoRaWAN and
  Meshtastic plain header fields. `success` when the CRC passes.
- **GFSK/FSK RC and telemetry links:** matched-area energy scan, local refinement, symbol rate, labels
  (FrSky-like, ELRS FLRC/GFSK-like).
- **BLE advertising:** access address, de-whitening, CRC-24 → PDU type and advertiser address.
- **ESB-like 2 Mb/s links (`esb.cpp`):** 2 Mb/s slicing behind a 0101 preamble, joint CRC-16/CCITT check of the
  packets sharing a sync word (span and init from the packets alone).
- **Frequency-label check:** `check_frequency_label` reports a wrong capture centre from CRC-valid BLE adverts.

## C++ API (`namespace si::rx::rc_link`, depends on core/ and io/ only)
```cpp
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/label_check.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc);          // LoRa + FSK + BLE, images suppressed
std::optional<si::Json> fix = check_frequency_label(em);     // nullopt = label consistent
```
Building blocks (`lora_scan`, `lora_refine`, `walk_packet`, `fsk_scan`, `refine_local`, `fsk_classify`,
`ble_try_decode`, `lora_phy::*`) are public for tests.

## Implementation notes
- dtypes and operation order as numpy (float32 spectra, `fft::np_*`, `np::cabs` / `np::cmul`); numpy temporary
  elision on ≥ 256 KiB arrays changes complex-multiply operand order, mirrored by `mul_conj(..., elided)`.
- `dsp_util` holds receiver-local helpers (`pymod`, `slice`, `upchirp`, `kmeans2`); candidates for core/.

## Parity / tests
- `tools/parity.py --receivers rc_link`: 16/16 inputs identical (9 rc_link signals, 141 energy signals).
- `tests/unit/test_rx_rc_link.cpp`: LoRa PHY (whitening, CRC, header checksum, symbol count), scan/refine/walk and
  FSK/BLE stages against Python vectors.

## 2026-10-22 market-gap batch (ported with Python, PLAN G2-G6)
`radiolink.cpp`, `hubsan.cpp`, `toys.cpp`, AFHDS 1 in `afhds2a.cpp`, LR1121 GFSK SYNC in `elrs.cpp`, variable XN297
payloads + the toy table in `xn297.cpp`, `sik::demod_bits_mid` / `two_level_threshold`, `lora_phy::mavlink_in`; the
FSK chain (`fsk.cpp`) gates, the split-packet carry and the RadioLink piece absorption, plugin fields and model rules
as in Python. `tests/unit/test_rc_more.cpp` (vectors `tests/vectors/rc_more/gen_rc_more.py`): hop tables, detect labels
and decoded packets, plugin fields and model blocks identical to Python; `test_drone_links.cpp` covers MAVLink over LoRa.

## 2026-10-20 (ported with Python)
- LoRa SF5 at 500 kHz in the scan; `lora_phy::decode_sf5_implicit` (implicit-header SX126x SF5 packets, length and
  CR found by the payload CRC-16, `cr_ambiguous`); the walk's repeated-symbol end rule scales with N (3 / 4 / 5).
- LoRa-only mode on 860–930 MHz for captures outside the RC bands (`detect(..., lora_only)`, `lora_scan` f_range,
  plugin `in_range` / `can_decode`). `tools/parity.py --receivers rc_link --inputs LIVE_LORA_SF5_890MHz,
  LIVE_PMR_860.4MHz,LIVE_ELRS_2G4_2440MHz`: identical (SF5 packet decoded by both).
- SF5 stage-1 threshold `kEtaSf5` = 40; `SliceCache` (per sub-band, `detect_lora` refines sub-band groups in
  parallel). Parity after both (with `real_434MHz_sf9`): 4/4 identical; C++ on `LIVE_PMR_860.4MHz` 13.1 -> 5.7 s.
- `tests/unit/test_rx_rc_link.cpp`: "SF5 implicit-header decode" on the live packet's symbols, an uncorrectable
  block and the ambiguous zero-tail case, against Python's output.

## 2026-10-17 (ported with Python)
- ESB-like 2 Mb/s links (`esb.cpp`, `esb_try` in fsk.cpp) and the `fsk_refine` level of a short packet in a long
  box (98th percentile when the 75th is at the floor). `tools/parity.py --inputs LIVE_ESB2M_2440MHz`: 15/15 signals
  identical, every decoded field (CRC span / init, frame hex, frame period) included.

## 2026-10-10 (ported with Python)
- LoRa sync symbols must be a sync word's nibbles x 8 (+-1).
- A preamble-only LoRa report needs the SFD search window past the end of the capture, or a seen SFD (z_sfd >= 4).

## 2026-10-23 (PLAN G10 / G12, ported with Python; research)
- `hott.{hpp,cpp}`: Graupner HoTT (CC2500 MSK, FEC: `fec_encode` / `fec_decode` hard-decision Viterbi, DN504 vector),
  `hitec.{hpp,cpp}`: Hitec Optima / Minima (sync = TX ID) and Futaba S-FHSS (`crc16_cc2400`, both inits),
  `adsl.{hpp,cpp}`: ADS-L SRD-860 M-band (`crc24`, XXTEA zero key) and the OGN tracker (`ogn::tea_key0`, Manchester-clean
  check), `lora_phy::fanet_parse` (FANET in `identify`). `fsk.cpp`: `FskMeas::cc2500_rc` / `adsl`, the new gates
  (FLARM-PHY pieces >= 40 chips), labels hott / hitec / sfhss / adsl / ogn_like, carry and covered pieces; plugin fields
  `cc2500_rc`, `adsl`, `ogn`, `fanet` and the model blocks. == python on `tests/vectors/g10_g14` (`test_g10_g14.cpp`).

## Signal gaps S7 (research)
- `core/ble_trackers.{hpp,cpp}`: the BLE tracker / device recognizers (python `core/ble_trackers.py`), called from
  `core/remote_id` `parse_ble_ad` (`ad.device`); the plugin adds `ble_devices` and the "Bluetooth LE tracker (...)"
  label.
- `bredr.{hpp,cpp}`: Bluetooth BR/EDR access codes (python `bredr.py`): `syncword(lap)`, `find(bits)`,
  `try_decode(y, ofs)` on 1 Mb/s hops that are not BLE / XN297 / nRF24 → "Bluetooth BR/EDR (LAP …)" and `bredr`
  fields. complex128 filter and discriminator as numpy (temporary-elision operand order). Vectors:
  `tests/vectors/signal_gaps/bredr*`, `ble_trackers*`.
