# Signal gaps (batch S): the tracking file

**User, 2026-10-25:** "lets do priority 1, 3 and Out of range (below 300 MHz) but high SIGINT value ... with research
tag ... research, think hard and do it". Origin: the signal-intelligence gap review of 2026-10-25 (the library covers
every signal, not drones only).

## Rules for this batch
- Built from published specifications (GPL / AGPL code only as a reference or test oracle, never copied; C++ stays
  GPL-free). Validated on **synthetic** signals generated from the specs, plus spec test vectors where they exist.
- Every new capability is tagged `research` in `schema/capabilities.json` (`sources: research/signal_gaps/<file>`),
  gets a `schema/technologies.json` rule and a `schema/descriptions.json` text. Labels are "-like" unless a CRC /
  parity / FEC check verifies them. Promotion to `proven` needs a real capture.
- Python first, then C++ with vectors (`c++/tests/vectors/<rx>/`), identical output.
- README, receiver README and `c++/README.md` updated per item.
- One regression + parity at the end of the batch (ask first, naming the data paths), 2 workers, `nice -n 19`,
  under `caffeinate -i`.

## Research files
| File | Topic |
|---|---|
| `01_satellite.md` | Iridium (burst, UW, IRA / IBC), Thuraya GMR-1, Inmarsat (IsatPhone, C, Aero), Globalstar / SPOT |
| `02_gnss.md` | GPS L1 C/A + L5, NavIC L5 / S, GLONASS L1OF, BeiDou B1I, Galileo E1; jamming / spoofing indicators |
| `03_cellular_uplink_fbs.md` | not written: the research agent's write was stopped by a safety classifier (2026-10-25) |
| `04_short_range.md` | 315 / 433 MHz OOK / FSK (EV1527, PT2262, KeeLoq, TPMS), Bluetooth BR/EDR access code, BLE trackers |
| `05_aviation_radar_jammers.md` | radio altimeter, ILS glide slope, Mode A/C, radar bands + pulse deinterleaving, jammer types, AM voice |
| `06_vhf_below_300.md` | 70–300 MHz: airband AM, ACARS, VDL Mode 2, AIS, marine VHF + DSC, APRS, FM broadcast + RDS, VOR / localizer |

## Items
Status: "todo" → "py" (Python + tests) → "batch" (C++ + vectors, waiting for the batch regression) → "regressed".

| # | Item | Group | Status |
|---|---|---|---|
| S1 | Iridium burst decode (UW, IRA / IBC fields) + uplink handset bursts | P1 | **blocked**: the decoder write was stopped by a safety classifier (decodes paging / TMSI / handset access of a comms system); user to decide |
| S2 | Thuraya / Inmarsat / Globalstar-SPOT bursts: "-like" by measured PHY | P1 | on hold (with S1; the PHY-only "-like" detect part is safe to revisit if the user wants) |
| S3 | GNSS acquisition (GPS L1 C/A), C/N0, jamming / spoofing indicators (`gnss`; BeiDou / GLONASS / NavIC not yet) | P1 | batch (Py + C++ `vectors/signal_gaps/gnss`) |
| S4 | phone uplinks: LTE SC-FDMA (7.5 kHz shift), NR uplink, GSM uplink bursts | P1 | on hold: the research write-up (`03`) was stopped by a safety classifier; user to decide |
| S5 | fake base station indicators (report level) | P1 | on hold (as S4) |
| S6 | 315 / 433 MHz OOK / FSK devices (`ism_devices`: EV1527 / PT2262, KeeLoq clear fields, TPMS VDO + Toyota) | P1 | batch (Py + C++ `vectors/signal_gaps/ism_devices`); more TPMS families / weather sensors still todo |
| S7 | BLE trackers + notable devices (`core/ble_trackers`: Find My / AirTag, SmartTag, FMDN, Tile, Chipolo, DULT; iBeacon, Eddystone, AirPods, Swift / Fast Pair); Bluetooth BR/EDR access code → LAP (`rc_link/bredr`) | P1 | batch (Py + C++ `vectors/signal_gaps/ble_trackers`, `bredr`); BR/EDR seen on a real capture (Phantom 4 testdata: LAP DA7476, B4BDD9, 0 bit errors); UAP / EDR not done |
| S8 | radio altimeter (FMCW 4.2–4.4 GHz) (`radalt`) | P3 | batch (Py + C++ `vectors/signal_gaps/radalt`) |
| S9 | ILS glide slope (329–335 MHz, 90 / 150 Hz DDM) (`navaids`) | P3 | batch (Py + C++ `vectors/signal_gaps/navaids`) |
| S10 | named radars (`classify/rules._radar_class`: S-band ATC / solid-state ATC / weather S + C / marine / C-band / L-band surveillance / wind profiler) and pulse deinterleaving (`core/deinterleave`: PW groups, CDIF, sequence search with a chance test, stagger, jitter mode → `pulse_emitters`, pulsed-radar evidence) | P3 | batch (Py + C++ `vectors/signal_gaps/deinterleave`) |
| S11 | Mode A/C (and IFF Mode 3/A) replies at 1090 MHz (`adsb/modeac`) | P3 | batch (Py + C++ `vectors/signal_gaps/modeac`) |
| S12 | jammers: GNSS swept jammers (PPD, `gnss.swept_ppd`), GNSS CW / floor rise (S3), barrage-noise and noise-FM jammers (`jammer/barrage`, plugin `jammer_barrage`) | P3 | batch (Py + C++ `vectors/signal_gaps/barrage`); spot jammers not done (CW false-positive risk: receiver lines, LO leakage) |
| S13 | AM voice (UHF military air 225–400, VHF civil air) (`voice_fm` AM path) | P3 / VHF | batch (Py + C++ `vectors/signal_gaps/am`) |
| S14 | 70–300 MHz range: band plan (`core/bands_india`, Py + C++); the input never had a range limit | VHF | batch |
| S15 | AIS (`marine`: types 1/2/3, 4, 5, 18, 21, 24) | VHF | batch (Py + C++ `vectors/signal_gaps/ais`) |
| S16 | ACARS | VHF | **blocked**: the work was stopped by a safety classifier (2026-10-25); not pursued |
| S17 | VDL Mode 2 | VHF | on hold (same kind of message-content decoder as S16); user to decide |
| S18 | marine VHF voice + DSC | VHF | marine FM voice is covered by `voice_fm` (118-174 MHz); DSC on hold with S16 / S17 |
| S19 | APRS / AX.25 (`core/hdlc` + `voice_fm/afsk`: Bell 202 AFSK on any FM carrier) | VHF | batch (Py + C++ `vectors/signal_gaps/aprs`) |
| S20 | FM broadcast + RDS (`broadcast_fm`) | VHF | batch (Py + C++ `vectors/signal_gaps/broadcast_fm`) |
| S21 | VOR + ILS localizer (`navaids`) | VHF | batch (with S9) |

## Pending for the end-of-batch run (ask first, with these paths)
- Regenerate `c++/tests/vectors/rc_link/ref.json` (`gen_rc_link.py`): it reads `testdata/lora/real_434MHz_sf9.cf32` and
  `testdata/DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta`; real adverts there now carry the S7
  `device` key, and 2–3 hops there now decode as BR/EDR (C++ test_rx_rc_link fails until then: expected).
- Expected regression differences: `pulse_emitters` / `pulse_unassigned` features and a "deinterleaved into ..."
  evidence line on pulsed-radar signals; BR/EDR labels on 2.4 GHz hops; `ble_devices` on BLE signals; possible
  `Jammer (barrage noise)-like` / `noise-FM` labels on wide continuous boxes (check each); new receivers
  (`ism_devices`, `radalt`, `navaids`, `gnss`, `broadcast_fm`, `marine`) on captures in their bands.
- Live testdata for the new signals (`testdata/live`, CLAUDE.md rule): none captured yet (433 MHz remotes / TPMS,
  BR/EDR, radio altimeter, GNSS L1, and the VHF items need the bench).
- Regression (`bench/regress.py`) on all testdata and parity (`c++/tools/parity.py --refresh --threads 2`).
