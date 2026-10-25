# Drone transmitters: the single tracking file

**This is the only plan to follow (user, 2026-10-22).** Read it first after a context reset. Every item lives here
with its status; update this file whenever an item changes state. Research details live in the numbered files of this
folder; this file only tracks them.

## Goal (user, 2026-10-21 / 22)
- Detect **any drone transmitter**: analog video, digital video, RC links (controller and receiver), telemetry, Remote
  ID / DroneID. Decode where possible and extract the **maximum information**: system, generation, **model**,
  channel, mode, antenna count, unit fingerprint (ppm), hop pattern, IDs and decoded fields.
- **Model level, for every family** (user, 2026-10-22): within a family, identify the model (e.g. Walksnail Avatar
  Pro / Mini 1S / Moonlight / GT; DJI O3 / O4 / O4 Pro; ELRS hardware; Crossfire Micro / Nano), or say "model
  unresolved" with the candidate set; never guess.
- **Low SNR** (synthetic SNR sweep per capability, threshold in the receiver README) and **milliseconds** (time every
  capability, Python and C++).

## Rules for this work
- Build from **research** (published specs; GPL code only as a reference / oracle, never copied; the C++ stays GPL-free).
  Validate on **synthetic** signals made from the specs.
- Tag every new capability `research` in `schema/capabilities.json` (technology or field capability, with
  `sources: research/drone_transmitters/<file>`) and add a `schema/technologies.json` rule for each new label.
  Labels are "-like" unless a CRC verifies them. Promotion to `proven` needs a real-transmitter capture (add `proof`).
- Python first, then the C++ port with unit-test vectors (`c++/tests/vectors/<rx>/`), identical output.
- Update `README.md`, the receiver README and `c++/README.md` for every item.
- **Batch workflow** (user, 2026-10-22): implement every item of the batch with its own unit tests and vectors; run
  the full regression (`bench/regress.py`) and parity (`c++/tools/parity.py --refresh --threads 2`) **once** at the
  end of the batch, under `nice -n 19`, 2 workers at most. Testdata may be used for regression and parity only.

## Research files (this folder)
| File | Topic |
|---|---|
| `01_analog_video.md` | NTSC / PAL / PAL-M / PAL-N timings, VTX RF, channel tables |
| `02_digital_video.md` | DJI, Walksnail, HDZero, Shark Byte, Autel, Siyi, Skydroid, COFDM, Wi-Fi FPV |
| `03_rc_links.md` | ELRS, Crossfire, Tracer, Ghost, FrSky, FlySky, DSM, Futaba, Hitec (+ `src_gpl_reference/`, oracle only) |
| `04_telemetry_remoteid.md` | SiK / RFD, Microhard, Doodle Labs, XBee, Remote ID, DroneID, ADS-B, FLARM (+ `src_sik/`, BSD-2) |
| `05_sik_mavlink_details.md` | Si443x packet, MAVLink framing / CRC_EXTRA |
| `06_dvbt_cofdm.md` | DVB-T, DVB-T2 P1, drone COFDM products (+ `src_etsi/`) |
| `07_elrs_details.md` | ELRS FLRC, OTA packets, SYNC, CRCs, IQ polarity |
| `08_crossfire_dsm.md` | Crossfire 150 Hz FSK, CYRF6936 DSM2 / DSMX |
| `09_ble_coded_droneid_legacy.md` | BLE LE Coded PHY, DroneID 8 vs 9 symbols |
| `10_cc2500_a7105_xn297_flarm_xbee.md` | FrSky D8 / D16, AFHDS 2A, XN297 / Bayang, FLARM, 802.15.4 / XBee |
| `12_rc_links_more.md` | mLRS, FrSky R9, S-FHSS / Hitec / HoTT / RadioLink, OpenLRS-ng, AFHDS 3 / FASST / Tracer / Ghost (unpublished), SX127x SF6 framing |
| `13_model_gaps_droneid.md` | FAA DOC manufacturer codes (~140 more) and serial prefixes -> model (~55 DJI models + other makers), DJI AeroScope product-type map (codes to 173), 8 vs 9 DroneID symbols by model, Walksnail header (GPL-only) |
| `14_market_gap_survey.md` | market gap survey 2026-10-22: VTXs, RC links, telemetry not yet covered, with priorities |
| `16_uat_adsl_fanet_ogn.md` | UAT 978 (DO-282B, RS constants), ADS-L SRD-860 (EASA Issue 2), FANET, OGN tracker (LDPC matrix GPL / unlicensed only), China GB 46750-2025 RID |
| `18_wifi_field_gaps.md` | Wi-Fi field gaps (batch W): what existed, drone-by-behaviour literature, hotspot names, Tello SDK, MAC randomisation |
| `17_hott_hitec_r9.md` | HoTT (CC2500 MSK + FEC, 16 hop tables), Hitec Optima / Minima, Futaba S-FHSS, FrSky R9, CC2500 FEC (DN504) |
| `11_model_level_identification.md` | per-family model discriminators: decision tables, models that look identical on air, cited sources |

## 1. Added (in the library)
Status "regressed" = covered by a full regression + parity run; "batch" = unit tests + vectors pass, waiting for the
end-of-batch regression.

| # | Capability | Py | C++ | Status |
|---|---|---|---|---|
| V1 | validation tag (`capabilities.json`, API `validation`) | yes | yes | regressed |
| A1–A4 | analog PAL-M / PAL-N synth, line-locked chroma standard, Lowband channel table | yes | yes | regressed |
| A5 | analog audio subcarrier FM demod: level, tone (real NTSC VTX: 239 Hz) | yes | yes | regressed |
| A6 | sync-phase late lock | closed, not reproduced (54/54 within 0.1 µs) | – | closed |
| A7 | analog SECAM (Db / Dr line alternation) | yes | yes | regressed |
| D1 | digital_vtx system catalog → `system_candidates` | yes | yes | regressed |
| D4 | DroneID 8 vs 9 symbols (`droneid_frame`), carrier fix | yes | yes | regressed |
| D5 | COFDM DVB-T (TPS) + DVB-T2 P1 (S1 / S2, margin gate) | yes | yes | regressed |
| D7 | Walksnail frame period (fold R²; real 7.162 ms) | yes | yes | regressed |
| R1 | RC protocol catalog → `protocol_candidates` | yes | yes | regressed |
| R2a | ELRS 2.4 GHz UID[5] bit 0 from IQ polarity (`elrs_identity`) | yes | yes | regressed |
| R2c | ELRS sub-GHz SYNC decode + self-verify (`elrs_sync`), implicit LoRa vs gr-lora_sdr | yes | yes | regressed |
| R6a | **FrSky D8 / D16 v1 / v2 / LBT**: sync + length + CRC → variant, TX ID (`frsky`) | yes | yes (`vectors/frsky`) | regressed |
| R7 | **FlySky AFHDS 2A**: ID code 0x5475C52A → "-like"; FEC + CRC under candidate conventions → type, TX / RX IDs (`afhds2a`) | yes | yes (`vectors/afhds2a`) | regressed |
| M0 | **model level, common output** `fields.model` {family, model, candidates, basis, evidence, unresolved_reason} (`core/drone_models`) | yes | yes (`vectors/drone_models`) | regressed |
| M1 | DJI DroneID product-type byte → model (codes 1–70; unknown stays unknown); frame layout → generation only | yes | yes (C++ has no DroneID payload decoder: generation path only) | regressed |
| M2 | Remote ID serial → CTA-2063-A manufacturer (`remote_id.manufacturer`, model block; BLE + Wi-Fi beacon / NAN / DSSS) | yes | yes (remote_id / ble / wifi_nan vectors regenerated) | regressed |
| M3–M4 | digital video model level: catalog systems → model lists, DJI O3 / O4 channel centres, `tx_antennas` (Walksnail 1 vs 2 antenna kits, O4 Lite), HDZero 17 MHz (`catalog.model_candidates`) | yes | yes (`vectors/digital_vtx` catalog cases) | regressed |
| M5 | RC / telemetry model level: ELRS chip class + domain, FrSky variant, AFHDS / Crossfire candidates, SiK RFD 3.x rate set (`rc_link/models.py`) | yes | yes (`vectors/rc_models`) | regressed |
| M8 | analog VTX model level: channel, audio subcarrier, pit-mode 5584 MHz, 1.2 / 1.3 GHz 2-channel units; otherwise unresolved | yes | yes (parity on analog testdata at batch end) | regressed |
| R8 | **XN297 / Bayang**: prefix 0x710F55, LFSR scrambler (matches the published table), Bayang checksum → address, TX ID, hop channels, sticks; CRC residual (`xn297`) | yes | yes (`vectors/xn297`) | regressed |
| T4 | **IEEE 802.15.4 / XBee**: non-coherent O-QPSK template demod, SFD, PHR, FCS, MAC header, Digi OUI → vendor + model block (`ieee802154`) | yes | yes (`vectors/ieee802154` incl. pipeline) | regressed |
| R2d | **ELRS hop recovery**: LCG block-shuffled sequence (matches the reference transliteration), UID2 / UID3 search, domain disambiguation, next channels (`elrs_sync.hop_recovery`) | yes | yes (`vectors/elrs_fhss` incl. pipeline) | regressed |
| M6 | DJI uplink hop timing → model candidates (`catalog.hop_model`; RFUAV table, low confidence; the hop-width definition is unconfirmed on a real Mini 4 Pro / Avata 2 uplink) | yes | yes (`vectors/digital_vtx` hop_model cases) | regressed |
| D3 | hop signatures → Herelink / Skydroid / Siyi / DJI OcuSync 2 / O4-range candidates (generalised `catalog.hop_model`) | yes | yes (`vectors/digital_vtx` hop_model cases) | regressed |
| L1c | XN297 full-bit soft demod + prefix-set decision level (~2.5 dB; 30/40 vs 10/40 at 5 dB in 8 MHz) | yes | yes (`vectors/xn297` low-SNR standalone cases) | regressed |
| T5 | **FLARM** detection + clear header: Manchester, sync F5 31 FA B6, radio ID + address type (only with a violation-free header), payload not decoded (`flarm`) | yes | yes (`vectors/flarm`) | regressed |
| M7a–d | FAA DOC manufacturer codes (142 in all), serial prefix → model (`serial_model`), DroneID product codes to 173 (DJI map, unconfirmed > 70), 8-symbol frame → Mavic Pro / Mavic 2 candidates | yes | yes (`vectors/drone_models`, regenerated remote_id / wifi_nan / ble / droneid vectors) | regressed |
| R2e | ELRS SF6 (200 Hz) SYNC: `decode_implicit` from SF6 (SX127x framing), gr-lora_sdr SF6 vectors byte-exact | yes | yes (`vectors/elrs_sync` SF6 cases) | regressed |
| R3 / T6 | **mLRS**: air-mode identification (kept apart from ELRS), SX127x SF6 frame decode + CRC → FrameSyncWord (`mlrs`) | yes | yes (`vectors/mlrs`) | regressed |
| R9 | TBS Crossfire 150 Hz: sync → packet bytes, grid channel / link | yes | yes | regressed |
| F1 | blind FSK frame: preamble length, sync candidate, payload bits (`features.frame`), known sync names, field `fsk_frame` | yes | yes (`vectors/fsk_frame`) | regressed |
| D8 | **HDZero** 27 MHz mode (measured 2026-10-22 on a user-confirmed HDZero Freestyle V2, `LIVE_HDZERO_5658MHz`): 1024 FFT at 27 MSps, CP 32 (1/32), 889 sc, 26 continual pilots, scattered pilots every 32 carriers shifting 8 per symbol, ~17-symbol preamble once per video frame (426 symbols = 16.66 ms, 60 fps), 1 TX, continuous. Detection on the HDZero channels by continual-pilot coherence vs neighbouring data carriers (down to about −6 dB link SNR on 32 ms; frame rate to −3 dB; 1 ms blocks work) → channel, carrier / ppm, frame rate, model block (every HDZero VTX: unresolved). Catalog + numerology table updated (`hdzero`) | yes | yes (`vectors/digital_vtx` hdzero cases) | regressed |
| G1 | **wfb-ng / OpenIPC / OpenHD** in the 802.11 parser: TA = BSSID = 57:42 + link_id + radio port, packet type, FEC block / fragment; default link domain named (`mac80211`) | yes | yes (`vectors/drone_links`) | regressed |
| G2 | **MAVLink over transparent LoRa** (`lora_phy.mavlink_in`) and **open Wi-Fi UDP / TCP** (`mac80211`), CRC_EXTRA-verified | yes | yes (`vectors/drone_links`) | regressed |
| G3 | **ExpressLRS LR1121 GFSK**: self-verifying SYNC packet whose UID bytes = the sync word (`elrs.gfsk_try_decode`); Gemini Xrossband leg pairing not done | yes | yes (`vectors/rc_more`) | regressed |
| G4 | **RadioLink** CC2500 (sync 0xABCD / 0xBADC / 0xABA5, PN9, CRC + checksum, hop-table check; `radiolink`) | yes | yes (`vectors/rc_more`) | regressed |
| G5 | **FlySky AFHDS 1** (in `afhds2a`) and **Hubsan** (`hubsan`) | yes | yes (`vectors/rc_more`) | regressed |
| G6 | **toy drones**: XN297 MJXQ / CX10 / JJRC345 / H8_3D / E016H / Potensic, nRF24 V2x2 / SymaX / X5C (`toys`) | yes | yes (`vectors/rc_more`) | regressed |
| G7 | **ESP-NOW + ELRS Backpack** (bind UID, MSP command, VTX band / channel; `mac80211`) | yes | yes (`vectors/drone_links`) | regressed |
| G8 | **wider blind-OFDM bands** 1.35–1.444 / 1.6–2.53 / 4.4–6.0 GHz, `band_context` + label suffix (1.4 GHz UAV band, L-band, government / professional bands); cellular band gate in the blind scan (LIVE_1840MHz digital_vtx 3.8 → 0.2 s Python, ~2.6 → ~0.65 s C++) | yes | yes (`vectors/digital_vtx`) | regressed |
| G9 | **3.3 GHz analog** range 3150–3490 MHz + channel groups 3G3A / 3G3B; the wide 4.9–6.0 GHz 80 / 96-channel tables are image-only: detected, no channel names | yes | yes (unit test) | regressed |
| G10 | **ADS-L SRD-860 M-band** (CRC-24, XXTEA zero key; UAS categories), **OGN tracker** (-like: no LDPC; Manchester-clean + parity + plausibility), **FANET** (LoRa 0xF1; India 866.2 MHz); O-band not done | yes | yes (`vectors/g10_g14`) | regressed |
| G12 | **Graupner HoTT** (CC2500 FEC Viterbi, DN504 vector), **Hitec** Optima / Minima (sync = TX ID), **Futaba S-FHSS** (CC2400 CRC, model hint); FrSky R9 not done (R4) | yes | yes (`vectors/g10_g14`) | regressed |
| G13 | **off-grid 802.11 sensitivity**: windowed (250 µs) CP trigger, gated on a whole-capture ratio ≥ 2; Mimosa in `LIVE_MIMOSA_4880MHz_weak` (new testdata) now decoded (21 PPDUs, 14 FCS-ok) | yes | yes (unit test; parity on the new input) | regressed |
| G14 | **UAT 978** (own Reed-Solomon `core/rs`, dump978 sample messages; emitter category 14 = UAV) and **China GB 46750-2025** RID in the F3411 element (DJI Mini 5 Pro frame) | yes | yes (`vectors/g10_g14`) | regressed |
| T1 | SiK telemetry + MAVLink | yes | yes | regressed |
| T2 | Remote ID over BLE legacy adverts | yes | yes | regressed |
| T2b | Remote ID over BLE Long Range (LE Coded S=8 / S=2, extended adverts) | yes | yes | regressed |
| T3 | Remote ID over Wi-Fi NAN | yes | yes | regressed |
| T7 | ADS-B emitter category B6 = UAV | yes | yes | regressed |
| T8 | DroneID payload state (decoded / opaque) | yes | no (no C++ DroneID payload decoder) | partial |
| L1 | coherent GFSK demod (BLE, Chase CRC): Remote ID 12 → 8 dB | yes | yes | regressed |

## 2. Next to add (researched; ready to implement), in this order
**Batch W: Wi-Fi field gaps + drone-vs-regular Wi-Fi (user 2026-10-25, research file `18_wifi_field_gaps.md`).**

| # | Item | Notes | Status |
|---|---|---|---|
| W1 | drone Wi-Fi **link profile** per BSS (`wifi/field.py`, both PHYs): AP, clients, PPDUs / airtime / bytes per direction, downlink gaps, video payloads; `drone_links` verdict from identity or behaviour (0.55 "-like", caveat), label suffix; `devices` finding `drone_link_behaviour`. Cross-block accumulation (video frame rate over ~1 s) not done: needs a real Wi-Fi drone capture first | Alipour-Fanid (TIFS), Alsoliman (ICDCN 2021) | regressed (Py + C++ `vectors/wifi_field`, `vectors/devices`) |
| W2 | **HE-SIG-A** (802.11ax): SU / MU / TB (ER SU not decoded), BSS colour, UL / DL, MCS, NSS, bandwidth, GI, TXOP, CRC-4; colour → BSSID (beacon HE Operation). **Real: 44 / 44 HE PPDUs of LIVE_WIFI_5755MHz pass** (proven); synthetic 100 % from 8 dB | 802.11ax-2021 27.3.11.7 | regressed (Py + C++) |
| W3 | drone / DIY hotspot SSIDs from vendor documents (RMTT-, HolyStoneFPV-, WiFi-720P-, ExpressLRS TX/RX, ArduPilot / PixRacer, CUAVWLINK, DroneBridge), "-like"; model level for the user's 12 Wi-Fi drones (Mavic Air, Mavic Mini, Spark RC, Phantom 3 Standard / 4K, Anafi, Yuneec CGO3, Ruko F11; file 18 §12 drones); GCS hotspots (Herelink / SIYI / Skydroid) have no fixed default name, cameras are not drones | file 18 | regressed (Py + C++) |
| W4 | open-link payloads: Tello SDK command / state (UDP 8889 / 8890; credentials withheld) → drone flag; video datagrams by NAL units (H.264 / H.265 Annex B, RTP dynamic PT + H.264) | file 18 | regressed (Py + C++) |
| W5 | single-stream 40 MHz HT and VHT 20 / 40 MHz payloads (fit in 56 MHz) | 802.11-2020 cl. 19 / 21 | **deferred**: no such PPDU in any testdata (every real VHT is 80 MHz or 2–3 streams); take up with a capture that has them |
| W6 | alerts: deauth flood and duplicate SSID / evil twin **already existed** (`receivers/devices`); added beacon flood (≥ 15 unregistered beaconing BSSIDs, ≥ 8 prefixes) | – | regressed (Py + C++) |
| W7 | MAC randomisation: probe signature + CFO linking **already existed** (`devices.link_randomized`); added sequence-counter continuity (`sequence_continues`, confidence 0.9), `seq_first` / `seq_last` per device | Vanhoef 2016 | regressed (Py + C++) |
| W8 | ESP-NOW Long Range, 802.11ah HaLow (865–868 MHz India) | research first | deferred (no PHY spec for ESP-LR; no captures) |

| # | Item | Source | Notes |
|---|---|---|---|
| G11 | CYRF6936 DSSS demod: DSM2 / DSMX / DSMR, DEVO, M-LINK | files 03 / 08 / 14 | **not pursued** (2026-10-23): follows R5 (DSM: passive timing detection only, already in the `protocols` catalog); the research write-up was stopped by a safety check |
| R4 | FrSky R9: not pursued | – | dropped |
| R6b | not pursued | – | dropped |
| R10 | not pursued | – | dropped |
| R11 | not pursued | – | dropped |
| R5 | DSM2 / DSMX: deferred (passive burst-timing detection only if taken up later; no identity work) | file 08 | deferred |
| R2f | ELRS cross-block tracking: not pursued | – | dropped |
| L1b | combining repeated Remote ID adverts: deferred. In one 34 ms block an aircraft sends ~1 advert per advertising channel, the 3 channels whiten differently, message types rotate and Location content changes, so identical adverts to combine are rare; worth it only across blocks (StreamAnalyzer) for static messages (Basic ID, Operator ID) | – | deferred |
| P1 | whole-pipeline latency to ms. **Measured (C++, 2 threads, one 2^21 block):** 10–16 s per block; front-end detection 0.15 s; receivers: digital_vtx 6–13 s, rc_link 4.8 s (2.4 GHz), dji ~3 s, wifi_ofdm(_scaled) 0.6–1.9 s. **Done:** Herelink PSS fold search only for carrier hypotheses in 2370–2520 MHz (Herelink is 2.4 GHz only): digital_vtx 12.7 → 8.9 s (5.8 GHz Walksnail), 7.6 → 3.5 s (5.8 GHz Wi-Fi); identical results at 2.4 GHz. **Cost map (C++ ms per 2^21 block, 2 threads):** 5.8 GHz: digital_vtx 3.5–8.9 s (blind_links ~60 %, blind_hops ~40 %, both dominated by slice IFFTs of lengths with large prime factors), dji ~2.9 s (band-gated already: 2.4 / 5.1 / 5.8 GHz only), wifi 0.6–1.9 s; 2.4 GHz: digital_vtx 6.2 s, rc_link 4.8 s, dji 2.8 s; 940 MHz: nr_ssb 5.7 s, rc_link 4.9 s, lte 3.4 s; 3.5 GHz: nr_ssb 6.5 s. **Finding:** the slice lengths cannot simply be rounded to FFT-friendly sizes: rate_for() ties the rate to an integer number of samples per Tu. **Done 2026-10-23 (exact: identical results):** NR SSB PSS correlation over next_fast_len(L) instead of the next power of two above L + 256 (the lags used never wrap): nr_ssb 6.8 → 3.65 s (n78), 10.1 → 5.4 s (1840), 5.8 → 3.1 s (940), PCI / MIB unchanged; digital_vtx blind band gate for cellular bands (G8): 1840 MHz 2.6 → 0.65 s; digital_vtx slice cache (6 most recent, read-only) + exact fma remainder in window_phases (replaces fmod, bit-identical on 57 M cases): 5805 MHz 8.75 → 7.66 s. **Measured limits:** digital_vtx slices have medium prime factors (rate_for), Bluestein is no faster; the rc_link cost is the LoRa stage-1 scan (every BW/4 sub-band x 8 SFs x 2 polarities x 3 BWs) and the fsk_scan local median filter; LTE already uses next_fast_len. Larger gains need algorithm changes that alter detection statistics (validate on testdata + synthetic SNR sweeps first). **Next targets, by gain:** nr_ssb (5.7–6.5 s wherever NR can exist), rc_link LoRa scan (~4.8 s), lte (3.4 s), digital_vtx blind slices (cache / reuse identical slices, batch IFFTs), dji (2.9 s) | – | in progress |

## 3. To research before adding
| # | Item | Why it is blocked |
|---|---|---|
| M7r | remaining model gaps: ICAO register (403), sub-250 g DJI models have no DOC serial range, product codes > 70 not confirmed: bench 2026-10-22 on 65 real recordings (file 13 §7): Mini 2 decodes (63 confirmed again); the Mini 4 Pro sends 8-symbol DroneID frames that do not decode (O4 encryption), so 93 stays unconfirmed; needs a decodable newer DJI (O3: Mini 3 Pro, Air 2S, Mavic 3, Avata) | sources blocked / need real decodes |
| D2 | Walksnail BPSK header symbols (one per antenna) + 6-bit frame counter | bench 2026-10-22 (`bench/walksnail_header/FINDINGS.md`): the 34 ms capture holds 4 frames; a ~9-symbol block per frame start is partly structured but not clean BPSK; needs a >= 0.5 s continuous Walksnail capture (64-frame counter cycle), ideally known 1- / 2-antenna units |
| R2b | ELRS FLRC sync word → UID[2..5] | the fixed 21-bit field is unpublished: needs a real FLRC capture |
| D6 / M4b | O4 Air Unit (Lite, 1 TX) vs Pro (2 TX) via tx_antennas | bench 2026-10-22 (`bench/o4_streams/FINDINGS.md`): the real Mini 4 Pro O4 bursts give kurtosis 1.70 with no Alamouti pairing (spatial multiplexing / DFT-spread / reference symbols): no verdict; needs the O4 reference-symbol map and known 1-TX / 2-TX O4 captures |
| T8b | DroneID 8-symbol payload layout | layout unverified |
| G3b | ELRS Gemini Xrossband: pair the simultaneous 900 MHz SF6 and 2.4 GHz SF7 legs | both legs are already seen separately; needs a capture to set the pairing rule |
| G6b | 250 kb/s XN297 toys (V911S, KF606, MJXQ E010 / Phoenix) and non-nRF toy chips (E012 / E015 HS6200, E129 / E010 r5 RF2500, E016H v2) | new PHY work (250 kb/s path, HS6200 scrambler / guard bytes, RF2500 4x FEC) |
| B1 | Skydroid (H12 / H16 / T10 / T12 / G12 / H30) and SIYI (MK15 / MK32 / HM30) video downlinks: numerology, TDD | unpublished; the dominant Indian agri links: need bench captures of the air units |
| B2 | DJI O4+ (Mini 5 Pro / Mavic 4 Pro) numerology; DJI Neo / Flip Wi-Fi signatures | unverified claims (15 kHz, 2048 FFT, ZC 232); need captures |
| B3 | analog-receiver LO leakage (RX5808 / RTC6715: CW on even MHz = channel − 479 MHz, scanning) as a ground-station detector | hypothesis from UNKNOWN_SIGNAL_4 (2026-10-22); needs a goggles-only capture |
| B4 | Ghost / Tracer, ACCESS, AFHDS 3, FASST, Parrot MARS, Skydio Connect MH, RubyFPV | unpublished / encrypted: timing or band "-like" labels after bench captures |
| D8b | HDZero 17 MHz narrow mode (RunCam Nano 90) and Shark Byte numerology | no capture: needs an HDZero VTX with the Nano 90 camera in 540p90 (expect 11.1 ms frames) and a Shark Byte TX |
| D8c | HDZero video picture from the JSCC subcarriers | the sample-to-subcarrier mapping is unpublished; research only (would need equalisation via the pilots, then the mapping) |

## 4. Validation log
- Batch W (2026-10-25, Wi-Fi field gaps): pytest 380 (+1 skip), si_tests 129/129 (23 253 assertions; `vectors/wifi_field`
  new, `vectors/wifi` and `vectors/devices` regenerated), modularity 0. Cost (Python): HE-SIG-A ~2 ms per HE PPDU, BSS
  link profiles 0.34 ms per 150 PPDUs. Regression vs g8to14 (`out/regress_wbatch_diff.log`, with P1 and the DroneID
  fixes): changed on 4 of 50 inputs, all intended: D_20251119-183758 (real DJI Enhanced Wi-Fi, 5 MHz) gains
  ", drone link-like (DJI-like)" (identity: SZ DJI OUI); HE PPDUs at their HE-SIG-A bandwidth (80 MHz) form their own
  "(HE)" signals on LIVE_WIFI_5755MHz / LIVE_OFDM_VTX_5805MHz / real_ht_mcs1_5240 (two unknown energy pieces at
  5753 / 5755 MHz become one at 5751 MHz). Parity (`--refresh --threads 2`): **50 / 50 identical**.
- Batch G1–G14 (2026-10-23): pytest 368 (+1 skip), si_tests 124/124, modularity 0. Regression G1–G7 vs pal1 changed 0/0
  on all 49 inputs; G8–G14 vs G1–G7 changed 0 on 49 of 50 inputs, the one change being the intended G8 band widening:
  D_20260211-161341 gains a 15 kHz-SCS CP-OFDM link at 2081.0 MHz (14 MHz, ~10 % duty, 7 dB; "-like", not the documented
  Microhard link at 2069 MHz; unknown, worth a bench check). New input LIVE_MIMOSA_4880MHz_weak (G13) decodes. Parity
  (`--refresh --threads 2`) **50 / 50 identical** (rc_link 44/44, digital_vtx 6/6, adsb 2/2). Earlier parity attempt
  stalled only because the Mac slept on battery: runs now under `caffeinate -i`.
- 2026-10-21: pytest 224/224; regression diff only fixed 5 false COFDM links; parity 47/47.
- 2026-10-22 (after L1 / D7 / T7 / A7 / A5): pytest 242/242, si_tests 91/91, regression only BLE gains, parity 47/47.
- 2026-10-22 evening, batch 2 (FrSky, AFHDS 2A, XN297, blind frame, model level M0–M5 / M8, + earlier batch items):
  pytest 315/315, si_tests 105/105. Regression vs l1b: no change from the new RC decoders or the model fields; the
  DVB-T2 P1 search (D5) gave 6 false links on real captures (ADS-B, weak analog VTX, clipped 940 / 899 MHz) → fixed
  with an S1 / S2 margin gate (≥ 0.4; real P1 ≥ 0.74 at −6 dB); re-run on those inputs: changed 0/0. Parity (`--refresh --threads 2`): **47 / 47 identical** (rc_link 44/44).
- 2026-10-22 late, batch 3 (T4 802.15.4, R2d ELRS hop recovery, M6 DJI hop timing): pytest 326/326, si_tests all pass;
  regression vs batch2: no change from the batch; the only differences are the DVB-T2 P1 false links removed by the
  margin gate (8 in total: also 3 on Analog_20260804 near 5 GHz and 1 on D_20251119-165139). Parity 47 / 47 identical.
- Batch 4 (P1 step 1: Herelink band gate): regression vs batch3 changed 0/0 on every input (Python wall time e.g.
  Wi-Fi 5755 33.9 → 16.1 s, Mimosa 4880 39.1 → 23.5 s; the 5805 Walksnail run rose 48 → 75 s under 2-worker load, 15 s
  alone); parity 47 / 47 identical. C++ digital_vtx 12.7 → 8.9 s (5805), 7.6 → 3.5 s (5755).
- User 2026-10-22: complete every PLAN item first; keep each new capability in milliseconds (time it); broader
  optimisation (P1) later, when needed for specific signals.
- Batch 5 (D3 hop signatures, L1c XN297 soft demod, T5 FLARM): pytest 329, si_tests all pass; regression vs batch4
  changed 0/0 on every input; parity 47 / 47 identical.
- After batch 6 (2026-10-22): real-recording check of DroneID codes (file 13 §7) -> fixes: 8-symbol candidates include
  the Mini 4 Pro (measured), and a decoded product name from the external decoder ("Mini 2") is used as the model
  (it was reported as an unknown code). Unit tests + vectors pass; regression on the 3 DJI inputs changed 0/0, parity
  47 / 47 identical.
- Analog PAL colour (2026-10-22, user-confirmed PAL colour camera, `LIVE_ANALOG_VTX_PAL_5740MHz`): the colour test
  missed the real PAL burst (energy −0.7 dB: porch / sync-edge leakage in the 1.8 µs gate, references aliased);
  fixed with a per-line detrend and an in-band reference grid (Py + C++) → PAL, coherence 0.875, 6.7 dB. Regression vs hdzero1
  changed 0/0 on all 48 inputs; parity **49 / 49 identical**.
- D8 HDZero (2026-10-22): pytest 343 (+1 skip), si_tests 113/113; regression vs batch6 changed 0/0 on all 47 inputs;
  parity **48 / 48 identical** (new input LIVE_HDZERO_5658MHz: digital_vtx 5/5).
- Batch 6 (M7a–d, R2e, R3/T6; R4 / R6b / R10 / R11 not pursued): pytest 337, si_tests all pass; regression vs batch5
  changed 0/0 on every input; parity 47 / 47 identical.
