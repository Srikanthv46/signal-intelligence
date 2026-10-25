# 14. Market gap survey: VTXs, RC links and telemetry the library does not cover (2026-10-22)

Three web surveys checked market products from 2020 to 2026 against the library's coverage (README signal table,
PLAN §1, research files 01–13). Nothing here is bench-verified.

- **Tags:** [S] source code read, [D] datasheet, spec or vendor page, [M-3p] third-party measurement, [U] uncertain or
  inferred.
- **GPL code:** wfb-ng, SoftRF, the DIY Multiprotocol module and dump978 (mutability fork) are oracles only.
  Implement from the published descriptions. dump978 from FlightAware is BSD-2.
- **Rule:** every item still goes bench → Python → C++.

## 1. Gaps in what our own code does (cheap, verified in the code)
| Gap | Where | Note |
|---|---|---|
| wfb-ng / OpenIPC / OpenHD found only by SSID | `wifi/droneid.py` | wfb-ng sends no beacons, so the SSID rule never fires (see §2) |
| Blind OFDM search only at 2.33–2.53 and 4.9–6.0 GHz | `digital_vtx` `FPV_VIDEO_RANGES` | misses the 1.35–1.44, 1.43–1.444 (Chinese 1.4 GHz UAV band), 1.6–2.33 and 4.4–4.9 GHz professional / UAV bands |
| Analog channel tables stop at 5325–5945 MHz, no 3.3 GHz | `analog_vtx/detector.py` | detection is generic; this is a labelling gap only |
| RC catalog-only protocols (DSM, R9, HoTT, Hitec, RadioLink, S-FHSS) | `rc_link/protocols.py` | rate / deviation / period matching only; no demodulation or decode |
| ELRS LR1121 GFSK (K modes) | `rc_link` | in the catalog, no demod or sync decode |
| Off-grid 802.11 scan misses weak, low-duty links | `wifi/ofdm_rx.offgrid_plan` | the Mimosa 4880 MHz link at ~6.6 dB / 20 % duty was not found (UNKNOWN_SIGNAL_4, 2026-10-22); a forced decode got 21 L-SIGs |

## 2. Video transmitters
| System | Bands | PHY / what is published | Passive ID | Priority |
|---|---|---|---|---|
| **wfb-ng family** (OpenIPC FPV, OpenHD, RunCam WiFiLink, EMAX Wyvern Link) | 5.18–5.885 GHz, 2.4 GHz | Raw 802.11n injection: data frames FC 08 01 to broadcast, TA = SA = **57:42:xx:xx:xx:xx**, last 4 bytes channel_id = (link_id << 8) + radio_port. Packet type 0x01 data, 0x02 session (about once per second). Defaults MCS1, 20 MHz, BCC, RS FEC 8/12. RTL8812EU also injects at **10 MHz** [S] | **Exact.** link_id (per installation), ports, MCS / BW / STBC / LDPC, FEC structure. Payload ChaCha20-Poly1305; OpenIPC default keys [U] | **High** |
| **Skydroid** H12 / H12 Pro / H16 / T10 / T12 / G12 / H30 (Indian agri drones) | 2402–2480; H16 + 5180–5240; G12 + 5745–5825 [D] | RC hops measured by RFUAV (in `hop_model`); H16 video 20 MHz [M-3p]. Video downlink PHY unpublished | needs a bench capture → numerology → "-like" | **High** (bench) |
| **SIYI** HM30 / MK15 / MK15E / MK32 / FT24 | HM30 5100–5825 (FCC 5745–5825); MK 2.4 GHz | HM30: 8 / 4 Mb/s, H.265; PHY unpublished | bench | **High** (bench) |
| **ELRS Backpack** (side channel) | 2.4 GHz, ESP-NOW (channel 1 [U]) | 802.11 action frames (category 127, Espressif OUI 18:fe:34 [U]); destination MAC = the **ELRS bind UID**; MSP payload (VRX band / channel, head tracking) [S] | pilot identity tying the ELRS link, the VTX channel and the goggles | **High** |
| DJI O4+ (Mini 5 Pro, Mavic 4 Pro) | 2.4 / 5.1 / 5.8 | one repo claims 15 kHz, 2048 FFT, ZC root 232, length 1201 [U]; 2T4R | bench | Medium |
| DJI Neo / Neo 2 / Flip phone mode | 802.11 2.4 / 5.8 | Wi-Fi | SSID / OUI rules; check whether a Wi-Fi DroneID IE is present | Medium |
| Chinese 1.4 GHz UAV links (Taisync ViULinx, generic "1.4G") | 1430–1444 MHz | 10 MHz OFDM, TDD, LDPC [D] | blind OFDM once the band is searched; band-aware label | Medium |
| Parrot ANAFI UKR (MARS) / USA MIL | MARS 1.8–5.9 GHz hopping, TX ≠ RX frequency, 12.6 ms beamforming; Wi-Fi 6; LoRa backup 863–870; USA MIL: Microhard pDDL1800 | encrypted [D] | evidence-only class; Microhard 1.8 GHz band entry | Medium |
| Skydio X10 (Connect SL / MH) | MH: 1790–1850, 2040–2110, 2200–2390, 2400–2500; SL 2.4 / 5 GHz 2T4R [D] | unknown (Doodle Labs Helix-like [U]) | blind OFDM in the government bands | Medium |
| Analog **3.3 GHz** (Kimpok, FT Systems; 3170–3480 MHz) and **wide 4.9–6.0 GHz** (Axisflying 80 ch, FT 96 ch, groups A B E F R D X L J U O H) | – | analog FM (generic) [D] | channel-table labels | Medium |
| RubyFPV | 2.3 / 2.4 / 2.5 / 5.8 GHz injection + 433 / 868 / 915 | own framing, undocumented | capture; flag 802.11 at 2.3 / 2.5 GHz as non-ISM injection-like | Medium |
| Autel SkyLink 3 900 MHz mode (FCC only), Hubsan SyncLeas 3, Potensic PixSync 4 (AR803x [U]), FIMI RokLink 3, Yuneec / toy Wi-Fi SSIDs, HDZero 17 MHz | – | – | captures when a unit turns up | Low |
| Persistent MPU5, Mobilicom SkyHopper, Doodle Labs Helix, Silvus, DTC | L / S / C government bands, 5–20 MHz OFDM, MIMO | encrypted | numerology catalog entries + the wider search band | Low–Medium |
| Cellular dongles (DJI Cellular Dongle 2, Elsight), Starlink (Ku, out of range), 6.7 / 7.2 GHz analog (out of range) | – | – | not attributable / not visible | – |

## 3. RC links
| System | PHY (published) | Passive ID | Priority |
|---|---|---|---|
| **ELRS K modes / LR1121** (K1000 900 MHz, K1000 / DK500 / DK250 2.4 GHz) | GFSK 300 kb/s, fdev 100 kHz, 16-bit preamble, 1000 µs, OTA4 / OTA8 [S] | same sync / UID decode as the other ELRS modes | **High** |
| **ELRS Gemini Xrossband** (dual LR1121) | simultaneous 900 MHz LoRa SF6 BW500 CR4/8 + 2.4 GHz LoRa SF7 BW812, 6.67 / 10 ms [S] | pairing the two legs → "ELRS Gemini Xrossband" model label | **High** |
| **RadioLink** (AT9S / AT10 / T8FB / T8S; R7FG, R6FG, R8EF, R8FM, R4FGM) | CC2500 GFSK 38.0 kb/s, dev 38.1 kHz, 399.9 kHz grid from 2402 MHz; sync **0xABCD** (Surface / Air), **0xBADC** (DumboRC), **0xABA5** (RC4G); 34-byte packet; 16 hops via LCG; 19.9 ms (RC4G 14.46 ms) [S] | TX ID, channels; sync → sub-protocol | **High** |
| **FlySky AFHDS 1** (FS-CT6B, FS-T6) | A7105; 16 hops from an 8-row table; 21-byte packet; 1510 µs [S] | TX ID | **High** |
| **Hubsan** H107 / H501 | A7105; one fixed channel of 12 (2410–2465 MHz); 16 bytes; 10 ms [S] | ID | High |
| **Toy drones, nRF24 / XN297 / BK2425** (SymaX, MJXq, WLtoys V2x2 / V911S, E01X / E016H / KF606 / JJRC, Potensic, CX10) | per-protocol address + hop tables [S] | add to `xn297.py`; Wi-Fi toy SSIDs (WiFi-UFO, WIFI_UAV, …) [U] | High |
| **Spektrum DSM2 / DSMX / DSMR, Walkera DEVO, Multiplex M-LINK** | CYRF6936 DSSS 1 Mcps, SOP codes, 16-byte packets, 11 / 22 ms (files 03 / 08) [S] | mfg ID → unit; one DSSS demod covers ~14 protocols | Medium |
| FrSky R9 ACCST | LoRa SF6 BW500 CR4/5, preamble 9; FCC 43 ch / EU 19 ch; 26 B (EU 14 B); 20 ms [S] | region from the grid; framing decode | Medium |
| HoTT, Hitec Optima / Minima, Corona | CC2500 MSK 250 kb/s + FEC (HoTT sync 0xD391); Hitec sync = TX ID [S] | TX ID | Medium–Low |
| ImmersionRC Ghost, TBS Tracer | Semtech 2.4 GHz LoRa / FLRC; frame rates published, SF / BW / sync not [U] | timing "-like" after a bench capture | Low |
| FrSky ACCESS / ACCESS LR / TD / TWIN, FlySky AFHDS 3 / ANT, Futaba FASST / FASSTest / T-FHSS, JR DMSS | not published / encrypted [U] | blind detection only | Low |

## 4. Telemetry, data links and identification
| System | PHY / framing (published) | Passive ID | Priority |
|---|---|---|---|
| **MAVLink in LoRa payloads** (Ebyte E22 / E32 / E220 transparent links, 433 / 865 MHz) | LoRa SX126x / SX127x; air rate → SF / BW [D]; sync probably 0x12 [U] | the existing LoRa demod + the SiK MAVLink parser → sysid, HEARTBEAT, position | **High** |
| **MAVLink over open Wi-Fi** (UDP 14550 / TCP 5760; ESP8266 bridges, DroneBridge) and **ESP-NOW** links | 802.11; ESP-NOW action frames [D] | MAVLink decode on open networks; ESP-NOW transmitter MAC. DroneBridge is AES-GCM | High |
| **ADS-L 4 SRD860** (EASA, Issue 2, 2025-12-01) | M-band 868.2 / 868.4: FLARM PHY (100 kchip/s Manchester), sync 0x72 0x4B. O-band 869.525: GFSK 38.4 kb/s sync 0xB4 0x2B, or GMSK 200 kb/s sync 0x2D 0xD4. Key 0 = XXTEA with a zero key, i.e. scrambling only [D] | full decode: address, category (UAS), position, **Remote ID payload** | High (EU) |
| OGN tracker (OGNTP) | FLARM PHY; sync AA 66 55 A5 96 99 96 5A; 20 B + LDPC parity; clear [S] | address, type, position | Medium |
| FANET / FANET+ | LoRa SF7 BW250 (EU), sync 0xF1 [S] | parser on our LoRa output: address, aircraft type 7 = UAV, position | Medium |
| UAT 978 (US) | CPFSK 1.041667 Mb/s; RS-coded (dump978 BSD-2) | ICAO, UAV category, position | Medium (global) |
| China RID GB 46750-2025 (from 2026-05-01) | BT5 / Wi-Fi, not F3411-compatible (opendroneid-core-c `libopendroneidcn`, Apache-2.0) [U] | serial, position | Medium (global) |
| Microhard P900 / pDDL900, RFD multipoint / AES, PilotAware, LoRa-APRS, Mode A/C, Iridium content, FLARM payload, cellular UE uplinks | mostly unpublished / encrypted / out of scope | FHSS fingerprints at best | Low |

**India:**
- No broadcast Remote ID standard is notified. NPNT / Digital Sky need no live link.
- Finished-drone imports have been banned since 2022, but components are allowed. The local fleet is therefore
  domestic airframes (ideaForge, Garuda, IoTechWorld) with Chinese radios, Skydroid and SIYI in particular.
- The 865–867 MHz and 433 MHz licence-free bands carry LoRa telemetry and Meshtastic.

## 5. Sources (main)
- wfb-ng: https://github.com/svpcom/wfb-ng (`src/wifibroadcast.hpp`) · OpenIPC: https://docs.openipc.org/use-cases/fpv/net-cards/rtl8812eu/
- ELRS: https://raw.githubusercontent.com/ExpressLRS/ExpressLRS/master/src/src/common.cpp · Backpack: https://github.com/ExpressLRS/Backpack
- DIY Multiprotocol (GPL, oracle): https://github.com/pascallanger/DIY-Multiprotocol-TX-Module (RadioLink, FlySky, Hubsan, DSM, Devo, HoTT, Hitec, FrSkyR9 files)
- ADS-L: https://www.easa.europa.eu/en/downloads/142821/en · SoftRF (GPL, oracle): https://github.com/lyusupov/SoftRF
- Skydroid FCC: https://fccid.io/2ATGZQZYZH12 · RFUAV: https://arxiv.org/html/2503.09033v1 · SIYI HM30: https://discuss.ardupilot.org/t/30km-1080p-60fps-150ms-siyi-hm30-full-hd-digital-video-link-telemetry-radio-system/75189
- Parrot ANAFI UKR white paper: https://www.parrot.com/assets/s3fs-public/document/ANAFI%20UKR%20white%20paper.pdf · Skydio X10: https://skydio.com/x10/technical-specs
- Taisync: https://www.unmannedsystemstechnology.com/company/taisync-technology/viulinx-pro/ · 3.3 GHz VTX: https://ft.systems/en/?p=3992 · wide 4.9–6 GHz: https://ft.systems/en/?p=1941
- dump978: https://github.com/flightaware/dump978 · Ebyte E22: https://mischianti.org/ebyte-lora-e22-device-for-arduino-esp32-or-esp8266-specs-and-basic-usage-1/
