# Wi-Fi: field gaps and drone-vs-regular Wi-Fi (research 2026-10-25, PLAN batch W)

**User question (2026-10-25):** what can a field operator get from Wi-Fi that the library lacks, and can it tell drones
that use Wi-Fi from ordinary Wi-Fi?

## What the library already had (checked in the code, not assumed)
- **PHY:** 802.11b DSSS / CCK; 802.11a/g; 802.11n HT single stream, 20 MHz, BCC / LDPC / STBC; VHT-SIG-A of payloads
  out of reach; 802.11 at 5 / 10 / 11 MHz clocks, off-grid and mirrored; per-PPDU transmitter attribution.
- **MAC:** full dissection, OUI vendor, security, WPS / P2P identity, EAPOL, ARP, IP / ports.
- **Drone identity:** Remote ID (F3411 beacon and NAN, DJI DroneID IE, GB 46750), wfb-ng, ESP-NOW / ELRS Backpack,
  MAVLink over open Wi-Fi, 9 vendor OUIs and 13 SSID patterns.
- **Report-level `receivers/devices`:** device tracks, randomised-MAC linking (probe signature + CFO), findings:
  Remote ID / signature drones, deauth floods, duplicate SSID / evil twin vs virtual BSSIDs, weak security, WPS,
  handshakes, client history. The first gap list (in the conversation) wrongly called deauth / evil-twin detection and
  probe fingerprinting missing; they existed. Batch W only extends them.

## Gaps taken up
| # | Gap | Done | Validation |
|---|---|---|---|
| W1 | drone Wi-Fi link by **behaviour** (encrypted / undecodable traffic) | `wifi/field.py` per signal: BSS profiles + identity / behaviour verdict; `devices` finding `drone_link_behaviour` | synthetic only (research) |
| W2 | 802.11ax **HE-SIG-A** (BSS colour, UL/DL, MCS, streams, TXOP) | `ofdm_rx.parse_he_sig_a` + plugin fields | **real**: 44 / 44 HE PPDUs of `LIVE_WIFI_5755MHz` pass CRC-4 (proven) |
| W3 | more drone / DIY hotspot names | `droneid.SSID_PATTERNS` | vendor docs (below); frames synthetic |
| W4 | open-link payloads | Tello SDK command / state, H.264 / H.265 / RTP video in `mac80211` | Tello SDK guide; frames synthetic |
| W6 | beacon flood | `devices/threats._beacon_flood` | synthetic |
| W7 | MAC randomisation: sequence counter | `devices/tracker.link_randomized` | synthetic |

## Deferred
- **W5** single-stream 40 MHz HT and VHT 20 / 40 MHz payloads: they fit a 56 MHz Pluto capture, but no testdata holds
  any (2026-10-25 regression outputs: every real VHT PPDU is 80 MHz or 2-3 streams; no 40 MHz HT). Take up when a
  capture has them (e.g. a Wi-Fi drone in 40 MHz mode).
- **W8** Espressif Long Range (ESP-NOW LR, 250 / 500 kb/s, not 802.11; the PHY is not published) and 802.11ah HaLow
  (865-868 MHz in India, S1G OFDM 31.25 kHz): research first, no captures.

## Behaviour of Wi-Fi drones (W1 sources)
- Alipour-Fanid, Dabaghchian, Wang, Zeng, "Machine Learning-Based Delay-Aware UAV Detection and Operation Mode
  Identification over Encrypted Wi-Fi Traffic", IEEE TIFS 2020 (arXiv 1905.06396): WPA2-encrypted drone traffic
  (DJI, Parrot Bebop, DBPower ...) identified from packet sizes and inter-arrival times only, 0.15-0.35 s, 85.7-95.2 %;
  operation modes (standby, hover, flying ...) 88.5-98.2 %.
- Alsoliman, Rigoni, Levorato, Pinotti, Tippenhauer, Conti, "COTS Drone Detection Using Video Streaming
  Characteristics", ICDCN 2021: up to 99 % from ~170 drone packets (820 ms); separates drones from similar video
  sources (security cameras).
- What one 34 ms block allows: the shape (one client, heavy steady downlink, small uplink), not the frame rate. Hence
  the 0.55 "-like" verdict with an explicit caveat; identity (Remote ID, Tello SDK, MAVLink, signatures) names a drone.

## Hotspot names (W3, vendor documents)
| SSID | Equipment | Source |
|---|---|---|
| `ExpressLRS TX` / `ExpressLRS RX` (password `expresslrs`) | ExpressLRS module / receiver in Wi-Fi setup mode | expresslrs.org Web UI / Wi-Fi updating guide |
| `ArduPilot`, `PixRacer` | MAVESP8266 Wi-Fi telemetry bridge | ardupilot.org ESP8266 telemetry page |
| `CUAVWLINKxxxx` | CUAV PW-Link | ardupilot.org CUAV PW-Link page |
| `DroneBridge for ESP32` | DroneBridge ESP32 | docs.px4.io ESP32 Wi-Fi module page |
| `RMTT-xxxxxx` | DJI RoboMaster TT (Tello Talent, with the open-source controller) | Tello SDK 3.0 user guide; DroneBlocks |
| `HolyStoneFPV-xxxx` | Holy Stone Wi-Fi FPV drones | Holy Stone HS710 manual |
| `WiFi-720P-xxxx` | toy Wi-Fi FPV drones (Eachine E58 class) | Deviation forum (single source) |
Not added: `RCUFO_xxxx` (one non-authoritative guide), GoPro `GPxxxxxxxx` (a camera, not a drone), Skydroid H16
`Android AP` (Android's generic default).

## Tello SDK (W4)
Ryze / DJI Tello SDK 2.0 / 3.0 user guide: Tello is an open AP at 192.168.10.1; commands are ASCII to UDP 8889
(`command`, `takeoff`, `land`, `rc a b c d`, `wifi ssid pass`, `battery?` ...; replies `ok` / `error` / values); state
is `key:value;` text to UDP 8890 (`pitch roll yaw vgx vgy vgz templ temph tof h bat baro time agx agy agz`, `mid` with a
mission pad); video is raw H.264 to UDP 11111. Toy drones stream H.264 over UDP / TCP on vendor ports (8800, 8888;
hackaday ProMark VR log, HTB forum), so video is recognised by its NAL units, not by port.

## MAC randomisation (W7)
Vanhoef, Matte, Cunche, Cardoso, Piessens, "Why MAC Address Randomization is not Enough", AsiaCCS 2016: probe-request
element fingerprints and sequence numbers that continue across address changes link randomised MACs.

## The 12 Wi-Fi drones of the user's list (2026-10-25)
User figure (source paper not found): "OFDM -> Wi-Fi 4": Holy Stone HS720E, HS100G; "Wi-Fi -> Wi-Fi 4": DJI Mavic Air,
Spark, Phantom 3 Standard, Mavic Mini 1, Yuneec Q500-HD, Parrot Anafi, Ruko F11 GIM, Hubsan X4 Air, DJI Tello, Skydio 2.
All are standard 802.11 (the phone apps join an ordinary AP), so the existing 802.11b/g/n receivers decode their
frames; DJI "Enhanced Wi-Fi" (Spark, Mavic Air, Mavic Mini) also runs at a 5 MHz clock (Department 13 via Bender),
which `ofdm_scaled` handles (`testdata/D_20251119-183758`: 5 MHz DJI link, SZ DJI OUIs), with DroneID in the beacon
vendor element (OUI 26:37:12, decoded). None of the 12 is on DJI's FAA Remote ID list.

| Drone | Wi-Fi link | Identification in the library |
|---|---|---|
| Holy Stone HS720E / HS100G | drone -> phone FPV (HS720 manual: camera 5500-5700 MHz; HS100: 2.4 GHz); RC 2.4 GHz non-Wi-Fi | SSID `HolyStoneFPV-` / `_` (manuals) -> "Holy Stone (Wi-Fi FPV: HS720 / HS100 series)" |
| Ruko F11 GIM | drone -> phone 5.8 GHz "enhanced Wi-Fi" 720p | SSID `RUKO-F11-xxxxxxxx-5G` (Ruko support page) -> "Ruko F11" |
| DJI Mavic Air | Enhanced Wi-Fi 2.4 / 5.8, drone <-> RC | SSID `MavicAir-` -> "DJI Mavic Air" (+ SZ DJI OUI 0.9); DroneID beacons |
| DJI Spark | 2.4 / 5.8 Wi-Fi to phone or Spark RC | `Spark-` -> DJI Spark, `Spark-RC` -> remote controller; DroneID beacons |
| DJI Phantom 3 Standard | 2.4 GHz Wi-Fi video drone -> RC (RC re-broadcasts); control 5.8 GHz non-Wi-Fi | `PHANTOM3_` (manual) -> "DJI Phantom 3 Standard / 4K (Wi-Fi video models)" |
| DJI Mavic Mini 1 | Enhanced Wi-Fi drone <-> RC (MT1SS5 5.8 only; MT1SD25 2.4 + 5.8) | `Mavic-Mini-` -> "DJI Mavic Mini" (SSID format not verified), SZ DJI OUIs, DroneID beacons; 5 MHz clock |
| Yuneec Q500 | CGO3 camera -> ST10 5.8 GHz Wi-Fi video; control 2.4 GHz non-Wi-Fi | `CGO3_` / `CGO3P_` (forum) -> "Yuneec CGO3 camera video link"; Yuneec OUI E0:B6:F5:8 (MA-M); Q500-HD / CGO2 name not found |
| Parrot Anafi | 802.11 a/b/g/n 2.4 / 5.8 drone <-> Skycontroller 3 | `ANAFI-` -> "Parrot Anafi"; Parrot OUIs |
| Hubsan X4 Air H501A | Wi-Fi video to phone; RC 2.4 GHz | Hubsan OUI 98:AA:FC:7 (MA-M); SSID not found |
| DJI Tello | 2.4 GHz 802.11n, open AP | `TELLO-` / `RMTT-`, Tello SDK command / state, H.264 on UDP 11111 |
| Skydio 2 | 2.4 / 5 GHz Wi-Fi drone <-> phone / controller | Skydio OUI 38:1D:14; SSID per-unit, format not found |
No captures of these drones are in `testdata`; the signatures and the behaviour verdict (W1) are tested on synthetic
frames. Sources: manuals (Holy Stone HS720 / HS100, DJI Phantom 3 Standard, Mavic Mini, Tello, Skydio support), Ruko
support page, Kismet `kismet_uav.conf` and DJI IE parser, Wireshark manuf, Bender arXiv 2207.10795, Yuneec / Mavic forums.
