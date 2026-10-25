# 11 — Model-level identification within each drone transmitter family (research 2026-10-22)

**Goal:** go from "family" (already detected; see 01–10) to "model / hardware class" using only features a passive
receiver can measure from IQ. Web research only; no local captures used. This file does not repeat the PHY facts in
01–10; it adds what separates models.

**Tags:** [spec] confirmed standard / datasheet; [vendor] manufacturer page or manual (retailer copies of it included);
[code] source code (GPL code = facts only, never copied into C++); [comm] community measurement / review / paper;
[U] unverified or our inference. "Pwr" = rated maximum; on air we can only bound it if the range is known, so power
is a **weak** feature everywhere.

**Reporting rule for the library:** when the measured features match more than one model, report
`family` + `model_candidates` (list) + `model: unresolved`. Never pick one model. Label any inferred model `-like`
and keep the `research` validation level (PLAN.md V1) until a real capture proves it.

**Features the library already measures (refer to these by name below):**
- `BW`: occupied bandwidth.
- `NUM`: OFDM numerology (`system_candidates`).
- `FP`: blind frame period.
- `TXA`: `tx_antennas` (STBC / SFBC / 1 stream).
- `PPM`: carrier ppm.
- `CF`: centre frequency / channel-table match.
- `DID`: DroneID 8 / 9 symbols (`droneid_frame`), plus the payload's product type when it decodes.
- `RID`: Remote ID fields.
- `HOP`: blind hop width / on-time / period (D3, todo).
- `RC`: rc_link rate / deviation / grid / `protocol_candidates` / ELRS domain and `elrs_sync`.
- `AUD`: analog audio subcarrier.

## 0. Identity carried in the clear (strongest: decoded, not inferred)

### 0.1 DJI DroneID product-type byte (`device_type`, after home lat in the v2 record)
Sources:
- RUB-SysSec/DroneSecurity `droneid_packet.py` `DRONEID_DRONE_TYPES`
  (https://raw.githubusercontent.com/RUB-SysSec/DroneSecurity/public_squash/src/droneid_packet.py), **AGPL-3.0, oracle only**.
- Kismet `dot11_ie_221_dji_droneid.h` `product_type_str()`
  (https://raw.githubusercontent.com/kismetwireless/kismet/master/dot11_parsers/dot11_ie_221_dji_droneid.h), **GPL-2+, oracle only**.

The two tables agree on every code they share [code]. The code -> name pairs are facts. Re-typing them as a data
table is OK; copying the source files is not.

| Code | Model | Code | Model | Code | Model |
|---|---|---|---|---|---|
| 1 | Inspire 1 | 20 | N2 | 41 | Mavic 2 |
| 2, 3 | Phantom 3 series | 21 | Spark | 44 | M200 V2 series |
| 4 | Phantom 3 Std | 23 | M600 Pro | 51 | Mavic 2 Enterprise |
| 5 | M100 | 24 | Mavic Air | 53 | Mavic Mini |
| 6 | ACEONE | 25 | M200 | 58 | Mavic Air 2 |
| 7 | WKM | 26 | Phantom 4 series | 59 | P4 Multispectral |
| 8 | NAZA | 27 | Phantom 4 Adv | 60 | M300 RTK |
| 9 | A2 | 28 | M210 | 61 | DJI FPV |
| 10 | A3 | 30 | M210 RTK | **63** | **Mini 2** (confirmed by real decodes: NDSS 2023 Fig. 8; anarkiwi sample) [comm] |
| 11 | Phantom 4 | 31 | A3 AG | 64 | Agras T10 |
| 12 | MG-1 | 32 | MG-2 | 65 | Agras T30 |
| 14 | M600 | 34 | MG-1A | 66 | Air 2S |
| 15 | Phantom 3 4K | 35 | Phantom 4 RTK | 67 | M30 (RUB only) |
| 16 | Mavic Pro | 36 | Phantom 4 Pro V2.0 | 68 | Mavic 3 (RUB) |
| 17 | Inspire 2 | 38 | MG-1P | 69 | Mavic 2 Enterprise Advanced |
| 18 | Phantom 4 Pro | 40 | MG-1P RTK | 70 | Mini SE |

- **Unassigned or unknown codes:** 13, 19, 22, 29, 33, 37, 39, 42, 43, 45–50, 52, 54–57, 62, and everything above 70.
  - **No public code** for Mini 3 / 3 Pro, Mini 4 Pro, Air 3, Avata / Avata 2, Mavic 3 Classic / Pro or Neo.
  - Report an unknown code as `product_type=<n> (unknown)`, never as the nearest name.
- **Conflicting table, do not use:** llorencroma/DJIDroneIDspoofer gives, for example, mini_2 = 0x0e "as seen on Aeroscope".
  This contradicts the real decodes (Mini 2 = 63) [U].
- **Encryption status:**
  - Mini 2, Mini 3 Pro, Air 2S and Mavic 3 are listed as "O2/O3 (unencrypted)"; O4 / Mini 5 as encrypted
    (https://github.com/alphafox02/antsdr_dji_droneid) [comm].
  - A Mini 4 Pro has been shown decoding (https://www.rtl-sdr.com/tech-minds-video-on-dji-drone-detection-on-the-antsdr-e200/) [comm].
  - EdgeSource claims encryption since January 2024 on Mavic 3 / Mini 4 Pro with the latest firmware
    (https://www.linkedin.com/pulse/security-risks-aeroscope-upgrade-module-edgesource-ofgbe) [U].
  - Conclusion: status depends on firmware. Keep `droneid_payload_state` = decoded / opaque.
- **The serial gives the model even when the product-type code is unknown.** See 0.2: DJI serials 1581F + model prefix.
- **8 vs 9 symbols (D4): conflict to resolve.** File 09 puts Mini 2 at 9 symbols and Mavic 2 and older at 8. One
  reading of NDSS says "Mini 2 and older omit symbol 1". Treat `droneid_frame` as a generation hint only, and let the
  decoded product type win.

### 0.2 Remote ID (F3411): UA type and CTA-2063-A serial -> manufacturer
- **UA type (Basic ID, 4 bits)** (opendroneid-core-c `opendroneid.h`, Apache-2.0,
  https://raw.githubusercontent.com/opendroneid/opendroneid-core-c/master/libopendroneid/opendroneid.h) [code mirroring spec]:
  - 0 none, 1 aeroplane, 2 helicopter / multirotor, 3 gyroplane, 4 hybrid VTOL;
  - 5 ornithopter, 6 glider, 7 kite, 8 free balloon, 9 captive balloon, 10 airship;
  - 11 free fall / parachute, 12 rocket, 13 tethered powered, 14 ground obstacle, 15 other.
- **ID type:** 1 = serial (CTA-2063-A), 2 = CAA registration, 3 = UTM UUID, 4 = session ID.
- **Serial format:**
  - 4-character manufacturer code (assigned by ICAO), then 1 length character (1–9, A–F = 1–15), then up to 15
    characters; the character set excludes O and I.
  - Sources [spec, via secondary]: https://datatracker.ietf.org/doc/html/draft-wiethuechter-drip-uas-sn-dns-02 and the
    Parrot staff breakdown https://forum.developer.parrot.com/t/anafi-ai-serial-number-dri-breakdown/22704.
  - The ICAO public register (https://www.icao.int/airnavigation/international-aviation-trust-framework-iatf/remote-id-number-registration)
    returned 403. **Fetch it manually.**

**Manufacturer codes:**

| Code | Maker | Src |
|---|---|---|
| 1581 | DJI (1581F…) | [vendor] https://support.dji.com/help/content?customId=en-us03400007747&spaceId=34&re=US&lang=en |
| 1748 | Autel (1748C…) | [comm] https://dronexl.co/2023/03/27/remote-id-is-your-drone-remote-id-compliant/ |
| 1588 | Parrot (1588E…) | [vendor staff] Parrot forum above |
| 1668 | Skydio (1668B…) | [comm] DroneXL above |
| 1596 | Dronetag | [vendor] https://help.dronetag.com/dronetag-dri/integration/usa-manufacturer/ |
| 1863 XAG, 2104 Talos, 1865 Holy Stone, 1910 Potensic, 1587 AgEagle/senseFly, 1792 uAvionix, 1817 Freefly, 1849 Spektrum, 1930 SwellPro, 1972 Lumenier, 1973 FrSky, 2071 NewBeeDrone, 1892 SJRC, 1875 Bwine, 1710 Atom UAV, 1902 Elsight | | [comm] compiled from FAA DOC list, https://jrupprechtlaw.com/remote-id-compliance-research/ |
| 1809 | Aerobits **or** idME (conflict) | [U] |

Not found: Yuneec, BetaFPV.

**DJI model from serial:**
- The RID serial equals the aircraft serial (20 characters). For 14-digit flight-controller serials it is
  1581F + FC serial + "0": Air 2S, M300, Agras T30/T40, P4P V2, Mavic 2, DJI FPV, Mavic Air 2.
  Sources: DJI FAQ above; https://dspalliance.org/dji-rid-doc-issues/ [vendor via 3rd party].
- **Model prefixes** [U, forum snippets only]:
  - Mavic 3 = 1581F45T (FAA DOC range);
  - Mavic 3 Classic = 1581F67P;
  - Mavic 3 Pro = 1581F67Q;
  - Mini 3 Pro = 1581F4XF or 1581F6L5;
  - Air 3 = 1581F6N8.
- **Authoritative per-model ranges:** the FAA DOC list https://uasdoc.faa.gov/listdocs. **Fetch it to build the table.**

**Library action:**
- Add `manufacturer` (from the MFR code) to `remote_id`.
- Add `model_candidates` from a serial-range table, built from the FAA DOC list.
- An add-on broadcast module (Dronetag, Aerobits, …) names the module maker, not the airframe. Report it as
  `rid_module`.

## 1. DJI (OcuSync 2 / O3 / O3+ / O4 / O4+, FPV, Wi-Fi)

### 1.1 Aircraft
EIRP limits are given as FCC / CE (dBm). Source: dji.com spec and support pages per model [vendor].

| Model | System | Bands (GHz) | EIRP 2.4 | EIRP 5.8 | TX chains | Stated BW / mode | Src |
|---|---|---|---|---|---|---|---|
| Mini 2, Mini 4K, Mini 2 SE | O2 | 2.4, 5.8 | 26/20 | 26/14 | not stated | 720p30 | dji.com/my/support/product/mini-2, /mini-2-se |
| Mini 3 | O2 | 2.4, 5.8 | 26/20 | 26/14 | 1T2R | 720p30 | dji.com/mini-3/specs |
| Mini SE | "Enhanced Wi-Fi" | 5.8 (SS5) or 2.4+5.8 (SD25) | <19 CE | <30 FCC | 1T | 720p30 | dji.com/my/support/product/mini-se |
| Mavic Air 2, Mavic 2 Pro/Zoom | OcuSync 2.0 | 2.4, 5.8 | 26/20 | 26/14 | not stated | 1080p30, 12 Mb/s | support pages |
| Mini 3 Pro | O3 | 2.4, 5.8 | 26/20 | 26/14 | 1T2R (DJI) vs 2T4R (retailers) [U] | **1.4/3/10/20/40 MHz** | dji.com/my/support/product/mini-3-pro |
| Air 2S | O3 | 2.4, 5.8 | **30**/20 | **30**/14 (SRRC **29**) | 2T4R | 1080p30 | dji.com/my/support/product/air-2s |
| Mavic 3 / Classic / Pro | O3+ | 2.4, 5.8 | 33/20 | 33/14 | 2T4R | 1080p60 | support / spec pages |
| Avata | O3 | **5.8 TX; 2.4 RX only** | – | 33/14 | 2T2R | 1080p100, ≤40 MHz | dji.com/my/support/product/avata |
| DJI FPV | DJI FPV (O1-class) | 2.4, 5.8 | 31.5/20 | 31.5/14 | 2T [U] | 810p120, ≤40 MHz, 50 Mb/s | dji.com/my/support/product/dji-fpv |
| Mini 4 Pro, Air 3, Air 3S | O4 | 2.4, **5.1**, 5.8 | 33/20 | 33/14 (5.1: 23 CE) | 2T4R | 1080p60 | spec pages |
| Mini 5 Pro, Mavic 4 Pro | O4+ | 2.4, 5.1, 5.8 | 33/20 | 33/14 | 2T4R | 1080p60 | spec pages |
| Avata 2 | O4 | 2.4, 5.1, 5.8 | 33/20 | 33/14 | 2T4R | **≤60 MHz**, 1080p100 | dji.com/avata-2/specs |
| **Neo** | O4 (+ Wi-Fi) | 2.4, 5.1, 5.8 | **26**/20 | **26**/14 | **1T2R** | ≤40 MHz | dji.com/neo/specs |
| **Flip** | O4 (+ Wi-Fi) | 2.4, 5.1, 5.8 | 30/20 | 33/14 | **1T at 2.4, 2T at 5.x** | – | dji.com/flip/specs |
| Matrice 30 | O3 Enterprise | 2.4, 5.8 | 33/20 | 33/14 | – | – | enterprise.dji.com/matrice-30/specs |
| Matrice 350 | O3 Enterprise | 2.4, **5.1**, 5.8 | 33/20 | 33/14 | 2T4R | – | enterprise.dji.com/matrice-350-rtk/specs |
| Matrice 4 | O4 Enterprise | 2.4, 5.1, 5.8 | 33/20 | 33/14 | 2T4R | – | enterprise.dji.com/matrice-4-series/specs |
| Phantom 4 Pro V2 | OcuSync 1-class | 2.4, 5.8 | – | – | – | 720p/1080p30 | support page |

**Published passive measurements:**

RFUAV Table 4, uplink hop width / on-time / interval / pattern period (https://arxiv.org/html/2503.09033v1) [comm, one unit each]:

| Model | Hop width | On-time | Interval | Pattern period | Video BW in the capture |
|---|---|---|---|---|---|
| DJI FPV | 5 MHz | 0.64 ms | 4 ms | 38.3 ms | 10 MHz |
| Avata 2 | 7 MHz | 0.41 ms | 2 ms | 20.6 ms | 10 MHz |
| Mavic 3 Pro | 4.78 MHz | 1.7 ms | 6 ms | 60 ms | 10 MHz |
| Mini 3 | 3.5 MHz | 0.56 ms | 5.96 ms | 40.0 ms | 10 MHz |
| Mini 4 Pro | 6.54 MHz | 0.404 ms | 2.5 ms | 24.05 ms | 10 MHz |

- Video BW was 10 MHz for every DJI model in the capture, so **video BW is not a model feature** (it is adaptive).

**DJI decision table:**

| Measured | Model candidates | Conf. | We measure? |
|---|---|---|---|
| DroneID decoded, product type known | that model | high | `DID` payload (external decoder) |
| DroneID / RID serial 1581F + prefix | model by FAA DOC range | high once table built | `RID` / `DID` |
| NUM = O2 15 kHz, 5.8 or 2.4 | O2 group: Mini 2 / 4K / 3 / 2 SE, Air 2, Mavic 2 → unresolved; HOP 3.5 MHz / 6 ms / 40 ms → Mini 3-like | med (HOP: one unit) | `NUM`; `HOP` todo |
| NUM = O3 (2048, SFBC 2T), no 5.1 | Air 2S, Mavic 3 family, Avata, M30, (Mini 3 Pro if 2T) → unresolved; HOP 4.78 / 1.7 / 6 ms → Mavic 3 Pro-like | med | `NUM`, `TXA`; `HOP` |
| O3 video on 5.8 GHz, **no 2.4 GHz downlink**, BW ≤ 40 MHz | Avata (or O3 Air Unit, §2) | med | `CF`, `BW` |
| NUM = O4 30 kHz | Mini 4 Pro / Air 3 / 3S / Avata 2 / Neo / Flip / M4 / O4 Air Unit | high (generation) | `NUM`, `FP` 2.5 ms |
| O4 + **1 TX** at 5.8 | Neo, O4 Air Unit (Lite), Flip at 2.4 only | med | `TXA` (D6 todo) |
| O4 + **60 MHz** BW | Avata 2, O4 Air Unit / Pro (with Goggles 3/N3) | med | `BW` |
| O4 HOP 7 MHz / 0.41 / 2 ms | Avata 2-like; 6.5 MHz / 0.40 / 2.5 ms → Mini 4 Pro-like | low–med (one unit each) | `HOP` |
| NUM = O4+ (15 kHz / 2048 / ZC 1201 [U]) | Mini 5 Pro, Mavic 4 Pro → unresolved | low | `NUM` (needs entry) |
| 40 MHz, 810p-class, DJI FPV channel plan | DJI FPV drone / FPV Air Unit / Vista (§2) | med | `BW`, `CF` |
| 802.11a at a 5 MHz clock + IE 26:37:12 | Spark / Mavic Air / Mini / Mini SE → product-type byte in the IE | high when decoded | wifi |

**Indistinguishable on air (report unresolved):**
1. Mini 2 / Mini 4K / Mini 2 SE (O2, same EIRP), unless HOP or DroneID differ.
2. Mavic 3 / 3 Classic / 3 Pro (O3+).
3. Mini 4 Pro / Air 3 / Air 3S / Mini 5 Pro / Mavic 4 Pro: same bands and 2T4R. HOP timing (RFUAV) is the only
   published separator, from one unit each.
4. Matrice 4 vs consumer O4 [U].

## 2. FPV digital air units: DJI FPV / O3 / O4, Walksnail Avatar, HDZero

### 2.1 DJI FPV-class air units
Sources: https://www.dji.com/au/support/product/o3-air-unit, https://www.dji.com/o4-air-unit/specs [spec];
https://oscarliang.com/caddx-vista-vtx [comm].

| Unit | TX | Channels at 20/10 MHz | 40 / 60 MHz | EIRP FCC | Notes |
|---|---|---|---|---|---|
| DJI FPV Air Unit (V1) | 2 | 5660/5695/5735/5770/5805/5878/5914 + CH8 5839 (CE: 5735/5770/5805/5839) | 40 MHz on CH1/2/3/8 | <30 | 720p120 |
| Caddx Vista, RunCam Link / Wasp | **1** | same | same | <30 | identical firmware |
| O3 Air Unit | 2T2R [spec], one dual-polarised antenna [U] | **5768.5 / 5804.5 / 5839.5** | 40: 5794.5 | <33 | 2.4 GHz RX only |
| O4 Air Unit | **1T1R** | **5768.5 / 5789.5 / 5814.5** | 40/60: 5794.5; plus 5.170–5.250 | **<30** | |
| O4 Air Unit Pro | **2T2R** | same as O4 | same | **<33** | |

### 2.2 Walksnail / Caddx Avatar HD
Source: caddxfpv.com product pages [spec].

| Model | TX ant. | EIRP FCC | Pwr steps | Src |
|---|---|---|---|---|
| Avatar HD V1 / Pro Kit | 2 | <30 | 25/200/700 (unlock 1000/1200) mW | https://caddxfpv.com/blogs/news/walksnail-avatar-system |
| Avatar HD Kit V2 | 2 (dual SKU) or 1 (single SKU) | <30 | same | https://caddxfpv.com/products/walksnail-avatar-hd-kit-v2-dual-antenna |
| Avatar GT / GT2 | GT 2, GT2 1 [vendor] | <30, max 33 (2 W) | + "2000mW_Auto" (fw 39.44.15) | https://caddxfpv.com/products/walksnail-avatar-gt-kit-only |
| Moonlight | 2 | <30 | same | https://caddxfpv.com/products/walksnail-moonlight-kit |
| Nano V3 | 1 | <27 (500 mW) | | https://caddxfpv.com/products/walksnail-avatar-hd-nano-kit-v3 |
| Mini 1S / 1S Lite | 1 | <25.5 (350–500 mW) | | https://caddxfpv.com/products/walksnail-avatar-hd-mini-1s-lite-kit |
| **Ascent Lite / Lite+** | 1 | <20 (100 mW) / 500 mW | **5 / 10 / 20 MHz / Auto**; 5.730–5.815 GHz | https://nextfpv.com/products/ascent-lite-kit [vendor] |

- **Channels and bandwidth** (https://oscarliang.com/fpv-channels, https://oscarliang.com/setup-avatar-fpv-system) [comm]:
  - Avatar HD uses the DJI FPV V1 channel plan.
  - **25 Mb/s "Standard" = 20 MHz** on CH1–8. **50 Mb/s "High" = 40 MHz** on 5695/5770/5878/5839.
  - **Race mode** (fw 34.40.15+): raceband R1–R8, 540p.
  - Goggles transmit on CH8 (5839) until they link.
- **Video frame rate (1080p60 / 100 / 120) does not change the RF link** [U]. Our `FP` 7.163 ms is the radio frame,
  not the video rate. Moonlight 4K is a recording mode only.

### 2.3 HDZero (Divimath DM6300)
Source: HDZero VTX manual v0.7 (https://flyingtech.co.uk/wp-content/uploads/2025/05/HDZero_VTX_UserManual_v0.7.pdf)
[vendor]; firmware https://github.com/hd-zero/hdzero-vtx (`dm6300.c`: BW_27M / BW_17M) [code].

- **PHY measured** (2026-10-22, Freestyle V2, 27 MHz mode): file 02 §3; the library names it by its pilots (`hdzero`).
- **One TX antenna on every model.**
- **Power ceilings:**
  - 25/200 mW on Whoop, Whoop Lite, Race V1/V2/V3, ECO, AIO5/AIO15;
  - Freestyle V1/V2 add 500 mW / 1 W (HAM unlock);
  - Shark Byte TX5S.1 200 mW, TX5M.1 500 mW.
- **17 MHz mode** only with the RunCam Nano 90 camera in a 540p60 menu setting [comm/U]; otherwise 27 MHz.
- **Channels:** R1–8 raceband, E1 5705, F1 5740, F2 5760, F4 5800, L1–8 5362–5621.
- **Possible synthesiser offsets:** a table read from the firmware suggests a few MHz off nominal on some channels [U].
  **Check on air:** it would be a firmware-generation fingerprint, not a model fingerprint.

### 2.4 FPV decision table

| Measured | Candidates | Conf. | We measure? |
|---|---|---|---|
| `NUM` walksnail_like + `TXA` = 2 | Avatar HD V1/Pro, V2 dual, GT, Moonlight | high (antenna count) | `NUM`, `TXA` |
| `NUM` walksnail_like + `TXA` = 1 | Mini 1S, Nano V3, V2 single, GT2 | med | `TXA` |
| walksnail_like at 40 MHz | any Avatar HD in High (50 Mb/s) mode | med | `BW` |
| walksnail_like on a raceband centre | Avatar HD in Race mode | med | `CF` |
| AR8030-like 5 / 10 MHz in 5.730–5.815 | Walksnail Ascent family | low–med | `BW`, `CF` |
| O3 numerology at 5804.5 / 5839.5 | O3 Air Unit, or O3 aircraft (Avata) | high (raster) | `CF`, `NUM` |
| O4 at 5789.5 / 5814.5, or 60 MHz, or 5.1 GHz | O4 Air Unit / Pro / O4 drone; `TXA` 1 → O4 Lite, 2 → Pro | high (generation), med (Lite vs Pro) | `CF`, `BW`, `TXA` (D6) |
| DJI FPV plan, O1 numerology, `TXA` 2 vs 1 | FPV Air Unit / DJI FPV drone vs Vista / Link | med | `TXA` |
| hdzero_like 27 MHz | any HDZero / Shark Byte → unresolved | – | `BW` |
| hdzero_like 17 MHz | HDZero + Nano 90 camera (Whoop / Lite / AIO class) | low–med | `BW` |

**Identical on air:**
- every HDZero VTX (power only);
- Vista = RunCam Link;
- Avatar Mini 1S = Nano V3 = V2 single-antenna (power only);
- Avatar V1 Pro = V2 dual = Moonlight; GT only above ~1.2 W.

## 3. Analog VTX
Sources:
- RTC6705 datasheet https://der-frickler.net/_media/modellbau/fpv/sky-rf/rtc6705-dst-001.pdf [spec];
- Betaflight `vtx_table.c` / `vtx_rtc6705.c` [code, GPL facts];
- TBS manuals and vendor pages as cited.

| Measured | Candidates | Conf. | We measure? |
|---|---|---|---|
| Carrier at **5584 MHz** with very low EIRP | TBS Unify Pro pit mode (https://team-blacksheep.com/products/prod:unify_pro_hv_sma) [vendor] | med | `CF` |
| 6.5 MHz audio subcarrier | VTX with a microphone: Unify Pro32 HV, Unify Evo, Rush Tank Solo / Max Solo, HGLRC Zeus (vendor pages, 01/agent sources) | low (many) | `AUD` (A5) |
| No audio subcarrier | Unify Pro32 Nano, many AIO / whoop VTX [U] | low | `AUD` |
| 6.0 + 6.5 MHz both | RTC6705-class with both audio channels [spec]; 1.2 GHz 9-channel modules [vendor] | low | `AUD` |
| 1080–1360 MHz grid (9 ch); only 1258 / 1280 | 1.2/1.3 GHz VTX; 2-channel locked VTX-1G3SE (https://www.getfpv.com/vtx-1g3se-2ch-video-transmitter-1-2ghz-1-3ghz.html) | med | `CF` |
| 3170–3480 MHz, 20 MHz grid | Kimpok-style 3.3 GHz (https://megadron.pl/en/products/kimpok-analog-transmitter-kp-t3335-33607.html) | med | `CF` (table to add) |
| 4.9–6 GHz extended / X band | Foxeer Reaper Extreme V3-class (https://rcdrone.top/products/foxeer-4-9g-6g-reaper-extreme-v3-2-5w-80ch-vtx) | low | `CF` |
| Carrier offset from integer MHz | the unit's crystal (Betaflight sets integer MHz; RTC6705 step 40 kHz [U]) | unit, not model | `PPM` |

**Not observable:** SmartAudio / Tramp / MSP control (wired UART), PAL / NTSC (set by the camera), the power
setting. Deviation differs between VTXs (C5VRX uses 10/12/16 MHz filter classes, https://github.com/Twotoz/C5VRX
[comm]), but **no per-model table exists**.

→ Report `analog_vtx, model unresolved` + `channel band` + `audio yes/no`.

## 4. RC links

### 4.1 ExpressLRS
Sources: `common.cpp` rate table, `FHSS.cpp` domains (https://github.com/ExpressLRS/ExpressLRS, GPL, facts only)
[code]; https://www.expresslrs.org/info/signal-health/ [vendor]. Domain grids and the sync channel are already in 07.

| Measured | Hardware class | Conf. | We measure? |
|---|---|---|---|
| FLRC 650 kb/s (F1000 / F500 / D250 / D500) | **SX1280** 2.4 GHz | high | rc_link FLRC class |
| GFSK 300 kb/s, 100 kHz dev at 2.4 GHz (K1000, DK250/500) | **LR1121** (LR2021 on master) | high | `RC` (rate 300k exists as elrs_gfsk_like) |
| 900 MHz LoRa **SF5** (200 Hz Full, 250 Hz) or 900 MHz GFSK K1000 | **LR1121** | high | `RC` (LoRa SF) |
| 900 + 2.4 GHz LoRa simultaneous, same timing (X100 / X150) | LR1121 dual-band TX | high | needs cross-band pairing |
| Two simultaneous copies in one band (Gemini) | dual-SX1280 / dual-LR1121 Gemini TX (BetaFPV SuperG, RadioMaster Nomad, …) | med; placement of the second copy [U] | todo |
| 900 LoRa SF6–SF9 (25–200 Hz) | SX127x **or** LR1121 → unresolved | – | `RC` |
| 2.4 LoRa SF5–SF8 | SX1280 **or** LR1121 → unresolved | – | `RC` |
| Channel count 40 / 20 / 13 / 4 / 3 / 8 / 80 | FCC915 / AU915 / EU868 / IN866 / AU-EU433 / US433-TH920 / ISM2G4 build | high | `RC` domain |

**Identical on air:** every module with the same chip, rate and domain:
- 2.4 GHz SX1280: Ranger / Micro / Nano, ES24TX, BetaFPV Micro / Nano, TX16S, Boxer, Zorro, Pocket internals;
- 900 MHz SX127x: ES900TX, BetaFPV 900, R9M flashed with ELRS.

They differ only by power, ppm and antenna. EU_CE_2400 LBT uses the same 80-channel grid as ISM2G4 [code], so it is
not separable. Unit identity: the UID bits from 07 (`elrs_identity`, `elrs_sync`) plus `PPM`.

### 4.2 TBS Crossfire / Tracer
Source: https://www.g3gg0.de/wordpress/fpv/fpv-analysis-of-tbs-crossfire/ [comm]; details in 08.

- Micro TX V2 (250 mW, 1 W after a firmware update), Nano TX (≤1 W), full TX (2 W [U]) and Tango 2 / Mambo internal
  are the **same protocol**; only power differs (https://oscarliang.com/tbs-crossfire-micro-tx-v2/,
  https://www.team-blacksheep.com/products/prod:xf_nano_tx [vendor]).
  → `crossfire, model unresolved`.
- **Plans that can be told apart:**
  - 868 / 915 (260 kHz grid at x.165 MHz, uplink 0–49, downlink +50);
  - **CE 863.093 + n·0.112 MHz**;
  - **AU 915.165–927.905**;
  - RACE (100 shared channels).
  These give the region build, not the model.
- **Tracer:** 2.4 GHz, 250 Hz, 25/100 mW (https://oscarliang.com/tbs-tracer/) [vendor]. Modulation unpublished
  (FLRC? [U]). Needs a capture.

### 4.3 FrSky, FlySky, Spektrum, Futaba, HoTT, Hitec
Sources: Betaflight `cc2500_frsky_*.c`; Multiprotocol `*_cc2500.ino`, `AFHDS2A_a7105.ino`, `DSM_cyrf6936.ino`;
DeviationTX `sfhss_cc2500.c` [code, GPL facts; derived rates assume a 26 MHz crystal]. Pages:
https://github.com/pascallanger/DIY-Multiprotocol-TX-Module, https://github.com/betaflight/betaflight,
https://github.com/DeviationTX/deviation.

| Measured (2.4 GHz unless noted) | Protocol / variant | Model meaning |
|---|---|---|
| CC2500 grid ~300 kHz from ~2404 MHz, 47 hops, ~31 kBd, 25 B | FrSky D8 | any D8 TX / multiprotocol; unresolved |
| same grid, ~70 kBd GFSK, 30 B, 9 / 18 ms (8 / 16 ch) | D16 ACCST v1 FCC | " |
| ~100 kBd 2-FSK, 35 B, grid from ~2405 | D16 v1 **EU-LBT** | region build |
| ~77 kBd GFSK, HW CRC on, extra bind byte | **ACCST v2** FCC (LBT: 2-FSK, 1500 µs turnaround) | firmware generation |
| ACCESS | encrypted, PHY unpublished [U] | – |
| 900 LoRa SF6 BW500 CR4/5, 20 ms, 43 ch (FCC) or 19 ch EU 250 kHz | FrSky R9 (as Multiprotocol implements it; stock = [U]) | R9M / R9M Lite / Lite Pro unresolved (power only) |
| A7105, 16 of 164 channels, 3.85 ms, TX 38 B | FlySky AFHDS 2A | FS-i6 = i6X = any AFHDS 2A TX |
| AFHDS 3 (PL18, NV14, FTr) | PHY unpublished [U] | – |
| CYRF6936 1 Mcps, 2 fixed channels | DSM2 | |
| CYRF6936, 23 channels | DSMX | 11 vs 22 ms frame = TX setting (not model) |
| CC2500 ~128 kBd, 30 ch at 1.5 MHz, 6.8 ms, 13 B | Futaba S-FHSS | |
| FASST / T-FHSS | PHY unpublished [U] | – |
| CC2500 MSK ~250 kBd, sync 0xD391, FEC, 75 ch, 10 ms, 50 B | Graupner HoTT | |
| CC2500 ~250 kBd, ~1.8 kHz dev, 21 ch, 22.5 ms | Hitec Optima / Minima | |

**Multiprotocol 4-in-1 vs genuine TX:**
- No published on-air fingerprint.
- Multiprotocol exposes a user fine-tune for the CC2500 frequency (and the A7105 LO), so a module may carry a
  deliberate carrier offset [code]. It is per unit, not a signature [U].
- Report `protocol, transmitter unresolved`.

### 4.4 mLRS
Source: https://github.com/olliw42/mLRS `common_conf.h` [code], https://ardupilot.org/sub/docs/common-mlrs-rc.html [vendor].

- **Fingerprint:** 91-byte frames in both directions.
- **2.4 GHz:** SX128x LoRa BW800, CR 4/5, preamble 12:
  - SF5 50 Hz (24 ch), SF6 31 Hz (18 ch), SF7 19 Hz (12 ch);
  - also FLRC 111 Hz.
- **900 MHz:**
  - SX126x SF5 BW500 = 31 Hz, **SX126x only**;
  - SF6 BW500 = 19 Hz, SX126x or SX127x;
  - **GFSK 100 kb/s, 50 kHz dev = FSK 50 Hz, SX126x only**.
- **Hop counts:** 915 FCC 25, 868 6, 866 IN 3, 433 2.
- A 31 Hz SF5 / FSK mode → SX126x hardware.

## 5. Telemetry / data links
Sources: ArduPilot SiK docs
(https://ardupilot.org/plane/docs/common-3dr-radio-advanced-configuration-and-technical-information.html) [vendor];
irlock RFD pages; Microhard, Doodle Labs, CubePilot, SIYI pages as cited [vendor].

| Measured | Candidates | Conf. |
|---|---|---|
| SiK GFSK, air rate in {2,4,8,16,19,24,32,48,64,96,128,192,250} kb/s, 433 / 915, grid (MAX−MIN)/(N+2) | any 100 mW SiK (3DR V1/V2, Holybro V3, mRo, HM-TRP) **or** original RFD900 → unresolved | high (family) |
| Air rate in {12, 56, 100, 125, 200, 224, 500, 750} kb/s | RFD900x / ux / 868x (3.x firmware) (https://irlock.com/products/rfd-900x-modem-bare) | med [U inference] |
| FHSS 868–870 MHz | RFD868x / 868ux | med |
| 4 or 8 MHz OFDM in 2402–2478, 2 streams | Microhard pMDDL2450 (2×2) (https://microhardcorp.com/pMDDL2450.php); single stream → pDDL2450 [U] | med |
| 5/10/20/40 MHz, 2 TX, 2.4 GHz | pMDDL2460 (https://www.microhardcorp.com/pMDDL2460.php) | low |
| 18/8/4 MHz at 3.2–4.8 GHz | pMDDL4000 | med |
| COFDM 26 or 40 MHz at 2.4 GHz | Doodle Labs RM2450-class (datasheet via digikey) | med |
| COFDM 3/5/10/20 MHz at 1625–2500 | Doodle Labs Helix / miniOEM / nanoOEM → unresolved | low |
| LTE-like 18 / 9 MHz, PSS 25/29/34 (Herelink) | Herelink v1.0 / v1.1 / Blue → unresolved (same Pinecone S1) (https://www.getfpv.com/cubepilot-herelink-blue-bundle-v1-1.html) | – |
| SIYI 5 GHz link | HM30 / MK15 / MK32 air units → unresolved; MK15E / MK32E are 2.4 GHz only | low |
| 2K COFDM 2 / 2.5 / 4 / 8 MHz, 300–1400 MHz, no DVB-T TPS | generic Chinese COFDM modules (many rebrands, identical) (https://www.asmag.com/product_detail.aspx?id=35994) | med |
| DVB-T TPS decodes (D5) | DVB-T-compliant broadcast-class TX (Vislink / IMT class [U]) | med |
| COFDM 1.25 / 1.5 / 1.75 MHz mesh | DTC MeshUltra (https://www.unmannedsystemstechnology.com/2021/07/dtc-unveils-new-software-defined-radio/) | low |

**Identical on air:** all 100 mW SiK boards; RFD900x vs ux vs x-US (region band edges only); Microhard
pMDDL2450 vs -SWP / -ENC; Doodle Labs form factors; Herelink versions; SIYI HM30 / MK15 / MK32.

## 6. Priorities for the library (what pays off)
1. **Decoded identity first:**
   - DroneID product-type table (§0.1) → `model` (high);
   - Remote ID MFR code → `manufacturer`;
   - serial range → `model_candidates`. Needs the FAA DOC list (https://uasdoc.faa.gov/listdocs) and the ICAO register.
2. **Generation by raster / numerology** (already measured):
   - O3 vs O4 by channel centre (5804.5 / 5839.5 vs 5789.5 / 5814.5);
   - 60 MHz or 5.1 GHz ⇒ O4;
   - Avata = 5.8-only O3 downlink.
3. **`tx_antennas`:**
   - O4 Lite vs Pro;
   - Walksnail 1-antenna vs 2-antenna kits;
   - Vista vs FPV Air Unit;
   - Neo / Flip at 2.4 GHz = 1T.
4. **ELRS chip class** from the rate mode; domain from the channel count.
5. **HOP timing (D3)** for DJI models within one generation. Only one published unit per model: low confidence until
   bench-verified.
6. Power ceilings, video modes, firmware menus: **not usable** passively (range unknown, payload-only).

## 7. Gaps / to fetch
- ICAO MFR code register (403).
- FAA DOC serial ranges.
- DroneID product codes for Mini 3, Mini 4 Pro, Air 3, Avata / Avata 2 and Neo (need real decodes).
- 8 vs 9 symbols per model (conflict in §0.1).
- FCC test reports (occupied BW, dwell) for Walksnail (Caddx), HDZero, Tracer, AFHDS 3, ACCESS, FASST.
- Placement of the ELRS Gemini second copy.
- Crossfire 50 Hz LoRa parameters.
- PDFs not readable via fetch (read with pdftotext):
  - RTC6705 datasheet, TBS Pro32 manual, Microhard brochures, RFD900x datasheet (URLs above);
  - Department 13 DroneID paper https://petapixel.com/assets/uploads/2022/08/Anatomy-of-DJI-Drone-ID-Implementation1.pdf.

## Licences
- **Facts only, oracle only (GPL / AGPL):** RUB-SysSec (AGPL), Kismet (GPL-2+), ExpressLRS, Betaflight,
  Multiprotocol, DeviationTX, mLRS (GPL).
- **Usable in C++:** opendroneid-core-c (Apache-2.0), proto17 (MIT), anarkiwi (Apache-2.0).
- Lookup tables of facts (codes → names, rate tables) are re-typed with a source citation, never copied as code.
