# 13 — Model-level gaps: Remote ID manufacturer codes, FAA DOC serial ranges, DroneID product codes, 8 vs 9 symbols, Walksnail header

Research agent, 2026-10-22. Web only; no local captures used. This file closes most of file 11 §7 and PLAN.md M7 / T8b /
D2, and does not repeat file 09 or file 11.

**Tags:** [spec] standard; [vendor] the manufacturer's own data (incl. what it files with the FAA); [code] source code
(GPL / AGPL = facts only); [comm] community / paper; [U] unverified or our inference.

## 0. Two new primary sources (fetch these again to refresh the tables)
1. **FAA UAS Declaration of Compliance (DOC), public JSON API** behind https://uasdoc.faa.gov/listdocs (the page is a
   React app; the API needs no login) [vendor, filed with the FAA]:
   - list: `https://uasdoc.faa.gov/api/v1/publicDOCRev?itemsPerPage=1000&pageIndex=0`. On 2026-10-22 it held 509
     public DOCs: 480 Remote ID (`docType=rid`, 390 "Unmanned Aircraft", 89 "Broadcast Module") and 29 operations
     over people.
   - one DOC: `https://uasdoc.faa.gov/api/v1/publicDOCRev/<RIDnnnnnnnnn>`. The field `publicRevision` carries
     `makeName`, `modelName`, `fccIdentifier`, `serialNumberStart`, `serialNumberEnd`, `docRevId`.
   - serials of one DOC: `https://uasdoc.faa.gov/api/v1/serialNumbers?findBy=docTrackingNumber&docTrackingNumber=<RID…>&docRevId=<rev>&itemsPerPage=1000&pageIndex=<n>`.
     Each item is either `type: "RANGE"` (`value` = `start-end`) or `type: "ANSI/CTA-2063-A"`, a single serial
     split into `mfrCode`, `lengthCode` and `mfrSerialNumber`.
   - The DOC app hard-codes the CTA-2063-A length code: `1`–`9` = 1–9 characters, `A`–`F` = 10–15 characters.
2. **DJI's own AeroScope type map**: `https://mydjiflight.dji.com/links/links/areoscope_type`. It redirects to a
   djicdn file named `areoscope_type_20260920`, a JSON `typeMap` {code: model}, 87 entries [vendor]. Bender cites
   this link as the source of the "AeroScope ID (value in the model field)" (arXiv 2207.10795 Table IV, ref [4]:
   https://arxiv.org/pdf/2207.10795).

## 1. CTA-2063-A manufacturer codes (beyond file 11)
- **ICAO register:** still not readable.
  - https://www.icao.int/operational-safety/doc-8643-aircraft-type-designators/manufacturers-codes returns 403,
    also with a browser user agent. The older `www2023.icao.int` host no longer resolves.
  - The AESA page that links to it gives no list
    (https://www.seguridadaerea.gob.es/en/ambitos/drones/requisitos-de-uas-drones/listado-de-fabricantes-y-modelos-de-uas).
  - ICAO assigns the 4-character codes on request by e-mail; no public list was found.
- **Substitute:** the codes **used in the serials and ranges each maker declared to the FAA** (DOC API above,
  2026-10-22). These are [vendor] self-declarations, not the ICAO register.
  - An OEM can declare another company's code. Examples: A2Z lists 1596 (Dronetag) serials; Camflite lists 1872
    (Veeniix) serials.
  - So a code names the **RID module or serial issuer**, not always the airframe maker.
- **Resolves the 1809 conflict:** both the "Aerobits" and the "idME" DOCs declare 1809. idME is Aerobits' Remote ID
  module line (https://www.unmannedsystemstechnology.com/case-study-demonstrating-the-remote-id-concept-using-the-aerobits-idme/),
  so **1809 = Aerobits**.
- **Unclear names:** the makes "UAS" (1774), "Broadcast" (1935) and "Broadcast Module" (2051) appear in the list as
  printed. Keep the code, report the maker as unknown.

`BM` = the DOC is for a broadcast (add-on) module.

| Code | Maker (DOC `makeName`) | Code | Maker | Code | Maker |
|---|---|---|---|---|---|
| 1595 | Sunflower Labs | 1801 | Drone Volt | 1906 | Aurelia (BM) |
| 1605, 1705 | Pyka | 1807 | Event38 | 1907 | **Anzu Robotics** (Raptor) |
| 1607 | Wingcopter | 1808 | Draganfly | 1908 | MicroFlight (BM) |
| 1620 | Aero Systems West | 1810 | Matternet | 1911 | Cider Copter |
| 1633, 1663 | Wing Aviation | 1812 | Ascent AeroSystems | 1914 | BRINC |
| 1634 | ACSL | 1813 | Harris Aerial | 1922 | TRIPLEFINE |
| 1641, AA23 | Applied Aeronautics | 1816 | Field / Explorer | 1924 | Arcsky |
| 1644 | Microdrones | 1826 | make "ELIOS 3 UAV" (Flyability [U]) | 1925 | DIY-D (BM) |
| 1656 | Quantum-Systems | 1836 | Coretronic Intelligent Robotics | 1926 | Skyfront |
| 1698 | RigiTech | 1839 | Teal Drones | 1928 | XM2 Labs |
| 1701 | Percepto | 1840 | Voliro | 1933 | Triad Drones |
| 1718 | Censys | 1843 | blueflite | 1934 | "C100" (make field) |
| 1720 | Leica Geosystems | 1860 | XiL | 1937 | PrecisionVision |
| 1728 | **Sony** (Airpeak S1) | 1862 | DeltaQuad | 1939 | SPECTA |
| 1744 | Wingtra | 1866 | SmartDrone | 1943 | LISHI |
| 1759 | ABZ Innovation | 1867 | Hoverfly | 1947 | Swoop Aero |
| 1767 | **Teledyne FLIR** | **1869** | **Ruko** (BM too) | 1949 | SKYFISH |
| 1774 | "UAS" (BM) | 1871 | Camflite | 1961 | Phoenix UAS (BM) |
| 1785 | Watts Innovations | 1872 | Veeniix (BM) | 1968 | Light Quadcopter UAV |
| 1786 | Amazon Prime Air | 1873 | Zipline | 1969 | SkySpecs |
| 1787 | BlueMark (BM); also Avidrone, Commaris | 1877 | Lucid Bots | **1971** | **FIMI** (X8 SE 2022, X8 Tele) |
| 1788 | Vision Aerial | 1879 | Flightwave | 1984 | Nokia |
| 1790 | Hylio | 1884 | Workhorse Aero | 1990 | FACEGLE |
| 1791 | Inspired Flight | 1886 | Flytrex | 1995 | Blackbird |
| 1795 | Zephyr Systems | 1889 | WISPR | 2003 | Holyton (BM too) |
| 1796 | Drone Defence (BM) | 1891 | Zing Drone Solutions (BM) | 2009 | Cinnamon Drones (BM) |
| 1798 | Sentera (BM) | 1893 | Pierce Aerospace (BM) | 2021 | SHANGHENG |
| 1809 | **Aerobits / idME** (BM) | 1894 | Drone-Clone Xperts | 2030 | Telegrid |
| | | 1895 | EXO | 2046 | Aurora Flight Sciences |
| | | 1898 | Fotokite | 2047 | Airobotics |
| | | 1900 | Tucok | 2049 | Aeroo |
| | | 1901 | Meteodrone | 2054 | Global Drone |
| | | 1904 | ShenZhen C-Fly | 2055 | Gleesfun (BM) |
| | | | | 2061 | Aerostar Dynamics |

More codes from the same data:
- 2068 Black Swift; 2078 CHCNAV; **2095 Zero Zero Robotics** (HOVERAir brand [U]); 2098 XPECE; 2105 Flyby Robotics;
  2108 Trace / Vesper UAS; 2110 DMR Technologies; **2111 TopXGun**; 2123 Dunlevy Consulting; 2124 Avol;
  2125 Exedy Drones; 2126 Caelone; 2129 Heilongjiang Huida; 2132 Angel Aerial Systems; **2137 Antigravity**;
  2146 Contixo; 2148 Sustainable Skylines; 2157 Ceres Air; 2160 Rainmaker; 2161 Tofect; 2164 Terra Drone;
  2165 / A2ZD A2Z; 2183 SubUAS; 2187 Amor Fati (Seneca); 2199 Voltair; 2218 Manna.
- Still **not found**: Yuneec, BetaFPV, EHang, Hubsan, Eachine. They have no public FAA RID DOC.

## 2. Model from serial: FAA DOC ranges
Source for every row: the DOC API (§0), URL `https://uasdoc.faa.gov/api/v1/publicDOCRev/<DOC>`, fetched 2026-10-22
[vendor]. Rows in **bold** were not in file 11.

### 2.1 DJI: every range is `1581F` + a 3-character model code + 12 characters
Each DJI range is `1581Fxxx000000000000` – `1581FxxxZZZZZZZZZZZZ` (older revisions end `FFFF…`). **The 3 characters after
`1581F` are the model.** They equal the old 14-digit flight-controller serial prefixes (e.g. Mavic 2 Zoom 0M6, Mavic 2
Pro 163, Air 2 1WN/3N3, FPV 37Q, Air 2S 3YT, M300 1ZN: Bender Table IV, https://arxiv.org/pdf/2207.10795) [vendor + comm].

| Model | Prefix | FCC ID | DOC |
|---|---|---|---|
| Phantom 4 RTK | **`1581F0UY`** | SS3-WM334R1801 | RID000000752 |
| Phantom 4 Pro V2.0 | **`1581F11V`** | SS3-WM331S1801 | RID000000759 |
| Mavic 2 Zoom | **`1581F0M6`** | SS3-L1Z1805 | RID000000758 |
| Mavic 2 Pro | **`1581F163`** | SS3-L1P1805 | RID000000757 |
| Mavic 2 Enterprise Advanced | **`1581F4GC`** | SS3-L1ZE1807 | RID000000753 |
| Mavic Air 2 (MA2UE1N / MA2UE3W) | **`1581F1WN` / `1581F3N3`** | SS3-MA2UE1N1911 | RID000000754 / 755 |
| DJI FPV | **`1581F37Q`** | SS3-FD1W4K2006 | RID000000756 |
| Air 2S | **`1581F3YT`** | SS3-DA2SUE12011 | RID000000007 |
| Mini 3 Pro | `1581F4XF` **and `1581F6L5`** (both in one DOC) | SS3-MT3M3VD | RID000000001 |
| Mini 3 | **`1581F5YH`** | SS3-MT3PD22 | RID000000062 |
| Mini 4 Pro | **`1581F6Z9`** | SS3-MT4MFVD23 | RID000000713 |
| Mini 5 Pro | **`1581F9DE`** | SS3-MT5MFND25 | RID000002581 |
| Air 3 | `1581F6N8` | SS3-EB3WBC2303 | RID000000714 |
| Air 3S | **`1581F895`** | SS3-CZ3SCL24 | RID000002046 |
| Avata | **`1581F4QW`** | SS3-QF2W4K22 | RID000000008 |
| Avata 2 | **`1581F6W8`** | SS3-QF3W4K23 | RID000001376 |
| Avata 360 (DVN3XT / DVN3NT) | **`1581FA8J` / `1581FBLK`** | SS3-DVN3NT | RID000002651 / 2650 |
| Mavic 3 / Mavic 3 Cine | `1581F45T` / **`1581F4QZ`** | SS3-L2P2104 | RID000000009 / 010 |
| Mavic 3 Classic | `1581F67P` | SS3-L2P2104 | RID000000041 |
| Mavic 3 Pro / Pro Cine | `1581F67Q` / **`1581F6MK`** | SS3-L2ES2212 | RID000000204 / 205 |
| Mavic 4 Pro (64 GB / 512 GB) | **`1581F8LQ` / `1581F986`** | SS3-L3AB2410 | RID000002375 / 2376 |
| Inspire 3 | **`1581F4Z4`** | SS3-T7402206 | RID000000192 |
| Mavic 3 Enterprise M3E / M3T / M3M | **`1581F5FH` / `5FJ` / `5FK`** | SS3-M3E2206 | RID000000019 / 020 / 057 |
| M3TA | **`1581FBEX`** | SS3-M3E2206 | RID000002584 |
| Matrice 3D / 3TD | **`1581F6QA` / `1581F6Q8`** | SS3-M3D2308 | RID000001166 / 1167 |
| M30 / M30T | **`1581F4BN` / `1581F5BK`** | SS3-M302110 | RID000000011 / 012 |
| M30 / M30T Dock version | **`1581F5BL` / `1581F5BM`** | SS3-M302110 | RID000000273 / 274 |
| M30T (EU & KR) | **`1581F62H`** | SS3-M302110 | RID000001136 |
| M300 RTK | **`1581F1ZN`** | SS3-M3001910 | RID000000021 |
| M350 RTK | **`1581F6GK`** | SS3-M3502301 | RID000000256 |
| Matrice 4E / 4T | **`1581F7FV` / `1581F7K3`** | SS3-JV312406 | RID000002198 / 2199 |
| Matrice 4D / 4TD | **`1581F8HH` / `1581F8HG`** | SS3-M4D2411 | RID000002246 / 2247 |
| Matrice 400 | **`1581F8DB`** | SS3-M4002412 | RID000002486 |
| FlyCart 30 / FlyCart 100 | **`1581F7C6` / `1581FAN4`** | SS3-E2MTR2302 / SS3-T100A2411 | RID000000760 / 2620 |
| Agras T10 **and** T30 | **`1581F446`** (shared) | SS3-T10A / SS3-T30A | RID000001163 / 038 |
| Agras T20P **and** T40 | **`1581F574`** (shared) | SS3-T40A2112 | RID000001181 / 039 |
| Agras T50 & T25 | **`1581F6BU`** | SS3-T50A2303 | RID000001594 |
| Agras T25P / T70P / T100 | **`1581F836` / `8ZX` / `8ZL`** | SS3-T60A2404 / T702412 / T100A2411 | RID000002587 / 2579 / 2576 |

- **Shared prefixes:** T10 / T30 (`446`) and T20P / T40 (`574`). Report both as candidates.
- **No DOC:** Mini 2, Mini 4K, Mini 2 SE, Mini SE and Neo. They are under 250 g, so no FAA Remote ID is required. Their
  serial is not a model key; use the DroneID product type (§3) instead.
- File 11's `1581F45T` (Mavic 3), `67P`, `67Q`, `4XF` / `6L5` and `6N8` are now confirmed [vendor].

### 2.2 Other makers (prefix = range start/end common prefix; "list" = common prefix of the first and last 1000 listed serials)

| Maker | Model → serial prefix | Kind | DOC |
|---|---|---|---|
| Autel | EVO II → `1748CH…`, e.g. `1748CHL78…` – `1748CHN29…` | range | RID000000030 |
| Autel | EVO II V3 → `1748FEV2H` | range | RID000001969 |
| Autel | EVO Lite **and** Lite+ → `1748CLTC` (shared, also MDXM / MDXM2) | range | RID000000032, 035, 026, 027 |
| Autel | EVO Lite (new) → `1748FEC1LTC`, `1748FEC2LTC`, `1748FEV5LTE` | range | RID000001086 |
| Autel | EVO Nano → `1748CNA…`; newer `1748FEC3NAC`, `1748FEC4NAC` | list / range | RID000000320, 1038 |
| Autel | EVO Nano+ → `1748…` (no common sub-prefix in the list) | list | RID000002429 |
| Autel | EVO Max → `1748FEV3HMA` | range | RID000000151 |
| Autel | MDX-1 → `1748FEV3HMM`, `…HMP`, `…HMQ` | range | RID000002344 |
| Autel | Alpha → `1748FEV4AL` | range | RID000001332 |
| Autel | Dragonfish Std / Lite / Pro → `1748CJD1` / `CJD2` / `CJD3`; newer `1748FDF1JD1` / `FDF2JD2` / `FDF3JD3` | range | RID000000048 / 046 / 047 |
| Parrot | ANAFI Ai → `1588E040452` | range (`…AA000001`–`…ZZ999999`) | RID000000067 |
| Parrot | ANAFI USA → `1588E040505` | range | RID000000068 |
| Parrot | ANAFI UKR GOV → `1588E040465` | range | RID000002663 |
| Skydio | X2E → `1668BE10JA00` | range | RID000000097 / 098 |
| Skydio | Skydio 2+ → `1668BS20` (`FB…`–`JB…`) | range | RID000000099 / 100 |
| Skydio | `X10E SR47P` (X10 family [U]) → `1668BR40` | list | RID000000749 / 750 |
| Skydio | R10 → `1668BR70` | list | RID000002670 |
| Skydio | Skydio 2 → `1668…` (no sub-prefix) | list | RID000000300 |
| Anzu | Raptor → `1907F7USE23B…`; Raptor T → `1907F7UEE23B…` | list | RID000001035 / 1036 |
| Potensic | ATOM 2 → `1910F14…` (also a `1910F…` catch-all); RID-916 module → `1910F916` | range | RID000002237 / 1125 |
| Holy Stone | HS600 `1865AHS600`, HS600D `1865CHS600D`, HS720G `1865CHS720G0`, HS920 `1865CHS920`; HS720 / 720E / 720R / 700E → `1865C00000[0-3]…` | list / range | RID000000327, 1175, 0411, 1150, 0472 / 462 / 402 / 458 |
| Ruko | F11 `18699F11`, F11PRO `1869AF11P`, F11PRO 2 `1869BF11P2`, F11GIM `1869BF11GM`, F11GIM2 `1869CF11GM2`, U11PRO `1869AU11P`, U11S `1869AU11S` | list | RID000000287, 288, 683, 521, 251, 286, 285 |
| SJRC | F11 4K PRO `1892AF110`, F11s 4K PRO `1892CF11S23`, F22 4K PRO `1892CF22A23`, F22s 4K PRO `1892CF22B23`, F5s PRO `1892CF5SA00`, F7 4K PRO `1892BF7A238` | list | RID000000414, 394, 492, 379, 473, 413 |
| Bwine | F7GB2 `1875BF7GB2`, F7GB3 `1875CF7GB3`, F7GIM2 `1875CF7GM2` | list / range | RID000000418, 2596, 1711 |
| Zero Zero Robotics | ZZ-H-1-003 `2095FH05`, ZZ-H-1-004 `2095FH04` | range | RID000002593 / 2594 |
| Sony | Airpeak S1 (ARS-S1) → `17287200…` | range + list | RID000000017 |
| Teledyne FLIR | SkyRanger / SkyRaider R70 / R80D → `1767CSRH` | list | RID000000191 |
| SwellPro | FD2 `1930FD02`, FD3 `1930FD03`, A-RID module `1930FARD` | range | RID000001501, 1762, 1321 |
| XAG | P100 `1863FPEG8218`, P100 Pro `1863FPEG8219`, P150 `1863FPEG8223`, P150 Max `1863F251`, P60 `1863F242` | range | RID000000162, 163, 2240, 2602, 2588 |
| Wing | Hummingbird → `16335…` (also `16635…`) | range / list | RID000001459, 0488 |
| Zipline / Amazon | P2 Zip, Sparrow `1873…`; MK27-2 `17868SN000…`, MK30 `17865…` | list / range | RID000002128, 0397, 0378, 1843 |

- **Caveat:** a "list" prefix is inferred from the first and last 1000 serials only [U]. Store the DOC number with each
  table row, so the tables can be rebuilt by script from the API.

## 3. DJI DroneID product-type codes beyond 70
- **Source:** DJI AeroScope `typeMap` (§0.2), [vendor, DJI-published data].
  - Codes 1–70 in this map equal the RUB / Kismet table in file 11 §0.1, code for code. That table includes the
    real-decode check Mini 2 = 63.
  - Bender (Table IV) states that the AeroScope ID is the value of the DroneID model field.
  - So the codes above 70 can be treated as DroneID `device_type` values: [vendor], **not yet confirmed by a real
    OcuSync decode** for each model.
- **Conflict settled:** the llorencroma spoofer table (`mini_3_pro` = 0x13, `mini_2` = 0x0e; repo has no licence,
  https://raw.githubusercontent.com/llorencroma/DJIDroneIDspoofer/master/Drone.py) contradicts this map. Ignore it.

| Code | Model | Code | Model | Code | Model |
|---|---|---|---|---|---|
| 72 | Agras T40 | 89 | M350 RTK | 104 | **Neo** |
| 73 | **Mini 3 Pro** (also Bender Table IV) | 90 | **Air 3** | 105 | **Air 3S** |
| 75 | **Avata** | 91 | Matrice 3D / 3TD | 108 | Mavic 4 Pro |
| 76 | Inspire 3 | 93 | **Mini 4 Pro** | 109 | Agras T100 |
| 77 | Mavic 3 Enterprise E / T / M | 94 | **Avata 2** | 112 | **Mini 5 Pro** |
| 78 | FlyCart 30 | 95 | Agras T60 | 124 | **Neo 2** |
| 82 | Agras T25 | 96 | Agras T25P | 127 | Avata 360 |
| 83 | Agras T50 | 99 | Matrice 4 T / E | 128 | Agras T100S / T70S |
| 84 | **Mavic 3 Pro** | 100 | Matrice 4D / 4TD | 130 | Matrice 5TD |
| 86 | **Mavic 3 Classic** | 102 | **Flip** | 132 | Agras T55 |
| 87 | **Mini 3** | 103 | M400 | 133 / 140 | Lito X1 / Lito 1 |
| 88 | **Mini 2 SE** | | | 166 / 167 / 173 | Agras T200 / FlyCart 200 / **Mini 4K** |

- **Still unassigned:** 71, 74, 79–81, 85, 92, 97, 98, 101, 106, 107, 110, 111, 113–123, 125, 126, 129, 131,
  134–139, 141–165, 168–172, > 173.
  - Mavic 3 Cine has no separate code (it sits under 68 Mavic 3) [U].
  - Bender's code 240 = Yuneec H480 (2022 map) is gone from the 2026 map.
- **Encrypted payloads:** for O4 models (Mini 4 Pro with newer firmware, Air 3, Neo, Flip, Mini 5 Pro …) the code
  is visible only if the payload is unencrypted (file 11 §0.1). The Inspire 3 (O3) is reported as encrypted too
  (README text of https://github.com/luyii-code-1/dji-ocusync-droneid-research, GPL-3 repo, facts only) [comm].
- **Wi-Fi DroneID:** the vendor IE carries the same model field (Kismet parser; Bender §II-B). That covers Mini SE,
  Spark, Mavic Air, and possibly Neo / Flip over Wi-Fi [U].

## 4. DroneID 8 vs 9 symbols, and the legacy payload
- **Resolved (NDSS 2023 §III-B, p. 6, https://www.ndss-symposium.org/wp-content/uploads/2023/02/ndss2023_f217_paper.pdf)** [comm]:
  - Quote: "some drones (Mavic 2 and older OcuSync drones) do not send the first symbol (Symbol 1 …), which results
    in a shorter frame duration of 576 µs. The other parameters remain unchanged."
  - The paper examined Mini 2, Mavic Air 2 and Mavic 2 Pro / Zoom, and states that all three use OcuSync 2.0.
  - **So 8 vs 9 is not the OcuSync version:** Mavic 2 (OcuSync 2.0) sends 8, Mini 2 / Air 2 (OcuSync 2.0, S1 SoC)
    send 9.
  - File 11's reading "Mini 2 and older omit symbol 1" is wrong.
- **RUB-SysSec receiver** (AGPL-3, facts only; https://github.com/RUB-SysSec/DroneSecurity, `src/droneid_receiver_offline.py`):
  the `--legacy` flag is described as "Support of legacy drones (Mavic Pro, Mavic 2)". **8-symbol = Mavic Pro
  (OcuSync 1) and Mavic 2 family** [code].
- **Bender** captured a Mavic Pro and notes "Sometimes, drone ID signals contain eight OFDM symbols". He gives no
  per-model split (https://arxiv.org/pdf/2207.10795 §III) [comm].
- **9-symbol, confirmed:** Mini 2, Mavic Air 2 (NDSS). Air 2S, Mavic 3 and Mini 3 Pro decode with 9-symbol tools
  (antsdr, proto17), so they are 9-symbol [U]. O4 keeps a 9-symbol shell for some encrypted packets; a 10-symbol
  variant has been discussed publicly (luyii repo README, GPL-3) [comm / U].
- **Library rule:**
  - `droneid_frame = 8` → candidates {Mavic Pro, Mavic Pro Platinum, Mavic 2 Pro / Zoom / Enterprise / Enterprise
    Advanced}. The Enterprise and Platinum members are [U]: same platform.
  - `droneid_frame = 9` → newer.
  - A decoded product type still wins.
- **Legacy (8-symbol) payload layout:** **no separate layout is published.**
  - RUB's decoder runs the same chain on both variants. In legacy mode it takes the bits of all 6 data symbols
    instead of dropping symbol 1. The rest is unchanged: Gold descrambler seed 0x12345678, Nc 1600, systematic
    stream from offset 4148, length 1412 (`src/qpsk.py`) [code].
  - RUB parses the result with **one** struct for every drone (`droneid_packet.py`). That is the 0x5B-byte record
    of NDSS Fig. 13:
    - `len_pack`@0, `zero_byte`@1, `version`@2, `sequence_num`@3 (u16), `state_info`@5 (u16);
    - `serial_num`@7 (16 B), `longitude`@0x17, `latitude`@0x1B (i32), `altitude`@0x1F, `height`@0x21;
    - `v_north` / `v_east` / `v_up`@0x23 / 25 / 27, `yaw`@0x29;
    - `gps_time`@0x2B (u64), `rc_latitude`@0x33, `rc_longitude`@0x37, `home_longitude`@0x3B, `home_latitude`@0x3F;
    - `device_type`@0x43, `uuid_len`@0x44, `uuid`@0x45 (20 B), `crc`@0x59 (u16).
  - So the 8-symbol frame carries the same 91-byte record with the same CRC16, as far as RUB's decodes of Mavic 2
    show [comm].
  - **Not published:** whether a Mavic Pro / Mavic 2 sends the version-1 record (pitch / roll instead of GPS time
    and pilot position: Bender Fig. 9, D13 white paper https://petapixel.com/assets/uploads/2022/08/Anatomy-of-DJI-Drone-ID-Implementation1.pdf).
    Parse it by the `version` byte [U]. T8b stays "layout = v2 struct, verify on a real Mavic 2".

## 5. Walksnail Avatar header (D2)
- **No non-GPL source found.** The only description is in two issues of the GPL-3 project v0l/waveshark, opened by
  the author on 2026-09-27 (https://github.com/v0l/waveshark/issues/227,
  https://github.com/v0l/waveshark/issues/226) [comm, from a GPL project; the issue text states facts, the code stays
  oracle only].
- **High-level facts** (issue text):
  - **Header:** OFDM symbols 0 and 1 of each downlink frame are BPSK, an independent stream per TX antenna.
  - Each stream is a linear code of rank about 13–16 across frames. Its span contains a **6-bit frame counter,
    period 64 frames = 458 ms** (64 × 7.163 ms).
  - The counter is **not plain binary**: "bits 0–3 count directly; bits 4 and 5 sit at other phases".
  - The remaining rank is unidentified (MCS / length suspected).
  - **The exact bit positions and equations are not in the issue text.** They are only in the GPL commit e2cd400
    (`decode::walksnail::Link`), so do not use them. Derive them on the bench instead: fold header bits over 64
    frames and solve the GF(2) rank.
- **Waveform (issue #226):**
  - 802.16e OFDMA 20 MHz profile: 22.4 MS/s, FFT 2048, CP 1/8, 1729 used carriers.
  - AMC bins of 8 data + 1 pilot carriers; STBC Matrix A over symbol pairs (fits our `tx_antennas`).
  - Pilot PRBS as in DVB-T at `k + 1745 − (s mod 32)`.
  - 802.16e LDPC (rates 1/2 … 5/6, n = 2304); four 12-symbol time-interleaver blocks; no per-frame scrambling.
  - The interleaver permutation is unknown, so the payload is not decodable.
- **Goggles uplink (issue #227):**
  - A 938 µs burst every 7.163 ms frame: 137 µs of unknown content, the downlink preamble, then 6 BPSK symbols on
    carriers +6 … +438 (about 4.7 MHz).
  - It is sent even with no VTX answering, so it **reveals a pilot's goggles**. The paired-VTX ID bits are unknown.
  - New detect-only feature: a pilot on the channel.

## 6. What changes in the library plan
- **M2 / M7:**
  - Add a `serial_prefix → model` table, generated by a script from the DOC API (§0) and stored with the DOC number:
    DJI 8-character prefixes; Autel / Skydio / Parrot as in §2.2.
  - Extend the manufacturer table with §1.
  - A prefix shared by two models (DJI T10 / T30, T20P / T40; Autel EVO Lite / Lite+) → `model_candidates`.
- **M1:** extend the DroneID product-type table with codes 72–173 (§3), with the source tag `vendor-aeroscope-map`.
  Codes 1–70 are unchanged.
- **D4:** `droneid_frame = 8` → Mavic Pro / Mavic 2 family candidates (§4). Close the conflict in file 11 §0.1.
- **D2:** the header facts are enough for a bench derivation (GF(2) rank over 64 frames). Nothing can be copied.
  Add a goggles-uplink burst detector (938 µs every 7.163 ms) as a separate item.

## 7. Bench check on real recordings (2026-10-22, user-approved, read-only)
Recordings: `../droneid-decoder-reference/DJI_Mini_Pro_4_drone_id/` and `data/drone_id/` (65 captures, 2 M samples each),
library `dji.detect` with DroneID decode, plus the external decoder's legacy mode for the 8-symbol bursts.
- **Mini 2:** 22 bursts in 22 captures decode (CRC ok), product type Mini 2 (code 63), one aircraft (serial prefix
  3NZCK1H0), 9-symbol frames (643 µs). Re-confirms code 63.
- **DJI Mini 4 Pro (folder label):** 33 bursts in 14 captures, all **8-symbol** (572 µs, first-symbol CP absent:
  0.03–0.10 vs 0.93 on the Mini 2). They do not decode in standard or legacy mode (outer CRC-24 fails): the O4 payload
  is encrypted or scrambled differently. So **code 93 (Mini 4 Pro) is not confirmed**, and the 8-symbol frame is
  **not** exclusive to the Mavic Pro / Mavic 2 family (§4): the library's candidate list now includes the Mini 4 Pro.
- **No code above 70 could be confirmed** from these recordings.
