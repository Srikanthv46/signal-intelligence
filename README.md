# spectrum_intel

Spectrum intelligence for IQ blocks (ADALM-Pluto and other SDRs, 300 MHz – 6 GHz, India). For every block the
pipeline:
- **detects** every signal from the IQ: bursty or continuous, narrowband or wideband;
- **localises** each one in absolute time and frequency;
- **identifies** it from the India band plan, confirmed by measured features;
- **decodes** it where a receiver exists;
- **explains** how each result was obtained, in the JSON report.

The library is modular. Each signal family is a self-contained **receiver** package that depends only on a small
shared `core`, so it can be lifted into another pipeline or removed from this one.

## Repository layout

| Folder | What |
|---|---|
| `python/` | the **Python** library (`spectrum_intel`) and its tests: the development reference, used by the signal-vault GUI |
| `c++/` | the **C++** library: the production build for edge devices, same results as Python (see `c++/README.md`) |
| `schema/` | the contract both share: response JSON Schema, technology taxonomy, receiver requirements |
| `testdata/` | captures with known contents (see `testdata/README.md`) |
| `out/` | run outputs (not source) |

**How the code evolves:** a new signal is detected and decoded on the bench first; once verified it goes into
`python/` (validated on `testdata/`), then into `c++/` (parity-checked against Python on `testdata/`). What is in
`python/` and `c++/` is tested and working; nothing unverified is added to either.

## Quick start

```
.venv/bin/spectrum-intel run <snapshot dir | .sigmf-meta | .sigmf-data | folder | sweep folder> \
        [--out out] [--workers 2] [--resume] [--no-png] \
        [--freq-range 2300:2600,4900:6000] [--enable analog_vtx,wifi] [--disable dji]
.venv/bin/spectrum-intel receivers                # list installed receivers
.venv/bin/spectrum-intel summary out/<sweep>      # rebuild the sweep summary only
```

**Install:** the `python/` folder is installed editable into `.venv` (`pip install -e python --config-settings
editable_mode=compat`), which is radioconda Python 3.10 with `--system-site-packages`, plus `py3gpp`, `pyModeS`,
`jsonschema`, `pytest` and `pyflakes`. `.venv/bin/python -m spectrum_intel.cli` does the same thing.

**Output location:** `<out>/<sweep or folder name>/<snapshot folder>/`, or `<out>/single/...` for a lone file.

| Scope | Files |
|---|---|
| Per snapshot | `<name>_report.json`, `<name>_annotated.png`, `<name>_vtx<MHz>_video.png` (one per demodulated analog VTX) |
| Per sweep | `summary.json`, `signals.csv`, `occupancy_map.png`, `run_log.json`, `_sweep/sweep_context.npz` |

**Python:** `from spectrum_intel.analyze import analyze, analyze_iq` takes a SigMF path or an in-memory block;
`spectrum_intel.report.writer.build(a)` then gives the report dict. For production use the API below.

**Memory:** keep `--workers` at 2 or fewer on an 8 GB machine. Each worker needs about 1 GB.
**Threads:** inside one block, receivers and independent hypotheses run on `SPECTRUM_INTEL_THREADS` threads
(default min(4, cores); 1 = strictly sequential). Results do not depend on the thread count.

## Production API: IQ in, JSON out

The acquisition pipeline hands over IQ samples, the sample rate, the centre frequency and the analog bandwidth;
the library returns one JSON response per block. Python and C++ return the same response
(`schema/response.schema.json`).

```python
from spectrum_intel.api import Engine, IQBlock, StreamAnalyzer

engine = Engine()                                     # loads the receivers once; reuse it
resp = engine.analyze(IQBlock(iq, sample_rate_hz=61.44e6, center_frequency_hz=2.437e9, bandwidth_hz=56e6,
                              start_time="2026-10-09T10:00:00Z", source_id="pluto-1"))
```
```cpp
#include <spectrum_intel/spectrum_intel.hpp>
si::Engine engine;
si::IQBlock b;  b.data = samples;  b.n_samples = n;  b.format = si::IqFormat::CI16;
b.sample_rate_hz = 61.44e6;  b.center_frequency_hz = 2.437e9;  b.bandwidth_hz = 56e6;
std::string json = engine.analyze_json(b);            // or engine.analyze(b) -> nlohmann::ordered_json
```
A C API (`c++/include/spectrum_intel/c_api.h`) serves other languages.

**Input:** complex samples, or interleaved I/Q as int16/int8/uint8/float (`iq_format`), not rescaled (raw ADC
counts keep the clipping check meaningful). The centre frequency is trusted as given. The bandwidth limits the
analysed band to the analog passband (unset = the whole sampled band). Optional: `start_time` (UTC of the first
sample: absolute times in the response), `gain_db`, `source_id`, `request_id`, and **`receiver_ppm`**, the SDR's
frequency error. A sweep calibrates the error itself; a single block cannot always do so. On the 928 MHz testdata
snapshot, the GSM receiver finds 8 BCCH carriers with `receiver_ppm=7.44` (the Pluto's TCXO) and none without it.
It also lets the pipeline flag lines locked to the receiver's own reference (`receiver_spur`).
So pass it when the acquisition side knows it (a factory or LTE/GPS calibration, or a GPS-disciplined clock). The engine never raises for bad input: the
response then has `status` `insufficient_input` (block too short, NaN samples, bad metadata) or `error`, with the
reason in `errors[]`.

**How many samples.** Measured on testdata at 61.44 MSps: below 16,384 samples the analysis cannot run; 16k–64k
samples (0.3–1 ms) give spurious detections; real decodes start near 131,072 (2 ms). So:
- **Minimum:** 65,536 samples. Smaller blocks are refused with `insufficient_input`. Blocks of 2048 or 4096
  samples must be accumulated first: `StreamAnalyzer` (Python and C++) collects contiguous blocks of one tuning and
  analyses each full block.
- **Recommended:** 2^21 = 2,097,152 samples (34 ms at 61.44 MSps). This covers one LTE/NR sync period with margin,
  one analog video field, and many Wi-Fi/BLE bursts. A `SHORT_BLOCK` warning flags shorter blocks.
- **Per receiver** (`schema/requirements.json`; the response's `coverage` block reports this per block):

| Receiver | Min. bandwidth | Min. duration | For full output | Why |
|---|---|---|---|---|
| wifi | 20 MHz | 0.2 ms | 102.4 ms | one PPDU; every AP beacons once per 102.4 ms |
| rc_link (BLE, LoRa, FSK) | 2 MHz | 0.4 ms | 0.5 s | a BLE advert ≤ 376 µs; a LoRa SF12 preamble ~0.3 s |
| dji | 15.36 MHz | 10 ms | 0.7 s | video/RC frames every few ms; DroneID about every 0.6 s |
| analog_vtx | 20 MHz | 5 ms | 40 ms | tens of video lines; a picture needs one field |
| jammer | 20 MHz | 3 ms | 30 ms | a sweep of ≥ 10 MHz; the fold needs ≥ 8 sweep periods (300 Hz–120 kHz sweeps) |
| voice_fm | 50 kHz | 10 ms | 0.4 s | CTCSS ≥ 0.3 s, DCS ≥ 0.35 s |
| digital_vtx (Herelink) | 1.92 MHz | 10 ms | 20 ms | PSS every 5 ms, folded over half-frames |
| digital_vtx (HDZero) | 24 MHz (the whole 23.44 MHz band inside the usable span) | 1 ms | 34 ms | pilot coherence over ≥ 24 symbols; the frame rate needs two per-frame preambles (2 frames at 60 fps) |
| lte | 1.92 MHz | 10 ms | 20 ms | PSS every 5 ms; SSS in subframes 0 and 5 |
| nr | 3.84 MHz (7.68 for 30 kHz SSB) | 20 ms | 40 ms | SSB period 20 ms |
| gsm | 0.5 MHz | 47 ms | 240 ms | FCCH every 46.15 ms; BCCH needs the 51-multiframe |
| adsb | 2 MHz (8 MSps preferred) | 0.12 ms | 10 s | a frame is 64/120 µs; a position needs an even/odd pair |

**Response** (`schema/response.schema.json`, version 1.3.0; 1.3.0 added `validation` (proven / research, see below); 1.2.0 added the `esb_like_link` technology; 1.1.0 added
the `dmr` / `nxdn` technologies, the `mcc` / `mnc` / `colour_code` identifiers, the `pmr_digital` receiver and the
`possible_lo_harmonic_image` flag). Units are SI and named in the key; times are seconds
from the first sample plus UTC when the start time is known:
- `status`, `engine` (implementation, version, receivers), `input` (what was analysed, including the analysed band);
- `quality`: ADC clipping, RMS level, noise floor, centre-frequency label check;
- `coverage`: per receiver, whether this block's band, bandwidth and length can carry what it looks for;
- `summary`: counts by decode status and family, technologies seen;
- `signals[]`, each with:
  - a stable `technology` code and `family` (`schema/technologies.json`: `lte`, `nr`, `gsm`, `wifi_ofdm`, `ble`,
    `lora`, `dji_droneid`, `analog_fpv_video`, `mode_s_adsb`, ... and `unknown`) and `category`
    (`emission`/`artifact`);
  - the human `label`, `confidence` and `inferred` (a "-like" rule label);
  - `description` (since schema 1.4.0): two or three plain sentences saying what the signal is and what the device does
    (drones are named only where the signal is drone-specific), e.g. "Video transmitter of an FPV drone (Walksnail Avatar class, ...): it sends
    the camera picture from the drone to the pilot's goggles. An FPV drone is switched on nearby, usually within a few
    kilometres." One text per technology, with label variants where a technology covers different things (Walksnail vs
    PMP 450 fixed wireless, ordinary Wi-Fi vs a Wi-Fi drone link, DJI video vs DroneID), in `schema/descriptions.json`
    (shared by Python and C++; a receiver artifact is described as one). Also in the report (`report.writer`) and in
    the signal-vault Signal Intelligence tab (label tooltip, detail panel);
  - `frequency` (centre/low/high/bandwidth/offset), `time` (start/end/duration, UTC, bursts), `snr_db`, `modulation`;
  - `decode`: status, reason, method, fields and steps;
  - `identifiers`: PCI, BSIC, ICAO, MAC, SSID, serial and so on, when decoded;
  - `evidence`, `band_allocation`, `flags`, and the `detection` and `features` provenance (left out with
    `compact`);
- `weak_detections[]`, `emitters[]`, `devices[]`, `findings[]`, `extensions`;
- `warnings[]` and `errors[]` with machine-readable codes (`SHORT_BLOCK`, `ADC_CLIPPING`,
  `BLOCK_TOO_SHORT_FOR_RECEIVER`, `RECEIVER_ERROR`, ...);
- `timing_ms`.

## Performance

Python, on the 12 SigMF captures in `testdata/` (Apple M1, 4 threads, measured while other jobs ran). Every report
is identical to the previous code's (`processing` timings aside):

| | before | after |
|---|---|---|
| total for 12 captures | 635 s | 339 s |
| NR-band capture (D_20260211-161016, 2060 MHz) | 162 s | 54 s |
| Herelink capture (two passes: label corrected) | 127 s | 73 s |
| typical 5.8 GHz drone capture (DJI Mini 2) | 25 s | 16.5 s |

What changed: the NR receiver caches its 3GPP sequences (py3gpp built every Gold sequence bit by bit in Python;
`core.seq3gpp.gold` is a bit-exact vectorised generator) and correlates its three PSS candidates in one batched FFT.
Receiver packages run concurrently, each on its own copy of the capture cache, so the results do not depend on
scheduling. NR raster hypotheses and the Herelink PSS-fold batches run on threads. The C++ library is the fast
path for edge devices (see `c++/README.md`).

C++ against Python, all receivers enabled, on the 25 testdata inputs. Both run 2 threads; the times are the pipeline
timings in each report, and the reports are identical (`c++/tools/parity.py`, 2026-10-10):

| | Python | C++ | speed-up |
|---|---|---|---|
| total for 25 inputs | 414 s* | 187 s | 2.2x |
| typical 5.8 GHz drone capture (2M samples, DJI Mini 2) | 17.2 s | 8.5 s | 2.0x |
| NR-band capture (D_20260211-161016) | 51.4 s | 22.9 s | 2.2x |
| Herelink capture (D_20251125-110453) | 88.3 s | 40.7 s | 2.2x |
| real 802.11n capture (real_ht_mcs1_5240) | 7.4 s | 3.5 s | 2.1x |
| live 868 MHz block (LTE B5/B8 + DMR/NXDN + ISM) | 14.9 s | 8.6 s | 1.7x |
| live TETRA dwell (390.78 MHz) | 3.6 s | 1.7 s | 2.1x |

\* D_20251119-182450 is counted at 61.5 s, its time when run alone; in the parity run the machine stalled
(1948 s), and C++ took 21.1 s on it. The Herelink capture is slower in Python than before (60 s) because the LTE
search now also tries bursty B40 candidates there (about 20 s) and finds no cell.

Most of the receiver time goes to three receivers: digital_vtx (blind OFDM numerology and the Herelink PSS fold),
nr_ssb (the raster and CFO hypotheses) and rc_link (the LoRa dechirp scan). Build C++ with `-DSI_RECEIVERS=` to keep
only the receivers a deployment needs.


## Validation tags: proven vs research
Every signal carries `validation` (the API signal and `decoded.validation` in the per-snapshot report), from
`schema/capabilities.json`:
- **`proven`**: the capability was tested against a real transmitter or capture. `proof` names it (testdata).
- **`research`**: the capability was built from published specifications (papers, standards, datasheets, open
  source used as reference) and checked on synthetic signals only. It has not yet been tested on a real
  transmitter. `sources` names the research file (`research/drone_transmitters/*.md`).
- **`research_fields`** lists decoded fields that come from research-level capabilities inside an otherwise proven
  receiver (e.g. `tx_antennas` in `digital_vtx`). `research_sources` gives where they came from.

**Promotion:** after testing against a real transmitter, set the entry's `level` to `proven` and add `proof` (capture
name + date) in `schema/capabilities.json`, plus a `testdata/README.md` row. No code changes. A technology missing
from the registry is reported as `research`. `python/tests/test_validation_tag.py` and the C++
`test_contract_validation.cpp` require every technology code to have an entry.

## Model-level identification (`fields.model`, research)
Within a family (DJI, Walksnail, HDZero, ELRS, …) a signal can carry `decoded.fields.model`, built in
`core/drone_models.py` (C++ `core/drone_models`) from `research/drone_transmitters/11_model_level_identification.md`:
- **Shape:** `{family, model, candidates, basis, evidence, unresolved_reason?}`. `model` is set only when exactly one
  model fits; otherwise it is null, `candidates` lists what fits and `unresolved_reason` says why (for example,
  several models share the same RF). An inferred model gets the "-like" suffix; `basis` is `decoded` (identity sent
  in the clear) or `inferred` (from RF features).
- **Decoded identity, DroneID:** the DJI DroneID product-type byte → model.
  - Codes 1–70 come from the public decoder tables and are confirmed by real decodes.
  - Codes up to 173 come from DJI's own AeroScope type map and are not yet confirmed per model.
  - Unknown codes stay unknown.
- **Decoded identity, Remote ID (BLE and Wi-Fi):** the serial → CTA-2063-A manufacturer code → manufacturer
  (`remote_id.manufacturer`). The 142 codes are file 11's plus the FAA DOC declarations. An add-on broadcast module
  names the module maker, not the airframe.
- **Serial → model:** the serial prefix → model, from the FAA DOC ranges: DJI `1581F` + model code for about 55
  models, plus Autel, Parrot, Skydio, Anzu, Potensic, Holy Stone, Ruko, SJRC, XAG and others. Shared prefixes
  (Agras T10 / T30) stay as candidates. Prefixes inferred from listed serials are "-like". DJI models under 250 g
  have no DOC.
- **Inferred:** a 9-symbol DroneID frame gives only the generation. An 8-symbol frame → candidates: the Mavic Pro /
  Mavic 2 family (NDSS 2023) and the DJI Mini 4 Pro, which sends 8-symbol frames whose payload does not decode
  (measured on 14 real recordings).
- **Inferred, digital video (`digital_vtx`):** the best conflict-free research-catalog systems map to their published
  model lists. The DJI O3 / O4 20 / 10 MHz channel centres (O3 5804.5 / 5839.5, O4 5789.5 / 5814.5 MHz) separate the
  generations. 2 TX antennas (STBC / SFBC) remove the single-antenna units: Walksnail Mini 1S / Nano V3 / V2 single /
  GT2 and the O4 Air Unit (Lite). The HDZero 17 MHz mode implies the RunCam Nano 90 camera (not yet measured; the
  27 MHz mode is, on a Freestyle V2). Units with identical RF
  stay unresolved: Walksnail V1 Pro = V2 dual = Moonlight (GT only above ~1.2 W); every HDZero VTX; Vista = RunCam
  Link; Mavic 3 / Classic / Pro; Mini 4 Pro / Air 3 / Air 3S.
- **RC and telemetry (`rc_link/models.py`):** ExpressLRS radio chip from the air mode (sub-GHz SF5 or 2.4 GHz
  300 kb/s GFSK ⇒ LR1121, FLRC ⇒ SX1280, other LoRa modes ⇒ SX127x / SX1280 or LR1121) plus the regulatory domain;
  FrSky variant = firmware generation / region build; the SiK air-rate set (the RFD 3.x rates ⇒ RFD900x / ux / 868x).
  The protocols whose transmitters differ only in power stay unresolved, with the candidates listed: Crossfire
  Micro / Nano / full TX, AFHDS 2A radios, FrSky radios vs multiprotocol modules, 100 mW SiK boards.
- **Analog VTX:** model unresolved (power, SmartAudio / Tramp control and the camera standard are not visible), with
  the channel and whether a microphone audio subcarrier is sent; 5584 MHz = TBS Unify Pro pit mode, 1258 / 1280 MHz =
  locked 2-channel 1.2 / 1.3 GHz VTXs.
- **Uplink / control hop timing (`digital_vtx` blind hops):** hop width, on-time and interval matched to published
  values. **DJI:** FPV, Avata 2, Mavic 3 Pro, Mini 3, Mini 4 Pro, the O4 range and OcuSync 2. **Others:** Herelink,
  Skydroid H12 / T10, Siyi MK15 / MK32 / FT24. One unit per model was published, so confidence is low and a single
  match is labelled "-like". The family names the matching brands. The real Walksnail uplink in testdata matches
  none.
- Research level: tables from published sources, checked on synthetic and table tests only.

## What is detected and decoded

**Status values:**
- **success:** a payload verified by a CRC or FCS, or for analog video, a demodulated picture.
- **partial:** the physical layer is synchronised and its identity or parameters are measured, but the payload is
  not decoded.
- **not_possible:** detected and classified only.

Every signal's `decoded.steps` and `detection` blocks in the JSON record exactly what was done.

| Signal | Receiver | Detected by | Decoded (fields) | Status |
|---|---|---|---|---|
| Wi-Fi 802.11a/g/n/ac/ax (OFDM) | `wifi` | L-STF/L-LTF preamble on the 20 MHz channel plan; **off the grid or outside the Wi-Fi bands** (fixed-wireless backhaul such as Mimosa at 4880 / 5755 / 5790 MHz): a 3.2 µs CP-correlation trigger on 20 MHz slices — whole-capture, or per 250 µs window for weak, low-duty links (PLAN G13: the Mimosa link at ~6.6 dB / 20 % duty in `LIVE_MIMOSA_4880MHz_weak`; run only where the whole-capture ratio ≥ 2, so noise costs nothing) — then blind 1 MHz centres around it (`ofdm_rx.offgrid_plan`) | L-SIG, HT-SIG/VHT-SIG-A; **HE-SIG-A (802.11ax, CRC-4): SU / MU / TB, BSS colour, uplink / downlink, MCS, streams, bandwidth, TXOP; colour → BSSID through a beacon's HE Operation element** (PLAN W2: all 44 real HE PPDUs of `LIVE_WIFI_5755MHz` pass; synthetic 100 % from 8 dB); non-HT MPDU (6–54 Mb/s) → FCS; **802.11n HT MCS 0–7, 20 MHz, BCC or LDPC, with or without STBC** (Alamouti from two antennas) → A-MPDU split on its delimiters, FCS per MPDU; VHT headers of payloads out of reach (80 MHz): direction, streams, MCS, airtime; full MAC dissection (below); **transmitter of every PPDU** (below). 2+ spatial streams and 40/80 MHz payloads are declined with the reason | success (FCS) / partial (payload needs >1 antenna) |
| 802.11 at a non-standard clock, frequency or orientation: 5/10 MHz channels, vendor data links (Microhard 11 MHz, DJI Wi-Fi 5 MHz), custom/off-grid channels (MikroTik), up/down-converted (mirrored) links | `wifi` (`wifi_ofdm_scaled`) | **any band**: bursty bands with sharp edges → blind numerology (CP lag = T_u, symbol line = T_s, guard 1/4 → clock 64/T_u, width fits 52–56 of 64 subcarriers) | the capture is time-scaled into a standard channel and the 802.11 receiver above runs on it, in both spectrum orientations: clock, mirroring, centre, every field above (MPDUs, FCS, MACs, vendor), transmitters | success (FCS) / partial |
| Wi-Fi 802.11b (DSSS/CCK) | `wifi` | Barker-11 phase-coherence detector, SFD-first acquisition at every sample phase | 1/2/5.5/11 Mb/s, long/short preamble: PLCP CRC-16, coherent Barker / equalised CCK PSDU, FCS (with FCS-guided correction and interference excision); full MAC dissection: per-BSS table (SSID, security, PHY, uptime, country, load, WPS), stations, one line per frame | success (FCS) / partial (header, or Barker-only) |
| Drone Remote ID (over Wi-Fi) | `wifi` | vendor element in any decoded 802.11 frame; **Wi-Fi NAN service discovery** (public action, Service Descriptor Attribute with service ID SHA-256("org.opendroneid.remoteid")[:6]; *research*) | **ASTM F3411 / ASD-STAN**: UAS ID, aircraft type, status, position, altitude, height, speed, track, operator ID and **operator position**; **DJI DroneID**: serial, model, position, home point, operator position — located by its CRC-16/Kermit; **China GB 46750-2025** packets in the same element (type byte 0xFF: serial, registration mark, UA class, aircraft / operator positions, altitudes, speed, track, Unix-ms timestamp; *research*, PLAN G14, checked against a published DJI Mini 5 Pro frame) | success (checksum / standard framing) |
| **Drone Remote ID over Bluetooth LE** (ASTM F3411 / EN 4709-002, legacy advertising and **BT5 Long Range**: LE Coded PHY S=8 / S=2, ADV_EXT_IND → AUX_ADV_IND) | `rc_link` + `core/remote_id` | CRC-valid BLE advert whose AD has Service Data UUID 0xFFFA, app code 0x0D | every F3411 message: UAS serial / ID type, UA type, status, drone lat/lon, altitudes, height, speed, track, vertical speed, accuracies, timestamp, **operator lat/lon / altitude**, operator ID, self-ID, area; plus the advert's AD fields (`ble_ad`: name, company, UUIDs, TX power); label "Drone Remote ID (Bluetooth LE, …)", technology `remote_id_ble` | success (CRC) — **research**: synthetic only, full decode ≥ 8 dB per advert (≥ 6 dB for a 376-bit BLE 4 advert) with the coherent GFSK detector + CRC Chase (`ble_coherent`, research); BT5 Long Range (`phy` "LE Coded S=8/S=2", `extended`: AdvA, ADI, AuxPtr) decodes to ~1 dB (synthetic) |
| Devices and findings | `devices` | no IQ: reads the decoded 802.11 frames of every signal | per-device inventory (vendor, role, SSID, activity, airtime, probed networks, clock offset), randomised-MAC linking (probe signature + clock offset, and since PLAN W7 a sequence-number counter continuing across the address change → confidence 0.9), and ranked findings: drones (Remote ID / signature, and **behaviour-only drone video links** from the Wi-Fi signals' `drone_links`), deauthentication floods, duplicate SSIDs, **beacon floods** (≥ 15 unregistered beaconing BSSIDs from ≥ 8 unrelated prefixes; PLAN W6, *research*), weak security | report `devices` / `findings` |
| Bluetooth LE advertising | `rc_link` | 1 Mb/s GFSK bursts on channels 37/38/39 (a BLE-wide 2.4 GHz burst too weak for a symbol rate still gets one coherent attempt) | access address, de-whitening, CRC-24 → PDU type, advertiser address; coherent detector + CRC Chase first (~4–5 dB more sensitive; measured modulation index; *research*), hard slicer as fallback; **capture frequency-label check**. A 1 Mb/s GFSK hop on any other channel is "Bluetooth (BLE/BR)-like" (BLE data channel or BR/EDR hop: the measured modulation index does not separate h 0.5 from 0.32) | success (CRC) / partial |
| **2 Mb/s GFSK packet links of the Enhanced ShockBurst family** (nRF24 / nRF52 ESB and clones: 2.4 GHz mice, keyboards, game controllers, dongles, toy RC) | `rc_link` (`esb.py`) | 2.4 GHz bursts 0.5–3 MHz wide sliced at 2 Mb/s; an alternating preamble of ≥ 14 bits starting in the first 32 bits | the packets of one link (same 16-bit sync word) are checked together: the **CRC-16/CCITT** end and init are found from the packets alone (no per-product table; ≥ 3 distinct packets on one CRC), giving CRC-verified frames (hex), sync word per node, node frame period, carrier offset on the 1 MHz raster. Payload meaning is vendor-specific and not interpreted | success (CRC) / partial |
| LoRa (433/865/890 MHz, any SX127x/SX126x link) | `rc_link` + `lora_phy` | blind dechirp scan (SF7–12 at 125 kHz, SF7–10 at 250 kHz, SF5–9 at 500 kHz) + SFD, in the RC bands and (LoRa only) on 860–930 MHz, where EU868 / FCC915 drone links and off-band telemetry radios sit inside cellular bands; sync symbols must be a sync word's nibbles × 8; a preamble-only report needs its SFD past the end of the capture or seen | explicit-header packets: length, coding rate, sync word, **payload bytes**, CRC-16; LoRaWAN (DevAddr, FCnt, FPort, join EUIs) and Meshtastic (sender, destination, packet id, hops) plain headers by sync word; **SF5 (SX126x) implicit-header packets with a CRC**: length and coding rate found by the CRC, payload bytes (live: 890.03 MHz telemetry link, `LIVE_LORA_SF5_890MHz`) | success (CRC) |
| ExpressLRS | `rc_link` | as above | SF, BW, preamble, payload length, ELRS grid and air-rate mode (live 2.4 GHz "100Hz Full", `LIVE_ELRS_2G4_2440MHz`); *research*: 2.4 GHz binding-phrase UID[5] bit 0 from the IQ polarity (`elrs_identity`); **sub-GHz SYNC packets decoded and self-verified** (implicit-header LoRa checked against gr-lora_sdr, CRC-14 seeded by the packet's own UID4/UID5 → hop index, nonce, rate, TLM ratio, UID bytes; `elrs_sync`, synthetic to 0 dB); *research*: **hop recovery** from a SYNC packet plus the link's other packets in the block (LCG + block-shuffled 256-entry sequence → UID2 / UID3 candidates, the regulatory domain among overlapping grids, the next channels when unique; `elrs_sync.hop_recovery`). SF6 (200 Hz) SYNC packets decode too (SX127x SF6 framing, checked against gr-lora_sdr). 2.4 GHz payloads: the SX1280 long-interleaving coding rate (CR_LI) is undocumented | partial / success (SYNC CRC) |
| mLRS | `rc_link` (`mlrs.py`) | LoRa bursts with preamble 12, normal IQ and 91-byte frames (2.4 GHz SF5–7 at 812.5 kHz; sub-GHz SF5 / SF6 at 500 kHz) | air mode → "mLRS-like" (kept apart from ExpressLRS by preamble, IQ polarity and frame length); the SX127x-framed SF6 mode (19 Hz "7x") is decoded: CRC-16/MCRF4XX over bytes 0–88 (byte order / crc1 coverage tried and named) → the per-binding FrameSyncWord (`mlrs`, model block with the radio chip class) | **research** (synthetic) |
| **TBS Crossfire 150 Hz** | `rc_link` (`crossfire.py`) | 85.1 kBd FSK (±42.48 kHz) in 850–935 MHz | sync word 0x3E6B3E after the preamble → 23 / 13-byte packets (raw + PN9 de-whitened), 260 kHz grid channel and up/down link (`crossfire`); label "TBS Crossfire-like"; the per-hop CRC init is unknown, so packets are not verified | **research** (synthetic, 12–20 dB) |
| **FLARM beacons** (detection + clear header) | `rc_link` (`flarm.py`) | 100 kchip/s, ±50 kHz FSK at 868.2 / 868.4 / 866.0 MHz | Manchester chips → sync F5 31 FA B6 → 24-bit radio ID and address type (FLARM / ICAO) from the clear header, reported only when the header has no Manchester violation (no CRC verifies it); the obfuscated payload is not decoded (`flarm`) | **research** (synthetic: ID 40/40 at 10 dB, 39/40 at 6 dB, never a wrong ID; 0.5 ms) |
| **ADS-L iConspicuity** (EASA SRD-860 M-band) | `rc_link` (`adsl.py`) | FLARM PHY (868.2 / 868.4 MHz, 100 kchip/s Manchester), sync F5 72 4B + length | CRC-24 (Mode-S polynomial), XXTEA zero-key descrambling (checked against a reference implementation) → address table / address, payload type, Traffic: category (11–13 = **UAS**), position, altitude, speed, climb, track, integrity (`adsl`); label "ADS-L iConspicuity (SRD-860)[, UAS]" | **research** (synthetic; PLAN G10; O-band not covered) |
| **OGN tracker** (OGNTP v1) | `rc_link` (`ogn.py`) | FLARM PHY, sync 0A F3 65 6C | zero Manchester violations over all 208 coded bits + header parity + plausible position (the LDPC matrix is only in unlicensed / GPL code: not used, label "-like") → address, address type, aircraft type (13 = UAV), TEA de-whitened position, altitude, speed (`ogn`) | **research** (synthetic; PLAN G10) |
| **FANET** (paraglider / UAV tracking on LoRa) | `rc_link` (`lora_phy.py`) | CRC-valid LoRa packet with sync 0xF1 (or a strictly valid tracking packet with 0x12); EU 868.2, **India 866.2 MHz** | manufacturer + device ID, message type, tracking (aircraft type 7 = UAV, position, altitude, speed, climb, heading), name, ground tracking (`fanet`) | **research** (byte-level tests; PLAN G10) |
| **FrSky D8 / D16 (ACCST)** | `rc_link` (`cc2500.py`) | 31 / 70 / 77 / 100 kb/s GFSK in 2.37–2.52 GHz | CC2500 sync 0xD391 / 0xD391D391 after the preamble → length byte and CRC: D8 hardware CRC-16, D16 software CRC (+ v2 hardware CRC) → variant (D8, D16 v1 / v2 / LBT v1) and TX ID (`frsky`); label "FrSky D16 v2 (ACCST)" ("-like" without a CRC pass) | **research** (synthetic: 18–20/20 packets at 8 dB, 14/20 at 6 dB; 0.2–0.4 ms) |
| **FlySky AFHDS 2A** | `rc_link` (`afhds2a.py`) | 500 kb/s 2-FSK (±186 kHz) in 2.37–2.52 GHz | A7105 ID code 0x5475C52A after the preamble (≤ 3 bit errors) → "FlySky AFHDS 2A-like"; (7,4) Hamming + CRC-16 decoded under candidate conventions (the datasheet leaves them unverified), reported only on a CRC pass with a known packet type → packet type, TX ID, RX ID, binding, and the convention that passed (`afhds2a`) | **research** (synthetic: decode 16–20/20 at 8 dB, ID-only lower; 0.3–0.8 ms) |
| **XN297 toy drones / Bayang** | `rc_link` (`xn297.py`) | 1 Mb/s GFSK in 2.37–2.52 GHz; also any packet-sized (100–600 µs, 0.3–1.6 MHz) constant-envelope 2.4 GHz burst, because the rate estimate fails at h ≈ 0.3 | prefix 0x710F55 → de-scrambled address (scrambler = 7-bit LFSR seeded 0x71, matches the published table) and bit-reversed payload; Bayang verified by its header + 8-bit sum → TX ID and hop channels (bind), sticks and flags (data); CRC residual reported (length constant unverified) (`xn297`) | **research** (synthetic, full-bit soft demod with the prefix-set level: 30/40 packets at ~12 dB in-band, 10/40 at ~10 dB; 0.3–0.7 ms) |
| **More toy drones (nRF24 / XN297)** | `rc_link` (`toys.py`) | the XN297 path above (any address, 10–21-byte payloads) and plain nRF24 at 1 Mb/s (fixed address by correlation, full-bit phase demod) | XN297: MJXQ (WLH08 / E010, X600 / X800), CX10 family, JJRC345, H8_3D / CG023 / H20, E016H, Potensic A20 by address + length + payload checksum; nRF24: V2x2 (WLtoys V202 / V222 / V262), SymaX, SymaX X5C by address + nRF CRC-16 (with or without the ESB control field) → protocol, address, TX ID; label "Toy drone RC link: …" | **research** (synthetic; 250 kb/s XN297 and HS6200 / RF2500 chips not covered) |
| **RadioLink** (T8FB / T8S air, RC4GS / RC6GS / RC8X surface, DumboRC, RC4G) | `rc_link` (`radiolink.py`) | 38 kBd GFSK (±38 kHz) on the 399.9 kHz CC2500 grid from 2402 MHz | doubled sync 0xABCD / 0xBADC / 0xABA5, PN9 dewhitening, CC2500 CRC-16 + checksum → sub-protocol, TX ID, 16 channels, hop slot checked against the TX ID's 16-channel hop table (LCG shuffles), receiver telemetry; the 9.3 ms packet's pieces join one emission | **research** (synthetic, decodes at 6 dB; AT9S / AT10 use an unpublished protocol) |
| **FlySky AFHDS 1** (FS-CT6B, FS-T6, 9X modules) | `rc_link` (`afhds2a.py`) | as AFHDS 2A (same A7105 ID code), 0.5 ms packets every 1.51 ms | 21-byte packet: header 0x55 / 0xAA, TX ID, 8 channels (µs); verified by a CRC convention or by header + channel ranges + the TX ID's hop table containing the measured channel → TX ID, hop table, bind | **research** (synthetic; CRC init unverified) |
| **Hubsan** (H107 family, H107 Plus, H301, H501 class) | `rc_link` (`hubsan.py`) | 100 kb/s 2-FSK (±186 kHz), 1.5–2.6 ms burst on one of 12 channels 2410–2465 MHz | session ID code after the preamble, 16 bytes summing to 0, model constants → ID code, sticks, bind / telemetry, model | **research** (synthetic, decodes at 6 dB) |
| **Graupner HoTT** | `rc_link` (`hott.py`) | 250 kBd MSK (±62.5 kHz), 0.6–4.4 ms burst on the 1.014 MHz grid from 2404.026 MHz | sync D391 ×2, 4×4 de-interleave + K=4 rate-½ Viterbi (CC2500 FEC; TI DN504 test vector reproduced), PN9, CRC-16 → TX ID, RX ID, hop table (one of 16) + hop index, 12 channels (µs), failsafe / bind, receiver telemetry (`cc2500_rc`) | **research** (synthetic; PLAN G12) |
| **Hitec Optima / Minima** | `rc_link` (`hitec.py`) | 250 kBd MSK on the 0.6 MHz grid from 2405.6 MHz | sync word = the TX ID (id3 id2 sent twice; bind D391), PN9, variable length, CRC-16, sync/ID cross-check → TX ID, 9 channels, repeat flag, bind hop table, telemetry (RX battery) | **research** (synthetic; PLAN G12) |
| **Futaba S-FHSS** | `rc_link` (`hitec.py`) | 128 kBd 2-FSK (±38 kHz), 0.4–1.9 ms bursts on the 1.5 MHz grid from 2404 MHz | sync D391, 13 bytes, CRC-16 in CC2400 mode (init undocumented: 0xFFFF / 0x0000 tried) → TX ID, model hint (T8J / T10J / TM-FH from bytes 3–4), hop index checked against the channel, hop code, 4 × 12-bit channels | **research** (synthetic; PLAN G12) |
| **ExpressLRS GFSK (LR1121 K / DK modes)** | `rc_link` (`elrs.py`) | 300 kb/s FSK (±100 kHz), ≤ 0.7 ms, at 2.4 GHz or 860–930 MHz | a SYNC packet whose CRC verifies with its own UID bytes and whose UID5 / UID4 equal the 16-bit sync word → UID4 / UID5, OTA version, whitening convention; model: LR1121-based module | **research** (synthetic; whitening bit order unverified, both tried) |
| **IEEE 802.15.4 (XBee / Zigbee telemetry)** | `rc_link` (`ieee802154.py`) | a 2.4 GHz burst on an 802.15.4 channel (2405 + 5k MHz), 1.2–4 MHz, ≥ 250 µs (tried before the ESB slicer) | non-coherent O-QPSK: discriminator signs per chip against 16 symbol templates made by modulation; preamble + SFD 0xA7, PHR, FCS (CRC-16 KERMIT) → MAC header: frame type, PAN IDs, short / 64-bit addresses; Digi OUI 00 13 A2 → "Digi XBee" (`ieee802154`, model block) | **research** (synthetic: 10/10 frames at ~7 dB in-band; 0.7 ms per frame) |
| **Any other FSK packet: blind frame structure** | `rc_link` (`fsk_frame.py`) | an undecoded FSK burst, 20 kBd–1.2 MBd, ≥ 10 dB, room for preamble + 32 bits | preamble length, 32-bit sync-word candidate after the alternating preamble, payload length (`features.frame`); per signal `fsk_frame`: sync candidates counted on 24 bits (a repeated one identifies the link), known sync words named as "-like" candidates (CC2500 FrSky, A7105 AFHDS 2A, Crossfire, SiK); nothing decoded or verified | **research** (synthetic, 18 dB) |
| GFSK/FSK RC / telemetry links | `rc_link` | matched-area energy, flank and constant-envelope tests | deviation, symbol rate, hop channels | partial |
| **RC / telemetry protocol candidates** (research catalog) | `rc_link` (`protocols.py`) | every FSK-family signal with a measured symbol rate | ranked `protocol_candidates` with matched / conflicting rate, deviation, band, packet period: FrSky D8 / D16 v1 / v2 / LBT, RadioLink, TBS Crossfire 150 Hz, mLRS FSK, FLARM, Futaba S-FHSS, Graupner HoTT, Hitec AFHSS, ExpressLRS GFSK / FLRC, FlySky AFHDS 2A, nRF24 / XN297 toys, OpenLRS-ng, XBee 900HP, SiK telemetry (13 rates, deviation ≈ rate) | **research** (synthetic / table tests) |
| **SiK telemetry radio + MAVLink** (3DR / Holybro / mRo / RFD900-compatible, 433 / 470 / 868 / 915 MHz) | `rc_link` (`sik.py`) + `core/mavlink` | sub-GHz GFSK burst at a SiK air rate with deviation ≈ rate; sync 0x2DD4 (either polarity) | **NETID**, air rate, whitening convention, radio CRC variant; **MAVLink** frames (v1 / v2, X.25 + CRC_EXTRA verified): vehicle type, autopilot (ArduPilot / PX4), armed, flight mode, **lat / lon / altitude / relative altitude / heading**, GPS fix / satellites, home position, ground / air speed, climb, throttle, battery V / %, attitude, status text, radio RSSI; label "SiK telemetry radio (MAVLink)" | success (MAVLink CRC) — **research**: synthetic only, full decode ≥ 9 dB; Golay (ECC) / AES links not decoded |
| **MAVLink over LoRa and open Wi-Fi** (Ebyte E22 / E32 / E220 transparent links, ESP8266 / DroneBridge bridges) | `rc_link` (`lora_phy.mavlink_in`), `wifi` (`mac80211`) | any explicit-header LoRa payload; UDP / TCP inside an unencrypted 802.11 frame (ports 14550 / 14555 / 14540 / 5760 named) | CRC_EXTRA-verified MAVLink frames → the SiK MAVLink summary (vehicle, autopilot, armed, mode, position, …); a vehicle sysid flags the frame as a drone | **research** (synthetic; a frame split over two LoRa packets is not reassembled) |
| **wfb-ng / OpenIPC / OpenHD digital FPV** (raw 802.11 injection) | `wifi` (`mac80211`) | every decoded 802.11 frame | broadcast data frames with TA = BSSID = 57:42 + channel_id → link_id (default link domain 0x7505D6 named), radio port and stream (video / MAVLink / tunnel / bind), packet type, FEC block / fragment; content encrypted (not decrypted) | **research** (synthetic frames; the 10 MHz RTL8812EU mode not tested) |
| **Drone Wi-Fi links by behaviour** (encrypted or not: toy / Wi-Fi drones, DJI Spark / Tello class, DIY links) | `wifi` (`field.py`, both PHYs) | every BSS in a decoded Wi-Fi signal, including PPDUs whose payload is out of reach (attributed by `links.py`) | per BSS (`wifi_links`): AP, clients, downlink / uplink PPDUs, airtime, bytes, big frames, downlink gaps, video payloads; `drone_links` verdict: **identity** (Remote ID, Tello SDK, MAVLink, wfb-ng, ESP-NOW, drone SSID / OUI) or **behaviour** (one client, downlink ≥ 3× uplink airtime, ≥ 5 % airtime or ≥ 1 Mb/s, small uplink frames → "drone / FPV video link-like", 0.55; 0.65 with a video payload). Label suffix ", drone link-like (…)" (PLAN W1; Alipour-Fanid TIFS 2020, Alsoliman ICDCN 2021) | **research** (synthetic; one block cannot separate an IP camera or a phone hotspot streaming video: the verdict says so) |
| **Open drone links' payloads** | `wifi` (`mac80211`) | UDP inside unencrypted 802.11 frames | **Tello SDK**: commands (UDP 8889, Wi-Fi credentials withheld), state (UDP 8890: battery, height, ToF, attitude, speeds, temperatures, flight time) → "DJI / Ryze Tello (Tello SDK)"; **video**: raw H.264 / H.265 Annex B NAL units, RTP with a dynamic payload type carrying H.264 (`ipv4.video`); vendor-documented hotspot names (RoboMaster TT, Holy Stone FPV, Eachine-class WiFi-720P, ExpressLRS TX / RX, ArduPilot / PixRacer MAVESP8266, CUAV PW-Link, DroneBridge) as platform signatures, and model-level names for Wi-Fi drones (DJI Spark / Spark RC, Mavic Air, Mavic Mini, Phantom 3 Standard / 4K, Tello; Parrot Anafi; Yuneec CGO3 camera link; Ruko F11; Holy Stone FPV; Skydio / Hubsan by OUI; research file 18) (PLAN W3 / W4) | **research** (synthetic frames) |
| **ESP-NOW and the ExpressLRS Backpack** | `wifi` (`mac80211`) | 802.11 vendor action frames (category 127, Espressif OUI) | ESP-NOW version and body; a CRC-valid MSP v2 frame between UID-derived MACs → **ELRS bind UID**, command, VTX band / channel (default Betaflight VTX table) | **research** (synthetic frames; ESP32 Long Range rate is not 802.11 and not decoded) |
| DJI OcuSync video / RC uplink | `dji` | 15/30 kHz CP-OFDM fold + ZC root (10/20 MHz mode from the ZC length); 15 kHz OFDM hop tracker | ZC root, bandwidth, hop pattern | partial |
| DJI DroneID | `dji` | 9-symbol (OcuSync 2/3) or 8-symbol (legacy, Mavic 2 / older; *research*) ZC burst | serial, GPS (external `droneid_api`, turbo + CRC; 9-symbol frames); `droneid_frame` = generation (*research*) | success |
| Analog FPV video (PAL/NTSC FM) | `analog_vtx` | line-rate fold of the sync tips (to about −12 dB SNR); needs ≥ 20 MHz of band; a "line clock" locked to the sample clock (fs/N within 2 ppm) with its sync tip on the LO and no field sync is rejected (capture-chain artifact) | standard, line rate/ppm, lines per field, sync tip, deviation, audio subcarrier (*research*: demodulated → level, dominant tone such as the propeller blade-pass, flatness), channel, **picture**; **colour standard NTSC / PAL / PAL-M / PAL-N / SECAM or monochrome** (`color_standard`, label PAL-M / PAL-N; line-locked burst phasors at each standard's subcarrier, PAL ±45° swing; SECAM from the line-alternating Db / Dr back-porch subcarrier (4.250 / 4.406 MHz, rest frequencies measured); subcarrier frequency = camera crystal; *research*: correct to 7 dB in-band SNR on synthetics, 2.3 ms; real NTSC (Analog_20260804) and real PAL (LIVE_ANALOG_VTX_PAL_5740MHz, colour camera) correct; the burst gate is detrended per line (porch level + sync-edge tail leaked over every frequency and hid the real PAL burst); Lowband LB channel table *research*; 3.3 GHz analog band 3150–3490 MHz searched and labelled with the two 8-channel groups 3G3A 3330–3470 / 3G3B 3170–3310 MHz, *research*, PLAN G9; the wide 4.9–6.0 GHz groups of 80/96-channel VTXs are detected but have no channel names: the vendor tables are published only as images); streaming video, one field per block (`FieldDemodulator`: power-weighted FM discriminator, line clock from every line's sync edge (H-PLL fit), sync-tip noise estimate, 3x3 median below ~5 dB SNR), **in colour** (NTSC / PAL: burst-locked chroma, burst gate tried at 3 offsets, PAL V switch + delay line, burst ACC, compressed-sync levels, burst and noise colour killers) when the VTX sends a colour burst | success (picture) |
| Analog FM voice: walkie-talkies (PMR446, LPD433, FRS), PMR | `voice_fm` | narrow peaks in 136–174 / 300–520 / 806–870 MHz; keyed, constant-envelope FM (envelope varies no more than the SNR explains, unless a sub-audible code is decoded), deviation ≤ 6 kHz; no 2400/4800 sym/s FSK line (digital PMR); the receiver's own clock lines skipped | channel (PMR446/LPD433/FRS plan), **DCS code** (Golay(23,12)-verified, with its on-air aliases) or **CTCSS tone**, deviation, keyed time, speech activity, **audio** (`<name>_voice<MHz>_audio.wav`; streaming `VoiceDemodulator`); clipping products of a strong handheld suppressed | success (DCS/CTCSS verified) / partial |
| TETRA downlink (380–400, 410–430, 450–470, 870–876, 915–921 MHz) | `pmr_digital` | 18 ksym/s line of the matched-filtered envelope; π/4-DQPSK; synchronisation / normal training sequences on the 510-bit slot grid | **BSCH of each synchronisation burst: MCC, MNC, colour code**, system code, TN/FN/MN, sharing mode, cell service level | success (CRC-CCITT) / partial (training sequences) / "TETRA-like carrier" (symbol-rate line only) |
| DMR (Tier II/III) | `pmr_digital` | 48-bit SYNC (BS/MS/direct, voice/data), correlation ≥ 0.9 and DMR deviation from the sync's own level fit | **colour code and data type** (slot type, Golay(20,8)), base station / mobile / direct, voice activity, BPTC(196,96) Hamming check, CSBK opcode / FID when its CRC passes | success (CSBK CRC) / partial |
| NXDN (4800 / 2400 sym/s) | `pmr_digital` | 20-bit frame sync word on the 192-symbol frame | symbol rate, channel width | partial |
| Herelink drone link | `digital_vtx` | LTE PSS fold, PSS purity ≥ 0.5 (rejects DJI sync), 5 ms TDD rhythm, ≥ 2400 MHz | carrier, PCI (SSS), CFO/ppm, control / video / uplink subframes | partial |
| **HDZero digital FPV video** (27 MHz mode) | `digital_vtx` (`hdzero.py`) | on every HDZero channel inside the capture (R1–R8, E1, F1, F2, F4, L1–L8): band cut at 1024 × 27 MSps from the block FFT, CP-fold timing + fractional CFO, then the 26 continual pilots must be coherent symbol to symbol (≥ 0.2 and ≥ 4× the data carriers beside them) over the integer carrier offsets; ~0.1 s per channel on a 2 M-sample block, 3 ms on 65 k | channel, carrier and ppm, 23.44 MHz / 889 subcarriers, per-frame preamble → video frame period / fps, SNR; model block (every HDZero VTX: identical on air, unresolved). Measured on a real HDZero Freestyle V2 (R1, 60 fps); detected to about −6 dB link SNR (32 ms), frame rate to −3 dB | partial (video is JSCC: no published format) |
| **Any CP-OFDM drone video / data link, of any make (also unknown ones)** | `digital_vtx` (`blind.py`) | no vendor knowledge: coherent lag-Tu autocorrelation of every 10/20 MHz slice (time-AGC'd, so Wi-Fi bursts do not dominate), confirmed by a symbol-rate line or per-window CP fold; Tu ≥ 10 µs (SCS ≤ 100 kHz); 1.35–1.444 / 1.6–2.53 / 4.4–6.0 GHz (PLAN G8: L-band, the Chinese 1430–1444 MHz UAV band and the 1.6–2.33 / 4.4–4.9 GHz government / professional bands, label suffix "in the … band" and `band_context`); slices and candidates whose band and IQ mirror are both cellular are not scanned (a capture inside B3 / B40 costs no blind scan) | Tu, CP, SCS, symbol period, occupied band (CP cross-spectrum, then in-band power while on vs off: multipath dips and Wi-Fi do not cut or widen it), bursts on local symbol grids (survives a timing re-anchor) extended / joined by in-band power, TDD frame period, duty cycle, SNR; the numerology table names the family (DJI OcuSync / O4-like, Lightbridge / WiMAX-like, Herelink-like, PMP 450-like fixed wireless, unknown). A link a protocol receiver identifies (DJI, Herelink) is attributed to it. **Transmit antennas** (`streams.py`, links ≥ 20 dB SNR, ≥ 200 subcarriers): per-subcarrier kurtosis and the Alamouti pair power correlation on runs of ≥ 20 full-band symbols → `2tx_stbc_like` / `2tx_sfbc_like` (label suffix `2-TX STBC-like`) / `1_stream_like` / unknown, in `decoded.fields.tx_antennas`. **System / model candidates** (`catalog.py`, *research*): ranked `system_candidates` (DJI O1 / OcuSync 1 / 2 / O3 / O4, O4 Air Unit Lite vs Pro, Walksnail dual- vs single-antenna air unit, AR8030 family, Lightbridge 2, Herelink, Autel, HDZero, PMP 450), each with what matched and what conflicted (Tu, CP, BW mode, band, channel raster, antennas, frame) | partial (`-like` / unknown) |
| **Hopping narrowband CP-OFDM RC uplinks, of any make** (DJI O4 30 kHz 4.4 MHz hops, OcuSync 2 / Autel 15 kHz 1.1 MHz hops, unknown ones) | `digital_vtx` (`blind.py` `blind_hops`) | windowed CP-fold scan at every numerology found in the capture and at the drone families' (15 / 30 kHz); each hit's edges from its envelope, band from the on-vs-off spectrum, kept only with a significant CP fold; hops of one link share dwell and width; missed hops re-tested at the hop cadence | per hop: time, band, channel, CP-fold z; per link: numerology, hop width / dwell, channels, hop period, the video link of the same numerology | partial (`-like` / unknown) |
| **COFDM video links with DVB-T framing** (UAV downlinks: Domo/DTC, IMT/Vislink, HawkEye, cheap "COFDM 2K" modules; DVB-T broadcast), 8 / 7 / 6 / 5 MHz and scaled 4 / 2.5 / 2 MHz, 200 MHz–6 GHz | `cofdm` | Welch block 0.951·B wide → CP fold at the scaled DVB-T clock (2K / 8K × GI) → continual pilots → TPS | mode, guard interval, bandwidth, exact carrier, **TPS: constellation, hierarchy, HP / LP code rate, frame number, cell-ID byte** (BCH-verified); label "COFDM video link (DVB-T 2K, 8 MHz, QPSK, HP 1/2, GI 1/8)", technology `cofdm_video_link`; **DVB-T2 P1 preamble** when one falls in the block: S1 (T2 SISO / MISO / Lite / FEF), FFT size and guard-interval class, mixed bit (*research*, synthetic to −3 dB) | success (TPS BCH) / partial — **research**: synthetic only, TPS decoded to 0 dB in-band SNR (2K, 8 / 7 MHz); 2 / 4 MHz and 8K need > 1 block for TPS |
| Drone emitters | `emitters` | association of the drone receivers' bursts | radios (hop pattern, frame clock, oscillator ppm), confidence, evidence | report `emitters` |
| LTE downlink (FDD/TDD) | `lte` | PSS 5 ms fold on the EARFCN raster of every LTE-band energy detection ≥ 3 dB, continuous or bursty (a lightly loaded or TDD downlink is bursty in energy); a bursty / weak candidate must give a decisive SSS | PCI (SSS, maximum-ratio combined with the PSS channel estimate smoothed over 225 kHz), duplex, CP, EARFCN, band (cells only on an LTE band raster), bandwidth (largest 3GPP bandwidth whose band is always-on and flat, 10–90th percentile within 8 dB), receiver ppm | partial (no MIB yet) |
| 5G NR (SSB) | `nr` | PSS on the GSCN raster | PCI (differential SSS), SSB index (DMRS), **MIB** (SFN, SCS, k_SSB, CORESET#0, barred) | success (CRC-24C) |
| GSM BCCH | `gsm` | FCCH tone on the ARFCN raster, searched within ±5 ppm of the receiver ppm (±30 ppm when unknown) on 0.4 ms STFT frames; every FCCH of a block must agree on the receiver error (±1.5 ppm of the median of the clean ones) | carrier, ARFCN, CFO, receiver ppm from the FCCH; SCH → BSIC, frame number when CRC-10 passes | partial / success |
| **Swept-carrier jammers** (one carrier swept periodically across a band, continuously: the common 2.4 / 5.8 GHz drone / Wi-Fi jammer) | `jammer` | anywhere, fs ≥ 20 MHz: cyclic line in the per-bin STFT occupancy (64 bins, ~1 µs rows) at the sweep rate; the STFT folded on the true period (shortest multiple of the line period with ≥ 80 % of the best phase structure; refused when a longer period has more: pulsed signals) shows one dominant carrier moving continuously (≤ 1 flyback, ≤ 2 band exits) over ≥ 10 MHz for ≥ 30 % of the period; video line rates refused | sweep rate / period (drift-corrected: ±10 ppm on fast sweeps), frequencies reached, side it leaves the capture through, shape (up-down / sawtooth), time in band, dominance, mean trajectory; owns its band (its energy pieces are absorbed). "Jammer (swept carrier)", technology `swept_jammer`, family `jamming`; in the radio-altimeter band "Swept carrier (FMCW-like)" (live: `LIVE_JAMMER_SWEEP_2460MHz`) | success (sweep measured; no payload exists) |
| ADS-B / Mode S 1090 | `adsb` | preamble correlation, PPM | every field of DF0/4/5/11/16/17/18/20/21 (identity, category, position, velocity, intent, integrity, capabilities, Comm-B air data, ACAS/radar reply fields), plus a per-aircraft summary with CPR position, wind/temperature and who interrogates it (pyModeS) | success (CRC-24) |
| **UAT 978 MHz ADS-B** (DO-282B downlink) | `adsb` (`uat.py`, plugin `uat`) | capture covering 978 MHz: CPFSK h 0.6 at 4 samples / bit, sync 0xEACDDA4E2 (≤ 4 bit errors) with the slicer centred on the sync's levels | own Reed-Solomon (GF(256) 0x187, fcr 120; parity checked against a reference implementation): Long RS(48,34) / Basic RS(30,18) as dump978 → address + qualifier, position, altitude, velocities, NIC, emitter category (**14 = UAV**), callsign; per-aircraft summary | **research** (dump978 sample messages + synthetic IQ; PLAN G14; no UAT in India) |

| **SSR Mode A/C replies** (and IFF Mode 3/A) | `adsb` (`modeac.py`) | 1090 MHz envelope: framing pulses 20.3 µs apart within 3 dB, every pulse between them on the 1.45 µs grid, SPI | Mode A code (octal) and, when valid, the Gillham Mode C altitude for every reply (the interrogation is not in view); 7500 / 7600 / 7700; repeated words are the verified ones (*research*, signal gaps S11) | partial |
| **AM voice** (VHF civil air 118–137 MHz, UHF military air 225–400 MHz) | `voice_fm` | per carrier: the envelope's 300–3400 Hz rms ≥ 4 % of the carrier in ≥ 10 % of 20 ms frames while the phase stays quiet, aeronautical bands only | service, 25 kHz channel, 121.5 / 243.0 MHz emergency, keyed time, **AM audio** (*research*, S13) | partial |
| **APRS / AX.25 packet radio** | `voice_fm` (`afsk.py`, `core/hdlc`) | Bell 202 AFSK (1200 / 2200 Hz) on any FM voice carrier | HDLC frames with CRC-16/X.25 → AX.25 UI addresses / path; APRS uncompressed, compressed and Mic-E positions, status, messages, objects (*research*, S19) | success (CRC) |
| **ILS glide slope / localizer, VOR** | `navaids` | carrier on the ICAO channel grid (glide path 329.15–335 MHz, localizer / VOR 108–118 MHz) | ILS: 90 / 150 Hz depths, DDM, SDM, side, paired channel, two-frequency systems; VOR: radial, AM depth, 9960 Hz FM deviation (*research*, S9 / S21) | partial |
| **FM broadcast + RDS** (87.6–107.9 MHz) | `broadcast_fm` | 100 kHz raster, ≥ 10 dB, deviation 15–110 kHz, 19 kHz pilot | stereo; RDS (≥ 45 ms): PI, PTY, TP / TA, PS name, RadioText, clock time from error-free (26,16) blocks (*research*, S20) | success (RDS) / partial |
| **AIS** (161.975 / 162.025 MHz) | `marine` | both channels, GMSK 9600 via a ±10 kHz discriminator | NRZI, HDLC, CRC-16/X.25 → types 1/2/3, 4, 5, 18, 21, 24: MMSI, flag (419 India), position, SOG / COG / heading, name, call sign, ship type, dimensions, destination; per-vessel list (*research*, S15) | success (CRC) |
| **GNSS: GPS L1 C/A, jamming and spoofing** | `gnss` | block covering 1575.42 MHz (≥ 2.1 MSps, ≥ 16 ms) | 32 PRNs: code phase, Doppler, C/N0; J/N and narrow lines near L1; spoofing indicators S1–S4; **swept GNSS jammers (personal privacy devices**, fs ≥ 30 MSps: period, span, centre, crosses L1) (*research*, S3 / S12) | partial |
| **Bluetooth LE trackers** (AirTag / Find My, SmartTag, Google FMDN, Tile, Chipolo, DULT) and notable devices (iBeacon, Eddystone, AirPods, Swift / Fast Pair) | `rc_link` (`core/ble_trackers`) | a CRC-valid BLE advert | vendor, kind, tracker flag, state / battery where the format carries them; label "Bluetooth LE tracker (…)"; no stable identity (rotating tokens) (*research*, S7) | success (CRC) |
| **Bluetooth BR/EDR (Classic)** | `rc_link` (`bredr.py`) | a 1 Mb/s GFSK hop that is not a BLE advert | 64-bit access code within 3 bits of the (64,30) BCH codeword of its own LAP → **LAP**, inquiry (GIAC / LIAC); UAP / NAP (vendor) are not on air in one block (*research*, S7) | partial |
| **315 / 433 MHz remotes and TPMS** | `ism_devices` | bursty carriers in 314–316 / 433.05–434.79 MHz; OOK pulses (envelope) or a two-level FSK discriminator | **TPMS** Continental / VDO and Toyota PMV-C210 (CRC-8): sensor ID, pressure, temperature; **EV1527 / PT2262** fixed codes (ID, data, tri-state; repeats as evidence); **KeeLoq** clear serial / buttons / battery (the hopping code is not reported); other OOK / FSK devices "-like" (*research*, S6) | success (TPMS CRC) / partial |
| **Radio altimeters** (4200–4400 MHz FMCW) | `radalt` | 128-point STFT ridge, straight crossings (R² ≥ 0.99) through the capture | sweep rate, triangle / sawtooth, turnarounds, sweep bandwidth, period (ITU-R M.2059 types) (*research*, S8) | partial |
| **Noise jammers** (barrage noise, noise-FM) | `jammer` (`barrage.py`, plugin `jammer_barrage`) | an unexplained energy box ≥ 1 MHz, continuous, ≥ 6 dB | flat Gaussian plateau (k = E\|x\|⁴/E²\|x\|² ≈ 2) or constant envelope with Gaussian-like IF; no cyclic prefix, no cyclic line, stationary; allocation match (`covers_band`) (*research*, S12) | partial |
**Wi-Fi MAC dissection** (`receivers/wifi/mac80211.py`, both PHYs):
- addresses and their roles, with the vendor from the IEEE registry; randomised addresses are flagged;
- for an AP: SSID (or hidden), security (WPA2/WPA3/OWE, ciphers, MFP), PHY generation (802.11b/g/n/ac/ax/be),
  uptime from the TSF, country, BSS load, HT/VHT operation, and WPS / Wi-Fi Direct device identity;
- for other frames: probe targets, auth/assoc status, deauth reasons, CCMP packet numbers, and the EAPOL
  handshake message;
- Remote ID and drone platform when a frame carries them (`receivers/wifi/droneid.py`);
- a one-line `summary` per frame.

**Wi-Fi transmitter attribution** (`receivers/wifi/links.py`): every 802.11 OFDM signal lists its `transmitters`
(MAC, vendor, PPDU count, formats, frame types, median received power and CFO, scrambler-counter behaviour), and
every PPDU its `ta` with `ta_evidence`, strongest first:
- **mac header:** an FCS-valid MPDU's TA;
- **response:** an ACK / Block Ack / CTS names in its RA the station that just transmitted (and an ACK/CTS is sent by
  whoever the previous frame addressed);
- **scrambler counter:** many chipsets (Atheros/Qualcomm, Microhard, DJI, the 5240 MHz station in the sweep) step the
  7-bit scrambler seed by one per PPDU. The seed is read from the SERVICE field (self-checked by 9 zero bits) even
  when the payload fails, and a gap in a station's counter that equals the number of unattributed PPDUs in it
  attributes them all — including 2-stream frames no single antenna can decode;
- `ta_likely` (weak): the one transmitter within 3 dB in received power, for what is still unattributed.

Everything else is **classified, not decoded**. The rules (`classify/rules.py`) combine the India NFAP-2022 band
plan with measured features: timing, envelope, cyclic lines, OFDM CP, FSK. This covers:
- UMTS-like, DVB-T2-like, PMR / TETRA-like (band and bandwidth only; `pmr_digital` reports the verified ones);
- radar and DME pulses (at least 2 pulses in the block: DME replies are pulse pairs, radars repeat). Pulses are
  measured on the signal's own band with a 1 µs envelope over the whole block (`features/measure.pulse_train`:
  pulse groups, width, PRI as the shortest interval cluster so staggered PRIs count, and a linear
  instantaneous-frequency fit for LFM chirps); the spectrogram's burst boxes are far longer than radar pulses. In
  2700–2900 MHz this gives "Pulsed radar (S-band ATC)-like" (live: 2888.57 MHz, staggered PRI ~1.17 ms, 1 µs + 69 µs
  ~0.9 MHz LFM pulse), and an LFM chirp is no longer read as OFDM. Measured only in radar / DME allocations
  (`core.bands_india.PULSED_BANDS`). Gates: bandwidth × pulse width ≥ 1 (shorter "pulses" are noise peaks), a DME
  label needs a pulse-pair group of 10–45 µs, nothing within 1.5 MHz of 1030/1090 MHz (SSR) is DME, and a PRI is
  claimed only when its interval cluster holds at least half of the intervals. Narrower detections of one wideband
  pulse train (same PRI within 1 %, the same pulse starts within 5 µs for at least half of the pulses, spans at most
  2 MHz apart) are merged into one signal together with the burst-only slivers on its pulses
  (`detect.merge.merge_pulse_trains`), and an autocorrelation lag longer than a quarter of the pulse width is not read
  as OFDM. In U-NII-2e (5470–5725 MHz, shared with Wi-Fi) a radar label needs a PRI: bursts without one are
  packets (live: Wi-Fi edge fragments at 5495.65 MHz were "Pulsed radar-like");
- **named radar classes** (research, signal gaps S10; `classify.rules._radar_class`): S-band ATC (solid-state when a
  long LFM pulse is present), S- / C-band weather, marine S-band, C-band, L-band surveillance, wind profiler, from the
  band plan, pulse width, PRI and LFM; the class only raises the confidence of a "-like" label;
- **pulse deinterleaving** (research, S10; `core/deinterleave`, feature `pulse_emitters`): a pulse train in a radar
  band is split into emitters by PRI (PW groups, cumulative difference histograms, sequence search with a chance test
  per chain, staggered PRIs as frames with their sub-intervals, a 10 % jitter mode); each emitter's PRI, stagger,
  jitter, pulse count and width go into the features and the pulsed-radar evidence (synthetic: interleaved constant,
  staggered and dual-PRF trains separated, no emitter from Poisson pulse times in 300 blocks);
- ISM short-range devices, CW carriers;
- Iridium-like bursts, only with a PSK signature (x² / x⁴ carrier line or the 25 kHz symbol line): a narrow burst
  without one is as likely a wandering clock-comb line.

An inferred label carries the "-like" suffix; below 0.5 confidence the label is `unknown`, with the evidence.

**Receiver artifacts are flagged, not reported as emissions:**
- LO leakage;
- IQ images;
- ADC clipping;
- receiver spurs, in one block (`core/spurs.py`, flag `receiver_spur`): a narrow line at an exact harmonic of the
  sample clock (k × fs: the Pluto's 430.08 MHz = 7 × 61.44 MHz), or, when the receiver ppm is given, a line on the
  integer-MHz grid of the receiver's own frequency scale (locked to its reference: an external line there would be
  displayed off by the receiver error, ~4 kHz at 550 MHz). The line is refined to ~1 Hz on the full-length FFT;
  receivers skip these lines (`ctx.cache["receiver_lines"]`). Over a sweep, also from the cross-tuning check;
- LO-harmonic responses, confirmed over a sweep (`summary.py`). In one block a switching mixer's 3rd/5th
  LO-harmonic images (no preselector on the Pluto) cannot be proven, so an unidentified signal outside the cellular
  bands that such an image of a cellular downlink would explain gets the flag `possible_lo_harmonic_image`, the
  candidate source bands (`features.lo_harmonic_image_sources`) and a CAUTION in its evidence; label and
  confidence are unchanged (a hint, not a verdict).

Detections below 6 dB without a decode go to `weak_detections`.

### Live LO-by-LO bench sweep, 300–2012 MHz (2026-10-10)

A Pluto (61.44 MSps, 56 MHz analog band, gain 40 dB, LO step 36 MHz) swept one LO at a time: the rx buffer was
destroyed before each retune, 4 buffers discarded (buffer 0 still held the previous LO; the spectrum settles from
buffer 2), 3 snapshots kept, analysed, and deleted before the next LO. Every detection was then checked on the IQ
(`bench/live_lo/`, findings in `bench/live_lo/FINDINGS.md`); the regression captures are in `testdata/live/`.

| Found | How | Library now |
|---|---|---|
| **TETRA** 390.03 + 391.53 MHz: **MCC 404 (India), MNC 4040, colour code 7** | π/4-DQPSK, sync bursts, BSCH CRC (52 + 19 of 78 in 1.75 s; 2 + 2 in one 34 ms block) | `pmr_digital`, success |
| **DMR** repeaters 859.82 (voice call) and 860.82 MHz (control), colour code 0; **NXDN** 859.35 MHz | SYNC / slot type / BPTC; NXDN frame sync | `pmr_digital`, partial (their CSBKs fail the standard CRC: vendor control channel) |
| **LTE B28** 760.5 MHz PCI 351, **LTE B5** 874.0 MHz PCI 169 | PSS/SSS | `lte` (missed before: bursty energy, weak SSS) |
| **5G NR n28** SSB 782.55 MHz, PCI 61 and PCI 98 (two cells) | MIB CRC | `nr` (unchanged) |
| Most other "signals" at 330–830 MHz: **LO-harmonic images** of the 1.8/2.1/2.3 GHz cellular bands (an LO step of +1 MHz moves them −4 MHz = 5th harmonic, +4 MHz = 3rd) | LO-dither test | `possible_lo_harmonic_image` caution; the fix is a preselector filter |
| Receiver / host lines: 307.2, 430.08, 921.6 MHz (k × 61.44), 536.000 / 552.000 MHz (reference-locked), Zynq 33.333 MHz harmonics (DDR 533.33 MHz), USB-device clock combs | exact frequency vs receiver ppm | `receiver_spur` (the first two kinds) |
| Capture self-noise: the noise floor jumps 4–5 dB ~12 ms into every buffer, modulated at the USB microframe rate (8 kHz) | buffer-time locking | acquisition side (USB cable/ferrite, antenna away from the host) |

Removed false labels: 6 "LoRa / ExpressLRS preamble SF6" (the capture self-noise), 2 "LoRa" inside LTE carriers,
"Analog FM voice" on clock / EMI lines (420.003, 430.08, 519.998 MHz) and "Analog FM voice + CTCSS" on DMR/NXDN.

**896–2012 MHz (same method, LO 904 … 1984).** The Pluto TCXO had drifted to +9.5–10 ppm (LTE- and FCCH-measured),
so the LOs were analysed with `receiver_ppm` 10. The library fixes that followed are in the last column; the
regression captures are `testdata/live/LIVE_940MHz`, `LIVE_1840MHz` and `LIVE_ADSB_1090MHz`.

| Found | How | Library today |
|---|---|---|
| **LTE B8** PCI 96 (940.8 MHz) and PCI 101 (949.0 MHz); **LTE B3** PCI 472 (1810.1), PCI 110 (1841.5), PCI 101 (1826.5, 20 MHz) | PSS/SSS | `lte`: all, with bandwidths (B3 PCI 101 since the single-block floor bridges wide carriers) |
| **GSM-900** BCCH carriers ARFCN 11/12/44–54/97–109 (13–16 per block, all at the same +9.6–9.8 ppm), **ARFCN 106 / 102 BSIC 47** (NCC 4, BCC 7); **GSM-1800** ARFCN 566 | FCCH + SCH CRC (bench `gsm_dl.py`, `gsm_fcch_cfo.py`) | `gsm`: FCCH searched over the CFO (ppm-independent, before: none with a 2.4 kHz ppm error) and checked for a common receiver error; SCH still mostly fails its CRC (no equaliser) |
| **GSM handset uplink** 899.6 MHz (ARFCN 48), training sequence 3 in 7/7 bursts with 0 errors (serving cell BCC 3) | GMSK demod + TSC match (bench `gsm_ul.py`) | not detected (`unknown`) |
| **ADS-B** (1 s dwell): 17 frames, aircraft 8015CB, 8015D1, 80175F, 80028D | Mode S CRC | `adsb`, success |
| **5G NR n78 PCI 972** (SSB 3449.28 MHz, MIB CRC 3/3) and **5 GHz Wi-Fi ch 36 BSSID c4:27:28:fa:d1:8f**, both decoded *through the LO's 3rd-harmonic image* at 1128–1184 / 1718–1758 MHz | conjugated IQ analysed at 3 × LO (bench `image_decode.py`) | n78 images flagged `possible_lo_harmonic_image`; Wi-Fi images not flagged |
| Clock combs: 27.008 MHz × 55–66 (1485–1783 MHz, 47 Hz wobble), 40 MHz × 40 / 48 (1600, 1920 MHz), DDR 1066.66 MHz, k × 61.44 MHz | exact harmonic numbers, LO dither | `unknown` / `receiver_spur` |

False labels found and removed: "Analog FPV video (PAL)" on a capture-locked line at the LO (1090 MHz dwell at 8 MSps,
below the receiver's 20 MHz minimum), "Iridium-like TDMA burst" on a 27 MHz comb line (1620.48 MHz: no PSK signature),
"DME/TACAN-like" on a narrow single blip (1201 MHz: one pulse).

**2012–4000 MHz (2026-10-15/16 night, LO 2020 … 3992, `receiver_ppm` 10).** Regression captures:
`testdata/live/LIVE_ISM_2452MHz`, `LIVE_SBAND_RADAR_2884MHz`, `LIVE_PULSED_3040MHz`.

| Found | How | Library today |
|---|---|---|
| **LTE B1** PCI 103/415/426, **B40** PCI 103/350/196, **B41** PCI 218 | PSS/SSS | `lte` |
| **5G NR n78 PCI 102** (SSB 3511.20 MHz) and **PCI 591** (SSB 3449.28 MHz), MIB CRC; PCI 531 (SSB 3331.20, PSS/SSS only) | SSB, PBCH | `nr`; wide n78 carriers are still split into an "NR-like" entry plus `unknown` pieces |
| Wi-Fi 802.11b/g/n SSIDs and BSSIDs, BLE adverts (CRC) on 37/38/39 | DSSS / OFDM, BLE CRC | `wifi`, `rc_link`; GFSK hops off the advertising channels are "Bluetooth (BLE/BR)-like" |
| **2 Mb/s GFSK "hopper" at 2426 / 2450 MHz** (2026-10-17, LO 2440): 94 µs packets every 1.000 ms, two nodes 145 µs apart (sync 0x31A9D906 / 0x31A9D886), 20–40 ms per channel; idle traffic repeating every 14 packets (3-bit sequence × 2 payloads); CRC-16/CCITT ending at bit 124, init 0xDD19 at the sync word; a second device at 2430.1 MHz (sync 0xC7B6A60A, 2 packets) | 2 Mb/s slicing, joint CRC search (`esb.py`) | "2.4 GHz GFSK link (ESB-like, 2 Mb/s)", `success` (`LIVE_ESB2M_2440MHz`: 44 of 47 packets CRC-valid; was 1–4 MHz `unknown` entries) |
| **S-band ATC radar** 2888.57 MHz: staggered PRI ~1.12–1.21 ms, 1 µs + 69 µs LFM (~0.9 MHz) | 1 µs envelope, PRI cluster, LFM fit | "Pulsed radar (S-band ATC)-like" |
| **Pulsed emitter 3036–3047 MHz**: 10–15 µs pulses, PRI 4.000 ms (250 Hz), in every dwell | dwells at 10 and 61.44 MSps; quiet control at 3100 MHz | one "Pulsed radar-like" signal (was 3 radar entries + slivers: `detect.merge.merge_pulse_trains`) |
| Pluto's own lines: 2666.67 MHz comb (5 × DDR 533.33, 1.92 MHz sidebands), 3200.0 MHz (6 × 533.33), a beat pair at LO ∓ 3 × (LO − 3266.67 MHz) near LO 3266.67 MHz, 3932.16 MHz (64 × 61.44) | fixed RF under LO steps / pair moving with the LO | `receiver_spur` for k × 61.44 only; the others are still `unknown` |

**One emitter, one signal (2026-10-16):** a wideband burst that the energy detector splits at its fading notches
(an 802.11 PPDU the receiver could not decode) is absorbed into the receiver signal whose span holds it when >= 3
burst-only pieces start and end together (+-20 us) and span >= 8 MHz (`analyze._absorb_coincident`; live: 5 pieces
at 5739–5754 MHz joined the 5.8 GHz Wi-Fi on ch 149, `testdata/live/LIVE_WIFI_5755MHz`). A narrow hop or RC
packet inside a Wi-Fi channel is a single piece and stays its own signal.

**Re-check of the known signals on fresh live captures (2026-10-16 night):** TETRA (a second network, **MCC 404 MNC
2024 colour code 1** on 395.159 / 395.409 / 395.905 MHz: 81 / 61 / 33 CRC-valid BSCH in 1.75 s,
`testdata/live/LIVE_TETRA_395.28MHz`; the MNC 4040 carriers were off), DMR colour code 0 and NXDN (856–862 MHz), NR
n28 PCI 75 (MIB), LTE B5 / B8 / B3 / B1 / B40 / B41 cells, GSM-900 FCCH on 13 ARFCNs, ADS-B (2 aircraft), Wi-Fi and BLE
were all found and decoded again. Removed false label: "SiK-like" on GSM-900 traffic bursts (957.0 / 958.0 MHz, 554 µs
GMSK on the 200 kHz raster: `rc_link` now leaves 1–8 GSM time slots on the GSM raster in a GSM band to the GSM side;
`testdata/live/LIVE_GSM_958MHz_fsk`).

**Swept jammer (signal-vault capture, 2026-10-19, confirmed by the user):** at 2460 MHz / gain 40 the Pluto's ADC
was saturated (31 % of the samples clipped) by one carrier swept up and down every 100.00 µs (10.000 kHz) from
~2474.5 MHz to below the capture's 2429.3 MHz edge; ADC clipping spread it over the whole band. The pipeline reported
nine "unknown" single-carrier boxes (the sweep rate read as a 10 / 20 kHz "symbol rate"). The `jammer` receiver now
reports one "Jammer (swept carrier)": 10000.25 Hz (long-baseline reference 10000.09–10000.21 Hz), seen
2435.3–2474.6 MHz and continuing below the capture, up-down, in band 88 % of the period; it absorbs the energy
pieces (including the one cut at the capture edge: an edge-cut box may now be absorbed by a receiver that owns its
band). One BLE ADV_IND (Apple Nearby Info) on channel 39 is decoded beside it. Python 0.1–0.3 s, C++ 40 ms per
block; parity identical (`testdata/live/LIVE_JAMMER_SWEEP_2460MHz`, `_jitter`). All 10 blocks of the session
(0.32 s apart) give one jammer each, 10000.13–10000.27 Hz, top 2469–2478 MHz; this needed two fixes found on them: a
fast transit smeared by the ~2 µs sweep jitter is not a flyback when the bins between the two peaks stay lit
(within 10 dB), and the trajectory is median-smoothed (a one-phase clipping product had set the top at 2484 MHz).

### Results on the 2026-10-04 Pluto sweep (713 snapshots)

Results are in `out/full/SWEEP_300-6000MHz_.../`. The sweep was last re-run for the 1.2, 2.4 and 5 GHz VTX bands
on 2026-10-05.

**Totals:** 4,090 signals at 6 dB SNR or more, plus 25,519 weak detections listed separately.
- 131 `success`, 242 `partial`, 3,717 `not_possible`.
- 349 of the signals are receiver artifacts: spurs, IQ images, LO-harmonic responses and LO leakage.
- 3,104 are `unknown`.

| Receiver | success | partial | What was found |
|---|---|---|---|
| Wi-Fi OFDM | 55 | 68 | FCS-valid frames, BSSID/SSID; HT/VHT headers |
| Wi-Fi 802.11b | 31 | 18 | beacons with SSID and DS channel (previous receiver; see below for the rebuilt one) |
| BLE advertising | 17 | 2 | CRC-24 valid adverts |
| 5G NR | 23 | 5 | PCI + MIB (n78 ×20, n28 ×8) |
| LTE | – | 75 | PCI, duplex, EARFCN, band (B1/B3/B5/B8/B28/B40/B41) |
| GSM | 1 | 65 | FCCH carriers / ARFCN; 1 SCH with BSIC |
| ADS-B | 4 | – | ICAO / callsign / altitude |
| RC links | – | 8 | single-burst "ExpressLRS FLRC/GFSK-like" at 2.4 GHz, confidence 0.6; **unconfirmed**, probably other 2.4 GHz traffic |
| Digital VTX | – | 1 | one unknown 44.7 kHz-SCS OFDM link at 2404.1 MHz, 16.6 dB SNR (unidentified) |
| Analog VTX, DJI, Herelink | – | – | none in this sweep (no drones were flying); all are validated on `testdata/` |

**ADS-B, live (2026-10-09):** 60 s at 1090 MHz from a balcony (8 MSps, gain 60; 33 s of samples) gave 129
CRC-valid frames from two aircraft, both fully summarised: IndiGo 6E311 (801648) climbing 9,700 → 11,700 ft out of
Hyderabad in a 24° left bank (31 positions, selected altitude 24,000 ft, ADS-B v2, NACp 10, TCAS), with Comm-B
replies to a Mode S radar giving IAS 222 kt, TAS 266 kt, Mach 0.408 and hence wind about 7 kt and SAT about +7 °C ±4,
and 18 ACAS replies; and 6E6294 (801903) at FL210 on track 028°. Indoors, the same set-up received nothing in 35 s.
In the sweep only 141036_221 holds a frame at its own LO; the other two earlier hits were stale buffers.

**Far DJI Mini 2 at 5770 MHz (2026-10-09):** on two detector-GUI captures (`detection-spectrogram-guis/data/5770_Ant2_Mini_2_*`,
drone 2 km away), the video downlink now matches the GUI's boxes: a 10 MHz link at 5760.49 MHz (171218) and a 20 MHz
link at 5780.48 MHz (171137). Before, the 10 MHz link hid behind its own failed 20 MHz candidate and was labelled
Herelink (DJI's 5 ms sync symbol passed the PSS purity check at 0.351 vs 0.35, now 0.5). Also, the 20 MHz link was
labelled 10 MHz because overlapping Wi-Fi cut the width measurement short; the mode now comes from the verified
ZC length. About 30 s per snapshot, about half of it in `digital_vtx`. On the 2.4 GHz set from the same flight
(72 of 155 snapshots run, about 62 s each with 2 workers), the controller uplink is found in most snapshots.

**Missed DJI links (2026-10-10):** two DJI video downlinks the `dji` receiver detected and then rejected are now
found. A Mini 2 at 500 m (D_20260515-192148, 5770 MHz): the carrier estimate was 2.35 MHz high because the scan's
best centre cut the band edge, so the ZC test failed; the estimate now re-centres until it stops moving (20 MHz
video at 5770.58 MHz, root 820, matching the labelled boxes). A 10 MHz link at 2466.5 MHz under Wi-Fi channel 11
(2484 MHz capture): its symbol grid moves between bursts, so the capture-level grid missed its ZC symbols; a
burst-by-burst re-timing fallback now finds it (root 116, 8 bursts). All testdata DJI results are unchanged; the
tv_2p4 video now survives +30 dB of added noise instead of +10 dB. C++ ported, parity identical.

**DJI bursts under overlapping signals (2026-10-10):** on D_20260515-192835 (Mini 2, 5770 MHz) Wi-Fi up to 14 dB
above the video cut its 7 labelled bursts into 18 pieces (70 % coverage). The `dji` receiver now tests each symbol
on the full band and on 4 sub-bands, bridges short gaps where another signal covers the whole band, takes the band
position from a per-bin median, and excludes in-band RC uplink hops: 7 bursts for 7 labels, none split, 93 %
coverage. Label scores elsewhere: tv_2p4 7/7 (97 %), D_20260515-192148 unchanged, and the saturated Mini 2 at 1 m
(tv_5p8_125115), never found before, now reports its video at 5770.54 MHz. No false alarm on noise or Herelink; the
Mini 2 still holds at +5 dB added noise. C++ ported, parity identical.

**Weak DJI O4 link on one antenna (2026-10-10):** D_20260609-151042 (FMComms5, 4 channels, drone 2 km away):
the video was found on ch1–3 but dropped on ch0, ~5 dB weaker, where the O4 repeated-reference test scored 6.35
< 8. The test now scores the two slot positions that carry the repeated pilot pair (about 2x the evidence on
every O4 link: ch0 9.5, ch3 26.5), low-SNR bursts are completed from the 2.5 ms TDD frame, and links with
evidence 5–8 are reported as `DJI OcuSync …-like` instead of dropped. All four channels now give the same 6
bursts. Note: the metadata says Mini 2 but the link is O4-type (30 kHz, repeated reference), unlike every other
Mini 2 capture (15 kHz, ZC sync).

**Any drone, found blind (2026-10-12):** a vendor-independent CP-OFDM link detector (`digital_vtx/blind.py`)
replaces the energy-based generic search. It needs no protocol: every drone video link in `testdata/` and `tmpdata/`
is found without its own receiver, with the right numerology, carrier and band — DJI OcuSync 2 (15 kHz, 10/20 MHz),
O4 (30 kHz, 10/20 MHz, 2.5 ms TDD), Herelink (LTE DL, 17.9 MHz at 2400.70 MHz in capture coordinates; its PSS gives 2400.63), and the O4 link on
D_20260609-151042 ch0 ~5 dB under the noise. Its sensitivity equals the DJI receiver's (Mini 2 and Mini 4 Pro to
+5 dB added noise, both lost at +10 dB). It also found a link the library had never reported: on the **DJI Phantom 4**
capture a 10.94 kHz CP-OFDM link (Tu 91.43 µs, CP 1/8: the 802.16e/WiMAX 10 MHz numerology), 10.3 MHz wide at
2426.1 MHz, in long TDD bursts (61 % duty) — very likely the Phantom 4's Lightbridge 2 video downlink, reported as
`OFDM drone video link-like (DJI Lightbridge / WiMAX-like numerology)`; and the 43.75 kHz fixed-wireless link under
the Mini 2 video on D_20260515-192835, labelled `OFDM fixed-wireless link-like (Cambium PMP 450-like)`. No link on
noise, Wi-Fi, the analog VTX, or (outside cellular bands) any of the live / sweep captures; every other testdata
report is unchanged. Links a protocol receiver identifies (DJI, Herelink) are attributed to it
(`explained_receiver_detections`). On D_20260609-151042 the `dji` receiver now also follows a link that re-times
mid-capture (the O4 grid jumps by 10.3 µs at 16.8 ms, and the frame phase with it): before, only the 6 bursts after
the jump were reported; now all 14 (ch3), and the O4 bandwidth is measured through the pilot dips (18.0 MHz, labelled
20 MHz; it was 7–11 MHz and labelled 10 MHz). C++ ported, parity identical.

**Complete bands, whole bursts and the hopping uplink (2026-10-13):** on D_20260609-151042 ch0 the O4 video was still
reported 14.1 MHz wide (true 18.1 MHz: 5761.5–5779.6 MHz, as the capture tool's labels on ch1), with bursts split in two
(6.93–7.68 + 8.25–8.60 ms) or starting late, and the drone's hopping RC uplink was not reported at all. The per-symbol CP
test is weak at ~0 dB; in-band power over 18 MHz is far stronger (every 36 µs symbol ~17 sigma). New `core/bursts.py`:
the band is the on-vs-off power ratio's best contiguous rectangle (multipath dips no longer end it), and bursts are the
in-band power runs (25th-percentile rank over sub-bands, so Wi-Fi over part of the band neither creates nor cuts a
burst; a much stronger signal over the band is "unknown" and bridged) that the receiver's CP evidence touches. Both
`dji` and the blind detector use it: all four channels now give the same 13 bursts, one per 2.5 ms frame, edges within
~50 µs of each other, band 5761.5–5779.6 MHz. The uplink: narrowband CP-OFDM hops with the video's numerology (30 kHz,
CP 72/1024), 4.4 MHz × 0.41 ms, a one-symbol preamble then 10 symbols, at 5747.5 / 5762.4 / 5792.5 MHz on a 2.505 ms
cadence that re-times with the video at 16.8 ms — the same hops on ch0–3; found blind by `blind_hops` (5–6 per
channel). The same detector finds DJI Mini 2's OcuSync 2 uplink (15 kHz, 1.1 MHz × 0.51 ms) on D_20260916 without the
DJI receiver (4 of 5 hops; the fifth sits on the LO tone), O4 uplink hops at ~2510 MHz on the 2484 MHz test captures,
and the PMP 450-like subscriber bursts at 5758.7 MHz on D_20260515-192148 (labelled fixed-wireless). No hops on noise,
Herelink, the Mini 2 / Mini 4 Pro / Phantom 4 video captures, Wi-Fi or the analog VTX. C++ ported, parity identical.

**802.11 at non-standard clocks and 802.11n LDPC/STBC (2026-10-07):** a sample of 80 sweep snapshots across
300–6000 MHz found one link the channel-plan scan cannot see: a MikroTik (Routerboard.com) station
6c:3b:6b:3c:d8:e5 sending 50-byte data frames every 0.35 ms on an off-grid 20 MHz channel at 5101.4 MHz
(snapshot 141319_273, 63 FCS-valid frames); no false link elsewhere (LTE, NR, GSM, DVB-T, standard Wi-Fi). The
11 real 802.11n MCS1 PPDUs at 5240 MHz (141325_432) that never decoded now give 8 FCS-valid MPDUs: the HT pilot
pattern rotates every symbol, and the decoder and its reference transmitter had both used a fixed one. The
sweep's `out/full` reports predate both.

**802.11b receiver rebuilt (2026-10-05):** compared on the fresh 2.4 GHz snapshots (same IQ):

| | PLCP headers | FCS-valid frames |
|---|---|---|
| previous receiver | 41 | 23 |
| this receiver | 74 | 32 |

No frame the old receiver decoded is lost. On the Phantom 4 test capture: 10 → 16 FCS-valid frames. It now
also decodes the 2 Mb/s, CCK and short-preamble frames the old one could not, down to 0–4 dB in-band SNR on
synthetic signals (`receivers/wifi/README.md`). The sweep's `out/full` reports predate this; re-run the 2.4 GHz
band to refresh them.

**Notes:**
- `receiver ppm` is +7.44, from the LTE PSS.
- These results were produced while the pipeline still corrected the sweep's stale-LO snapshots (3 of the 4 repeats
  after each retune held the previous step's IQ). That correction has been removed: the library now trusts the
  capture metadata, so this sweep should be recaptured with the buffers flushed after each retune before it is
  re-analysed.
- Reports that were not re-run since the 2026-10-05 restructuring still name the old research modules in
  `decoded.method` (for example `det_rc_link`); their results are otherwise identical. A regression check on 18
  captures of all types showed 0 differing reports.

## Viewing results in signal-vault

signal-vault (`../signal-vault`, the acquisition app) has a **Signal Intelligence** tab that shows what this library
detects and decodes, so every change here can be checked by eye:
- the analysed block's spectrogram, with a box per signal: green decoded, yellow partial, red not decoded,
  grey artifact;
- a table of signals, with label, status, frequency, bandwidth, time, SNR, confidence and what was decoded;
- the full JSON of the selected signal: detection method, evidence, `decoded.steps` and fields;
- **▶ View feed** for an analog FPV video signal, which opens a window with the demodulated live picture, in colour
  when the VTX sends a colour burst (the window says `colour` or `mono`);
- **🔊 Listen** for an FM voice signal (walkie-talkie, PMR): its demodulated audio on the speakers, live from every
  new block or played through the capture and the following ones at the same frequency.

Sources:
- **Live:** "Analyse the live stream" or "Analyse next block" on whatever the radio delivers.
- **Replay:** Device `replay`, or `sv gui --replay <capture or folder>`. Recorded SigMF files (a capture, a recording
  session or a whole sweep) are fed through the engine at real time, as if they came from a radio.
- **Files:** analyse a capture; list a folder; or open an existing `*_report.json` from `spectrum-intel run`. With
  a reports folder set, a listed capture shows its existing report instead of being re-analysed.

```
.venv/bin/sv gui                                   # this repo's .venv has signal-vault installed (editable)
.venv/bin/sv gui --replay testdata/Analog_20260804-161815
```

The tab only draws. signal-vault runs the library in a worker process through `signal_vault/intel/`. It always runs
the code in `python/` as it is on disk: the package is installed editable from `python/` (or, if not installed,
`../spectrum_intelligence/python` is put on the path), and the bridge fingerprints the library's source before
each analysis. After an edit it starts a fresh worker and re-imports the video/audio demodulators, so no GUI
restart is needed. It uses this library's public API:
- `spectrum_intel.analyze.analyze_iq(iq, fs, fc, ...)` and `analyze(path, ...)`;
- `report.writer.build(analysis)`, which gives the report dict;
- `decode.registry.receiver_names()` and `plugins(enable)`;
- `receivers.analog_vtx.FieldDemodulator.from_fields(signal["decoded"]["fields"]).process(iq, fs, fc)`
  (`VideoFrame.image` luma; `VideoFrame.rgb` colour or None);
- `receivers.voice_fm.VoiceDemodulator.from_fields(signal["decoded"]["fields"], audio_fs).process(iq, fs, fc)`,
  offered for every signal whose `decoded.fields.audio` is `"nbfm"`.

**Keep these stable.** Anything new that a receiver puts in the report (labels, `decoded.fields`, `detection`) shows
up in the tab without any change to signal-vault.

## Architecture

```
python/spectrum_intel/
  api/             production API: Engine, IQBlock, StreamAnalyzer, response builder; schema -> ../../../schema
  core/            shared, numpy/scipy only: channelizers, spur finder, 15 kHz OFDM scan, Emission record,
                   SigMF loaders, 3GPP sequences, India band plan, thread helpers (parallel.py)
  receivers/       one package per signal family (see table); DSP imports only core
    <name>/__init__.py   public API          <name>/plugin.py   pipeline adapter (PLUGINS)
    <name>/README.md     what / how / API / validation
  decode/          plugin API (base.py), discovery (registry.py), Emission adapter
  detect/          energy detection: STFT, noise floor, CFAR, merge, artifacts
  features/        timing, envelope, cyclic, OFDM, FSK features
  classify/        India band rules for undecoded signals
  analyze.py       one snapshot: detect -> features -> plugins -> classify
  sweep.py         sweep context: noise template from neighbouring tunings, receiver ppm
  report/ viz/     JSON writer, annotated spectrogram
  summary.py       sweep summary, spur / LO-harmonic checks, occupancy map
  research/        standalone research CLIs (spectrum_survey, drone_survey); not used by the pipeline
  cli.py
python/tests/      unit, modularity, API and receiver tests
c++/               the same library in C++ (c++/README.md, c++/PORTING.md)
schema/            response.schema.json, technologies.json, requirements.json (shared by both)
```

The pipeline for one snapshot (`analyze.py`):
1. **Detect:**
   - coarse STFT (60 kHz × 16.7 µs) and fine STFT (3.75 kHz × 267 µs);
   - noise floor from the sweep-neighbour template, or a morphological opening, with a per-frame gain g(t);
   - continuous detector, plus hysteresis CFAR for bursts;
   - merge, valley split and burst grouping;
   - flags for LO leakage, IQ images and clipping.
2. **Features** for every energy detection.
3. **Plugins:**
   - `capture_signals`: receiver-defined signals, which absorb the energy detections they explain. Receiver
     packages run concurrently, each on its own copy of the capture cache; a plugin that reads another package's
     cache entries names it in `after` (`nr` and `gsm` read the receiver ppm `lte` measures);
   - `arbitrate`;
   - per-signal `decode` for whatever is left;
   - `postprocess` (emitter association);
   - `frequency_label_check`: if a plugin proves the centre-frequency label wrong, the snapshot is re-analysed at
     the corrected frequency.
4. **Classify** the undecoded signals, then the **report** and **PNG**.

### Taking a capability out (example: analog VTX)
- **In another pipeline:** `pip install -e spectrum_intelligence/python`, then
  `from spectrum_intel.receivers.analog_vtx import detect, video_frame`. This loads only `core` and that receiver;
  `tests/test_modularity.py` checks this for every receiver.
- **As a copy:** take `python/spectrum_intel/core/` and `python/spectrum_intel/receivers/analog_vtx/`, and drop
  `plugin.py` if the pipeline isn't needed. In C++: `-DSI_RECEIVERS=...` builds a binary without it.
- **Out of this pipeline:** `--disable analog_vtx` for one run, or delete the folder. The registry finds receivers
  at run time, and nothing else imports them.

### Adding a receiver
1. Verify it on the bench first. Then create `python/spectrum_intel/receivers/<name>/` with the DSP modules
   (importing only `spectrum_intel.core`), an
   `__init__.py` re-exporting the API, a `README.md`, and a `plugin.py`.
2. In `plugin.py`, subclass `decode.base.Decoder`, or `decode.emission_adapter.EmissionReceiver` if `detect()`
   returns Emission records, and set `PLUGINS = [...]`.
3. Add the name to `decode.registry.ORDER` if the priority matters.
4. Add its sample-rate and duration needs to `schema/requirements.json`, and a technology code to
   `schema/technologies.json` if it defines a new kind of signal.
5. Run `pytest tests` (including `test_modularity`), then update this README's signal table.
6. Port it to `c++/src/receivers/<name>/` following `c++/PORTING.md`, and run the parity check
   (`c++/tools/parity.py --receivers <name>`).

## Tests

```
cd python && ../.venv/bin/python -m pytest -q tests     # unit + API + modularity + emitter association (~2 min)
cd python/tests/receivers && ../../../.venv/bin/python test_det_digital_vtx.py   # receiver validations (tests/receivers/README.md)
c++/build/si_tests && .venv/bin/python c++/tools/parity.py --cpp c++/build/si_analyze   # C++ unit tests, C++ = Python
```
Test captures with known contents are in `testdata/` (see `testdata/README.md`); `testdata/live/` holds the
2026-10-10 live Pluto captures (TETRA, DMR/NXDN, LTE B5/B28, receiver spurs, LO-harmonic images). For a change that
touches many receivers, `bench/regress.py run <tag>` / `diff <a> <b>` compares the Python reports on all testdata.

## Known limitations / next steps
- **Jammers:** swept-carrier jammers (one sweep per block, ≥ ~+4 dB per MHz over the noise and not much weaker
  than the random traffic in its bins), GNSS swept jammers (`gnss`), and barrage-noise / noise-FM jammers (research,
  synthetic only; a Gaussian plateau < 10 MHz is called a jammer only on an allocation's edges, since a WCDMA carrier
  looks the same). Spot and hopping jammers are not identified, and an FMCW radar outside the radio-altimeter band
  would be called a jammer.
- **Signal-gaps batch (2026-10-25, all *research*):** every new capability (Mode A/C, AM voice, APRS, ILS / VOR,
  FM + RDS, AIS, GNSS, BLE trackers, BR/EDR LAP, 315 / 433 MHz remotes + TPMS, radio altimeters, noise jammers, radar
  classes and deinterleaving) is validated on synthetic signals and published test vectors only; none has a real
  capture in `testdata/live` yet. The band plan (`core/bands_india`) now starts at 70 MHz, but the VHF items (ILS
  localizer, VOR, FM, AIS, APRS, civil air band) need a receiver that tunes below 300 MHz. Not done: satellite-phone bursts
  (Iridium, Thuraya, Inmarsat), ACARS / VDL2 / DSC, phone uplinks and fake-base-station indicators (blocked or on hold,
  see `research/signal_gaps/PLAN.md`), spot jammers, more TPMS families, the BR/EDR UAP and an EDR detector.
- **Licences of the Python dependencies:** the Python package uses pyModeS (GPL-3) for ADS-B and py3gpp (AGPL-3)
  for the NR polar decoder. Shipping the Python package therefore carries those licences. The C++ library has
  neither: its ADS-B and NR code is written from the standards (DO-260B / ICAO, 3GPP TS 38.211/38.212), and those
  packages serve only as test oracles. Its third-party code is pocketfft (BSD-3), nlohmann/json (MIT) and doctest
  (MIT).
- **Input assumptions:** the library trusts the capture metadata. Each snapshot's centre frequency must be the LO
  it was captured at, so acquisition must flush the receiver buffers after every retune. Input problems are fixed
  at the input side, not in this library.
- **Acquisition (Pluto sweeps):**
  - Clipping in the cellular bands.
  - 32.5 ms snapshots are too short for a full GSM multiframe, LTE MIB repeats and voice; a TETRA BSCH is decoded
    only when a synchronisation burst falls in the block, an NXDN frame sync only when its 40 ms frame does.
  - No preselector: strong 1.8–2.6 GHz cellular downlinks appear at 330–830 MHz through the mixer's 3rd and 5th
    LO harmonics. One block only gets a `possible_lo_harmonic_image` caution (which also lands on some real
    signals); a filter per sub-band, or a second tuning, settles it.
  - Capture self-noise from the USB link (floor +4–5 dB ~12 ms into each buffer, 8 kHz microframe modulation, combs
    near 420 and 481–492 MHz); fix at the acquisition side. Pluto Zynq clock harmonics (33.333 MHz, DDR 533.33 MHz)
    are not on the sample-clock / reference grid and are not flagged.
  - The LO-harmonic images reach above 1 GHz too: n78 (3.4–3.6 GHz) at 1128–1220 MHz and an in-room 5 GHz Wi-Fi
    AP at 1708–1758 MHz (3rd harmonic). The band-plan caution only knows cellular sources, so the Wi-Fi images are
    not flagged.
- **Analog video colour:** only the streaming `FieldDemodulator` (the GUI feed) decodes colour; the report's
  `<name>_vtx<MHz>_video.png` is still luma (a diagnostic stack of both fields in time order, so it starts mid-field).
  Colour is checked on synthetic bars, one real NTSC VTX and one real PAL VTX (`LIVE_ANALOG_VTX_PAL_5740MHz`, 2026-10-22:
  field-aligned colour picture, burst field 5.2). At low SNR the chroma is noisy before the luma is (FM noise rises with frequency); the colour
  killer then falls back to grey. On a full-screen colour-bar test pattern the detector can mislocate the VTX
  (polarity / sync tip; seen on synthetic bars only), which spoils the luma as well. The streaming demodulator's line-rate
  refit (±40 ppm around the detector's rate) can run to its search edge on a noisy block; the sync depth is then
  under-read and the colour picture comes out too bright.
- **Single-block noise floor:** without a sweep noise model the floor is a 3 MHz morphological opening; a wide
  always-on carrier (≥ 10 dB above the band's quiet level) is bridged by interpolating the floor from the quiet bins
  on either side (live: LTE B3 PCI 101, 20 MHz, now detected and decoded). A wide carrier less than 10 dB up, or a
  band with no quiet bins, is still read as floor.
- **GSM:** the SCH still needs an MLSE equaliser (training sequences match 57–64 of 64 but most CRCs fail); the
  uplink (handset bursts) is not detected.
- **Not decoded yet:**
  - LTE MIB/SIB1, which would give the operator (PLMN);
  - GSM SCH for most carriers, which needs an MLSE equaliser;
  - NR SIB1;
  - UMTS scrambling code;
  - DVB-T2 L1;
  - TETRA MAC/MLE broadcast (SYSINFO: main carrier, location area); DMR link control / voice (AMBE+2); NXDN
    LICH/SACCH;
  - Wi-Fi payloads that need more than one antenna or more bandwidth: VHT/HE, HT with 2+ spatial streams, 40 MHz
    HT (each is reported with the reason, not skipped silently; their transmitter is still attributed when the
    response or scrambler-counter evidence allows). Single-stream 40 MHz HT and VHT 20 / 40 MHz would fit a Pluto
    capture but no testdata holds any (2026-10-25: every real VHT PPDU is 80 MHz or 2-3 streams), so PLAN W5 waits
    for such a capture; HE (802.11ax) payloads are not decoded (their HE-SIG-A is).
  - Drone-over-Wi-Fi by behaviour (PLAN W1) looks at one block: an IP camera or a phone hotspot streaming video to
    one client has the same shape, and the video frame rate needs ~1 s of blocks; the behaviour verdict stays
    "-like" at 0.55 and says so. Identity (Remote ID, Tello SDK, MAVLink, signatures) is what names a drone.
- **DJI:** RC uplink hops inside the video band are missed (D_20260916-171218: 2 of 4 hops, at 5764 MHz). A burst
  edge next to a much stronger signal over the whole band is uncertain: it can run ~0.2 ms into it (D_20260609 ch1
  14.13 ms, ch3 3.46 ms); a link too weak to show in in-band power (≲ −3 dB, Mini 2 + 5 dB noise) keeps the CP-evidence
  bursts, which can still fragment. The saturated Mini 2 at 1 m (tv_5p8_125115) is found at 0 dB but not reliably with added
  noise. On D_20260515-192148 one of 3
  labelled uplink hops is missed and the two found are about 2 MHz above their label boxes (unexplained).
- **Blind OFDM links (`digital_vtx`):** a PHY alone cannot tell a drone from a fixed-wireless link of the same
  numerology; the label says what the numerology resembles (`-like`), and an unknown numerology is `unknown OFDM
  (TDD) link`. Not found blind: non-OFDM or CP-less video (analog FM has its own receiver; HDZero's 1/32 CP is too short
  for the blind scan, so the `hdzero` path finds it by its pilots instead, 27 MHz mode only), links with Tu < 10 µs (802.11-derived links are the Wi-Fi receivers'), non-OFDM hopping control links
  (FHSS FSK / LoRa: `rc_link`), hops of a numerology that is neither in the capture's links nor a drone family's
  (15 / 30 kHz), a hop hidden under a much stronger burst (D_20260609 ch0 29.6 ms), links weaker than ~+5 dB below the DJI Mini 2 test capture, and links outside
  2.33–2.53 / 4.9–6.0 GHz (the plugin's range). The Lightbridge identity of the Phantom 4 link is inferred, not
  verified. The frame period of very long bursts (Phantom 4) is not meaningful.
- **Transmit antennas (`digital_vtx` `tx_antennas`, 2026-10-21):** evidence, not identification. Validated on
  synthetic links and on real 802.11 with SIG-stated modes (Microhard STBC: kurtosis 1.60, pair delta −0.33;
  FCS-verified single-stream frames: 1.07–1.27, ~0), and found on one real drone link (LIVE_OFDM_VTX_5805MHz,
  Walksnail-like: `2tx_stbc_like`, delta −0.18). It has **not** been checked against Walksnail units of known
  antenna count, so it does not name a model. One stream per subcarrier also covers several antennas with cyclic delay
  or switching. Spatial multiplexing, varying allocations and non-QAM data give `unknown`. Links below 20 dB SNR,
  under 200 subcarriers, or without 20 contiguous full-band data symbols (DJI O4's interleaved comb-2 reference
  symbols break the runs) give `unknown` with the reason.
- **802.11b:** most FCS failures left on real captures are collisions with a co-channel OFDM frame, which no
  single-receiver excision can undo. PPDUs that began before the capture are reported as DSSS-like, not decoded.
- **Wi-Fi at 5 GHz:** most data traffic is 80 MHz VHT, which does not fit in a 61.44 MSps Pluto capture (every VHT
  PPDU in the live 2026-10-05 captures was 80 MHz). What is decodable is the non-HT management and control traffic,
  which carries the MAC identities, plus the VHT headers. A 16-/64-QAM bit-mapping bug that failed every 24–54 Mb/s
  frame was fixed on 2026-10-05 (live 5 GHz: 39 → 70 FCS-valid frames); the receiver now decodes every frame
  gr-ieee802-11 decodes on the same IQ, plus 16 more. Capture 5 GHz at about 40 dB gain: 50 dB rails the ADC near
  an AP and costs the 16-QAM frames.
- **802.11n (HT) payloads:** single-stream BCC is confirmed on real frames (8 of 11 MCS1 PPDUs at 5240 MHz in an
  ADC-clipping band pass FCS) and LDPC + STBC on the Microhard link (MCS7, 3 PPDUs with CRC-valid MPDUs; most of its
  PPDUs arrive at about 16 dB SINR, below what MCS7 needs, so the rest are attributed by evidence, not decoded).
  2.4 GHz client traffic at MCS 3–7 has not been tested on real IQ yet.
- **802.11 at non-standard clocks:** needs some burst activity (≥ 150 µs) with sharp spectral edges, a 1/4 guard
  interval (802.11's long guard) and a clock of 2.5–26 MHz; links with only short-guard HT data and no legacy
  preamble, 40 MHz-type links (two 20 MHz halves) and continuous (non-bursty) OFDM are not searched. Validated on
  synthetic 5/10/11 MHz links (mirrored or not), the real Microhard (11 MHz, mirrored) and DJI (5 MHz) links and an
  off-grid MikroTik channel in the sweep; it adds about 1 s per 61.44 MSps snapshot.
- **Remote ID is self-declared**, broadcast in clear and unauthenticated: the values are what the aircraft chose
  to say, and an aircraft that does not broadcast it is invisible to this. Platform signatures (vendor OUI, SSID)
  identify equipment, never a flight.
- **Devices/findings:** range from signal level only — locating an emitter needs several receivers or a
  directional antenna. A 33 ms snapshot sees roughly one beacon in three, so device counts are lower bounds and
  the finding thresholds are per capture, not per second.
- **ESB-like 2 Mb/s links (`rc_link`):** 2 Mb/s only (1 Mb/s and 250 kb/s ESB are not sliced), CRC-16 only (an
  8-bit CRC is not searched), and a link needs ≥ 3 packets of distinct content in the block for a CRC verdict
  (fewer stay `partial`). The vendor and the payload fields are not identified (the bench link is idle traffic of
  an unknown product). Packets inside a Wi-Fi channel are partly missed by the narrowband scan. One block sees only
  the channels used during its 34 ms, so the hop set needs several blocks. The energy detector still adds an
  `unknown` entry around a channel shared with BLE adverts (2426 MHz), as for BLE alone.
- **LoRa:** payloads are decoded for explicit-header packets, and for SF5 implicit-header packets with a CRC (the
  CRC fixes length and coding rate; with a last block of zero nibbles the CR stays ambiguous: `cr_ambiguous`).
  Explicit-header SF5 and all SF6 frames are not verified on hardware and not decoded. Other implicit-header links
  (ExpressLRS, most RC) carry no length or coding rate on air and stay at SF/BW/length. ExpressLRS 2.4 GHz uses the
  SX1280 long-interleaving coding rates (CR_LI 4/6, 4/8): their interleaver is not published, and the standard
  chain fails even on the first block, so 2.4 GHz ELRS payloads need a reverse-engineering capture (many packets of
  varying content) first. SF5 adds 0.1–0.3 s per 34 ms sub-GHz block to the LoRa scan (Python). Captures that touch no RC band but
  860–930 MHz now run the LoRa scan on that range only (no FSK: GSM/LTE uplink bursts would look like RC links),
  (rc_link-only, interleaved A/B: unchanged on `LIVE_868MHz` (8.1–8.9 s Python); `LIVE_LORA_SF5_890MHz` 2.3 s
  Python; the 1.75 s `LIVE_PMR_860.4MHz` dwell, newly scanned, ~7 s Python / C++ whole report 5.7 s). SF5 stage 1
  uses eta 40 (24 let ~170 DMR/NXDN false candidates per dwell into refinement) and refinement shares each
  sub-band's full-capture slice (200 -> 87 slices on the dwell, bit-identical). LoRaWAN and Meshtastic payloads are encrypted; only
  their plain framing is read. SF12 needs about 0.3 s of capture to see its preamble, so it is invisible in the
  sweep's 32.5 ms snapshots; record 433/865 MHz at a low sample rate (2.5 MSps) for about a second instead.
- **Pluto at 2.5 MSps (2026-10-05, 433.92 MHz):** every 0.8 s block holds the same artifacts at the same offsets:
  wideband noise gated 23 ms on / 12 ms off for the first 265 ms, a 65 ms stretch about 12 dB lower, then comb
  lines. They are tied to the block, not to the air, so they come from acquisition (cause not yet found; note that
  the RF bandwidth stays at 61.44 MHz while the sample rate is 2.5 MSps) and give thousands of weak narrowband detections per block, none reported as a signal.
- **Walkie-talkies (`voice_fm`):** analog FM only. Scrambled voice (inversion, rolling-code) and digital voice
  (DMR, dPMR, NXDN, P25, TETRA) are not turned into audio: a scrambled carrier is still found, with its channel and
  squelch code, but plays unintelligible, and digital voice falls to the band rules. DCS and CTCSS need ≥ 0.35 s
  and ≥ 0.3 s keyed, so 32.5 ms sweep snapshots give the carrier, channel and deviation only. Validated on one real
  PMR446 radio (channel 8, D053N, 25 blocks) that was keyed with no speech and close enough to clip the ADC; the
  code is not yet confirmed against the radio's display, and a capture with speech is still needed. Keep a handheld
  a few metres from the Pluto and use 0–10 dB gain.
- **ADS-B:** a position needs an even/odd pair within 10 s (or the receiver position in SigMF `core:geolocation`),
  so the sweep's 32.5 ms snapshots give identity and altitude only; record 1090 MHz for seconds (8 MSps) for
  tracks. Comm-B registers (wind, temperature, IAS/Mach) appear only while a Mode S radar interrogates the aircraft.
  Derived wind uses the magnetic heading as true (India's declination is 0–2°) and is good to about ±5–10 kt; SAT
  carries the quantisation error it reports (±3–5 °C). Reception needs sky view: indoors gave no frames. The
  country is the ICAO address block and the operator the callsign prefix: what the transponder is set to say.
- **Unknowns:** many "unknown" detections are fragments of known emitters. Cross-snapshot grouping and an ML
  classifier for the remaining unknowns are planned.
