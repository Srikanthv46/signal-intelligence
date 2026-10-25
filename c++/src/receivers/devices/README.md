# devices (C++) — device tracking and findings over decoded 802.11 frames

Port of `python/spectrum_intel/receivers/devices` (plugin.py, tracker.py, threats.py). Not a receiver: it
decodes no IQ. It reads the frames other receivers decoded through the generic report contract
(`signal.decoded.fields.ppdus`, each PPDU with `mac` or, for an A-MPDU, `macs`) and imports no receiver.

## What it does
- **Frames**: every FCS-verified PPDU (`fcs_ok` not false) of every receiver signal becomes one frame per MPDU;
  802.11b records (`start_s`/`end_s`/`cfo_hz`) and OFDM records (`start_us`/`dur_us`/`cfo_khz`) are normalised;
  the source is the first word of the signal's decoder method.
- **Tracker**: one device per transmitting address (TA/SA/addr2, or the BSSID of an AP-only subtype), with
  role, vendor, first/last seen, airtime, bytes, SNR/CFO medians, channels, frame types, rates, retry and
  power-save fractions, packet-number max, peers, probed SSIDs, probe-request fingerprint, drone / Remote ID,
  AP fields (SSID, security, uptime, country, BSS load, WPS); stations are tied to the BSS they talk to.
  Strongest first.
- **Randomised MAC linking**: same probe fingerprint and CFO medians within 2.5 kHz.
- **Findings** (most serious first): drone Remote ID / drone signature, deauth bursts, duplicate SSID /
  multi-BSSID, open / weak / legacy-cipher networks, active WPS, WPA key exchange, MAC rotation, client probe
  history, self-named devices, inventory.
- Plugin `devices` (package `devices`): capture-level, defines no signals; `postprocess` adds `devices`,
  `findings`, `device_summary` (and `device_groups`) to the report and `device_macs` to each signal's decoded
  fields; `finalize_report` adds `signal_ids` to devices and findings.

## C++ API (`namespace si::rx::devices`, depends on core/ only)
```cpp
#include "receivers/devices/devices.hpp"
std::vector<Frame> fr = frames_of(signal.decoded.fields, "receivers.wifi");
Json devs = tracks(fr);                     // tracker.hpp
Json groups = link_randomized(devs);
Json found = findings(devs, fr, groups);    // threats.hpp
Json extras = device_report(signals);       // the whole postprocess (sets device_macs)
finalize(report);                           // the finalize_report hook
```
`pyjson.hpp` (`detail::`) holds the Python dict semantics the port needs: truthiness, `a or b`,
`collections.Counter` (insertion order, stable `most_common`, `1 == 1.0` keys), and the key strings
`json.dumps` writes for numeric dict keys (`{"24.0": 5}` for a float rate, `{"36": 1}` for an int channel).

## Exactness notes
- Values copied from the frames keep their JSON type, so `round()` of an int SNR or CFO median stays an int
  and dict keys and fingerprint strings read as Python writes them.
- A signal whose decoder method is empty makes Python's `method.split()[0]` raise; the C++ throws the same
  `IndexError: list index out of range` (recorded as `devices_postprocess_error`).

## Tests and parity
- `tests/unit/test_rx_devices.cpp` against `tests/vectors/devices/cases.json` (written by `gen_devices.py`
  from the Python package): 4 synthetic cases covering every tracker path and every rule; extras,
  per-signal `device_macs` and `finalize_report` must be identical.
- `parity.py --receivers devices`: 16/16 inputs identical (no decoded frames without Wi-Fi, so no extras).
- `parity.py --receivers devices,wifi`: 16/16 inputs identical; 8 of them carry devices and findings (2-12 devices,
  2-31 frames each), compared as report-level extras. wifi_dsss 14/14, wifi_ofdm 12/12, wifi_ofdm_scaled 3/3.
