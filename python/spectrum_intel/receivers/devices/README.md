# devices — who is transmitting, what they are doing, and what deserves attention

This receiver demodulates nothing. It reads the 802.11 frames the other receivers decoded and turns
them into the two things an operator actually asks for: **a list of the emitters present** and **a
list of what is worth looking at**. One frame identifies a radio; the same address seen across a
capture gives a pattern of life, and that is what encryption does not hide.

It imports no receiver. Frames reach it through the generic report contract
(`signal.decoded["fields"]["ppdus"]`), so 802.11b and OFDM feed it the same way and anything else that
yields addressed frames (BLE, later) can too.

## What it produces

**`report["devices"]`** — one row per transmitter heard. An address only ever seen as a *destination*
is not a device we heard; it is recorded as a peer.
- identity: `mac`, `vendor` (IEEE registry), `address` (global, or randomised/locally administered),
  `role` (access point or station), `device_name` (WPS / Wi-Fi Direct);
- presence: `first_seen_s`, `last_seen_s`, `n_frames`, `snr_db_median`, `snr_db_max`, `channels`;
- activity: `airtime_ms`, `bytes_seen`, `frame_types`, `rates_mbps`, `retry_fraction`,
  `power_save_fraction`, `encryption_packet_counter_max` (how much the device has sent since it got
  its key), `talks_to`;
- for an access point: `ssid` (or `<hidden>`), `security`, `ap_uptime`, `country`, `bss_load`, `wps`;
- for a station: `associated_bssid`, `probed_ssids` (networks it has joined before), and
  `probe_fingerprint`;
- `drone` / `remote_id` when the frame carried them;
- `seq_first` / `seq_last` — the 802.11 sequence number of its first and last frame (PLAN W7);
- `cfo_hz_median` — the transmitter's crystal offset, which is a property of the hardware.

**`report["device_groups"]`** — randomised MACs that are probably one handset. Phones rotate their
address, so counting addresses overstates how many devices are present. Two things survive a MAC
change: the probe-request element signature (fixed by the driver and OS) and the transmitter's clock
offset. Agreement on both is evidence, not proof, so each group carries a confidence and its reason.
Many drivers also keep one 12-bit sequence counter across the address change (Vanhoef et al., AsiaCCS 2016): when
the next address's first sequence number continues the previous one's last within 64, the group says
`sequence_continues` and its confidence rises to 0.9 (research, PLAN W7).

**`report["findings"]`** — ranked critical / high / medium / low / info. Every finding says what was
observed, what it means, and how sure it is.

| kind | severity | fires on |
|---|---|---|
| `drone_remote_id` | high | an aircraft broadcasting ASTM F3411 or DJI Remote ID — carries the UAS ID, its position and **the operator's position** |
| `drone_signature` | medium | a drone vendor OUI or a drone-shaped SSID, with no Remote ID: identifies equipment, not a flight |
| `drone_link_behaviour` | medium | a Wi-Fi signal's `drone_links` entry with basis "behaviour" (receivers/wifi/field.py): one client receiving a steady heavy downlink and sending small frames, with no identity; says that an IP camera or a phone hotspot looks the same (research, PLAN W1) |
| `beacon_flood` | high | ≥ 15 beaconing access points without an IEEE OUI, from ≥ 8 unrelated address prefixes: a beacon-flooding tool (virtual BSSIDs of one router share a prefix and do not fire it; research, PLAN W6) |
| `deauth_flood` | critical / high | ≥ 25 (or broadcast-addressed) deauthenticate/disassociate frames from one transmitter; ≥ 8 is "high" |
| `duplicate_ssid` | high / medium | one SSID on several radios with different security or vendor — an impersonating AP looks like this, and so does a mesh; the telling sign is a copy offering weaker security |
| `multi_bssid` | info | the same SSID on BSSIDs that differ only in the locally-administered bit and a byte: one router's virtual BSSs, **not** an attack |
| `weak_encryption` | medium | WEP / TKIP / WPA1 with no modern cipher |
| `legacy_cipher_offered` | low | TKIP offered *alongside* CCMP — a compatibility mode, not a broken network |
| `open_network` | low | no link-layer encryption |
| `wps_active` | medium | the AP is in WPS enrolment mode |
| `key_exchange` | info | a WPA handshake: marks exactly when a device associated |
| `mac_rotation` | low | several addresses that look like one device |
| `client_history` | info | the networks a station is searching for — where it has been |
| `named_device` | info | a device that broadcasts its own name |
| `inventory` | info | the count and the named networks |

Rules that cannot separate an attack from ordinary behaviour say so and sit at a lower severity. The
three cases that would otherwise cry wolf are tested for explicitly: virtual BSSIDs on one radio, TKIP
beside CCMP, and a single roaming deauthentication.

## Across a sweep

Snapshots are analysed in separate processes, so `summary.py` is the only place an emitter can be
followed from one snapshot to the next. It merges devices and findings by MAC into
`summary.json`: `devices`, `device_summary`, `findings`, each with `n_snapshots` and `seen_at_fc_mhz`.
A device heard in many snapshots is a fixture of the area; one heard in a single snapshot was passing
through or is at the edge of range — which is itself worth knowing, so it is reported, not filtered.

## API (stdlib only — no numpy, no receiver imports)
```python
from spectrum_intel.receivers.devices import frame, tracks, link_randomized, findings
frames = [frame(mac_dict, t0=..., dur_s=..., snr_db=..., cfo_hz=..., channel=6)]
devs   = tracks(frames)                 # one row per transmitter, strongest first
groups = link_randomized(devs)          # randomised MACs that are probably one device
found  = findings(devs, frames, groups, drone_links)  # ranked; drone_links: the Wi-Fi signals' fields.drone_links
```

## Limitations
- **Range, not position.** Signal level gives a rough distance only. Locating an emitter needs several
  receivers or a directional antenna; one Pluto cannot do it.
- **Only what was decoded.** A device whose frames all failed FCS is absent. 802.11b and non-HT OFDM
  decode well; 2+ stream and LDPC payloads do not, so busy MIMO clients are under-represented.
- **One capture is 33 ms.** A beacon interval is about 102 ms, so a snapshot sees roughly one beacon in
  three. Counts of devices are lower bounds, and `n_frames` is not a traffic rate.
- **Thresholds are per capture**, not per second. The deauthentication thresholds were chosen for a
  33 ms snapshot; a longer dwell would need them scaled.
- MAC-rotation linking needs probe requests. A station that only sends data frames cannot be linked
  this way.

## Validation
`tests/receivers/test_det_devices.py` (pytest, ~1 s): roles and fields; MAC-rotation linking by
fingerprint and clock, including a device that must *not* be linked; a deauthentication flood raised
and a single one not; virtual BSSIDs distinguished from an evil twin; security findings proportionate;
a drone finding carrying the operator position; findings sorted and self-describing; and the plugin
driven purely through the generic report contract for both 802.11b and OFDM field names.
