# adsb — Mode S / ADS-B 1090 MHz

**Detects and decodes** Mode S replies and squitters, decodes every field each one carries, and summarises each
aircraft. A frame is kept only when its parity checks:
- DF17/18 need CRC residual 0.
- DF11 needs a residual below 128, which is the interrogator code.
- DF0/4/5/16/20/21 carry the address in their parity, so they are kept only when that address is confirmed by a
  DF11/17/18 frame in the same capture (the dump1090 rule).

Status is `success`.

## Per frame (`messages.py`)
| Format | Fields |
|---|---|
| DF17/18 TC1–4 | callsign, emitter category (light … heavy, rotorcraft, UAV, glider, …) |
| DF17/18 TC5–8 | surface speed and track, CPR |
| DF17/18 TC9–18 / 20–22 | baro / GNSS altitude, CPR (odd flag, lat/lon), surveillance status (alert, SPI), NIC-B |
| DF17/18 TC19 | ground speed + track, or IAS/TAS + heading; vertical rate and source, NACv, GNSS−baro difference |
| DF17/18 TC28 | emergency state (general, medical, min fuel, no comms, unlawful interference, downed) + squawk; TCAS RA broadcast |
| DF17/18 TC29 | selected altitude (MCP/FCU or FMS), baro setting, selected heading, autopilot/VNAV/LNAV/alt-hold/approach, NACp, SIL, NICbaro, TCAS operational |
| DF17/18 TC31 | ADS-B version, NIC-A, NACp, GVA, SIL (per hour/sample), SDA, TCAS/1090ES-IN/UAT-IN/ARV/TS capabilities, single antenna, heading reference, IDENT, TCAS RA active |
| DF17 / DF18 | transponder capability (CA) / control field (ADS-B non-transponder, TIS-B, ADS-R) |
| DF11 | capability, interrogator: II code (0 = acquisition squitter or II-0 radar) or SI code |
| DF0 / DF16 (ACAS) | vertical status, ACAS sensitivity level, reply info (ACAS type or max airspeed), cross-link, altitude; DF16 MV and RA reports |
| DF4/5/20/21 (radar) | flight status (airborne/ground, alert, SPI), downlink request, interrogator identifier (IIS), altitude or squawk |
| DF20/21 Comm-B | register inferred by pyModeS: BDS17 (capabilities), 20 (callsign), 40 (selected alt, baro), 44/45 (MRAR: wind, temperature, pressure, humidity, turbulence, icing, wind shear), 50 (roll, true track, GS, track rate, TAS), 60 (magnetic heading, IAS, Mach, baro and inertial vertical rate) |

Every frame also has its time, raw hex, parity residual and SNR (pulse amplitude over the envelope median).

## Per aircraft (`tracks.py`)
`tracks(msgs, ref=None)` groups frames by ICAO address. It needs `t` on each frame, in seconds; the plugin sets it
from the SigMF `core:datetime`. Each aircraft's summary covers:
- **Identity:** country (ICAO address block), callsign, operator and IATA flight number (`tables.OPERATORS`),
  category, transponder level, ADS-B version, squawk and emergency state.
- **Position:**
  - The first fix needs an even/odd pair within 10 s (global CPR).
  - Later frames decode locally against the previous fix.
  - Given the receiver position (`ref`, or SigMF `core:geolocation`), a single frame also decodes; this adds
    range, bearing and elevation.
  - A fix that implies more than 1000 kt from the previous one is flagged and not used.
- **Kinematics:** altitude range and rate, climbing/level/descending, ground-speed range, turn over the
  observation, GNSS−baro difference.
- **Intent and integrity:** selected altitude/heading, baro setting, autopilot modes; NACp (as a position
  accuracy), NIC, SIL, GVA, NACv, SDA.
- **Capabilities:** TCAS, 1090ES-IN, UAT-IN, ARV, TS, ACAS type from DF0.
- **Comm-B:** ambiguous BDS50/BDS60 replies are resolved against the aircraft's own ADS-B speed, track and
  altitude (pyModeS `is50or60`).
- **Met:**
  - Wind is the ground vector minus the air vector. The air vector uses BDS50 TAS and BDS60 heading. That heading
    is magnetic, used as true; India's declination is 0–2°.
  - Static air temperature is (TAS / (38.97·Mach))², with its quantisation uncertainty and the ISA deviation.
  - BDS44, when present, is reported as sent.
- **Interrogation environment:** ACAS replies (other aircraft's TCAS), Mode S radar replies, beam dwells and the
  radar revisit period, radar IIS codes, and the DF11 II/SI codes.
- **Signal and flags:** SNR range, frame rate, and flags for changed callsign or squawk, emergency squawks,
  implausible positions, GNSS/baro disagreement beyond what temperature explains, a Comm-B callsign that differs
  from the ADS-B one, and ACAS RAs.

## How
- The 1090 MHz band (±2 MHz) is resampled to 8 MSps, and the magnitude envelope is taken.
- An 8 µs preamble correlation (pulses at 0/1/3.5/4.5 µs) is followed by PPM bit decisions.
- CRC-24 uses the polynomial 0xFFF409.

## API (depends on `spectrum_intel.core` + pyModeS)
```python
from spectrum_intel.receivers.adsb import demod, fields, tracks, FS
msgs = fields(demod(np.abs(y_at_FS), thr))
for m in msgs: m["t"] = t0 + m["t_s"]
aircraft = tracks(msgs, ref=(lat, lon, alt_m))      # ref optional
```

## Validation
- `tests/test_core.py`: CRC on a public DF17 frame (KLM1023); a single bit error is detected.
- `tests/test_adsb.py`: on real frames from a live capture (2026-10-09, `testdata/adsb/`):
  - every format's fields;
  - 2 aircraft with 31 + 3 positions;
  - BDS50/60 resolution, wind and temperature;
  - a local decode against a receiver position;
  - a flagged position jump;
  - the plugin end to end on real IQ.
- The positions are self-consistent: a smooth climbing left turn out of Hyderabad, matching the BDS50 roll
  (−24°) and the TC19 track. They have not been checked against an external tracker.

## Capturing
The sweep's 32.5 ms snapshots catch about 1 frame per aircraft, too few for a position. Record 1090 MHz at
8 MSps for seconds, with a view of the sky (indoors gave nothing), at gain about 60.

## UAT 978 MHz (`uat.py`, plugin `uat`; research, PLAN G14)
- DO-282B ADS-B downlink: CPFSK 1.041667 Mb/s, h 0.6. The 978 MHz channel at 4 samples / bit; the phase change over one
  bit; sync 0xEACDDA4E2 (≤ 4 bit errors) at the best of 4 timing phases, the slicer centred on the sync's 1 / 0 levels
  (absorbs the CFO).
- Reed-Solomon over GF(256), polynomial 0x187, first root 120 (`core/rs.py`, our own Berlekamp-Massey / Chien / Forney;
  dump978's RS library is LGPL and only its parity output was used as an oracle): Long RS(48,34) when the payload type
  is not 0, else Basic RS(30,18), as FlightAware dump978 (BSD-2).
- Fields: payload type, address qualifier, address, position, altitude, NIC, air / ground, velocities; Mode Status:
  emitter category (14 = unmanned aerial vehicle → `uav`), callsign (base 40), emergency, MOPS. One capture-level signal
  "UAT 978 MHz ADS-B[, UAV]" with per-aircraft summary. Ground uplinks (FIS-B) are not decoded.
- Validation: dump978's sample messages (positions, callsign N5130E) and synthetic IQ (≈7 dB with RS corrections).
  There is no UAT traffic in India.

## SSR Mode A/C replies (`modeac.py`; research, signal gaps S11)

**Validation: research.** ICAO Annex 10 Vol IV 3.1.1.6 via `research/signal_gaps/05_aviation_radar_jammers.md` section
3; the Gillham table checked against pyModeS (test oracle only) on all 2048 words; synthetic replies only.
- **Reply:** framing pulses F1 / F2 20.3 µs apart, 13 positions on the 1.45 µs grid (C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2
  B4 D4), SPI 24.65 µs after F1; 0.45 µs pulses ≥ 12 dB over the envelope median; F1 / F2 within 3 dB; every pulse
  between them on the grid (else garbled: overlapping replies or Mode S data); X empty; pulses ≥ 9 dB weaker than the
  framing pulses are ignored.
- **Fields** (`mode_ac`): per distinct word the Mode A `code` (octal ABCD), `altitude_ft_if_mode_c` (when the bits
  are a valid Gillham code), count, SPI, `special_code` (7500 / 7600 / 7700) and `repeated`. The interrogation is
  usually not in view, so every reply is read both ways. Military IFF Mode 3 equals Mode A.
- **Signals:** with Mode S frames, `mode_ac` is added to the Mode S signal; without any, a signal
  "SSR Mode A/C replies[ (squawk 7700 emergency)]" (0.8 when a word repeats, else "-like" 0.6), status `partial`:
  replies carry no parity, so repetition is the only verification.
