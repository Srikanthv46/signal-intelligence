"""
protocols.py -- which RC / telemetry / identification protocol an FSK-family burst most resembles (research catalog).

Research, not identification: rates, deviations, bands and frame periods from published firmware and datasheets
(research/drone_transmitters/03_rc_links.md, 04_telemetry_remoteid.md); GPL sources are used as reference only (no
code). Validated on synthetic bursts only. Every candidate lists what matched and what conflicted.

Entry: (name, symbol rate Hz, deviation Hz or None, bands MHz, packet period ms or None, note). "dev_ratio" entries
(SiK) match on deviation / rate instead of an absolute deviation.
"""
from __future__ import annotations

S24 = ((2400.0, 2484.0),)
SUB = ((420.0, 480.0), (860.0, 930.0))
P = [
    ("FrSky D8 (ACCST, CC2500)", 31.0e3, 31.7e3, S24, 9.0, "47 ch on a 299.927 kHz grid"),
    ("RadioLink (CC2500)", 38.0e3, 38.0e3, S24, 20.0, "16 ch, 400 kHz grid from 2402 MHz"),
    ("FrSky D16 / ACCST v1 (FCC)", 70.0e3, 57.0e3, S24, 9.0, "47 ch, 299.927 kHz grid"),
    ("FrSky ACCST v2 (D16 v2)", 77.0e3, 57.0e3, S24, 9.0, "hop table from the 16-bit TX ID"),
    ("FrSky D16 LBT (EU)", 100.0e3, 70.0e3, S24, 9.0, "2-FSK, base 2405 MHz"),
    ("TBS Crossfire 150 Hz (FSK)", 85.1e3, 42.48e3, ((860.0, 930.0),), 6.666, "260 kHz grid, downlink +50 channels"),
    ("mLRS 50 Hz FSK", 100.0e3, 50.0e3, SUB, 20.0, "91-byte frames, whitening"),
    ("FLARM-compatible (2-FSK, Manchester)", 100.0e3, 50.0e3, ((860.0, 870.0), (902.0, 928.0)), None,
     "India 866.0 MHz; ~6 ms every 0.6-1.4 s; aircraft type 13 = UAV"),
    ("Futaba S-FHSS (CC2500)", 128.1e3, 38.1e3, S24, 6.8, "30 ch, 1.5 MHz hop grid, packets in pairs"),
    ("Graupner HoTT (MSK, CC2500)", 250.0e3, 62.5e3, S24, 10.0, "MSK + FEC, 338 kHz grid"),
    ("Hitec AFHSS / Minima (MSK, CC2500)", 250.0e3, 62.5e3, S24, 22.5, "sync word = TX ID, 300 kHz grid"),
    ("ExpressLRS LR1121 GFSK", 300.0e3, 100.0e3, S24 + ((860.0, 930.0),), 1.0, "1000 Hz / DVDA"),
    ("FlySky AFHDS 2A (A7105)", 500.0e3, None, S24, 3.85, "16 ch from the TX ID, 0.5 MHz grid"),
    ("ExpressLRS / mLRS / TBS Tracer FLRC (SX1280)", 650.0e3, None, S24, None, "ELRS FLRC sync word = UID[2..5]"),
    ("nRF24 / XN297 toy RC (Bayang, V2x2, ...)", 1.0e6, 160.0e3, S24, None, "Bayang 2 ms, V2x2 4 ms"),
    ("OpenLRS-ng 57.6 kb/s", 57.6e3, 28.75e3, SUB, None, "Si443x, rf_magic bind ID"),
    ("OpenLRS-ng 125 kb/s", 125.0e3, 125.0e3, SUB, None, "Si443x"),
    ("XBee-PRO 900HP 200 kb/s", 200.0e3, None, ((902.0, 928.0),), None, "FSK FHSS, 64 ch"),
]
# SiK (Si443x GFSK, deviation = 625 Hz x reg 0x72; deviation ~ rate): rate kb/s -> deviation kHz (radio_443x.c)
SIK = {2: 1.875, 4: 3.75, 8: 8.125, 16: 16.25, 19: 18.75, 24: 23.75, 32: 31.875, 48: 48.125, 64: 63.75,
       96: 96.25, 128: 128.125, 192: 158.75, 250: 158.75}


def _in(f, bands):
    return f is not None and any(lo <= f <= hi for lo, hi in bands)


def candidates(sym_rate, dev_hz, f_mhz, interval_ms=None, top=3):
    """Ranked protocols for a measured FSK burst (median over a signal's bursts). Rate must agree within 5 %;
    deviation within 25 % (when the protocol fixes it); band and packet period add matches or conflicts."""
    if not sym_rate:
        return []
    out = []

    def add(name, m, c, note):
        out.append(dict(protocol=name, score=len(m) - len(c), matched=m, conflicts=c, note=note))
    for name, R, dev, bands, per, note in P:
        if abs(sym_rate - R) > 0.05 * R:
            continue
        m, c = [f"rate {R / 1e3:g} kb/s"], []
        if dev is not None and dev_hz:
            (m if abs(dev_hz - dev) <= 0.25 * dev else c).append(f"deviation {dev / 1e3:g} kHz")
        if f_mhz is not None:
            (m if _in(f_mhz, bands) else c).append("band")
        if per is not None and interval_ms is not None:
            (m if abs(interval_ms - per) <= 0.05 * per else c).append(f"packet period {per:g} ms")
        add(name, m, c, note)
    for kb, dk in SIK.items():                                 # SiK / RFD900 (SiK-compatible) telemetry
        if abs(sym_rate - kb * 1e3) > 0.05 * kb * 1e3:
            continue
        m, c = [f"rate {kb} kb/s"], []
        if dev_hz:
            (m if abs(dev_hz - dk * 1e3) <= 0.25 * dk * 1e3 else c).append(f"deviation {dk:g} kHz (h ~ 2)")
        if f_mhz is not None:
            (m if _in(f_mhz, ((433.0, 435.0), (470.0, 471.0), (868.0, 870.0), (902.0, 928.0))) else c).append("band")
        add(f"SiK telemetry radio {kb} kb/s (3DR / Holybro / mRo / RFD900 compatible)", m, c,
            "GFSK, PN9 whitening, NETID header; TDM window ~114 ms at 64 kb/s")
    out.sort(key=lambda d: (-d["score"], len(d["conflicts"])))
    return out[:top]
