"""Reference 802.11b DSSS/CCK transmitter (IEEE 802.11-2020 Clause 15/16), for validating the receiver.

    ppdu_chips(mpdu, rate_mbps, preamble)  -> complex chips at 11 Mchip/s (FCS appended)
    beacon(bssid, ssid, channel, ...)      -> MPDU bytes without FCS

Bits go out LSB first per field; the scrambler is x^7+x^4+1 (seed 0x6C long, 0x1B short).
"""
from __future__ import annotations

import numpy as np

BARKER = np.array([1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1], float)
SIGNAL = {1.0: 0x0A, 2.0: 0x14, 5.5: 0x37, 11.0: 0x6E}
DQPSK = {(0, 0): 0.0, (0, 1): np.pi / 2, (1, 1): np.pi, (1, 0): 3 * np.pi / 2}
QPSK = {(0, 0): 0.0, (0, 1): np.pi / 2, (1, 0): np.pi, (1, 1): 3 * np.pi / 2}


def _bits(value, n):
    return [(value >> i) & 1 for i in range(n)]


def _byte_bits(data):
    return [b for byte in data for b in _bits(byte, 8)]


def _scramble(bits, seed):
    reg = [(seed >> i) & 1 for i in range(7)]          # reg[k] = s(n-1-k)
    out = []
    for b in bits:
        s = b ^ reg[3] ^ reg[6]
        out.append(s)
        reg = [s] + reg[:6]
    return out


def crc16(bits):
    crc = 0xFFFF
    for b in bits:
        fb = ((crc >> 15) & 1) ^ b
        crc = (crc << 1) & 0xFFFF
        if fb:
            crc ^= 0x1021
    crc ^= 0xFFFF
    return [(crc >> (15 - i)) & 1 for i in range(16)]


def crc32(data):
    crc = 0xFFFFFFFF
    for byte in data:
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (0xEDB88320 & -(crc & 1))
    return crc ^ 0xFFFFFFFF


def plcp_header(n_bytes, rate):
    """SIGNAL, SERVICE, LENGTH (us) and CRC-16 bits for an n_bytes PSDU."""
    service = 0x04                                      # locked clocks
    if rate in (5.5, 11.0):
        length = int(np.ceil(8 * n_bytes / rate))
        if rate == 11.0 and length * 11 - 8 * n_bytes >= 8:
            service |= 0x80                             # length extension
    else:
        length = int(8 * n_bytes / rate)
    hdr = _bits(SIGNAL[rate], 8) + _bits(service, 8) + _bits(length, 16)
    return hdr + crc16(hdr)


def ppdu_chips(mpdu, rate=1.0, preamble="long"):
    """One PPDU as 11 Mchip/s complex chips. `mpdu` excludes the FCS, which is appended."""
    psdu = bytes(mpdu) + crc32(bytes(mpdu)).to_bytes(4, "little")
    if preamble == "long":
        pre = [1] * 128 + _bits(0xF3A0, 16)
        seed, hdr_rate = 0x6C, 1.0
    else:
        pre = [0] * 56 + _bits(0x05CF, 16)
        seed, hdr_rate = 0x1B, 2.0
    hdr = plcp_header(len(psdu), rate)
    data = _byte_bits(psdu)
    s = _scramble(pre + hdr + data, seed)
    n_pre = len(pre)
    syms, ph = [], 0.0
    for b in s[:n_pre]:                                 # SYNC + SFD: DBPSK
        ph += np.pi * b
        syms.append(np.exp(1j * ph) * BARKER)
    hb = s[n_pre:n_pre + 48]
    if hdr_rate == 1.0:
        for b in hb:
            ph += np.pi * b
            syms.append(np.exp(1j * ph) * BARKER)
    else:
        for i in range(0, 48, 2):
            ph += DQPSK[(hb[i], hb[i + 1])]
            syms.append(np.exp(1j * ph) * BARKER)
    db = s[n_pre + 48:]
    if rate == 1.0:
        for b in db:
            ph += np.pi * b
            syms.append(np.exp(1j * ph) * BARKER)
    elif rate == 2.0:
        for i in range(0, len(db), 2):
            ph += DQPSK[(db[i], db[i + 1])]
            syms.append(np.exp(1j * ph) * BARKER)
    else:
        bps = 8 if rate == 11.0 else 4
        for k, i in enumerate(range(0, len(db), bps)):
            d = db[i:i + bps]
            ph += DQPSK[(d[0], d[1])] + (np.pi if k % 2 else 0.0)
            if rate == 11.0:
                p2, p3, p4 = QPSK[(d[2], d[3])], QPSK[(d[4], d[5])], QPSK[(d[6], d[7])]
            else:
                p2, p3, p4 = d[2] * np.pi + np.pi / 2, 0.0, d[3] * np.pi
            syms.append(np.exp(1j * ph) * np.array([
                np.exp(1j * (p2 + p3 + p4)), np.exp(1j * (p3 + p4)), np.exp(1j * (p2 + p4)), -np.exp(1j * p4),
                np.exp(1j * (p2 + p3)), np.exp(1j * p3), -np.exp(1j * p2), 1.0]))
    return np.concatenate(syms)


def beacon(bssid="02:11:22:33:44:55", ssid="test", channel=6, seq=1, tsf_us=123456789, interval_tu=100,
           privacy=True):
    """A minimal beacon: SSID, supported rates, DS parameter set, and RSN (WPA2-PSK, CCMP) when privacy."""
    mac = bytes(int(x, 16) for x in bssid.split(":"))
    hdr = bytes([0x80, 0x00, 0, 0]) + b"\xff" * 6 + mac + mac + (seq << 4).to_bytes(2, "little")
    cap = 0x0001 | 0x0020 | (0x0010 if privacy else 0)  # ESS, short preamble, privacy
    body = tsf_us.to_bytes(8, "little") + interval_tu.to_bytes(2, "little") + cap.to_bytes(2, "little")
    s = ssid.encode()
    body += bytes([0, len(s)]) + s + bytes([1, 4, 0x82, 0x84, 0x8B, 0x96]) + bytes([3, 1, channel])
    if privacy:
        rsn = (bytes([1, 0]) + bytes([0x00, 0x0F, 0xAC, 4]) + bytes([1, 0, 0x00, 0x0F, 0xAC, 4])
               + bytes([1, 0, 0x00, 0x0F, 0xAC, 2]) + bytes([0, 0]))
        body += bytes([48, len(rsn)]) + rsn
    return hdr + body


def data_frame(bssid="02:11:22:33:44:55", sta="02:aa:bb:cc:dd:ee", n_payload=100, seq=7, seed=0):
    """A protected QoS-data frame from the AP (FromDS): CCMP header (PN = seq), random ciphertext."""
    a = lambda m: bytes(int(x, 16) for x in m.split(":"))  # noqa: E731
    hdr = bytes([0x88, 0x42, 0x2C, 0x00]) + a(sta) + a(bssid) + a(bssid) + (seq << 4).to_bytes(2, "little")
    hdr += bytes([0, 0])                                # QoS control
    pn = (seq + 1).to_bytes(6, "little")
    ccmp = bytes([pn[0], pn[1], 0, 0x20]) + pn[2:6]     # Ext IV, key ID 0
    rng = np.random.default_rng(seed)
    return hdr + ccmp + bytes(rng.integers(0, 256, n_payload, dtype=np.uint8))


# ------------------------------------------------------------------------------------------- channel model
CHIP_RATE = 11e6


def _rc_pulse(t, beta=0.5):
    t = np.asarray(t, float)
    den = 1 - (2 * beta * t) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.sinc(t) * np.cos(np.pi * beta * t) / den
    h[np.abs(den) < 1e-8] = np.pi / 4 * np.sinc(1 / (2 * beta))
    return h


def render(chips, n_out, fs, t0_chips=0.0, ppm=0.0, taps=(), span=6):
    """Chips as a continuous raised-cosine waveform sampled at fs, with the transmitter's chip clock off by
    ppm and multipath taps [(delay_ns, complex gain)] added to the direct path."""
    n = np.arange(n_out)
    out = np.zeros(n_out, complex)
    for delay_ns, g in [(0.0, 1.0)] + list(taps):
        t = n / fs * CHIP_RATE * (1 + ppm * 1e-6) - t0_chips - delay_ns * 1e-9 * CHIP_RATE
        k0 = np.floor(t).astype(int)
        for d in range(-span, span + 1):
            k = k0 + d
            ok = (k >= 0) & (k < len(chips))
            out[ok] += g * chips[k[ok]] * _rc_pulse(t[ok] - k[ok])
    return out


def scene(ppdus, n, fs=61.44e6, fc=2437e6, noise=1.0, seed=0):
    """A capture holding PPDUs. Each is a dict: mpdu, rate, preamble, f_hz (carrier), t_s, snr_db (in the
    22 MHz channel), and optional cfo_hz, ppm (chip clock; default locked to cfo), taps."""
    rng = np.random.default_rng(seed)
    x = (rng.normal(size=n) + 1j * rng.normal(size=n)) * np.sqrt(noise / 2)
    for p in ppdus:
        chips = ppdu_chips(p["mpdu"], p.get("rate", 1.0), p.get("preamble", "long"))
        i0 = int(p["t_s"] * fs)
        m = min(int(len(chips) / CHIP_RATE * fs) + 200, n - i0)
        cfo = p.get("cfo_hz", 0.0)
        ppm = p.get("ppm", cfo / p["f_hz"] * 1e6)
        y = render(chips, m, fs, rng.uniform(0, 1), ppm, p.get("taps", ()))
        a = np.sqrt(10 ** (p["snr_db"] / 10) * noise * 22e6 / fs / np.mean(np.abs(y[:min(m, 5 * len(chips))]) ** 2))
        f = p["f_hz"] - fc + cfo
        x[i0:i0 + m] += a * y * np.exp(2j * np.pi * f * (np.arange(m) + i0) / fs + 1j * rng.uniform(0, 2 * np.pi))
    return x.astype(np.complex64)
