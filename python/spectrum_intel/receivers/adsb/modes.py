"""Mode S / ADS-B 1090 MHz: envelope, preamble correlation, PPM bits, CRC-24."""
from __future__ import annotations

import numpy as np

from .messages import fields  # noqa: F401  (re-export: fields used to live here)


FS = 8e6
SPU = 8                                   # samples per us
F_ADSB = 1090e6


def _crc24(bits):
    g = 0xFFF409
    reg = 0
    for b in bits:
        top = (reg >> 23) & 1
        reg = ((reg << 1) & 0xFFFFFF) | int(b)
        if top:
            reg ^= g
    for _ in range(24):
        top = (reg >> 23) & 1
        reg = (reg << 1) & 0xFFFFFF
        if top:
            reg ^= g
    return reg


def _residual(bits):
    """Mode S parity: CRC of the data bits xor the transmitted parity field."""
    n = len(bits)
    data, par = bits[:n - 24], bits[n - 24:]
    p = 0
    for b in par:
        p = (p << 1) | int(b)
    return _crc24(data) ^ p


def demod(env, thr):
    pre = np.zeros(8 * SPU)
    for t in (0, 1, 3.5, 4.5):
        pre[int(t * SPU):int(t * SPU) + SPU // 2] = 1
    pre = pre - pre.mean()
    c = np.correlate(env, pre, "valid")
    cand = np.where((c > thr) & (c == np.maximum.reduce([np.roll(c, k) for k in range(-3, 4)])))[0]
    out, busy = [], -1
    noise = float(np.median(env)) or 1e-12
    for s in cand:
        if s < busy:
            continue
        d0 = s + 8 * SPU
        if d0 + 112 * SPU > len(env):
            break
        a = env[d0:d0 + 112 * SPU].reshape(112, SPU)
        bits = (a[:, :SPU // 2].sum(1) > a[:, SPU // 2:].sum(1)).astype(np.uint8)
        df = int("".join(map(str, bits[:5])), 2)
        n = 112 if df >= 16 else 56
        bits = bits[:n]
        res = _residual(bits)
        # DF17/18: CRC exact. DF11: parity = interrogator id (< 128). DF0/4/5/16/20/21:
        # parity = aircraft address -> only kept if that address is confirmed elsewhere
        if df in (17, 18):
            ok = res == 0
        elif df == 11:
            ok = res < 128
        else:
            ok = df in (0, 4, 5, 16, 20, 21)
        if not ok:
            continue
        hexmsg = "%0*X" % (n // 4, int("".join(map(str, bits)), 2))
        on = np.maximum(a[:n, :SPU // 2].mean(1), a[:n, SPU // 2:].mean(1)).mean()    # pulse amplitude
        out.append({"t_s": s / FS, "df": df, "bits": n, "hex": hexmsg, "parity_residual": f"{res:06X}",
                    "snr_db": round(float(20 * np.log10(on / noise)), 1),
                    "aa": "%06X" % int("".join(map(str, bits[8:32])), 2) if df in (11, 17, 18) else None})
        busy = d0 + n * SPU
    # address/parity frames are accepted only for addresses confirmed by DF11/17/18 (dump1090 rule)
    known = {m["aa"] for m in out if m["aa"]}
    for m in out:
        if m["df"] not in (11, 17, 18):
            m["confirmed"] = m["parity_residual"] in known
            m["aa"] = m["parity_residual"]
        else:
            m["confirmed"] = True
    return [m for m in out if m["confirmed"]]
