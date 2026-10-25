"""4-level FSK land-mobile radio: DMR (ETSI TS 102 361-1) and NXDN, from the FM discriminator output.

    discriminator(y, fs)        -> instantaneous frequency (Hz), median removed
    symbol_line(d, fs)          -> (rate, ratio): 2400 / 4800 symbol-rate line of the discriminator output
    dmr(d, fs)                  -> DMR bursts: sync type, colour code / data type (Golay), BPTC check, CSBK CRC
    nxdn(d, fs, rate)           -> NXDN frame-sync words

Dibits (TS 102 361-1 10.1.2; NXDN uses the same map): +3 -> 01, +1 -> 00, -1 -> 10, -3 -> 11.
DMR bursts are demodulated one at a time: RRC (0.2) matched filter, timing and level/offset taken from the
burst's own 24-symbol sync (outer symbols only), then sliced. Depends on numpy only.
"""
from __future__ import annotations

import numpy as np

RATES = (2400.0, 4800.0)
SYNC_ERR_COARSE = 6
SYNC_ERR = 4
SYNC_CORR_MIN = 0.9                      # matched-filtered sync vs its ideal symbols (real bursts: >= 0.98)
DMR_DEV_UNIT = 648.0                     # Hz per symbol level (+-1 -> +-648 Hz, +-3 -> +-1944 Hz), 10.1.2
DMR_GAIN_RANGE = (0.5, 1.5)              # fitted level vs DMR_DEV_UNIT
# TS 102 361-1 9.1.1 SYNC patterns (48 bits)
DMR_SYNC = {"bs_voice": 0x755FD7DF75F7, "bs_data": 0xDFF57D75DF5D, "ms_voice": 0x7F7D5DD57DFD,
            "ms_data": 0xD5D7F77FD757, "ms_rc": 0x77D55F7DFD77, "direct1_voice": 0x5D577F7757FF,
            "direct1_data": 0xF7FDD5DDFD55, "direct2_voice": 0x7DFFD5F55D5F, "direct2_data": 0xD7557F5FF7F5}
NXDN_FSW = 0xCDF59                       # 20-bit frame sync word (10 symbols)
NXDN_FRAME_BITS = 384                    # 192 symbols
DATA_TYPES = {0: "PI header", 1: "voice LC header", 2: "terminator with LC", 3: "CSBK", 4: "MBC header",
              5: "MBC continuation", 6: "data header", 7: "rate 1/2 data", 8: "rate 3/4 data", 9: "idle",
              10: "rate 1 data", 11: "unified single block data"}
CRC_MASK = {0: 0x6969, 3: 0xA5A5, 4: 0xAAAA, 6: 0xCCCC, 11: 0x3333}      # B.3.12 (CRC-CCITT data types)
BPTC_TYPES = {0, 1, 2, 3, 4, 5, 6, 7, 9, 11}                              # BPTC(196,96) payloads (not 8, 10)


def _bits(v: int, n: int) -> np.ndarray:
    return np.array([(v >> (n - 1 - i)) & 1 for i in range(n)], np.int8)


_SYM = {(0, 1): 3.0, (0, 0): 1.0, (1, 0): -1.0, (1, 1): -3.0}
_DIBIT = {3: (0, 1), 1: (0, 0), -1: (1, 0), -3: (1, 1)}


def _symbols(v: int, n_bits: int) -> np.ndarray:
    b = _bits(v, n_bits)
    return np.array([_SYM[(b[2 * i], b[2 * i + 1])] for i in range(n_bits // 2)], np.float32)


def discriminator(y: np.ndarray, fs: float) -> np.ndarray:
    d = np.angle(y[1:] * np.conj(y[:-1])) * (fs / (2 * np.pi))
    return d - np.median(d)


def symbol_line(d: np.ndarray, fs: float) -> tuple[float | None, float]:
    """Strongest 4FSK symbol-rate line in the squared derivative of the discriminator output (4 kHz low-passed by
    a moving average), over the median of its +-40..400 Hz neighbourhood."""
    k = max(1, int(round(fs / 9600)))
    dl = np.convolve(d, np.ones(k) / k, "same")
    e = np.diff(dl) ** 2
    if len(e) < 64:
        return None, 0.0
    E = np.abs(np.fft.rfft(e - e.mean()))
    fe = np.fft.rfftfreq(len(e), 1 / fs)
    best = (None, 0.0)
    for r in RATES:
        on, bg = np.abs(fe - r) <= 15, (np.abs(fe - r) > 40) & (np.abs(fe - r) < 400)
        if on.any() and bg.any():
            q = float(E[on].max() / (np.median(E[bg]) + 1e-30))
            if q > best[1]:
                best = (r, q)
    return best


def _rrc(beta: float, sps: float, span: int = 6) -> np.ndarray:
    t = np.arange(-int(span * sps), int(span * sps) + 1) / sps
    h = np.sinc(t) * np.cos(np.pi * beta * t) / (1 - (2 * beta * t) ** 2 + 1e-9)
    return h / h.sum()


def _coarse_bits(d: np.ndarray, sps: float) -> tuple[np.ndarray, float]:
    """Boxcar-integrated, one timing phase for the whole block, levels from the 75th percentile -> bits.
    Only used to locate sync words; returns (bits, phase in samples)."""
    w = max(1, int(round(0.6 * sps)))
    z = np.convolve(d, np.ones(w) / w, "same")
    best = None
    for ph in np.arange(0, sps, sps / 10):
        idx = (ph + sps * np.arange(int((len(z) - ph) / sps) - 1)).astype(int)
        if len(idx) < 24:
            continue
        v = z[idx]
        a = np.percentile(np.abs(v), 75) + 1e-9
        q = np.clip(np.round((v / a * 3 + 3) / 2) * 2 - 3, -3, 3)
        err = float(np.mean((v / a * 3 - q) ** 2))
        if best is None or err < best[0]:
            best = (err, v, a, ph)
    if best is None:
        return np.zeros(0, np.int8), 0.0
    _, v, a, ph = best
    q = np.clip(np.round((v / a * 3 + 3) / 2) * 2 - 3, -3, 3).astype(int)
    return np.array([b for s in q for b in _DIBIT[int(s)]], np.int8), float(ph)


def _find(bits: np.ndarray, pattern: int, n: int, max_err: int, with_err: bool = False) -> list:
    if len(bits) < n:
        return []
    c = np.correlate(2 * bits.astype(np.float32) - 1, 2 * _bits(pattern, n).astype(np.float32) - 1, "valid")
    e = np.rint((n - c) / 2).astype(int)
    hits = np.flatnonzero(e <= max_err)
    return [(int(i), int(e[i])) for i in hits] if with_err else [int(i) for i in hits]


# ------------------------------------------------------------------------------------------- DMR FEC
_GOLAY_G = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0],
                     [0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
                     [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1],
                     [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1],
                     [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],
                     [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1],
                     [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0],
                     [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1]], np.int8)  # B.3.1 table B.11
_GOLAY_CW = (np.array([[(v >> (7 - i)) & 1 for i in range(8)] for v in range(256)], np.int8) @ _GOLAY_G) % 2
_H15_P = np.array([[1, 0, 0, 1], [1, 1, 0, 1], [1, 1, 1, 1], [1, 1, 1, 0], [0, 1, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1],
                   [1, 0, 1, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]], np.int8)              # B.3.4 B.15
_H13_P = np.array([[1, 1, 1, 1], [1, 1, 1, 0], [0, 1, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 1], [1, 1, 0, 0],
                   [0, 1, 1, 0], [0, 0, 1, 1]], np.int8)                                         # B.3.4 B.14
_DEINT = np.array([(181 * i) % 196 for i in range(196)])                 # B.1.1: matrix bit i sent at 181 i mod 196


def golay_20_8(r: np.ndarray) -> tuple[int | None, int]:
    """Nearest Golay(20,8) codeword (minimum distance 8: up to 3 errors corrected). Returns (data byte, errors);
    None when the distance is > 3."""
    dist = np.sum(_GOLAY_CW != r[None, :], axis=1)
    v = int(np.argmin(dist))
    return (v if dist[v] <= 3 else None), int(dist[v])


def _hamming_fix(word: np.ndarray, P: np.ndarray) -> bool:
    """Single-error correction in place of a systematic Hamming codeword [data | data P]. True if valid after."""
    k = P.shape[0]
    syn = (word[:k] @ P + word[k:]) % 2
    if not syn.any():
        return True
    cols = np.vstack([P, np.eye(P.shape[1], dtype=np.int8)])            # syndrome of an error at each position
    hit = np.flatnonzero((cols == syn).all(axis=1))
    if len(hit) == 1:
        word[hit[0]] ^= 1
        return True
    return False


def bptc_196_96(r: np.ndarray) -> tuple[np.ndarray, int, int]:
    """De-interleave and decode BPTC(196,96) (B.1.1): Hamming(13,9) columns then (15,11) rows, single-error
    correction each. Returns (96 info bits I(95)..I(0), rows valid of 9, columns valid of 15)."""
    m = r[_DEINT].astype(np.int8)
    M = m[1:].reshape(13, 15).copy()
    cols = sum(_hamming_fix(M[:, j], _H13_P) for j in range(15))
    rows = sum(_hamming_fix(M[i], _H15_P) for i in range(9))
    return np.concatenate([M[0, 3:11], M[1:9, :11].ravel()]), int(rows), int(cols)


def crc_ccitt(bits: np.ndarray) -> int:
    """B.3.8: CRC-CCITT, initial remainder 0, inverted."""
    r = 0
    for b in bits:
        fb = ((r >> 15) & 1) ^ int(b)
        r = (r << 1) & 0xFFFF
        if fb:
            r ^= 0x1021
    return r ^ 0xFFFF


def _octets(b: np.ndarray) -> bytes:
    return bytes(int("".join(map(str, b[8 * i:8 * i + 8])), 2) for i in range(len(b) // 8))


# ------------------------------------------------------------------------------------------- DMR
def dmr(d: np.ndarray, fs: float, sr: float = 4800.0) -> list[dict]:
    """Every DMR burst whose 48-bit SYNC is in the block (one record per burst)."""
    sps = fs / sr
    bits, _ = _coarse_bits(d, sps)
    z = np.convolve(d, _rrc(0.2, sps), "same")
    out, seen = [], set()
    for name, pat in DMR_SYNC.items():
        ref = _symbols(pat, 48)
        for h in _find(bits, pat, 48, SYNC_ERR_COARSE):
            s0 = (h // 2) * sps
            best = None
            for dt in np.arange(-1.5 * sps, 1.5 * sps, sps / 20):          # timing from the sync itself
                idx = (s0 + dt + sps * np.arange(24)).astype(int)
                if idx.min() < 0 or idx.max() >= len(z):
                    continue
                v = z[idx] - z[idx].mean()
                c = float(np.dot(v, ref) / (np.linalg.norm(v) * np.linalg.norm(ref) + 1e-9))
                if best is None or c > best[0]:
                    best = (c, dt)
            if best is None:
                continue
            c, dt = best
            if c < SYNC_CORR_MIN:
                continue
            t_sync = s0 + dt
            key = int(round(t_sync / sps))
            if any(abs(key - k) < 12 for k in seen):
                continue
            idx = (t_sync + sps * np.arange(24)).astype(int)
            gain, dc = np.polyfit(ref.astype(np.float64), z[idx], 1)
            # the level fit would scale anything to the pattern: DMR's own deviation must come out of it
            if not DMR_GAIN_RANGE[0] * DMR_DEV_UNIT <= gain <= DMR_GAIN_RANGE[1] * DMR_DEV_UNIT:
                continue
            rec = {"sync": name, "t_s": round(float(t_sync / fs), 6), "sync_corr": round(c, 3),
                   "deviation_hz_per_level": round(float(gain), 1)}
            # whole burst: 54 symbols, SYNC (24), 54 symbols
            idx = (t_sync - 54 * sps + sps * np.arange(132)).astype(int)
            if idx.min() >= 0 and idx.max() < len(z):
                sym = (z[idx] - dc) / gain
                q = np.clip(np.round((sym + 3) / 2) * 2 - 3, -3, 3).astype(int)
                bb = np.array([b for s in q for b in _DIBIT[int(s)]], np.int8)
                rec["sync_bit_errors"] = int(np.sum(bb[108:156] != _bits(pat, 48)))
                if rec["sync_bit_errors"] > SYNC_ERR:
                    continue
                rec["symbol_mse"] = round(float(np.mean((sym - q) ** 2)), 3)
                if name.endswith("data") or name == "ms_rc":
                    _data_burst(bb, rec)
            elif name.endswith("data"):                      # slot type alone (5 ms either side of the sync)
                idx = (t_sync - 5 * sps + sps * np.arange(34)).astype(int)
                if idx.min() >= 0 and idx.max() < len(z):
                    sym = (z[idx] - dc) / gain
                    q = np.clip(np.round((sym + 3) / 2) * 2 - 3, -3, 3).astype(int)
                    bb = np.array([b for s in q for b in _DIBIT[int(s)]], np.int8)
                    _slot_type(np.concatenate([bb[:10], bb[58:68]]), rec)
            seen.add(key)
            out.append(rec)
    return sorted(out, key=lambda r: r["t_s"])


def _slot_type(st: np.ndarray, rec: dict) -> None:
    v, nerr = golay_20_8(st)
    if v is not None:
        rec.update(colour_code=v >> 4, data_type=DATA_TYPES.get(v & 15, f"reserved {v & 15}"),
                   data_type_code=v & 15, slot_type_errors=nerr)


def _data_burst(bb: np.ndarray, rec: dict) -> None:
    _slot_type(np.concatenate([bb[98:108], bb[156:166]]), rec)
    dt = rec.get("data_type_code")
    if dt not in BPTC_TYPES:
        return
    info, rows, cols = bptc_196_96(np.concatenate([bb[:98], bb[166:264]]))
    rec.update(bptc_rows_valid=rows, bptc_cols_valid=cols)
    if dt in CRC_MASK and rows == 9 and cols == 15:
        crc_rx = int("".join(map(str, info[80:96])), 2) ^ CRC_MASK[dt]
        rec["crc_ok"] = crc_ccitt(info[:80]) == crc_rx
        if rec["crc_ok"] and dt == 3:
            o = _octets(info[:80])
            rec.update(csbk_last_block=o[0] >> 7, csbk_protect=(o[0] >> 6) & 1, csbko=o[0] & 0x3F,
                       feature_set_id=o[1], csbk_data_hex=o[2:10].hex())


# ------------------------------------------------------------------------------------------- NXDN
def nxdn(d: np.ndarray, fs: float, sr: float) -> list[tuple[int, int]]:
    """(sample position, bit errors) of NXDN frame-sync words (<= 1 bit error of 20)."""
    sps = fs / sr
    bits, ph = _coarse_bits(d, sps)
    return [(int(ph + (h // 2) * sps), e) for h, e in _find(bits, NXDN_FSW, 20, 1, with_err=True)]


def frame_consistent(pos_samples: list[int], fs: float, sr: float) -> int:
    """Largest number of frame-sync hits on one 192-symbol frame grid (within 2 symbols)."""
    if not pos_samples:
        return 0
    fr = 192 * fs / sr
    tol = 2 * fs / sr
    best = 0
    for p in pos_samples:
        best = max(best, sum(abs(((q - p) / fr) - round((q - p) / fr)) * fr <= tol for q in pos_samples))
    return best
