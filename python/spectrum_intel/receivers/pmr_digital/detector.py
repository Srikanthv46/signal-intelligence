"""Digital land-mobile radio detector: TETRA, DMR, NXDN.

    detect(ch, ctx=None) -> list[Emission]      families "tetra", "dmr", "nxdn"

1. Candidates: the capture spectrum (Channelizer FFT) in 1 kHz bins; peaks >= 8 dB over the median floor, at most
   30 kHz wide, inside the land-mobile ranges.
2. TETRA (380-400, 410-430, 450-470, 870-876, 915-921 MHz): the 25 kHz channel at 144 kSps (8 samples/symbol);
   the 18 ksym/s line of the matched-filtered envelope must stand out; pi/4-DQPSK demodulation; synchronisation
   and normal training sequences; BSCH of every synchronisation burst (CRC) -> MCC, MNC, colour code.
3. 4FSK (136-174, 300-520, 806-870 MHz): the 12.5 kHz channel at 48 kSps (constant envelope), FM discriminator;
   DMR SYNC words (4800 sym/s) -> slot type (Golay), BPTC(196,96), CSBK CRC; else NXDN frame sync (4800 or
   2400 sym/s). The 2400/4800 symbol-rate line is reported, not required.
A carrier is reported only with a verified structure (training/sync sequences); a TETRA symbol-rate line alone
gives a "-like" TETRA carrier. Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np
from scipy import signal as ss

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.core.emission import Emission

from . import fsk4, tetra

TETRA_RANGES = [(380e6, 400e6), (410e6, 430e6), (450e6, 470e6), (870e6, 876e6), (915e6, 921e6)]
FSK_RANGES = [(136e6, 174e6), (300e6, 520e6), (806e6, 870e6)]
RES = 1000.0
PEAK_DB = 8.0
MAX_WIDTH = 30e3
TETRA_LINE = 6.0                 # 18 ksym/s envelope line over its neighbourhood: worth demodulating
TETRA_LIKE_LINE = 9.0            # ... and enough on its own (no training sequence found) for a "-like" carrier
MAX_ENV_CV = 0.35                # 4FSK is constant-envelope
FSK_FS = 48000.0


def _in(f, ranges):
    return any(a <= f <= b for a, b in ranges)


def candidates(ch: Channelizer, band: tuple[float, float], max_n: int = 24) -> list[dict]:
    fs, fc, L = ch.fs, ch.fc, ch.L
    lo, hi = band
    rng = [r for r in TETRA_RANGES + FSK_RANGES if r[0] < hi and r[1] > lo]
    if not rng:
        return []
    g = max(1, int(round(RES * L / fs)))
    order = np.argsort(ch.f)
    P = np.abs(ch.X[order]) ** 2
    f = ch.f[order]
    n = len(P) // g
    Pg, fg = P[:n * g].reshape(n, g).mean(1), f[:n * g].reshape(n, g).mean(1)
    fa = fc + fg
    inb = (fa >= lo) & (fa <= hi)
    if inb.sum() < 20:
        return []
    floor = float(np.median(Pg[inb]))
    L_db = 10 * np.log10(Pg / floor + 1e-30)
    pk, pr = ss.find_peaks(np.where(inb, L_db, 0.0), height=PEAK_DB, prominence=PEAK_DB,
                           width=(None, MAX_WIDTH / RES), rel_height=0.5, distance=max(1, int(6e3 / RES)))
    out = []
    for i, w in zip(pk, pr["widths"]):
        f_abs = float(fa[i])
        if not any(a <= f_abs <= b for a, b in rng):
            continue
        a, b = max(0, i - 12), min(n, i + 13)                # centroid over +-12 kHz (a TETRA channel is flat)
        p = np.maximum(Pg[a:b] - floor, 0.0)
        fcen = float(np.sum(fa[a:b] * p) / np.sum(p)) if p.sum() > 0 else f_abs
        out.append({"f_hz": fcen, "peak_db": float(L_db[i]), "width_hz": float(w * RES), "floor": floor / (fs * L)})
    out.sort(key=lambda c: -c["peak_db"])
    keep = []
    for c in out:                                            # one candidate per 12.5 kHz
        if all(abs(c["f_hz"] - k["f_hz"]) >= 12.5e3 for k in keep):
            keep.append(c)
    return keep[:max_n]


def _snr(y: np.ndarray, floor_density: float, bw: float) -> float | None:
    noise = floor_density * bw
    p = float(np.mean(np.abs(y) ** 2))
    return round(10 * np.log10(max(p - noise, 1e-30) / noise), 1) if noise > 0 and p > noise else None


def _tetra(ch, c, dur):
    y, fs = ch.slice(c["f_hz"], 12.5e3, tetra.SR * tetra.SPS)
    line = tetra.symbol_line(y, fs)
    if line < TETRA_LINE:
        return None
    best = None
    for q in (0, -1, 1):                                     # quarter-turn ambiguity of the carrier estimate
        bits, soft, err, foff = tetra.demod(y, q)
        tr = tetra.find_training(bits)
        sync, normal = [h for h, _ in tr["sync"]], [h for h, _ in tr["normal"]]
        bsch = []
        for h in sync:
            if h >= tetra.BSCH_BITS:
                ok, pdu = tetra.decode_bsch(soft[h - tetra.BSCH_BITS:h])
                if ok:
                    bsch.append({**pdu, "t_s": round(h / 2 / tetra.SR, 6)})
        n_norm = tetra.slot_consistent(normal)
        # the 22-bit normal training sequences also hit by chance in a long capture: they count only when most
        # of them share one slot phase; the 38-bit sync sequence practically never does
        norm_ok = n_norm if n_norm >= max(3, 0.5 * len(normal)) else 0
        score = (len(bsch), tetra.slot_consistent(sync), norm_ok)
        if best is None or score > best[0]:
            best = (score, err, foff, sync, normal, bsch, norm_ok)
    (_, n_sync_slot, norm_ok), err, foff, sync, normal, bsch, norm_ok = best
    n_slot = max(n_sync_slot, norm_ok)
    framed = bool(bsch) or n_sync_slot >= 2 or norm_ok >= 3
    f_car = c["f_hz"] + foff
    feats = {"carrier_hz": round(f_car), "symbol_rate_line": round(line, 1), "phase_error_rad": round(err, 3),
             "n_sync_bursts": len(sync),
             "n_normal_training": len(normal), "slot_consistent_hits": n_slot, "n_bsch_crc_ok": len(bsch),
             "modulation": "pi/4-DQPSK 18 ksym/s", "snr_db": _snr(y, c["floor"], 25e3)}
    if bsch:
        ids: dict = {}                                       # identity -> count, in order of first appearance
        for b in bsch:
            k = (b["mcc"], b["mnc"], b["colour_code"])
            ids[k] = ids.get(k, 0) + 1
        top = max(ids, key=ids.get)                          # most frequent; ties: the first seen
        first = next(b for b in bsch if (b["mcc"], b["mnc"], b["colour_code"]) == top)
        feats.update({k: first[k] for k in ("mcc", "mnc", "colour_code", "system_code", "sharing_mode",
                                            "ts_reserved_frames", "u_plane_dtx", "frame18_extension",
                                            "cell_service_level", "late_entry", "neighbour_cell_broadcast")})
        feats["bsch"] = [{k: b[k] for k in ("t_s", "timeslot", "frame", "multiframe")} for b in bsch[:20]]
        feats["bsch_identities"] = len(ids)
        label, conf = "tetra_downlink_bsch", 0.97
    elif framed:
        label, conf = "tetra_downlink", 0.85
    elif line >= TETRA_LIKE_LINE:
        label, conf = "tetra_like_carrier", 0.55
    else:
        return None
    return Emission(family="tetra", label=label, start_us=0.0, end_us=dur * 1e6, f_lo_mhz=(f_car - 12.5e3) / 1e6,
                    f_hi_mhz=(f_car + 12.5e3) / 1e6, center_mhz=f_car / 1e6, snr_db=feats["snr_db"],
                    confidence=conf, decoded={"crc_ok": True} if bsch else None, features=feats)


def _fsk(ch, c, dur):
    y, fs = ch.slice(c["f_hz"], 6.25e3, FSK_FS)
    env = np.abs(y)
    if float(np.std(env) / (np.mean(env) + 1e-30)) > MAX_ENV_CV:
        return None                                           # not constant envelope: not FSK
    d = fsk4.discriminator(y, fs)
    # the symbol-rate line is reported, not required: in a 34 ms block it is too weak to gate on (a DMR burst at
    # 11 dB gives x2.6), while the sync words themselves are the verification
    rate, line = fsk4.symbol_line(d, fs)
    snr = _snr(y, c["floor"], 12.5e3)
    base = {"symbol_rate_line": round(line, 1), "snr_db": snr,
            "deviation_hz_p99": round(float(np.percentile(np.abs(d), 99)))}
    bursts = fsk4.dmr(d, fs, 4800.0)
    if bursts:
        return [_dmr_emission(c, dur, bursts, {**base, "symbol_rate": 4800.0})]
    for sr in (4800.0, 2400.0):
        hits = fsk4.nxdn(d, fs, sr)
        pos = [p for p, _ in hits]
        n_frame = fsk4.frame_consistent(pos, fs, sr)
        # a 34 ms block holds at most one 40 ms NXDN frame: one exact sync word, or several on the frame grid
        if n_frame >= 2 or any(e == 0 for _, e in hits):
            half = 3.125e3 if sr == 2400 else 6.25e3
            return [Emission(family="nxdn", label="nxdn", start_us=pos[0] / fs * 1e6, end_us=dur * 1e6,
                             f_lo_mhz=(c["f_hz"] - half) / 1e6, f_hi_mhz=(c["f_hz"] + half) / 1e6,
                             center_mhz=c["f_hz"] / 1e6, snr_db=snr, confidence=0.85,
                             features={**base, "symbol_rate": sr, "n_frame_sync": len(pos),
                                       "frame_consistent_hits": n_frame, "modulation": f"4FSK {sr:g} sym/s",
                                       "channel_bw_khz": 6.25 if sr == 2400 else 12.5})]
    return None


def _dmr_emission(c, dur, bursts, base):
    kinds = {}
    for b in bursts:
        kinds[b["sync"]] = kinds.get(b["sync"], 0) + 1
    cc = [b["colour_code"] for b in bursts if "colour_code" in b]
    dts = {}
    for b in bursts:
        if "data_type" in b:
            dts[b["data_type"]] = dts.get(b["data_type"], 0) + 1
    csbk = [b for b in bursts if b.get("crc_ok") and "csbko" in b]
    crc_checked = [b for b in bursts if "crc_ok" in b]
    bs = sum(v for k, v in kinds.items() if k.startswith("bs_"))
    ms = sum(v for k, v in kinds.items() if k.startswith("ms_"))
    direct = sum(v for k, v in kinds.items() if k.startswith("direct"))
    role = "bs" if bs >= max(ms, direct) else ("ms" if ms >= direct else "direct")
    feats = {**base, "modulation": "4FSK 4800 sym/s", "channel_bw_khz": 12.5, "n_bursts": len(bursts),
             "sync_types": kinds, "source": {"bs": "base station / repeater", "ms": "mobile station",
                                            "direct": "direct mode (TDMA)"}[role],
             "voice_bursts": sum(v for k, v in kinds.items() if k.endswith("voice")),
             "data_types": dts, "bursts": [{k: v for k, v in b.items()} for b in bursts[:12]]}
    if cc:
        feats["colour_code"] = max(sorted(set(cc)), key=cc.count)     # most frequent; ties: the lowest
        feats["colour_code_votes"] = cc.count(feats["colour_code"])
    if crc_checked:
        feats["crc_checked"] = len(crc_checked)
        feats["crc_ok"] = sum(bool(b["crc_ok"]) for b in crc_checked)
    if csbk:
        feats["csbk"] = [{k: b[k] for k in ("csbko", "feature_set_id", "csbk_last_block", "csbk_data_hex")}
                         for b in csbk[:8]]
    t0 = min(b["t_s"] for b in bursts)
    return Emission(family="dmr", label=f"dmr_{role}", start_us=t0 * 1e6, end_us=dur * 1e6,
                    f_lo_mhz=(c["f_hz"] - 6.25e3) / 1e6, f_hi_mhz=(c["f_hz"] + 6.25e3) / 1e6,
                    center_mhz=c["f_hz"] / 1e6, snr_db=base["snr_db"], confidence=0.9 if cc else 0.8,
                    decoded={"crc_ok": True} if csbk else None, features=feats)


def detect(ch: Channelizer, ctx: dict | None = None) -> list[Emission]:
    """ctx: optional {"band": (lo, hi) analysed band, "exclude_hz": receiver lines to skip}."""
    ctx = ctx or {}
    band = ctx.get("band") or (ch.fc - 0.45 * ch.fs, ch.fc + 0.45 * ch.fs)
    excl = ctx.get("exclude_hz") or []
    dur = ch.L / ch.fs
    out: list[Emission] = []
    for c in candidates(ch, band):
        if any(abs(c["f_hz"] - f) < 3e3 for f in excl):
            continue
        e = _tetra(ch, c, dur) if _in(c["f_hz"], TETRA_RANGES) else None
        if e is not None:
            out.append(e)
            continue
        if _in(c["f_hz"], FSK_RANGES):
            out += _fsk(ch, c, dur) or []
    return out
