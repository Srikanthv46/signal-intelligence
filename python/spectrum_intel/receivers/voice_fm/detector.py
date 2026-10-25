"""Narrowband FM voice (walkie-talkies, PMR, land-mobile) detector.

    detect(x, fs, fc, ctx=None) -> list[Emission]   family "voice_fm"

1. Candidates: the capture's spectrum averaged into 1 kHz bins; runs >= 10 dB over the
   median floor and <= 40 kHz wide, inside the voice ranges (136-174, 300-520, 806-870 MHz).
2. Per candidate, the +-10 kHz channel at 48 kSps: 5 ms frame power -> keyed time
   (>= 10 dB over the floor), constant envelope during it (FM), FM discriminator ->
   carrier frequency, peak deviation (<= 6 kHz: narrowband), voice-band shape.
3. Narrowband *data* (2/4-FSK: POCSAG, DMR, dPMR) has a platykurtic voice band (two or
   four levels) and is left to the band rules; an unmodulated carrier is left to the CW rule.
4. Sub-audible signalling on the longest keyed segment: DCS (Golay-verified) or CTCSS.
5. Channel plan match (PMR446 / LPD433 / FRS) with the receiver ppm if the pipeline has it.
6. Clipping products and IQ images of a stronger carrier (k x its offset from the LO,
   k = -1, +-3, +-5, +-7, 9) are dropped and noted on the carrier.
Depends on numpy/scipy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np
from scipy import signal as ss

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.core.emission import Emission

from . import afsk
from .channels import match_channel
from .subaudio import subaudible

VOICE_RANGES = [(118e6, 174e6), (225e6, 520e6), (806e6, 870e6)]
# AM voice (research, signal gaps S13): aeronautical bands (ICAO Annex 10 Vol V; India NFAP-2022)
AIR_BANDS = [(118e6, 137e6, "VHF air band"), (225e6, 400e6, "UHF military air band")]
AIR_EMERGENCY = {121.5e6: "aeronautical emergency 121.5 MHz", 243.0e6: "military emergency 243.0 MHz"}
AM_FRAME = 0.05               # a 20 ms frame with speech: envelope speech-band rms >= 5 % of the carrier
AM_MIN_RMS = 0.04
AM_SPEECH_FRAC = 0.31         # share of the slice's envelope noise (+-10 kHz) inside 300-3400 Hz
RES = 1000.0                  # candidate search resolution, Hz
MID = 48000.0                 # channel analysis rate
FRAME_S = 0.005
MAX_DEV = 6000.0
KEYED_DB = 10.0
FSK_RATES = (2400.0, 4800.0)  # 4FSK digital voice/data symbol rates (DMR, P25, NXDN 9600/4800, dPMR)
FSK_LINE = 6.0                # symbol-rate line over its neighbourhood: data, not speech
ENV_CV_FLOOR = 0.12           # FM envelope cv allowed whatever the SNR
ENV_CV_EXCESS = 1.6           # ... and up to this multiple of what the noise alone gives


def _runs(m: np.ndarray) -> list[tuple[int, int]]:
    d = np.diff(np.r_[0, m.astype(np.int8), 0])
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def candidates(ch: Channelizer, max_n: int = 12) -> list[dict]:
    fs, fc, L = ch.fs, ch.fc, ch.L
    g = max(1, int(round(RES * L / fs)))
    order = np.argsort(ch.f)
    P = np.abs(ch.X[order]) ** 2
    f = ch.f[order]
    n = len(P) // g
    Pg, fg = P[:n * g].reshape(n, g).mean(1), f[:n * g].reshape(n, g).mean(1)
    fa = fc + fg
    inr = np.zeros(n, bool)
    for a, b in VOICE_RANGES:
        inr |= (fa >= a) & (fa <= b)
    inr &= np.abs(fg) < 0.45 * fs
    if inr.sum() < 20:
        return []
    floor = float(np.median(Pg[np.abs(fg) < 0.45 * fs]))
    L_db = 10 * np.log10(Pg / floor + 1e-30)
    # narrow peaks: >= 10 dB over the floor and over their surroundings, <= 40 kHz wide at half
    # prominence (a strong or clipped carrier's skirt can stay above the floor far wider than that)
    pk, pr = ss.find_peaks(L_db, height=KEYED_DB, prominence=KEYED_DB, width=(None, 40e3 / RES),
                           rel_height=0.5, distance=max(1, int(3e3 / RES)))
    out = []
    for i, w in zip(pk, pr["widths"]):
        if not inr[i]:
            continue
        a, b = max(0, i - 4), min(n, i + 5)
        p = Pg[a:b]
        out.append({"f_hz": float(fc + np.sum(fg[a:b] * p) / np.sum(p)), "peak_db": float(L_db[i]),
                    "width_hz": float(w * RES), "noise_density": floor / (fs * L)})
    # noise_density: noise mean|y|^2 in B Hz = floor_bin * B / (fs * L) (Parseval)
    return sorted(out, key=lambda c: -c["peak_db"])[:max_n]


def _kurtosis(v: np.ndarray) -> float:
    v = v - v.mean()
    s = np.mean(v ** 2)
    return float(np.mean(v ** 4) / (s ** 2 + 1e-30))


def measure(ch: Channelizer, c: dict) -> dict | None:
    """Keyed time, envelope, deviation, voice band and sub-audio of one candidate."""
    y, fs_y = ch.slice(c["f_hz"], 10e3, MID)
    n_fr = max(1, int(FRAME_S * fs_y))
    nf = len(y) // n_fr
    if nf < 4:
        return None
    p = np.mean(np.abs(y[:nf * n_fr].reshape(nf, n_fr)) ** 2, 1)
    noise = c["noise_density"] * 20e3                             # the +-10 kHz slice's noise power
    keyed = p > noise * 10 ** (KEYED_DB / 10)
    # keyed segments; gaps under 30 ms bridged (syllable dips never drop an FM carrier, fades can)
    k = keyed.copy()
    for a, b in _runs(~keyed):
        if 0 < a and b < nf and (b - a) * FRAME_S < 0.03:
            k[a:b] = True
    segs = [(a, b) for a, b in _runs(k) if (b - a) * FRAME_S >= 0.01]
    dur = len(y) / fs_y
    keyed_s = sum(b - a for a, b in segs) * FRAME_S
    if not segs or keyed_s < min(0.04, 0.5 * dur):
        return None
    idx = np.concatenate([np.arange(a * n_fr, b * n_fr) for a, b in segs])
    env = np.abs(y[idx])
    cv = float(np.std(env) / (np.mean(env) + 1e-30))
    d = np.angle(y[1:] * np.conj(y[:-1])) * (fs_y / (2 * np.pi))
    d = np.r_[d[:1], d]
    lp = ss.butter(4, 4000, fs=fs_y, output="sos")
    dl = ss.sosfiltfilt(lp, d)
    dk = dl[idx]
    mean_off = float(np.mean(dk))
    dev = float(np.percentile(np.abs(dk - mean_off), 99.5))
    bp = ss.butter(4, [300, 3400], "bandpass", fs=fs_y, output="sos")
    v = ss.sosfiltfilt(bp, d)[idx]
    nv = len(v) // int(0.02 * fs_y)
    vf = v[:nv * int(0.02 * fs_y)].reshape(nv, -1) if nv else np.zeros((0, 1))
    vrms = np.sqrt(np.mean(vf ** 2, 1)) if nv else np.zeros(0)
    # AM: the envelope's speech band relative to the carrier (an FM carrier's envelope only carries noise)
    ea = np.abs(y)
    ae = ss.sosfiltfilt(bp, ea - float(np.mean(env)))[idx] / (float(np.mean(env)) + 1e-30)
    af = ae[:nv * int(0.02 * fs_y)].reshape(nv, -1) if nv else np.zeros((0, 1))
    arms = np.sqrt(np.mean(af ** 2, 1)) if nv else np.zeros(0)
    a, b = max(segs, key=lambda s: s[1] - s[0])
    sub = subaudible(d[a * n_fr:b * n_fr], fs_y)
    fsk_line = _symbol_line(dl[a * n_fr:b * n_fr], fs_y)
    snr = 10 * np.log10(np.median(p[keyed]) / noise) if keyed.any() else None
    return {"f_hz": c["f_hz"] + mean_off, "segs_s": [(a * FRAME_S, b * FRAME_S) for a, b in segs],
            "keyed_s": keyed_s, "duration_s": dur, "envelope_cv": cv, "deviation_hz": dev,
            "voice_rms_hz": float(np.sqrt(np.mean(v ** 2))), "voice_kurtosis": _kurtosis(v),
            "voice_activity": float(np.mean(vrms > 250)) if nv else 0.0, "subaudible": sub, "fsk_line": fsk_line,
            "am_voice_rms": float(np.sqrt(np.mean(ae ** 2))), "am_activity": float(np.mean(arms > AM_FRAME)) if nv else 0.0,
            "packets": afsk.decode(d, fs_y),
            "snr_db": None if snr is None else float(snr), "peak_db": c["peak_db"]}


def _symbol_line(dl: np.ndarray, fs: float) -> dict:
    """Strongest FSK symbol-rate line in the squared derivative of the (4 kHz low-passed) discriminator output,
    over the median of its +-40..400 Hz neighbourhood. Digital FSK has it; speech and sub-audible tones do not."""
    if len(dl) < int(0.02 * fs):
        return {"rate_hz": None, "ratio": 0.0}
    e = np.diff(dl) ** 2
    E = np.abs(np.fft.rfft(e - e.mean()))
    fe = np.fft.rfftfreq(len(e), 1 / fs)
    best = {"rate_hz": None, "ratio": 0.0}
    for r in FSK_RATES:
        on, bg = np.abs(fe - r) <= 15, (np.abs(fe - r) > 40) & (np.abs(fe - r) < 400)
        if on.any() and bg.any():
            q = float(E[on].max() / (np.median(E[bg]) + 1e-30))
            if q > best["ratio"]:
                best = {"rate_hz": r, "ratio": q}
    return best


def _noise_cv(snr_db: float) -> float:
    """Envelope coefficient of variation of a constant-envelope carrier at this SNR (noise alone): ~1/sqrt(2 SNR)."""
    return 1.0 / np.sqrt(2 * 10 ** (snr_db / 10))


def _artifact_of(m: dict, stronger: list[dict], fc: float) -> dict | None:
    """A product of a stronger carrier at k x its offset from the LO. ADC clipping of a
    carrier at offset f puts copies at (4j+1) f (-3f, 5f, -7f, 9f ...) with |k| x its
    deviation, only ~10 dB weaker when clipping is hard; the IQ image (k = -1) has the same
    deviation, so it must also be >= 15 dB weaker."""
    fo = m["f_hz"] - fc
    for s in stronger:
        so = s["f_hz"] - fc
        if abs(so) < 1e3:
            continue
        for k in (-1, 3, -3, 5, -5, 7, -7, 9, -9, 11):
            if abs(fo - k * so) >= 300 + 60 * abs(k):
                continue
            if k == -1:
                ok = s["peak_db"] - m["peak_db"] >= 15
            else:
                r = m["deviation_hz"] / max(s["deviation_hz"], 1.0)
                ok = 0.6 * abs(k) <= r <= 1.4 * abs(k)
            if ok:
                return {"of_mhz": round(s["f_hz"] / 1e6, 5), "k": k}
    return None


def detect(x: np.ndarray, fs: float, fc: float, ctx: dict | None = None) -> list[Emission]:
    """ctx: optional {"ch": core.dsp.Channelizer of x, "ppm": receiver frequency error,
    "exclude_hz": the receiver's own clock lines (skipped)}."""
    ctx = ctx or {}
    if not any(fc - fs / 2 < b and fc + fs / 2 > a for a, b in VOICE_RANGES):
        return []
    ch = ctx.get("ch") or Channelizer(np.asarray(x, np.complex64), fs, fc)
    ppm = float(ctx.get("ppm") or 0.0)
    excl = ctx.get("exclude_hz") or []
    cands = [c for c in candidates(ch) if not any(abs(c["f_hz"] - f) < 2e3 for f in excl)]
    measured = [m for m in (measure(ch, c) for c in cands) if m is not None]   # strongest first
    kept: list[dict] = []
    for i, m in enumerate(measured):
        art = _artifact_of(m, measured[:i], fc)
        if art is None:
            kept.append(m)
            continue
        host = next(s for s in measured[:i] if round(s["f_hz"] / 1e6, 5) == art["of_mhz"])
        host.setdefault("artifacts", []).append({"f_mhz": round(m["f_hz"] / 1e6, 5), **art})
    return [e for e in (_emission(m, ppm) for m in kept) if e is not None]


def _noise_fm_rms(snr_db: float) -> float:
    """Speech-band (300-3400 Hz) rms of the discriminator output for an unmodulated carrier at this SNR: phase noise
    1 / (2 SNR) spread over the +-10 kHz slice, frequency noise density f^2 times that."""
    nd = 1.0 / (2 * 10 ** (snr_db / 10)) / 20e3
    return float(np.sqrt(2 * nd * (3400.0 ** 3 - 300.0 ** 3) / 3))


def air_channel(f_hz: float, ppm: float = 0.0) -> dict:
    """Aeronautical service and nearest 25 kHz channel (8.33 kHz channel names are not resolved)."""
    f = f_hz * (1 + ppm * 1e-6)
    out = {}
    for a, b, name in AIR_BANDS:
        if a <= f < b:
            out["service"] = name
            nom = round(f / 25e3) * 25e3
            out["channel_mhz"] = round(nom / 1e6, 4)
            out["offset_hz"] = round(f - nom)
            break
    for fe, name in AIR_EMERGENCY.items():
        if abs(f - fe) <= 5e3:
            out["emergency"] = name
    return out


def _am_emission(m: dict, ppm: float) -> Emission | None:
    """Full-carrier AM voice (research): a carrier whose envelope carries speech while its phase stays quiet."""
    snr = m["snr_db"]
    air = air_channel(m["f_hz"], ppm)
    if not air or snr is None or snr < 10 or m["fsk_line"]["ratio"] >= FSK_LINE:
        return None                                        # AM voice only in the aeronautical bands (OOK remotes
                                                           # and pulsed data elsewhere also move the envelope)
    noise_am = _noise_cv(snr) * np.sqrt(AM_SPEECH_FRAC)
    if m["am_voice_rms"] < max(AM_MIN_RMS, 3 * noise_am) or m["am_activity"] < 0.1:
        return None
    if m["voice_rms_hz"] > 3 * _noise_fm_rms(snr) + 150:
        return None                                        # the phase moves: FM (or noise), not AM
    conf, why = 0.7, ["carrier with speech-band envelope modulation",
                      f"envelope speech rms {m['am_voice_rms'] * 100:.0f} % of the carrier "
                      f"(noise alone {noise_am * 100:.1f} %)",
                      f"phase quiet: discriminator speech-band rms {m['voice_rms_hz']:.0f} Hz",
                      f"speech in {m['am_activity'] * 100:.0f} % of 20 ms frames", f"keyed {m['keyed_s']:.2f} s"]
    why.append(air.get("emergency") or f"{air['service']}, channel {air['channel_mhz']} MHz "
                                       f"({air['offset_hz']:+d} Hz)")
    if m["am_activity"] >= 0.2:
        conf += 0.1
    feats = {"center_mhz": round(m["f_hz"] / 1e6, 6), "modulation_rms": round(m["am_voice_rms"], 3),
             "channel_bw_khz": 25.0, "keyed_s": round(m["keyed_s"], 3),
             "keyed_fraction": round(m["keyed_s"] / m["duration_s"], 3),
             "envelope_cv": round(m["envelope_cv"], 3), "voice_activity": round(m["am_activity"], 3),
             "phase_rms_hz": round(m["voice_rms_hz"]), "snr_db": round(snr, 1), "evidence": why, **air}
    if m.get("artifacts"):
        feats["suppressed_artifacts"] = m["artifacts"]
    segs = m["segs_s"]
    return Emission(family="voice_fm", label="am_voice", start_us=segs[0][0] * 1e6, end_us=segs[-1][1] * 1e6,
                    f_lo_mhz=(m["f_hz"] - 4e3) / 1e6, f_hi_mhz=(m["f_hz"] + 4e3) / 1e6,
                    center_mhz=m["f_hz"] / 1e6, snr_db=feats["snr_db"], confidence=round(min(conf, 0.9), 2),
                    decoded=None, features={**feats, "segments_s": [[round(a, 3), round(b, 3)] for a, b in segs]})


def _packet_emission(m: dict) -> Emission | None:
    """AX.25 / APRS packets (CRC-valid) carried by the FM channel (research, signal gaps S19)."""
    pk = m.get("packets") or []
    if not pk:
        return None
    aprs = any(p.get("aprs") for p in pk)
    srcs = sorted({p["src"] for p in pk})
    feats = {"center_mhz": round(m["f_hz"] / 1e6, 6), "deviation_hz": round(m["deviation_hz"]),
             "channel_bw_khz": 12.5 if m["deviation_hz"] <= 3200 else 25.0, "keyed_s": round(m["keyed_s"], 3),
             "snr_db": None if m["snr_db"] is None else round(m["snr_db"], 1), "stations": srcs, "packets": pk,
             "evidence": [f"{len(pk)} AX.25 UI frame(s) with a valid CRC-16/X.25 (Bell 202 AFSK 1200 bit/s)",
                          f"stations {srcs}"]}
    segs = m["segs_s"]
    half = 6.25e3 if m["deviation_hz"] <= 3200 else 12.5e3
    return Emission(family="voice_fm", label="aprs" if aprs else "ax25_packet", start_us=segs[0][0] * 1e6,
                    end_us=segs[-1][1] * 1e6, f_lo_mhz=(m["f_hz"] - half) / 1e6, f_hi_mhz=(m["f_hz"] + half) / 1e6,
                    center_mhz=m["f_hz"] / 1e6, snr_db=feats["snr_db"], confidence=0.97,
                    decoded={"crc_ok": True, "packets": len(pk)},
                    features={**feats, "segments_s": [[round(a, 3), round(b, 3)] for a, b in segs]})


def _emission(m: dict, ppm: float) -> Emission | None:
    pk = _packet_emission(m)
    if pk is not None:
        return pk
    am = _am_emission(m, ppm)
    if am is not None:
        return am
    sub = m["subaudible"] or {}
    if m["envelope_cv"] > 0.35 or m["deviation_hz"] > MAX_DEV:
        return None                                        # not a constant-envelope narrowband FM carrier
    if not sub and m["snr_db"] is not None and m["envelope_cv"] > max(ENV_CV_FLOOR, ENV_CV_EXCESS * _noise_cv(m["snr_db"])):
        return None                                        # envelope varies beyond noise: two tones / AM, not FM
    modulated = bool(sub) or m["voice_activity"] >= 0.05 or m["deviation_hz"] >= 300
    if not modulated:
        return None                                        # unmodulated carrier: CW, not voice
    if not sub and m["voice_kurtosis"] < 2.2 and m["voice_rms_hz"] > 500:
        return None                                        # two/four-level FSK data, not speech
    if m["fsk_line"]["ratio"] >= FSK_LINE:
        return None                                        # FSK symbol-rate line: digital voice/data, not analog
    chan = match_channel(m["f_hz"], ppm)
    conf, why = 0.5, ["constant-envelope narrowband FM carrier",
                      f"peak deviation {m['deviation_hz']:.0f} Hz", f"keyed {m['keyed_s']:.2f} s"]
    if sub.get("type") == "DCS":
        conf += 0.3; why.append(f"DCS {sub['dcs_code']} ({sub['dcs_words_verified']} Golay-verified words)")
    elif sub.get("ctcss_hz"):
        conf += 0.25; why.append(f"CTCSS {sub['ctcss_hz']} Hz (+{sub['ctcss_snr_db']} dB)")
    elif sub:
        conf += 0.1; why.append(f"sub-audible tone {sub['ctcss_measured_hz']} Hz (non-standard)")
    if m["voice_activity"] >= 0.1:
        conf += 0.1; why.append(f"speech-band activity {m['voice_activity'] * 100:.0f}% of 20 ms frames")
    if chan:
        conf += 0.1; why.append(f"{chan['service']} channel {chan['channel']} ({chan['offset_hz']:+d} Hz)")
    if m["keyed_s"] >= 0.3:
        conf += 0.05
    conf = min(conf, 0.97)
    chan_bw = (chan or {}).get("channel_bw_khz") or (12.5 if m["deviation_hz"] <= 3200 else 25.0)
    half = min(chan_bw * 1e3 / 2, m["deviation_hz"] + 3400)
    feats = {"center_mhz": round(m["f_hz"] / 1e6, 6), "deviation_hz": round(m["deviation_hz"]),
             "channel_bw_khz": chan_bw, "keyed_s": round(m["keyed_s"], 3),
             "keyed_fraction": round(m["keyed_s"] / m["duration_s"], 3),
             "envelope_cv": round(m["envelope_cv"], 3), "voice_activity": round(m["voice_activity"], 3),
             "voice_rms_hz": round(m["voice_rms_hz"]), "voice_kurtosis": round(m["voice_kurtosis"], 2),
             "snr_db": None if m["snr_db"] is None else round(m["snr_db"], 1), "evidence": why,
             **(chan or {}), **sub}
    if m.get("artifacts"):
        feats["suppressed_artifacts"] = m["artifacts"]
    label = "nbfm_voice" + (f"_{chan['service'].lower().replace('/', '_')}_ch{chan['channel']}" if chan else "")
    segs = m["segs_s"]
    return Emission(family="voice_fm", label=label, start_us=segs[0][0] * 1e6, end_us=segs[-1][1] * 1e6,
                    f_lo_mhz=(m["f_hz"] - half) / 1e6, f_hi_mhz=(m["f_hz"] + half) / 1e6,
                    center_mhz=m["f_hz"] / 1e6, snr_db=feats["snr_db"], confidence=round(conf, 2),
                    decoded=sub or None, features={**feats, "segments_s": [[round(a, 3), round(b, 3)] for a, b in segs]})
