#!/usr/bin/env python3
"""Bench: capture one flushed Pluto buffer straight into testdata as SigMF ci16_le (raw ADC counts).
Same flush as lo_step: rx buffer destroyed before the retune, --discard buffers dropped.

    .venv/bin/python bench/live_lo/capture_testdata.py --freq 868 --name LIVE_868MHz [--rate 61.44e6 --bw 56e6
        --n 2097152 --gain 40] [--out testdata/live]
    --trigger 2888.57 [--max-tries 400]: keep reading buffers until one holds >= --min-pulses pulses (1 us envelope
    > median + 12 dB within +-1 MHz of that frequency), e.g. a rotating radar's beam pass.
"""
import argparse, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lo_step import open_pluto, FULL_SCALE

REPO = Path(__file__).resolve().parents[2]

def n_pulses(x, fs, f_off, half=1e6, thr_db=12.0):
    """Pulses (1 us envelope above median + thr_db, gaps < 3 us merged) within f_off +- half."""
    X = np.fft.fft(x.astype(np.complex64)); f = np.fft.fftfreq(len(x), 1 / fs)
    X[np.abs(f - f_off) > half] = 0
    e = np.abs(np.fft.ifft(X)) ** 2
    m = max(1, int(fs * 1e-6)); e = e[:len(e) // m * m].reshape(-1, m).mean(1)
    on = 10 * np.log10(e + 1e-20) > 10 * np.log10(np.median(e) + 1e-20) + thr_db
    idx = np.flatnonzero(on)
    return 0 if len(idx) == 0 else int(1 + np.sum(np.diff(idx) > 3))


def whole_burst(x, fs, f_off, half, min_ms, thr_db=12.0):
    """Longest burst (ms) within f_off +- half that starts and ends inside the buffer, or 0."""
    X = np.fft.fft(x.astype(np.complex64)); f = np.fft.fftfreq(len(x), 1 / fs)
    X[np.abs(f - f_off) > half] = 0
    e = np.abs(np.fft.ifft(X)) ** 2
    m = max(1, int(fs * 100e-6)); e = e[:len(e) // m * m].reshape(-1, m).mean(1)
    on = 10 * np.log10(e + 1e-20) > 10 * np.log10(np.median(e) + 1e-20) + thr_db
    d = np.diff(np.r_[0, on.astype(int), 0]); st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    ok = [(b - a) * 0.1 for a, b in zip(st, en) if a > 0 and b < len(on) and (b - a) * 0.1 >= min_ms]
    return max(ok, default=0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--rate", type=float, default=61.44e6); ap.add_argument("--bw", type=float, default=56e6)
    ap.add_argument("--n", type=int, default=1 << 21); ap.add_argument("--gain", type=float, default=40)
    ap.add_argument("--discard", type=int, default=4); ap.add_argument("--uri", default="usb:1.2.5")
    ap.add_argument("--out", default=str(REPO / "testdata" / "live")); ap.add_argument("--note", default="")
    ap.add_argument("--trigger", type=float, help="MHz: wait for a pulsed signal there")
    ap.add_argument("--min-pulses", type=int, default=10); ap.add_argument("--max-tries", type=int, default=400)
    ap.add_argument("--burst-ms", type=float, help="with --trigger: wait for a whole burst at least this long")
    ap.add_argument("--half-mhz", type=float, default=0.3)
    ap.add_argument("--ppm", type=float, help="receiver ppm to record as bench:receiver_ppm")
    a = ap.parse_args()
    sdr = open_pluto(a.uri, a.rate, a.bw, a.n)
    try:
        sdr.rx_destroy_buffer(); sdr.rx_lo = int(a.freq * 1e6); sdr.rx_hardwaregain_chan0 = a.gain
        for _ in range(a.discard):
            sdr.rx()
        t0 = datetime.now(timezone.utc); x = sdr.rx()
        for k in range(a.max_tries if a.trigger else 0):
            if a.burst_ms:
                ms = whole_burst(x, a.rate, (a.trigger - a.freq) * 1e6, a.half_mhz * 1e6, a.burst_ms)
                if ms:
                    print(f"trigger: whole {ms:.1f} ms burst in buffer {k}")
                    break
            else:
                n = n_pulses(x, a.rate, (a.trigger - a.freq) * 1e6)
                if n >= a.min_pulses:
                    print(f"trigger: {n} pulses in buffer {k}")
                    break
            t0 = datetime.now(timezone.utc); x = sdr.rx()
        else:
            if a.trigger:
                raise SystemExit(f"no trigger in {a.max_tries} buffers")
    finally:
        sdr.rx_destroy_buffer(); del sdr
    iq = np.empty((len(x), 2), "<i2"); iq[:, 0], iq[:, 1] = x.real, x.imag; del x
    clip = float(np.mean(np.abs(iq) >= FULL_SCALE - 1))
    d = Path(a.out) / a.name; d.mkdir(parents=True, exist_ok=True)
    raw = iq.tobytes(); (d / f"{a.name}.sigmf-data").write_bytes(raw)
    meta = {"global": {"core:datatype": "ci16_le", "core:num_channels": 1, "core:recorder": "spectrum_intelligence bench",
                       "core:sample_rate": a.rate, "core:sha512": hashlib.sha512(raw).hexdigest(), "core:version": "1.2.0",
                       "core:description": a.note, "signalvault:device": "pluto", "signalvault:gain_db": a.gain,
                       "signalvault:gain_mode": "manual", "signalvault:rf_bandwidth_hz": a.bw,
                       "bench:flush": f"rx buffer destroyed before retune, {a.discard} buffers discarded",
                       **({"bench:receiver_ppm": a.ppm} if a.ppm is not None else {})},
            "captures": [{"core:datetime": t0.isoformat(), "core:frequency": a.freq * 1e6, "core:sample_start": 0}],
            "annotations": []}
    (d / f"{a.name}.sigmf-meta").write_text(json.dumps(meta, indent=4))
    print(f"{d}  {a.n / a.rate * 1e3:.1f} ms  clip {clip:.1e}  {len(raw) / 1e6:.1f} MB")

if __name__ == "__main__":
    main()
