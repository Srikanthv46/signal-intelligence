#!/usr/bin/env python3
"""Bench: capture one LO from the Pluto, keep only post-retune snapshots, analyse them.

One invocation = one LO. The radio buffer is destroyed before the retune, the first --discard buffers after
it are read and thrown away (their spectra are logged so we can see when the new LO is settled), then
--keep snapshots are captured --gap seconds apart. If the ADC clips, the gain is lowered and the LO redone.

Each kept snapshot is analysed with spectrum_intel.api.Engine. The compact responses go to
<out>/<MHz>/snapNN.json; the IQ goes to <scratch>/<MHz>/snapNN.npy (ci16) only for follow-up research on this
LO and is deleted with --clean once we move on.

    .venv/bin/python bench/live_lo/lo_step.py --freq 300 [--gain 40] [--keep 3] [--discard 4]
    .venv/bin/python bench/live_lo/lo_step.py --clean 300      # drop that LO's IQ
"""
from __future__ import annotations

import argparse
import gc
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
SCRATCH = Path("/private/tmp/claude-501/-Users-srikanth-Documents-modules-spectrum-intelligence/"
               "c94d801a-125c-43e1-8ae4-81b504e1452f/scratchpad/lo_iq")
FULL_SCALE = 2048          # Pluto: 12-bit ADC, sign-extended into int16
CLIP_FRAC = 1e-3           # more than 0.1 % of samples at the rails -> back the gain off


def open_pluto(uri, rate, bw, n):
    sys.path.insert(0, "/Users/srikanth/Documents/modules/signal-vault/src")
    from signal_vault.drivers.libiio import load_adi
    adi = load_adi()
    sdr = adi.Pluto(uri=uri)
    sdr.sample_rate = int(rate)
    sdr.rx_rf_bandwidth = int(bw)
    sdr.rx_buffer_size = int(n)
    sdr.gain_control_mode_chan0 = "manual"
    return sdr


def psd_db(x, nfft=1024):
    m = len(x) // nfft
    X = np.fft.fftshift(np.fft.fft(x[: m * nfft].reshape(m, nfft) * np.hanning(nfft), axis=1), axes=1)
    return 10 * np.log10(np.mean(np.abs(X) ** 2, axis=0) + 1e-12)


def capture(sdr, f_hz, gain, discard, keep, gap):
    sdr.rx_destroy_buffer()                 # no buffer from the previous LO survives the retune
    sdr.rx_lo = int(f_hz)
    sdr.rx_hardwaregain_chan0 = gain
    time.sleep(0.05)                        # LO lock
    log, prev = [], None
    for i in range(discard):
        x = sdr.rx()
        p = psd_db(x)
        log.append(dict(buf=i, rms=float(np.sqrt(np.mean(np.abs(x) ** 2))),
                        psd_change_db=None if prev is None else float(np.median(np.abs(p - prev)))))
        prev = p
        del x
    snaps = []
    for k in range(keep):
        if k:
            time.sleep(gap)
        t0 = datetime.now(timezone.utc)
        x = sdr.rx()
        p = psd_db(x)
        clip = float(np.mean((np.abs(x.real) >= FULL_SCALE - 1) | (np.abs(x.imag) >= FULL_SCALE - 1)))
        log.append(dict(buf=f"keep{k}", rms=float(np.sqrt(np.mean(np.abs(x) ** 2))), clip=clip,
                        psd_change_db=float(np.median(np.abs(p - prev)))))
        prev = p
        iq16 = np.empty((len(x), 2), np.int16)
        iq16[:, 0], iq16[:, 1] = x.real, x.imag
        del x
        snaps.append((t0, iq16, clip))
    return snaps, log


def summarize(resp):
    rows = []
    for s in resp.get("signals", []):
        f, d = s.get("frequency", {}), s.get("decode", {})
        rows.append(dict(
            f_mhz=round(f.get("center_hz", 0) / 1e6, 3), bw_khz=round((f.get("bandwidth_hz") or 0) / 1e3, 1),
            tech=s.get("technology"), label=s.get("label"), conf=s.get("confidence"),
            snr=s.get("snr_db"), dec=d.get("status"), ids=s.get("identifiers") or None,
            reason=d.get("reason") if d.get("status") not in ("decoded",) else None))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, help="LO in MHz")
    ap.add_argument("--rate", type=float, default=61.44e6)
    ap.add_argument("--bw", type=float, default=56e6)
    ap.add_argument("--n", type=int, default=1 << 21)
    ap.add_argument("--gain", type=float, default=40)
    ap.add_argument("--min-gain", type=float, default=0)
    ap.add_argument("--discard", type=int, default=4)
    ap.add_argument("--keep", type=int, default=3)
    ap.add_argument("--gap", type=float, default=0.3)
    ap.add_argument("--ppm", type=float, default=7.4)
    ap.add_argument("--uri", default="usb:1.2.5")
    ap.add_argument("--dither", type=float, default=1.0, help="MHz; also capture LO+dither for the harmonic-image "
                    "map (0 = off)")
    ap.add_argument("--no-analyze", action="store_true")
    ap.add_argument("--clean", type=float, metavar="MHZ", help="delete the stored IQ of that LO and exit")
    a = ap.parse_args()

    if a.clean is not None:
        for d in (SCRATCH / f"{a.clean:g}", SCRATCH / f"{a.clean + a.dither:g}"):
            shutil.rmtree(d, ignore_errors=True)
            print(f"removed {d}")
        return

    tag = f"{a.freq:g}"
    out, scr = OUT / tag, SCRATCH / tag
    out.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(scr, ignore_errors=True)
    scr.mkdir(parents=True)

    sdr = open_pluto(a.uri, a.rate, a.bw, a.n)
    gain = a.gain
    try:
        while True:
            snaps, log = capture(sdr, a.freq * 1e6, gain, a.discard, a.keep, a.gap)
            worst = max(c for _, _, c in snaps)
            if worst <= CLIP_FRAC or gain - 10 < a.min_gain:
                break
            print(f"gain {gain:g} dB clips {worst:.2%} -> {gain - 10:g} dB")
            del snaps
            gain -= 10
        dsnap = None
        if a.dither:
            dsnap, _ = capture(sdr, (a.freq + a.dither) * 1e6, gain, a.discard, 1, 0)
    finally:
        sdr.rx_destroy_buffer()
        del sdr
    gc.collect()

    meta = dict(freq_mhz=a.freq, rate=a.rate, bw=a.bw, n=a.n, gain_db=gain, ppm=a.ppm, buffers=log,
                snapshots=[dict(file=str(scr / f"snap{k:02d}.npy"), utc=t.isoformat(), clip=c)
                           for k, (t, _, c) in enumerate(snaps)])
    for k, (_, iq16, _) in enumerate(snaps):
        np.save(scr / f"snap{k:02d}.npy", iq16)
    (out / "capture.json").write_text(json.dumps(meta, indent=1))
    if dsnap:
        dtag = f"{a.freq + a.dither:g}"
        (OUT / dtag).mkdir(parents=True, exist_ok=True); (SCRATCH / dtag).mkdir(parents=True, exist_ok=True)
        np.save(SCRATCH / dtag / "snap00.npy", dsnap[0][1])
        (OUT / dtag / "capture.json").write_text(json.dumps(dict(meta, freq_mhz=a.freq + a.dither, dither_of=a.freq,
            buffers=None, snapshots=[dict(file=str(SCRATCH / dtag / "snap00.npy"), utc=dsnap[0][0].isoformat(),
                                          clip=dsnap[0][2])])))
        del dsnap
        import subprocess
        subprocess.run([sys.executable, str(HERE / "harmonic_map.py"), tag, dtag])
    print(json.dumps(dict(gain_db=gain, buffers=log), indent=None))

    if a.no_analyze:
        return
    from spectrum_intel.api import Engine, IQBlock
    eng = Engine(detail="compact")
    for k, (t0, iq16, _) in enumerate(snaps):
        x = iq16[:, 0].astype(np.float32) + 1j * iq16[:, 1].astype(np.float32)
        snaps[k] = None                       # release the int16 copy as we go
        t = time.time()
        resp = eng.analyze(IQBlock(x, sample_rate_hz=a.rate, center_frequency_hz=a.freq * 1e6, bandwidth_hz=a.bw,
                                   start_time=t0.isoformat().replace("+00:00", "Z"), gain_db=gain,
                                   receiver_ppm=a.ppm, source_id="pluto-bench"))
        del x
        (out / f"snap{k:02d}.json").write_text(json.dumps(resp))
        print(f"\n== snap{k:02d}  {time.time() - t:.1f}s  status={resp.get('status')}  "
              f"warnings={[w.get('code') for w in resp.get('warnings', [])]}")
        for r in summarize(resp):
            print("  ", {k2: v for k2, v in r.items() if v is not None})
        del resp
        gc.collect()


if __name__ == "__main__":
    main()
