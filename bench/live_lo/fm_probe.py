#!/usr/bin/env python3
"""Bench: FM-demodulate one carrier of a dwell capture and characterise the baseband (audio/data).

    .venv/bin/python bench/live_lo/fm_probe.py c420 420.0029 [--bw 12.5e3] [--wav out.wav]
Prints channel SNR, deviation, the audio spectrum in bands and its strongest lines, the power envelope
(keyed on/off?), and a crude FSK test (two-tone or two-level baseband histogram).
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
from scipy.io import wavfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lo_step import SCRATCH

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("name"); ap.add_argument("f", type=float)
    ap.add_argument("--bw", type=float, default=12.5e3); ap.add_argument("--wav")
    a = ap.parse_args()
    m = json.loads((SCRATCH / "dwell" / f"{a.name}.json").read_text()); s = np.load(SCRATCH / "dwell" / f"{a.name}.npy")
    x = s[:, 0].astype(np.float32) + 1j * s[:, 1].astype(np.float32); fs = m["rate"]; fc = m["freq_mhz"] * 1e6
    y = x * np.exp(-2j * np.pi * (a.f * 1e6 - fc) / fs * np.arange(len(x))).astype(np.complex64); del x
    y = ss.resample_poly(y, 1, int(fs // 48e3)); fo = fs / int(fs // 48e3)
    Y = np.abs(np.fft.fftshift(np.fft.fft(y))) ** 2; fy = np.fft.fftshift(np.fft.fftfreq(len(y), 1 / fo))
    snr = 10 * np.log10(Y[np.abs(fy) < a.bw / 2].mean() / Y[np.abs(np.abs(fy) - 18e3) < 3e3].mean())
    occ = fy[(Y > Y.max() / 100) & (np.abs(fy) < 20e3)]
    yf = ss.lfilter(*ss.butter(5, a.bw / 2 / (fo / 2)), y)
    env = np.abs(yf); b = int(5e-3 * fo); e = 20 * np.log10(env[:len(env) // b * b].reshape(-1, b).mean(1) + 1e-9)
    d = np.angle(yf[1:] * np.conj(yf[:-1])) * fo / (2 * np.pi)
    A = np.abs(np.fft.rfft((d - d.mean()) * np.hanning(len(d)))); fa = np.fft.rfftfreq(len(d), 1 / fo)
    pw = lambda lo, hi: 10 * np.log10(np.mean(A[(fa > lo) & (fa < hi)] ** 2))
    s_ = (fa > 50) & (fa < 6000); An = A[s_] / np.median(A[s_]); top = np.argsort(An)[-8:][::-1]
    print(f"{a.f} MHz: channel SNR {snr:.1f} dB, -20 dB occupied {occ.min()/1e3:+.1f}..{occ.max()/1e3:+.1f} kHz, "
          f"dev rms {np.std(d)/1e3:.2f} kHz, mean offset {np.mean(d):+.0f} Hz")
    print(f"  envelope (5 ms) dB: min {e.min():.1f} p10 {np.percentile(e,10):.1f} med {np.median(e):.1f} max {e.max():.1f}")
    print(f"  audio band dB: 30-250 {pw(30,250):.0f}  300-1000 {pw(300,1000):.0f}  1000-2500 {pw(1000,2500):.0f}  3000-6000 {pw(3000,6000):.0f}")
    print(f"  strongest baseband lines (Hz, x median): {[(round(fa[s_][i], 1), round(An[i])) for i in top]}")
    dl = ss.lfilter(*ss.butter(4, 3e3 / (fo / 2)), d)
    h, edges = np.histogram(dl, bins=30); print("  baseband histogram (lowpassed):", h.tolist(), f"range {edges[0]:.0f}..{edges[-1]:.0f} Hz")
    if a.wav:
        au = ss.lfilter(*ss.butter(4, [300 / (fo / 2), 3400 / (fo / 2)], "band"), d)
        wavfile.write(a.wav, int(fo), (au / (np.abs(au).max() + 1e-9) * 30000).astype(np.int16)); print("  wav:", a.wav)

if __name__ == "__main__":
    main()
