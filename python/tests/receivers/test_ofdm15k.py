#!/usr/bin/env python3
"""Synthetic self-test for core.ofdm15k (the non-Wi-Fi parts of the research spectrum_survey).

Injects LTE-like 15 kHz OFDM hops (72 subcarriers = 1.08 MHz, 7-symbol 0.5 ms
slots, normal CP) and Wi-Fi-numerology OFDM bursts (64-pt @ 20 MSps, 0.8 us GI)
into noise at 61.44 MSps, then checks:
  * every hop found (blind pass + hop tracker), centre within 30 kHz, start within 20 us
  * Wi-Fi-numerology OFDM never reported as 15 kHz OFDM
  * blind numerology labels the Wi-Fi-like bursts wifi_like_ofdm
  * pure noise -> 0 detections
Expect: all PASS.
"""
import numpy as np
from spectrum_intel.core import ofdm15k as SS, spectral as SP
FS, FC = 61.44e6, 5770e6
rng = np.random.default_rng(7)


def lte_slot(n_sc=72, nfft=4096):
    """One 0.5 ms slot at 61.44 MSps: CP 320 then 6x CP 288, N = 4096 (15 kHz)."""
    out = []
    for i in range(7):
        X = np.zeros(nfft, complex)
        k = np.r_[1:n_sc // 2 + 1, -n_sc // 2:0]                  # DC unused, like LTE DL
        X[k % nfft] = (rng.choice([-1, 1], n_sc) + 1j * rng.choice([-1, 1], n_sc)) / np.sqrt(2)
        s = np.fft.ifft(X) * nfft / np.sqrt(n_sc)                  # unit power
        cp = 320 if i == 0 else 288
        out.append(np.r_[s[-cp:], s])
    return np.concatenate(out)


def wifi_like(n_sym=40):
    """802.11a-numerology OFDM (52 of 64 subcarriers @ 20 MSps, 16-sample GI), resampled to 61.44 MSps."""
    from scipy.signal import resample_poly
    syms = []
    for _ in range(n_sym):
        X = np.zeros(64, complex)
        k = np.r_[1:27, -26:0]
        X[k % 64] = rng.choice([-1, 1], 52) + 1j * rng.choice([-1, 1], 52)
        s = np.fft.ifft(X)
        syms.append(np.r_[s[-16:], s])
    return resample_poly(np.concatenate(syms), 384, 125)


def build(n=2_000_000, snr_db=12):
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    hops = [(5.6e-3, 5774.5e6), (9.6e-3, 5796.6e6), (15.6e-3, 5796.6e6), (19.6e-3, 5774.5e6),
            (25.6e-3, 5752.3e6), (29.6e-3, 5781.1e6)]
    g = np.sqrt(10 ** (snr_db / 10) * 1.08e6 / FS)                  # in-band SNR = snr_db
    t = np.arange(n)
    for t0, f in hops:
        s = lte_slot(); a = int(t0 * FS)
        x[a:a + len(s)] += g * s * np.exp(2j * np.pi * (f - FC) / FS * t[a:a + len(s)])
    wifi = [(2.0e-3, 5745e6), (12.0e-3, 5785e6), (22.0e-3, 5760e6)]
    for t0, f in wifi:
        s = wifi_like(); a = int(t0 * FS)
        s = s / np.sqrt(np.mean(np.abs(s) ** 2)) * 10 ** (15 / 20)
        x[a:a + len(s)] += s * np.exp(2j * np.pi * (f - FC) / FS * t[a:a + len(s)])
    return x, hops, wifi


def run_15k(x):
    ch = SP.Channelizer(x, FS)
    hits, stats = SS.ofdm15k_scan(ch, FC)
    ob = SS.merge_ofdm15k([b for b in (SS.refine_ofdm15k(ch, h) for h in hits) if b is not None])
    ob = SS.merge_ofdm15k(ob + SS.track_hops(ch, FC, ob, stats, len(x) / FS))
    return ch, ob


def main():
    ok = True
    x, hops, wifi = build()
    ch, ob = run_15k(x)
    for t0, f in hops:
        m = [b for b in ob if abs(b["start_s"] - t0) < 20e-6 and abs(FC + b["f_off"] - f) < 30e3]
        good = len(m) == 1 and SS.label_ofdm15k(m[0])[0] == "ocusync_uplink_like"
        print(f"{'PASS' if good else 'FAIL'} hop {t0*1e3:5.1f} ms {f/1e6:.1f} MHz"
              + (f"  -> {(FC + m[0]['f_off'])/1e6:.4f} MHz, start {m[0]['start_s']*1e3:.3f} ms, "
                 f"bw {m[0]['bw_hz']/1e6:.2f} MHz{' (tracked)' if m[0].get('tracked') else ''}" if m else ""))
        ok &= good
    extra = [b for b in ob if not any(abs(b["start_s"] - t0) < 20e-6 for t0, _ in hops)]
    print(f"{'PASS' if not extra else 'FAIL'} false 15 kHz detections: {len(extra)}")
    ok &= not extra
    for t0, f in wifi:
        nm = SS.ofdm_numerology(ch, f - FC, 9e6, t0 + 20e-6, t0 + 150e-6)
        lab = SS.classify_numerology(nm[0])[0] if nm else None
        good = lab == "wifi_like_ofdm"
        print(f"{'PASS' if good else 'FAIL'} numerology of Wi-Fi-like burst at {f/1e6:.0f} MHz: {lab}"
              + (f" (Tu {nm[0]*1e6:.2f} us)" if nm else ""))
        ok &= good
    noise = (rng.standard_normal(2_000_000) + 1j * rng.standard_normal(2_000_000)) / np.sqrt(2)
    _, ob_n = run_15k(noise)
    print(f"{'PASS' if not ob_n else 'FAIL'} pure noise -> {len(ob_n)} 15 kHz detections")
    ok &= not ob_n
    print("\nALL PASS" if ok else "\nSOME FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
