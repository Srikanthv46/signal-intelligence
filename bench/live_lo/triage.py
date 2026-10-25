#!/usr/bin/env python3
"""Bench: triage every Engine detection of one LO into SELF / IMAGE / REAL.

Tests (all from this LO's stored captures; see FINDINGS.md for why):
  window : detection lies wholly inside the capture-locked self-noise window (t < 12.4 ms of the buffer);
  clock  : raw centre (no ppm correction) within 1.5 kHz of n x 61.44 MHz or n x 40 MHz -> our own clocks;
  image  : LO+delta dither capture -> display shift 0 (real) / +4d (H3) / -4d (H5) / +8d (H7);
  seen   : in how many of the kept snapshots it was detected.

    .venv/bin/python bench/live_lo/triage.py 472 [--all]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlap import mean_psd
from harmonic_map import HYP

HERE = Path(__file__).resolve().parent
SELF_T = 12.4e-3
CLOCKS = (61.44e6, 40e6)
PLUTO_DDR = 533.333e6           # Zynq DDR3 clock (own oscillator: ~+4 ppm vs our display), with ±k x 1.92 MHz lines
# site EMI clocks measured on this bench (FINDINGS.md): (fundamental, offset vs our raw display in ppm)
EMI = ((12e6, 7.0), (40e6, -4.5))

def on_clock(f):
    if abs(f - PLUTO_DDR) < 10e3 or (abs(f - PLUTO_DDR) < 12e6 and abs(((f - PLUTO_DDR) / 1.92e6) % 1 - 0.5) > 0.49):
        return True
    if abs(f - round(f / 1e6) * 1e6) < 300:  # an external line here would be ~4 kHz off (our +7.4 ppm): own reference
        return True
    ps = 100e6 / 3 * (1 + 5e-6)              # Zynq PS clock 33.333 MHz harmonics (DDR = 16x)
    if abs(f - round(f / ps) * ps) < 1.5e3:
        return True
    return any(abs(f - round(f / c) * c) < 1.5e3 for c in CLOCKS)

def on_emi(f):
    for c, ppm in EMI:
        n = round(f / c)
        if abs(f - n * c * (1 + ppm * 1e-6)) < 1.5e3:
            return f"EMI {n}x{c/1e6:g}MHz"
    return None

def capture_lock(cap, fc, bw):
    """In-band power, first 12 ms minus after 13 ms, averaged over snapshots (dB). Air does not know our buffer."""
    fs, lo = cap["rate"], cap["freq_mhz"] * 1e6; m = 4096; d, cs = [], []
    for s in cap["snapshots"]:
        a = np.load(s["file"], mmap_mode="r")
        def psd(i0, i1):
            x = a[i0:i1, 0].astype(np.float32) + 1j * a[i0:i1, 1].astype(np.float32)
            X = np.abs(np.fft.fftshift(np.fft.fft(x[:len(x) // m * m].reshape(-1, m) * np.hanning(m), axis=1), axes=1)) ** 2
            f = (np.arange(m) - m // 2) * fs / m + lo; sel = np.abs(f - fc) <= max(bw / 2, 2 * fs / m)
            return 10 * np.log10(X[:, sel].mean(0) + 1e-12)
        p1, p2 = psd(int(0.6e-3 * fs), int(12e-3 * fs)), psd(int(13e-3 * fs), int(24.4e-3 * fs))
        d.append(10 * np.log10(np.mean(10 ** (p1 / 10)) / np.mean(10 ** (p2 / 10))))
        if len(p1) >= 8 and max(p1.std(), p2.std()) > 1.0:
            cs.append(np.corrcoef(p1, p2)[0, 1])
    return float(np.mean(d)), (float(np.mean(cs)) if cs else None)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    out = HERE / "out" / a.tag
    snaps = sorted(out.glob("snap*.json"))
    cap = json.loads((out / "capture.json").read_text()); lo = cap["freq_mhz"] * 1e6
    dtag = next((p.name for p in (HERE / "out").iterdir() if (p / "capture.json").exists()
                 and json.loads((p / "capture.json").read_text()).get("dither_of") == cap["freq_mhz"]), None)
    fa, ca, pa = mean_psd(a.tag, m=1 << 14); pa -= np.median(pa)
    fine_f, _, fine_p = mean_psd(a.tag, m=1 << 17); fine_f = fine_f + ca
    fine_p = fine_p - np.median(fine_p)
    if dtag:
        fb, cb, pb = mean_psd(dtag, m=1 << 14); pb -= np.median(pb); delta = cb - ca
        fbf, _, pbf = mean_psd(dtag, m=1 << 17); pbf -= np.median(pbf)
    groups = []
    for k, p in enumerate(snaps):
        for s in json.loads(p.read_text())["signals"]:
            f = s["frequency"]; t = s.get("time", {})
            fc, bw = f["center_hz"], max(f.get("bandwidth_hz") or 0, 10e3)
            g = next((g for g in groups if abs(g["fc"] - fc) < max(bw, g["bw"]) / 2 + 10e3), None)
            rec = dict(snap=k, t0=t.get("start_s"), t1=t.get("end_s"), label=s["label"], tech=s["technology"],
                       snr=s.get("snr_db"), dec=s["decode"]["status"])
            if g is None:
                groups.append(dict(fc=fc, bw=bw, recs=[rec]))
            else:
                g["recs"].append(rec)
    print(f"LO {a.tag} MHz: {len(groups)} detection groups (dither {dtag})")
    rows = []
    for g in sorted(groups, key=lambda g: g["fc"]):
        fc, bw, recs = g["fc"], g["bw"], g["recs"]
        seen = len({r["snap"] for r in recs})
        in_win = all(r["t1"] is not None and r["t1"] <= SELF_T for r in recs)
        sel = np.abs(fine_f - fc) < max(bw / 2, 5e3); raw = fine_f[sel][np.argmax(fine_p[sel])] if sel.any() else fc
        clock = on_clock(raw)
        img = None
        if dtag:
            narrow = bw < 100e3                      # narrow: fine PSD, small window (keeps broadband images out)
            W = 40e3 if narrow else max(1.5 * bw, 300e3); step = 500.0 if narrow else 10e3
            FA, PA, FB, PB = ((fine_f, fine_p, fbf + cb, pbf) if narrow else (ca + fa, pa, cb + fb, pb))
            grid = np.arange(fc - W, fc + W, step)
            A = np.interp(grid, FA, PA)
            if A.std() > 0.3:
                sc = {k: np.corrcoef(A, np.interp(grid + k * delta, FB, PB))[0, 1] for k in HYP}
                k = max(sc, key=sc.get); n, sgn = HYP[k]
                if n != 1 and sc[k] > 0.5 and sc[k] - sc[0] > 0.15:
                    d = fc - ca; img = (n, (n * ca + d) if sgn > 0 else (n * ca - d), sc[k])
        emi = on_emi(raw) if bw < 100e3 else None
        lock, shape = capture_lock(cap, fc, bw)
        locked = abs(lock) >= 6 or (shape is not None and shape < 0.3)
        verdict = ("SELF-window" if in_win else "SELF-clock" if clock else
                   f"SELF-locked({lock:+.0f}dB,r={shape if shape is None else round(shape, 2)})" if locked else emi if emi else
                   f"IMAGE H{img[0]} of {img[1]/1e6:.1f}" if img else "REAL")
        best = max(recs, key=lambda r: (r["dec"] == "decoded", r["snr"] or -99))
        rows.append((verdict, fc, bw, seen, best))
        if a.all or verdict == "REAL":
            print(f"  {fc/1e6:10.4f} MHz {bw/1e3:7.1f} kHz  seen {seen}/{len(snaps)}  snr {best['snr']}  {verdict:22s} "
                  f"{best['label']} [{best['dec']}]  raw {raw/1e6:.5f}")
    from collections import Counter
    print("  verdicts:", dict(Counter(r[0].split(' ')[0] for r in rows)))

if __name__ == "__main__":
    main()
