"""Reference for tests/unit/test_rx_dji.cpp: receivers.dji building blocks on small synthetic
inputs, a synthetic DroneID burst (DroneID never fires on testdata), and detect() on the three
DJI testdata captures (OcuSync 2 video with ZC slots, O4 30 kHz video, RC uplink hops).

Run: nice -n 19 .venv/bin/python c++/tests/vectors/dji/gen_dji.py"""
import json
from pathlib import Path

import numpy as np

from spectrum_intel.io import sigmf
from spectrum_intel.receivers.dji import detector as D

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
rng = np.random.default_rng(7)


def cplx(v):
    return [[float(z.real), float(z.imag)] for z in v]


def em_json(e):
    def val(v):
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        if isinstance(v, (np.integer, int)):
            return int(v)
        if isinstance(v, (np.floating, float)):
            return float(v)
        if isinstance(v, (list, tuple)):
            return [val(u) for u in v]
        return v
    return dict(family=e.family, label=e.label, start_us=float(e.start_us), end_us=float(e.end_us),
                f_lo_mhz=float(e.f_lo_mhz), f_hi_mhz=float(e.f_hi_mhz), center_mhz=float(e.center_mhz),
                snr_db=None if e.snr_db is None else float(e.snr_db), confidence=float(e.confidence),
                features={k: val(v) for k, v in e.features.items()})


ref = {}
# ---- ZC helpers
z = D.zc_freq(137, 601)
ref["zc_freq_137_601"] = cplx(np.r_[z[:4], z[298:303], z[-3:]])
X73 = np.roll(D.zc_freq(50, 73), 5) * np.exp(1j * 0.7) + 0.3 * (rng.standard_normal(73) + 1j * rng.standard_normal(73))
cc, r = D.zc_root_search(X73, 73)
ref["zc_search_73"] = dict(X=cplx(X73), cc=cc, root=r)
cc1, r1 = D.zc_root_search(X73, 73, [50])
ref["zc_search_73_fixed"] = dict(cc=cc1, root=r1)
# DD test on a 1024-bin symbol carrying root 192 (L=601) with a timing offset
N, L = 1024, 601
Xs = np.zeros(N, complex)
Xs[N // 2 - 300:N // 2 + 301] = D.zc_freq(192, L) * np.exp(-2j * np.pi * 3 * np.arange(L) / N)
Xs += 0.2 * (rng.standard_normal(N) + 1j * rng.standard_normal(N))
acc, mag = D._zc_dd(Xs, L)
ref["zc_dd"] = dict(X=cplx(Xs), acc=[float(acc.real), float(acc.imag)], mag=float(mag), root=D._root_from_dd(acc, L))
# ---- CP metric / fold / robust z / segmentation on a synthetic grid
P = D.SLOT_S * 15.36e6 / 7
n = 40000
y = 0.5 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
for k in range(8, 30):                                     # symbols with a CP copy, 1 burst
    s = int(round(300 + k * P))
    y[s:s + 72] += 2 * y[s + 1024:s + 1096] / 2 + 1.5 * np.exp(1j * rng.uniform(0, 6.3, 72))
    y[s + 1024:s + 1096] = y[s:s + 72]
m, cs = D._cp_metric(y, 1024, 72)
fz = D._robust_z(D._fold(m, P))
runs, t, on, off = D._segment(m, int(np.argmax(fz)), P)
(HERE / "cp_y.c128").write_bytes(y.astype(np.complex128).tobytes())
ref["cp"] = dict(n=n, m_head=m[:5].tolist(), m_sum=float(np.sum(m)), cs_last=[float(cs[-1].real), float(cs[-1].imag)],
                 fold_argmax=int(np.argmax(fz)), fold_max=float(fz.max()), runs=[list(map(int, r_)) for r_ in runs],
                 t=t.tolist()[:5], on=on.tolist(), off=off.tolist())
# ---- frame period
ref["frame_period"] = [dict(s=s_, T=(None if D._frame_period(s_) is None else float(D._frame_period(s_))))
                       for s_ in ([0.1, 4.1, 10.1, 14.1, 20.1], [0.0, 1.0, 2.0, 3.0], [0.2, 2.7, 5.2, 12.7, 15.2],
                                  [0.0, 1.5, 7.0, 8.5, 9.0, 16.0, 20.0])]
# ---- synthetic DroneID burst at 15.36 MSps (9 symbols, ZC 600 / 147 in symbols 4 and 6)
fs = 15.36e6
fc = 2414.5e6
sym = []
for i in range(9):
    F = np.zeros(1024, complex)
    if i in (3, 5):
        F[512 - 300:512 + 301] = D.zc_freq(600 if i == 3 else 147, 601)
    else:
        F[512 - 300:512 + 301] = np.exp(1j * np.pi / 2 * rng.integers(0, 4, 601)) / np.sqrt(1)
        F[512] = 0
    td = np.fft.ifft(np.fft.ifftshift(F)) * 32
    sym.append(np.r_[td[-D.DID_CP[i]:], td])
burst = np.concatenate(sym)
x = 0.05 * (rng.standard_normal(65536) + 1j * rng.standard_normal(65536))
s0 = 20000
x[s0:s0 + len(burst)] += burst * np.exp(2j * np.pi * (1.1e6 + 2.3e3) * np.arange(len(burst)) / fs)
x = x.astype(np.complex64)
(HERE / "droneid_cf32.bin").write_bytes(x.tobytes())
ch = D._FastChannelizer(D.C.Channelizer(x.astype(np.complex128), fs, []))
sc = D._droneid_cp_score(ch.slice(-0.18e6, 5.0e6, fs))
ref["droneid"] = dict(fs=fs, fc=fc, score_argmax=int(np.argmax(sc)), score_max=float(sc.max()),
                      emissions=[em_json(e) for e in D.detect(x, fs, fc, {"decode_droneid": False})])
# ---- the three DJI testdata captures
caps = {}
for rel in ("DJI_Mini_2_20251119-165139/D_20251119-165139", "DJI_Mini_Pro_4_20251119-172214/D_20251119-172214",
            "D_20260916-171343/D_20260916-171343"):
    cap = sigmf.load(REPO / "testdata" / (rel + ".sigmf-meta"))
    caps[rel] = [em_json(e) for e in D.detect(cap.iq, cap.fs, cap.fc, {"decode_droneid": False})]
ref["captures"] = caps
(HERE / "ref.json").write_text(json.dumps(ref))
print({k: len(v) for k, v in caps.items()}, len(ref["droneid"]["emissions"]))
