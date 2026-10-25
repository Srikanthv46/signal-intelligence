"""Reference for tests/unit/test_rx_digital_vtx.cpp: receivers.digital_vtx building blocks and
detect() on a small synthetic capture (20 ms at 7.68 MSps): a 5 ms-periodic LTE PSS (root 29) at
+1.2 MHz (Herelink path) and 2.5 ms-periodic bursts of a 2 MHz CP-OFDM link (128-point, Tu 16.7 us,
CP 1/8) at -2 MHz (generic blind-link path, blind.py), in complex Gaussian noise.

Run: nice -n 19 .venv/bin/python c++/tests/vectors/digital_vtx/gen_digital_vtx.py"""
import json
from pathlib import Path

import numpy as np

from scipy.signal import resample_poly

from spectrum_intel.receivers.digital_vtx import blind as BL, catalog as CAT, detector as D, hdzero as HZ, streams as ST

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(11)
FS, FC, L = 7.68e6, 2.44e9, 153600


def cplx(v):
    return [[float(z.real), float(z.imag)] for z in v]


def val(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (np.integer, int)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return float(v)
    if isinstance(v, (list, tuple)):
        return [val(x) for x in v]
    return v


def em_json(e):
    return dict(family=e.family, label=e.label, start_us=val(e.start_us), end_us=val(e.end_us),
                f_lo_mhz=val(e.f_lo_mhz), f_hi_mhz=val(e.f_hi_mhz), center_mhz=val(e.center_mhz),
                snr_db=val(e.snr_db), confidence=val(e.confidence),
                features={k: val(v) for k, v in e.features.items()})


def capture():
    n = np.arange(L)
    x = (rng.standard_normal(L) + 1j * rng.standard_normal(L)) * 0.05
    pss = D._pss_time(29, FS)                                  # 512 samples, unit-ish power
    pss = pss / np.sqrt(np.mean(np.abs(pss) ** 2)) * 0.6
    for k in range(4):                                          # every 5 ms
        a = 3000 + k * 38400
        x[a:a + len(pss)] += pss * np.exp(2j * np.pi * 1.2e6 * n[a:a + len(pss)] / FS)
    # CP-OFDM bursts: 128-point IFFT (60 kHz SCS at 7.68 MSps), 16 used subcarriers each side (~2 MHz),
    # CP 16; a 0.86 ms burst every 2.5 ms
    used = np.r_[1:17, 112:128]
    for b0 in range(9000, L - 7000, 19200):
        syms = []
        for _ in range(46):
            S = np.zeros(128, complex)
            S[used] = (rng.choice([-1, 1], used.size) + 1j * rng.choice([-1, 1], used.size)) / np.sqrt(2)
            s = np.fft.ifft(S) * 6
            syms.append(np.r_[s[-16:], s])
        burst = np.concatenate(syms)
        burst = np.repeat(burst, 1)
        m = n[b0:b0 + len(burst)]
        x[b0:b0 + len(burst)] += burst * np.exp(-2j * np.pi * 2.0e6 * m / FS)
    return x.astype(np.complex64)


x = capture()
(HERE / "synth.cf32").write_bytes(x.tobytes())
ref = {"fs": FS, "fc": FC}
ref["pss_time_29_192"] = cplx(D._pss_time(29, 1.92e6))
ref["pss_time_34_384"] = cplx(D._pss_time(34, 3.84e6))
X = np.fft.fft(x.astype(np.complex128))
y, ofs = D._slice(X, FS, 1.2345e6, 0.5e6, 0.96e6)
ref["slice"] = dict(y=cplx(y[:64]), n=len(y), ofs=ofs, sum=cplx([y.sum()]))
sr = D.pss_fold_search(X, FS, 3.0e6)
ref["fold"] = dict(offs=val(list(sr["offs"])), Z={str(u): val(list(z)) for u, z in sr["Z"].items()},
                   T={str(u): val(list(t)) for u, t in sr["T"].items()}, ofs=sr["ofs"], M=sr["M"],
                   P=sr["P"], nt=sr["nt"], global_median=sr["global_median"])
c = D._pss_candidates(sr, 7.0)[0]
t_inst = [(c["tau"] + k * sr["P"]) / sr["ofs"] for k in range(6)]
t_inst = [t for t in t_inst if t + 70e-6 < L / FS]
f_c, u, s = D._carrier_search(X, FS, c["f_off"], t_inst)
ref["carrier"] = dict(f0=c["f_off"], t_inst=t_inst, f_c=float(f_c), root=int(u), score=float(s))
sl = lambda f, h, r: D._slice(X, FS, f, h, r)
usable = FS / 2
ref["rate_for"] = [list(BL.rate_for(t, b)) for t, b in ((66.667e-6, 18e6), (33.3335e-6, 9e6), (16.6667e-6, 1e6))]
hits = sorted(BL.scan(sl, usable), key=lambda h: (-h["z"], h["W"], h["off"], h["N"]))
ref["scan"] = [{k: val(v) for k, v in h.items()} for h in hits[:12]]
ref["n_scan"] = len(hits)
h = next(h for h in hits if h["tu"] >= BL.TU_MIN)
cf = BL.confirm(sl, h["off"], h["W"], h["N"], h["ofs"])
ref["confirm"] = dict(off=h["off"], W=h["W"], N=h["N"], ofs=h["ofs"], result=val(list(cf)) if cf else None)
tu, ts = cf[0], cf[1]
runs, tk, zr, ts2, ng, frac = BL.segment(sl, h["off"], h["W"], tu, ts)
ref["segment"] = dict(runs=val([list(r) for r in runs]), n_tk=len(tk), tk=val(list(tk[:40])), zr=val(list(zr[:40])),
                      zr_sum=float(np.sum(zr)), ts=float(ts2), ngrid=ng, frac=float(frac))
T, con = BL.frame_period(tk, zr, ts2)
ref["frame_period"] = [val(T), val(con)]
ref["frame_complete"] = val([list(r) for r in BL.frame_complete(runs, tk, zr, ts2, T)]) if T else None
fr, c = BL.band_cp(sl, h["off"], h["W"], tu, ts2, tk[zr > 2])
ref["band_cp"] = dict(n=len(c), fr=val(list(fr[::8])), c=val(list(c[::8])))
ref["measure_band"] = val(list(BL.measure_band(sl, h["off"], h["W"], tu, ts2, tk[zr > 2], usable)))
noise = float(D.C.noise_floor_psd(x)["floor_power_per_sample"])
ref["links"] = [{k: val(v) if not isinstance(v, list) else [val(list(r)) for r in v] for k, v in lk.items()}
                for lk in BL.blind_links(sl, FS, FC, usable, noise)]
ref["classify"] = {f"{t}|{b}": val(list(D.classify_numerology(t, None, b if b > 0 else None)))
                   for t, b in ((12.8, 0), (66.6, 0), (91.4, 9.6), (91.4, 18.0), (33.3, 9.0), (22.86, 0), (16.7, 0), (37.9, 0))}
ctx = {"report_wifi": False}
ref["detect"] = [em_json(e) for e in D.detect(x, FS, FC, ctx)]
# transmit streams (streams.tx_streams) on two small link slices: 512-point, Tu 91.43 us, CP 1/8, 448 subcarriers,
# 3 bursts of 26 symbols, 2x 6-tap channels, 30 dB SNR: one stream (64-QAM) and Alamouti STBC (16-QAM)
def link_slice(kind, M, g):
    N, CP, OFS = 512, 64, 5.6e6
    used = np.r_[-224:0, 1:225]
    def qam(n):
        m = int(np.sqrt(M)); lv = np.arange(-(m - 1), m, 2, dtype=float)
        return (g.choice(lv, n) + 1j * g.choice(lv, n)) / np.sqrt(np.mean(lv ** 2) * 2)
    def chan():
        h = (g.standard_normal(6) + 1j * g.standard_normal(6)) * np.exp(-np.arange(6) / 2)
        return np.fft.fft(h / np.sqrt(np.sum(np.abs(h) ** 2)) / np.sqrt(2), N)[used % N]
    H1, H2 = chan(), chan()
    out, runs, t = [np.zeros(900, complex)], [], 900
    for _ in range(3):
        S = qam(26 * used.size).reshape(26, -1)
        if kind == "stbc":
            F = np.empty_like(S)
            F[0::2] = H1 * S[0::2] + H2 * S[1::2]
            F[1::2] = -H1 * np.conj(S[1::2]) + H2 * np.conj(S[0::2])
        else:
            F = H1 * S
        G = np.zeros((26, N), complex)
        G[:, used % N] = F
        xx = np.fft.ifft(G, axis=1)
        xx = np.hstack([xx[:, -CP:], xx]).ravel()
        runs.append((t / OFS, (t + len(xx)) / OFS, 26))
        out += [xx, np.zeros(900, complex)]
        t += len(xx) + 900
    yy = np.concatenate(out)
    pw = np.mean(np.abs(yy[np.abs(yy) > 0]) ** 2)
    yy = yy + (g.standard_normal(len(yy)) + 1j * g.standard_normal(len(yy))) * np.sqrt(pw / 2e3)
    yy = (yy * np.exp(2j * np.pi * 2.3e3 / OFS * np.arange(len(yy)))).astype(np.complex64)
    return yy, runs, OFS, N / OFS, (N + CP) / OFS, 448 * OFS / N


g2 = np.random.default_rng(23)
ref["tx_streams"] = {}
for kind, M in (("1tx", 64), ("stbc", 16)):
    yy, runs, ofs_, tu_, ts_, bw_ = link_slice(kind, M, g2)
    (HERE / f"streams_{kind}.cf32").write_bytes(yy.tobytes())
    ref["tx_streams"][kind] = dict(ofs=ofs_, tu=tu_, ts=ts_, bw=bw_, runs=[list(r) for r in runs], dur=len(yy) / ofs_,
                                   out=ST.tx_streams(yy.astype(np.complex128), ofs_, tu_, ts_, bw_, runs, 30.0,
                                                     len(yy) / ofs_))
    print(kind, ref["tx_streams"][kind]["out"])
# research catalog (catalog.candidates) on measured-link tuples
# + the model level (catalog.model_candidates, research PLAN M3 / M4) on each result
ref["catalog"] = [dict(args=list(a), out=CAT.candidates(*a), model=CAT.model_candidates(CAT.candidates(*a), a[3], a[2], a[4]))
                  for a in (
    (91.432, 0.1253, 19.2, 5805.07, 2, 6.673), (33.33, 0.0703, 18.0, 5789.5, 2, 2.5), (33.33, 0.0703, 18.0, 5789.5, None, 2.5),
    (66.67, 0.0703, 18.0, 5804.5, 2, None), (91.43, 0.125, 9.6, 2426.6, None, None), (60.1, 0.125, 2.4, 2440.0, None, 2.5),
    (91.43, 0.125, 18.9, 5805.0, None, None), (66.67, 0.0703, 36.0, 5839.5, 2, None), (33.33, 0.0703, 36.0, 5200.0, 2, 2.5),
    (50.0, None, 17.0, 5806.0, None, None), (50.0, None, 27.0, 5806.0, None, None), (66.67, 0.0703, 18.0, 5768.5, 2, None),
    (37.926, 0.03125, 23.44, 5658.04, None, 16.661), (37.93, 0.0312, 17.2, 5399.0, 1, None))]
ref["hop_model"] = [[a, CAT.hop_model(*a)] for a in ((6.6, 410.0, 2.48), (7.0, 410.0, 2.0), (4.8, 940.8, 7.17),
                                                    (1.14, 513.0, 4.0), (5.0, 640.0, 4.0), (None, 1.0, 1.0),
                                                    (1.1, 510.0, None), (2.96, 520.0, 5.16), (6.06, 250.0, 2.97),
                                                    (3.32, 1550.0, 4.99), (4.6, 1610.0, 9.0), (5.5, 410.0, 2.2))]
# HDZero 27 MHz mode (hdzero.py): slice plans, a synthetic link (R1 + 20 kHz, ~10 dB, 30.72 MSps, 90 symbols with an
# 8-symbol preamble at symbol 30), the channel search and detect() on it, and the model block
g3 = np.random.default_rng(31)
FSZ, FCZ = 30.72e6, 5657e6
yz = resample_poly(HZ.synth(90, g3), 256, 225)
yz = yz * np.exp(2j * np.pi * (1e6 + 20e3) * np.arange(len(yz)) / FSZ)
pz = np.mean(np.abs(yz) ** 2) * FSZ / (889 * HZ.SCS)
xz = (yz + (g3.standard_normal(len(yz)) + 1j * g3.standard_normal(len(yz))) * np.sqrt(pz / 10 / 2)).astype(np.complex64)
(HERE / "hdzero.cf32").write_bytes(xz.tobytes())
XZ = np.fft.fft(xz.astype(np.complex128))
logz = []
lz = HZ.hdzero_links(lambda f, h, r: D._slice(XZ, FSZ, f, h, r), FSZ, FCZ, FSZ / 2, len(xz), log=logz.append)
ref["hdzero"] = dict(fs=FSZ, fc=FCZ, n=len(xz),
                     plans=[[n, f, *HZ.slice_plan(n, f)] for n, f in ((2097152, 61.44e6), (65536, 61.44e6),
                            (2000000, 61.44e6), (1998848, 61.44e6), (2000000, 56e6), (len(xz), FSZ), (100003, 20e6))],
                     links=[{k: val(v) for k, v in d.items()} for d in lz], log=[{k: val(v) for k, v in d.items()} for d in logz],
                     detect=[em_json(e) for e in D.detect(xz, FSZ, FCZ, {"report_wifi": False})],
                     model=CAT.hdzero_model("27 MHz (normal)"))
print("hdzero", ref["hdzero"]["log"], [(e["family"], e["label"]) for e in ref["hdzero"]["detect"]])
(HERE / "ref.json").write_text(json.dumps(ref))
print(len(ref["detect"]), "emissions;", {(e["family"], e["label"]) for e in ref["detect"]}, "links",
      [(l["tu"], l["cen"], l["bw"], len(l["runs"]), l["frame_s"]) for l in ref["links"]])
