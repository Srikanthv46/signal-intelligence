"""Reference vectors for tests/unit/test_rx_rc_link.cpp, written by the Python rc_link receiver.

Run: .venv/bin/python c++/tests/vectors/rc_link/gen_rc_link.py
Inputs are the testdata captures (read by the C++ test from SI_TESTDATA); only Python's
intermediate results are stored here (raw little-endian arrays + JSON).
"""
import json
from pathlib import Path

import numpy as np
from scipy import fft as sfft

from spectrum_intel.core.emission import load_sigmf
from spectrum_intel.receivers.rc_link import detector as D
from spectrum_intel.receivers.rc_link import lora_phy as P
from spectrum_intel.receivers.rc_link.label_check import check_frequency_label

OUT = Path(__file__).parent
REPO = Path(__file__).resolve().parents[4]


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


ref = {}
# ------------------------------------------------------------------ LoRa: real SF9 packet at 434 MHz
x = np.fromfile(REPO / "testdata/lora/real_434MHz_sf9.cf32", np.complex64)
fs, fc = 500e3, 433.92e6
X = sfft.fft(x)
D._upchirp(512).tofile(OUT / "upchirp_512.bin")
y, ofs, fa = D._slice(X, fs, 87500.0, 62500.0, 125e3)
y.tofile(OUT / "slice_plain.bin")
y2, _, _ = D._slice(X, fs, np.float32(87512.34), 62500.0, 125e3, t_adv=0.37)
y2.tofile(OUT / "slice_f32_tadv.bin")
ref["slice"] = dict(ofs=ofs, f_act=fa, n=len(y))
zu, zd, _, _ = D._streams(X, fs, 87500.0, 125e3, 512, 1, t_adv=0.37)
zu.tofile(OUT / "streams_zup.bin")
ss = []
for st in (0, 100, 4096, 20000, len(zu) - 512):
    S, pk, r, s2 = D._sym_spec(zu, st, 512)
    ss.append(dict(start=st, pk=float(pk), r=float(r), s2=float(s2), argmax=int(np.argmax(S)), S_sum=float(np.sum(S.astype(np.float64)))))
ref["sym_spec"] = ss
rng = np.random.default_rng(3)
v = (rng.random(6) * 512).astype(np.float32).astype(float).tolist(); w = (rng.random(6) * 50).astype(np.float32).astype(float).tolist()
ref["circ_mean"] = dict(v=v, w=w, out=float(D._circ_mean(v, 512, w)))
cands, _ = D.lora_scan(X, fs, fc, 125e3, [7, 8, 9, 10, 11, 12])
ref["lora_scan_n"] = len(cands)
cands = D._nms_candidates(cands, 125e3)
ref["lora_nms"] = clean(cands)
# capture the walk inputs of every refined candidate
walk_in = []
orig_walk = D._walk_packet
def spy(z_up, z_dn, b_sfd, N, cfo, pol, **kw):
    walk_in.append((z_up, z_dn, float(b_sfd), int(N), float(cfo), int(pol)))
    return orig_walk(z_up, z_dn, b_sfd, N, cfo, pol, **kw)
D._walk_packet = spy
refs = []
for c in cands:
    walk_in.clear()
    w = D.lora_refine(X, fs, fc, c)
    if w is None:
        refs.append(None)
        continue
    pr = w if w["partial"] else w["prem"]
    refs.append(dict(partial=bool(w["partial"]), n_pre_est=int(pr["n_pre_est"]), pre_start=int(pr["pre_start"]),
                     pre_end=int(pr["pre_end"]), k_u=float(pr["k_u"]),
                     walk=None if w["partial"] else clean({k: v for k, v in w.items() if k != "prem"})))
D._walk_packet = orig_walk
ref["lora_refine"] = refs
# one walk input pair (the last walk of the last full candidate) for a direct walk_packet test
for c in cands:
    walk_in.clear()
    D._walk_packet = spy
    w = D.lora_refine(X, fs, fc, c)
    D._walk_packet = orig_walk
    if w is not None and not w["partial"]:
        zu, zd, b, N, cfo, pol = walk_in[-1]
        zu.tofile(OUT / "walk_zup.bin"); zd.tofile(OUT / "walk_zdn.bin")
        ref["walk_direct"] = dict(n=len(zu), b=b, N=N, cfo=cfo, pol=pol, out=clean(orig_walk(zu, zd, b, N, cfo, pol)))
        break
ems = D.detect_lora(X, fs, fc)
ref["lora_emissions"] = [clean(e.to_dict()) for e in ems]
# ------------------------------------------------------------------ lora_phy
ref["lora_phy"] = []
for e in ems:
    if e.decoded:
        ref["lora_phy"].append(dict(sf=e.features["sf"], bw=e.features["bw_hz"], decoded=clean(e.decoded)))
ref["whitening"] = P.whitening_sequence(16).tolist()
ref["crc16"] = int(P.crc16(b"123456789"))
# ------------------------------------------------------------------ FSK + BLE: Phantom 4 capture
xp, fsp, fcp, _ = load_sigmf(str(REPO / "testdata/DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta"))
xp = xp.astype(np.complex64)
fcands = D.fsk_scan(xp, fsp, fcp)
ref["fsk_cands"] = clean(fcands)
meas = []
for c in fcands[:40]:
    for r in D._refine_local(xp, fsp, fcp, c):
        meas.append(clean({k: v for k, v in r.items()}))
ref["fsk_meas"] = meas
fe = D.detect_fsk(xp, None, fsp, fcp)
ref["fsk_emissions"] = [clean(e.to_dict()) for e in fe]
ref["label_check"] = clean(check_frequency_label(fe))
json.dump(ref, open(OUT / "ref.json", "w"), indent=0)
print("ok", len(cands), len(fcands), len(meas), len(fe))
