#!/usr/bin/env python3
"""Synthetic reference for BLE LE Coded PHY (BT5 Long Range Remote ID; rc_link/ble_coded; research, PLAN T2b).
  cases.json  viterbi on noisy coded bits, parse_extended, and per burst (8 MSps, bursts.cf32): find_coded, decode,
              ble_coded_decode;
  capture.cf32 / capture.json  a 15.36 MSps capture (fc 2406 MHz) with an ADV_EXT_IND on ch 37 (S=8) and an
              AUX_ADV_IND carrying an F3411 message pack on data ch 3 (S=8) at 3 dB, and rc_link detect's
              decoded dicts on it (the whole pipeline: 500 kBd-apparent gate, whole-box decode).

    nice -n 19 .venv/bin/python -I c++/tests/vectors/ble_coded/gen_ble_coded.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
for p in ("python", "python/tests", "python/tests/receivers"):
    sys.path.insert(0, str(REPO / p))
from test_det_droneid import odid_ie_body                                    # noqa: E402
from spectrum_intel.receivers.rc_link import ble_coded as BC, detector as D, gfsk_coherent as GC, synth as S   # noqa

AA = np.array([(D.BLE_ADV_AA >> i) & 1 for i in range(32)], np.uint8)
ADDR = bytes([0x21, 0x43, 0x65, 0x87, 0xA9, 0xCB])
rng = np.random.default_rng(9)
ref = {}

b = list(rng.integers(0, 2, 120)) + [0, 0, 0]
c = 2.0 * np.array(BC.conv_encode(b)) - 1 + 0.9 * rng.standard_normal(2 * len(b))
ref["viterbi"] = dict(llr=c.tolist(), out="".join(map(str, BC.viterbi(c.reshape(-1, 2)).tolist())))
ref["known0"] = "".join(map(str, BC.known_symbols(AA, 0).tolist()))
pe = S.ext_adv_payload(ADDR, b"\x02\x01\x06", did=0x456, sid=3, aux=(8, 166, 2))
ref["parse_extended"] = dict(hex=pe.hex(), out={k: (v.hex() if isinstance(v, bytes) else v)
                                                for k, v in BC.parse_extended(pe).items()})

body = bytes([0]) + odid_ie_body()[1:]
aux = S.ext_adv_payload(ADDR, S.remote_id_adv_data(body), did=0x456, sid=3)
prim = S.ext_adv_payload(ADDR, b"", did=0x456, sid=3, aux=(3, 166, 2))
cases, chunks, at = [], [], 0
for snr, payload, ch, ci in ((8, prim, 37, 0), (2, aux, 3, 0), (6, aux, 3, 1)):
    iq, _ = S.ble_coded_burst(payload, ch=ch, ci=ci, cfo_hz=rng.uniform(-30e3, 30e3), fs=8e6, rng=rng)
    x = np.r_[np.zeros(300), iq, np.zeros(300)]
    n0 = 1 / (1e6 * 10 ** (snr / 10))
    x = x + (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x))) * np.sqrt(n0 * 8e6 / 2)
    X = np.fft.fft(x); f = np.fft.fftfreq(len(x), 1 / 8e6); X[np.abs(f) > 0.75e6] = 0
    y = np.fft.ifft(X).astype(np.complex64)
    fc_ = BC.find_coded(GC.diff_products(y), AA)
    r = BC.decode(y, fc_[0][0], AA) if fc_ else None
    cases.append(dict(offset=at, n=len(y), ch=ch, find=[[int(i), float(v)] for i, v in fc_],
                      ci=r["ci"] if r else None, aa_err=r["aa_bit_errors"] if r else None, h=r["h"] if r else None,
                      bits="".join(map(str, r["bits"].tolist())) if r and "bits" in r else None,
                      out=D.ble_coded_decode(y, 8e6, ch, AA)))
    chunks.append(y); at += len(y)
    print(snr, ch, ci, "decoded", bool(cases[-1]["out"]))
np.concatenate(chunks).astype(np.complex64).tofile(OUT / "bursts.cf32")
ref["cases"] = cases

# whole pipeline on a 15.36 MSps capture
FS, FC = 15.36e6, 2406e6
x = S.noise_capture(1 << 18, rng)
ip, _ = S.ble_coded_burst(prim, ch=37, cfo_hz=9e3, fs=FS, rng=rng)
ia, _ = S.ble_coded_burst(aux, ch=3, cfo_hz=9e3, fs=FS, rng=rng)
y, _ = S.embed(x, FS, FC, [dict(t0_us=1000, f_hz=2402e6, iq=ip, bw_hz=1e6),
                           dict(t0_us=5000, f_hz=2410e6, iq=ia, bw_hz=1e6)], 3)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = D.detect(y, FS, FC)
ref["capture"] = dict(fs=FS, fc=FC, decoded=[e.decoded for e in em if (e.decoded or {}).get("phy")])
print("pipeline:", [(d["pdu_name"], d["channel_index_rf"]) for d in ref["capture"]["decoded"]])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
