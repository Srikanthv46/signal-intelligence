#!/usr/bin/env python3
"""Synthetic reference for the coherent BLE GFSK detector + CRC Chase (rc_link/gfsk_coherent; research, PLAN L1).
  bursts.cf32   BLE 4 Remote ID adverts at 8 MSps (one after another; offsets in cases.json), 6-10 dB in-band SNR,
                nominal / off-nominal modulation index, with and without drift;
  cases.json    per burst: find_access_address, demodulate (bits, h, cfo), ble_coherent_decode and ble_try_decode;
                plus CRC syndrome / Chase cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/ble_coherent/gen_ble_coherent.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
sys.path.insert(0, str(REPO / "python"))
sys.path.insert(0, str(REPO / "python" / "tests"))
sys.path.insert(0, str(REPO / "python" / "tests" / "receivers"))
from test_gfsk_coherent import ADDR, AA, basic_id_advert, location_advert, _burst     # noqa: E402
from spectrum_intel.receivers.rc_link import detector as D, gfsk_coherent as GC, synth as S   # noqa: E402

rng = np.random.default_rng(42)
cases, chunks, at = [], [], 0
for snr, dev, drift, adv in ((10, 250e3, 0, basic_id_advert(0)), (8, 250e3, 0, location_advert(1)),
                             (6, 250e3, 0, basic_id_advert(2)), (8, 228e3, 100, location_advert(3)),
                             (7, 270e3, 0, basic_id_advert(4))):
    bits = np.array(S.ble_adv_bits(ADDR, adv, ch=37), np.uint8)
    y = _burst(snr, rng, bits, dev=dev, drift=drift).astype(np.complex64)
    aa = GC.find_access_address(GC.diff_products(y), AA)
    i0 = aa[0][0]
    nb = min(8, i0 // 8)
    dm = GC.demodulate(y, i0 - nb * 8, nb + 32 + 16 + 8 * 255 + 24)
    c = dict(offset=at, n=len(y), snr=snr, dev=dev, drift=drift, aa=[[int(i), float(v)] for i, v in aa],
             s0=int(i0 - nb * 8), bits="".join(map(str, dm[0].tolist())), h=dm[2], cfo_hz=dm[3], dt=dm[4],
             truth_errors=int(np.count_nonzero(dm[0][:len(bits) - 2] != bits[:min(len(dm[0]), len(bits) - 2)])),
             coherent=D.ble_coherent_decode(y, 8e6, 37, AA), try_decode=D.ble_try_decode(y, 8e6, 2402e6))
    cases.append(c)
    chunks.append(y)
    at += len(y)
    print(snr, dev, drift, "h %.3f" % dm[2], "errors", c["truth_errors"], "decoded", bool(c["coherent"]),
          "fixes", (c["coherent"] or {}).get("bit_corrections"))
np.concatenate(chunks).astype(np.complex64).tofile(OUT / "bursts.cf32")

# CRC syndrome / Chase
bits = np.array(S.ble_adv_bits(ADDR, location_advert(), ch=37), np.uint8)
d = D._ble_dewhiten(bits[40:], 37)
n_pdu = 16 + 8 * int(sum(int(d[8 + i]) << i for i in range(8)))
syn, contrib = GC.crc_syndrome_parts(d, n_pdu)
rel = np.linspace(2, 3, 40 + n_pdu + 24)
chase = []
for flips in ([], [100], [100, 101], [30, 31, 200, 201], [5, 6, 77, 78, 150, 151], [12]):
    e = d.copy()
    for q in flips:
        e[q] ^= 1
    r = rel.copy()
    for q in flips:
        r[40 + q] = 0.01 * (1 + q % 7)
    got = GC.chase(e, n_pdu, r, 40, protect=range(8, 16))
    chase.append(dict(bits="".join(map(str, e[:n_pdu + 24].tolist())), rel=r.tolist(),
                      syndrome=GC.crc_syndrome_parts(e, n_pdu)[0],
                      out=None if got is None else ["".join(map(str, got[0][:n_pdu + 24].tolist())), got[1]]))
(OUT / "cases.json").write_text(json.dumps(dict(cases=cases, n_pdu=n_pdu, contrib=contrib.tolist(), chase=chase),
                                           indent=1, default=float))
